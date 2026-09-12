"""Did AlphaFold already see a homolog of each target at training time?

The deposit-date gate only guarantees that THIS entry was not in the training
set. It says nothing about homologs. If a close relative was already in the PDB
before 2021-09-30 (the AlphaFold 3 training cutoff), the model has effectively
seen the fold, and predicting the monomer becomes template recall rather than
prediction.

Method: for every protein chain, run an MMseqs2 sequence search against the PDB
restricted to entries RELEASED ON OR BEFORE 2021-09-30 -- i.e. exactly the set
AlphaFold 3 could have trained on -- and record the best sequence identity.

Interpretation:
  >= 95%   essentially the same protein was already known: monomer fold is
           memorised, only the assembly/ligand/PTM aspect is a real test
  70-95%   close homolog: fold certain, details (loops, side chains) open
  30-70%   distant homolog: fold likely right, accuracy genuinely uncertain
  < 30%    no useful homolog: a real fold-prediction test (rare)

Note this is measured per CHAIN. A known monomer fold does not mean the COMPLEX
is known -- which matters here, because eight of the nine classes are about
interactions, not monomer folds. The pair-level check is done separately below.
"""
import json
import os
import time
import urllib.request

OP = urllib.request.build_opener(urllib.request.ProxyHandler({}))
SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
TRAIN_CUTOFF = "2021-09-30"


def post(payload, timeout=300):
    req = urllib.request.Request(
        SEARCH, data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json"})
    return json.load(OP.open(req, timeout=timeout))


def homolog_search(seq, cutoff_date=TRAIN_CUTOFF, identity=0.15, rows=25):
    """Best pre-cutoff PDB match for this sequence."""
    q = {
        "query": {
            "type": "group", "logical_operator": "and",
            "nodes": [
                {"type": "terminal", "service": "sequence",
                 "parameters": {"evalue_cutoff": 10, "identity_cutoff": identity,
                                "sequence_type": "protein", "value": seq}},
                {"type": "terminal", "service": "text",
                 "parameters": {"attribute": "rcsb_accession_info.initial_release_date",
                                "operator": "less_or_equal", "value": cutoff_date}},
            ]},
        "return_type": "polymer_entity",
        "request_options": {"paginate": {"start": 0, "rows": rows},
                            "results_content_type": ["experimental"],
                            "scoring_strategy": "sequence"},
    }
    d = post(q)
    hits = []
    for x in d.get("result_set", []) or []:
        best = None
        for m in x.get("services", []) or []:
            for n in m.get("nodes", []) or []:
                for mm in n.get("match_context", []) or []:
                    si = mm.get("sequence_identity")
                    ev = mm.get("evalue")
                    if si is not None:
                        if best is None or si > best[0]:
                            best = (si, ev, mm.get("alignment_length"))
        if best:
            hits.append({"entity": x["identifier"], "identity": best[0],
                         "evalue": best[1], "alen": best[2]})
    hits.sort(key=lambda h: -h["identity"])
    return d.get("total_count", 0), hits


recs = json.load(open("targets41_full.json"))

cache = {}
if os.path.exists("homolog.json"):
    cache = json.load(open("homolog.json"))

jobs = []
for r in recs:
    for c in r["chains"]:
        if c["type"] != "Protein":
            continue
        # search with the sequence that is actually resolved in the structure
        seq = c.get("construct_seq") or c["seq"]
        if len(seq) < 12:
            continue          # too short for a meaningful sequence search
        key = f"{r['id']}_{c['entity']}"
        jobs.append((key, r, c, seq))

print(f"protein chains to search: {len(jobs)} (cached {len(cache)})", flush=True)

for key, r, c, seq in jobs:
    if key in cache:
        continue
    for a in range(4):
        try:
            total, hits = homolog_search(seq)
            cache[key] = {"total": total, "hits": hits[:8], "len": len(seq),
                          "desc": c["desc"], "pdb": r["id"], "entity": c["entity"],
                          "no": r["no"], "category": r["category"]}
            top = hits[0]["identity"] * 100 if hits else 0
            print(f"  #{r['no']:2d} {key:12s} len={len(seq):4d} "
                  f"pre-cutoff hits={total:5d} best={top:5.1f}%  {c['desc'][:40]}",
                  flush=True)
            break
        except Exception as ex:
            print(f"  retry {key} {str(ex)[:80]}", flush=True)
            time.sleep(5)
    else:
        cache[key] = None
    json.dump(cache, open("homolog.json", "w"))
    time.sleep(0.4)

print("done:", len(cache))
