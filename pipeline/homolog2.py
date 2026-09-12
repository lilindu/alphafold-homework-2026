"""Was a homolog of each target already available to AlphaFold?

Two distinct leakage pathways, with two different date bounds:

  TRAINING window  -- released <= 2021-09-30 (the AF3 training cutoff)
      "has the model memorised this fold from training data?"

  TEMPLATE window  -- released <= 2025-02-03 (the latest template date a student
      can select on the Server; the default is 2021-09-30)
      "could a homolog be handed to the model as a structural template at
       inference time?"

The template window is the stricter test and strictly contains the training
window. Structures released between the two dates were never trained on, but
their coordinates can still be fed in directly as templates -- a more direct form
of leakage than memorisation.

Both are searched server-side per chain (RCSB sequence search = MMseqs2) so each
answer is the true best hit in that window, not a sample of a larger result set.
"""
import json
import os
import time
import urllib.request

OP = urllib.request.build_opener(urllib.request.ProxyHandler({}))
SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"

TRAIN_CUTOFF = "2021-09-30"     # AF3 training cutoff (= default template cutoff)
TEMPLATE_MAX = "2025-02-03"     # latest template date selectable on the Server


def post(payload, timeout=300):
    req = urllib.request.Request(
        SEARCH, data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json"})
    return json.load(OP.open(req, timeout=timeout))


def best_hit(seq, upto, rows=20):
    q = {
        "query": {
            "type": "group", "logical_operator": "and",
            "nodes": [
                {"type": "terminal", "service": "sequence",
                 "parameters": {"evalue_cutoff": 1, "identity_cutoff": 0.2,
                                "sequence_type": "protein", "value": seq}},
                {"type": "terminal", "service": "text",
                 "parameters": {
                     "attribute": "rcsb_accession_info.initial_release_date",
                     "operator": "less_or_equal", "value": upto}},
            ]},
        "return_type": "polymer_entity",
        "request_options": {"paginate": {"start": 0, "rows": rows},
                            "results_content_type": ["experimental"],
                            "scoring_strategy": "sequence"},
    }
    d = post(q)
    best = None
    for x in d.get("result_set", []) or []:
        for s in x.get("services", []) or []:
            for n in s.get("nodes", []) or []:
                for mm in n.get("match_context", []) or []:
                    si = mm.get("sequence_identity")
                    if si is None:
                        continue
                    cand = {"entity": x["identifier"], "identity": si,
                            "evalue": mm.get("evalue"),
                            "alen": mm.get("alignment_length")}
                    if best is None or si > best["identity"]:
                        best = cand
    return d.get("total_count", 0), best


recs = json.load(open("targets41_full.json"))
cache = json.load(open("homolog2.json")) if os.path.exists("homolog2.json") else {}

targets = []
for r in recs:
    for c in r["chains"]:
        if c["type"] != "Protein":
            continue
        seq = c.get("construct_seq") or c["seq"]     # what is actually resolved
        if len(seq) < 12:
            continue
        targets.append((f"{r['id']}_{c['entity']}", r, c, seq))

print(f"chains to search: {len(targets)}  cached: {len(cache)}", flush=True)

for key, r, c, seq in targets:
    if key in cache:
        continue
    rec = {"no": r["no"], "pdb": r["id"], "entity": c["entity"],
           "category": r["category"], "desc": c["desc"], "len": len(seq)}
    ok = True
    for tag, date in (("train", TRAIN_CUTOFF), ("template", TEMPLATE_MAX)):
        for a in range(4):
            try:
                tot, b = best_hit(seq, date)
                rec[tag] = {"total": tot, "best": b}
                break
            except Exception as ex:
                if a == 3:
                    ok = False
                time.sleep(4)
        time.sleep(0.3)
    if not ok:
        print(f"  FAILED {key}", flush=True)
    cache[key] = rec
    tr = rec.get("train", {}).get("best")
    tp = rec.get("template", {}).get("best")
    print(f"  #{r['no']:2d} {key:12s} len={len(seq):4d}  "
          f"训练窗口 {(tr['identity']*100 if tr else 0):5.1f}%  "
          f"模板窗口 {(tp['identity']*100 if tp else 0):5.1f}%   {c['desc'][:34]}",
          flush=True)
    json.dump(cache, open("homolog2.json", "w"))

print("done:", len(cache))
