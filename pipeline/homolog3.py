"""How much did AlphaFold already know about each target?

The v2 search API does not return `match_context`, so per-hit sequence identity
cannot be read out (my previous attempt parsed a field that is never present and
reported 0.0% for every chain -- including myoglobin, which has 1131 pre-2021
homologs; the number was meaningless).

`total_count` IS reliable, so identity is obtained by probing the
`identity_cutoff` parameter: ask "are there any pre-cutoff hits at >=95%?", then
>=80%, and so on. The first threshold that returns hits brackets the best
homolog. This uses only counts, which the API definitely reports.

Two windows, because there are two leakage pathways:
  TRAIN    released <= 2021-09-30  -- could be memorised during training
  TEMPLATE released <= 2025-02-03  -- can be fed in as a structural template at
           inference time (the latest date a student can select). Strictly
           contains the training window, so it is the stricter test.
"""
import json
import os
import time
import urllib.request

OP = urllib.request.build_opener(urllib.request.ProxyHandler({}))
S = "https://search.rcsb.org/rcsbsearch/v2/query"
TRAIN = "2021-09-30"
TEMPLATE = "2025-02-03"
BINS = [0.95, 0.80, 0.60, 0.40, 0.25]


def post(p, t=180):
    """RCSB replies 204 with an EMPTY body when nothing matches, which json.load
    cannot parse. An empty body means zero hits -- real information, not an
    error -- so it is returned as an explicit zero-count result."""
    r = urllib.request.Request(S, data=json.dumps(p).encode(),
                               headers={"Content-Type": "application/json"})
    resp = OP.open(r, timeout=t)
    body = resp.read()
    if not body.strip():
        return {"total_count": 0}
    return json.loads(body)


def count(seq, upto, ident):
    q = {"query": {"type": "group", "logical_operator": "and", "nodes": [
            {"type": "terminal", "service": "sequence",
             "parameters": {"evalue_cutoff": 1, "identity_cutoff": ident,
                            "sequence_type": "protein", "value": seq}},
            {"type": "terminal", "service": "text",
             "parameters": {"attribute": "rcsb_accession_info.initial_release_date",
                            "operator": "less_or_equal", "value": upto}}]},
         "return_type": "polymer_entity",
         "request_options": {"paginate": {"start": 0, "rows": 1},
                             "results_content_type": ["experimental"]}}
    for a in range(4):
        try:
            return post(q).get("total_count", 0)
        except Exception:
            time.sleep(3 + 2 * a)
    return None


def probe(seq, upto):
    """-> (best_identity_bin, count_at_that_bin, count_at_25pct)"""
    for b in BINS:
        n = count(seq, upto, b)
        time.sleep(0.25)
        if n is None:
            return None, None, None
        if n > 0:
            wide = count(seq, upto, 0.25) if b != 0.25 else n
            time.sleep(0.25)
            return b, n, wide
    return 0.0, 0, 0


recs = json.load(open("targets41_full.json"))
cache = json.load(open("homolog3.json")) if os.path.exists("homolog3.json") else {}

jobs = []
for r in recs:
    for c in r["chains"]:
        if c["type"] != "Protein":
            continue
        seq = c.get("construct_seq") or c["seq"]
        if len(seq) < 12:
            continue
        jobs.append((f"{r['id']}_{c['entity']}", r, c, seq))

print(f"chains: {len(jobs)}  cached: {len(cache)}", flush=True)
for key, r, c, seq in jobs:
    v = cache.get(key)
    if v and v.get("train_n") is not None and v.get("tmpl_n") is not None:
        continue
    tb, tn, tw = probe(seq, TRAIN)
    pb, pn, pw = probe(seq, TEMPLATE)
    cache[key] = {"no": r["no"], "pdb": r["id"], "entity": c["entity"],
                  "category": r["category"], "desc": c["desc"], "len": len(seq),
                  "train_bin": tb, "train_n": tn, "train_wide": tw,
                  "tmpl_bin": pb, "tmpl_n": pn, "tmpl_wide": pw}
    json.dump(cache, open("homolog3.json", "w"))
    print(f"  #{r['no']:2d} {key:12s} len={len(seq):4d}  "
          f"train>={(tb or 0)*100:4.0f}% (n={tn}, wide={tw})  "
          f"tmpl>={(pb or 0)*100:4.0f}% (n={pn})   {c['desc'][:30]}", flush=True)
print("done", len(cache))
