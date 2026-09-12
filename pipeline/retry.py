"""Retry the searches that failed, with longer timeouts and more patience.

A failed search was being silently dropped by the rollup, which then reported a
structure's homology from whatever chain DID succeed -- so 14-3-3 sigma (278
pre-2021 homologs at >=95%) was reported as having none, because only its 13-aa
partner peptide's result survived. Failures must be resolved, not skipped.
"""
import json, time, urllib.request

OP = urllib.request.build_opener(urllib.request.ProxyHandler({}))
S = "https://search.rcsb.org/rcsbsearch/v2/query"
BINS = [0.95, 0.60, 0.30]


def post(p, t=240):
    r = urllib.request.Request(S, data=json.dumps(p).encode(),
                               headers={"Content-Type": "application/json"})
    resp = OP.open(r, timeout=t)
    b = resp.read()
    return {"total_count": 0} if not b.strip() else json.loads(b)


def count(seq, ident):
    q = {"query": {"type": "group", "logical_operator": "and", "nodes": [
            {"type": "terminal", "service": "sequence",
             "parameters": {"evalue_cutoff": 1, "identity_cutoff": ident,
                            "sequence_type": "protein", "value": seq}},
            {"type": "terminal", "service": "text",
             "parameters": {"attribute": "rcsb_accession_info.initial_release_date",
                            "operator": "less_or_equal", "value": "2021-09-30"}}]},
         "return_type": "polymer_entity",
         "request_options": {"paginate": {"start": 0, "rows": 1},
                             "results_content_type": ["experimental"]}}
    last = None
    for a in range(6):
        try:
            return post(q).get("total_count", 0)
        except Exception as e:
            last = e
            time.sleep(4 + 4 * a)
    print("   STILL FAILING:", type(last).__name__, str(last)[:70], flush=True)
    return None


sc = json.load(open("screen.json"))
todo = [s for s, v in sc.items() if v.get("bin") is None]
print("retrying", len(todo), flush=True)
for i, seq in enumerate(todo, 1):
    got = None
    for b in BINS:
        n = count(seq, b)
        if n is None:
            got = None
            break
        time.sleep(0.4)
        if n > 0:
            got = (b, n)
            break
    else:
        got = (0.0, 0)
    if got:
        sc[seq] = {"bin": got[0], "n": got[1], "len": len(seq)}
        print(f"  {i}/{len(todo)} len={len(seq):4d} -> bin={got[0]} n={got[1]}", flush=True)
    else:
        sc[seq] = {"bin": None, "n": None, "len": len(seq), "failed": True}
        print(f"  {i}/{len(todo)} len={len(seq):4d} -> unresolved", flush=True)
    json.dump(sc, open("screen.json", "w"))
print("done", flush=True)
