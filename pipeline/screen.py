"""Homolog screening over candidate pools, to see whether cleaner (lower-homology)
alternatives exist before re-selecting.

Scope agreed: PTM (all), membrane (all), hetero (100), ligand (100).

Two cost reductions, both justified by the 41-target data already collected:

  * TRAINING WINDOW ONLY (released <= 2021-09-30). For the 41 targets the
    training and template windows returned the same identity bin in almost every
    case, so searching both would double the runtime for almost no extra signal.
    The template window is re-checked later, only on the structures actually
    selected.

  * DEDUPLICATE BY SEQUENCE. The pools contain the same protein in several
    entries; identical sequences get one search.

Empty HTTP body = zero hits (RCSB returns 204), which is real information.
"""
import json
import os
import re
import time
import urllib.request

OP = urllib.request.build_opener(urllib.request.ProxyHandler({}))
S = "https://search.rcsb.org/rcsbsearch/v2/query"
TRAIN = "2021-09-30"
BINS = [0.95, 0.60, 0.30]          # coarse ladder: enough to rank candidates

BUILTIN_LIG = {"ADP", "ATP", "AMP", "GTP", "GDP", "FAD", "NAD", "NAP", "NDP",
               "HEM", "HEC", "PLM", "OLA", "MYR", "CIT", "CLA", "CHL", "BCL", "BCB"}
MEMB_ANNOT = {"PDBTM", "MemProtMD", "mpstruc"}
AB_RE = re.compile(r"heavy chain|light chain|\bfab\b|\bigg\b|antibody|immunoglobulin"
                   r"|\bscfv\b|\bnanobody\b|\bvhh\b|single[- ]domain antibody|sybody"
                   r"|\bmab\b|\bfv fragment\b|\bvariable domain\b", re.I)
NB_RE = re.compile(r"nanobody|\bvhh\b|single[- ]domain antibody|sybody|camelid", re.I)
TCR_RE = re.compile(r"\bT[- ]cell receptor\b|\bTCR\b|beta-2 microglobulin|\bMHC\b"
                    r"|\bHLA\b|histocompatibility|\bH-2[A-Z]", re.I)


def post(p, t=180):
    r = urllib.request.Request(S, data=json.dumps(p).encode(),
                               headers={"Content-Type": "application/json"})
    resp = OP.open(r, timeout=t)
    b = resp.read()
    return {"total_count": 0} if not b.strip() else json.loads(b)


def count(seq, ident, upto=TRAIN):
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
            time.sleep(3 + 3 * a)
    return None


def probe(seq):
    for b in BINS:
        n = count(seq, b)
        time.sleep(0.2)
        if n is None:
            return None, None
        if n > 0:
            return b, n
    return 0.0, 0


# ---------- rebuild the category pools ----------
valid = json.load(open("census_valid.json"))


def P(r):
    return [c for c in r["chains"] if c["type"] == "Protein"]


def NUC(r):
    return [c for c in r["chains"] if c["type"] in ("DNA", "RNA")]


def blob(r):
    return r["title"] + " || " + " | ".join(c["desc"] for c in r["chains"])


def memb(r):
    return any(set(c["annot"]) & MEMB_ANNOT for c in r["chains"])


def has_ab_chain(r):
    if AB_RE.search(r["title"]) or NB_RE.search(r["title"]):
        return True
    return any(AB_RE.search(c["desc"]) or NB_RE.search(c["desc"]) for c in r["chains"])


def cofac(r):
    return [l for l in r["ligands"] if l["code"] in BUILTIN_LIG]


def otherlig(r):
    return [l for l in r["ligands"] if l["code"] not in BUILTIN_LIG
            and (l["mw"] or 0) >= 150 and (l["atoms"] or 0) >= 10]


def clean(r):
    return not r["branched"] and not r["sugars"]


def resol(r):
    return r["resolution"] or 9.9


def uniprots(r):
    s = set()
    for c in r["chains"]:
        s |= set(c["uniprots"])
    return s


POOLS = {}
POOLS["ptm"] = [r for r in valid.values()
                if r["ptms"] and len(P(r)) >= 2 and not has_ab_chain(r) and clean(r)]
POOLS["memb"] = [r for r in valid.values()
                 if memb(r) and P(r) and not has_ab_chain(r) and clean(r)
                 and 250 <= sum(c["len"] * c["copies"] for c in r["chains"]) <= 2400]
POOLS["hetero"] = [r for r in valid.values()
                   if len(P(r)) >= 2 and not NUC(r) and not has_ab_chain(r)
                   and not memb(r) and not TCR_RE.search(blob(r))
                   and not cofac(r) and not otherlig(r)
                   and all(c["len"] >= 50 for c in P(r)) and len(uniprots(r)) >= 2
                   and clean(r) and resol(r) <= 2.6]
POOLS["ligand"] = [r for r in valid.values()
                   if otherlig(r) and not cofac(r) and not NUC(r)
                   and not has_ab_chain(r) and not memb(r)
                   and clean(r) and resol(r) <= 2.2]

LIMITS = {"ptm": None, "memb": None, "hetero": 100, "ligand": 100}
for k in POOLS:
    POOLS[k].sort(key=resol)
    if LIMITS[k]:
        POOLS[k] = POOLS[k][:LIMITS[k]]
    print(f"{k:8s} screening {len(POOLS[k])} structures", flush=True)

# unique sequences to search
seqmap = {}
for cat, rs in POOLS.items():
    for r in rs:
        for c in P(r):
            if len(c["seq"]) < 12:
                continue
            seqmap.setdefault(c["seq"], []).append((cat, r["id"], c["entity"]))
print(f"unique protein sequences to search: {len(seqmap)}", flush=True)

cache = json.load(open("screen.json")) if os.path.exists("screen.json") else {}
todo = [s for s in seqmap if s not in cache]
print(f"cached {len(cache)}, to search {len(todo)}", flush=True)

t0 = time.time()
for i, seq in enumerate(todo, 1):
    b, n = probe(seq)
    cache[seq] = {"bin": b, "n": n, "len": len(seq)}
    if i % 25 == 0:
        json.dump(cache, open("screen.json", "w"))
        el = time.time() - t0
        print(f"  {i}/{len(todo)}  elapsed {el/60:.1f}min  "
              f"eta {(el/i)*(len(todo)-i)/60:.1f}min", flush=True)
json.dump(cache, open("screen.json", "w"))
json.dump({s: v for s, v in seqmap.items()}, open("screen_map.json", "w"))
print("done", len(cache), flush=True)
