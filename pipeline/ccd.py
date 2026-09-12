"""Stage 3: for every ligand seen in the pools, fetch its CCD initial release
date. A code is usable on AlphaFold Server only if it already existed in CCD
version 2024_10_28."""
import json
import time
import urllib.request

OPENER = urllib.request.build_opener(urllib.request.ProxyHandler({}))
GRAPHQL = "https://data.rcsb.org/graphql"

Q = """
query($ids:[String!]!){
  chem_comps(comp_ids:$ids){
    rcsb_id
    chem_comp { id name formula type formula_weight }
    rcsb_chem_comp_info {
      initial_release_date
      atom_count
      atom_count_heavy
ive: bond_count
    }
  }
}
"""

# the 'ive:' alias above is a typo guard -- rebuild cleanly:
Q = """
query($ids:[String!]!){
  chem_comps(comp_ids:$ids){
    rcsb_id
    chem_comp { id name formula type formula_weight }
    rcsb_chem_comp_info {
      initial_release_date
      atom_count
      atom_count_heavy
      bond_count
    }
  }
}
"""


def post(payload, timeout=180):
    req = urllib.request.Request(
        GRAPHQL, data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json"})
    return json.load(OPENER.open(req, timeout=timeout))


meta = json.load(open("meta.json"))
codes = set()
for e in meta.values():
    for ne in e.get("nonpolymer_entities") or []:
        cc = (ne.get("nonpolymer_comp") or {}).get("chem_comp") or {}
        if cc.get("id"):
            codes.add(cc["id"])

codes = sorted(codes)
print("unique ligand codes:", len(codes))

out = {}
B = 50
for i in range(0, len(codes), B):
    chunk = codes[i:i + B]
    for attempt in range(3):
        try:
            d = post({"query": Q, "variables": {"ids": chunk}})
            for c in d["data"]["chem_comps"] or []:
                out[c["rcsb_id"]] = c
            break
        except Exception as ex:
            print("retry", i, ex)
            time.sleep(3)

json.dump(out, open("ccd.json", "w"))
print("fetched", len(out))

CUT = "2024-10-28"
late = [k for k, v in out.items()
        if (v.get("rcsb_chem_comp_info") or {}).get("initial_release_date", "")[:10] > CUT]
print(f"ligands NOT in CCD {CUT} (unusable on Server): {len(late)}")
print(sorted(late)[:60])
