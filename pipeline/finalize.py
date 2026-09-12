"""Finalize the 41: correct copy numbers to the biological assembly, drop
non-biological additives, and emit the assignment pack.

Two defects fixed here:

1. Copy numbers came from the asymmetric unit. classify.py counted
   auth_asym_ids, which is the ASU, so 10VC appeared as 6 copies of each chain
   (crystal packing) rather than the 1:1 complex a student should submit.
   Copies are now taken from biological assembly 1, scaled up when the instance
   list is truncated by symmetry operators.

2. Detergents and buffer salts were being emitted as ligands. LMT (lauryl
   maltoside) is the detergent used to solubilise a membrane protein and BO3 is
   a borate buffer ion -- neither belongs in the model. Cholesterol (CLR) is
   kept, since for a membrane protein it is a genuine structural lipid.
"""
import json
import os
import time
import urllib.request
from collections import Counter

OP = urllib.request.build_opener(urllib.request.ProxyHandler({}))
G = "https://data.rcsb.org/graphql"
OUT = "/Users/lilindu/alphafold-homework-2026"

# additives that must never be entered: detergents, cryoprotectants, buffers
EXTRA_IGNORE = {
    "LMT", "LDA", "LMN", "DDQ", "C8E", "BOG", "OGA", "BNG", "MC3", "DMU", "UMQ",
    "F09", "LMU", "PCW", "P6G", "TWT", "SDS", "TRT", "TX4", "BO3", "B3P", "BEZ",
    "NHE", "CXS", "MPO", "EPE", "TAM", "BTB", "MRD", "SGM", "PGE", "1PE", "2PE",
    "7PE", "12P", "15P", "XPE", "PE4", "PE8", "PEG", "PG4", "PG0", "GOL", "EDO",
    "MPD", "DMS", "IPA", "MOH", "EOH", "ACT", "FMT", "SO4", "PO4", "NO3", "AZI",
    "SCN", "IMD", "TRS", "MES", "CAC", "OCT", "D10", "DD9", "UND", "HEZ", "PGO",
    "PDO", "BME", "DTT", "DTU", "DTV", "TCE", "BU3", "NH4", "UNX", "UNL", "PIN",
}
BUILTIN_LIG = {"ADP", "ATP", "AMP", "GTP", "GDP", "FAD", "NAD", "NAP", "NDP",
               "HEM", "HEC", "PLM", "OLA", "MYR", "CIT", "CLA", "CHL", "BCL", "BCB"}
BUILTIN_ION = {"MG", "ZN", "CL", "CA", "NA", "MN", "K", "FE", "CU", "CO"}
ENTKEY = {"Protein": "proteinChain", "DNA": "dnaSequence", "RNA": "rnaSequence"}

QA = """
query($ids:[String!]!){
 entries(entry_ids:$ids){
  rcsb_id
  assemblies{
   rcsb_id
   rcsb_assembly_info{polymer_entity_instance_count}
   polymer_entity_instances{
    rcsb_polymer_entity_instance_container_identifiers{entity_id}}}}}
"""


def post(p, t=300):
    r = urllib.request.Request(G, data=json.dumps(p).encode(),
                               headers={"Content-Type": "application/json"})
    return json.load(OP.open(r, timeout=t))


alloc = json.load(open("alloc41.json"))
ids = sorted({it["rec"]["id"] for v in alloc.values() for it in v})
print("targets:", len(ids))

raw = {}
for i in range(0, len(ids), 20):
    ch = ids[i:i + 20]
    for a in range(4):
        try:
            d = post({"query": QA, "variables": {"ids": ch}})
            if "data" not in d:
                raise RuntimeError(str(d.get("errors"))[:150])
            for e in d["data"]["entries"] or []:
                raw[e["rcsb_id"]] = e.get("assemblies") or []
            break
        except Exception as ex:
            print("retry", str(ex)[:90], flush=True)
            time.sleep(4)

stoich = {}
for eid, asms in raw.items():
    if not asms:
        stoich[eid] = None
        continue
    a = asms[0]
    rep = (a.get("rcsb_assembly_info") or {}).get("polymer_entity_instance_count") or 0
    c = Counter(
        x["rcsb_polymer_entity_instance_container_identifiers"]["entity_id"]
        for x in (a.get("polymer_entity_instances") or [])
        if x.get("rcsb_polymer_entity_instance_container_identifiers"))
    listed = sum(c.values())
    per = dict(c)
    note = ""
    if listed and rep and listed != rep:
        if rep % listed == 0:
            f = rep // listed
            per = {k: v * f for k, v in per.items()}
            note = f"按对称操作展开 ×{f}"
        else:
            note = f"PDB 报告 {rep} 条链,实例表列出 {listed} 条,请自行核对拷贝数"
    stoich[eid] = {"per": per, "reported": rep, "listed": listed, "note": note}
json.dump(stoich, open("stoich41.json", "w"), indent=1)

changed = []
records = []
for cls in ["班1", "班2", "班3", "班4"]:
    for it in alloc[cls]:
        r = dict(it["rec"])
        eid = r["id"]
        st = stoich.get(eid)

        chains = []
        for c in r["chains"]:
            c = dict(c)
            if st and st["per"]:
                newn = st["per"].get(c["entity"])
                if newn and newn != c["copies"]:
                    changed.append((eid, c["entity"], c["copies"], newn))
                    c["copies"] = newn
            chains.append(c)

        ligs, ions, dropped = [], [], []
        for l in r["ligands"]:
            if l["code"] in EXTRA_IGNORE:
                dropped.append((l["code"], l["count"]))
            else:
                ligs.append(dict(l))
        for code, n in r["ions"]:
            if code in BUILTIN_ION:
                ions.append({"code": code, "count": n})
            else:
                dropped.append((code, n))

        tokens = sum(c["len"] * c["copies"] for c in chains)
        tokens += sum((l["atoms"] or 0) * l["count"] for l in ligs)
        tokens += sum(i["count"] for i in ions)

        records.append({
            "class": cls, "category": it["category"], "label": it["label"],
            "note": it["note"], "id": eid, "title": r["title"],
            "resolution": r["resolution"], "method": r["method"],
            "deposit": r["deposit"], "release": r["release"],
            "chains": chains, "ligands": ligs, "ions": ions,
            "omitted": list(r["omitted"]) + dropped,
            "ptms": r["ptms"], "tokens": tokens,
            "assembly_note": (st or {}).get("note", ""),
        })

print(f"\ncopy-number corrections (ASU -> biological assembly): {len(changed)}")
for eid, ent, old, new in changed:
    print(f"   {eid} entity {ent}: {old} -> {new}")

over = [r for r in records if r["tokens"] > 5000]
print(f"\ntoken range: {min(r['tokens'] for r in records)}-{max(r['tokens'] for r in records)}"
      f"  over 5000: {len(over)}")

for c in ["班1", "班2", "班3", "班4"]:
    n = sum(1 for r in records if r["class"] == c)
    print(f"  {c}: {n} 题")
print("TOTAL:", len(records))

os.makedirs(OUT, exist_ok=True)
json.dump(records, open(os.path.join(OUT, "targets41.json"), "w"),
          ensure_ascii=False, indent=1)
print("wrote targets41.json")
