"""Finalize the 7-class, 41-target pack into the standalone repo.

Applies, in order, every correction established earlier in this project:

  1. allocation fix -- the 13th antibody goes to 班4 (11 题), not 班1
  2. copy numbers from BIOLOGICAL ASSEMBLY 1, not the asymmetric unit
  3. ghost copies dropped (an assembly member with <50% of its residues resolved
     is lattice occupancy, not a real subunit)
  4. protein input switched to FULL-LENGTH UniProt sequence
  5. PTM positions remapped to full-length numbering and verified against the
     parent amino acid
  6. detergents / cryoprotectants / buffer salts excluded from the ligand list
  7. per-chain observation statistics recorded, so the handout can state how many
     residues exist in the model but not in the experiment
"""
import json
import os
import time
import urllib.request

OP = urllib.request.build_opener(urllib.request.ProxyHandler({}))
G = "https://data.rcsb.org/graphql"
REPO = "/Users/lilindu/alphafold_homework_41"

BUILTIN_LIG = {"ADP", "ATP", "AMP", "GTP", "GDP", "FAD", "NAD", "NAP", "NDP",
               "HEM", "HEC", "PLM", "OLA", "MYR", "CIT", "CLA", "CHL", "BCL", "BCB"}
BUILTIN_ION = {"MG", "ZN", "CL", "CA", "NA", "MN", "K", "FE", "CU", "CO"}
STD = set("ACDEFGHIKLMNPQRSTVWY")
TOKEN_LIMIT = 5000
MIN_OBS = 0.50

EXTRA_IGNORE = {
    "LMT", "LDA", "LMN", "DDQ", "C8E", "BOG", "OGA", "BNG", "MC3", "DMU", "UMQ",
    "F09", "LMU", "PCW", "P6G", "TWT", "SDS", "TRT", "TX4", "BO3", "B3P", "BEZ",
    "NHE", "CXS", "MPO", "EPE", "TAM", "BTB", "MRD", "SGM", "PGE", "1PE", "2PE",
    "7PE", "12P", "15P", "XPE", "PE4", "PE8", "PEG", "PG4", "PG0", "GOL", "EDO",
    "MPD", "DMS", "IPA", "MOH", "EOH", "ACT", "FMT", "SO4", "PO4", "NO3", "AZI",
    "SCN", "IMD", "TRS", "MES", "CAC", "OCT", "D10", "DD9", "UND", "HEZ", "PGO",
    "PDO", "BME", "DTT", "DTU", "DTV", "TCE", "BU3", "NH4", "UNX", "UNL", "PIN",
    "CLR", "Y01", "CHS", "PLC", "POV", "PEE", "PGT", "3PH", "DGA", "LPP",
}
PTM_PARENT = {"SEP": "S", "TPO": "T", "PTR": "Y", "NEP": "H", "HIP": "H",
              "ALY": "K", "MLY": "K", "M3L": "K", "MLZ": "K", "KCR": "K",
              "YHA": "K", "CIR": "R", "2MR": "R", "AGM": "R", "HYP": "P",
              "HY3": "P", "LYZ": "K", "AHB": "N", "MCS": "C", "P1L": "C",
              "SNC": "C", "SNN": "N", "TRF": "W"}

Q = """
query($ids:[String!]!){
 entries(entry_ids:$ids){
  rcsb_id
  assemblies{
   rcsb_id
   rcsb_assembly_info{polymer_entity_instance_count}
   polymer_entity_instances{
    rcsb_polymer_entity_instance_container_identifiers{entity_id}}}
  polymer_entities{
   rcsb_id
   entity_poly{pdbx_seq_one_letter_code_can rcsb_entity_polymer_type}
   rcsb_polymer_entity{pdbx_description}
   rcsb_polymer_entity_align{
    reference_database_name reference_database_accession
    aligned_regions{entity_beg_seq_id ref_beg_seq_id length}}
   uniprots{rcsb_id rcsb_uniprot_protein{sequence}}
   polymer_entity_instances{
    rcsb_polymer_entity_instance_container_identifiers{auth_asym_id}
    rcsb_polymer_instance_feature{
     type feature_positions{beg_seq_id end_seq_id}}}}}}
"""


def post(p, t=300):
    r = urllib.request.Request(G, data=json.dumps(p).encode(),
                               headers={"Content-Type": "application/json"})
    return json.load(OP.open(r, timeout=t))


alloc = json.load(open("alloc7.json"))
BAN = ["班1", "班2", "班3", "班4"]

# ---------- 1. allocation fix: extra antibody to 班4 ----------
flat = []
for b in BAN:
    for it in alloc[b]:
        flat.append((b, it))
by_cat = {}
for b, it in flat:
    by_cat.setdefault(it["category"], []).append(it)

fixed = {b: [] for b in BAN}
for key in ["hetero", "dna", "rna", "ligand", "ptm"]:
    for j, it in enumerate(by_cat[key][:4]):
        fixed[BAN[j]].append(it)
for j, it in enumerate(by_cat["memb"][:8]):
    fixed[BAN[j % 4]].append(it)
abs_ = by_cat["ab"]
for j, it in enumerate(abs_[:12]):
    fixed[BAN[j % 4]].append(it)
if len(abs_) > 12:
    fixed["班4"].append(abs_[12])
for b in BAN:
    print(f"{b}: {len(fixed[b])} 题  " +
          str({k: sum(1 for x in fixed[b] if x['category'] == k)
               for k in ['hetero','dna','rna','ligand','ptm','memb','ab']}))

ids = sorted({it["rec"]["id"] for b in BAN for it in fixed[b]})
print("targets:", len(ids))

# ---------- fetch structural metadata ----------
cache = json.load(open("fin7_raw.json")) if os.path.exists("fin7_raw.json") else {}
todo = [i for i in ids if i not in cache]
for i in range(0, len(todo), 10):
    ch = todo[i:i + 10]
    for a in range(4):
        try:
            d = post({"query": Q, "variables": {"ids": ch}})
            if "data" not in d:
                raise RuntimeError(str(d.get("errors"))[:180])
            for e in d["data"]["entries"] or []:
                cache[e["rcsb_id"]] = e
            break
        except Exception as ex:
            print("  retry", str(ex)[:90], flush=True)
            time.sleep(4)
json.dump(cache, open("fin7_raw.json", "w"))
print("metadata cached:", len(cache))
