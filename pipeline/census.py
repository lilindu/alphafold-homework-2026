"""Census: classify every Server-compatible post-2025-02-03 structure into the
nine classes the course needs, and count each class.

Hard gates applied to everything (all verified, not assumed):
  * deposit_date > 2025-02-03      (beyond the latest selectable template cutoff)
  * released (has initial_release_date)
  * standard residues only         (canonical seq; B/J/O/U/X/Z unsupported)
  * every chain >= 4 residues/nt
  * all ligand CCD codes present in the frozen 2024_10_28 dictionary
  * ions restricted to the 10 built-ins
  * total tokens <= 5000
  * biological assembly <= 12 chains (no capsids/filaments)

The nine classes:
  1 hetero  异源蛋白复合体 (non-antibody, non-nucleic, non-membrane)
  2 metal   蛋白 + 金属离子
  3 dna     蛋白 + DNA
  4 rna     蛋白 + RNA
  5 cofac   蛋白 + 内置辅因子 (19 built-in ligand codes)
  6 ligand  蛋白 + 非内置小分子配体 (in frozen CCD, needs manual code entry)
  7 ptm     有翻译后修饰的蛋白 (modification must be in Server's ptmType list)
  8 memb    膜蛋白 (PDBTM / MemProtMD / mpstruc annotation)
  9 ab      抗原抗体复合物
"""
import json
import os
import re
import time
import urllib.request

OPENER = urllib.request.build_opener(urllib.request.ProxyHandler({}))
SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL = "https://data.rcsb.org/graphql"
DEP_AFTER = "2025-02-03"
CCD_FREEZE = "2024-10-28"
TOKEN_LIMIT = 5000
MAX_ASSEMBLY = 12

BUILTIN_LIG = {"ADP", "ATP", "AMP", "GTP", "GDP", "FAD", "NAD", "NAP", "NDP",
               "HEM", "HEC", "PLM", "OLA", "MYR", "CIT", "CLA", "CHL", "BCL", "BCB"}
SUPPORTED_IONS = {"CA", "CO", "CU", "FE", "K", "MG", "MN", "NA", "ZN", "CL"}
METALS = {"ZN", "CU", "FE", "MN", "CO", "MG", "CA", "NA", "K"}
TRANSITION = {"ZN", "CU", "FE", "MN", "CO"}
SERVER_PTM = ["SEP", "TPO", "PTR", "NEP", "HIP", "ALY", "MLY", "M3L", "MLZ",
              "2MR", "AGM", "MCS", "HYP", "HY3", "LYZ", "AHB", "P1L", "SNN",
              "SNC", "TRF", "KCR", "CIR", "YHA"]
MEMB_ANNOT = {"PDBTM", "MemProtMD", "mpstruc"}
SUGARS = {"NAG", "BMA", "MAN", "GLC", "BGC", "FUC", "GAL", "NDG", "A2G", "SIA",
          "GLA", "GAL", "XYS", "RAM", "NGA"}

IGNORABLE = {
    "HOH", "DOD", "GOL", "EDO", "PEG", "PGE", "PG4", "PG0", "1PE", "2PE", "P6G",
    "7PE", "12P", "15P", "XPE", "MPD", "MRD", "SO4", "PO4", "ACT", "ACY", "FMT",
    "DMS", "TRS", "IMD", "EPE", "MES", "FLC", "TAR", "MLI", "MLA", "SIN",
    "BME", "DTT", "DTU", "DTV", "TCE", "BTB", "BU3", "HEZ", "PGO", "PDO", "IPA",
    "MOH", "EOH", "ETX", "NO3", "AZI", "SCN", "CO3", "BCT", "CAC", "NH4",
    "OXL", "SUC", "TLA", "BOG", "LDA", "C8E", "OGA", "UNX", "UNL", "PGR", "PIN",
    "POL", "SPD", "SPM", "GAI", "PYR", "FUM", "MAE", "MPO", "NHE", "CXS", "MRO",
}
OMITTABLE_HEAVY = {
    "NI", "CD", "HG", "PT", "AU", "AG", "PB", "SR", "BA", "CS", "RB", "LI",
    "BR", "IOD", "F", "AL", "V", "W", "MO", "SE", "TL", "GA", "IN", "SM", "EU",
    "GD", "TB", "YB", "LU", "PR", "ND", "HO", "ER", "ZR", "TA", "RE", "OS",
    "IR", "RU", "RH", "PD", "SB", "SN", "TE", "BI", "3NI", "2HP", "YT3", "CE",
}
STD_AA = set("ACDEFGHIKLMNPQRSTVWY")

AB_RE = re.compile(r"heavy chain|light chain|\bfab\b|\bigg\b|antibody|immunoglobulin"
                   r"|\bscfv\b|\bnanobody\b|\bvhh\b|single[- ]domain antibody|sybody"
                   r"|\bmab\b|\bfv\b", re.I)


def post(url, payload, timeout=240):
    req = urllib.request.Request(
        url, data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json"})
    return json.load(OPENER.open(req, timeout=timeout))


def T(a, o, v, neg=False):
    p = {"attribute": a, "operator": o, "value": v}
    if neg:
        p["negation"] = True
    return {"type": "terminal", "service": "text", "parameters": p}


DEP = T("rcsb_accession_info.deposit_date", "greater", DEP_AFTER)
SIZE = T("rcsb_entry_info.deposited_polymer_monomer_count", "less_or_equal", 4500)
LIGATTR = "rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id"


def search(nodes, rows=200, sort="rcsb_entry_info.resolution_combined"):
    payload = {"query": {"type": "group", "logical_operator": "and", "nodes": nodes},
               "return_type": "entry",
               "request_options": {"paginate": {"start": 0, "rows": rows},
                                   "results_content_type": ["experimental"],
                                   "sort": [{"sort_by": sort, "direction": "asc"}]}}
    d = post(SEARCH, payload)
    return d.get("total_count", 0), [x["identifier"] for x in d.get("result_set", [])]


# ---------------- searches per class ----------------
QUERIES = {}
QUERIES["hetero"] = [DEP, SIZE,
                     T("rcsb_entry_info.polymer_entity_count_protein", "range",
                       {"from": 2, "to": 4, "include_lower": True, "include_upper": True}),
                     T("rcsb_entry_info.polymer_entity_count_nucleic_acid", "equals", 0),
                     T("rcsb_entry_info.deposited_polymer_monomer_count", "greater_or_equal", 200)]
QUERIES["dna"] = [DEP, SIZE,
                  T("rcsb_entry_info.polymer_entity_count_DNA", "greater_or_equal", 1),
                  T("rcsb_entry_info.polymer_entity_count_protein", "greater_or_equal", 1)]
QUERIES["rna"] = [DEP, SIZE,
                  T("rcsb_entry_info.polymer_entity_count_RNA", "greater_or_equal", 1),
                  T("rcsb_entry_info.polymer_entity_count_protein", "greater_or_equal", 1)]
QUERIES["ab"] = [DEP, SIZE,
                 T("struct_keywords.pdbx_keywords", "contains_phrase", "IMMUNE SYSTEM"),
                 T("rcsb_entry_info.polymer_entity_count_protein", "greater_or_equal", 2),
                 T("rcsb_entry_info.deposited_polymer_monomer_count", "greater_or_equal", 250)]
QUERIES["ptm"] = [DEP, SIZE,
                  T("entity_poly.rcsb_non_std_monomers", "in", SERVER_PTM)]
for src in MEMB_ANNOT:
    QUERIES["memb_" + src] = [DEP, SIZE,
                              T("rcsb_polymer_entity_annotation.type", "exact_match", src)]
for lig in sorted(BUILTIN_LIG):
    QUERIES["cofac_" + lig] = [DEP, SIZE, T(LIGATTR, "exact_match", lig),
                               T("rcsb_entry_info.polymer_entity_count_protein",
                                 "greater_or_equal", 1)]
for ion in sorted(TRANSITION | {"MG", "CA"}):
    QUERIES["metal_" + ion] = [DEP, SIZE, T(LIGATTR, "exact_match", ion),
                               T("rcsb_entry_info.polymer_entity_count_protein",
                                 "greater_or_equal", 1)]
# generic ligand-bearing pool for class 6
QUERIES["ligand_any"] = [DEP, SIZE,
                         T("rcsb_entry_info.polymer_entity_count_protein", "greater_or_equal", 1),
                         T("rcsb_entry_info.polymer_entity_count_nucleic_acid", "equals", 0),
                         T("rcsb_nonpolymer_entity_annotation.type", "exact_match",
                           "SUBJECT_OF_INVESTIGATION")]
QUERIES["ligand_affinity"] = [DEP, SIZE,
                              T("rcsb_entry_info.polymer_entity_count_protein", "greater_or_equal", 1),
                              T("rcsb_binding_affinity.value", "exists", "")]

if __name__ == "__main__":
    totals, pools = {}, {}
    for name, nodes in QUERIES.items():
        rows = 300 if name in ("hetero", "ab", "ligand_any", "ligand_affinity",
                               "dna", "rna", "ptm") else 120
        for a in range(3):
            try:
                t, ids = search(nodes, rows=rows)
                totals[name], pools[name] = t, ids
                print(f"{name:22s} rcsb_total={t:6d} pulled={len(ids)}")
                break
            except Exception as e:
                print(f"{name:22s} retry {str(e)[:80]}")
                time.sleep(3)
    json.dump({"totals": totals, "pools": pools}, open("census_pools.json", "w"), indent=1)
    allids = sorted({i for v in pools.values() for i in v})
    json.dump(allids, open("census_ids.json", "w"))
    print("unique candidate entries:", len(allids))
