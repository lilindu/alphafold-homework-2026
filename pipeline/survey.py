"""Survey RCSB for AlphaFold-Server homework candidates.

Hard constraints (from AlphaFold Server FAQ + user requirements):
  - deposit_date > 2025-02-03  (beyond the latest selectable template cutoff,
    so no template pathway can leak the answer)
  - total polymer residues <= 4500 (5000-token job limit, leaving ligand room)
  - any ligand's CCD code must exist in CCD version 2024_10_28
    -> chem_comp initial_release_date <= 2024-10-28
  - ions/PTMs limited to the built-in list (checked separately)
"""
import json
import urllib.request

OPENER = urllib.request.build_opener(urllib.request.ProxyHandler({}))
SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"

DEPOSIT_AFTER = "2025-02-03"


def post(url, payload, timeout=180):
    req = urllib.request.Request(
        url, data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json"})
    return json.load(OPENER.open(req, timeout=timeout))


def T(attribute, operator, value, negation=False):
    p = {"attribute": attribute, "operator": operator, "value": value}
    if negation:
        p["negation"] = True
    return {"type": "terminal", "service": "text", "parameters": p}


def FT(value):
    return {"type": "terminal", "service": "full_text",
            "parameters": {"value": value}}


def search(nodes, rows=40, sort_attr="rcsb_entry_info.resolution_combined",
           direction="asc"):
    payload = {
        "query": {"type": "group", "logical_operator": "and", "nodes": nodes},
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": rows},
            "results_content_type": ["experimental"],
            "sort": [{"sort_by": sort_attr, "direction": direction}],
        },
    }
    d = post(SEARCH, payload)
    return d.get("total_count", 0), [x["identifier"] for x in d.get("result_set", [])]


DEPOSITED = T("rcsb_accession_info.deposit_date", "greater", DEPOSIT_AFTER)
SMALL = T("rcsb_entry_info.deposited_polymer_monomer_count", "less_or_equal", 4500)


def size_range(lo, hi):
    return [T("rcsb_entry_info.deposited_polymer_monomer_count", "greater_or_equal", lo),
            T("rcsb_entry_info.deposited_polymer_monomer_count", "less_or_equal", hi)]


def xray_good(max_res=2.5):
    return [T("rcsb_entry_info.resolution_combined", "less_or_equal", max_res)]


CATEGORIES = {}

# --- 1. single-chain protein, one protein entity, no nucleic acid ---
CATEGORIES["mono_protein"] = [
    DEPOSITED,
    T("rcsb_entry_info.polymer_entity_count_protein", "equals", 1),
    T("rcsb_entry_info.polymer_entity_count_nucleic_acid", "equals", 0),
    T("rcsb_entry_info.deposited_polymer_entity_instance_count", "equals", 1),
    *size_range(120, 500),
    *xray_good(2.0),
]

# --- 2. homodimer: 1 protein entity, 2 instances ---
CATEGORIES["homodimer"] = [
    DEPOSITED,
    T("rcsb_entry_info.polymer_entity_count_protein", "equals", 1),
    T("rcsb_entry_info.polymer_entity_count_nucleic_acid", "equals", 0),
    T("rcsb_entry_info.deposited_polymer_entity_instance_count", "equals", 2),
    *size_range(200, 1000),
    *xray_good(2.3),
]

# --- 3. heteromeric protein complex: >=2 distinct protein entities ---
CATEGORIES["hetero_complex"] = [
    DEPOSITED,
    T("rcsb_entry_info.polymer_entity_count_protein", "equals", 2),
    T("rcsb_entry_info.polymer_entity_count_nucleic_acid", "equals", 0),
    *size_range(200, 1200),
    *xray_good(2.8),
]

# --- 4. protein-DNA ---
CATEGORIES["protein_dna"] = [
    DEPOSITED,
    T("rcsb_entry_info.polymer_entity_count_DNA", "greater_or_equal", 1),
    T("rcsb_entry_info.polymer_entity_count_protein", "greater_or_equal", 1),
    *size_range(100, 1500),
]

# --- 5. protein-RNA ---
CATEGORIES["protein_rna"] = [
    DEPOSITED,
    T("rcsb_entry_info.polymer_entity_count_RNA", "greater_or_equal", 1),
    T("rcsb_entry_info.polymer_entity_count_protein", "greater_or_equal", 1),
    *size_range(100, 1500),
]

# --- 6. antibody / Fab / nanobody complexes (immune system keyword) ---
CATEGORIES["antibody"] = [
    DEPOSITED,
    T("struct_keywords.pdbx_keywords", "contains_phrase", "IMMUNE SYSTEM"),
    T("rcsb_entry_info.polymer_entity_count_protein", "greater_or_equal", 3),
    *size_range(400, 1500),
    *xray_good(3.0),
]

# --- 7. membrane proteins via annotation ---
CATEGORIES["membrane"] = [
    DEPOSITED,
    T("rcsb_polymer_entity_annotation.type", "exact_match", "PDBTM"),
    *size_range(200, 2000),
]

for name, nodes in CATEGORIES.items():
    try:
        total, ids = search(nodes, rows=25)
        print(f"{name}: total={total}")
        print("   ", " ".join(ids))
    except Exception as e:
        print(f"{name}: ERROR {e}")
