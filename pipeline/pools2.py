"""Stage 6: targeted searches to fill gaps — built-in-cofactor entries and
genuine membrane proteins (via PDBTM/MemProtMD/mpstruc annotations)."""
import json
import urllib.request

OPENER = urllib.request.build_opener(urllib.request.ProxyHandler({}))
SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
DEP_AFTER = "2025-02-03"

BUILTIN_LIG = ["ATP", "ADP", "AMP", "GTP", "GDP", "FAD", "NAD", "NAP", "NDP",
               "HEM", "HEC", "PLM", "OLA", "MYR", "CIT", "CLA", "CHL", "BCL", "BCB"]


def post(payload, timeout=180):
    req = urllib.request.Request(
        SEARCH, data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json"})
    return json.load(OPENER.open(req, timeout=timeout))


def T(a, o, v):
    return {"type": "terminal", "service": "text",
            "parameters": {"attribute": a, "operator": o, "value": v}}


DEP = T("rcsb_accession_info.deposit_date", "greater", DEP_AFTER)
LIGATTR = "rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id"


def search(nodes, rows=50, sort_attr="rcsb_entry_info.resolution_combined"):
    payload = {
        "query": {"type": "group", "logical_operator": "and", "nodes": nodes},
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": rows},
            "results_content_type": ["experimental"],
            "sort": [{"sort_by": sort_attr, "direction": "asc"}],
        },
    }
    d = post(payload)
    return d.get("total_count", 0), [x["identifier"] for x in d.get("result_set", [])]


pools = {}

# ---- built-in cofactor holo structures, monomeric-ish ----
for lig in BUILTIN_LIG:
    nodes = [DEP, T(LIGATTR, "exact_match", lig),
             T("rcsb_entry_info.polymer_entity_count_nucleic_acid", "equals", 0),
             T("rcsb_entry_info.deposited_polymer_monomer_count", "greater_or_equal", 150),
             T("rcsb_entry_info.deposited_polymer_monomer_count", "less_or_equal", 1200),
             T("rcsb_entry_info.resolution_combined", "less_or_equal", 2.4)]
    try:
        total, ids = search(nodes, rows=30)
        pools["lig_" + lig] = ids
        print(f"lig_{lig:4s} total={total:5d} kept={len(ids)}")
    except Exception as e:
        print("lig_" + lig, "ERR", str(e)[:100])

# ---- genuine membrane proteins ----
for src in ["PDBTM", "MemProtMD", "mpstruc"]:
    nodes = [DEP, T("rcsb_polymer_entity_annotation.type", "exact_match", src),
             T("rcsb_entry_info.deposited_polymer_monomer_count", "greater_or_equal", 200),
             T("rcsb_entry_info.deposited_polymer_monomer_count", "less_or_equal", 2600)]
    try:
        total, ids = search(nodes, rows=50)
        pools["memb_" + src] = ids
        print(f"memb_{src:10s} total={total:5d} kept={len(ids)}")
    except Exception as e:
        print("memb_" + src, "ERR", str(e)[:100])

# ---- GPCR / channel / transporter by name, multi-subunit signalling ----
for kw, label in [("G protein-coupled receptor", "gpcr"), ("ion channel", "channel"),
                  ("transporter", "transporter"), ("ABC transporter", "abc")]:
    nodes = [DEP, T("rcsb_polymer_entity.pdbx_description", "contains_phrase", kw),
             T("rcsb_entry_info.deposited_polymer_monomer_count", "greater_or_equal", 250),
             T("rcsb_entry_info.deposited_polymer_monomer_count", "less_or_equal", 2600)]
    try:
        total, ids = search(nodes, rows=40)
        pools["kw_" + label] = ids
        print(f"kw_{label:12s} total={total:5d} kept={len(ids)}")
    except Exception as e:
        print("kw_" + label, "ERR", str(e)[:100])

# ---- extra homodimers (pool was thin: only 7 valid) ----
nodes = [DEP,
         T("rcsb_entry_info.polymer_entity_count_protein", "equals", 1),
         T("rcsb_entry_info.polymer_entity_count_nucleic_acid", "equals", 0),
         T("rcsb_entry_info.deposited_polymer_entity_instance_count", "equals", 2),
         T("rcsb_entry_info.deposited_polymer_monomer_count", "greater_or_equal", 220),
         T("rcsb_entry_info.deposited_polymer_monomer_count", "less_or_equal", 900),
         T("rcsb_entry_info.resolution_combined", "less_or_equal", 2.6)]
total, ids = search(nodes, rows=100)
pools["homodimer2"] = ids
print(f"homodimer2 total={total} kept={len(ids)}")

# ---- extra protein-RNA (want variety beyond CRISPR/LGP2) ----
nodes = [DEP,
         T("rcsb_entry_info.polymer_entity_count_RNA", "greater_or_equal", 1),
         T("rcsb_entry_info.polymer_entity_count_protein", "greater_or_equal", 1),
         T("rcsb_entry_info.deposited_polymer_monomer_count", "greater_or_equal", 120),
         T("rcsb_entry_info.deposited_polymer_monomer_count", "less_or_equal", 1400)]
total, ids = search(nodes, rows=100, sort_attr="rcsb_accession_info.deposit_date")
pools["protein_rna2"] = ids
print(f"protein_rna2 total={total} kept={len(ids)}")

# ---- nucleic-acid-only, more variety ----
nodes = [DEP,
         T("rcsb_entry_info.polymer_entity_count_protein", "equals", 0),
         T("rcsb_entry_info.polymer_entity_count_nucleic_acid", "greater_or_equal", 1),
         T("rcsb_entry_info.deposited_polymer_monomer_count", "greater_or_equal", 25),
         T("rcsb_entry_info.deposited_polymer_monomer_count", "less_or_equal", 300)]
total, ids = search(nodes, rows=100, sort_attr="rcsb_accession_info.deposit_date")
pools["nucleic_only2"] = ids
print(f"nucleic_only2 total={total} kept={len(ids)}")

# ---- enzyme + substrate-analog heterocomplexes, non-immune ----
nodes = [DEP,
         T("rcsb_entry_info.polymer_entity_count_protein", "equals", 2),
         T("rcsb_entry_info.polymer_entity_count_nucleic_acid", "equals", 0),
         T("struct_keywords.pdbx_keywords", "contains_phrase", "IMMUNE SYSTEM"),
         T("rcsb_entry_info.deposited_polymer_monomer_count", "greater_or_equal", 200),
         T("rcsb_entry_info.deposited_polymer_monomer_count", "less_or_equal", 1100)]
# note: this one is *immune* on purpose, to enrich nanobody pool
total, ids = search(nodes, rows=60)
pools["nanobody2"] = ids
print(f"nanobody2 total={total} kept={len(ids)}")

old = json.load(open("pools.json"))
old.update(pools)
json.dump(old, open("pools.json", "w"), indent=1)

known = set(json.load(open("ids.json")))
newids = sorted({i for v in pools.values() for i in v} - known)
json.dump(newids, open("ids_new.json", "w"))
print("new ids to fetch:", len(newids))
