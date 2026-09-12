"""Stage 1: build candidate pools from RCSB search, save to pools.json."""
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


def T(attr, op, val, negation=False):
    p = {"attribute": attr, "operator": op, "value": val}
    if negation:
        p["negation"] = True
    return {"type": "terminal", "service": "text", "parameters": p}


def search(nodes, rows=50, sort_attr="rcsb_entry_info.resolution_combined",
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


DEP = T("rcsb_accession_info.deposit_date", "greater", DEPOSIT_AFTER)


def sz(lo, hi):
    return [T("rcsb_entry_info.deposited_polymer_monomer_count", "greater_or_equal", lo),
            T("rcsb_entry_info.deposited_polymer_monomer_count", "less_or_equal", hi)]


def res(mx):
    return [T("rcsb_entry_info.resolution_combined", "less_or_equal", mx)]


CATS = {
    # ---------- easy tier ----------
    "mono_protein": [DEP,
                     T("rcsb_entry_info.polymer_entity_count_protein", "equals", 1),
                     T("rcsb_entry_info.polymer_entity_count_nucleic_acid", "equals", 0),
                     T("rcsb_entry_info.deposited_polymer_entity_instance_count", "equals", 1),
                     *sz(120, 450), *res(2.0)],
    "mono_with_cofactor": [DEP,
                           T("rcsb_entry_info.polymer_entity_count_protein", "equals", 1),
                           T("rcsb_entry_info.polymer_entity_count_nucleic_acid", "equals", 0),
                           T("rcsb_entry_info.deposited_polymer_entity_instance_count", "equals", 1),
                           T("rcsb_entry_info.deposited_nonpolymer_entity_instance_count", "greater_or_equal", 1),
                           *sz(150, 500), *res(2.2)],
    "homodimer": [DEP,
                  T("rcsb_entry_info.polymer_entity_count_protein", "equals", 1),
                  T("rcsb_entry_info.polymer_entity_count_nucleic_acid", "equals", 0),
                  T("rcsb_entry_info.deposited_polymer_entity_instance_count", "equals", 2),
                  *sz(200, 900), *res(2.3)],
    "homotetramer": [DEP,
                     T("rcsb_entry_info.polymer_entity_count_protein", "equals", 1),
                     T("rcsb_entry_info.polymer_entity_count_nucleic_acid", "equals", 0),
                     T("rcsb_entry_info.deposited_polymer_entity_instance_count", "equals", 4),
                     *sz(300, 1600), *res(2.5)],
    # ---------- medium tier ----------
    "hetero2": [DEP,
                T("rcsb_entry_info.polymer_entity_count_protein", "equals", 2),
                T("rcsb_entry_info.polymer_entity_count_nucleic_acid", "equals", 0),
                *sz(200, 1100), *res(2.8)],
    "hetero3": [DEP,
                T("rcsb_entry_info.polymer_entity_count_protein", "equals", 3),
                T("rcsb_entry_info.polymer_entity_count_nucleic_acid", "equals", 0),
                *sz(300, 1600), *res(3.2)],
    "peptide_complex": [DEP,
                        T("rcsb_entry_info.polymer_entity_count_protein", "equals", 2),
                        T("rcsb_entry_info.polymer_entity_count_nucleic_acid", "equals", 0),
                        *sz(150, 700), *res(2.5)],
    "protein_dna": [DEP,
                    T("rcsb_entry_info.polymer_entity_count_DNA", "greater_or_equal", 1),
                    T("rcsb_entry_info.polymer_entity_count_protein", "greater_or_equal", 1),
                    *sz(100, 1400)],
    "protein_rna": [DEP,
                    T("rcsb_entry_info.polymer_entity_count_RNA", "greater_or_equal", 1),
                    T("rcsb_entry_info.polymer_entity_count_protein", "greater_or_equal", 1),
                    *sz(100, 1400)],
    # ---------- hard tier ----------
    "antibody": [DEP,
                 T("struct_keywords.pdbx_keywords", "contains_phrase", "IMMUNE SYSTEM"),
                 T("rcsb_entry_info.polymer_entity_count_protein", "greater_or_equal", 3),
                 *sz(400, 1400), *res(3.2)],
    "nanobody": [DEP,
                 T("struct_keywords.pdbx_keywords", "contains_phrase", "IMMUNE SYSTEM"),
                 T("rcsb_entry_info.polymer_entity_count_protein", "equals", 2),
                 *sz(250, 900), *res(3.0)],
    "membrane_kw": [DEP,
                    T("struct_keywords.pdbx_keywords", "contains_phrase", "MEMBRANE PROTEIN"),
                    *sz(200, 2200)],
    "transport_kw": [DEP,
                     T("struct_keywords.pdbx_keywords", "contains_phrase", "TRANSPORT PROTEIN"),
                     *sz(200, 2200)],
    "signaling_gpcr": [DEP,
                       T("struct_keywords.pdbx_keywords", "contains_phrase", "SIGNALING PROTEIN"),
                       T("rcsb_entry_info.polymer_entity_count_protein", "greater_or_equal", 3),
                       *sz(600, 2500)],
    "large_assembly": [DEP,
                       T("rcsb_entry_info.polymer_entity_count_protein", "greater_or_equal", 5),
                       T("rcsb_entry_info.polymer_entity_count_nucleic_acid", "equals", 0),
                       *sz(600, 3000)],
    "nucleic_only": [DEP,
                     T("rcsb_entry_info.polymer_entity_count_protein", "equals", 0),
                     T("rcsb_entry_info.polymer_entity_count_nucleic_acid", "greater_or_equal", 1),
                     *sz(20, 250)],
    "de_novo": [DEP,
                T("rcsb_entity_source_organism.ncbi_scientific_name", "exact_match",
                  "synthetic construct"),
                T("rcsb_entry_info.polymer_entity_count_nucleic_acid", "equals", 0),
                *sz(100, 800), *res(2.5)],
    "virus_protein": [DEP,
                      T("struct_keywords.pdbx_keywords", "contains_phrase", "VIRAL PROTEIN"),
                      *sz(150, 1500)],
}

pools = {}
for name, nodes in CATS.items():
    try:
        total, ids = search(nodes, rows=50)
        pools[name] = ids
        print(f"{name}: total={total} kept={len(ids)}")
    except Exception as e:
        print(f"{name}: ERROR {e}")

with open("pools.json", "w") as f:
    json.dump(pools, f, indent=1)

allids = sorted({i for v in pools.values() for i in v})
print("unique entries:", len(allids))
with open("ids.json", "w") as f:
    json.dump(allids, f)
