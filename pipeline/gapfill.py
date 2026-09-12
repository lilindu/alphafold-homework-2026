"""Stage 9: fill the thin oligomeric pools using assembly-level symmetry
annotation (rcsb_struct_symmetry.oligomeric_state), which is the correct
criterion -- ASU chain counts are not the biological unit."""
import json
import urllib.request

OPENER = urllib.request.build_opener(urllib.request.ProxyHandler({}))
SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"


def post(p, timeout=180):
    r = urllib.request.Request(SEARCH, data=json.dumps(p).encode(),
                               headers={"Content-Type": "application/json"})
    return json.load(OPENER.open(r, timeout=timeout))


def T(a, o, v):
    return {"type": "terminal", "service": "text",
            "parameters": {"attribute": a, "operator": o, "value": v}}


DEP = T("rcsb_accession_info.deposit_date", "greater", "2025-02-03")
SYM = "rcsb_struct_symmetry.oligomeric_state"


def run(nodes, rows=150, sort="rcsb_entry_info.resolution_combined"):
    p = {"query": {"type": "group", "logical_operator": "and", "nodes": nodes},
         "return_type": "entry",
         "request_options": {"paginate": {"start": 0, "rows": rows},
                             "results_content_type": ["experimental"],
                             "sort": [{"sort_by": sort, "direction": "asc"}]}}
    d = post(p)
    return d.get("total_count", 0), [x["identifier"] for x in d.get("result_set", [])]


pools = json.load(open("pools.json"))
new = {}

# homodimers by assembly symmetry
t, ids = run([DEP, T(SYM, "exact_match", "Homo 2-mer"),
              T("rcsb_entry_info.polymer_entity_count_protein", "equals", 1),
              T("rcsb_entry_info.polymer_entity_count_nucleic_acid", "equals", 0),
              T("rcsb_entry_info.deposited_polymer_monomer_count", "greater_or_equal", 230),
              T("rcsb_entry_info.deposited_polymer_monomer_count", "less_or_equal", 900),
              T("rcsb_entry_info.resolution_combined", "less_or_equal", 2.2)], rows=150)
new["homodimer_sym"] = ids
print(f"Homo 2-mer: total={t} kept={len(ids)}")

# homo 4-mer and larger
for state in ["Homo 4-mer", "Homo 3-mer", "Homo 6-mer", "Homo 8-mer"]:
    t, ids = run([DEP, T(SYM, "exact_match", state),
                  T("rcsb_entry_info.polymer_entity_count_protein", "equals", 1),
                  T("rcsb_entry_info.polymer_entity_count_nucleic_acid", "equals", 0),
                  T("rcsb_entry_info.deposited_polymer_monomer_count", "greater_or_equal", 300),
                  T("rcsb_entry_info.deposited_polymer_monomer_count", "less_or_equal", 1800),
                  T("rcsb_entry_info.resolution_combined", "less_or_equal", 2.4)], rows=80)
    key = "oligo_" + state.replace(" ", "_")
    new[key] = ids
    print(f"{state}: total={t} kept={len(ids)}")

# monomer by assembly symmetry, to strengthen A1
t, ids = run([DEP, T(SYM, "exact_match", "Monomer"),
              T("rcsb_entry_info.polymer_entity_count_protein", "equals", 1),
              T("rcsb_entry_info.polymer_entity_count_nucleic_acid", "equals", 0),
              T("rcsb_entry_info.deposited_polymer_monomer_count", "greater_or_equal", 130),
              T("rcsb_entry_info.deposited_polymer_monomer_count", "less_or_equal", 420),
              T("rcsb_entry_info.resolution_combined", "less_or_equal", 1.9)], rows=120)
new["mono_sym"] = ids
print(f"Monomer: total={t} kept={len(ids)}")

pools.update(new)
json.dump(pools, open("pools.json", "w"), indent=1)

known = set(json.load(open("meta_all.json")))
todo = sorted({i for v in new.values() for i in v} - known)
json.dump(todo, open("ids_gap.json", "w"))
print("new ids needing metadata:", len(todo))
