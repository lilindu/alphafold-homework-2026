"""Census phase 1: bulk metadata for all census candidates (checkpointed)."""
import json, os, time, urllib.request

OP = urllib.request.build_opener(urllib.request.ProxyHandler({}))
G = "https://data.rcsb.org/graphql"

Q = """
query($ids:[String!]!){
 entries(entry_ids:$ids){
  rcsb_id
  struct{title}
  rcsb_accession_info{deposit_date initial_release_date}
  exptl{method}
  rcsb_entry_info{resolution_combined deposited_polymer_monomer_count
    polymer_entity_count_protein polymer_entity_count_DNA polymer_entity_count_RNA
    branched_entity_count}
  struct_keywords{pdbx_keywords}
  polymer_entities{
   rcsb_id
   entity_poly{pdbx_seq_one_letter_code pdbx_seq_one_letter_code_can rcsb_entity_polymer_type}
   rcsb_polymer_entity_container_identifiers{auth_asym_ids}
   rcsb_polymer_entity{pdbx_description}
   rcsb_entity_source_organism{ncbi_scientific_name}
   rcsb_polymer_entity_annotation{type}
   uniprots{rcsb_id}}
  nonpolymer_entities{
   rcsb_id
   nonpolymer_comp{chem_comp{id name formula_weight}}
   rcsb_nonpolymer_entity_container_identifiers{auth_asym_ids}}
  branched_entities{
   rcsb_id
   rcsb_branched_entity{pdbx_description}
   rcsb_branched_entity_container_identifiers{chem_comp_monomers auth_asym_ids}}}}
"""


def post(p, t=300):
    r = urllib.request.Request(G, data=json.dumps(p).encode(),
                               headers={"Content-Type": "application/json"})
    return json.load(OP.open(r, timeout=t))


ids = json.load(open("census_ids.json"))
out = json.load(open("cmeta.json")) if os.path.exists("cmeta.json") else {}
todo = [i for i in ids if i not in out]
print(f"census ids {len(ids)}, cached {len(out)}, fetching {len(todo)}", flush=True)

B = 25
for i in range(0, len(todo), B):
    ch = todo[i:i + B]
    for a in range(4):
        try:
            d = post({"query": Q, "variables": {"ids": ch}})
            if "data" not in d:
                raise RuntimeError(str(d.get("errors"))[:200])
            for e in d["data"]["entries"] or []:
                out[e["rcsb_id"]] = e
            break
        except Exception as ex:
            print("retry", i, str(ex)[:90], flush=True)
            time.sleep(4)
    if (i // B) % 20 == 0:
        json.dump(out, open("cmeta.json", "w"))
        print(f"  {i + len(ch)}/{len(todo)}", flush=True)

json.dump(out, open("cmeta.json", "w"))
print("cmeta entries:", len(out))

# top up the CCD cache for any ligand code we have not dated yet
ccd = json.load(open("ccd.json"))
codes = set()
for e in out.values():
    for ne in e.get("nonpolymer_entities") or []:
        cc = (ne.get("nonpolymer_comp") or {}).get("chem_comp") or {}
        if cc.get("id"):
            codes.add(cc["id"])
todo2 = sorted(codes - set(ccd))
print("new ligand codes to date:", len(todo2), flush=True)
QC = ("query($ids:[String!]!){chem_comps(comp_ids:$ids){rcsb_id "
      "chem_comp{id name formula formula_weight} "
      "rcsb_chem_comp_info{initial_release_date atom_count atom_count_heavy}}}")
for i in range(0, len(todo2), 60):
    ch = todo2[i:i + 60]
    for a in range(4):
        try:
            d = post({"query": QC, "variables": {"ids": ch}})
            for c in d["data"]["chem_comps"] or []:
                ccd[c["rcsb_id"]] = c
            break
        except Exception as ex:
            print("retry ccd", str(ex)[:80], flush=True)
            time.sleep(4)
json.dump(ccd, open("ccd.json", "w"))
print("ccd cache:", len(ccd))
