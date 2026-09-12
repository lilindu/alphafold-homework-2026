"""Stage 2: fetch entry metadata + ligand/ion inventory via RCSB GraphQL."""
import json
import time
import urllib.request

OPENER = urllib.request.build_opener(urllib.request.ProxyHandler({}))
GRAPHQL = "https://data.rcsb.org/graphql"

Q = """
query($ids:[String!]!){
  entries(entry_ids:$ids){
    rcsb_id
    struct { title }
    rcsb_accession_info { deposit_date initial_release_date }
    exptl { method }
    rcsb_entry_info {
      resolution_combined
      deposited_polymer_monomer_count
      polymer_entity_count_protein
      polymer_entity_count_DNA
      polymer_entity_count_RNA
      deposited_polymer_entity_instance_count
      deposited_nonpolymer_entity_instance_count
      selected_polymer_entity_types
    }
    struct_keywords { pdbx_keywords }
    polymer_entities {
      rcsb_id
      entity_poly { pdbx_seq_one_letter_code_can rcsb_entity_polymer_type }
      rcsb_polymer_entity_container_identifiers { auth_asym_ids }
      rcsb_polymer_entity { pdbx_description }
      rcsb_entity_source_organism { ncbi_scientific_name }
    }
    nonpolymer_entities {
      rcsb_id
      nonpolymer_comp {
        chem_comp { id name formula type formula_weight }
        rcsb_chem_comp_descriptor { SMILES }
      }
      rcsb_nonpolymer_entity_container_identifiers { auth_asym_ids }
    }
  }
}
"""


def post(payload, timeout=180):
    req = urllib.request.Request(
        GRAPHQL, data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json"})
    return json.load(OPENER.open(req, timeout=timeout))


ids = json.load(open("ids.json"))
out = {}
B = 40
for i in range(0, len(ids), B):
    chunk = ids[i:i + B]
    for attempt in range(3):
        try:
            d = post({"query": Q, "variables": {"ids": chunk}})
            for e in d["data"]["entries"] or []:
                out[e["rcsb_id"]] = e
            break
        except Exception as ex:
            print("retry", i, ex)
            time.sleep(3)
    print(f"{i + len(chunk)}/{len(ids)}", flush=True)

json.dump(out, open("meta.json", "w"))
print("fetched", len(out))
