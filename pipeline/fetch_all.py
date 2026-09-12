"""Fetch metadata for every id in pools.json into meta_all.json (incremental)."""
import json
import os
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
    }
    struct_keywords { pdbx_keywords }
    polymer_entities {
      rcsb_id
      entity_poly { pdbx_seq_one_letter_code_can rcsb_entity_polymer_type }
      rcsb_polymer_entity_container_identifiers { auth_asym_ids }
      rcsb_polymer_entity { pdbx_description }
      rcsb_entity_source_organism { ncbi_scientific_name }
      rcsb_polymer_entity_annotation { type name }
      uniprots { rcsb_id }
    }
    nonpolymer_entities {
      rcsb_id
      nonpolymer_comp { chem_comp { id name formula formula_weight } }
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


pools = json.load(open("pools.json"))
want = sorted({i for v in pools.values() for i in v})

out = {}
if os.path.exists("meta_all.json"):
    out = json.load(open("meta_all.json"))

todo = [i for i in want if i not in out]
print(f"target {len(want)}, have {len(out)}, fetching {len(todo)}")

B = 40
for i in range(0, len(todo), B):
    chunk = todo[i:i + B]
    for attempt in range(4):
        try:
            d = post({"query": Q, "variables": {"ids": chunk}})
            for e in d["data"]["entries"] or []:
                out[e["rcsb_id"]] = e
            break
        except Exception as ex:
            print("retry", i, str(ex)[:90])
            time.sleep(4)
    if (i // B) % 5 == 0:
        json.dump(out, open("meta_all.json", "w"))
        print(f"  {i + len(chunk)}/{len(todo)}", flush=True)

json.dump(out, open("meta_all.json", "w"))
print("total meta:", len(out))
