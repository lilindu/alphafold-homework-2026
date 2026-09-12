"""Fetch biological-assembly composition for every candidate.

The asymmetric unit is NOT the biological unit: 9QDV has two chains in the ASU
but two separate monomeric assemblies. Oligomeric-state groups must be filtered
on assembly composition instead.
"""
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
    assemblies {
      rcsb_id
      rcsb_assembly_info {
        polymer_entity_instance_count
        polymer_entity_count
        polymer_entity_instance_count_protein
        polymer_entity_instance_count_DNA
        polymer_entity_instance_count_RNA
      }
    }
  }
}
"""


def post(payload, timeout=180):
    req = urllib.request.Request(
        GRAPHQL, data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json"})
    return json.load(OPENER.open(req, timeout=timeout))


valid = json.load(open("valid.json"))
want = sorted(valid)
out = json.load(open("asm.json")) if os.path.exists("asm.json") else {}
todo = [i for i in want if i not in out]
print(f"valid {len(want)}, have {len(out)}, fetching {len(todo)}")

B = 40
for i in range(0, len(todo), B):
    chunk = todo[i:i + B]
    for a in range(4):
        try:
            d = post({"query": Q, "variables": {"ids": chunk}})
            if "data" not in d:
                raise RuntimeError(str(d)[:200])
            for e in d["data"]["entries"] or []:
                out[e["rcsb_id"]] = e.get("assemblies") or []
            break
        except Exception as ex:
            print("retry", i, str(ex)[:110])
            time.sleep(4)
    if (i // B) % 5 == 0:
        json.dump(out, open("asm.json", "w"))
        print(f"  {i + len(chunk)}/{len(todo)}", flush=True)

json.dump(out, open("asm.json", "w"))
print("assemblies cached:", len(out))
