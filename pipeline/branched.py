"""Stage 14: fetch BRANCHED entities (oligosaccharides) for all candidates.

PDB deposits glycans/oligosaccharides as a third entity type, alongside polymer
and nonpolymer. My original metadata query fetched only polymer + nonpolymer, so
sugars deposited as branched entities were invisible: 9R6D was selected into the
"no ligand" tier while actually carrying a bound trehalose.

On AlphaFold Server a glycan is NOT a free ligand -- it is attached to a residue
through the proteinChain `glycans` field (residues + position). So we need both
the monomer composition and the attachment site.
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
    rcsb_entry_info { branched_entity_count }
    branched_entities {
      rcsb_id
      pdbx_entity_branch { type }
      rcsb_branched_entity { pdbx_description }
      rcsb_branched_entity_container_identifiers { chem_comp_monomers auth_asym_ids }
      pdbx_entity_branch_descriptor { type descriptor }
      branched_entity_instances {
        rcsb_id
        rcsb_branched_entity_instance_container_identifiers { auth_asym_id asym_id }
        rcsb_branched_struct_conn {
          role
          connect_type
          connect_partner { label_asym_id label_comp_id label_seq_id label_atom_id }
          connect_target { label_asym_id label_comp_id label_seq_id auth_seq_id label_atom_id }
        }
      }
    }
  }
}
"""


def post(payload, timeout=240):
    req = urllib.request.Request(
        GRAPHQL, data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json"})
    return json.load(OPENER.open(req, timeout=timeout))


valid = json.load(open("valid.json"))
want = sorted(valid)
out = json.load(open("branched.json")) if os.path.exists("branched.json") else {}
todo = [i for i in want if i not in out]
print(f"valid {len(want)}, cached {len(out)}, fetching {len(todo)}")

B = 25
for i in range(0, len(todo), B):
    chunk = todo[i:i + B]
    for a in range(4):
        try:
            d = post({"query": Q, "variables": {"ids": chunk}})
            if "data" not in d:
                raise RuntimeError(str(d.get("errors"))[:200])
            for e in d["data"]["entries"] or []:
                out[e["rcsb_id"]] = {
                    "count": (e.get("rcsb_entry_info") or {}).get("branched_entity_count") or 0,
                    "entities": e.get("branched_entities") or [],
                }
            break
        except Exception as ex:
            print("retry", i, str(ex)[:110])
            time.sleep(4)
    if (i // B) % 10 == 0:
        json.dump(out, open("branched.json", "w"))
        print(f"  {i + len(chunk)}/{len(todo)}", flush=True)

json.dump(out, open("branched.json", "w"))
withsugar = {k: v for k, v in out.items() if v["count"]}
print("cached:", len(out), "| entries with branched glycans:", len(withsugar))
