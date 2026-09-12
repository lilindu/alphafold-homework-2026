"""Stage 10: derive the exact per-entity stoichiometry of biological assembly 1.

Why this matters: the asymmetric unit is often several copies of the complex.
9U19's ASU holds protein x2 + two DNA strands x2 (i.e. two copies of a 1:1:1
complex); entering that verbatim would double every chain. Assembly 1's
per-entity counts give the right input spec.

Consistency check: if the per-instance list does not sum to the reported
polymer_entity_instance_count, the assembly was built with symmetry operators
and the list is incomplete -- those entries are flagged for manual handling
rather than silently trusted.
"""
import json
import time
import urllib.request
from collections import Counter

OPENER = urllib.request.build_opener(urllib.request.ProxyHandler({}))
GRAPHQL = "https://data.rcsb.org/graphql"

Q = """
query($ids:[String!]!){
  entries(entry_ids:$ids){
    rcsb_id
    assemblies {
      rcsb_id
      rcsb_assembly_info { polymer_entity_instance_count }
      polymer_entity_instances {
        rcsb_polymer_entity_instance_container_identifiers { entity_id auth_asym_id }
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


chosen = json.load(open("chosen45.json"))
ids = [r["id"] for r in chosen]

raw = {}
B = 20
for i in range(0, len(ids), B):
    chunk = ids[i:i + B]
    for a in range(4):
        try:
            d = post({"query": Q, "variables": {"ids": chunk}})
            if "data" not in d:
                raise RuntimeError(str(d)[:200])
            for e in d["data"]["entries"] or []:
                raw[e["rcsb_id"]] = e.get("assemblies") or []
            break
        except Exception as ex:
            print("retry", i, str(ex)[:100])
            time.sleep(4)

json.dump(raw, open("stoich_raw.json", "w"))

out = {}
for eid, asms in raw.items():
    if not asms:
        out[eid] = {"ok": False, "reason": "no assembly"}
        continue
    a = asms[0]
    reported = (a.get("rcsb_assembly_info") or {}).get("polymer_entity_instance_count") or 0
    insts = a.get("polymer_entity_instances") or []
    c = Counter(
        x["rcsb_polymer_entity_instance_container_identifiers"]["entity_id"]
        for x in insts
        if x.get("rcsb_polymer_entity_instance_container_identifiers"))
    listed = sum(c.values())
    out[eid] = {
        "ok": listed == reported and listed > 0,
        "assembly": a.get("rcsb_id"),
        "reported": reported,
        "listed": listed,
        "per_entity": dict(c),
        "n_assemblies": len(asms),
    }

json.dump(out, open("stoich.json", "w"), indent=1)

bad = {k: v for k, v in out.items() if not v["ok"]}
print(f"consistent: {len(out) - len(bad)}/{len(out)}")
for k, v in bad.items():
    print("  FLAG", k, v)
