"""Stage 16: screen all 45 targets for MODIFIED RESIDUES.

Why this is needed: every earlier stage read `pdbx_seq_one_letter_code_can`, the
*canonical* sequence, in which modified residues are mapped back to their parent
amino acid (SEP->S, TPO->T, PTR->Y, MSE->M ...). So a phosphopeptide passes a
"standard 20 residues" filter with the phosphates silently erased. #22 9T9W is
"diphosphorylated I-kappa-B-alpha degron peptide" -- the modification IS the
recognition event, and the generated job file drops it.

The non-canonical `pdbx_seq_one_letter_code` writes modified residues as
parenthesised CCD codes, e.g. "...GL(SEP)DES(TPO)PQ...", so it recovers both the
identity and the 1-based position that the Server's `modifications` field wants.
"""
import json
import re
import time
import urllib.request

OPENER = urllib.request.build_opener(urllib.request.ProxyHandler({}))
GRAPHQL = "https://data.rcsb.org/graphql"

# ptmType codes the Server accepts (server/README.md), protein chains
SERVER_PTM = {
    "SEP", "TPO", "PTR", "NEP", "HIP", "ALY", "MLY", "M3L", "MLZ", "2MR", "AGM",
    "MCS", "HYP", "HY3", "LYZ", "AHB", "P1L", "SNN", "SNC", "TRF", "KCR", "CIR",
    "YHA",
}
SERVER_DNA_MOD = {"5CM", "C34", "5HC", "6OG", "6MA", "1CC", "8OG", "5FC", "3DR"}
SERVER_RNA_MOD = {"PSU", "5MC", "OMC", "4OC", "5MU", "OMU", "UR3", "A2M", "MA6",
                  "6MZ", "2MG", "OMG", "7MG", "RSQ"}
# selenomethionine is a crystallographic substitute for Met, not a real PTM:
# entering plain M is correct, so it must not be reported as a missing feature
CRYST_SUBSTITUTE = {"MSE", "MLE"}

Q = """
query($ids:[String!]!){
  entries(entry_ids:$ids){
    rcsb_id
    polymer_entities {
      rcsb_id
      entity_poly {
        pdbx_seq_one_letter_code
        pdbx_seq_one_letter_code_can
        rcsb_entity_polymer_type
      }
      rcsb_polymer_entity { pdbx_description }
      rcsb_polymer_entity_container_identifiers {
        auth_asym_ids
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


def parse_modified(raw):
    """Return [(position_1based, code)] from a non-canonical one-letter sequence."""
    if not raw:
        return []
    s = raw.replace("\n", "").strip()
    out = []
    pos = 0
    i = 0
    while i < len(s):
        if s[i] == "(":
            j = s.find(")", i)
            if j == -1:
                break
            pos += 1
            out.append((pos, s[i + 1:j]))
            i = j + 1
        else:
            pos += 1
            i += 1
    return out


targets = json.load(open(
    "/Users/lilindu/pdb-structures-for-homework/alphafold_homework/targets.json"))
ids = [r["id"] for r in targets]

raw = {}
B = 20
for i in range(0, len(ids), B):
    chunk = ids[i:i + B]
    for a in range(4):
        try:
            d = post({"query": Q, "variables": {"ids": chunk}})
            if "data" not in d:
                raise RuntimeError(str(d.get("errors"))[:200])
            for e in d["data"]["entries"] or []:
                raw[e["rcsb_id"]] = e.get("polymer_entities") or []
            break
        except Exception as ex:
            print("retry", i, str(ex)[:110])
            time.sleep(4)

json.dump(raw, open("ptm_raw.json", "w"))

report = {}
for rec in targets:
    eid = rec["id"]
    found = []
    for pe in raw.get(eid) or []:
        ep = pe.get("entity_poly") or {}
        ptype = ep.get("rcsb_entity_polymer_type") or "?"
        mods = parse_modified(ep.get("pdbx_seq_one_letter_code"))
        if not mods:
            continue
        ent = pe["rcsb_id"].split("_")[-1]
        desc = ((pe.get("rcsb_polymer_entity") or {}).get("pdbx_description") or "")
        for pos, code in mods:
            if code in CRYST_SUBSTITUTE:
                status = "substitute"
            elif ptype == "Protein":
                status = "supported" if code in SERVER_PTM else "UNSUPPORTED"
            elif ptype == "DNA":
                status = "supported" if code in SERVER_DNA_MOD else "UNSUPPORTED"
            elif ptype == "RNA":
                status = "supported" if code in SERVER_RNA_MOD else "UNSUPPORTED"
            else:
                status = "UNSUPPORTED"
            found.append({"entity": ent, "type": ptype, "pos": pos,
                          "code": code, "status": status, "desc": desc})
    if found:
        report[eid] = found

json.dump(report, open("ptm_report.json", "w"), indent=1)

print(f"targets carrying modified residues: {len(report)}/45\n")
for rec in targets:
    f = report.get(rec["id"])
    if not f:
        continue
    real = [x for x in f if x["status"] != "substitute"]
    tag = "  <-- REAL PTM, currently missing from job file" if real else "  (MSE only)"
    print(f"#{rec['no']:2d} {rec['id']} [{rec['group']}]{tag}")
    print(f"     {rec['title'][:100]}")
    from collections import Counter
    for st in ("supported", "UNSUPPORTED", "substitute"):
        items = [x for x in f if x["status"] == st]
        if not items:
            continue
        c = Counter(x["code"] for x in items)
        detail = ", ".join(f"{k}×{v}" for k, v in c.items())
        if st == "substitute":
            print(f"     {st}: {detail}")
        else:
            sites = "; ".join(f"{x['code']}@{x['pos']}(entity {x['entity']},{x['type']})"
                              for x in items[:8])
            print(f"     {st}: {detail}  ->  {sites}")
