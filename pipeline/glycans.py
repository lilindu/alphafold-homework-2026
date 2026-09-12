"""Stage 15: work out, per target, exactly how each sugar must be entered.

Three distinct cases, and they need different instructions:

  1. Branched entity covalently linked to Asn/Ser/Thr (role=N-Glycosylation or
     O-Glycosylation in rcsb_branched_struct_conn) -> enter via the proteinChain
     `glycans` field, as residues + position. Server supports only
     BGC/BMA/GLC/MAN/NAG on Asn (plus FUC on Ser/Thr), max 8 residues, and
     cannot be told which atoms form the bond.
  2. Branched entity NOT linked to protein (a free oligosaccharide sitting in a
     binding site, e.g. trehalose in 9R6D) -> there is no way to enter this. The
     `glycans` field requires an attachment position; a free disaccharide is
     neither a supported ion nor enterable as a plain CCD ligand pair.
  3. Single sugar deposited as a nonpolymer entity -> if covalently bonded to a
     residue, same as case 1; if free, it is an ordinary CCD ligand.

The attachment position from PDB is label_seq_id (1-based into the full entity
sequence), which is what the Server's `position` field wants, but the deposited
construct may start partway into the sequence -- so we verify the residue at that
index really is N/S/T before emitting an instruction.
"""
import json
import time
import urllib.request

OPENER = urllib.request.build_opener(urllib.request.ProxyHandler({}))
GRAPHQL = "https://data.rcsb.org/graphql"

SERVER_GLYCAN_ON_ASN = {"BGC", "BMA", "GLC", "MAN", "NAG"}
SERVER_GLYCAN_ON_ST = {"BGC", "BMA", "FUC", "GLC", "MAN", "NAG"}
MAX_GLYCAN_RESIDUES = 8

Q = """
query($ids:[String!]!){
  entries(entry_ids:$ids){
    rcsb_id
    nonpolymer_entities {
      rcsb_id
      nonpolymer_comp { chem_comp { id } }
      nonpolymer_entity_instances {
        rcsb_id
        rcsb_nonpolymer_entity_instance_container_identifiers { auth_asym_id }
        rcsb_nonpolymer_struct_conn {
          role connect_type
          connect_partner { label_comp_id label_seq_id label_atom_id }
          connect_target { label_comp_id label_seq_id auth_seq_id label_atom_id }
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


targets = json.load(open(
    "/Users/lilindu/pdb-structures-for-homework/alphafold_homework/targets.json"))
branched = json.load(open("branched.json"))
ids = [r["id"] for r in targets]

# nonpolymer covalent links (to spot covalently bound single sugars / inhibitors)
npconn = {}
B = 20
for i in range(0, len(ids), B):
    chunk = ids[i:i + B]
    for a in range(4):
        try:
            d = post({"query": Q, "variables": {"ids": chunk}})
            if "data" not in d:
                raise RuntimeError(str(d.get("errors"))[:200])
            for e in d["data"]["entries"] or []:
                npconn[e["rcsb_id"]] = e.get("nonpolymer_entities") or []
            break
        except Exception as ex:
            print("retry", i, str(ex)[:110])
            time.sleep(4)
json.dump(npconn, open("npconn.json", "w"))


def seq_of(rec, entity):
    for c in rec["chains"]:
        if c["entity"] == entity:
            return c["seq"]
    return ""


def protein_entities(rec):
    return [c for c in rec["chains"] if c["type"] == "Protein"]


report = {}
for rec in targets:
    eid = rec["id"]
    info = {"glycan_attached": [], "glycan_free": [], "covalent_nonpolymer": []}

    # --- branched entities ---
    for e in (branched.get(eid) or {}).get("entities", []) or []:
        ci = e["rcsb_branched_entity_container_identifiers"]
        monomers = ci.get("chem_comp_monomers") or []
        desc = (e.get("rcsb_branched_entity") or {}).get("pdbx_description") or ""
        linked = None
        n_res = 0
        for bi in e.get("branched_entity_instances") or []:
            conns = bi.get("rcsb_branched_struct_conn") or []
            # count sugar residues by distinct auth_seq_id on the target side
            seqs = set()
            for sc in conns:
                tg = sc.get("connect_target") or {}
                if tg.get("auth_seq_id") is not None:
                    seqs.add(tg["auth_seq_id"])
                role = (sc.get("role") or "")
                if "Glycosylation" in role:
                    p = sc.get("connect_partner") or {}
                    linked = {"role": role,
                              "residue": p.get("label_comp_id"),
                              "seq_id": p.get("label_seq_id")}
            n_res = max(n_res, len(seqs))
        entry = {"monomers": monomers, "desc": desc, "n_residues": n_res}
        if linked:
            entry.update(linked)
            # verify the attachment index against the entity sequence
            ok_site, site_note = False, ""
            for c in protein_entities(rec):
                s = c["seq"]
                idx = linked.get("seq_id")
                if idx and 1 <= idx <= len(s):
                    aa = s[idx - 1]
                    exp = {"ASN": "N", "SER": "S", "THR": "T"}.get(linked["residue"] or "")
                    if exp and aa == exp:
                        ok_site = True
                        entry["position"] = idx
                        entry["attach_aa"] = aa
                        entry["attach_entity"] = c["entity"]
                        break
                    site_note = f"index {idx} of entity {c['entity']} is {aa}, expected {exp}"
            entry["site_verified"] = ok_site
            if not ok_site:
                entry["site_note"] = site_note or "attachment index not found in any protein entity"
            allowed = (SERVER_GLYCAN_ON_ASN if entry.get("attach_aa") == "N"
                       else SERVER_GLYCAN_ON_ST)
            entry["monomers_supported"] = all(m in allowed for m in monomers)
            entry["size_ok"] = n_res <= MAX_GLYCAN_RESIDUES
            info["glycan_attached"].append(entry)
        else:
            info["glycan_free"].append(entry)

    # --- covalently bonded nonpolymer components ---
    for ne in npconn.get(eid) or []:
        code = ((ne.get("nonpolymer_comp") or {}).get("chem_comp") or {}).get("id")
        for ni in ne.get("nonpolymer_entity_instances") or []:
            for sc in ni.get("rcsb_nonpolymer_struct_conn") or []:
                if sc.get("connect_type") == "covalent bond":
                    p = sc.get("connect_partner") or {}
                    t = sc.get("connect_target") or {}
                    info["covalent_nonpolymer"].append({
                        "code": code, "role": sc.get("role"),
                        "partner": p.get("label_comp_id"),
                        "partner_seq": p.get("label_seq_id"),
                        "target": t.get("label_comp_id")})
    report[eid] = info

json.dump(report, open("glycan_report.json", "w"), indent=1)

for eid, info in report.items():
    if any(info.values()):
        print("===", eid)
        for g in info["glycan_attached"]:
            print(f"  ATTACHED {g['monomers']} n={g['n_residues']} on {g.get('residue')}"
                  f" pos={g.get('position')} verified={g.get('site_verified')}"
                  f" supported={g.get('monomers_supported')} size_ok={g.get('size_ok')}"
                  + (f" [{g.get('site_note')}]" if g.get("site_note") else ""))
        for g in info["glycan_free"]:
            print(f"  FREE OLIGO {g['monomers']} n={g['n_residues']} -- cannot be entered")
        seen = set()
        for c in info["covalent_nonpolymer"]:
            k = (c["code"], c["partner"], c["partner_seq"])
            if k in seen:
                continue
            seen.add(k)
            print(f"  COVALENT LIGAND {c['code']} -> {c['partner']}{c['partner_seq']}"
                  f" role={c['role']}")
