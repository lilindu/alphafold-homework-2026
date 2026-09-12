"""Stage 11: resolve final input spec per target and emit the assignment pack."""
import json
import os
import re
from collections import Counter

OUT = "/Users/lilindu/pdb-structures-for-homework/alphafold_homework"

chosen = json.load(open("chosen45.json"))
meta = json.load(open("meta_all.json"))
ccd = json.load(open("ccd.json"))
stoich = json.load(open("stoich.json"))
raw = json.load(open("stoich_raw.json"))

BUILTIN_LIG = {"ADP", "ATP", "AMP", "GTP", "GDP", "FAD", "NAD", "NAP", "NDP",
               "HEM", "HEC", "PLM", "OLA", "MYR", "CIT", "CLA", "CHL", "BCL", "BCB"}
BUILTIN_ION = {"MG", "ZN", "CL", "CA", "NA", "MN", "K", "FE", "CU", "CO"}

TIER_NAME = {"A": "较易", "B": "中等", "C": "有挑战"}


def entity_map(eid):
    """entity_id -> (type, seq, desc, org) for one PDB entry."""
    out = {}
    for pe in meta[eid].get("polymer_entities") or []:
        ent = pe["rcsb_id"].split("_")[-1]
        ep = pe.get("entity_poly") or {}
        seq = (ep.get("pdbx_seq_one_letter_code_can") or "").replace("\n", "").strip().upper()
        orgs = [o.get("ncbi_scientific_name")
                for o in (pe.get("rcsb_entity_source_organism") or [])
                if o.get("ncbi_scientific_name")]
        out[ent] = {
            "type": ep.get("rcsb_entity_polymer_type") or "?",
            "seq": seq,
            "desc": ((pe.get("rcsb_polymer_entity") or {}).get("pdbx_description") or "").strip(),
            "org": orgs[0] if orgs else "",
        }
    return out


def resolve_counts(eid):
    """Per-entity copy numbers in biological assembly 1.

    When the instance list is truncated by symmetry operators, scale it up to the
    reported instance count. Returns (counts, note) or (None, reason).
    """
    st = stoich.get(eid) or {}
    per = {k: v for k, v in (st.get("per_entity") or {}).items()}
    if not per:
        return None, "no assembly data"
    reported, listed = st.get("reported") or 0, st.get("listed") or 0
    if listed == reported:
        return per, ""
    if listed and reported and reported % listed == 0:
        f = reported // listed
        return {k: v * f for k, v in per.items()}, f"按对称操作展开 ×{f}"
    return per, f"注意:PDB 报告 {reported} 条链,实例表只列出 {listed} 条,请自行核对拷贝数"


COMP = {"A": "T", "T": "A", "G": "C", "C": "G"}


def is_revcomp(a, b):
    return len(a) == len(b) and all(COMP.get(x) == y for x, y in zip(a, b[::-1]))


def ccd_atoms(code):
    info = (ccd.get(code) or {}).get("rcsb_chem_comp_info") or {}
    return info.get("atom_count_heavy") or info.get("atom_count") or 0


def ccd_name(code):
    return ((ccd.get(code) or {}).get("chem_comp") or {}).get("name", "").strip()


records = []
for r in chosen:
    eid = r["id"]
    ents = entity_map(eid)
    counts, note = resolve_counts(eid)
    if counts is None:
        counts = {k: 1 for k in ents}
        note = "无装配信息,按每种链 1 份处理"

    chains, tokens = [], 0
    for ent, n in sorted(counts.items(), key=lambda kv: -len(ents.get(kv[0], {}).get("seq", ""))):
        e = ents.get(ent)
        if not e or not e["seq"]:
            continue
        chains.append({"entity": ent, "type": e["type"], "seq": e["seq"],
                       "len": len(e["seq"]), "copies": n,
                       "desc": e["desc"], "org": e["org"]})
        tokens += len(e["seq"]) * n

    # double-stranded DNA detection
    dna = [c for c in chains if c["type"] == "DNA"]
    ds_pairs = []
    for i in range(len(dna)):
        for j in range(i + 1, len(dna)):
            if is_revcomp(dna[i]["seq"], dna[j]["seq"]):
                ds_pairs.append((dna[i]["entity"], dna[j]["entity"]))

    ligands, ions = [], []
    for code, n, atoms, name in r["ligands"]:
        ligands.append({"code": code, "count": n, "atoms": atoms or ccd_atoms(code),
                        "name": name or ccd_name(code),
                        "builtin": code in BUILTIN_LIG})
        tokens += (atoms or ccd_atoms(code)) * n
    for code, n in r["ions"]:
        ions.append({"code": code, "count": n, "builtin": code in BUILTIN_ION})
        tokens += n

    records.append({
        "id": eid, "group": r["group"], "group_label": r["group_label"],
        "tier": r["tier"], "group_note": r["group_note"],
        "title": r["title"], "resolution": r["resolution"], "method": r["method"],
        "deposit": r["deposit"], "release": r["release"],
        "chains": chains, "ligands": ligands, "ions": ions,
        "omitted": r["omitted"], "tokens": tokens,
        "assembly_note": note, "ds_pairs": ds_pairs,
        "n_prot": sum(1 for c in chains if c["type"] == "Protein"),
        "n_nuc": sum(1 for c in chains if c["type"] in ("DNA", "RNA")),
    })

order = {"A": 0, "B": 1, "C": 2}
records.sort(key=lambda x: (order[x["tier"]], x["group"], x["id"]))
for n, rec in enumerate(records, 1):
    rec["no"] = n

over = [r for r in records if r["tokens"] > 5000]
print("targets:", len(records), "| token max:", max(r["tokens"] for r in records),
      "| over limit:", len(over))
print("tier spread:", Counter(r["tier"] for r in records))
nonbuiltin = sorted({l["code"] for r in records for l in r["ligands"] if not l["builtin"]})
print("ligands needing manual CCD entry:", len(nonbuiltin))
print(" ", " ".join(nonbuiltin))
badion = sorted({i["code"] for r in records for i in r["ions"] if not i["builtin"]})
print("unsupported ions still present (must be a bug if non-empty):", badion)

os.makedirs(OUT, exist_ok=True)
json.dump(records, open(os.path.join(OUT, "targets.json"), "w"),
          indent=1, ensure_ascii=False)
print("wrote", os.path.join(OUT, "targets.json"))
