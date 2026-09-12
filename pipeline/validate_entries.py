"""Stage 4: validate every candidate against AlphaFold Server's real limits and
emit a per-category shortlist.

Server constraints encoded here (all from alphafoldserver.com/faq):
  * 5,000 token limit: 1/residue, 1/nucleotide, 1/ligand atom, 1/ion
  * protein input: standard 20 letters only (B J O U X Z unsupported)
  * each polymer chain >= 4 residues/nucleotides
  * ions limited to the built-in 10 (CCD: CA CO CU FE K MG MN NA ZN CL)
  * ligands: any CCD code, but dictionary is frozen at version 2024_10_28
  * no covalent-link specification, no waters, no hydrogens
"""
import json
import re

CCD_FREEZE = "2024-10-28"
DEPOSIT_AFTER = "2025-02-03"
TOKEN_LIMIT = 5000

meta = json.load(open("meta_all.json"))
ccd = json.load(open("ccd.json"))
pools = json.load(open("pools.json"))

# ions the Server exposes natively
SUPPORTED_IONS = {"CA", "CO", "CU", "FE", "K", "MG", "MN", "NA", "ZN", "CL"}

# crystallisation additives / cryoprotectants: a student would simply not enter
# these, so their presence must not disqualify an entry
IGNORABLE = {
    "HOH", "DOD", "GOL", "EDO", "PEG", "PGE", "PG4", "PG0", "1PE", "2PE", "P6G",
    "7PE", "12P", "15P", "XPE", "MPD", "MRD", "SO4", "PO4", "ACT", "ACY", "FMT",
    "DMS", "TRS", "IMD", "EPE", "MES", "CIT", "FLC", "TAR", "MLI", "MLA", "SIN",
    "BME", "DTT", "DTU", "DTV", "TCE", "BTB", "BU3", "HEZ", "PGO", "PDO", "IPA",
    "MOH", "EOH", "ETX", "NO3", "AZI", "SCN", "CO3", "BCT", "CAC", "CAD", "NH4",
    "OXL", "GLC", "SUC", "TLA", "MAN", "BOG", "LDA", "C8E", "OGA", "UNX", "UNL",
    "PGR", "PIN", "POL", "SPD", "SPM", "GAI", "ARS", "PYR", "FUM", "MAE",
}

# non-supported ions/heavies frequently present incidentally; if the entry is
# otherwise good, a student just omits them (noted in the sheet)
OMITTABLE_HEAVY = {
    "NI", "CD", "HG", "PT", "AU", "AG", "PB", "SR", "BA", "CS", "RB", "LI",
    "BR", "IOD", "F", "AL", "V", "W", "MO", "SE", "TL", "GA", "IN", "SM", "EU",
    "GD", "TB", "YB", "LU", "PR", "ND", "HO", "ER", "ZR", "TA", "RE", "OS",
    "IR", "RU", "RH", "PD", "SB", "SN", "TE", "BI", "U1", "Y1", "3NI", "2HP",
}

STD_AA = set("ACDEFGHIKLMNPQRSTVWY")


def ccd_ok(code):
    """Is this CCD code present in the frozen 2024_10_28 dictionary?"""
    rec = ccd.get(code)
    if not rec:
        return False
    d = (rec.get("rcsb_chem_comp_info") or {}).get("initial_release_date")
    if not d:
        return False
    return d[:10] <= CCD_FREEZE


def ccd_atoms(code):
    rec = ccd.get(code) or {}
    info = rec.get("rcsb_chem_comp_info") or {}
    return info.get("atom_count_heavy") or info.get("atom_count") or 0


def ccd_name(code):
    rec = ccd.get(code) or {}
    return ((rec.get("chem_comp") or {}).get("name") or "").strip()


def ccd_weight(code):
    rec = ccd.get(code) or {}
    return (rec.get("chem_comp") or {}).get("formula_weight") or 0


def classify(entry):
    """Return (verdict, info dict) for one PDB entry."""
    eid = entry["rcsb_id"]
    info = {"id": eid}
    ai = entry.get("rcsb_accession_info") or {}
    ei = entry.get("rcsb_entry_info") or {}

    dep = (ai.get("deposit_date") or "")[:10]
    rel = (ai.get("initial_release_date") or "")[:10]
    info["deposit"] = dep
    info["release"] = rel
    if dep <= DEPOSIT_AFTER:
        return "reject:deposit", info
    if not rel:
        return "reject:unreleased", info

    methods = [m.get("method") for m in (entry.get("exptl") or [])]
    info["method"] = "/".join([m for m in methods if m]) or "?"
    resl = ei.get("resolution_combined") or []
    info["resolution"] = resl[0] if resl else None
    info["title"] = ((entry.get("struct") or {}).get("title") or "").strip()
    info["keywords"] = ((entry.get("struct_keywords") or {}).get("pdbx_keywords") or "").strip()

    # ---- polymers ----
    chains = []
    tokens = 0
    for pe in entry.get("polymer_entities") or []:
        ep = pe.get("entity_poly") or {}
        seq = (ep.get("pdbx_seq_one_letter_code_can") or "").replace("\n", "").strip().upper()
        ptype = ep.get("rcsb_entity_polymer_type") or "?"
        auth = ((pe.get("rcsb_polymer_entity_container_identifiers") or {}).get("auth_asym_ids") or [])
        desc = ((pe.get("rcsb_polymer_entity") or {}).get("pdbx_description") or "").strip()
        orgs = [o.get("ncbi_scientific_name") for o in (pe.get("rcsb_entity_source_organism") or []) if o.get("ncbi_scientific_name")]
        n_copies = len(auth) or 1
        if not seq or len(seq) < 4:
            return "reject:shortchain", info
        if ptype == "Protein":
            bad = set(seq) - STD_AA
            if bad:
                return "reject:nonstdaa:" + "".join(sorted(bad)), info
        elif ptype == "DNA":
            if set(seq) - set("ACGT"):
                return "reject:nonstddna", info
        elif ptype == "RNA":
            if set(seq) - set("ACGU"):
                return "reject:nonstdrna", info
        else:
            return "reject:polytype:" + ptype, info
        tokens += len(seq) * n_copies
        chains.append({"type": ptype, "len": len(seq), "copies": n_copies,
                       "desc": desc, "org": orgs[:1], "seq": seq})
    if not chains:
        return "reject:nopolymer", info
    info["chains"] = chains
    info["n_prot_entities"] = ei.get("polymer_entity_count_protein") or 0
    info["n_dna"] = ei.get("polymer_entity_count_DNA") or 0
    info["n_rna"] = ei.get("polymer_entity_count_RNA") or 0
    info["polymer_tokens"] = tokens

    # ---- ligands ----
    keep_lig, keep_ion, omit, blocked = [], [], [], []
    for ne in entry.get("nonpolymer_entities") or []:
        comp = (ne.get("nonpolymer_comp") or {}).get("chem_comp") or {}
        code = comp.get("id")
        if not code:
            continue
        n = len(((ne.get("rcsb_nonpolymer_entity_container_identifiers") or {}).get("auth_asym_ids") or [])) or 1
        if code in IGNORABLE:
            omit.append((code, n))
            continue
        if code in SUPPORTED_IONS:
            keep_ion.append((code, n))
            tokens += n
            continue
        if code in OMITTABLE_HEAVY:
            omit.append((code, n))
            continue
        # a genuine ligand: must be in the frozen dictionary
        if not ccd_ok(code):
            blocked.append((code, n))
            continue
        na = ccd_atoms(code)
        keep_lig.append((code, n, na, ccd_name(code)))
        tokens += na * n

    info["ligands"] = keep_lig
    info["ions"] = keep_ion
    info["omitted"] = omit
    info["blocked"] = blocked
    info["tokens"] = tokens

    if blocked:
        return "reject:ccd_frozen", info
    if tokens > TOKEN_LIMIT:
        return "reject:tokens", info
    return "ok", info


results, rejects = {}, {}
for eid, e in meta.items():
    verdict, info = classify(e)
    if verdict == "ok":
        results[eid] = info
    else:
        rejects[eid] = (verdict, info)

print(f"validated OK: {len(results)} / {len(meta)}")
from collections import Counter
print(Counter(v.split(':', 1)[1] for v, _ in rejects.values()).most_common())

json.dump(results, open("valid.json", "w"), indent=1)
json.dump({k: v[0] for k, v in rejects.items()}, open("rejects.json", "w"), indent=1)

# report per pool category how many survived
print()
for cat, ids in pools.items():
    kept = [i for i in ids if i in results]
    print(f"{cat:18s} {len(kept):3d}/{len(ids)}")
