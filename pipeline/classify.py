"""Census phase 2: apply Server gates, then classify into the nine course classes
and count each class.

A structure may satisfy several classes; the census reports both the raw
per-class counts (how many candidates exist for that class) and a primary
assignment used for reporting, so the totals are interpretable.
"""
import json
import re
from collections import Counter, defaultdict

CCD_FREEZE = "2024-10-28"
DEP_AFTER = "2025-02-03"
TOKEN_LIMIT = 5000
MAX_ASSEMBLY = 12

BUILTIN_LIG = {"ADP", "ATP", "AMP", "GTP", "GDP", "FAD", "NAD", "NAP", "NDP",
               "HEM", "HEC", "PLM", "OLA", "MYR", "CIT", "CLA", "CHL", "BCL", "BCB"}
SUPPORTED_IONS = {"CA", "CO", "CU", "FE", "K", "MG", "MN", "NA", "ZN", "CL"}
# ions that count as a genuine "protein + metal" teaching case
TEACHING_METALS = {"ZN", "CU", "FE", "MN", "CO", "MG", "CA"}
STRUCTURAL_METALS = {"ZN", "CU", "FE", "MN", "CO"}

SERVER_PTM = {"SEP", "TPO", "PTR", "NEP", "HIP", "ALY", "MLY", "M3L", "MLZ",
              "2MR", "AGM", "MCS", "HYP", "HY3", "LYZ", "AHB", "P1L", "SNN",
              "SNC", "TRF", "KCR", "CIR", "YHA"}
CRYST_SUBSTITUTE = {"MSE", "MLE"}
MEMB_ANNOT = {"PDBTM", "MemProtMD", "mpstruc"}

SUGARS = {"NAG", "BMA", "MAN", "GLC", "BGC", "FUC", "GAL", "NDG", "A2G", "SIA",
          "GLA", "XYS", "RAM", "NGA", "BGL", "GCU", "ADA", "MBG"}

IGNORABLE = {
    "HOH", "DOD", "GOL", "EDO", "PEG", "PGE", "PG4", "PG0", "1PE", "2PE", "P6G",
    "7PE", "12P", "15P", "XPE", "MPD", "MRD", "SO4", "PO4", "ACT", "ACY", "FMT",
    "DMS", "TRS", "IMD", "EPE", "MES", "FLC", "TAR", "MLI", "MLA", "SIN",
    "BME", "DTT", "DTU", "DTV", "TCE", "BTB", "BU3", "HEZ", "PGO", "PDO", "IPA",
    "MOH", "EOH", "ETX", "NO3", "AZI", "SCN", "CO3", "BCT", "CAC", "NH4",
    "OXL", "SUC", "TLA", "BOG", "LDA", "C8E", "OGA", "UNX", "UNL", "PGR", "PIN",
    "POL", "SPD", "SPM", "GAI", "PYR", "FUM", "MAE", "MPO", "NHE", "CXS", "MRO",
    "SGM", "OCT", "D10", "DD9", "UND", "DPO", "TZZ", "ACE", "NH2", "PE4", "PE8",
}
OMITTABLE_HEAVY = {
    "NI", "CD", "HG", "PT", "AU", "AG", "PB", "SR", "BA", "CS", "RB", "LI",
    "BR", "IOD", "F", "AL", "V", "W", "MO", "SE", "TL", "GA", "IN", "SM", "EU",
    "GD", "TB", "YB", "LU", "PR", "ND", "HO", "ER", "ZR", "TA", "RE", "OS",
    "IR", "RU", "RH", "PD", "SB", "SN", "TE", "BI", "3NI", "2HP", "YT3", "CE",
    "XE", "KR", "AR", "HO3", "4MO", "6MO",
}
STD_AA = set("ACDEFGHIKLMNPQRSTVWY")

AB_RE = re.compile(r"heavy chain|light chain|\bfab\b|\bigg\b|antibody|immunoglobulin"
                   r"|\bscfv\b|\bnanobody\b|\bvhh\b|single[- ]domain antibody|sybody"
                   r"|\bmab\b|\bfv fragment\b", re.I)
NB_RE = re.compile(r"nanobody|\bvhh\b|single[- ]domain antibody|sybody|camelid", re.I)

meta = json.load(open("cmeta.json"))
ccd = json.load(open("ccd.json"))
asm = json.load(open("casm.json"))


def ccd_ok(code):
    r = ccd.get(code)
    if not r:
        return False
    d = (r.get("rcsb_chem_comp_info") or {}).get("initial_release_date")
    return bool(d) and d[:10] <= CCD_FREEZE


def ccd_atoms(code):
    i = (ccd.get(code) or {}).get("rcsb_chem_comp_info") or {}
    return i.get("atom_count_heavy") or i.get("atom_count") or 0


def ccd_name(code):
    return (((ccd.get(code) or {}).get("chem_comp") or {}).get("name") or "").strip()


def ccd_mw(code):
    return (((ccd.get(code) or {}).get("chem_comp") or {}).get("formula_weight") or 0)


def parse_mods(raw):
    if not raw:
        return []
    s = raw.replace("\n", "").strip()
    out, pos, i = [], 0, 0
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


def bio(eid):
    a = asm.get(eid) or []
    best = None
    for x in a:
        ai = x.get("rcsb_assembly_info") or {}
        n = ai.get("polymer_entity_instance_count") or 0
        if best is None or n > best["inst"]:
            best = {"inst": n, "ents": ai.get("polymer_entity_count") or 0,
                    "prot": ai.get("polymer_entity_instance_count_protein") or 0,
                    "dna": ai.get("polymer_entity_instance_count_DNA") or 0,
                    "rna": ai.get("polymer_entity_instance_count_RNA") or 0}
    return best


def evaluate(eid, e):
    """Apply the Server gates. Return dict or (None, reason)."""
    ai = e.get("rcsb_accession_info") or {}
    dep = (ai.get("deposit_date") or "")[:10]
    rel = (ai.get("initial_release_date") or "")[:10]
    if not dep or dep <= DEP_AFTER:
        return None, "deposit"
    if not rel:
        return None, "unreleased"

    ei = e.get("rcsb_entry_info") or {}
    b = bio(eid)
    if not b or b["inst"] <= 0:
        return None, "no_assembly"
    if b["inst"] > MAX_ASSEMBLY:
        return None, "assembly_too_big"

    chains, tokens = [], 0
    ptms, bad_ptm = [], []
    for pe in e.get("polymer_entities") or []:
        ep = pe.get("entity_poly") or {}
        seq = (ep.get("pdbx_seq_one_letter_code_can") or "").replace("\n", "").strip().upper()
        ptype = ep.get("rcsb_entity_polymer_type") or "?"
        if not seq or len(seq) < 4:
            return None, "shortchain"
        if ptype == "Protein":
            if set(seq) - STD_AA:
                return None, "nonstd_aa"
        elif ptype == "DNA":
            if set(seq) - set("ACGT"):
                return None, "nonstd_dna"
        elif ptype == "RNA":
            if set(seq) - set("ACGU"):
                return None, "nonstd_rna"
        else:
            return None, "polytype"
        ent = pe["rcsb_id"].split("_")[-1]
        auth = ((pe.get("rcsb_polymer_entity_container_identifiers") or {}).get("auth_asym_ids") or [])
        desc = ((pe.get("rcsb_polymer_entity") or {}).get("pdbx_description") or "").strip()
        orgs = [o.get("ncbi_scientific_name") for o in (pe.get("rcsb_entity_source_organism") or [])
                if o.get("ncbi_scientific_name")]
        ann = {a.get("type") for a in (pe.get("rcsb_polymer_entity_annotation") or []) if a.get("type")}
        ups = {u.get("rcsb_id") for u in (pe.get("uniprots") or []) if u.get("rcsb_id")}
        n = len(auth) or 1
        tokens += len(seq) * n
        # modified residues, from the non-canonical sequence
        for pos, code in parse_mods(ep.get("pdbx_seq_one_letter_code")):
            if code in CRYST_SUBSTITUTE:
                continue
            if ptype == "Protein" and len(code) >= 3 and code not in ("DA", "DT", "DG", "DC"):
                (ptms if code in SERVER_PTM else bad_ptm).append(
                    {"entity": ent, "pos": pos, "code": code})
        chains.append({"entity": ent, "type": ptype, "seq": seq, "len": len(seq),
                       "copies": n, "desc": desc, "org": orgs[0] if orgs else "",
                       "annot": sorted(ann), "uniprots": sorted(ups)})
    if not chains:
        return None, "nopolymer"

    ligs, ions, omit, blocked, sugars = [], [], [], [], []
    for ne in e.get("nonpolymer_entities") or []:
        cc = (ne.get("nonpolymer_comp") or {}).get("chem_comp") or {}
        code = cc.get("id")
        if not code:
            continue
        n = len(((ne.get("rcsb_nonpolymer_entity_container_identifiers") or {}).get("auth_asym_ids") or [])) or 1
        if code in IGNORABLE:
            omit.append((code, n)); continue
        if code in SUPPORTED_IONS:
            ions.append((code, n)); tokens += n; continue
        if code in OMITTABLE_HEAVY:
            omit.append((code, n)); continue
        if code in SUGARS:
            sugars.append((code, n)); continue
        if not ccd_ok(code):
            blocked.append((code, n)); continue
        a = ccd_atoms(code)
        ligs.append({"code": code, "count": n, "atoms": a, "name": ccd_name(code),
                     "mw": ccd_mw(code)})
        tokens += a * n
    if blocked:
        return None, "ccd_frozen"

    branched = []
    for be in e.get("branched_entities") or []:
        ci = be.get("rcsb_branched_entity_container_identifiers") or {}
        branched.append({"monomers": ci.get("chem_comp_monomers") or [],
                         "desc": ((be.get("rcsb_branched_entity") or {}).get("pdbx_description") or "")})

    if tokens > TOKEN_LIMIT:
        return None, "tokens"

    resl = ei.get("resolution_combined") or []
    return {
        "id": eid,
        "title": ((e.get("struct") or {}).get("title") or "").strip(),
        "deposit": dep, "release": rel,
        "method": "/".join(m.get("method") for m in (e.get("exptl") or []) if m.get("method")),
        "resolution": resl[0] if resl else None,
        "keywords": ((e.get("struct_keywords") or {}).get("pdbx_keywords") or "").strip(),
        "chains": chains, "ligands": ligs, "ions": ions, "omitted": omit,
        "sugars": sugars, "branched": branched,
        "ptms": ptms, "unsupported_ptms": bad_ptm,
        "tokens": tokens, "bio": b,
    }, None


ok, rej = {}, Counter()
for eid, e in meta.items():
    r, why = evaluate(eid, e)
    if r:
        ok[eid] = r
    else:
        rej[why] += 1

print(f"passed Server gates: {len(ok)} / {len(meta)}")
print("rejections:", rej.most_common(), "\n")
json.dump(ok, open("census_valid.json", "w"))
