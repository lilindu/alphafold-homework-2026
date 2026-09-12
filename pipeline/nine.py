"""Census phase 3: classify the gate-passing structures into the nine classes
and report counts.

Class definitions (as agreed):
  1 hetero  异源蛋白复合体 -- >=2 distinct protein entities, non-antibody,
            non-nucleic, non-membrane
  2 metal   蛋白 + 金属离子 (built-in ions only)
  3 dna     蛋白 + DNA
  4 rna     蛋白 + RNA
  5 cofac   蛋白 + 内置辅因子 (one of the 19 built-in ligand codes)
  6 ligand  蛋白 + 非内置小分子配体 (in frozen CCD; entered via CCD code)
  7 ptm     翻译后修饰 (modification in the Server's ptmType list)
  8 memb    膜蛋白 (PDBTM / MemProtMD / mpstruc annotation)
  9 ab      抗原抗体复合物 (>=2 antibody chains + >=1 antigen chain)

Counts are reported two ways, because they answer different questions:
  eligible  -- how many gate-passing structures satisfy that class at all
                (this is the pool depth available for picking)
  exclusive -- how many satisfy exactly one class (cleanest teaching examples)
"""
import json
import re
from collections import Counter, defaultdict

valid = json.load(open("census_valid.json"))

BUILTIN_LIG = {"ADP", "ATP", "AMP", "GTP", "GDP", "FAD", "NAD", "NAP", "NDP",
               "HEM", "HEC", "PLM", "OLA", "MYR", "CIT", "CLA", "CHL", "BCL", "BCB"}
TEACHING_METALS = {"ZN", "CU", "FE", "MN", "CO", "MG", "CA"}
STRUCTURAL_METALS = {"ZN", "CU", "FE", "MN", "CO"}
MEMB_ANNOT = {"PDBTM", "MemProtMD", "mpstruc"}

AB_RE = re.compile(r"heavy chain|light chain|\bfab\b|\bigg\b|antibody|immunoglobulin"
                   r"|\bscfv\b|\bnanobody\b|\bvhh\b|single[- ]domain antibody|sybody"
                   r"|\bmab\b|\bfv fragment\b", re.I)
NB_RE = re.compile(r"nanobody|\bvhh\b|single[- ]domain antibody|sybody|camelid", re.I)
TCR_RE = re.compile(r"\bT[- ]cell receptor\b|\bTCR\b|beta-2 microglobulin|\bMHC\b"
                    r"|\bHLA\b|histocompatibility|\bH-2[A-Z]", re.I)


def prot_chains(r):
    return [c for c in r["chains"] if c["type"] == "Protein"]


def nuc_chains(r):
    return [c for c in r["chains"] if c["type"] in ("DNA", "RNA")]


def has_dna(r):
    return any(c["type"] == "DNA" for c in r["chains"])


def has_rna(r):
    return any(c["type"] == "RNA" for c in r["chains"])


def blob(r):
    return r["title"] + " || " + " | ".join(c["desc"] for c in r["chains"])


def is_membrane(r):
    return any(set(c["annot"]) & MEMB_ANNOT for c in r["chains"])


def ab_chains(r):
    return [c for c in prot_chains(r)
            if AB_RE.search(c["desc"]) and not NB_RE.search(c["desc"])]


def nb_chains(r):
    return [c for c in prot_chains(r) if NB_RE.search(c["desc"])]


def antigen_chains(r):
    return [c for c in prot_chains(r) if not AB_RE.search(c["desc"])]


def is_antibody_complex(r):
    ab, nb, ag = ab_chains(r), nb_chains(r), antigen_chains(r)
    if TCR_RE.search(blob(r)):
        return False
    if len(ab) >= 2 and len(ag) >= 1:
        return True
    if len(nb) >= 1 and len(ag) >= 1:
        return True
    return False


def builtin_cofactors(r):
    return [l for l in r["ligands"] if l["code"] in BUILTIN_LIG]


def other_ligands(r):
    """Non-built-in, drug-like: needs a manual CCD code. Exclude tiny fragments
    so this class means a real small-molecule complex."""
    return [l for l in r["ligands"]
            if l["code"] not in BUILTIN_LIG and (l["mw"] or 0) >= 150
            and (l["atoms"] or 0) >= 10]


def teaching_metals(r):
    return [(c, n) for c, n in r["ions"] if c in TEACHING_METALS]


def structural_metals(r):
    return [(c, n) for c, n in r["ions"] if c in STRUCTURAL_METALS]


def classes_of(r):
    out = set()
    nprot_ent = len(prot_chains(r))
    memb = is_membrane(r)
    ab = is_antibody_complex(r)

    if ab:
        out.add("ab")
    if memb:
        out.add("memb")
    if has_dna(r):
        out.add("dna")
    if has_rna(r):
        out.add("rna")
    if r["ptms"]:
        out.add("ptm")
    if builtin_cofactors(r):
        out.add("cofac")
    if other_ligands(r):
        out.add("ligand")
    if teaching_metals(r):
        out.add("metal")
    # class 1 is deliberately narrow
    if (nprot_ent >= 2 and not nuc_chains(r) and not ab and not memb
            and not TCR_RE.search(blob(r))):
        out.add("hetero")
    return out


tagged = {}
for eid, r in valid.items():
    cs = classes_of(r)
    if cs:
        tagged[eid] = cs

elig = Counter()
for cs in tagged.values():
    for c in cs:
        elig[c] += 1

excl = Counter()
for cs in tagged.values():
    if len(cs) == 1:
        excl[next(iter(cs))] += 1

LABEL = {
    "hetero": "1 异源蛋白复合体(非抗体/非核酸/非膜)",
    "metal":  "2 蛋白 + 金属离子",
    "dna":    "3 蛋白 + DNA",
    "rna":    "4 蛋白 + RNA",
    "cofac":  "5 蛋白 + 内置辅因子",
    "ligand": "6 蛋白 + 非内置小分子配体",
    "ptm":    "7 有翻译后修饰的蛋白",
    "memb":   "8 膜蛋白",
    "ab":     "9 抗原抗体复合物",
}
ORDER = ["hetero", "metal", "dna", "rna", "cofac", "ligand", "ptm", "memb", "ab"]

print(f"gate-passing structures: {len(valid)}")
print(f"of which fall into >=1 of the nine classes: {len(tagged)}\n")
print(f"{'类别':42s} {'可用池':>7s} {'仅属此类':>9s}   需要")
need = {c: 4 for c in ORDER}
need["ab"] = 9
for c in ORDER:
    flag = "" if elig[c] >= need[c] else "   <-- 不足!"
    print(f"{LABEL[c]:42s} {elig[c]:7d} {excl[c]:9d}   {need[c]}{flag}")

print("\n--- class 7 (PTM) detail: which modifications are available ---")
ptmc = Counter()
for eid, r in valid.items():
    for p in r["ptms"]:
        ptmc[p["code"]] += 1
for code, n in ptmc.most_common():
    print(f"   {code}: {n} entries")

print("\n--- class 2 (metal) detail ---")
mc = Counter()
for eid, r in valid.items():
    for c, n in teaching_metals(r):
        mc[c] += 1
print("   ", dict(mc.most_common()))
sm = sum(1 for r in valid.values() if structural_metals(r))
print(f"    with a structural/catalytic transition metal (Zn/Cu/Fe/Mn/Co): {sm}")

print("\n--- class 5 (built-in cofactor) detail ---")
cc = Counter()
for r in valid.values():
    for l in builtin_cofactors(r):
        cc[l["code"]] += 1
print("   ", dict(cc.most_common()))

print("\n--- class 9 (antibody) detail ---")
fab = sum(1 for r in valid.values() if is_antibody_complex(r) and len(ab_chains(r)) >= 2)
nb = sum(1 for r in valid.values() if is_antibody_complex(r) and nb_chains(r))
print(f"    Fab/IgG-type: {fab}    nanobody/VHH-type: {nb}")

json.dump({k: sorted(v) for k, v in tagged.items()}, open("census_classes.json", "w"))
print("\nwrote census_classes.json")
