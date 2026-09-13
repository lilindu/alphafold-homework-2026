"""Final selection: 41 targets over SEVEN classes (metal and cofactor dropped,
their slots given to membrane and antibody).

Per 班: hetero 1 + DNA 1 + RNA 1 + ligand 1 + PTM 1 + membrane 2 + antibody 3
        = 10, and 班4 gets one extra antibody = 11.
Totals: hetero 4, DNA 4, RNA 4, ligand 4, PTM 4, membrane 8, antibody 13 = 41.

Homology screening is folded in where it was run (PTM, membrane, hetero, ligand):
those classes prefer candidates whose chains have the FEWEST pre-2021 homologs,
so the prediction is a real test rather than template recall. DNA/RNA/antibody
were not screened -- for DNA/RNA the nucleic partner has no training-homolog
notion, and antibody frameworks are unavoidably >=80% (measured), so screening
them would not change any choice.

PTM keeps its interface requirement: the modified residue must make a
non-covalent, multi-point contact across a protein-protein interface
(2.2 A <= min dist < 3.0 A, >= 20 atoms within 4 A). Among those, pick the
lowest homology -- ranked by homolog COUNT, since 7 homologs at >=95% is a much
weaker leak than 278.
"""
import json
import re

valid = json.load(open("census_valid.json"))
iface = json.load(open("iface.json"))
screen = json.load(open("screen.json"))

MEMB_ANNOT = {"PDBTM", "MemProtMD", "mpstruc"}
BUILTIN_LIG = {"ADP", "ATP", "AMP", "GTP", "GDP", "FAD", "NAD", "NAP", "NDP",
               "HEM", "HEC", "PLM", "OLA", "MYR", "CIT", "CLA", "CHL", "BCL", "BCB"}
AB_RE = re.compile(r"heavy chain|light chain|\bfab\b|\bigg\b|antibody|immunoglobulin"
                   r"|\bscfv\b|\bnanobody\b|\bvhh\b|single[- ]domain antibody|sybody"
                   r"|\bmab\b|\bfv fragment\b|\bvariable domain\b", re.I)
NB_RE = re.compile(r"nanobody|\bvhh\b|single[- ]domain antibody|sybody|camelid", re.I)
TCR_RE = re.compile(r"\bT[- ]cell receptor\b|\bTCR\b|beta-2 microglobulin|\bMHC\b"
                    r"|\bHLA\b|histocompatibility|\bH-2[A-Z]", re.I)


def P(r):
    return [c for c in r["chains"] if c["type"] == "Protein"]


def NUC(r):
    return [c for c in r["chains"] if c["type"] in ("DNA", "RNA")]


def blob(r):
    return r["title"] + " || " + " | ".join(c["desc"] for c in r["chains"])


def memb(r):
    return any(set(c["annot"]) & MEMB_ANNOT for c in r["chains"])


def has_ab_chain(r):
    if AB_RE.search(r["title"]) or NB_RE.search(r["title"]):
        return True
    return any(AB_RE.search(c["desc"]) or NB_RE.search(c["desc"]) for c in r["chains"])


def abch(r):
    return [c for c in P(r) if AB_RE.search(c["desc"]) and not NB_RE.search(c["desc"])]


def nbch(r):
    return [c for c in P(r) if NB_RE.search(c["desc"])]


def agch(r):
    return [c for c in P(r)
            if not AB_RE.search(c["desc"]) and not NB_RE.search(c["desc"])]


def is_ab_complex(r):
    if TCR_RE.search(blob(r)) or not agch(r):
        return False
    return len(abch(r)) >= 2 or len(nbch(r)) >= 1


def cofac(r):
    return [l for l in r["ligands"] if l["code"] in BUILTIN_LIG]


def otherlig(r):
    return [l for l in r["ligands"] if l["code"] not in BUILTIN_LIG
            and (l["mw"] or 0) >= 150 and (l["atoms"] or 0) >= 10]


def clean(r):
    return not r["branched"] and not r["sugars"]


def resol(r):
    return r["resolution"] or 9.9


def uniprots(r):
    s = set()
    for c in r["chains"]:
        s |= set(c["uniprots"])
    return s


def homology(r):
    """(worst_bin, worst_count) over screened protein chains; None if unscreened."""
    best_bin, best_n = None, None
    for c in P(r):
        v = screen.get(c["seq"])
        if not v or v.get("bin") is None:
            continue
        b, n = v["bin"], v["n"]
        if best_bin is None or b > best_bin or (b == best_bin and n > best_n):
            best_bin, best_n = b, n
    return best_bin, best_n


def hom_rank(r):
    b, n = homology(r)
    if b is None:
        return (9, 9)                 # unscreened sorts last
    return (b, n or 0)


def ptm_iface(pid):
    d = iface.get(pid)
    if not d or not d.get("ptms"):
        return None
    best = None
    for p in d["ptms"]:
        if p["min_dist"] is None:
            continue
        if best is None or p["n_close"] > best[1]:
            best = (p["min_dist"], p["n_close"], p["code"])
    return best


def ptm_ok(r):
    it = ptm_iface(r["id"])
    if not it:
        return False
    mind, ncl, _ = it
    return 2.2 <= mind < 3.0 and ncl >= 20 and len({c["seq"] for c in P(r)}) >= 2


STOP = {"crystal", "structure", "cryo", "complex", "with", "bound", "from", "the",
        "and", "form", "state", "resolution", "angstrom", "human", "mutant",
        "variant", "wild", "type", "apo", "holo", "domain", "protein", "analysis",
        "structural", "class", "conformation"}


def topic(r):
    t = re.sub(r"[^a-z0-9 ]", " ", r["title"].lower())
    return set(w for w in t.split() if len(w) > 3 and w not in STOP)


def longest(r):
    return max(r["chains"], key=lambda c: c["len"])


def similar(a, b):
    ua, ub = uniprots(a), uniprots(b)
    if ua and ub and (ua & ub):
        return True
    la, lb = longest(a), longest(b)
    if abs(la["len"] - lb["len"]) <= max(3, 0.05 * max(la["len"], lb["len"])):
        n = min(50, la["len"], lb["len"])
        if n >= 25 and la["seq"][:n] == lb["seq"][:n]:
            return True
    return len(topic(a) & topic(b)) >= 2


CLASSES = [
    dict(key="hetero", n=4, label="异源蛋白复合体",
         note="两条以上不同蛋白链的复合体(不含抗体链、不含核酸、非膜蛋白)。"
              "考察跨链共进化信号能否把界面摆对",
         test=lambda r: (len(P(r)) >= 2 and not NUC(r) and not has_ab_chain(r)
                         and not memb(r) and not TCR_RE.search(blob(r))
                         and not cofac(r) and not otherlig(r)
                         and all(c["len"] >= 50 for c in P(r))
                         and len(uniprots(r)) >= 2 and clean(r) and resol(r) <= 2.6),
         order="homology"),
    dict(key="dna", n=4, label="蛋白 + DNA",
         note="只含标准 A/C/G/T;双链须分别输入两条互补链",
         test=lambda r: (any(c["type"] == "DNA" for c in r["chains"])
                         and not any(c["type"] == "RNA" for c in r["chains"])
                         and P(r) and not has_ab_chain(r)
                         and sum(c["len"] * c["copies"] for c in r["chains"]
                                 if c["type"] == "DNA") >= 10
                         and clean(r) and resol(r) <= 2.8),
         order="diverse"),
    dict(key="rna", n=4, label="蛋白 + RNA",
         note="RNA 构象自由度大,是公认较难的一类",
         test=lambda r: (any(c["type"] == "RNA" for c in r["chains"])
                         and P(r) and not has_ab_chain(r)
                         and sum(c["len"] * c["copies"] for c in r["chains"]
                                 if c["type"] == "RNA") >= 8
                         and clean(r) and resol(r) <= 3.2),
         order="diverse"),
    dict(key="ligand", n=4, label="蛋白 + 小分子配体",
         note="配体不在 Server 的 19 种内置辅因子里,走任意 CCD 代码那条路;"
              "代码已核实在冻结的 CCD 2024_10_28 字典中",
         test=lambda r: (otherlig(r) and not cofac(r) and not NUC(r)
                         and not has_ab_chain(r) and not memb(r)
                         and clean(r) and resol(r) <= 2.2),
         order="homology"),
    dict(key="ptm", n=4, label="有翻译后修饰的蛋白",
         note="修饰残基经坐标实测介导蛋白-蛋白互作(最近距离 2.2–3.0 Å,"
              "4 Å 内接触原子 ≥ 20)。该类受体多为反复研究的识别模块,"
              "同源体无法避免,已在同源体数量上取最低",
         test=lambda r: (r["ptms"] and len(P(r)) >= 2 and not has_ab_chain(r)
                         and clean(r) and ptm_ok(r)),
         order="homology"),
    dict(key="memb", n=8, label="膜蛋白",
         note="带 PDBTM / MemProtMD / mpstruc 跨膜注释。Server 不建模膜平面,"
              "跨膜螺旋排布与构象态最易错;本类优先选训练窗口内查不到同源体的",
         test=lambda r: (memb(r) and P(r) and not has_ab_chain(r) and clean(r)
                         and 250 <= sum(c["len"] * c["copies"] for c in r["chains"]) <= 2400),
         order="homology"),
]

AB = dict(key="ab", label="抗原抗体复合物",
          note="抗体 CDR 环构象 + 表位定位。抗体骨架本身高度保守(实测 ≥80%),"
               "无法规避,真正未知的是表位识别;建议跑多个 seed 按 ipTM 排序")


def pick(pool, n, order):
    cands = list(pool)
    if order == "homology":
        cands.sort(key=lambda r: (hom_rank(r), resol(r)))
    else:
        cands.sort(key=lambda r: r["tokens"])
        k = len(cands)
        bands = [cands[(k * b) // n:(k * (b + 1)) // n] or cands for b in range(n)]
        out = []
        for band in bands:
            for r in sorted(band, key=resol):
                if not any(x["id"] == r["id"] for x in out) and \
                        not any(similar(r, x) for x in out):
                    out.append(r)
                    break
        cands = out + [c for c in cands if c not in out]
    chosen = []
    for r in cands:
        if len(chosen) >= n:
            break
        if any(x["id"] == r["id"] for x in chosen):
            continue
        if any(similar(r, x) for x in chosen):
            continue
        chosen.append(r)
    return chosen[:n]


used, result = set(), {}
for spec in CLASSES:
    pool = [r for r in valid.values() if r["id"] not in used and spec["test"](r)]
    got = pick(pool, spec["n"], spec["order"])
    for r in got:
        used.add(r["id"])
    result[spec["key"]] = got
    print(f"{spec['key']:8s} pool={len(pool):5d} picked={len(got)}/{spec['n']}")
    for r in got:
        b, nn = homology(r)
        lab = ("未查" if b is None else "无" if b == 0 else
               f">={int(b*100)}%")
        print(f"     {r['id']}  同源{lab:>6s}({nn if nn is not None else '-'})"
              f"  res={str(resol(r))[:4]:4s} tok={r['tokens']:5d}  {r['title'][:52]}")

# ---- antibodies: 13, split between peptide and folded-domain epitopes ----
def epitope(r):
    ag = agch(r)
    return "peptide" if ag and max(c["len"] for c in ag) <= 30 else "domain"


abpool = [r for r in valid.values() if r["id"] not in used and is_ab_complex(r)
          and clean(r) and resol(r) <= 3.0
          and 300 <= r["tokens"] <= 1600]
pep = [r for r in abpool if epitope(r) == "peptide"]
dom = [r for r in abpool if epitope(r) == "domain"]
nb = [r for r in abpool if nbch(r)]
sel_nb = pick(nb, 3, "diverse")
for r in sel_nb:
    used.add(r["id"])
pep = [r for r in pep if r["id"] not in used]
dom = [r for r in dom if r["id"] not in used]
sel_pep = pick(pep, 5, "diverse")
for r in sel_pep:
    used.add(r["id"])
dom = [r for r in dom if r["id"] not in used]
sel_dom = pick(dom, 5, "diverse")
for r in sel_dom:
    used.add(r["id"])
abs_ = sel_nb + sel_pep + sel_dom
print(f"ab       pool={len(abpool):5d} picked={len(abs_)}/13"
      f"  (纳米抗体 {len(sel_nb)}, 线性肽表位 {len(sel_pep)}, 折叠域表位 {len(sel_dom)})")
for r in abs_:
    print(f"     {r['id']}  {epitope(r):8s} res={str(resol(r))[:4]:4s} tok={r['tokens']:5d}"
          f"  {r['title'][:52]}")
result["ab"] = abs_

# ---- allocate to the four 班 ----
BAN = ["班1", "班2", "班3", "班4"]
alloc = {b: [] for b in BAN}
SPEC = {s["key"]: s for s in CLASSES}
SPEC["ab"] = AB

for i, key in enumerate(["hetero", "dna", "rna", "ligand", "ptm"]):
    for j, r in enumerate(result[key]):
        alloc[BAN[(i + j) % 4]].append((SPEC[key], r))
# membrane: 2 per 班
for j, r in enumerate(result["memb"]):
    alloc[BAN[j % 4]].append((SPEC["memb"], r))
# antibody: 3 per 班 + 1 extra to 班4; spread the types
order = []
for k in range(5):
    for grp in (sel_pep, sel_dom):
        if k < len(grp):
            order.append(grp[k])
order = sel_nb[:1] + order + sel_nb[1:]
order = [r for i, r in enumerate(order) if r["id"] not in {x["id"] for x in order[:i]}]
for j, r in enumerate(order):
    alloc[BAN[j % 4]].append((SPEC["ab"], r))

print()
tot = 0
for b in BAN:
    items = alloc[b]
    tot += len(items)
    cats = Counter = {}
    for s, r in items:
        cats[s["key"]] = cats.get(s["key"], 0) + 1
    print(f"{b}: {len(items):2d} 题  {cats}")
print("TOTAL:", tot)

out = {b: [{"category": s["key"], "label": s["label"], "note": s["note"], "rec": r}
           for s, r in alloc[b]] for b in BAN}
json.dump(out, open("alloc7.json", "w"), ensure_ascii=False, indent=1)
print("wrote alloc7.json")
