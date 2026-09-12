"""Final selection: 41 structures = 4 classes x (8 categories + 2~3 antibodies).

班1/2/3 = 10 题 (8 类各 1 + 抗体 2)
班4     = 11 题 (8 类各 1 + 抗体 3)

Fairness mechanism: within each category the 4 picks are difficulty-matched
(minimal token spread), then assigned to the four 班 by a rotating offset, so no
班 systematically receives the largest example in every category.
"""
import json
import re
from collections import Counter

valid = json.load(open("census_valid.json"))
iface = json.load(open("iface.json"))

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


def P(r):
    return [c for c in r["chains"] if c["type"] == "Protein"]


def NUC(r):
    return [c for c in r["chains"] if c["type"] in ("DNA", "RNA")]


def blob(r):
    return r["title"] + " || " + " | ".join(c["desc"] for c in r["chains"])


def memb(r):
    return any(set(c["annot"]) & MEMB_ANNOT for c in r["chains"])


def abch(r):
    return [c for c in P(r) if AB_RE.search(c["desc"]) and not NB_RE.search(c["desc"])]


def nbch(r):
    return [c for c in P(r) if NB_RE.search(c["desc"])]


def agch(r):
    return [c for c in P(r) if not AB_RE.search(c["desc"])]


def is_ab(r):
    if TCR_RE.search(blob(r)):
        return False
    return (len(abch(r)) >= 2 or len(nbch(r)) >= 1) and len(agch(r)) >= 1


def cofac(r):
    return [l for l in r["ligands"] if l["code"] in BUILTIN_LIG]


def otherlig(r):
    return [l for l in r["ligands"] if l["code"] not in BUILTIN_LIG
            and (l["mw"] or 0) >= 150 and (l["atoms"] or 0) >= 10]


def metals(r):
    return [(c, n) for c, n in r["ions"] if c in TEACHING_METALS]


def instances(r):
    return sum(c["copies"] for c in r["chains"])


def resol(r):
    return r["resolution"] or 9.9


def clean(r):
    """No un-enterable features: free oligosaccharides or loose sugars."""
    return not r["branched"] and not r["sugars"]


def upkey(r):
    ups = set()
    for c in r["chains"]:
        ups |= set(c["uniprots"])
    if ups:
        return ("up", tuple(sorted(ups)))
    big = max(r["chains"], key=lambda c: c["len"])
    return ("seq", big["seq"][:40])


STOP = {"crystal", "structure", "cryo", "complex", "with", "bound", "from", "the",
        "and", "form", "state", "resolution", "angstrom", "human", "mutant",
        "variant", "wild", "type", "apo", "holo", "domain", "protein", "analysis"}


def topic(r):
    t = re.sub(r"[^a-z0-9 ]", " ", r["title"].lower())
    return set(w for w in t.split() if len(w) > 3 and w not in STOP)


# ---------- PTM interface evidence ----------
def ptm_iface(pid):
    """Best interface contact for a PTM in this entry: (min_dist, n_contacts)."""
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


# ---------- category definitions ----------
CATS = [
    dict(key="hetero", label="异源蛋白复合体",
         note="两条以上不同蛋白链形成的复合体(非抗体、不含核酸、非膜蛋白)。"
              "考察跨链共进化信号能否把界面摆对",
         test=lambda r: (len(P(r)) >= 2 and not NUC(r) and not is_ab(r) and not memb(r)
                         and not TCR_RE.search(blob(r))
                         and not cofac(r) and not otherlig(r)
                         and all(c["len"] >= 50 for c in P(r))
                         and clean(r) and resol(r) <= 2.6),
         span=(250, 1100)),
    dict(key="metal", label="蛋白 + 金属离子",
         note="金属为 Server 内置的 10 种离子之一,优先选 Zn/Cu/Fe/Mn/Co 等有明确配位几何的过渡金属。"
              "看模型能否把配位残基摆到正确距离",
         test=lambda r: (metals(r) and any(c in STRUCTURAL_METALS for c, _ in r["ions"])
                         and not NUC(r) and not is_ab(r) and not memb(r)
                         and not cofac(r) and not otherlig(r)
                         and clean(r) and resol(r) <= 2.0),
         span=(150, 700)),
    dict(key="dna", label="蛋白 + DNA",
         note="只含标准 A/C/G/T;双链需分别输入两条互补链(或用 “+ Reverse complement”)",
         test=lambda r: (any(c["type"] == "DNA" for c in r["chains"])
                         and not any(c["type"] == "RNA" for c in r["chains"])
                         and P(r) and not is_ab(r)
                         and sum(c["len"] * c["copies"] for c in r["chains"]
                                 if c["type"] == "DNA") >= 10
                         and clean(r) and resol(r) <= 2.8),
         span=(150, 900)),
    dict(key="rna", label="蛋白 + RNA",
         note="RNA 构象自由度大,是公认较难的一类",
         test=lambda r: (any(c["type"] == "RNA" for c in r["chains"])
                         and P(r) and not is_ab(r)
                         and sum(c["len"] * c["copies"] for c in r["chains"]
                                 if c["type"] == "RNA") >= 8
                         and clean(r) and resol(r) <= 3.0),
         span=(200, 1100)),
    dict(key="cofac", label="蛋白 + 内置辅因子",
         note="配体在 Server 下拉菜单的 19 种内置辅因子之列,无需手输 CCD 代码",
         test=lambda r: (cofac(r) and not otherlig(r) and not NUC(r)
                         and not is_ab(r) and not memb(r)
                         and clean(r) and resol(r) <= 2.2),
         span=(180, 900)),
    dict(key="ligand", label="蛋白 + 非内置小分子配体",
         note="配体不在内置清单里,须在 request builder 里加 “CCD Code” 条目手动输入;"
              "代码已核实存在于冻结的 2024_10_28 字典中",
         test=lambda r: (otherlig(r) and not cofac(r) and not NUC(r)
                         and not is_ab(r) and not memb(r)
                         and clean(r) and resol(r) <= 2.2),
         span=(200, 900)),
    dict(key="ptm", label="有翻译后修饰的蛋白",
         note="修饰残基经坐标实测位于蛋白-蛋白界面上(重原子最近距离 < 3 Å,"
              "4 Å 内接触原子数 ≥ 20),即修饰本身就是结合决定簇",
         test=lambda r: (r["ptms"] and len(P(r)) >= 2 and not is_ab(r)
                         and (ptm_iface(r["id"]) or (99, 0))[0] < 3.0
                         and (ptm_iface(r["id"]) or (99, 0))[1] >= 20
                         and clean(r)),
         span=(200, 1400)),
    dict(key="memb", label="膜蛋白",
         note="有 PDBTM / MemProtMD / mpstruc 跨膜注释。Server 不建模膜平面(FAQ 明示),"
              "跨膜螺旋排布与构象态最易出错",
         test=lambda r: (memb(r) and P(r) and not is_ab(r) and clean(r)),
         span=(250, 2000)),
]

AB_CAT = dict(key="ab", label="抗原抗体复合物",
              note="抗体 CDR 环构象 + 表位定位,文献公认难点。建议跑多个 seed,"
                   "按 ranking_score / chain_pair_iptm 选最优模型",
              span=(400, 1200))


def pick_matched(pool, n, span):
    """Pick n difficulty-matched, non-redundant entries: choose the window (in
    token order) with the smallest spread that can still be filled."""
    lo, hi = span
    cands = [r for r in pool if lo <= r["tokens"] <= hi]
    cands.sort(key=lambda r: r["tokens"])
    best = None
    for i in range(len(cands)):
        chosen, ups, tops = [], set(), []
        for j in range(i, len(cands)):
            r = cands[j]
            uk = upkey(r)
            tp = topic(r)
            if uk in ups:
                continue
            if any(len(tp & t) >= 3 for t in tops):
                continue
            chosen.append(r)
            ups.add(uk)
            tops.append(tp)
            if len(chosen) == n:
                break
        if len(chosen) == n:
            spread = chosen[-1]["tokens"] / max(1, chosen[0]["tokens"])
            score = (spread, sum(resol(x) for x in chosen) / n)
            if best is None or score < best[0]:
                best = (score, chosen)
    return best[1] if best else []


used = set()
result = {}

for cat in CATS:
    pool = [r for r in valid.values()
            if r["id"] not in used and cat["test"](r)]
    got = pick_matched(pool, 4, cat["span"])
    for r in got:
        used.add(r["id"])
    result[cat["key"]] = got
    toks = [r["tokens"] for r in got]
    print(f"{cat['key']:8s} pool={len(pool):5d} picked={len(got)} "
          f"tokens={toks} spread={max(toks)/min(toks):.2f}x" if got
          else f"{cat['key']:8s} pool={len(pool):5d} picked=0  <-- FAILED")

# ---------- antibodies: 8 Fab-type + 1 nanobody ----------
fab_pool = [r for r in valid.values()
            if r["id"] not in used and is_ab(r) and len(abch(r)) >= 2
            and not nbch(r) and clean(r) and resol(r) <= 3.0]
nb_pool = [r for r in valid.values()
           if r["id"] not in used and is_ab(r) and nbch(r) and clean(r)
           and resol(r) <= 3.0]
fabs = pick_matched(fab_pool, 8, AB_CAT["span"])
for r in fabs:
    used.add(r["id"])
nbs = pick_matched(nb_pool, 1, (250, 900))
for r in nbs:
    used.add(r["id"])
result["ab_fab"] = fabs
result["ab_nb"] = nbs
print(f"ab_fab   pool={len(fab_pool):5d} picked={len(fabs)} "
      f"tokens={[r['tokens'] for r in fabs]}")
print(f"ab_nb    pool={len(nb_pool):5d} picked={len(nbs)} "
      f"tokens={[r['tokens'] for r in nbs]}")

# ---------- assign to the four 班 ----------
CLASSES = ["班1", "班2", "班3", "班4"]
alloc = {c: [] for c in CLASSES}

for i, cat in enumerate(CATS):
    got = sorted(result[cat["key"]], key=lambda r: r["tokens"])
    for j, r in enumerate(got):
        cls = CLASSES[(i + j) % 4]          # rotate so no 班 always gets the biggest
        alloc[cls].append((cat, r))

fabs_sorted = sorted(result["ab_fab"], key=lambda r: r["tokens"])
for j, r in enumerate(fabs_sorted):
    alloc[CLASSES[j % 4]].append((AB_CAT, r))
for r in result["ab_nb"]:
    alloc["班4"].append((AB_CAT, r))

print()
total = 0
for c in CLASSES:
    items = alloc[c]
    total += len(items)
    tk = sum(r["tokens"] for _, r in items)
    print(f"{c}: {len(items)} 题, token 合计 {tk}, 平均 {tk // len(items)}")
print("TOTAL:", total)

out = {}
for c in CLASSES:
    out[c] = [{"category": cat["key"], "label": cat["label"],
               "note": cat["note"], "rec": r} for cat, r in alloc[c]]
json.dump(out, open("alloc41.json", "w"), ensure_ascii=False, indent=1)
print("wrote alloc41.json")
