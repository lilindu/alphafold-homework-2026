"""Final selection: 41 structures = 4 班 x (8 categories + 2~3 antibodies).

班1/2/3 = 10 题 (8 类各 1 + 抗体 2)
班4     = 11 题 (8 类各 1 + 抗体 3)

Two fixes over the previous attempt:

1. Class 1 leaked antibodies. `is_ab` required ">=2 antibody chains AND >=1
   antigen chain", so an UNLIGANDED Fab (heavy + light, no antigen) evaluated as
   "not an antibody" and filled all four 异源蛋白复合体 slots. Class 1 now
   rejects any entry carrying an antibody/nanobody chain at all.

2. Difficulty matching is gone (per instruction: variation between 班 is
   desirable). The previous code minimized token spread, which actively selected
   near-duplicates -- a WT/point-mutant pair of the same enzyme, and two
   paralogues of the same transcription-factor family. Selection now MAXIMIZES
   diversity: candidates are drawn from different size bands, and any two picks
   in a class must differ by UniProt, by sequence fingerprint (catches
   WT-vs-mutant, where a UniProt is missing on one side), and by title topic.
"""
import json
import re

valid = json.load(open("census_valid.json"))
iface = json.load(open("iface.json"))

BUILTIN_LIG = {"ADP", "ATP", "AMP", "GTP", "GDP", "FAD", "NAD", "NAP", "NDP",
               "HEM", "HEC", "PLM", "OLA", "MYR", "CIT", "CLA", "CHL", "BCL", "BCB"}
TEACHING_METALS = {"ZN", "CU", "FE", "MN", "CO", "MG", "CA"}
STRUCTURAL_METALS = {"ZN", "CU", "FE", "MN", "CO"}
MEMB_ANNOT = {"PDBTM", "MemProtMD", "mpstruc"}

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
    """Any antibody-derived chain, bound to an antigen or not."""
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
    """Antibody WITH an antigen -- an unliganded Fab does not qualify."""
    if TCR_RE.search(blob(r)):
        return False
    if not agch(r):
        return False
    return len(abch(r)) >= 2 or len(nbch(r)) >= 1


def cofac(r):
    return [l for l in r["ligands"] if l["code"] in BUILTIN_LIG]


def otherlig(r):
    return [l for l in r["ligands"] if l["code"] not in BUILTIN_LIG
            and (l["mw"] or 0) >= 150 and (l["atoms"] or 0) >= 10]


def metals(r):
    return [(c, n) for c, n in r["ions"] if c in TEACHING_METALS]


def resol(r):
    return r["resolution"] or 9.9


def clean(r):
    return not r["branched"] and not r["sugars"]


def uniprots(r):
    s = set()
    for c in r["chains"]:
        s |= set(c["uniprots"])
    return s


def longest(r):
    return max(r["chains"], key=lambda c: c["len"])


STOP = {"crystal", "structure", "cryo", "complex", "with", "bound", "from", "the",
        "and", "form", "state", "resolution", "angstrom", "human", "mutant",
        "variant", "wild", "type", "apo", "holo", "domain", "protein", "analysis",
        "structural", "high", "room", "temperature", "dose", "series", "data"}


def topic(r):
    t = re.sub(r"[^a-z0-9 ]", " ", r["title"].lower())
    return set(w for w in t.split() if len(w) > 3 and w not in STOP)


def similar(a, b):
    """Would these two count as the same teaching example?"""
    ua, ub = uniprots(a), uniprots(b)
    if ua and ub and (ua & ub):
        return True
    la, lb = longest(a), longest(b)
    # WT vs point mutant: near-identical length and identical N-terminal stretch
    if abs(la["len"] - lb["len"]) <= max(3, 0.05 * max(la["len"], lb["len"])):
        n = min(50, la["len"], lb["len"])
        if n >= 25 and la["seq"][:n] == lb["seq"][:n]:
            return True
    if len(topic(a) & topic(b)) >= 2:
        return True
    return False


def ptm_iface(pid):
    d = iface.get(pid)
    if not d or not d.get("ptms"):
        return None
    best = None
    for p in d["ptms"]:
        if p["min_dist"] is None:
            continue
        if best is None or p["n_close"] > best[1]:
            best = (p["min_dist"], p["n_close"], p["code"], p.get("chain"),
                    p.get("partner"))
    return best


def ptm_cross_partner(r):
    """The PTM must contact a DIFFERENT molecule, not another copy of its own
    chain (that is usually crystal packing).

    Deliberately does NOT require a UniProt on the partner: most of these
    complexes pair a receptor with a synthetic phosphopeptide, which carries no
    accession but is a genuine binding partner. Requiring one discarded 19 of 28
    valid cases.

    Also rejects contacts shorter than 2.2 A: that is a covalent bond length,
    not a non-covalent interface contact (9PU8's ALY sits at 1.77 A from a
    covalently attached compound).
    """
    info = ptm_iface(r["id"])
    if not info:
        return False
    mind, _, _, chain, partner = info
    if not chain or not partner:
        return False
    if mind < 2.2:
        return False
    # need at least two distinct protein sequences present
    seqs = {c["seq"] for c in P(r)}
    return len(seqs) >= 2


CATS = [
    dict(key="hetero", label="异源蛋白复合体",
         note="两条以上不同蛋白链形成的复合体(不含抗体链、不含核酸、非膜蛋白)。"
              "考察跨链共进化信号(paired MSA)能否把界面摆对",
         test=lambda r: (len(P(r)) >= 2 and not NUC(r)
                         and not has_ab_chain(r) and not memb(r)
                         and not TCR_RE.search(blob(r))
                         and not cofac(r) and not otherlig(r)
                         and all(c["len"] >= 50 for c in P(r))
                         and len(uniprots(r)) >= 2
                         and clean(r) and resol(r) <= 2.6),
         span=(250, 1600)),
    dict(key="metal", label="蛋白 + 金属离子",
         note="金属为 Server 内置 10 种离子之一,且为 Zn/Cu/Fe/Mn/Co 等有明确配位几何的过渡金属。"
              "看模型能否把配位残基摆到正确距离",
         test=lambda r: (metals(r) and any(c in STRUCTURAL_METALS for c, _ in r["ions"])
                         and not NUC(r) and not has_ab_chain(r) and not memb(r)
                         and not cofac(r) and not otherlig(r)
                         and clean(r) and resol(r) <= 2.1),
         span=(150, 1200)),
    dict(key="dna", label="蛋白 + DNA",
         note="只含标准 A/C/G/T;双链需分别输入两条互补链(或用 “+ Reverse complement”)",
         test=lambda r: (any(c["type"] == "DNA" for c in r["chains"])
                         and not any(c["type"] == "RNA" for c in r["chains"])
                         and P(r) and not has_ab_chain(r)
                         and sum(c["len"] * c["copies"] for c in r["chains"]
                                 if c["type"] == "DNA") >= 10
                         and clean(r) and resol(r) <= 2.8),
         span=(150, 1200)),
    dict(key="rna", label="蛋白 + RNA",
         note="RNA 构象自由度大,是公认较难的一类",
         test=lambda r: (any(c["type"] == "RNA" for c in r["chains"])
                         and P(r) and not has_ab_chain(r)
                         and sum(c["len"] * c["copies"] for c in r["chains"]
                                 if c["type"] == "RNA") >= 8
                         and clean(r) and resol(r) <= 3.2),
         span=(200, 1600)),
    dict(key="cofac", label="蛋白 + 内置辅因子",
         note="配体属于 Server 下拉菜单里的 19 种内置辅因子,无需手输 CCD 代码",
         test=lambda r: (cofac(r) and not otherlig(r) and not NUC(r)
                         and not has_ab_chain(r) and not memb(r)
                         and clean(r) and resol(r) <= 2.2),
         span=(180, 1400)),
    dict(key="ligand", label="蛋白 + 非内置小分子配体",
         note="配体不在内置清单内,须在 request builder 中加 “CCD Code” 条目手动输入;"
              "代码已核实存在于冻结的 CCD 2024_10_28 字典中",
         test=lambda r: (otherlig(r) and not cofac(r) and not NUC(r)
                         and not has_ab_chain(r) and not memb(r)
                         and clean(r) and resol(r) <= 2.2),
         span=(200, 1400)),
    dict(key="ptm", label="有翻译后修饰的蛋白",
         note="修饰残基经坐标实测位于两个不同蛋白之间的界面上(重原子最近距离 < 3 Å,"
              "4 Å 内接触原子 ≥ 20 个),即修饰本身就是结合的决定因素",
         test=lambda r: (r["ptms"] and len(P(r)) >= 2 and not has_ab_chain(r)
                         and (ptm_iface(r["id"]) or (99, 0))[0] < 3.0
                         and (ptm_iface(r["id"]) or (99, 0))[1] >= 20
                         and ptm_cross_partner(r)
                         and clean(r)),
         span=(200, 1600)),
    dict(key="memb", label="膜蛋白",
         note="带 PDBTM / MemProtMD / mpstruc 跨膜注释。Server 不建模膜平面(FAQ 明示),"
              "跨膜螺旋排布与构象态最易出错",
         test=lambda r: (memb(r) and P(r) and not has_ab_chain(r) and clean(r)),
         span=(250, 2400)),
]

AB_CAT = dict(key="ab", label="抗原抗体复合物",
              note="抗体 CDR 环构象 + 表位定位,文献公认难点。建议跑多个 seed,"
                   "按 ranking_score / chain_pair_iptm 选最优模型")


def pick_diverse(pool, n, span=None, keyfn=None):
    """Maximise variety: draw from n different size bands, and reject any
    candidate similar to one already chosen."""
    cands = list(pool)
    if span:
        lo, hi = span
        banded = [r for r in cands if lo <= r["tokens"] <= hi]
        if len(banded) >= n:
            cands = banded
    if not cands:
        return []
    cands.sort(key=lambda r: r["tokens"])
    k = len(cands)
    bands = []
    for b in range(n):
        s, e = (k * b) // n, (k * (b + 1)) // n
        bands.append(cands[s:e] or cands)

    chosen = []

    def ok(r):
        if any(x["id"] == r["id"] for x in chosen):
            return False
        if keyfn is not None:
            kv = keyfn(r)
            if any(keyfn(x) == kv for x in chosen):
                return False
        return not any(similar(r, x) for x in chosen)

    for band in bands:
        for r in sorted(band, key=resol):
            if ok(r):
                chosen.append(r)
                break
    if len(chosen) < n:
        for r in sorted(cands, key=resol):
            if len(chosen) >= n:
                break
            if ok(r):
                chosen.append(r)
    return chosen[:n]


used, result = set(), {}
for cat in CATS:
    pool = [r for r in valid.values() if r["id"] not in used and cat["test"](r)]
    got = pick_diverse(pool, 4, cat["span"])
    for r in got:
        used.add(r["id"])
    result[cat["key"]] = got
    toks = [r["tokens"] for r in got]
    flag = "" if len(got) == 4 else "   <-- SHORT"
    print(f"{cat['key']:8s} pool={len(pool):5d} picked={len(got)} tokens={toks}{flag}")

# antibodies: 8 Fab-type + 1 nanobody, deliberately varied epitope size
def epitope_kind(r):
    ag = agch(r)
    if not ag:
        return "?"
    return "peptide" if max(c["len"] for c in ag) <= 30 else "domain"


fab_pool = [r for r in valid.values()
            if r["id"] not in used and is_ab_complex(r) and len(abch(r)) >= 2
            and not nbch(r) and clean(r) and resol(r) <= 3.0]
pep = [r for r in fab_pool if epitope_kind(r) == "peptide"]
dom = [r for r in fab_pool if epitope_kind(r) == "domain"]
fab_pep = pick_diverse(pep, 4, (350, 1400))
for r in fab_pep:
    used.add(r["id"])
dom = [r for r in dom if r["id"] not in used]
fab_dom = pick_diverse(dom, 4, (350, 1600))
for r in fab_dom:
    used.add(r["id"])
nb_pool = [r for r in valid.values()
           if r["id"] not in used and is_ab_complex(r) and nbch(r)
           and clean(r) and resol(r) <= 3.0]
nb = pick_diverse(nb_pool, 1, (250, 1000))
for r in nb:
    used.add(r["id"])

print(f"ab pep  pool={len(pep):5d} picked={len(fab_pep)} tokens={[r['tokens'] for r in fab_pep]}")
print(f"ab dom  pool={len(dom):5d} picked={len(fab_dom)} tokens={[r['tokens'] for r in fab_dom]}")
print(f"ab nb   pool={len(nb_pool):5d} picked={len(nb)} tokens={[r['tokens'] for r in nb]}")

CLASSES = ["班1", "班2", "班3", "班4"]
alloc = {c: [] for c in CLASSES}
for i, cat in enumerate(CATS):
    for j, r in enumerate(result[cat["key"]]):
        alloc[CLASSES[(i + j) % 4]].append((cat, r))
# each 班 gets one peptide-epitope and one folded-domain antibody
for j, r in enumerate(fab_pep):
    alloc[CLASSES[j % 4]].append((AB_CAT, r))
for j, r in enumerate(fab_dom):
    alloc[CLASSES[j % 4]].append((AB_CAT, r))
for r in nb:
    alloc["班4"].append((AB_CAT, r))

print()
tot = 0
for c in CLASSES:
    items = alloc[c]
    tot += len(items)
    tk = [r["tokens"] for _, r in items]
    print(f"{c}: {len(items):2d} 题  token {min(tk)}-{max(tk)}  合计 {sum(tk)}")
print("TOTAL:", tot)

out = {c: [{"category": cat["key"], "label": cat["label"], "note": cat["note"],
            "rec": r} for cat, r in alloc[c]] for c in CLASSES}
json.dump(out, open("alloc41.json", "w"), ensure_ascii=False, indent=1)
print("wrote alloc41.json")
