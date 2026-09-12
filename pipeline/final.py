"""Stage 5: pick 45 non-redundant targets across difficulty tiers."""
import json
import re
from collections import defaultdict

valid = json.load(open("valid.json"))
pools = json.load(open("pools.json"))

BUILTIN_LIG = {"ADP", "ATP", "AMP", "GTP", "GDP", "FAD", "NAD", "NAP", "NDP",
               "HEM", "HEC", "PLM", "OLA", "MYR", "CIT", "CLA", "CHL", "BCL", "BCB"}

AB_WORDS = re.compile(r"heavy chain|light chain|\bfab\b|\bigg\b|antibody|immunoglobulin|\bscfv\b|\bfv\b", re.I)
NB_WORDS = re.compile(r"nanobody|\bvhh\b|single[- ]domain antibody|sybody", re.I)
MEMB_WORDS = re.compile(r"channel|transporter|receptor|permease|ATPase|pump|porin|translocase|symporter|antiporter|exchanger|rhodopsin|integrase|flippase", re.I)


def longest_chain(info, ptype="Protein"):
    cs = [c for c in info["chains"] if c["type"] == ptype]
    return max(cs, key=lambda c: c["len"]) if cs else None


def descs(info):
    return " | ".join(c["desc"] for c in info["chains"])


def famkey(info):
    """Crude family key for de-duplication."""
    c = longest_chain(info) or info["chains"][0]
    return (c["seq"][:45], c["len"] // 25)


def titlekey(info):
    t = re.sub(r"[^a-z0-9 ]", " ", info["title"].lower())
    t = re.sub(r"\b(crystal|structure|cryo|em|of|the|a|an|in|complex|with|bound|to|at|and|form|state|resolution|angstrom|x|ray)\b", " ", t)
    return tuple(sorted(set(w for w in t.split() if len(w) > 3))[:6])


def is_ab(info):
    ch = info["chains"]
    ab = [c for c in ch if AB_WORDS.search(c["desc"])]
    non = [c for c in ch if not AB_WORDS.search(c["desc"]) and not NB_WORDS.search(c["desc"]) and c["type"] == "Protein"]
    return len(ab) >= 2 and len(non) >= 1


def is_nb(info):
    ch = info["chains"]
    nb = [c for c in ch if NB_WORDS.search(c["desc"])]
    non = [c for c in ch if not NB_WORDS.search(c["desc"]) and c["type"] == "Protein"]
    return len(nb) >= 1 and len(non) >= 1


def n_real_lig(info):
    return len(info["ligands"])


def has_builtin_only(info):
    return all(c in BUILTIN_LIG for c, *_ in info["ligands"])


def prot_entities(info):
    return sum(1 for c in info["chains"] if c["type"] == "Protein")


def nuc_entities(info):
    return sum(1 for c in info["chains"] if c["type"] in ("DNA", "RNA"))


def total_res(info):
    return sum(c["len"] * c["copies"] for c in info["chains"])


def cand(cat):
    return [valid[i] for i in pools.get(cat, []) if i in valid]


# ------------------------------------------------------------------
CATEGORY_SPECS = []


def spec(key, label, tier, note, source_cats, test, count, sort=None):
    CATEGORY_SPECS.append(dict(key=key, label=label, tier=tier, note=note,
                               cats=source_cats, test=test, count=count,
                               sort=sort or (lambda i: (i["resolution"] or 9))))


# ---------- Tier A: expected easier ----------
spec("A-mono", "单域单体蛋白 (single-domain monomer)", "A",
     "单链、无配体、高分辨率;MSA 深、折叠常见,是 AlphaFold 的强项",
     ["mono_protein"],
     lambda i: (prot_entities(i) == 1 and nuc_entities(i) == 0
                and len(i["chains"]) == 1 and i["chains"][0]["copies"] == 1
                and n_real_lig(i) == 0 and 130 <= total_res(i) <= 420
                and (i["resolution"] or 9) <= 1.8),
     5)

spec("A-homodimer", "同源二聚体 (homodimer)", "A",
     "单一序列 ×2,界面对称,AlphaFold-Multimer 一般表现良好",
     ["homodimer", "homotetramer"],
     lambda i: (prot_entities(i) == 1 and nuc_entities(i) == 0
                and i["chains"][0]["copies"] == 2
                and 220 <= total_res(i) <= 900
                and (i["resolution"] or 9) <= 2.2),
     3)

spec("A-oligomer", "同源四聚体/多聚体 (homo-oligomer)", "A",
     "单一序列多拷贝,考察对称寡聚体的组装",
     ["homotetramer"],
     lambda i: (prot_entities(i) == 1 and nuc_entities(i) == 0
                and i["chains"][0]["copies"] >= 4
                and 320 <= total_res(i) <= 1600),
     2)

spec("A-cofactor", "单体 + 内置辅因子/离子", "A",
     "配体全部在 Server 下拉列表内(ATP/HEM/FAD/NAD 等或常见离子),无需手输 CCD",
     ["mono_with_cofactor", "mono_protein"],
     lambda i: (prot_entities(i) == 1 and nuc_entities(i) == 0
                and len(i["chains"]) == 1
                and 1 <= n_real_lig(i) <= 2 and has_builtin_only(i)
                and 150 <= total_res(i) <= 600
                and (i["resolution"] or 9) <= 2.2),
     5)

# ---------- Tier B: medium ----------
spec("B-hetero2", "异源二聚体 (2 条不同链)", "B",
     "两个不同蛋白的界面,需要跨链共进化信号",
     ["hetero2", "peptide_complex"],
     lambda i: (prot_entities(i) == 2 and nuc_entities(i) == 0
                and not is_ab(i) and not is_nb(i)
                and all(c["len"] >= 40 for c in i["chains"])
                and 220 <= total_res(i) <= 1100
                and (i["resolution"] or 9) <= 2.6),
     4)

spec("B-hetero3", "异源三元复合物 (3 条不同链)", "B",
     "三个不同亚基,界面数量增加,ipTM 更容易掉",
     ["hetero3"],
     lambda i: (prot_entities(i) == 3 and nuc_entities(i) == 0
                and not is_ab(i) and not is_nb(i)
                and 300 <= total_res(i) <= 1600),
     3)

spec("B-peptide", "蛋白 + 短肽 (peptide in groove)", "B",
     "短肽链 (<30 aa) 结合位点;pTM 对短链偏低,须看 pLDDT/PAE",
     ["peptide_complex", "hetero2"],
     lambda i: (prot_entities(i) == 2 and nuc_entities(i) == 0
                and not is_ab(i) and not is_nb(i)
                and min(c["len"] for c in i["chains"]) <= 30
                and 160 <= total_res(i) <= 800),
     2)

spec("B-dna", "蛋白–DNA 复合物", "B",
     "标准 A/C/G/T;双链需分别输入两条互补链",
     ["protein_dna"],
     lambda i: (nuc_entities(i) >= 1 and prot_entities(i) >= 1
                and all(c["type"] != "RNA" for c in i["chains"])
                and 120 <= total_res(i) <= 1200),
     4)

spec("B-rna", "蛋白–RNA 复合物", "B",
     "RNA 构象自由度大,是公认较难的一类",
     ["protein_rna"],
     lambda i: (any(c["type"] == "RNA" for c in i["chains"])
                and prot_entities(i) >= 1
                and 120 <= total_res(i) <= 1200),
     2)

# ---------- Tier C: hard ----------
spec("C-ab", "抗体/Fab–抗原复合物", "C",
     "CDR 环 + 表位定位;文献公认难点,务必跑多个 seed",
     ["antibody"],
     lambda i: (is_ab(i) and 400 <= total_res(i) <= 1300),
     5)

spec("C-nb", "纳米抗体/VHH–抗原复合物", "C",
     "单域抗体,链数少但表位仍难定位",
     ["nanobody", "antibody"],
     lambda i: (is_nb(i) and 250 <= total_res(i) <= 900),
     2)

spec("C-membrane", "膜蛋白 (通道/转运体/受体)", "C",
     "Server 不知道膜平面,跨膜螺旋排布易错",
     ["membrane_kw", "transport_kw"],
     lambda i: (MEMB_WORDS.search(descs(i)) is not None
                and prot_entities(i) >= 1
                and 250 <= total_res(i) <= 2000),
     3)

spec("C-assembly", "大型多亚基复合物 (≥5 种链)", "C",
     "亚基多、界面多,token 接近上限",
     ["large_assembly", "signaling_gpcr"],
     lambda i: (prot_entities(i) >= 5 and nuc_entities(i) == 0
                and 700 <= total_res(i) <= 3200),
     3)

spec("C-rnafold", "纯核酸结构 (RNA/DNA,无蛋白)", "C",
     "无蛋白、无 MSA 支撑,AlphaFold 3 在此类最弱",
     ["nucleic_only"],
     lambda i: (prot_entities(i) == 0 and nuc_entities(i) >= 1
                and 25 <= total_res(i) <= 250),
     2)

# ------------------------------------------------------------------
chosen, used_fam, used_title, used_id = [], set(), set(), set()

for sp in CATEGORY_SPECS:
    seen = {}
    for c in sp["cats"]:
        for info in cand(c):
            seen[info["id"]] = info
    pool = [i for i in seen.values() if sp["test"](i)]
    pool.sort(key=sp["sort"])
    picked = 0
    for info in pool:
        if picked >= sp["count"]:
            break
        if info["id"] in used_id:
            continue
        fk, tk = famkey(info), titlekey(info)
        if fk in used_fam or (tk and tk in used_title):
            continue
        used_id.add(info["id"]); used_fam.add(fk); used_title.add(tk)
        rec = dict(info)
        rec["group"] = sp["key"]; rec["group_label"] = sp["label"]
        rec["tier"] = sp["tier"]; rec["group_note"] = sp["note"]
        chosen.append(rec); picked += 1
    print(f"{sp['key']:14s} pool={len(pool):4d} picked={picked}/{sp['count']}")

print("\nTOTAL:", len(chosen))
json.dump(chosen, open("chosen.json", "w"), indent=1)

for r in chosen:
    lig = ",".join(f"{c}x{n}" for c, n, *_ in r["ligands"]) or "-"
    ion = ",".join(f"{c}x{n}" for c, n in r["ions"]) or "-"
    print(f"{r['group']:12s} {r['id']}  res={r['resolution']}  tok={r['tokens']:5d} "
          f"prot={prot_entities(r)} nuc={nuc_entities(r)} lig={lig} ion={ion}")
    print(f"      {r['title'][:110]}")
