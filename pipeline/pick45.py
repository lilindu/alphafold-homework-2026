"""Stage 7: pick exactly 45 non-redundant targets across three difficulty tiers.

Fixes over the first pass:
  * antibody/immunoglobulin entries no longer leak into non-antibody groups
    (checks title as well as chain descriptions)
  * de-duplication by UniProt accession set + family key, so near-identical
    MHC/TCR/FABP entries cannot fill several slots
  * membrane group requires a real membrane annotation (PDBTM/MemProtMD/mpstruc)
    rather than a keyword that also matches soluble ectodomains
  * cofactor group draws on the dedicated built-in-ligand pools
"""
import json
import re

valid = json.load(open("valid.json"))
meta = json.load(open("meta_all.json"))
pools = json.load(open("pools.json"))

BUILTIN_LIG = {"ADP", "ATP", "AMP", "GTP", "GDP", "FAD", "NAD", "NAP", "NDP",
               "HEM", "HEC", "PLM", "OLA", "MYR", "CIT", "CLA", "CHL", "BCL", "BCB"}

AB = re.compile(r"heavy chain|light chain|\bfab\b|\bf\(ab\|\bigg\b|antibody|immunoglobulin"
                r"|\bscfv\b|\bnanobody\b|\bvhh\b|single[- ]domain antibody|sybody"
                r"|\bmab\b|\bfv fragment\b", re.I)
NB = re.compile(r"nanobody|\bvhh\b|single[- ]domain antibody|sybody|camelid", re.I)
TCR = re.compile(r"\bT[- ]cell receptor\b|\bTCR\b|beta-2 microglobulin|\bMHC\b|\bHLA\b"
                 r"|histocompatibility", re.I)
MEMB_ANNOT = {"PDBTM", "MemProtMD", "mpstruc"}


def m(info):
    return meta[info["id"]]


def title(info):
    return info.get("title", "")


def descs(info):
    return " | ".join(c["desc"] for c in info["chains"])


def blob(info):
    return title(info) + " || " + descs(info)


def prot_ents(info):
    return sum(1 for c in info["chains"] if c["type"] == "Protein")


def nuc_ents(info):
    return sum(1 for c in info["chains"] if c["type"] in ("DNA", "RNA"))


def total_res(info):
    return sum(c["len"] * c["copies"] for c in info["chains"])


def uniprots(info):
    acc = set()
    for pe in m(info).get("polymer_entities") or []:
        for u in pe.get("uniprots") or []:
            if u.get("rcsb_id"):
                acc.add(u["rcsb_id"])
    return acc


def annots(info):
    out = set()
    for pe in m(info).get("polymer_entities") or []:
        for a in pe.get("rcsb_polymer_entity_annotation") or []:
            if a.get("type"):
                out.add(a["type"])
    return out


def is_membrane(info):
    return bool(annots(info) & MEMB_ANNOT)


def has_ab(info):
    return bool(AB.search(blob(info)))


def is_nb(info):
    if not NB.search(blob(info)):
        return False
    nb = [c for c in info["chains"] if NB.search(c["desc"])]
    other = [c for c in info["chains"]
             if c["type"] == "Protein" and not NB.search(c["desc"])]
    return len(nb) >= 1 and len(other) >= 1


def is_ab_complex(info):
    ch = info["chains"]
    abc = [c for c in ch if AB.search(c["desc"]) and not NB.search(c["desc"])]
    ag = [c for c in ch if c["type"] == "Protein" and not AB.search(c["desc"])]
    return len(abc) >= 2 and len(ag) >= 1


def is_tcr_mhc(info):
    return bool(TCR.search(blob(info)))


def ligs(info):
    return info["ligands"]


def builtin_only(info):
    return bool(ligs(info)) and all(c in BUILTIN_LIG for c, *_ in ligs(info))


def resol(info):
    return info["resolution"] or 9.9


def famkey(info):
    up = uniprots(info)
    if up:
        return ("up", tuple(sorted(up)))
    c = max(info["chains"], key=lambda x: x["len"])
    return ("seq", c["seq"][:40])


def wordkey(info):
    t = re.sub(r"[^a-z0-9 ]", " ", title(info).lower())
    drop = {"crystal", "structure", "cryo", "em", "of", "the", "a", "an", "in",
            "complex", "with", "bound", "to", "at", "and", "form", "state",
            "resolution", "angstrom", "ray", "x", "from", "its", "by", "for",
            "human", "mutant", "variant", "wild", "type", "apo", "holo"}
    ws = [w for w in t.split() if len(w) > 3 and w not in drop]
    return tuple(sorted(set(ws))[:5])


def pool(*names):
    seen = {}
    for n in names:
        for i in pools.get(n, []):
            if i in valid:
                seen[i] = valid[i]
    return list(seen.values())


LIGPOOLS = [f"lig_{x}" for x in
            ["ATP", "ADP", "AMP", "GTP", "GDP", "FAD", "NAD", "NAP", "NDP",
             "HEM", "HEC", "PLM", "OLA", "MYR", "CIT"]]

SPECS = [
    # ---------------- Tier A ----------------
    dict(key="A1", n=5, tier="A",
         label="单域单体蛋白",
         note="单链、无配体、超高分辨率。MSA 深、折叠常见,预期 pLDDT>90",
         src=["mono_protein"],
         test=lambda i: (prot_ents(i) == 1 and nuc_ents(i) == 0
                         and len(i["chains"]) == 1 and i["chains"][0]["copies"] == 1
                         and not ligs(i) and not has_ab(i)
                         and 130 <= total_res(i) <= 420 and resol(i) <= 1.9)),
    dict(key="A2", n=4, tier="A",
         label="同源二聚体",
         note="单一序列 ×2,界面对称,AlphaFold-Multimer 的强项",
         src=["homodimer", "homodimer2"],
         test=lambda i: (prot_ents(i) == 1 and nuc_ents(i) == 0
                         and i["chains"][0]["copies"] == 2
                         and not has_ab(i) and not is_tcr_mhc(i)
                         and 230 <= total_res(i) <= 900 and resol(i) <= 2.2)),
    dict(key="A3", n=2, tier="A",
         label="同源四聚体及以上",
         note="单一序列多拷贝,考察对称寡聚体组装",
         src=["homotetramer"],
         test=lambda i: (prot_ents(i) == 1 and nuc_ents(i) == 0
                         and i["chains"][0]["copies"] >= 4
                         and not has_ab(i)
                         and 320 <= total_res(i) <= 1600 and resol(i) <= 2.4)),
    dict(key="A4", n=4, tier="A",
         label="单体 + 内置辅因子",
         note="配体全部在 Server 下拉列表内(ATP/HEM/FAD/NAD…),学生无需手输 CCD",
         src=LIGPOOLS + ["mono_with_cofactor"],
         test=lambda i: (prot_ents(i) == 1 and nuc_ents(i) == 0
                         and len(i["chains"]) == 1 and i["chains"][0]["copies"] <= 2
                         and builtin_only(i) and len(ligs(i)) <= 2
                         and not has_ab(i)
                         and 150 <= total_res(i) <= 700 and resol(i) <= 2.2)),
    # ---------------- Tier B ----------------
    dict(key="B1", n=4, tier="B",
         label="异源二聚体",
         note="两个不同蛋白的界面,依赖跨链共进化信号",
         src=["hetero2", "peptide_complex"],
         test=lambda i: (prot_ents(i) == 2 and nuc_ents(i) == 0
                         and not has_ab(i) and not is_tcr_mhc(i)
                         and min(c["len"] for c in i["chains"]) >= 45
                         and 230 <= total_res(i) <= 1100 and resol(i) <= 2.6)),
    dict(key="B2", n=2, tier="B",
         label="异源三元复合物",
         note="三个不同亚基,界面增多,ipTM 更易下降",
         src=["hetero3"],
         test=lambda i: (prot_ents(i) == 3 and nuc_ents(i) == 0
                         and not has_ab(i) and not is_tcr_mhc(i)
                         and 300 <= total_res(i) <= 1600 and resol(i) <= 3.0)),
    dict(key="B3", n=2, tier="B",
         label="蛋白 + 短肽",
         note="短肽 (<30 aa) 结合槽。pTM 对短链系统性偏低,须改看 pLDDT/PAE",
         src=["peptide_complex", "hetero2"],
         test=lambda i: (prot_ents(i) == 2 and nuc_ents(i) == 0
                         and not has_ab(i) and not is_tcr_mhc(i)
                         and min(c["len"] for c in i["chains"]) <= 30
                         and 160 <= total_res(i) <= 800 and resol(i) <= 2.5)),
    dict(key="B4", n=4, tier="B",
         label="蛋白–DNA 复合物",
         note="仅标准 A/C/G/T。双链需按 FAQ 分别输入两条互补链",
         src=["protein_dna"],
         test=lambda i: (nuc_ents(i) >= 1 and prot_ents(i) >= 1
                         and all(c["type"] != "RNA" for c in i["chains"])
                         and not has_ab(i)
                         and 120 <= total_res(i) <= 1200)),
    dict(key="B5", n=4, tier="B",
         label="蛋白–RNA 复合物",
         note="RNA 构象自由度大,公认较难的一类",
         src=["protein_rna", "protein_rna2"],
         test=lambda i: (any(c["type"] == "RNA" for c in i["chains"])
                         and prot_ents(i) >= 1 and not has_ab(i)
                         and 120 <= total_res(i) <= 1300)),
    # ---------------- Tier C ----------------
    dict(key="C1", n=5, tier="C",
         label="抗体 / Fab–抗原复合物",
         note="CDR 环构象 + 表位定位,文献公认难点。务必跑多个 seed 并按 ipTM 排序",
         src=["antibody", "nanobody2"],
         test=lambda i: (is_ab_complex(i) and not is_nb(i)
                         and 400 <= total_res(i) <= 1300 and resol(i) <= 3.2)),
    dict(key="C2", n=2, tier="C",
         label="纳米抗体 / VHH–抗原复合物",
         note="单域抗体,链数少但表位仍难定位",
         src=["nanobody", "nanobody2", "antibody"],
         test=lambda i: (is_nb(i) and 250 <= total_res(i) <= 900 and resol(i) <= 3.0)),
    dict(key="C3", n=3, tier="C",
         label="膜蛋白(通道 / 转运体 / 受体)",
         note="Server 不知道膜平面(FAQ 明示),跨膜螺旋排布与构象态易错",
         src=["memb_PDBTM", "memb_MemProtMD", "memb_mpstruc",
              "kw_gpcr", "kw_channel", "kw_transporter", "kw_abc", "membrane_kw"],
         test=lambda i: (is_membrane(i) and prot_ents(i) >= 1
                         and not has_ab(i)
                         and 250 <= total_res(i) <= 2200)),
    dict(key="C4", n=2, tier="C",
         label="大型多亚基复合物",
         note="≥5 种不同链,界面多、token 接近上限",
         src=["large_assembly", "signaling_gpcr"],
         test=lambda i: (prot_ents(i) >= 5 and nuc_ents(i) == 0
                         and not has_ab(i)
                         and 700 <= total_res(i) <= 3200)),
    dict(key="C5", n=2, tier="C",
         label="纯核酸结构(无蛋白)",
         note="没有蛋白 MSA 支撑,AlphaFold 3 在此类最弱",
         src=["nucleic_only", "nucleic_only2"],
         test=lambda i: (prot_ents(i) == 0 and nuc_ents(i) >= 1
                         and 25 <= total_res(i) <= 300)),
]

assert sum(s["n"] for s in SPECS) == 45, sum(s["n"] for s in SPECS)

chosen = []
used_id, used_fam, used_word = set(), set(), set()

for sp in SPECS:
    cands = [i for i in pool(*sp["src"]) if sp["test"](i)]
    cands.sort(key=lambda i: (resol(i), -total_res(i)))
    got = 0
    for info in cands:
        if got >= sp["n"]:
            break
        if info["id"] in used_id:
            continue
        fk, wk = famkey(info), wordkey(info)
        if fk in used_fam:
            continue
        if wk and wk in used_word:
            continue
        used_id.add(info["id"])
        used_fam.add(fk)
        if wk:
            used_word.add(wk)
        rec = dict(info)
        rec.update(group=sp["key"], group_label=sp["label"], tier=sp["tier"],
                   group_note=sp["note"])
        chosen.append(rec)
        got += 1
    print(f"{sp['key']} {sp['label'][:20]:22s} pool={len(cands):4d} picked={got}/{sp['n']}")

print("\nTOTAL", len(chosen))
json.dump(chosen, open("chosen45.json", "w"), indent=1)
