"""Stage 8: final selection of 45, now filtering oligomeric groups on the
BIOLOGICAL ASSEMBLY rather than the asymmetric unit, and forcing subject
diversity so no group is filled with near-duplicates.
"""
import json
import re

valid = json.load(open("valid.json"))
meta = json.load(open("meta_all.json"))
pools = json.load(open("pools.json"))
asm = json.load(open("asm.json"))

BUILTIN_LIG = {"ADP", "ATP", "AMP", "GTP", "GDP", "FAD", "NAD", "NAP", "NDP",
               "HEM", "HEC", "PLM", "OLA", "MYR", "CIT", "CLA", "CHL", "BCL", "BCB"}

AB = re.compile(r"heavy chain|light chain|\bfab\b|\bigg\b|antibody|immunoglobulin"
                r"|\bscfv\b|\bnanobody\b|\bvhh\b|single[- ]domain antibody|sybody"
                r"|\bmab\b", re.I)
NB = re.compile(r"nanobody|\bvhh\b|single[- ]domain antibody|sybody|camelid", re.I)
TCRMHC = re.compile(r"\bT[- ]cell receptor\b|\bTCR\b|beta-2 microglobulin|\bMHC\b"
                    r"|\bHLA\b|histocompatibility|\bH-2[A-Z]\b", re.I)
MEMB_ANNOT = {"PDBTM", "MemProtMD", "mpstruc"}


def M(i):
    return meta[i["id"]]


def blob(i):
    return i.get("title", "") + " || " + " | ".join(c["desc"] for c in i["chains"])


def prot_ents(i):
    return sum(1 for c in i["chains"] if c["type"] == "Protein")


def nuc_ents(i):
    return sum(1 for c in i["chains"] if c["type"] in ("DNA", "RNA"))


def total_res(i):
    return sum(c["len"] * c["copies"] for c in i["chains"])


def resol(i):
    return i["resolution"] or 9.9


def uniprots(i):
    s = set()
    for pe in M(i).get("polymer_entities") or []:
        for u in pe.get("uniprots") or []:
            if u.get("rcsb_id"):
                s.add(u["rcsb_id"])
    return s


def annots(i):
    s = set()
    for pe in M(i).get("polymer_entities") or []:
        for a in pe.get("rcsb_polymer_entity_annotation") or []:
            if a.get("type"):
                s.add(a["type"])
    return s


def assemblies(i):
    out = []
    for a in asm.get(i["id"]) or []:
        ai = a.get("rcsb_assembly_info") or {}
        out.append({
            "inst": ai.get("polymer_entity_instance_count") or 0,
            "ents": ai.get("polymer_entity_count") or 0,
            "prot": ai.get("polymer_entity_instance_count_protein") or 0,
            "dna": ai.get("polymer_entity_instance_count_DNA") or 0,
            "rna": ai.get("polymer_entity_instance_count_RNA") or 0,
        })
    return out


def bio(i):
    """Largest biological assembly."""
    a = assemblies(i)
    return max(a, key=lambda x: x["inst"]) if a else None


def bio_is(i, prot_inst=None, prot_ents_n=None, min_inst=None, max_inst=None):
    b = bio(i)
    if not b:
        return False
    if prot_inst is not None and b["prot"] != prot_inst:
        return False
    if prot_ents_n is not None and b["ents"] != prot_ents_n:
        return False
    if min_inst is not None and b["inst"] < min_inst:
        return False
    if max_inst is not None and b["inst"] > max_inst:
        return False
    return True


def has_ab(i):
    return bool(AB.search(blob(i)))


def is_nb(i):
    if not NB.search(blob(i)):
        return False
    nb = [c for c in i["chains"] if NB.search(c["desc"])]
    oth = [c for c in i["chains"] if c["type"] == "Protein" and not NB.search(c["desc"])]
    return len(nb) >= 1 and len(oth) >= 1


def is_ab_complex(i):
    ab = [c for c in i["chains"] if AB.search(c["desc"]) and not NB.search(c["desc"])]
    ag = [c for c in i["chains"] if c["type"] == "Protein" and not AB.search(c["desc"])]
    return len(ab) >= 2 and len(ag) >= 1


def is_tcr(i):
    return bool(TCRMHC.search(blob(i)))


def is_memb(i):
    return bool(annots(i) & MEMB_ANNOT)


def ligs(i):
    return i["ligands"]


def builtin_only(i):
    return bool(ligs(i)) and all(c in BUILTIN_LIG for c, *_ in ligs(i))


def famkey(i):
    up = uniprots(i)
    if up:
        return ("up", tuple(sorted(up)))
    c = max(i["chains"], key=lambda x: x["len"])
    return ("seq", c["seq"][:40])


STOP = {"crystal", "structure", "cryo", "of", "the", "a", "an", "in", "complex",
        "with", "bound", "to", "at", "and", "form", "state", "resolution",
        "angstrom", "ray", "from", "its", "by", "for", "human", "mutant",
        "variant", "wild", "type", "apo", "holo", "double", "domain"}


def topic(i):
    t = re.sub(r"[^a-z0-9 ]", " ", i.get("title", "").lower())
    ws = [w for w in t.split() if len(w) > 3 and w not in STOP]
    return set(ws[:8])


# Reject huge biological assemblies (viral capsids, filaments): the deposited
# biological unit is hundreds of chains, so the experimental "answer" is not the
# few-chain protomer a student would actually submit, and it blows the 5,000
# token limit. 9TBH/9TBI (Hepatitis A capsid, 180 chains) are caught by this.
MAX_ASSEMBLY_INSTANCES = 12


def assembly_sane(i):
    b = bio(i)
    return bool(b) and 0 < b["inst"] <= MAX_ASSEMBLY_INSTANCES


def pool(*names):
    seen = {}
    for n in names:
        for x in pools.get(n, []):
            if x in valid:
                seen[x] = valid[x]
    return list(seen.values())


LIGP = [f"lig_{x}" for x in ["ATP", "ADP", "AMP", "GTP", "GDP", "FAD", "NAD",
                             "NAP", "NDP", "HEM", "HEC", "PLM", "OLA", "MYR", "CIT"]]

SPECS = [
    dict(key="A1", n=5, tier="A", label="单域单体蛋白",
         note="单链、无配体、超高分辨率;MSA 深、折叠常见,预期 pLDDT > 90。生物学单元确认为单体",
         src=["mono_protein", "mono_sym"],
         test=lambda i: (prot_ents(i) == 1 and nuc_ents(i) == 0
                         and len(i["chains"]) == 1 and i["chains"][0]["copies"] == 1
                         and not ligs(i) and not has_ab(i) and not is_tcr(i)
                         and bio_is(i, prot_inst=1)
                         and 130 <= total_res(i) <= 420 and resol(i) <= 1.9)),
    dict(key="A2", n=4, tier="A", label="同源二聚体",
         note="生物学单元为二聚体(非晶体学堆积),单一序列 ×2,界面对称",
         src=["homodimer", "homodimer2", "homodimer_sym"],
         test=lambda i: (prot_ents(i) == 1 and nuc_ents(i) == 0
                         and not has_ab(i) and not is_tcr(i)
                         and bio_is(i, prot_inst=2, prot_ents_n=1)
                         and 230 <= total_res(i) <= 900 and resol(i) <= 2.2)),
    dict(key="A3", n=2, tier="A", label="同源四聚体及以上",
         note="生物学单元 ≥4 个相同亚基,考察对称寡聚体组装",
         src=["homotetramer", "oligo_Homo_4-mer", "oligo_Homo_6-mer", "oligo_Homo_8-mer", "oligo_Homo_3-mer"],
         test=lambda i: (prot_ents(i) == 1 and nuc_ents(i) == 0
                         and not has_ab(i) and not is_tcr(i)
                         and bio_is(i, prot_ents_n=1, min_inst=4)
                         and 320 <= total_res(i) <= 1600 and resol(i) <= 2.4)),
    dict(key="A4", n=4, tier="A", label="单体 + 内置辅因子",
         uniq=lambda i: sorted(c for c, *_ in i["ligands"])[0] if i["ligands"] else None,
         note="配体全部在 Server 下拉菜单内(ATP/HEM/FAD/GDP…),无需手动输入 CCD 代码",
         src=LIGP + ["mono_with_cofactor"],
         test=lambda i: (prot_ents(i) == 1 and nuc_ents(i) == 0
                         and len(i["chains"]) == 1
                         and builtin_only(i) and len(ligs(i)) <= 2
                         and not has_ab(i) and not is_tcr(i)
                         and 150 <= total_res(i) <= 700 and resol(i) <= 2.2)),
    dict(key="B1", n=4, tier="B", label="异源二聚体",
         note="两个不同蛋白的界面,依赖跨链共进化信号(paired MSA)",
         src=["hetero2", "peptide_complex"],
         test=lambda i: (prot_ents(i) == 2 and nuc_ents(i) == 0
                         and not has_ab(i) and not is_tcr(i)
                         and min(c["len"] for c in i["chains"]) >= 45
                         and 230 <= total_res(i) <= 1100 and resol(i) <= 2.6)),
    dict(key="B2", n=2, tier="B", label="异源三元复合物",
         note="三种不同亚基,界面数量增加,ipTM 更容易下降",
         src=["hetero3"],
         test=lambda i: (prot_ents(i) == 3 and nuc_ents(i) == 0
                         and not has_ab(i) and not is_tcr(i)
                         and 300 <= total_res(i) <= 1600 and resol(i) <= 3.0)),
    dict(key="B3", n=2, tier="B", label="蛋白 + 短肽",
         note="短肽(≤30 aa)结合槽。pTM 对短链系统性偏低(FAQ 明示),须改看 pLDDT/PAE",
         src=["peptide_complex", "hetero2"],
         test=lambda i: (prot_ents(i) == 2 and nuc_ents(i) == 0
                         and not has_ab(i) and not is_tcr(i)
                         and min(c["len"] for c in i["chains"]) <= 30
                         and 160 <= total_res(i) <= 800 and resol(i) <= 2.5)),
    dict(key="B4", n=4, tier="B", label="蛋白–DNA 复合物",
         note="仅标准 A/C/G/T;双链须分别输入两条互补链(或用 reverse-complement 选项)",
         src=["protein_dna"],
         test=lambda i: (nuc_ents(i) >= 1 and prot_ents(i) >= 1
                         and all(c["type"] != "RNA" for c in i["chains"])
                         and not has_ab(i) and not is_tcr(i)
                         and sum(c["len"] * c["copies"] for c in i["chains"]
                                 if c["type"] == "DNA") >= 8
                         and 120 <= total_res(i) <= 1200)),
    dict(key="B5", n=4, tier="B", label="蛋白–RNA 复合物",
         note="RNA 构象自由度大,是公认较难的一类",
         src=["protein_rna", "protein_rna2"],
         test=lambda i: (any(c["type"] == "RNA" for c in i["chains"])
                         and prot_ents(i) >= 1 and not has_ab(i) and not is_tcr(i)
                         and sum(c["len"] * c["copies"] for c in i["chains"]
                                 if c["type"] == "RNA") >= 6
                         and 120 <= total_res(i) <= 1300)),
    dict(key="C1", n=5, tier="C", label="抗体 / Fab–抗原复合物",
         note="CDR 环构象 + 表位定位,公认难点。务必跑多个 seed,并按 ipTM / chain_pair_iptm 排序",
         src=["antibody", "nanobody2"],
         test=lambda i: (is_ab_complex(i) and not is_nb(i) and not is_tcr(i)
                         and 400 <= total_res(i) <= 1300 and resol(i) <= 3.2)),
    dict(key="C2", n=2, tier="C", label="纳米抗体 / VHH–抗原复合物",
         note="单域抗体,链数少但表位仍难定位",
         src=["nanobody", "nanobody2", "antibody"],
         test=lambda i: (is_nb(i) and not is_tcr(i)
                         and 250 <= total_res(i) <= 900 and resol(i) <= 3.0)),
    dict(key="C3", n=3, tier="C", label="膜蛋白(通道 / 转运体 / 受体)",
         note="有 PDBTM/MemProtMD/mpstruc 跨膜注释。Server 不建模膜平面(FAQ 明示),跨膜螺旋排布与构象态易错",
         src=["memb_PDBTM", "memb_MemProtMD", "memb_mpstruc", "kw_gpcr",
              "kw_channel", "kw_transporter", "kw_abc", "membrane_kw", "transport_kw"],
         test=lambda i: (is_memb(i) and prot_ents(i) >= 1
                         and not has_ab(i) and not is_tcr(i)
                         and 250 <= total_res(i) <= 2200)),
    dict(key="C4", n=2, tier="C", label="TCR–pMHC 复合物",
         note="T 细胞受体识别 MHC-肽;5 条链、界面浅,是免疫结构预测的经典难题",
         src=["large_assembly", "signaling_gpcr", "hetero3", "antibody"],
         test=lambda i: (is_tcr(i) and prot_ents(i) >= 4 and nuc_ents(i) == 0
                         and 500 <= total_res(i) <= 1400)),
    dict(key="C5", n=2, tier="C", label="纯核酸结构(无蛋白)",
         note="没有蛋白 MSA 支撑,AlphaFold 3 在此类表现最弱",
         src=["nucleic_only", "nucleic_only2"],
         test=lambda i: (prot_ents(i) == 0 and nuc_ents(i) >= 1
                         and 25 <= total_res(i) <= 300)),
]
assert sum(s["n"] for s in SPECS) == 45

chosen = []
used_id, used_fam = set(), set()
used_topics = []


def topic_clash(t, limit=3):
    return any(len(t & u) >= limit for u in used_topics)


for sp in SPECS:
    cands = [i for i in pool(*sp["src"]) if assembly_sane(i) and sp["test"](i)]
    if sp["key"] in ("A1", "A2", "A3", "B1", "B2", "B3"):
        cands.sort(key=lambda i: (len(ligs(i)), resol(i), -total_res(i)))
    else:
        cands.sort(key=lambda i: (resol(i), -total_res(i)))
    got = 0
    used_uniq = set()
    for info in cands:
        if got >= sp["n"]:
            break
        if info["id"] in used_id:
            continue
        fk = famkey(info)
        if fk in used_fam:
            continue
        tp = topic(info)
        if topic_clash(tp):
            continue
        uf = sp.get("uniq")
        if uf is not None:
            uv = uf(info)
            if uv is None or uv in used_uniq:
                continue
            used_uniq.add(uv)
        used_id.add(info["id"])
        used_fam.add(fk)
        used_topics.append(tp)
        rec = dict(info)
        b = bio(info)
        rec["bio_assembly"] = b
        rec.update(group=sp["key"], group_label=sp["label"], tier=sp["tier"],
                   group_note=sp["note"])
        chosen.append(rec)
        got += 1
    status = "OK " if got == sp["n"] else "SHORT"
    print(f"{status} {sp['key']} {sp['label'][:22]:24s} pool={len(cands):4d} picked={got}/{sp['n']}")

print("\nTOTAL", len(chosen))
json.dump(chosen, open("chosen45.json", "w"), indent=1)
