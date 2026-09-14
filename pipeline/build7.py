"""Build the final 41-target pack into the standalone repo.

Applies every correction established during this project. Each is a real defect
found earlier, not a hypothetical:

  1. copy numbers from BIOLOGICAL ASSEMBLY 1, not the asymmetric unit
     (10VC's ASU held 6 copies of a 1:1 complex)
  2. ghost copies dropped -- an assembly member with <50% of its residues
     resolved is lattice occupancy, not a subunit (30ZS declared 3 trypsins,
     two of which had 3 residues each)
  3. protein input = FULL-LENGTH UniProt sequence (user's decision: reflects the
     real situation, where only the full sequence is in hand)
  4. PTM positions remapped to full-length numbering, then verified to land on
     the correct parent amino acid
  5. detergents / cryoprotectants / buffer salts excluded
  6. per-chain observation stats recorded, so students know how many predicted
     residues have no experimental counterpart to compare against
"""
import json
import os
from collections import Counter

REPO = "/Users/lilindu/alphafold-homework-2026"
BAN = ["班1", "班2", "班3", "班4"]

# alloc7.json still uses the old 班N keys, but targets.json and the job files use
# the class names 甲乙丙丁 and a per-class prefix. The rename is strictly a rename:
# no target changed classes. The mapping below was re-derived by aligning the ID
# SETS (not from memory) -- old 班4 is 甲班, the 11-student class, because the
# 13th antibody lands there.
RENAME = {"班4": "甲班", "班1": "乙班", "班2": "丙班", "班3": "丁班"}
PREFIX = {"甲班": "jia", "乙班": "yi", "丙班": "bing", "丁班": "ding"}
CLS_ORDER = ["甲班", "乙班", "丙班", "丁班"]

BUILTIN_LIG = {"ADP", "ATP", "AMP", "GTP", "GDP", "FAD", "NAD", "NAP", "NDP",
               "HEM", "HEC", "PLM", "OLA", "MYR", "CIT", "CLA", "CHL", "BCL", "BCB"}
BUILTIN_ION = {"MG", "ZN", "CL", "CA", "NA", "MN", "K", "FE", "CU", "CO"}
STD = set("ACDEFGHIKLMNPQRSTVWY")
TOKEN_LIMIT = 5000
MIN_OBS = 0.50

EXTRA_IGNORE = {
    "LMT", "LDA", "LMN", "DDQ", "C8E", "BOG", "OGA", "BNG", "MC3", "DMU", "UMQ",
    "F09", "LMU", "PCW", "P6G", "TWT", "SDS", "TRT", "TX4", "BO3", "B3P", "BEZ",
    "NHE", "CXS", "MPO", "EPE", "TAM", "BTB", "MRD", "SGM", "PGE", "1PE", "2PE",
    "7PE", "12P", "15P", "XPE", "PE4", "PE8", "PEG", "PG4", "PG0", "GOL", "EDO",
    "MPD", "DMS", "IPA", "MOH", "EOH", "ACT", "FMT", "SO4", "PO4", "NO3", "AZI",
    "SCN", "IMD", "TRS", "MES", "CAC", "OCT", "D10", "DD9", "UND", "HEZ", "PGO",
    "PDO", "BME", "DTT", "DTU", "DTV", "TCE", "BU3", "NH4", "UNX", "UNL", "PIN",
    "CLR", "Y01", "CHS", "PLC", "POV", "PEE", "PGT", "3PH", "DGA", "LPP", "OLC",
}
PTM_PARENT = {"SEP": "S", "TPO": "T", "PTR": "Y", "NEP": "H", "HIP": "H",
              "ALY": "K", "MLY": "K", "M3L": "K", "MLZ": "K", "KCR": "K",
              "YHA": "K", "CIR": "R", "2MR": "R", "AGM": "R", "HYP": "P",
              "HY3": "P", "LYZ": "K", "AHB": "N", "MCS": "C", "P1L": "C",
              "SNC": "C", "SNN": "N", "TRF": "W"}

raw = json.load(open("fin7_raw.json"))
alloc = json.load(open("alloc7.json"))
screen = json.load(open("screen.json"))
iface = json.load(open("iface.json"))
# 每题对应论文,来自 evidence/publications.json(随 repo 版本化,不是临时缓存)。
# 缺引用直接报错:作业的一个硬要求是「41 个结构全部已有论文或预印本」,
# 让它在生成阶段就卡住,而不是等到 handout 里出现空行才发现。
PUB = json.load(open(os.path.join(REPO, "evidence", "publications.json")))

# ---- re-apply the allocation fix (13th antibody -> 班4) ----
by_cat = {}
for b in BAN:
    for it in alloc[b]:
        by_cat.setdefault(it["category"], []).append(it)
fixed = {b: [] for b in BAN}
for key in ["hetero", "dna", "rna", "ligand", "ptm"]:
    for j, it in enumerate(by_cat[key][:4]):
        fixed[BAN[j]].append(it)
for j, it in enumerate(by_cat["memb"][:8]):
    fixed[BAN[j % 4]].append(it)
abs_ = by_cat["ab"]
for j, it in enumerate(abs_[:12]):
    fixed[BAN[j % 4]].append(it)
if len(abs_) > 12:
    fixed["班4"].append(abs_[12])


def parse_entities(eid):
    e = raw[eid]
    out = {}
    for pe in e.get("polymer_entities") or []:
        ent = pe["rcsb_id"].split("_")[-1]
        ep = pe.get("entity_poly") or {}
        construct = (ep.get("pdbx_seq_one_letter_code_can") or "").replace("\n", "").strip().upper()
        ptype = ep.get("rcsb_entity_polymer_type") or "?"
        desc = ((pe.get("rcsb_polymer_entity") or {}).get("pdbx_description") or "").strip()
        ups = [{"acc": u.get("rcsb_id"),
                "seq": ((u.get("rcsb_uniprot_protein") or {}).get("sequence") or "").upper()}
               for u in (pe.get("uniprots") or [])]
        ups = [u for u in ups if u["seq"]]
        regions, acc = [], None
        for al in pe.get("rcsb_polymer_entity_align") or []:
            if (al.get("reference_database_name") or "").upper().startswith("UNIPROT"):
                acc = al.get("reference_database_accession")
                for x in al.get("aligned_regions") or []:
                    if all(x.get(k) is not None for k in
                           ("entity_beg_seq_id", "ref_beg_seq_id", "length")):
                        regions.append((x["entity_beg_seq_id"], x["ref_beg_seq_id"], x["length"]))
        insts = []
        for pi in pe.get("polymer_entity_instances") or []:
            ch = ((pi.get("rcsb_polymer_entity_instance_container_identifiers") or {})
                  .get("auth_asym_id"))
            miss, ranges = 0, []
            for f in pi.get("rcsb_polymer_instance_feature") or []:
                if f.get("type") == "UNOBSERVED_RESIDUE_XYZ":
                    for fp in f.get("feature_positions") or []:
                        b0 = fp.get("beg_seq_id")
                        e0 = fp.get("end_seq_id") or b0
                        if b0:
                            miss += e0 - b0 + 1
                            ranges.append([b0, e0])
            obs = len(construct) - miss
            insts.append({"chain": ch, "observed": obs, "ranges": ranges,
                          "frac": (obs / len(construct)) if construct else 0.0})
        insts.sort(key=lambda x: -x["frac"])
        out[ent] = {"construct": construct, "type": ptype, "desc": desc,
                    "uniprots": ups, "regions": regions, "align_acc": acc,
                    "instances": insts}
    return out


def assembly_counts(eid):
    e = raw[eid]
    asms = e.get("assemblies") or []
    if not asms:
        return {}, ""
    a = asms[0]
    rep = (a.get("rcsb_assembly_info") or {}).get("polymer_entity_instance_count") or 0
    c = Counter(x["rcsb_polymer_entity_instance_container_identifiers"]["entity_id"]
                for x in (a.get("polymer_entity_instances") or [])
                if x.get("rcsb_polymer_entity_instance_container_identifiers"))
    listed = sum(c.values())
    per, note = dict(c), ""
    if listed and rep and listed != rep:
        if rep % listed == 0:
            f = rep // listed
            per = {k: v * f for k, v in per.items()}
            note = f"生物学装配按对称操作展开 ×{f}"
        else:
            note = f"PDB 报告 {rep} 条链、实例表列出 {listed} 条,拷贝数请自行核对"
    return per, note


def remap(pos, regions):
    for ebeg, rbeg, length in regions:
        if ebeg <= pos < ebeg + length:
            return rbeg + (pos - ebeg)
    return None


def homology_of(seq):
    v = screen.get(seq)
    if not v or v.get("bin") is None:
        return None
    return {"bin": v["bin"], "n": v["n"]}


records, issues = [], []
for b in BAN:
    for it in fixed[b]:
        r = it["rec"]
        eid = r["id"]
        ents = parse_entities(eid)
        counts, asm_note = assembly_counts(eid)

        chains, ghosts = [], []
        for ent, d in ents.items():
            n = counts.get(ent, 1)
            # A ghost copy is an instance that EXISTS in the deposited coordinates
            # but is mostly unresolved. Only then is a copy not real.
            #
            # Do NOT compare the assembly count against the number of ASU
            # instances: they live in different spaces. The assembly count is
            # symmetry-expanded, so a 1-instance ASU legitimately yields 2 copies
            # in a dimeric assembly. Comparing them made 30TL (14-3-3 sigma, a
            # genuine dimer) and 24UX (homodimer) look like they had ghosts, and
            # wrongly demoted both to monomers -- with an empty ghost list, which
            # is what exposed the error.
            insts = d["instances"]
            bad = [i for i in insts if i["frac"] < MIN_OBS]
            good = [i for i in insts if i["frac"] >= MIN_OBS]
            if bad and good:
                # scale the assembly count by the fraction of instances that are real
                n_new = max(1, round(n * len(good) / len(insts)))
                if n_new != n:
                    ghosts.append((ent, n, n_new,
                                   [(g["chain"], g["observed"], g["frac"]) for g in bad]))
                    n = n_new

            seq, source, upacc, regions = d["construct"], "construct", None, d["regions"]
            clen = len(d["construct"])
            if d["type"] == "Protein" and len(d["uniprots"]) == 1 and regions:
                full = d["uniprots"][0]["seq"]
                if not (set(full) - STD):
                    seq, source = full, "uniprot"
                    upacc = d["uniprots"][0]["acc"]
            best = d["instances"][0] if d["instances"] else None
            chains.append({
                "entity": ent, "type": d["type"], "desc": d["desc"],
                "seq": seq, "len": len(seq), "copies": n, "source": source,
                "uniprot": upacc, "regions": regions,
                "construct_len": clen,
                "observed": best["observed"] if best else None,
                "frac_observed": round(best["frac"], 3) if best else None,
                "unobserved_ranges": best["ranges"] if best else [],
                "homology": homology_of(d["construct"]),
                "reason_construct": (
                    None if source == "uniprot" else
                    "合成肽 / 人工设计蛋白,无天然全长序列" if not d["uniprots"] else
                    "融合构建体,无单一全长形式" if len(d["uniprots"]) > 1 else
                    "缺少 UniProt 比对信息"),
            })

        # PTM remap
        ptms = []
        for p in r["ptms"]:
            ch = next((c for c in chains if c["entity"] == p["entity"]), None)
            if not ch:
                continue
            pos = p["pos"]
            if ch["source"] == "uniprot":
                np_ = remap(pos, ch["regions"])
                if np_ is None or not (1 <= np_ <= ch["len"]):
                    issues.append((eid, f"修饰 {p['code']}@{pos} 无法映射到全长编号"))
                    continue
                pos = np_
            want = PTM_PARENT.get(p["code"])
            aa = ch["seq"][pos - 1]
            if want and aa != want:
                issues.append((eid, f"修饰 {p['code']} 落在 {aa}{pos},应为 {want}"))
                continue
            ptms.append({"entity": p["entity"], "code": p["code"], "pos": pos,
                         "pos_construct": p["pos"], "parent": aa})

        if eid not in PUB:
            raise SystemExit(f"evidence/publications.json 里没有 {eid} 的引用")

        ligs, ions, dropped = [], [], []
        for l in r["ligands"]:
            (dropped if l["code"] in EXTRA_IGNORE else ligs).append(l)
        for code, n in r["ions"]:
            (ions if code in BUILTIN_ION else dropped).append(
                {"code": code, "count": n} if code in BUILTIN_ION else {"code": code, "count": n})
        ions = [x for x in ions if isinstance(x, dict)]

        tokens = sum(c["len"] * c["copies"] for c in chains)
        tokens += sum((l["atoms"] or 0) * l["count"] for l in ligs)
        tokens += sum(i["count"] for i in ions)

        it_iface = iface.get(eid)
        ev = None
        if it_iface and it_iface.get("ptms"):
            bestp = max((p for p in it_iface["ptms"] if p["min_dist"] is not None),
                        key=lambda p: p["n_close"], default=None)
            if bestp:
                ev = {"min_dist": bestp["min_dist"], "n_close": bestp["n_close"],
                      "code": bestp["code"]}

        records.append({
            "class": RENAME[b], "category": it["category"], "label": it["label"],
            "note": it["note"], "id": eid, "title": r["title"],
            "publication": PUB[eid],
            "resolution": r["resolution"], "method": r["method"],
            "deposit": r["deposit"], "release": r["release"],
            "chains": chains, "ligands": ligs, "ions": ions,
            "dropped": [(d["code"], d["count"]) if isinstance(d, dict)
                        else (d["code"], d["count"]) for d in dropped],
            "omitted_pdb": list(r["omitted"]),
            "ptms": ptms, "ptm_interface": ev,
            "tokens": tokens, "assembly_note": asm_note,
            "ghosts": ghosts,
        })

CAT_ORDER = ["hetero", "dna", "rna", "ligand", "ptm", "memb", "ab"]
records.sort(key=lambda r: (CLS_ORDER.index(r["class"]),
                            CAT_ORDER.index(r["category"]), r["id"]))
_seen = {}
for r in records:
    _seen[r["class"]] = _seen.get(r["class"], 0) + 1
    r["no"] = _seen[r["class"]]          # 班内序号,不是全局连续编号
    r["prefix"] = PREFIX[r["class"]]

print("targets:", len(records))
for b in CLS_ORDER:
    print(f"  {b}: {sum(1 for r in records if r['class'] == b)} 题")
print("token range:", min(r["tokens"] for r in records), "-", max(r["tokens"] for r in records))
print("over limit:", sum(1 for r in records if r["tokens"] > TOKEN_LIMIT))
print("ghost-copy fixes:", sum(1 for r in records if r["ghosts"]))
for r in records:
    for ent, was, now, gl in r["ghosts"]:
        print(f"   {r['id']} e{ent}: copies {was}->{now}, 空拷贝 {gl}")
if issues:
    print("ISSUES:")
    for a, b2 in issues:
        print("  ", a, b2)

os.makedirs(REPO, exist_ok=True)
json.dump(records, open(os.path.join(REPO, "targets.json"), "w"),
          ensure_ascii=False, indent=1)
print("wrote", os.path.join(REPO, "targets.json"))
