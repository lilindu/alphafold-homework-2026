"""Stage 12: emit the assignment pack (job JSONs, FASTAs, master sheet, roster)."""
import json
import os
import re

OUT = "/Users/lilindu/pdb-structures-for-homework/alphafold_homework"
recs = json.load(open(os.path.join(OUT, "targets.json")))

BUILTIN_LIG = {"ADP", "ATP", "AMP", "GTP", "GDP", "FAD", "NAD", "NAP", "NDP",
               "HEM", "HEC", "PLM", "OLA", "MYR", "CIT", "CLA", "CHL", "BCL", "BCB"}
# free sugars must not be entered as ligands: on the Server a glycan is attached
# to a residue via the glycans field, so a loose NAG would be chemically wrong
SUGARS = {"NAG", "BMA", "MAN", "GLC", "BGC", "FUC", "GAL", "NDG", "A2G", "SIA"}

TIER_NAME = {"A": "较易", "B": "中等", "C": "有挑战"}
KEY = {"proteinChain": "sequence", "dnaSequence": "sequence", "rnaSequence": "sequence"}
ENTKEY = {"Protein": "proteinChain", "DNA": "dnaSequence", "RNA": "rnaSequence"}


def slug(s, n=40):
    s = re.sub(r"[^A-Za-z0-9]+", "_", s).strip("_")
    return s[:n]


def split_ligands(rec):
    """-> (enter_now, manual_ccd, skip) buckets."""
    enter, manual, skip = [], [], []
    for l in rec["ligands"]:
        if l["code"] in SUGARS:
            skip.append(l)
        elif l["code"] in BUILTIN_LIG:
            enter.append(l)
        else:
            manual.append(l)
    return enter, manual, skip


def job_json(rec):
    seqs = []
    for c in rec["chains"]:
        et = ENTKEY.get(c["type"])
        if not et:
            continue
        seqs.append({et: {"sequence": c["seq"], "count": c["copies"]}})
    enter, manual, _ = split_ligands(rec)
    for l in enter:
        seqs.append({"ligand": {"ligand": "CCD_" + l["code"], "count": l["count"]}})
    for l in manual:
        seqs.append({"ligand": {"ligand": l["code"], "count": l["count"]}})
    for i in rec["ions"]:
        seqs.append({"ion": {"ion": i["code"], "count": i["count"]}})
    return [{"name": f"{rec['no']:02d}_{rec['id']}", "modelSeeds": [],
             "sequences": seqs, "dialect": "alphafoldserver", "version": 1}]


def input_lines(rec):
    """Human-readable input spec for the handout."""
    out = []
    for c in rec["chains"]:
        label = {"Protein": "蛋白质链", "DNA": "DNA 链", "RNA": "RNA 链"}.get(c["type"], c["type"])
        cp = f" ×{c['copies']}" if c["copies"] > 1 else ""
        d = c["desc"] or "(未注明)"
        out.append(f"{label}{cp}:{c['len']} nt/aa — {d}")
    enter, manual, skip = split_ligands(rec)
    for l in enter:
        out.append(f"配体(下拉菜单内):{l['code']} ×{l['count']} — {l['name'][:60]}")
    for l in manual:
        out.append(f"配体(需手输 CCD 代码):{l['code']} ×{l['count']} — {l['name'][:60]}")
    for i in rec["ions"]:
        out.append(f"离子:{i['code']} ×{i['count']}")
    for l in skip:
        out.append(f"不要输入:{l['code']}(糖基,应作为糖链挂在残基上,不能当游离配体)")
    return out


def cautions(rec):
    """Per-target warnings a student needs before submitting."""
    cs = []
    enter, manual, skip = split_ligands(rec)
    if manual:
        cs.append("在 request builder 里加 “CCD Code” 条目,手动输入代码:"
                  + "、".join(l["code"] for l in manual))
    if skip:
        cs.append("PDB 中的糖基(" + "、".join(l["code"] for l in skip)
                  + ")不作为游离配体输入;若要建糖基化,须用 glycans 字段挂到指定残基上")
    if rec["ds_pairs"]:
        cs.append("这些 DNA 链两两互补,是同一段双链的两条股,必须都输入(或用 “+ Reverse complement”)")
    if rec["assembly_note"]:
        cs.append(rec["assembly_note"])
    if rec["omitted"]:
        codes = "、".join(f"{c}×{n}" for c, n in rec["omitted"][:8])
        cs.append(f"结晶助剂/非支持重原子已剔除,不要输入:{codes}")
    if rec["n_prot"] == 0:
        cs.append("纯核酸,没有蛋白 MSA 支撑,预期置信度低,这本身就是要观察的现象")
    if any(c["len"] < 20 for c in rec["chains"] if c["type"] == "Protein"):
        cs.append("含极短肽链(<20 aa):pTM 会系统性偏低(FAQ 明示),评估请以 pLDDT/PAE 为主")
    if rec["tier"] == "C":
        cs.append("属于挑战组:建议跑 3-5 个不同 seed,按 ranking_score / ipTM 选最优模型")
    return cs


def metric_hint(rec):
    """Guidance keyed on chain INSTANCES, not distinct entities: a homodimer has
    one entity but two chains, so it still has an interface to score."""
    inst = sum(c["copies"] for c in rec["chains"])
    homo = len(rec["chains"]) == 1 and inst > 1
    if rec["n_prot"] == 0:
        return "无蛋白:看 pLDDT 与整体 RMSD;pTM 对短核酸不可靠"
    if inst == 1:
        return "单链:主看 pLDDT + Cα RMSD / TM-score"
    if homo:
        return ("同源多聚体:先看 ipTM / chain_pair_iptm 判断亚基间排布,"
                "再看 pLDDT;单个亚基对得上但装配错位是这一组的典型失败模式")
    return "多链:主看 ipTM 与 chain_pair_iptm(界面),再看 pLDDT(各链自身)"


os.makedirs(OUT, exist_ok=True)
jobs_dir = os.path.join(OUT, "job_files")
seq_dir = os.path.join(OUT, "sequences")
os.makedirs(jobs_dir, exist_ok=True)
os.makedirs(seq_dir, exist_ok=True)

all_jobs = []
for rec in recs:
    j = job_json(rec)
    all_jobs.extend(j)
    stem = f"{rec['no']:02d}_{rec['id']}"
    with open(os.path.join(jobs_dir, stem + ".json"), "w") as f:
        json.dump(j, f, indent=1)
    with open(os.path.join(seq_dir, stem + ".fasta"), "w") as f:
        for c in rec["chains"]:
            hdr = f">{rec['id']}_entity{c['entity']}|{c['type']}|copies={c['copies']}|{c['desc']}"
            f.write(hdr + "\n")
            for k in range(0, len(c["seq"]), 60):
                f.write(c["seq"][k:k + 60] + "\n")

with open(os.path.join(OUT, "all_45_jobs.json"), "w") as f:
    json.dump(all_jobs, f, indent=1)

print("job files:", len(os.listdir(jobs_dir)))
print("fasta files:", len(os.listdir(seq_dir)))
json.dump({r["id"]: {"inputs": input_lines(r), "cautions": cautions(r),
                     "metric": metric_hint(r)} for r in recs},
          open("/tmp/afhw/_handout_bits.json", "w"),
          indent=1, ensure_ascii=False)
print("wrote handout fragments")


