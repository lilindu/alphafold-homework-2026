"""Emit the assignment pack: per-target job JSON, FASTA, and the master sheet."""
import json
import os

OUT = "/Users/lilindu/alphafold-homework-2026"
recs = json.load(open(os.path.join(OUT, "targets41.json")))

BUILTIN_LIG = {"ADP", "ATP", "AMP", "GTP", "GDP", "FAD", "NAD", "NAP", "NDP",
               "HEM", "HEC", "PLM", "OLA", "MYR", "CIT", "CLA", "CHL", "BCL", "BCB"}
ENTKEY = {"Protein": "proteinChain", "DNA": "dnaSequence", "RNA": "rnaSequence"}
COMP = {"A": "T", "T": "A", "G": "C", "C": "G"}
CLASSES = ["班1", "班2", "班3", "班4"]

CAT_ORDER = ["hetero", "metal", "dna", "rna", "cofac", "ligand", "ptm", "memb", "ab"]
CAT_NUM = {k: i + 1 for i, k in enumerate(CAT_ORDER)}


def revcomp_pairs(rec):
    dna = [c for c in rec["chains"] if c["type"] == "DNA"]
    out = []
    for i in range(len(dna)):
        for j in range(i + 1, len(dna)):
            a, b = dna[i]["seq"], dna[j]["seq"]
            if len(a) == len(b) and all(COMP.get(x) == y for x, y in zip(a, b[::-1])):
                out.append((dna[i]["entity"], dna[j]["entity"]))
    return out


def job_json(rec, no):
    seqs = []
    for c in rec["chains"]:
        et = ENTKEY.get(c["type"])
        if not et:
            continue
        d = {"sequence": c["seq"], "count": c["copies"]}
        if c["type"] == "Protein":
            mods = [p for p in rec["ptms"] if p["entity"] == c["entity"]]
            if mods:
                d["modifications"] = [
                    {"ptmType": "CCD_" + p["code"], "ptmPosition": p["pos"]}
                    for p in mods]
        seqs.append({et: d})
    for l in rec["ligands"]:
        code = ("CCD_" + l["code"]) if l["code"] in BUILTIN_LIG else l["code"]
        seqs.append({"ligand": {"ligand": code, "count": l["count"]}})
    for i in rec["ions"]:
        seqs.append({"ion": {"ion": i["code"], "count": i["count"]}})
    return [{"name": f"{no:02d}_{rec['id']}", "modelSeeds": [], "sequences": seqs,
             "dialect": "alphafoldserver", "version": 1}]


def inputs(rec):
    out = []
    for c in rec["chains"]:
        lab = {"Protein": "蛋白质链", "DNA": "DNA 链", "RNA": "RNA 链"}[c["type"]]
        cp = f" ×{c['copies']}" if c["copies"] > 1 else ""
        out.append(f"{lab}{cp}:{c['len']} {'aa' if c['type']=='Protein' else 'nt'}"
                   f" — {c['desc'] or '(未注明)'}")
        mods = [p for p in rec["ptms"] if p["entity"] == c["entity"]]
        for p in mods:
            out.append(f"　└ 翻译后修饰:{p['code']} 加在第 {p['pos']} 位")
    for l in rec["ligands"]:
        tag = "下拉菜单内" if l["code"] in BUILTIN_LIG else "需手输 CCD 代码"
        out.append(f"配体({tag}):{l['code']} ×{l['count']} — {l['name'][:56]}")
    for i in rec["ions"]:
        out.append(f"离子:{i['code']} ×{i['count']}")
    return out


def cautions(rec):
    cs = []
    manual = [l["code"] for l in rec["ligands"] if l["code"] not in BUILTIN_LIG]
    if manual:
        cs.append("在 request builder 里加 “CCD Code” 条目,手动输入:" + "、".join(manual))
    if rec["ptms"]:
        cs.append("修饰要用 proteinChain 的 modifications 字段(ptmType + ptmPosition),"
                  "不是当作配体输入")
    if revcomp_pairs(rec):
        cs.append("其中两条 DNA 链互为反向互补,是同一段双链的两条股,必须都输入"
                  "(或用 “+ Reverse complement”)")
    if rec["assembly_note"]:
        cs.append(rec["assembly_note"])
    if rec["omitted"]:
        cs.append("以下为结晶助剂/去垢剂/不支持的重原子,已剔除,不要输入:"
                  + "、".join(f"{c}×{n}" for c, n in rec["omitted"][:9]))
    short = [c for c in rec["chains"] if c["type"] == "Protein" and c["len"] < 20]
    if short:
        cs.append("含极短肽链(<20 aa):pTM 对短链系统性偏低(FAQ 明示),"
                  "评估以 pLDDT / PAE 为主")
    if rec["category"] == "ab":
        cs.append("建议跑 3–5 个不同 seed,按 ranking_score / chain_pair_iptm 选最优模型")
    return cs


def metric(rec):
    inst = sum(c["copies"] for c in rec["chains"])
    if inst == 1:
        return "单链:pLDDT + Cα RMSD / TM-score"
    return "多链:ipTM 与 chain_pair_iptm 看界面,pLDDT 看各链自身"


jobs_dir = os.path.join(OUT, "job_files")
seq_dir = os.path.join(OUT, "sequences")
os.makedirs(jobs_dir, exist_ok=True)
os.makedirs(seq_dir, exist_ok=True)

by_class = {c: [] for c in CLASSES}
for r in recs:
    by_class[r["class"]].append(r)
for c in CLASSES:
    by_class[c].sort(key=lambda r: (CAT_NUM[r["category"]], r["id"]))

no = 0
allj = []
for c in CLASSES:
    for r in by_class[c]:
        no += 1
        r["no"] = no
        j = job_json(r, no)
        allj.extend(j)
        stem = f"{no:02d}_{r['id']}"
        json.dump(j, open(os.path.join(jobs_dir, stem + ".json"), "w"), indent=1)
        with open(os.path.join(seq_dir, stem + ".fasta"), "w") as f:
            for ch in r["chains"]:
                f.write(f">{r['id']}_entity{ch['entity']}|{ch['type']}"
                        f"|copies={ch['copies']}|{ch['desc']}\n")
                for k in range(0, len(ch["seq"]), 60):
                    f.write(ch["seq"][k:k + 60] + "\n")
json.dump(allj, open(os.path.join(OUT, "all_41_jobs.json"), "w"), indent=1)
json.dump(recs, open(os.path.join(OUT, "targets41.json"), "w"),
          ensure_ascii=False, indent=1)
print("job files:", len(os.listdir(jobs_dir)), "| fasta:", len(os.listdir(seq_dir)))
json.dump({r["id"]: {"inputs": inputs(r), "cautions": cautions(r), "metric": metric(r)}
           for r in recs}, open("/tmp/afhw/bits41.json", "w"), ensure_ascii=False, indent=1)
print("wrote handout fragments")
