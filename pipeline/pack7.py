"""Emit job files, FASTAs and the master sheet into the standalone repo."""
import json
import os

REPO = "/Users/lilindu/alphafold_homework_41"
recs = json.load(open(os.path.join(REPO, "targets.json")))

BUILTIN_LIG = {"ADP", "ATP", "AMP", "GTP", "GDP", "FAD", "NAD", "NAP", "NDP",
               "HEM", "HEC", "PLM", "OLA", "MYR", "CIT", "CLA", "CHL", "BCL", "BCB"}
ENTKEY = {"Protein": "proteinChain", "DNA": "dnaSequence", "RNA": "rnaSequence"}
COMP = {"A": "T", "T": "A", "G": "C", "C": "G"}
BAN = ["班1", "班2", "班3", "班4"]
CAT_ORDER = ["hetero", "dna", "rna", "ligand", "ptm", "memb", "ab"]
CAT_NUM = {k: i + 1 for i, k in enumerate(CAT_ORDER)}
CAT_LABEL = {"hetero": "异源蛋白复合体", "dna": "蛋白 + DNA", "rna": "蛋白 + RNA",
             "ligand": "蛋白 + 小分子配体", "ptm": "有翻译后修饰的蛋白",
             "memb": "膜蛋白", "ab": "抗原抗体复合物"}


def hom_label(h):
    if not h:
        return "未筛查"
    b, n = h["bin"], h["n"]
    if b == 0:
        return "无同源体"
    pct = ">=95%" if b >= 0.95 else ">=60%" if b >= 0.60 else ">=30%"
    return f"{pct}({n} 个)"


def job_json(rec):
    seqs = []
    for c in rec["chains"]:
        et = ENTKEY.get(c["type"])
        if not et:
            continue
        d = {"sequence": c["seq"], "count": c["copies"]}
        if c["type"] == "Protein":
            mods = [p for p in rec["ptms"] if p["entity"] == c["entity"]]
            if mods:
                d["modifications"] = [{"ptmType": "CCD_" + p["code"],
                                       "ptmPosition": p["pos"]} for p in mods]
        seqs.append({et: d})
    for l in rec["ligands"]:
        code = ("CCD_" + l["code"]) if l["code"] in BUILTIN_LIG else l["code"]
        seqs.append({"ligand": {"ligand": code, "count": l["count"]}})
    for i in rec["ions"]:
        seqs.append({"ion": {"ion": i["code"], "count": i["count"]}})
    return [{"name": f"{rec['no']:02d}_{rec['id']}", "modelSeeds": [],
             "sequences": seqs, "dialect": "alphafoldserver", "version": 1}]


def revcomp_pairs(rec):
    dna = [c for c in rec["chains"] if c["type"] == "DNA"]
    out = []
    for i in range(len(dna)):
        for j in range(i + 1, len(dna)):
            a, b = dna[i]["seq"], dna[j]["seq"]
            if len(a) == len(b) and all(COMP.get(x) == y for x, y in zip(a, b[::-1])):
                out.append((dna[i]["entity"], dna[j]["entity"]))
    return out


def inputs(rec):
    out = []
    for c in rec["chains"]:
        lab = {"Protein": "蛋白质链", "DNA": "DNA 链", "RNA": "RNA 链"}[c["type"]]
        cp = f" ×{c['copies']}" if c["copies"] > 1 else ""
        unit = "aa" if c["type"] == "Protein" else "nt"
        line = f"{lab}{cp}:{c['len']} {unit} — {c['desc'] or '(未注明)'}"
        out.append(line)
        if c["type"] == "Protein":
            if c["source"] == "uniprot":
                obs = c.get("observed")
                extra = (f",实验结构里只解出 {obs} 个残基"
                         if obs is not None and obs < c["len"] else "")
                out.append(f"　└ 全长序列 {c['uniprot']}(构建体 {c['construct_len']} aa"
                           f"{extra})")
            else:
                out.append(f"　└ 用沉积序列:{c['reason_construct']}")
            h = c.get("homology")
            if h is not None:
                out.append(f"　└ 训练截止前同源体:{hom_label(h)}")
            for p in [x for x in rec["ptms"] if x["entity"] == c["entity"]]:
                out.append(f"　└ 翻译后修饰:{p['code']} 加在第 {p['pos']} 位"
                           f"({p['parent']},构建体编号 {p['pos_construct']})")
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
        cs.append("修饰用 proteinChain 的 modifications 字段(ptmType + ptmPosition),"
                  "位置已换算为全长编号")
    if revcomp_pairs(rec):
        cs.append("两条 DNA 链互为反向互补,是同一段双链的两股,必须都输入")
    if rec["assembly_note"]:
        cs.append(rec["assembly_note"])
    drop = rec.get("dropped") or []
    if drop or rec.get("omitted_pdb"):
        allc = list(drop) + list(rec.get("omitted_pdb") or [])
        seen, uniq = set(), []
        for code, n in allc:
            if code not in seen:
                seen.add(code)
                uniq.append(f"{code}×{n}")
        cs.append("去垢剂/结晶助剂/不支持的重原子已剔除,不要输入:" + "、".join(uniq[:10]))
    partial = [c for c in rec["chains"]
               if c["type"] == "Protein" and c.get("frac_observed") is not None
               and c["frac_observed"] < 0.95]
    if partial:
        d = "、".join(f"{c['desc'][:16]} {int(c['frac_observed']*100)}%" for c in partial[:3])
        cs.append(f"以下链在实验结构里并未全部解出({d}),比较时只对齐两者共有的残基")
    short = [c for c in rec["chains"] if c["type"] == "Protein" and c["len"] < 20]
    if short:
        cs.append("含极短肽链(<20 aa):pTM 对短链系统性偏低(FAQ 明示),"
                  "评估以 pLDDT / PAE 为主")
    if rec["category"] == "ab":
        cs.append("建议跑 3–5 个不同 seed,按 ranking_score / chain_pair_iptm 选最优模型")
    if rec["category"] == "memb":
        cs.append("Server 不知道膜平面,跨膜螺旋的相对排布是本题最可能出错的地方")
    return cs


def metric(rec):
    inst = sum(c["copies"] for c in rec["chains"])
    if inst == 1:
        return "单链:pLDDT + Cα RMSD / TM-score"
    return "多链:ipTM 与 chain_pair_iptm 看界面,pLDDT 看各链自身"


jobs = os.path.join(REPO, "job_files")
seqs = os.path.join(REPO, "sequences")
os.makedirs(jobs, exist_ok=True)
os.makedirs(seqs, exist_ok=True)

allj = []
for r in recs:
    j = job_json(r)
    allj.extend(j)
    stem = f"{r['no']:02d}_{r['id']}"
    json.dump(j, open(os.path.join(jobs, stem + ".json"), "w"), indent=1)
    with open(os.path.join(seqs, stem + ".fasta"), "w") as f:
        for c in r["chains"]:
            src = c["source"] if c["type"] == "Protein" else "pdb"
            f.write(f">{r['id']}_entity{c['entity']}|{c['type']}|copies={c['copies']}"
                    f"|source={src}|{c['desc']}\n")
            for k in range(0, len(c["seq"]), 60):
                f.write(c["seq"][k:k + 60] + "\n")
json.dump(allj, open(os.path.join(REPO, "all_41_jobs.json"), "w"), indent=1)
json.dump({r["id"]: {"inputs": inputs(r), "cautions": cautions(r), "metric": metric(r)}
           for r in recs}, open("/tmp/afhw/bits7.json", "w"),
          ensure_ascii=False, indent=1)
print("job files:", len(os.listdir(jobs)), "| fasta:", len(os.listdir(seqs)))
