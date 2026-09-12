"""Stage 13: write the master assignment sheet + per-student handouts."""
import json
import os

OUT = "/Users/lilindu/pdb-structures-for-homework/alphafold_homework"
recs = json.load(open(os.path.join(OUT, "targets.json")))
bits = json.load(open("/tmp/afhw/_handout_bits.json"))

TIER_NAME = {"A": "较易", "B": "中等", "C": "有挑战"}
TIER_DESC = {
    "A": "预期 AlphaFold 3 表现好。用来建立基线:看清高置信度预测长什么样。",
    "B": "预期中等。界面、核酸引入了真实难度,置信度指标开始分化。",
    "C": "预期困难。文献中公认的弱项,用来看清模型的边界在哪里。",
}


def fmt_chains(rec):
    parts = []
    for c in rec["chains"]:
        t = {"Protein": "蛋白", "DNA": "DNA", "RNA": "RNA"}.get(c["type"], c["type"])
        cp = f"×{c['copies']}" if c["copies"] > 1 else ""
        parts.append(f"{t}{c['len']}{cp}")
    return " + ".join(parts)


lines = []
A = lines.append

A("# AlphaFold 3 结构预测作业 — 45 个目标")
A("")
A("每位学生分配 1 个 PDB 结构:在 AlphaFold Server 上做预测,再与实验结构比较。")
A("")
A("## 这批结构是怎么挑的")
A("")
A("三条硬性条件,逐条用 RCSB 检索接口核验过,不是凭印象挑的:")
A("")
A("1. **全部在 2025-02-03 之后 deposit 到 PDB。** AlphaFold 3 的训练数据截止于")
A("   2021-09-30,而 Server 允许学生把模板截止日期最晚调到 2025-02-03。选在这条线")
A("   之后,意味着无论学生怎么设模板选项,都不可能把答案本身当模板调进来。")
A("2. **配体都能真正输进 Server。** Server 的 CCD 字典锁定在 version 2024_10_28,")
A("   2024-10-28 之后新登记的化学组分代码会被判为 invalid。候选池里 40% 的配体正是")
A("   因此被排除的 —— 这是「deposit 时间越新」与「配体可用」之间的真实冲突。")
A("   离子只支持内置 10 种,已逐个核对。")
A("3. **token 数在 5,000 上限内。** 按 FAQ 的算法(蛋白 1/残基、核酸 1/碱基、")
A("   配体 1/原子、离子 1/个)逐个算过,最大的一个是 1,682。")
A("")
A("另外做了三件容易被忽略的事:")
A("")
A("- **拷贝数按生物学装配算,不是按晶体学不对称单位。** 这两者常常不同 ——")
A("  比如有些条目不对称单位里有两条链,但生物学单元其实是单体,照抄就会多建一条链。")
A("- **剔除了超大装配。** 病毒衣壳一类的生物学单元有上百条链,「正确答案」不是学生")
A("  会提交的那几条链,不适合作题。")
A("- **去冗余。** 同一 UniProt、同一主题的近似结构不会重复占用名额。")
A("")
A("## 难度分层")
A("")
A("| 组 | 类型 | 数量 | 预期 |")
A("|---|---|---|---|")
seen = []
for r in recs:
    k = (r["tier"], r["group"], r["group_label"])
    if k not in seen:
        seen.append(k)
for tier, g, label in seen:
    n = sum(1 for r in recs if r["group"] == g)
    A(f"| {g} | {label} | {n} | {TIER_NAME[tier]} |")
A("")
for t in ("A", "B", "C"):
    n = sum(1 for r in recs if r["tier"] == t)
    A(f"- **{TIER_NAME[t]}({n} 个)**:{TIER_DESC[t]}")
A("")
A("## 目标清单")
A("")
A("| # | PDB | 组 | 难度 | 分辨率 | token | 组成 | 名称 |")
A("|---|---|---|---|---|---|---|---|")
for r in recs:
    res = f"{r['resolution']:.2f}" if r["resolution"] else "—"
    title = r["title"].replace("|", "/")
    if len(title) > 58:
        title = title[:57] + "…"
    A(f"| {r['no']} | [{r['id']}](https://www.rcsb.org/structure/{r['id']}) "
      f"| {r['group']} | {TIER_NAME[r['tier']]} | {res} | {r['tokens']} "
      f"| {fmt_chains(r)} | {title} |")
A("")

A("## 每个目标的输入清单")
A("")
A("`job_files/` 里有对应的 JSON,可直接用 Server 的 “Upload JSON” 导入,省去手工输序列。")
A("`sequences/` 是同样内容的 FASTA。下面列出的是同一份规格的人读版本。")
A("")
for r in recs:
    b = bits[r["id"]]
    res = f"{r['resolution']:.2f} Å" if r["resolution"] else "—"
    A(f"### {r['no']}. {r['id']} — {r['group_label']}({TIER_NAME[r['tier']]})")
    A("")
    A(f"{r['title']}")
    A("")
    A(f"- 方法 / 分辨率:{r['method']} / {res}")
    A(f"- deposit:{r['deposit']} · release:{r['release']} · 估算 token:{r['tokens']}")
    A(f"- 为什么归到这组:{r['group_note']}")
    A(f"- 评估重点:{b['metric']}")
    A("")
    A("输入:")
    for line in b["inputs"]:
        A(f"- {line}")
    if b["cautions"]:
        A("")
        A("注意:")
        for c in b["cautions"]:
            A(f"- {c}")
    A("")
    A(f"- 文件:`job_files/{r['no']:02d}_{r['id']}.json` · "
      f"`sequences/{r['no']:02d}_{r['id']}.fasta`")
    A("")

with open(os.path.join(OUT, "README.md"), "w") as f:
    f.write("\n".join(lines))
print("wrote README.md:", len(lines), "lines")
