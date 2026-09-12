"""Write the master assignment sheet for the 41-structure, four-班 pack."""
import json
import os

OUT = "/Users/lilindu/pdb-structures-for-homework/alphafold_homework_41"
recs = json.load(open(os.path.join(OUT, "targets41.json")))
bits = json.load(open("/tmp/afhw/bits41.json"))

CLASSES = ["班1", "班2", "班3", "班4"]
CAT_ORDER = ["hetero", "metal", "dna", "rna", "cofac", "ligand", "ptm", "memb", "ab"]
CAT_NUM = {k: i + 1 for i, k in enumerate(CAT_ORDER)}
CAT_LABEL = {
    "hetero": "异源蛋白复合体", "metal": "蛋白 + 金属离子", "dna": "蛋白 + DNA",
    "rna": "蛋白 + RNA", "cofac": "蛋白 + 内置辅因子",
    "ligand": "蛋白 + 非内置小分子配体", "ptm": "有翻译后修饰的蛋白",
    "memb": "膜蛋白", "ab": "抗原抗体复合物",
}

by_class = {c: [] for c in CLASSES}
for r in recs:
    by_class[r["class"]].append(r)
for c in CLASSES:
    by_class[c].sort(key=lambda r: r["no"])

L = []
A = L.append

A("# AlphaFold 3 结构预测作业 — 41 个目标(四个班)")
A("")
A("每位学生 1 个 PDB 结构:在 AlphaFold Server 上预测,再与实验结构比较。")
A("班1 / 班2 / 班3 各 10 题,班4 共 11 题。")
A("")
A("## 每班的九类构成")
A("")
A("| 类别 | 班1 | 班2 | 班3 | 班4 |")
A("|---|---|---|---|---|")
for k in CAT_ORDER:
    row = []
    for c in CLASSES:
        ids = [r["id"] for r in by_class[c] if r["category"] == k]
        row.append(" / ".join(ids) if ids else "—")
    A(f"| {CAT_NUM[k]}. {CAT_LABEL[k]} | " + " | ".join(row) + " |")
A(f"| **合计** | **{len(by_class['班1'])}** | **{len(by_class['班2'])}** "
  f"| **{len(by_class['班3'])}** | **{len(by_class['班4'])}** |")
A("")
A("第 9 类每班 2 个(班4 为 3 个,多出的一个是纳米抗体/VHH)。")
A("四个班同类的例子取自不同体系,难度不刻意对齐 —— 班际差异本身有比较价值。")
A("")

A("## 这批结构是怎么筛的")
A("")
A("从 PDB 检索出 3126 个候选,逐项用接口核验后 2082 个通过 Server 的全部限制,")
A("再按九类划分选出 41 个。四道硬门槛:")
A("")
A("1. **deposit 日期晚于 2025-02-03。** AlphaFold 3 训练数据截止 2021-09-30,")
A("   而 Server 允许把模板截止日期最晚调到 2025-02-03。选在这条线之后,学生无论")
A("   怎么设模板选项都不可能把答案本身当模板调进来。")
A("2. **配体必须能真正输进去。** Server 的 CCD 字典冻结在 version 2024_10_28,")
A("   之后新登记的代码会被判 invalid。这一条挡掉了 744 个结构 —— 是本次筛选中")
A("   最大的单一淘汰原因,也是「deposit 越新」与「配体可用」之间的真实冲突。")
A("3. **token ≤ 5000。** 按 FAQ 算法(蛋白 1/残基、核酸 1/碱基、配体 1/原子、")
A("   离子 1/个)逐个计算,本批最大 2060。")
A("4. **生物学装配 ≤ 12 条链。** 病毒衣壳一类的装配有上百条链,「正确答案」不是")
A("   学生会提交的那几条链。这一条挡掉 128 个。")
A("")
A("另外三件容易出错、已逐个核对的事:")
A("")
A("- **拷贝数按生物学装配算,不是按晶体学不对称单位。** 本批 22 处需要修正:")
A("  例如 10VC 不对称单位里每种链有 6 份(晶体堆积),生物学单元其实是 1:1。")
A("  照抄不对称单位会让学生多建几条链。")
A("- **去垢剂、结晶助剂不作为配体输入。** 膜蛋白里的 LMT(十二烷基麦芽糖苷)、")
A("  硼酸缓冲液 BO3 之类已剔除,每题都列出了「不要输入」清单。胆固醇(CLR)保留,")
A("  因为对膜蛋白而言它是真实的结构脂。")
A("- **翻译后修饰用 modifications 字段,不是配体。** 已核对每个修饰的位点残基")
A("  确实是对应的母体氨基酸(SEP 落在 Ser、TPO 落在 Thr、ALY 落在 Lys)。")
A("")
A("### 第 7 类(翻译后修饰)的额外筛选")
A("")
A("按您的要求,只选修饰对结构有明确影响的 —— 判据不是标题措辞,而是从坐标实测:")
A("修饰残基重原子到另一条链的最近距离 < 3 Å,且 4 Å 内接触原子 ≥ 20 个。")
A("四个例子的实测值:")
A("")
A("| PDB | 修饰 | 最近距离 | 4Å 内接触 | 体系 |")
A("|---|---|---|---|---|")
PTM_EVID = {
    "9S1S": ("SEP", "2.67 Å", 56, "CDC14A 磷酸酶识别磷酸化底物"),
    "9X8S": ("ALY", "2.66 Å", 57, "GAS41 YEATS 结构域读取乙酰化赖氨酸"),
    "30TL": ("TPO", "2.63 Å", 34, "14-3-3σ 结合 ERα 磷酸肽"),
    "9T9W": ("SEP", "2.61 Å", 32, "β-TrCP 识别双磷酸化 IκBα degron"),
}
for r in recs:
    if r["category"] == "ptm":
        e = PTM_EVID.get(r["id"])
        if e:
            A(f"| {r['id']} | {e[0]} | {e[1]} | {e[3+0-1] if False else e[2]} | {e[3]} |")
A("")
A("即修饰基团本身被对方蛋白多点抓住,是结合的决定因素,而不是恰好靠近。")
A("学生可以把修饰去掉再跑一次,直接看界面预测是否变差。")
A("")

A("## 分班清单")
A("")
for c in CLASSES:
    A(f"### {c}({len(by_class[c])} 题)")
    A("")
    A("| # | PDB | 类别 | 分辨率 | token | 组成 | 名称 |")
    A("|---|---|---|---|---|---|---|")
    for r in by_class[c]:
        comp = " + ".join(
            f"{{'Protein':'蛋白','DNA':'DNA','RNA':'RNA'}}[ch['type']]" for ch in [])
        parts = []
        for ch in r["chains"]:
            t = {"Protein": "蛋白", "DNA": "DNA", "RNA": "RNA"}[ch["type"]]
            cp = f"×{ch['copies']}" if ch["copies"] > 1 else ""
            parts.append(f"{t}{ch['len']}{cp}")
        comp = " + ".join(parts)
        res = f"{r['resolution']:.2f}" if r["resolution"] else "—"
        title = r["title"].replace("|", "/")
        if len(title) > 46:
            title = title[:45] + "…"
        A(f"| {r['no']} | [{r['id']}](https://www.rcsb.org/structure/{r['id']}) "
          f"| {CAT_NUM[r['category']]}. {CAT_LABEL[r['category']]} | {res} "
          f"| {r['tokens']} | {comp} | {title} |")
    A("")

A("## 每题的输入规格")
A("")
A("`job_files/` 下有对应 JSON,可用 Server 的 “Upload JSON” 直接导入,免去手工输序列;")
A("`sequences/` 是同样内容的 FASTA。下面是同一份规格的人读版本。")
A("")
for c in CLASSES:
    A(f"### {c}")
    A("")
    for r in by_class[c]:
        b = bits[r["id"]]
        res = f"{r['resolution']:.2f} Å" if r["resolution"] else "—"
        A(f"#### {r['no']}. {r['id']} — 第 {CAT_NUM[r['category']]} 类"
          f" {CAT_LABEL[r['category']]}")
        A("")
        A(r["title"])
        A("")
        A(f"- {r['method']} / {res} · deposit {r['deposit']} · release {r['release']}"
          f" · 估算 token {r['tokens']}")
        A(f"- 这一类考察什么:{r['note']}")
        A(f"- 评估重点:{b['metric']}")
        A("")
        A("输入:")
        for line in b["inputs"]:
            A(f"- {line}")
        if b["cautions"]:
            A("")
            A("注意:")
            for x in b["cautions"]:
                A(f"- {x}")
        A("")
        A(f"文件:`job_files/{r['no']:02d}_{r['id']}.json` · "
          f"`sequences/{r['no']:02d}_{r['id']}.fasta`")
        A("")

with open(os.path.join(OUT, "README.md"), "w") as f:
    f.write("\n".join(L))
print("wrote README.md:", len(L), "lines")
