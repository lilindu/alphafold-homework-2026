"""Write README.md for the standalone 41-target repo."""
import json
import os
import re

REPO = "/Users/lilindu/alphafold-homework-2026"
recs = json.load(open(os.path.join(REPO, "targets.json")))
bits = json.load(open("/tmp/afhw/bits7.json"))

BAN = ["甲班", "乙班", "丙班", "丁班"]
CAT_ORDER = ["hetero", "dna", "rna", "ligand", "ptm", "memb", "ab"]
CAT_NUM = {k: i + 1 for i, k in enumerate(CAT_ORDER)}
CAT_LABEL = {"hetero": "异源蛋白复合体", "dna": "蛋白 + DNA", "rna": "蛋白 + RNA",
             "ligand": "蛋白 + 小分子配体", "ptm": "有翻译后修饰的蛋白",
             "memb": "膜蛋白", "ab": "抗原抗体复合物"}
PER_BAN = {"hetero": 1, "dna": 1, "rna": 1, "ligand": 1, "ptm": 1, "memb": 2, "ab": 3}

# ---------------------------------------------------------------- 中文标点
# 本文档用全角句读(。、——),但 , : ; ( ) 一直是半角,两者混排。这里统一:
# 只在半角标点的**左右任一侧是中文**时换成全角,于是英文术语、URL、代码路径、
# 千分位(15,831)都不受影响。三类东西先摘出来原样保留:
#   `行内代码`  |  Markdown 链接语法 ](…)  |  裸 URL
# 括号按**配对**决定,避免出现全角开、半角闭。
HALF2FULL = {",": "，", ":": "：", ";": "；", "(": "（", ")": "）"}
CJK = re.compile(r"[　-〿一-鿿＀-￯]")
PROTECT = re.compile(r"`[^`]*`|\]\([^)]*\)|https?://\S+")


def _cn_once(s):
    pairs, stack = {}, []
    for j, ch in enumerate(s):
        if ch == "(":
            stack.append(j)
        elif ch == ")" and stack:
            pairs[stack.pop()] = j
    full = set()
    for i, j in pairs.items():
        if (CJK.search(s[i-1] if i else "") or CJK.search(s[j+1] if j + 1 < len(s) else "")
                or CJK.search(s[i+1:j])):
            full |= {i, j}
    out = []
    for j, ch in enumerate(s):
        if ch in HALF2FULL:
            prev = s[j-1] if j else ""
            nxt = s[j+1] if j + 1 < len(s) else ""
            if ch in "()":
                if j in full:
                    ch = HALF2FULL[ch]
            elif CJK.search(prev) or CJK.search(nxt):
                ch = HALF2FULL[ch]
        out.append(ch)
    return "".join(out)


def _cn(s):
    """迭代到稳定:括号转成全角后,紧跟其后的 , : ; 的邻居变了,需要重判一次
    (如 `…(非内置 CCD):6FN` 里的冒号,首轮看到的左邻还是半角 `)`)。"""
    for _ in range(4):
        t = _cn_once(s)
        if t == s:
            return s
        s = t
    return s


def cn_punct(line):
    out, last = [], 0
    for m in PROTECT.finditer(line):
        out.append(_cn(line[last:m.start()]))
        out.append(m.group(0))
        last = m.end()
    out.append(_cn(line[last:]))
    return "".join(out)



def hom_label(h):
    if not h:
        return "未筛查"
    b, n = h["bin"], h["n"]
    if b == 0:
        return "无"
    p = "≥95%" if b >= 0.95 else "≥60%" if b >= 0.60 else "≥30%"
    return f"{p} ({n})"


def worst_hom(r):
    hs = [c["homology"] for c in r["chains"] if c.get("homology")]
    if not hs:
        return None
    return max(hs, key=lambda h: (h["bin"], h["n"] or 0))


def hom_stats(cat):
    """(无同源的题数, 该类题数, 其余题各自的同源体数,升序)。

    正文里关于同源体的叙述一律从这个函数取数。原先这些数字是写死的,换题后
    就与下面的分班清单表自相矛盾(配体类换入 9P9O 后仍写着「2 题仅 1 个同源体」;
    PTM 类换入 9QNG 后仍写着「按数量取最低」)。
    """
    zero, tot, rest = 0, 0, []
    for r in recs:
        if r["category"] != cat:
            continue
        tot += 1
        h = worst_hom(r)
        if h is None:
            continue                     # 未筛查(DNA/RNA/抗体三类不做筛查)
        if h["bin"] == 0:
            zero += 1
        else:
            rest.append(h["n"] or 0)
    return zero, tot, sorted(rest)


def comp(r):
    """All four input kinds, in the order they appear in the job JSON."""
    parts = []
    for c in r["chains"]:
        t = {"Protein": "蛋白", "DNA": "DNA", "RNA": "RNA"}[c["type"]]
        cp = f"×{c['copies']}" if c["copies"] > 1 else ""
        parts.append(f"{t}{c['len']}{cp}")
    for p in r["ptms"]:
        parts.append(f"**{p['code']}**@{p['pos']}")
    for l in r["ligands"]:
        cp = f"×{l['count']}" if l["count"] > 1 else ""
        parts.append(f"配体 {l['code']}{cp}")
    for i in r["ions"]:
        cp = f"×{i['count']}" if i["count"] > 1 else ""
        parts.append(f"离子 {i['code']}{cp}")
    return " + ".join(parts)


by = {b: [r for r in recs if r["class"] == b] for b in BAN}
for b in BAN:
    by[b].sort(key=lambda r: r["no"])

L = []


def A(line):
    L.append(cn_punct(line))       # 统一中文标点(见上面的说明)


A("# AlphaFold 3 结构预测作业(41 题,四个班)")
A("")
A("每位学生一个 PDB 结构:在 [AlphaFold Server](https://alphafoldserver.com) 上做预测,")
A("再与实验解出的结构比较。甲班 11 题,乙班/丙班/丁班 各 10 题。")
A("")
_nj = sum(1 for r in recs if (r.get("publication") or {}).get("type") == "journal")
_np = len(recs) - _nj
A("## 七类构成")
A("")
A("| 类别 | 每班 | 甲班 | 乙班 | 丙班 | 丁班 |")
A("|---|---|---|---|---|---|")
for k in CAT_ORDER:
    row = []
    for b in BAN:
        ids = [r["id"] for r in by[b] if r["category"] == k]
        row.append(" ".join(ids) if ids else "—")
    n = PER_BAN[k]
    extra = "3(甲班 4)" if k == "ab" else str(n)
    A(f"| {CAT_NUM[k]}. {CAT_LABEL[k]} | {extra} | " + " | ".join(row) + " |")
A(f"| **合计** | | **{len(by['甲班'])}** | **{len(by['乙班'])}** "
  f"| **{len(by['丙班'])}** | **{len(by['丁班'])}** |")
A("")

A("## 这些结构满足以下要求")
A("")
A("1. **deposit 晚于 2025-02-03。** AlphaFold 3 训练数据截止 2021-09-30,而 Server 允许把")
A("   模板截止日期最晚调到 2025-02-03;选在这条线之后,学生无论怎么设模板选项,都不可能")
A("   把答案本身当模板调进来。")
A("2. **每条链 ≥ 4 个残基/碱基,且只含标准残基/碱基。** 更短的链、以及未知残基 `X`，")
A("   Server 都拒收。")
A("3. **配体都在冻结的 CCD 字典(version 2024_10_28)里。** 之后新登记的代码会被判 invalid；")
A("   这一条淘汰了 762 个候选(占被剔除总数 73%)—— 结构越新,配体越可能是新登记的。")
A("4. **token ≤ 5,000。** 按 FAQ 算法:蛋白 1/残基、核酸 1/碱基、配体 1/原子、离子 1/个。")
A("5. **生物学单元 ≤ 12 条链。** 上百条链的生物学单元里,「正确答案」不是学生会提交的那几条链。")
A(f"6. **全部已有正式论文或预印本**({_nj} 篇正式论文 + {_np} 篇预印本),每题都给出引用。")
A("")
A("从 PDB 中 deposit 晚于 2025-02-03 的 15,831 个实验结构出发,按类定向检索并合并去重得")
A("3,126 个,**用前 5 个要求过滤后剩 2,082 个**,再从中选出这 41 个 —— 第 6 条是对这 41 个")
A("逐一核对过的,不是筛选条件。")
A("")
A("**另外逐个核对过四件事:**")
A("")
A("- **拷贝数按生物学单元算,不是晶体学不对称单位。** 两者常不同:曾发现某条目不对称")
A("  单位里每种链 6 份(晶体堆积),生物学单元其实是 1:1。")
A("- **空拷贝已剔除。** 装配成员若只有极少数残基有坐标,那是晶格占位而非真实亚基")
A("  (曾发现某条目声明 3 份胰蛋白酶,其中两份各只有 3 个残基)。")
A("- **去垢剂、结晶助剂、缓冲盐不作为配体输入。** 每题都列出「不要输入」清单。")
A("- **翻译后修饰用 `modifications` 字段,不是配体。** 已核对每个修饰落在正确的母体")
A("  氨基酸上(SEP 在 Ser、TPO 在 Thr、ALY 在 Lys)。")
A("")
A("### 一个需要先知道的现实:多数结构在训练截止前已有同源体")
A("")
A("deposit 日期只保证「这个条目」不在训练集里,不保证同源蛋白不在。我们对四类做过同源")
A("筛查(检索训练截止 2021-09-30 前发布的 PDB 结构):")
A("")
A("| 类别 | 筛查结构数 | 所有链同源度 <60% |")
A("|---|---|---|")
A("| 1. 异源复合体 | 17 | 6 |")
A("| 4. 小分子配体 | 100 | 13 |")
A("| 5. 翻译后修饰 | 34 | **0** |")
A("| 6. 膜蛋白 | 74 | **23** |")
A("")
_zm, _tm, _ = hom_stats("memb")
_zh, _th, _ = hom_stats("hetero")
A(f"据此:**异源复合体 {_zh}/{_th} 题、膜蛋白 {_zm}/{_tm} 题选用训练截止前查不到同源体的"
  f"结构**,是真正的从头预测测试。")
_zl, _tl, _rl = hom_stats("ligand")
if _rl:
    A(f"配体类 {_zl}/{_tl} 题无同源,其余 {len(_rl)} 题的同源体数为 "
      + "、".join(f"{n} 个" for n in _rl) + "。")
A("")
_, _tp, _rp = hom_stats("ptm")
if _rp:
    A(f"**翻译后修饰类找不到干净的,这是结构性的**:被 Server 支持的修饰主要出现在 14-3-3、"
      f"β-TrCP 这类反复研究的识别模块上 —— 本作业 {len(_rp)} 道 PTM 题都已有同源体,"
      f"数量从 {_rp[0]} 个到 {_rp[-1]} 个,泄漏程度差很多。")
A("")
A("**DNA / RNA / 抗体三类未做筛查**,理由:核酸端没有「训练同源」概念;抗体骨架实测")
A("必然 ≥80%(免疫球蛋白高度保守),无法规避 —— 但抗原端往往查不到同源体,而表位识别")
A("正是这类要考的东西。")
A("")
A("**这对完成作业的要求:** 同源度高的题请分开回答两个问题 —— （1）单体折叠预测得怎样")
A("(若同源度 ≥95%,这一项接近查表)；（2）相互作用预测得怎样(界面/配体/修饰)。")
A("后者才是该题的实际考点,也是 AlphaFold 目前真正的弱项。")
A("")

A("## 输入用全长序列,但实验结构往往只是一部分")
A("")
A("PDB 里其实有三条不同的序列,常被混为一谈:")
A("")
A("| | 含义 | 性质 |")
A("|---|---|---|")
A("| A 全长 | UniProt 规范序列,完整基因产物 | — |")
A("| B 构建体 | 实验者实际放进管子里的东西 | 主动截短:为结晶只表达某个结构域,可能带标签 |")
A("| C 观测到的 | 真正有坐标的残基 | 被动损失:无序区、柔性末端看不到 |")
A("")
A("**本作业输入 A(全长)**,因为研究者在真实场景里手上只有全长序列。表达标签因此自动")
A("消失(UniProt 序列不含 His-tag)。代价是:预测模型里会有实验结构中根本不存在的部分。")
A("")
A("**所以比较时必须只对齐两者共有的残基。** 每题都标出了该链的全长长度、构建体长度和")
A("实际解出的残基数;若直接拿全长模型和实验结构算整体 RMSD,数字没有意义。")
A("")

A("## 分班清单")
A("")
for b in BAN:
    A(f"### {b}({len(by[b])} 题)")
    A("")
    A("| # | PDB | 类别 | 分辨率 | token | 组成 | 训练前同源 | 名称 |")
    A("|---|---|---|---|---|---|---|---|")
    for r in by[b]:
        res = f"{r['resolution']:.2f}" if r["resolution"] else "—"
        t = r["title"].replace("|", "/")
        A(f"| {r['no']} | [{r['id']}](https://www.rcsb.org/structure/{r['id']}) "
          f"| {CAT_NUM[r['category']]} | {res} | {r['tokens']} | {comp(r)} "
          f"| {hom_label(worst_hom(r))} | {t} |")
    A("")

A("## 逐题输入规格")
A("")
A("`job_files/` 下有对应 JSON,可用 Server 的 “Upload JSON” 直接导入;`sequences/` 是")
A("同样内容的 FASTA。下面是人读版本。")
A("")
A("**上传 JSON 后不需要手输任何东西** —— 序列、配体的 CCD 代码、离子、翻译后修饰都会")
A("自动填进 request builder(点 “Open draft” 即可核对)。下面标「非内置 CCD」的配体只是")
A("说它不在下拉菜单的 19 种内置辅因子里,不是说要自己敲。")
A("")
for b in BAN:
    A(f"### {b}")
    A("")
    for r in by[b]:
        d = bits[r["id"]]
        res = f"{r['resolution']:.2f} Å" if r["resolution"] else "—"
        A(f"#### {b} {r['no']}. {r['id']} — 第 {CAT_NUM[r['category']]} 类 "
          f"{CAT_LABEL[r['category']]}")
        A("")
        A(r["title"])
        A("")
        pub = r.get("publication") or {}
        if pub:
            kind = "预印本" if pub.get("type") == "preprint" else "论文"
            link = (f"https://doi.org/{pub['doi']}" if pub.get("doi")
                    else f"https://www.rcsb.org/structure/{r['id']}")
            A(f"- {kind}:[*{pub['journal']}* {pub['year']}]({link})")
        A(f"- {r['method']} / {res} · deposit {r['deposit']} · release {r['release']} "
          f"· 估算 token {r['tokens']}")
        A(f"- 这一类考察什么:{r['note']}")
        A(f"- 评估重点:{d['metric']}")
        if r.get("ptm_interface"):
            e = r["ptm_interface"]
            A(f"- 修饰介导互作的实测证据:{e['code']} 到对方链最近 {e['min_dist']} Å，"
              f"4 Å 内接触原子 {e['n_close']} 个")
        A("")
        A("输入:")
        for line in d["inputs"]:
            A(f"- {line}")
        if d["cautions"]:
            A("")
            A("注意:")
            for c in d["cautions"]:
                A(f"- {c}")
        A("")
        A(f"文件:`job_files/{r['prefix']}_{r['no']}_{r['id']}.json` · "
          f"`sequences/{r['prefix']}_{r['no']}_{r['id']}.fasta`")
        A("")

A("## 目录")
A("")
A("```")
A("README.md          本文件")
A("targets.json       41 题的结构化数据")
A("job_files/         41 个可直接 Upload JSON 导入的作业文件")
A("sequences/         41 个 FASTA")
A("all_41_jobs.json   一次性导入全部 41 题")
A("pipeline/          筛选流程的全部脚本(可复现)")
A("evidence/          筛选依据与出处:同源筛查、界面实测、装配拷贝数等原始结果,")
A("                   以及逐题的论文引用(publications.json)")
A("```")
A("")
A("## 数据来源")
A("")
A("结构数据来自 [RCSB PDB](https://www.rcsb.org)(CC0)。")
A("AlphaFold Server 的限制依据其 [FAQ](https://alphafoldserver.com/faq) 与")
A("[Release Updates](https://alphafoldserver.com/release-updates);")
A("非内置 CCD 配体输入自 2026-08-19 起开放。")
A("模型见 Abramson et al., *Nature* 630:493–500 (2024)。")

with open(os.path.join(REPO, "README.md"), "w") as f:
    f.write("\n".join(L))
print("wrote README.md:", len(L), "lines")
