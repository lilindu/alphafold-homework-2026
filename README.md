# AlphaFold 3 结构预测作业(41 题,四个班)

每位学生一个 PDB 结构:在 [AlphaFold Server](https://alphafoldserver.com) 上做预测,
再与实验解出的结构比较。甲班 11 题,乙班/丙班/丁班 各 10 题。

## 七类构成

| 类别 | 每班 | 甲班 | 乙班 | 丙班 | 丁班 |
|---|---|---|---|---|---|
| 1. 异源蛋白复合体 | 1 | 9V8I | 11KP | 9SCD | 9SDC |
| 2. 蛋白 + DNA | 1 | 36BT | 9NYI | 9U19 | 9Y1J |
| 3. 蛋白 + RNA | 1 | 9QS9 | 9SIT | 9VHE | 9TEL |
| 4. 蛋白 + 小分子配体 | 1 | 24UX | 9PV1 | 29LA | 21ZG |
| 5. 有翻译后修饰的蛋白 | 1 | 30TL | 9S1S | 9T9W | 9X8S |
| 6. 膜蛋白 | 2 | 9O0I 9VEB | 9MAU 9PS4 | 9M0S 9UET | 9M2H 9N93 |
| 7. 抗原抗体复合物 | 3(甲班 4) | 9NCY 9ULL 9ULM 9Y3I | 9NJY 9P4C 9S9E | 9NKZ 9VDY 9ZRO | 9IA3 9NW4 9Q1L |
| **合计** | | **11** | **10** | **10** | **10** |

四个班同类的结构取自不同体系,难度不刻意对齐 —— 班际差异本身有比较价值。

## 筛选过程

### 漏斗

| 步骤 | 数量 | 说明 |
|---|---|---|
| PDB 中 deposit 晚于 2025-02-03 的实验结构 | 15,831 | 全库基数 |
| 按类定向检索后合并去重 | 3,126 | |
| − 序列含未知/非标准残基或碱基 | −157 | 主体是未知残基 `X` |
| − 有链短于 4 个残基/碱基 | −30 | Server 要求每条链 ≥ 4 |
| − **配体不在冻结的 CCD 字典里** | **−762** | 最大的单一淘汰原因 |
| − 生物学装配超过 12 条链 | −93 | 病毒衣壳、纤维等 |
| − token 超过 5,000 | −2 | |
| **通过全部限制** | **2,082** | 从中选出 41 个 |

### 四道硬门槛

1. **deposit 晚于 2025-02-03。** AlphaFold 3 训练数据截止 2021-09-30,而 Server 允许
   把模板截止日期最晚调到 2025-02-03。选在这条线之后,学生无论怎么设模板选项,都不
   可能把答案本身当模板调进来。
2. **配体必须真能输进去。** Server 的 CCD 字典冻结在 version 2024_10_28,之后新登记
   的代码会被判 invalid。这一条淘汰了 762 个(占被剔除总数 73%):结构越新,配体越
   可能是新登记的 —— 这是「deposit 要新」与「配体要可用」之间的真实冲突。
3. **token ≤ 5,000。** 按 FAQ 算法:蛋白 1/残基、核酸 1/碱基、配体 1/原子、离子 1/个。
4. **生物学装配 ≤ 12 条链。** 装配上百条链时,「正确答案」不是学生会提交的那几条链。

### 逐个核对过的四件事

- **拷贝数按生物学装配算,不是晶体学不对称单位。** 两者常不同:曾发现某条目不对称
  单位里每种链 6 份(晶体堆积),生物学单元其实是 1:1。
- **空拷贝已剔除。** 装配成员若只有极少数残基有坐标,那是晶格占位而非真实亚基
  (曾发现某条目声明 3 份胰蛋白酶,其中两份各只有 3 个残基)。
- **去垢剂、结晶助剂、缓冲盐不作为配体输入。** 每题都列出「不要输入」清单。
- **翻译后修饰用 `modifications` 字段,不是配体。** 已核对每个修饰落在正确的母体
  氨基酸上(SEP 在 Ser、TPO 在 Thr、ALY 在 Lys)。

## 两个必须先讲清的问题

### 一、输入用全长序列,但实验结构往往只是一部分

PDB 里其实有三条不同的序列,常被混为一谈:

| | 含义 | 性质 |
|---|---|---|
| A 全长 | UniProt 规范序列,完整基因产物 | — |
| B 构建体 | 实验者实际放进管子里的东西 | 主动截短:为结晶只表达某个结构域,可能带标签 |
| C 观测到的 | 真正有坐标的残基 | 被动损失:无序区、柔性末端看不到 |

**本作业输入 A(全长)**,因为研究者在真实场景里手上只有全长序列。表达标签因此自动
消失(UniProt 序列不含 His-tag)。代价是:预测模型里会有实验结构中根本不存在的部分。

**所以比较时必须只对齐两者共有的残基。** 每题都标出了该链的全长长度、构建体长度和
实际解出的残基数;若直接拿全长模型和实验结构算整体 RMSD,数字没有意义。

### 二、多数结构在训练截止前已有同源体

deposit 日期只保证「这个条目」不在训练集里,不保证同源蛋白不在。我们对四类做了
同源筛查(检索训练截止 2021-09-30 前发布的 PDB 结构):

| 类别 | 筛查结构数 | 所有链同源度 <60% |
|---|---|---|
| 膜蛋白 | 74 | **23** |
| 小分子配体 | 100 | 13 |
| 异源复合体 | 17 | 6 |
| 翻译后修饰 | 34 | **0** |

据此:**膜蛋白 8 题和异源复合体 4 题全部选用训练窗口内查不到同源体的结构**,是真正
的从头预测测试。配体类 2 题无同源、2 题仅 1 个同源体。

**翻译后修饰类找不到干净的,这是结构性的**:被 Server 支持的修饰主要出现在 14-3-3、
β-TrCP 这类反复研究的识别模块上。已按同源体**数量**取最低(3 个 vs 278 个,同属
≥95% 档但泄漏程度差很多)。

**DNA / RNA / 抗体三类未做筛查**,理由:核酸端没有「训练同源」概念;抗体骨架实测
必然 ≥80%(免疫球蛋白高度保守),无法规避 —— 但抗原端往往查不到同源体,而表位识别
正是这类要考的东西。

**这对写报告的要求:** 同源度高的题请分开回答两个问题 —— (1) 单体折叠预测得怎样
(若同源度 ≥95%,这一项接近查表);(2) 相互作用预测得怎样(界面/配体/修饰)。
后者才是该题的实际考点,也是 AlphaFold 目前真正的弱项。

## 分班清单

### 甲班(11 题)

| # | PDB | 类别 | 分辨率 | token | 组成 | 训练前同源 | 名称 |
|---|---|---|---|---|---|---|---|
| 1 | [9V8I](https://www.rcsb.org/structure/9V8I) | 1 | 1.70 | 644 | 蛋白534 + 蛋白110 | 无 | Crystal structure for YxiD-YxxD toxin-immunity protein complex from Bacillus subtilis 6633. |
| 2 | [36BT](https://www.rcsb.org/structure/36BT) | 2 | 1.85 | 2942 | DNA38 + 蛋白1435 + 蛋白1435 + 配体 6FN + 离子 MG | 未筛查 | HIV-1 reverse transcriptase in complex with DNAddG Aptamer and unincorporated ISL-triphosphate |
| 3 | [9QS9](https://www.rcsb.org/structure/9QS9) | 3 | 1.44 | 336 | RNA6 + 蛋白162×2 + 离子 MG×6 | 未筛查 | Structure and mechanism of the broad spectrum CRISPR-associated ring nuclease Crn4 |
| 4 | [24UX](https://www.rcsb.org/structure/24UX) | 4 | 1.20 | 668 | 蛋白321×2 + 配体 SAH | 无 | Crystal structure of FPP-methyltransferase PcFPPMT from Pseudomonas chlororaphis O6 in complex with SAH |
| 5 | [30TL](https://www.rcsb.org/structure/30TL) | 5 | 1.20 | 509 | 蛋白5×2 + 蛋白248×2 + **TPO**@4 + 离子 CA + 离子 CL + 离子 MG | ≥95% (278) | 14-3-3sigma protein binding to ERalpha-weak peptide (AAA mutation) |
| 6 | [9O0I](https://www.rcsb.org/structure/9O0I) | 6 | 3.86 | 825 | 蛋白315×2 + 蛋白195 | 无 | Cryo-EM structure of Local KwaA-KwaB complex |
| 7 | [9VEB](https://www.rcsb.org/structure/9VEB) | 6 | 3.56 | 1086 | 蛋白543×2 | 无 | Cryo-EM structure of the aspartate:alanine antiporter AspT mutant L60C |
| 8 | [9NCY](https://www.rcsb.org/structure/9NCY) | 7 | 1.63 | 840 | 蛋白230 + 蛋白213 + 蛋白397 | 未筛查 | Fab1392 in complex with the C-terminal alpha-TSR domain of the P. falciparum circumsporozoite protein |
| 9 | [9ULL](https://www.rcsb.org/structure/9ULL) | 7 | 1.63 | 810 | 蛋白360 + 蛋白219 + 蛋白231 | 未筛查 | Mogamulizumab in complex with CCR4 N-terminus peptide (S14-S24) |
| 10 | [9ULM](https://www.rcsb.org/structure/9ULM) | 7 | 2.01 | 940 | 蛋白360 + 蛋白130 + 蛋白219 + 蛋白231 | 未筛查 | Mogamulizumab in complex with CCR4 N-terminus peptide (N2-C29) |
| 11 | [9Y3I](https://www.rcsb.org/structure/9Y3I) | 7 | 1.80 | 1120 | 蛋白301×2 + 蛋白132×2 + 蛋白127×2 | 未筛查 | Crystal Structure of PA14 Cif Bound to Nanobodies VHH108 and VHH219 |

### 乙班(10 题)

| # | PDB | 类别 | 分辨率 | token | 组成 | 训练前同源 | 名称 |
|---|---|---|---|---|---|---|---|
| 1 | [11KP](https://www.rcsb.org/structure/11KP) | 1 | 1.48 | 1182 | 蛋白837 + 蛋白345 | 无 | Crystal structure of the Caenorhabditis elegans telomeric POT-1-TEBP-1 complex interface |
| 2 | [9NYI](https://www.rcsb.org/structure/9NYI) | 2 | 1.98 | 1508 | DNA6×4 + 蛋白371×4 | 未筛查 | Structure of HalA in complex with oligodeoxyadenylate |
| 3 | [9SIT](https://www.rcsb.org/structure/9SIT) | 3 | 1.75 | 2942 | RNA60 + DNA46 + 蛋白335×6 + DNA46 + 蛋白181 + 蛋白255 + 蛋白344 | 未筛查 | Type I-F_HNH variant Cascade bound to dsDNA, HNH domain in inwards position |
| 4 | [9PV1](https://www.rcsb.org/structure/9PV1) | 4 | 1.20 | 674 | 蛋白310×2 + 配体 AKG×2 + 配体 BTN×2 + 配体 FE2 + 离子 CL | 无 | Biotin halogenase BtnX, anaerobic structure with Fe(II), biotin, alpha-ketoglutarate, chloride |
| 5 | [9S1S](https://www.rcsb.org/structure/9S1S) | 5 | 1.62 | 1212 | 蛋白6 + 蛋白603×2 + **SEP**@2 | ≥60% (3) | Crystal structure of C278S mutant of mouse CDC14A in complex with a model phosphopeptide |
| 6 | [9MAU](https://www.rcsb.org/structure/9MAU) | 6 | 2.87 | 563 | 蛋白563 | 无 | Cryo-EM structure of human OAT1 in the apo state |
| 7 | [9PS4](https://www.rcsb.org/structure/9PS4) | 6 | 3.29 | 1755 | 蛋白585×3 | 无 | Cryo-EM structure of NCLX without calcium (class 1) |
| 8 | [9NJY](https://www.rcsb.org/structure/9NJY) | 7 | 1.58 | 1428 | 蛋白989 + 蛋白225 + 蛋白214 | 未筛查 | Terminal two domains of ClfA002 with bound Fab of AZD7745 |
| 9 | [9P4C](https://www.rcsb.org/structure/9P4C) | 7 | 1.52 | 1049 | 蛋白622 + 蛋白211 + 蛋白216 | 未筛查 | Crystal structure of Mesothelin C-terminal peptide-RO4 Fab complex |
| 10 | [9S9E](https://www.rcsb.org/structure/9S9E) | 7 | 1.31 | 370 | 蛋白125×2 + 蛋白60×2 | 未筛查 | Co-crystal of broadly neutralizing biparatopic monomeric VHH in complex with cardiotoxin (P01468) naja pallida |

### 丙班(10 题)

| # | PDB | 类别 | 分辨率 | token | 组成 | 训练前同源 | 名称 |
|---|---|---|---|---|---|---|---|
| 1 | [9SCD](https://www.rcsb.org/structure/9SCD) | 1 | 1.52 | 1051 | 蛋白710 + 蛋白341 | 无 | Structure of S. pombe PNUTS (565 - 644) bound to Swd2.2 - crystal form 2 |
| 2 | [9U19](https://www.rcsb.org/structure/9U19) | 2 | 1.18 | 396 | 蛋白371 + DNA12 + DNA12 + 离子 MG | 未筛查 | Crystal structure of the NKX2.1 homeodomain in complex with a 12-bp DNA duplex containing a CACG motif variant |
| 3 | [9VHE](https://www.rcsb.org/structure/9VHE) | 3 | 2.50 | 2960 | DNA74 + RNA65 + 蛋白216 + 蛋白550×4 + 蛋白311 + 配体 ATP×3 + 离子 MG | 未筛查 | cryoEM structure of retron-Eco7 complex |
| 4 | [29LA](https://www.rcsb.org/structure/29LA) | 4 | 1.13 | 289 | 蛋白275 + 配体 4NC + 离子 K + 离子 CL + 离子 MG | ≥30% (1) | L-DOPA extradiol dioxygenase from Beta vulgaris in complex with 4-nitrocatechol |
| 5 | [9T9W](https://www.rcsb.org/structure/9T9W) | 5 | 1.16 | 922 | 蛋白317 + 蛋白605 + **SEP**@32 + **SEP**@36 | ≥95% (7) | Crystal structure of beta-TrCP bound by diphosphorylated I-kappa-B-alpha degron peptide |
| 6 | [9M0S](https://www.rcsb.org/structure/9M0S) | 6 | 3.50 | 600 | 蛋白549 + 配体 ACO | 无 | Acetyl-CoA-bound SLC33A1 in a cytoplasm-facing conformation |
| 7 | [9UET](https://www.rcsb.org/structure/9UET) | 6 | 3.68 | 814 | 蛋白406×2 + 离子 MG×2 | 无 | Cryo-EM structure of human choline-phosphotransferase 1 |
| 8 | [9NKZ](https://www.rcsb.org/structure/9NKZ) | 7 | 1.48 | 844 | 蛋白230 + 蛋白217 + 蛋白397 | 未筛查 | Crystal structure of Fab MAM01 in complex with NANP6 peptide from circumsporozoite protein |
| 9 | [9VDY](https://www.rcsb.org/structure/9VDY) | 7 | 2.28 | 1523 | 蛋白218 + 蛋白232 + 蛋白1073 | 未筛查 | hA5-6 Fab bound to SFTSV glycoprotein Gn |
| 10 | [9ZRO](https://www.rcsb.org/structure/9ZRO) | 7 | 1.40 | 3871 | 蛋白3433 + 蛋白214 + 蛋白224 | 未筛查 | Neutralizing W037 Fab antibody fragment in complex with West Nile Virus EDIII |

### 丁班(10 题)

| # | PDB | 类别 | 分辨率 | token | 组成 | 训练前同源 | 名称 |
|---|---|---|---|---|---|---|---|
| 1 | [9SDC](https://www.rcsb.org/structure/9SDC) | 1 | 1.70 | 644 | 蛋白77×4 + 蛋白84×4 | 无 | RelSI toxin-antitoxin complex |
| 2 | [9Y1J](https://www.rcsb.org/structure/9Y1J) | 2 | 1.55 | 398 | DNA5 + DNA10 + DNA16 + 蛋白335 + 配体 F2A + 离子 NA + 离子 MG | 未筛查 | S180R human DNA polymerase beta, Ternary complex dT:dAmpCpp |
| 3 | [9TEL](https://www.rcsb.org/structure/9TEL) | 3 | 1.44 | 728 | RNA10 + RNA10 + 蛋白674 + 配体 ADP + 配体 ALF + 离子 MG + 离子 ZN | 未筛查 | Structure of chicken LGP2 bound to 10-mer RNA mismatched duplex that mimics the influenza B virus vRNA promoter (panhandle) and to ADP-AlF4-Mg. |
| 4 | [21ZG](https://www.rcsb.org/structure/21ZG) | 4 | 1.40 | 350 | 蛋白338 + 配体 DHB + 离子 FE | ≥30% (1) | Crystal structure of the petrobactin-binding protein FatB from Bacillus cereus complexed with ferric siderophore mimic, Fe(3,4-DHB)2 |
| 5 | [9X8S](https://www.rcsb.org/structure/9X8S) | 5 | 1.70 | 922 | 蛋白7×2 + 蛋白227×4 + **ALY**@4 | ≥95% (7) | Crystal structure of the human GAS41 YEATS domain in complex with an acetylated YFV capsid peptide (K4ac) |
| 6 | [9M2H](https://www.rcsb.org/structure/9M2H) | 6 | 3.40 | 497 | 蛋白485 + 配体 3C4 | 无 | Structure of the auxin importer AUX1 in Arabidopsis thaliana in the CHPAA-bound state |
| 7 | [9N93](https://www.rcsb.org/structure/9N93) | 6 | 2.95 | 807 | 蛋白807 | 无 | Human TMEM63A mutant V53M lipid-open state |
| 8 | [9IA3](https://www.rcsb.org/structure/9IA3) | 7 | 1.11 | 634 | 蛋白191 + 蛋白228 + 蛋白215 | 未筛查 | Bc8.108 Fab bound to preS2 peptide |
| 9 | [9NW4](https://www.rcsb.org/structure/9NW4) | 7 | 1.82 | 827 | 蛋白229 + 蛋白220 + 蛋白378 | 未筛查 | Structure of CISV1 antibody bound to PvCSP repeat peptide |
| 10 | [9Q1L](https://www.rcsb.org/structure/9Q1L) | 7 | 1.56 | 1038 | 蛋白214 + 蛋白230 + 蛋白593 + 离子 CL | 未筛查 | Crystal structure of the walnut allergen Jug r 2 bound to the human-derived Fab 6D12 |

## 逐题输入规格

`job_files/` 下有对应 JSON,可用 Server 的 “Upload JSON” 直接导入;`sequences/` 是
同样内容的 FASTA。下面是人读版本。

**上传 JSON 后不需要手输任何东西** —— 序列、配体的 CCD 代码、离子、翻译后修饰都会
自动填进 request builder(点 “Open draft” 即可核对)。下面标「任意 CCD」的配体只是
说它不在下拉菜单的 19 种内置辅因子里,不是说要自己敲。

### 甲班

#### 甲班 1. 9V8I — 第 1 类 异源蛋白复合体

Crystal structure for YxiD-YxxD toxin-immunity protein complex from Bacillus subtilis 6633.

- X-RAY DIFFRACTION / 1.70 Å · deposit 2025-05-29 · release 2025-12-24 · 估算 token 644
- 这一类考察什么:两条以上不同蛋白链的复合体(不含抗体链、不含核酸、非膜蛋白)。考察跨链共进化信号能否把界面摆对
- 评估重点:多链:ipTM 与 chain_pair_iptm 看界面,pLDDT 看各链自身

输入:
- 蛋白质链:534 aa — YxiDCTD
- 　└ 全长序列 E0TXK9(构建体 141 aa,实验结构里只解出 107 个残基)
- 　└ 训练截止前同源体:无同源体
- 蛋白质链:110 aa — YxxD
- 　└ 全长序列 A0A9Q4H945(构建体 110 aa,实验结构里只解出 108 个残基)
- 　└ 训练截止前同源体:无同源体

注意:
- 以下链在实验结构里并未全部解出(YxiDCTD 75%),比较时只对齐两者共有的残基

文件:`job_files/jia_1_9V8I.json` · `sequences/jia_1_9V8I.fasta`

#### 甲班 2. 36BT — 第 2 类 蛋白 + DNA

HIV-1 reverse transcriptase in complex with DNAddG Aptamer and unincorporated ISL-triphosphate

- ELECTRON MICROSCOPY / 1.85 Å · deposit 2026-05-29 · release 2026-09-02 · 估算 token 2942
- 这一类考察什么:只含标准 A/C/G/T;双链须分别输入两条互补链
- 评估重点:多链:ipTM 与 chain_pair_iptm 看界面,pLDDT 看各链自身

输入:
- DNA 链:38 nt — DNAddG (38-MER)
- 蛋白质链:1435 aa — p51 RT
- 　└ 全长序列 P04585(构建体 455 aa,实验结构里只解出 375 个残基)
- 蛋白质链:1435 aa — Reverse transcriptase/ribonuclease H
- 　└ 全长序列 P04585(构建体 562 aa,实验结构里只解出 530 个残基)
- 配体(任意 CCD):6FN ×1 — 2'-deoxy-4'-ethynyl-2-fluoroadenosine 5'-(tetrahydrogen 
- 离子:MG ×1

注意:
- 以下配体不在下拉菜单里,而是任意 CCD 代码:6FN。上传 job JSON 会自动填入 “CCD Code” 栏;只有从零手工搭建输入时才需自己敲
- 以下链在实验结构里并未全部解出(p51 RT 82%、Reverse transcri 94%),比较时只对齐两者共有的残基

文件:`job_files/jia_2_36BT.json` · `sequences/jia_2_36BT.fasta`

#### 甲班 3. 9QS9 — 第 3 类 蛋白 + RNA

Structure and mechanism of the broad spectrum CRISPR-associated ring nuclease Crn4

- X-RAY DIFFRACTION / 1.44 Å · deposit 2025-04-04 · release 2025-12-03 · 估算 token 336
- 这一类考察什么:RNA 构象自由度大,是公认较难的一类
- 评估重点:多链:ipTM 与 chain_pair_iptm 看界面,pLDDT 看各链自身

输入:
- RNA 链:6 nt — Cyclic polyA
- 蛋白质链 ×2:162 aa — Crn4 H15A
- 　└ 全长序列 A0A1H0EUC2(构建体 149 aa,实验结构里只解出 125 个残基)
- 离子:MG ×6

注意:
- 以下链在实验结构里并未全部解出(Crn4 H15A 83%),比较时只对齐两者共有的残基

文件:`job_files/jia_3_9QS9.json` · `sequences/jia_3_9QS9.fasta`

#### 甲班 4. 24UX — 第 4 类 蛋白 + 小分子配体

Crystal structure of FPP-methyltransferase PcFPPMT from Pseudomonas chlororaphis O6 in complex with SAH

- X-RAY DIFFRACTION / 1.20 Å · deposit 2026-03-22 · release 2026-08-12 · 估算 token 668
- 这一类考察什么:配体不在 Server 的 19 种内置辅因子里,走任意 CCD 代码那条路;代码已核实在冻结的 CCD 2024_10_28 字典中
- 评估重点:多链:ipTM 与 chain_pair_iptm 看界面,pLDDT 看各链自身

输入:
- 蛋白质链 ×2:321 aa — Methyltransferase domain protein
- 　└ 全长序列 A0AB33WVX4(构建体 311 aa,实验结构里只解出 292 个残基)
- 　└ 训练截止前同源体:无同源体
- 配体(任意 CCD):SAH ×1 — S-ADENOSYL-L-HOMOCYSTEINE

注意:
- 以下配体不在下拉菜单里,而是任意 CCD 代码:SAH。上传 job JSON 会自动填入 “CCD Code” 栏;只有从零手工搭建输入时才需自己敲
- 生物学装配按对称操作展开 ×2
- 去垢剂/结晶助剂/不支持的重原子已剔除,不要输入:GOL×1、IMD×1
- 以下链在实验结构里并未全部解出(Methyltransferas 93%),比较时只对齐两者共有的残基

文件:`job_files/jia_4_24UX.json` · `sequences/jia_4_24UX.fasta`

#### 甲班 5. 30TL — 第 5 类 有翻译后修饰的蛋白

14-3-3sigma protein binding to ERalpha-weak peptide (AAA mutation)

- X-RAY DIFFRACTION / 1.20 Å · deposit 2026-05-13 · release 2026-08-26 · 估算 token 509
- 这一类考察什么:修饰残基经坐标实测介导蛋白-蛋白互作(最近距离 2.2–3.0 Å,4 Å 内接触原子 ≥ 20)。该类受体多为反复研究的识别模块,同源体无法避免,已在同源体数量上取最低
- 评估重点:多链:ipTM 与 chain_pair_iptm 看界面,pLDDT 看各链自身
- 修饰介导互作的实测证据:TPO 到对方链最近 2.63 Å,4 Å 内接触原子 34 个

输入:
- 蛋白质链 ×2:5 aa — ERa peptide pT594 weak mutant
- 　└ 用沉积序列:合成肽 / 人工设计蛋白,无天然全长序列
- 　└ 翻译后修饰:TPO 加在第 4 位(T,构建体编号 4)
- 蛋白质链 ×2:248 aa — 14-3-3 protein sigma
- 　└ 全长序列 P31947(构建体 236 aa,实验结构里只解出 236 个残基)
- 　└ 训练截止前同源体:>=95%(278 个)
- 离子:CA ×1
- 离子:CL ×1
- 离子:MG ×1

注意:
- 修饰用 proteinChain 的 modifications 字段(ptmType + ptmPosition),位置已换算为全长编号
- 生物学装配按对称操作展开 ×2
- 含极短肽链(ERa peptide pT594 weak mutant 5 aa):FAQ 明示 pTM 对短于 16 残基的链系统性偏低,该值接近 0 不代表预测失败。评估以 pLDDT / PAE 为主;要看界面就取 chain_pair_iptm 里「受体链 × 该肽链」那一格,不要用整体 ipTM

文件:`job_files/jia_5_30TL.json` · `sequences/jia_5_30TL.fasta`

#### 甲班 6. 9O0I — 第 6 类 膜蛋白

Cryo-EM structure of Local KwaA-KwaB complex

- ELECTRON MICROSCOPY / 3.86 Å · deposit 2025-04-02 · release 2025-08-06 · 估算 token 825
- 这一类考察什么:带 PDBTM / MemProtMD / mpstruc 跨膜注释。Server 不建模膜平面,跨膜螺旋排布与构象态最易错;本类优先选训练窗口内查不到同源体的
- 评估重点:多链:ipTM 与 chain_pair_iptm 看界面,pLDDT 看各链自身

输入:
- 蛋白质链 ×2:315 aa — Kiwa protein KwaB
- 　└ 全长序列 P0DW46(构建体 321 aa)
- 　└ 训练截止前同源体:无同源体
- 蛋白质链:195 aa — Kiwa protein KwaA
- 　└ 全长序列 P0DW45(构建体 203 aa)
- 　└ 训练截止前同源体:无同源体

注意:
- Server 不知道膜平面,跨膜螺旋的相对排布是本题最可能出错的地方

文件:`job_files/jia_6_9O0I.json` · `sequences/jia_6_9O0I.fasta`

#### 甲班 7. 9VEB — 第 6 类 膜蛋白

Cryo-EM structure of the aspartate:alanine antiporter AspT mutant L60C

- ELECTRON MICROSCOPY / 3.56 Å · deposit 2025-06-09 · release 2025-08-06 · 估算 token 1086
- 这一类考察什么:带 PDBTM / MemProtMD / mpstruc 跨膜注释。Server 不建模膜平面,跨膜螺旋排布与构象态最易错;本类优先选训练窗口内查不到同源体的
- 评估重点:多链:ipTM 与 chain_pair_iptm 看界面,pLDDT 看各链自身

输入:
- 蛋白质链 ×2:543 aa — Aspartate/alanine antiporter
- 　└ 全长序列 Q8L3K8(构建体 549 aa)
- 　└ 训练截止前同源体:无同源体

注意:
- Server 不知道膜平面,跨膜螺旋的相对排布是本题最可能出错的地方

文件:`job_files/jia_7_9VEB.json` · `sequences/jia_7_9VEB.fasta`

#### 甲班 8. 9NCY — 第 7 类 抗原抗体复合物

Fab1392 in complex with the C-terminal alpha-TSR domain of the P. falciparum circumsporozoite protein

- X-RAY DIFFRACTION / 1.63 Å · deposit 2025-02-17 · release 2026-03-11 · 估算 token 840
- 这一类考察什么:抗体 CDR 环构象 + 表位定位。抗体骨架本身高度保守(实测 ≥80%),无法规避,真正未知的是表位识别;建议跑多个 seed 按 ipTM 排序
- 评估重点:多链:ipTM 与 chain_pair_iptm 看界面,pLDDT 看各链自身

输入:
- 蛋白质链:230 aa — Fab1392 heavy chain
- 　└ 用沉积序列:合成肽 / 人工设计蛋白,无天然全长序列
- 蛋白质链:213 aa — Fab1392 light chain
- 　└ 用沉积序列:合成肽 / 人工设计蛋白,无天然全长序列
- 蛋白质链:397 aa — Circumsporozoite protein
- 　└ 全长序列 Q7K740(构建体 73 aa,实验结构里只解出 65 个残基)

注意:
- 去垢剂/结晶助剂/不支持的重原子已剔除,不要输入:SO4×2、GOL×3、PEG×3
- 以下链在实验结构里并未全部解出(Circumsporozoite 89%),比较时只对齐两者共有的残基
- 建议跑 3–5 个不同 seed,按 ranking_score / chain_pair_iptm 选最优模型

文件:`job_files/jia_8_9NCY.json` · `sequences/jia_8_9NCY.fasta`

#### 甲班 9. 9ULL — 第 7 类 抗原抗体复合物

Mogamulizumab in complex with CCR4 N-terminus peptide (S14-S24)

- X-RAY DIFFRACTION/SOLUTION SCATTERING / 1.63 Å · deposit 2025-04-20 · release 2026-04-22 · 估算 token 810
- 这一类考察什么:抗体 CDR 环构象 + 表位定位。抗体骨架本身高度保守(实测 ≥80%),无法规避,真正未知的是表位识别;建议跑多个 seed 按 ipTM 排序
- 评估重点:多链:ipTM 与 chain_pair_iptm 看界面,pLDDT 看各链自身

输入:
- 蛋白质链:360 aa — C-C chemokine receptor type 4
- 　└ 全长序列 P51679(构建体 11 aa,实验结构里只解出 10 个残基)
- 蛋白质链:219 aa — light chain
- 　└ 用沉积序列:合成肽 / 人工设计蛋白,无天然全长序列
- 蛋白质链:231 aa — heavy chain
- 　└ 用沉积序列:合成肽 / 人工设计蛋白,无天然全长序列

注意:
- 以下链在实验结构里并未全部解出(C-C chemokine re 90%、heavy chain 93%),比较时只对齐两者共有的残基
- 建议跑 3–5 个不同 seed,按 ranking_score / chain_pair_iptm 选最优模型

文件:`job_files/jia_9_9ULL.json` · `sequences/jia_9_9ULL.fasta`

#### 甲班 10. 9ULM — 第 7 类 抗原抗体复合物

Mogamulizumab in complex with CCR4 N-terminus peptide (N2-C29)

- X-RAY DIFFRACTION / 2.01 Å · deposit 2025-04-20 · release 2026-04-22 · 估算 token 940
- 这一类考察什么:抗体 CDR 环构象 + 表位定位。抗体骨架本身高度保守(实测 ≥80%),无法规避,真正未知的是表位识别;建议跑多个 seed 按 ipTM 排序
- 评估重点:多链:ipTM 与 chain_pair_iptm 看界面,pLDDT 看各链自身

输入:
- 蛋白质链:360 aa — C-C chemokine receptor type 4
- 　└ 全长序列 P51679(构建体 28 aa,实验结构里只解出 11 个残基)
- 蛋白质链:130 aa — anti-kappa VHH
- 　└ 用沉积序列:合成肽 / 人工设计蛋白,无天然全长序列
- 蛋白质链:219 aa — light chain
- 　└ 用沉积序列:合成肽 / 人工设计蛋白,无天然全长序列
- 蛋白质链:231 aa — heavy chain
- 　└ 用沉积序列:合成肽 / 人工设计蛋白,无天然全长序列

注意:
- 以下链在实验结构里并未全部解出(C-C chemokine re 39%、anti-kappa VHH 92%、heavy chain 92%),比较时只对齐两者共有的残基
- 建议跑 3–5 个不同 seed,按 ranking_score / chain_pair_iptm 选最优模型

文件:`job_files/jia_10_9ULM.json` · `sequences/jia_10_9ULM.fasta`

#### 甲班 11. 9Y3I — 第 7 类 抗原抗体复合物

Crystal Structure of PA14 Cif Bound to Nanobodies VHH108 and VHH219

- X-RAY DIFFRACTION / 1.80 Å · deposit 2025-09-02 · release 2026-09-09 · 估算 token 1120
- 这一类考察什么:抗体 CDR 环构象 + 表位定位。抗体骨架本身高度保守(实测 ≥80%),无法规避,真正未知的是表位识别;建议跑多个 seed 按 ipTM 排序
- 评估重点:多链:ipTM 与 chain_pair_iptm 看界面,pLDDT 看各链自身

输入:
- 蛋白质链 ×2:301 aa — CFTR inhibitory factor
- 　└ 全长序列 A0A0M3KL26(构建体 295 aa,实验结构里只解出 294 个残基)
- 蛋白质链 ×2:132 aa — Nanobody VHH108
- 　└ 用沉积序列:合成肽 / 人工设计蛋白,无天然全长序列
- 蛋白质链 ×2:127 aa — Nanobody VHH219
- 　└ 用沉积序列:合成肽 / 人工设计蛋白,无天然全长序列

注意:
- 去垢剂/结晶助剂/不支持的重原子已剔除,不要输入:FMT×2
- 以下链在实验结构里并未全部解出(Nanobody VHH108 93%、Nanobody VHH219 93%),比较时只对齐两者共有的残基
- 建议跑 3–5 个不同 seed,按 ranking_score / chain_pair_iptm 选最优模型

文件:`job_files/jia_11_9Y3I.json` · `sequences/jia_11_9Y3I.fasta`

### 乙班

#### 乙班 1. 11KP — 第 1 类 异源蛋白复合体

Crystal structure of the Caenorhabditis elegans telomeric POT-1-TEBP-1 complex interface

- X-RAY DIFFRACTION / 1.48 Å · deposit 2026-03-02 · release 2026-08-12 · 估算 token 1182
- 这一类考察什么:两条以上不同蛋白链的复合体(不含抗体链、不含核酸、非膜蛋白)。考察跨链共进化信号能否把界面摆对
- 评估重点:多链:ipTM 与 chain_pair_iptm 看界面,pLDDT 看各链自身

输入:
- 蛋白质链:837 aa — Double-strand telomeric DNA-binding proteins 1
- 　└ 全长序列 O62329(构建体 111 aa,实验结构里只解出 110 个残基)
- 　└ 训练截止前同源体:无同源体
- 蛋白质链:345 aa — Protection of telomeres homolog 1
- 　└ 全长序列 A0A7R7JK61(构建体 172 aa,实验结构里只解出 165 个残基)
- 　└ 训练截止前同源体:无同源体

文件:`job_files/yi_1_11KP.json` · `sequences/yi_1_11KP.fasta`

#### 乙班 2. 9NYI — 第 2 类 蛋白 + DNA

Structure of HalA in complex with oligodeoxyadenylate

- ELECTRON MICROSCOPY / 1.98 Å · deposit 2025-03-27 · release 2025-05-07 · 估算 token 1508
- 这一类考察什么:只含标准 A/C/G/T;双链须分别输入两条互补链
- 评估重点:多链:ipTM 与 chain_pair_iptm 看界面,pLDDT 看各链自身

输入:
- DNA 链 ×4:6 nt — DNA (5'-D(*AP*AP*AP*AP*AP*A)-3')
- 蛋白质链 ×4:371 aa — Structure of HalA in complex with oligodeoxyadenylate
- 　└ 用沉积序列:合成肽 / 人工设计蛋白,无天然全长序列

注意:
- 以下链在实验结构里并未全部解出(Structure of Hal 89%),比较时只对齐两者共有的残基

文件:`job_files/yi_2_9NYI.json` · `sequences/yi_2_9NYI.fasta`

#### 乙班 3. 9SIT — 第 3 类 蛋白 + RNA

Type I-F_HNH variant Cascade bound to dsDNA, HNH domain in inwards position

- ELECTRON MICROSCOPY / 1.75 Å · deposit 2025-08-29 · release 2026-02-18 · 估算 token 2942
- 这一类考察什么:RNA 构象自由度大,是公认较难的一类
- 评估重点:多链:ipTM 与 chain_pair_iptm 看界面,pLDDT 看各链自身

输入:
- RNA 链:60 nt — crRNA
- DNA 链:46 nt — Target strand
- 蛋白质链 ×6:335 aa — Cas7f
- 　└ 全长序列 A0AAX7FM28(构建体 335 aa,实验结构里只解出 324 个残基)
- DNA 链:46 nt — Non-target strand
- 蛋白质链:181 aa — Cas6f
- 　└ 全长序列 A0AAX7FM27(构建体 181 aa,实验结构里只解出 178 个残基)
- 蛋白质链:255 aa — Cas5f
- 　└ 全长序列 A0AAX7FM22(构建体 255 aa,实验结构里只解出 252 个残基)
- 蛋白质链:344 aa — Cas8f fusion with HNH
- 　└ 全长序列 A0AAX7FM29(构建体 344 aa,实验结构里只解出 340 个残基)

注意:
- 两条 DNA 链互为反向互补,是同一段双链的两股,必须都输入

文件:`job_files/yi_3_9SIT.json` · `sequences/yi_3_9SIT.fasta`

#### 乙班 4. 9PV1 — 第 4 类 蛋白 + 小分子配体

Biotin halogenase BtnX, anaerobic structure with Fe(II), biotin, alpha-ketoglutarate, chloride

- X-RAY DIFFRACTION / 1.20 Å · deposit 2025-07-31 · release 2026-06-24 · 估算 token 674
- 这一类考察什么:配体不在 Server 的 19 种内置辅因子里,走任意 CCD 代码那条路;代码已核实在冻结的 CCD 2024_10_28 字典中
- 评估重点:多链:ipTM 与 chain_pair_iptm 看界面,pLDDT 看各链自身

输入:
- 蛋白质链 ×2:310 aa — Biotin halogenase BtnX
- 　└ 全长序列 A8LT50(构建体 312 aa)
- 　└ 训练截止前同源体:无同源体
- 配体(任意 CCD):AKG ×2 — 2-OXOGLUTARIC ACID
- 配体(任意 CCD):BTN ×2 — BIOTIN
- 配体(任意 CCD):FE2 ×1 — FE (II) ION
- 离子:CL ×1

注意:
- 以下配体不在下拉菜单里,而是任意 CCD 代码:AKG、BTN、FE2。上传 job JSON 会自动填入 “CCD Code” 栏;只有从零手工搭建输入时才需自己敲
- 去垢剂/结晶助剂/不支持的重原子已剔除,不要输入:EDO×2、PEG×1、SO4×1

文件:`job_files/yi_4_9PV1.json` · `sequences/yi_4_9PV1.fasta`

#### 乙班 5. 9S1S — 第 5 类 有翻译后修饰的蛋白

Crystal structure of C278S mutant of mouse CDC14A in complex with a model phosphopeptide

- X-RAY DIFFRACTION / 1.62 Å · deposit 2025-07-21 · release 2025-12-10 · 估算 token 1212
- 这一类考察什么:修饰残基经坐标实测介导蛋白-蛋白互作(最近距离 2.2–3.0 Å,4 Å 内接触原子 ≥ 20)。该类受体多为反复研究的识别模块,同源体无法避免,已在同源体数量上取最低
- 评估重点:多链:ipTM 与 chain_pair_iptm 看界面,pLDDT 看各链自身
- 修饰介导互作的实测证据:SEP 到对方链最近 2.67 Å,4 Å 内接触原子 56 个

输入:
- 蛋白质链:6 aa — Model phosphohexapeptide
- 　└ 用沉积序列:合成肽 / 人工设计蛋白,无天然全长序列
- 　└ 翻译后修饰:SEP 加在第 2 位(S,构建体编号 2)
- 蛋白质链 ×2:603 aa — Dual specificity protein phosphatase CDC14A
- 　└ 全长序列 Q6GQT0(构建体 344 aa,实验结构里只解出 339 个残基)
- 　└ 训练截止前同源体:>=60%(3 个)

注意:
- 修饰用 proteinChain 的 modifications 字段(ptmType + ptmPosition),位置已换算为全长编号
- 去垢剂/结晶助剂/不支持的重原子已剔除,不要输入:GOL×2、SO4×2
- 以下链在实验结构里并未全部解出(Model phosphohex 66%),比较时只对齐两者共有的残基
- 含极短肽链(Model phosphohexapeptide 6 aa):FAQ 明示 pTM 对短于 16 残基的链系统性偏低,该值接近 0 不代表预测失败。评估以 pLDDT / PAE 为主;要看界面就取 chain_pair_iptm 里「受体链 × 该肽链」那一格,不要用整体 ipTM

文件:`job_files/yi_5_9S1S.json` · `sequences/yi_5_9S1S.fasta`

#### 乙班 6. 9MAU — 第 6 类 膜蛋白

Cryo-EM structure of human OAT1 in the apo state

- ELECTRON MICROSCOPY / 2.87 Å · deposit 2025-03-14 · release 2025-06-18 · 估算 token 563
- 这一类考察什么:带 PDBTM / MemProtMD / mpstruc 跨膜注释。Server 不建模膜平面,跨膜螺旋排布与构象态最易错;本类优先选训练窗口内查不到同源体的
- 评估重点:单链:pLDDT + Cα RMSD / TM-score

输入:
- 蛋白质链:563 aa — Isoform 2 of Solute carrier family 22 member 6
- 　└ 全长序列 Q4U2R8(构建体 572 aa,实验结构里只解出 506 个残基)
- 　└ 训练截止前同源体:无同源体

注意:
- 以下链在实验结构里并未全部解出(Isoform 2 of Sol 88%),比较时只对齐两者共有的残基
- Server 不知道膜平面,跨膜螺旋的相对排布是本题最可能出错的地方

文件:`job_files/yi_6_9MAU.json` · `sequences/yi_6_9MAU.fasta`

#### 乙班 7. 9PS4 — 第 6 类 膜蛋白

Cryo-EM structure of NCLX without calcium (class 1)

- ELECTRON MICROSCOPY / 3.29 Å · deposit 2025-07-24 · release 2025-09-10 · 估算 token 1755
- 这一类考察什么:带 PDBTM / MemProtMD / mpstruc 跨膜注释。Server 不建模膜平面,跨膜螺旋排布与构象态最易错;本类优先选训练窗口内查不到同源体的
- 评估重点:多链:ipTM 与 chain_pair_iptm 看界面,pLDDT 看各链自身

输入:
- 蛋白质链 ×3:585 aa — NCLX
- 　└ 全长序列 Q6AXS0(构建体 585 aa,实验结构里只解出 490 个残基)
- 　└ 训练截止前同源体:无同源体

注意:
- 以下链在实验结构里并未全部解出(NCLX 83%),比较时只对齐两者共有的残基
- Server 不知道膜平面,跨膜螺旋的相对排布是本题最可能出错的地方

文件:`job_files/yi_7_9PS4.json` · `sequences/yi_7_9PS4.fasta`

#### 乙班 8. 9NJY — 第 7 类 抗原抗体复合物

Terminal two domains of ClfA002 with bound Fab of AZD7745

- X-RAY DIFFRACTION / 1.58 Å · deposit 2025-02-28 · release 2025-07-30 · 估算 token 1428
- 这一类考察什么:抗体 CDR 环构象 + 表位定位。抗体骨架本身高度保守(实测 ≥80%),无法规避,真正未知的是表位识别;建议跑多个 seed 按 ipTM 排序
- 评估重点:多链:ipTM 与 chain_pair_iptm 看界面,pLDDT 看各链自身

输入:
- 蛋白质链:989 aa — Clumping factor A
- 　└ 全长序列 Q99VJ4(构建体 302 aa,实验结构里只解出 302 个残基)
- 蛋白质链:225 aa — Human antibody, heavy chain fragment, antigen binding
- 　└ 用沉积序列:合成肽 / 人工设计蛋白,无天然全长序列
- 蛋白质链:214 aa — Human antibody, light chain, antigen binding
- 　└ 用沉积序列:合成肽 / 人工设计蛋白,无天然全长序列

注意:
- 去垢剂/结晶助剂/不支持的重原子已剔除,不要输入:GOL×2
- 建议跑 3–5 个不同 seed,按 ranking_score / chain_pair_iptm 选最优模型

文件:`job_files/yi_8_9NJY.json` · `sequences/yi_8_9NJY.fasta`

#### 乙班 9. 9P4C — 第 7 类 抗原抗体复合物

Crystal structure of Mesothelin C-terminal peptide-RO4 Fab complex

- X-RAY DIFFRACTION / 1.52 Å · deposit 2025-06-16 · release 2025-11-12 · 估算 token 1049
- 这一类考察什么:抗体 CDR 环构象 + 表位定位。抗体骨架本身高度保守(实测 ≥80%),无法规避,真正未知的是表位识别;建议跑多个 seed 按 ipTM 排序
- 评估重点:多链:ipTM 与 chain_pair_iptm 看界面,pLDDT 看各链自身

输入:
- 蛋白质链:622 aa — Mesothelin, cleaved form
- 　└ 全长序列 Q13421(构建体 17 aa,实验结构里只解出 17 个残基)
- 蛋白质链:211 aa — Heavy chain of the RO4 fab fragment
- 　└ 用沉积序列:合成肽 / 人工设计蛋白,无天然全长序列
- 蛋白质链:216 aa — Light chain of the RO4 Fab fragment
- 　└ 用沉积序列:合成肽 / 人工设计蛋白,无天然全长序列

注意:
- 去垢剂/结晶助剂/不支持的重原子已剔除,不要输入:GOL×1
- 建议跑 3–5 个不同 seed,按 ranking_score / chain_pair_iptm 选最优模型

文件:`job_files/yi_9_9P4C.json` · `sequences/yi_9_9P4C.fasta`

#### 乙班 10. 9S9E — 第 7 类 抗原抗体复合物

Co-crystal of broadly neutralizing biparatopic monomeric VHH in complex with cardiotoxin (P01468) naja pallida

- X-RAY DIFFRACTION / 1.31 Å · deposit 2025-08-06 · release 2026-09-02 · 估算 token 370
- 这一类考察什么:抗体 CDR 环构象 + 表位定位。抗体骨架本身高度保守(实测 ≥80%),无法规避,真正未知的是表位识别;建议跑多个 seed 按 ipTM 排序
- 评估重点:多链:ipTM 与 chain_pair_iptm 看界面,pLDDT 看各链自身

输入:
- 蛋白质链 ×2:125 aa — Variable Domain of Heavy-Chain only Antibody (VHH) TPL0870_01_G09_Wt
- 　└ 用沉积序列:合成肽 / 人工设计蛋白,无天然全长序列
- 蛋白质链 ×2:60 aa — Cytotoxin 1
- 　└ 全长序列 P01468(构建体 60 aa)

注意:
- 去垢剂/结晶助剂/不支持的重原子已剔除,不要输入:ACT×3
- 以下链在实验结构里并未全部解出(Variable Domain  94%),比较时只对齐两者共有的残基
- 建议跑 3–5 个不同 seed,按 ranking_score / chain_pair_iptm 选最优模型

文件:`job_files/yi_10_9S9E.json` · `sequences/yi_10_9S9E.fasta`

### 丙班

#### 丙班 1. 9SCD — 第 1 类 异源蛋白复合体

Structure of S. pombe PNUTS (565 - 644) bound to Swd2.2 - crystal form 2

- X-RAY DIFFRACTION / 1.52 Å · deposit 2025-08-10 · release 2026-08-05 · 估算 token 1051
- 这一类考察什么:两条以上不同蛋白链的复合体(不含抗体链、不含核酸、非膜蛋白)。考察跨链共进化信号能否把界面摆对
- 评估重点:多链:ipTM 与 chain_pair_iptm 看界面,pLDDT 看各链自身

输入:
- 蛋白质链:710 aa — Cleavage and polyadenylation factor complex subunit C74.02c
- 　└ 全长序列 O74535(构建体 87 aa,实验结构里只解出 79 个残基)
- 　└ 训练截止前同源体:无同源体
- 蛋白质链:341 aa — Uncharacterized WD repeat-containing protein C824.04
- 　└ 全长序列 Q9UT39(构建体 341 aa,实验结构里只解出 325 个残基)
- 　└ 训练截止前同源体:无同源体

注意:
- 以下链在实验结构里并未全部解出(Cleavage and pol 90%),比较时只对齐两者共有的残基

文件:`job_files/bing_1_9SCD.json` · `sequences/bing_1_9SCD.fasta`

#### 丙班 2. 9U19 — 第 2 类 蛋白 + DNA

Crystal structure of the NKX2.1 homeodomain in complex with a 12-bp DNA duplex containing a CACG motif variant

- X-RAY DIFFRACTION / 1.18 Å · deposit 2026-01-29 · release 2026-03-11 · 估算 token 396
- 这一类考察什么:只含标准 A/C/G/T;双链须分别输入两条互补链
- 评估重点:多链:ipTM 与 chain_pair_iptm 看界面,pLDDT 看各链自身

输入:
- 蛋白质链:371 aa — Homeobox protein Nkx-2.1
- 　└ 全长序列 P43699(构建体 59 aa,实验结构里只解出 59 个残基)
- DNA 链:12 nt — 12bp reverse complementary
- DNA 链:12 nt — CACG-containing 12bp forward strand
- 离子:MG ×1

注意:
- 两条 DNA 链互为反向互补,是同一段双链的两股,必须都输入

文件:`job_files/bing_2_9U19.json` · `sequences/bing_2_9U19.fasta`

#### 丙班 3. 9VHE — 第 3 类 蛋白 + RNA

cryoEM structure of retron-Eco7 complex

- ELECTRON MICROSCOPY / 2.50 Å · deposit 2025-06-16 · release 2025-12-31 · 估算 token 2960
- 这一类考察什么:RNA 构象自由度大,是公认较难的一类
- 评估重点:多链:ipTM 与 chain_pair_iptm 看界面,pLDDT 看各链自身

输入:
- DNA 链:74 nt — msdDNA
- RNA 链:65 nt — msrRNA
- 蛋白质链:216 aa — Retron Ec78 putative HNH endonuclease
- 　└ 全长序列 P0DV92(构建体 216 aa,实验结构里只解出 210 个残基)
- 蛋白质链 ×4:550 aa — Retron Ec78 probable ATPase
- 　└ 全长序列 P0DV91(构建体 550 aa,实验结构里只解出 447 个残基)
- 蛋白质链:311 aa — Retron Ec78 reverse transcriptase
- 　└ 全长序列 Q46666(构建体 311 aa,实验结构里只解出 273 个残基)
- 配体(下拉菜单内置):ATP ×3 — ADENOSINE-5'-TRIPHOSPHATE
- 离子:MG ×1

注意:
- 以下链在实验结构里并未全部解出(Retron Ec78 prob 81%、Retron Ec78 reve 87%),比较时只对齐两者共有的残基

文件:`job_files/bing_3_9VHE.json` · `sequences/bing_3_9VHE.fasta`

#### 丙班 4. 29LA — 第 4 类 蛋白 + 小分子配体

L-DOPA extradiol dioxygenase from Beta vulgaris in complex with 4-nitrocatechol

- X-RAY DIFFRACTION / 1.13 Å · deposit 2026-03-19 · release 2026-08-12 · 估算 token 289
- 这一类考察什么:配体不在 Server 的 19 种内置辅因子里,走任意 CCD 代码那条路;代码已核实在冻结的 CCD 2024_10_28 字典中
- 评估重点:单链:pLDDT + Cα RMSD / TM-score

输入:
- 蛋白质链:275 aa — 4,5-DOPA dioxygenase extradiol 1
- 　└ 全长序列 I3PFJ9(构建体 278 aa,实验结构里只解出 265 个残基)
- 　└ 训练截止前同源体:>=30%(1 个)
- 配体(任意 CCD):4NC ×1 — 4-NITROCATECHOL
- 离子:K ×1
- 离子:CL ×1
- 离子:MG ×1

注意:
- 以下配体不在下拉菜单里,而是任意 CCD 代码:4NC。上传 job JSON 会自动填入 “CCD Code” 栏;只有从零手工搭建输入时才需自己敲
- 去垢剂/结晶助剂/不支持的重原子已剔除,不要输入:PEG×1、NI×1、EDO×1

文件:`job_files/bing_4_29LA.json` · `sequences/bing_4_29LA.fasta`

#### 丙班 5. 9T9W — 第 5 类 有翻译后修饰的蛋白

Crystal structure of beta-TrCP bound by diphosphorylated I-kappa-B-alpha degron peptide

- X-RAY DIFFRACTION / 1.16 Å · deposit 2025-11-17 · release 2026-04-08 · 估算 token 922
- 这一类考察什么:修饰残基经坐标实测介导蛋白-蛋白互作(最近距离 2.2–3.0 Å,4 Å 内接触原子 ≥ 20)。该类受体多为反复研究的识别模块,同源体无法避免,已在同源体数量上取最低
- 评估重点:多链:ipTM 与 chain_pair_iptm 看界面,pLDDT 看各链自身
- 修饰介导互作的实测证据:SEP 到对方链最近 2.61 Å,4 Å 内接触原子 32 个

输入:
- 蛋白质链:317 aa — NF-kappa-B inhibitor alpha
- 　└ 全长序列 P25963(构建体 12 aa,实验结构里只解出 10 个残基)
- 　└ 训练截止前同源体:无同源体
- 　└ 翻译后修饰:SEP 加在第 32 位(S,构建体编号 5)
- 　└ 翻译后修饰:SEP 加在第 36 位(S,构建体编号 9)
- 蛋白质链:605 aa — F-box/WD repeat-containing protein 1A
- 　└ 全长序列 Q9Y297(构建体 365 aa,实验结构里只解出 355 个残基)
- 　└ 训练截止前同源体:>=95%(7 个)

注意:
- 修饰用 proteinChain 的 modifications 字段(ptmType + ptmPosition),位置已换算为全长编号
- 去垢剂/结晶助剂/不支持的重原子已剔除,不要输入:EDO×1
- 以下链在实验结构里并未全部解出(NF-kappa-B inhib 83%),比较时只对齐两者共有的残基

文件:`job_files/bing_5_9T9W.json` · `sequences/bing_5_9T9W.fasta`

#### 丙班 6. 9M0S — 第 6 类 膜蛋白

Acetyl-CoA-bound SLC33A1 in a cytoplasm-facing conformation

- ELECTRON MICROSCOPY / 3.50 Å · deposit 2025-02-25 · release 2025-04-23 · 估算 token 600
- 这一类考察什么:带 PDBTM / MemProtMD / mpstruc 跨膜注释。Server 不建模膜平面,跨膜螺旋排布与构象态最易错;本类优先选训练窗口内查不到同源体的
- 评估重点:单链:pLDDT + Cα RMSD / TM-score

输入:
- 蛋白质链:549 aa — Acetyl-coenzyme A transporter 1
- 　└ 全长序列 O00400(构建体 560 aa,实验结构里只解出 424 个残基)
- 　└ 训练截止前同源体:无同源体
- 配体(任意 CCD):ACO ×1 — ACETYL COENZYME *A

注意:
- 以下配体不在下拉菜单里,而是任意 CCD 代码:ACO。上传 job JSON 会自动填入 “CCD Code” 栏;只有从零手工搭建输入时才需自己敲
- 以下链在实验结构里并未全部解出(Acetyl-coenzyme  75%),比较时只对齐两者共有的残基
- Server 不知道膜平面,跨膜螺旋的相对排布是本题最可能出错的地方

文件:`job_files/bing_6_9M0S.json` · `sequences/bing_6_9M0S.fasta`

#### 丙班 7. 9UET — 第 6 类 膜蛋白

Cryo-EM structure of human choline-phosphotransferase 1

- ELECTRON MICROSCOPY / 3.68 Å · deposit 2025-04-09 · release 2025-06-18 · 估算 token 814
- 这一类考察什么:带 PDBTM / MemProtMD / mpstruc 跨膜注释。Server 不建模膜平面,跨膜螺旋排布与构象态最易错;本类优先选训练窗口内查不到同源体的
- 评估重点:多链:ipTM 与 chain_pair_iptm 看界面,pLDDT 看各链自身

输入:
- 蛋白质链 ×2:406 aa — Cholinephosphotransferase 1
- 　└ 全长序列 Q8WUD6(构建体 406 aa,实验结构里只解出 368 个残基)
- 　└ 训练截止前同源体:无同源体
- 离子:MG ×2

注意:
- 去垢剂/结晶助剂/不支持的重原子已剔除,不要输入:POV×2
- 以下链在实验结构里并未全部解出(Cholinephosphotr 90%),比较时只对齐两者共有的残基
- Server 不知道膜平面,跨膜螺旋的相对排布是本题最可能出错的地方

文件:`job_files/bing_7_9UET.json` · `sequences/bing_7_9UET.fasta`

#### 丙班 8. 9NKZ — 第 7 类 抗原抗体复合物

Crystal structure of Fab MAM01 in complex with NANP6 peptide from circumsporozoite protein

- X-RAY DIFFRACTION / 1.48 Å · deposit 2025-03-02 · release 2026-03-04 · 估算 token 844
- 这一类考察什么:抗体 CDR 环构象 + 表位定位。抗体骨架本身高度保守(实测 ≥80%),无法规避,真正未知的是表位识别;建议跑多个 seed 按 ipTM 排序
- 评估重点:多链:ipTM 与 chain_pair_iptm 看界面,pLDDT 看各链自身

输入:
- 蛋白质链:230 aa — Heavy Chain of Fab MAM01
- 　└ 用沉积序列:合成肽 / 人工设计蛋白,无天然全长序列
- 蛋白质链:217 aa — Light Chain of Fab MAM01
- 　└ 用沉积序列:合成肽 / 人工设计蛋白,无天然全长序列
- 蛋白质链:397 aa — Circumsporozoite protein
- 　└ 全长序列 Q7K740(构建体 24 aa,实验结构里只解出 15 个残基)

注意:
- 以下链在实验结构里并未全部解出(Circumsporozoite 62%),比较时只对齐两者共有的残基
- 建议跑 3–5 个不同 seed,按 ranking_score / chain_pair_iptm 选最优模型

文件:`job_files/bing_8_9NKZ.json` · `sequences/bing_8_9NKZ.fasta`

#### 丙班 9. 9VDY — 第 7 类 抗原抗体复合物

hA5-6 Fab bound to SFTSV glycoprotein Gn

- X-RAY DIFFRACTION / 2.28 Å · deposit 2025-06-09 · release 2026-01-07 · 估算 token 1523
- 这一类考察什么:抗体 CDR 环构象 + 表位定位。抗体骨架本身高度保守(实测 ≥80%),无法规避,真正未知的是表位识别;建议跑多个 seed 按 ipTM 排序
- 评估重点:多链:ipTM 与 chain_pair_iptm 看界面,pLDDT 看各链自身

输入:
- 蛋白质链:218 aa — hA5-6 Fab light chain
- 　└ 用沉积序列:合成肽 / 人工设计蛋白,无天然全长序列
- 蛋白质链:232 aa — hA5-6 Fab heavy chain
- 　└ 用沉积序列:合成肽 / 人工设计蛋白,无天然全长序列
- 蛋白质链:1073 aa — Envelopment polyprotein
- 　└ 全长序列 W5VWE0(构建体 338 aa,实验结构里只解出 314 个残基)

注意:
- 以下链在实验结构里并未全部解出(hA5-6 Fab heavy  92%、Envelopment poly 92%),比较时只对齐两者共有的残基
- 建议跑 3–5 个不同 seed,按 ranking_score / chain_pair_iptm 选最优模型

文件:`job_files/bing_9_9VDY.json` · `sequences/bing_9_9VDY.fasta`

#### 丙班 10. 9ZRO — 第 7 类 抗原抗体复合物

Neutralizing W037 Fab antibody fragment in complex with West Nile Virus EDIII

- X-RAY DIFFRACTION / 1.40 Å · deposit 2025-12-20 · release 2026-07-29 · 估算 token 3871
- 这一类考察什么:抗体 CDR 环构象 + 表位定位。抗体骨架本身高度保守(实测 ≥80%),无法规避,真正未知的是表位识别;建议跑多个 seed 按 ipTM 排序
- 评估重点:多链:ipTM 与 chain_pair_iptm 看界面,pLDDT 看各链自身

输入:
- 蛋白质链:3433 aa — Envelope protein E
- 　└ 全长序列 Q9Q6P4(构建体 102 aa,实验结构里只解出 99 个残基)
- 蛋白质链:214 aa — W037 Fab Light Chain
- 　└ 用沉积序列:合成肽 / 人工设计蛋白,无天然全长序列
- 蛋白质链:224 aa — W037 Fab Heavy Chain
- 　└ 用沉积序列:合成肽 / 人工设计蛋白,无天然全长序列

注意:
- 建议跑 3–5 个不同 seed,按 ranking_score / chain_pair_iptm 选最优模型

文件:`job_files/bing_10_9ZRO.json` · `sequences/bing_10_9ZRO.fasta`

### 丁班

#### 丁班 1. 9SDC — 第 1 类 异源蛋白复合体

RelSI toxin-antitoxin complex

- X-RAY DIFFRACTION / 1.70 Å · deposit 2025-08-13 · release 2026-07-22 · 估算 token 644
- 这一类考察什么:两条以上不同蛋白链的复合体(不含抗体链、不含核酸、非膜蛋白)。考察跨链共进化信号能否把界面摆对
- 评估重点:多链:ipTM 与 chain_pair_iptm 看界面,pLDDT 看各链自身

输入:
- 蛋白质链 ×4:77 aa — Toxin
- 　└ 全长序列 I6X520(构建体 77 aa)
- 　└ 训练截止前同源体:无同源体
- 蛋白质链 ×4:84 aa — RelI
- 　└ 全长序列 I6Y9Z5(构建体 84 aa,实验结构里只解出 77 个残基)
- 　└ 训练截止前同源体:无同源体

注意:
- 以下链在实验结构里并未全部解出(RelI 91%),比较时只对齐两者共有的残基

文件:`job_files/ding_1_9SDC.json` · `sequences/ding_1_9SDC.fasta`

#### 丁班 2. 9Y1J — 第 2 类 蛋白 + DNA

S180R human DNA polymerase beta, Ternary complex dT:dAmpCpp

- X-RAY DIFFRACTION / 1.55 Å · deposit 2025-08-29 · release 2026-01-28 · 估算 token 398
- 这一类考察什么:只含标准 A/C/G/T;双链须分别输入两条互补链
- 评估重点:多链:ipTM 与 chain_pair_iptm 看界面,pLDDT 看各链自身

输入:
- DNA 链:5 nt — Downstream
- DNA 链:10 nt — Primer strand
- DNA 链:16 nt — Template strand
- 蛋白质链:335 aa — DNA polymerase beta
- 　└ 全长序列 P06746(构建体 335 aa,实验结构里只解出 326 个残基)
- 配体(任意 CCD):F2A ×1 — 2'-deoxy-5'-O-[(S)-hydroxy{[(S)-hydroxy(phosphonooxy)pho
- 离子:NA ×1
- 离子:MG ×1

注意:
- 以下配体不在下拉菜单里,而是任意 CCD 代码:F2A。上传 job JSON 会自动填入 “CCD Code” 栏;只有从零手工搭建输入时才需自己敲

文件:`job_files/ding_2_9Y1J.json` · `sequences/ding_2_9Y1J.fasta`

#### 丁班 3. 9TEL — 第 3 类 蛋白 + RNA

Structure of chicken LGP2 bound to 10-mer RNA mismatched duplex that mimics the influenza B virus vRNA promoter (panhandle) and to ADP-AlF4-Mg.

- X-RAY DIFFRACTION / 1.44 Å · deposit 2025-11-25 · release 2025-12-31 · 估算 token 728
- 这一类考察什么:RNA 构象自由度大,是公认较难的一类
- 评估重点:多链:ipTM 与 chain_pair_iptm 看界面,pLDDT 看各链自身

输入:
- RNA 链:10 nt — RNA (5'-R(P*AP*GP*UP*AP*GP*UP*AP*AP*CP*A)-3')
- RNA 链:10 nt — RNA (5'-R(*UP*GP*CP*UP*UP*CP*UP*GP*CP*U)-3')
- 蛋白质链:674 aa — RNA helicase
- 　└ 全长序列 G0YYQ5(构建体 673 aa,实验结构里只解出 661 个残基)
- 配体(下拉菜单内置):ADP ×1 — ADENOSINE-5'-DIPHOSPHATE
- 配体(任意 CCD):ALF ×1 — TETRAFLUOROALUMINATE ION
- 离子:MG ×1
- 离子:ZN ×1

注意:
- 以下配体不在下拉菜单里,而是任意 CCD 代码:ALF。上传 job JSON 会自动填入 “CCD Code” 栏;只有从零手工搭建输入时才需自己敲

文件:`job_files/ding_3_9TEL.json` · `sequences/ding_3_9TEL.fasta`

#### 丁班 4. 21ZG — 第 4 类 蛋白 + 小分子配体

Crystal structure of the petrobactin-binding protein FatB from Bacillus cereus complexed with ferric siderophore mimic, Fe(3,4-DHB)2

- X-RAY DIFFRACTION / 1.40 Å · deposit 2026-01-04 · release 2026-04-22 · 估算 token 350
- 这一类考察什么:配体不在 Server 的 19 种内置辅因子里,走任意 CCD 代码那条路;代码已核实在冻结的 CCD 2024_10_28 字典中
- 评估重点:单链:pLDDT + Cα RMSD / TM-score

输入:
- 蛋白质链:338 aa — Ferric anguibactin-binding protein
- 　└ 全长序列 Q815N5(构建体 302 aa,实验结构里只解出 295 个残基)
- 　└ 训练截止前同源体:>=30%(1 个)
- 配体(任意 CCD):DHB ×1 — 3,4-DIHYDROXYBENZOIC ACID
- 离子:FE ×1

注意:
- 以下配体不在下拉菜单里,而是任意 CCD 代码:DHB。上传 job JSON 会自动填入 “CCD Code” 栏;只有从零手工搭建输入时才需自己敲
- 去垢剂/结晶助剂/不支持的重原子已剔除,不要输入:IPA×1、EDO×1

文件:`job_files/ding_4_21ZG.json` · `sequences/ding_4_21ZG.fasta`

#### 丁班 5. 9X8S — 第 5 类 有翻译后修饰的蛋白

Crystal structure of the human GAS41 YEATS domain in complex with an acetylated YFV capsid peptide (K4ac)

- X-RAY DIFFRACTION / 1.70 Å · deposit 2025-10-20 · release 2025-11-26 · 估算 token 922
- 这一类考察什么:修饰残基经坐标实测介导蛋白-蛋白互作(最近距离 2.2–3.0 Å,4 Å 内接触原子 ≥ 20)。该类受体多为反复研究的识别模块,同源体无法避免,已在同源体数量上取最低
- 评估重点:多链:ipTM 与 chain_pair_iptm 看界面,pLDDT 看各链自身
- 修饰介导互作的实测证据:ALY 到对方链最近 2.66 Å,4 Å 内接触原子 57 个

输入:
- 蛋白质链 ×2:7 aa — Yellow Fever Virus Capsid Protein
- 　└ 用沉积序列:合成肽 / 人工设计蛋白,无天然全长序列
- 　└ 翻译后修饰:ALY 加在第 4 位(K,构建体编号 4)
- 蛋白质链 ×4:227 aa — YEATS domain-containing protein 4
- 　└ 全长序列 O95619(构建体 133 aa,实验结构里只解出 133 个残基)
- 　└ 训练截止前同源体:>=95%(7 个)

注意:
- 修饰用 proteinChain 的 modifications 字段(ptmType + ptmPosition),位置已换算为全长编号
- 含极短肽链(Yellow Fever Virus Capsid Protein 7 aa):FAQ 明示 pTM 对短于 16 残基的链系统性偏低,该值接近 0 不代表预测失败。评估以 pLDDT / PAE 为主;要看界面就取 chain_pair_iptm 里「受体链 × 该肽链」那一格,不要用整体 ipTM

文件:`job_files/ding_5_9X8S.json` · `sequences/ding_5_9X8S.fasta`

#### 丁班 6. 9M2H — 第 6 类 膜蛋白

Structure of the auxin importer AUX1 in Arabidopsis thaliana in the CHPAA-bound state

- ELECTRON MICROSCOPY / 3.40 Å · deposit 2025-02-27 · release 2025-05-28 · 估算 token 497
- 这一类考察什么:带 PDBTM / MemProtMD / mpstruc 跨膜注释。Server 不建模膜平面,跨膜螺旋排布与构象态最易错;本类优先选训练窗口内查不到同源体的
- 评估重点:单链:pLDDT + Cα RMSD / TM-score

输入:
- 蛋白质链:485 aa — Auxin transporter protein 1
- 　└ 全长序列 Q96247(构建体 485 aa,实验结构里只解出 431 个残基)
- 　└ 训练截止前同源体:无同源体
- 配体(任意 CCD):3C4 ×1 — (3-CHLORO-4-HYDROXYPHENYL)ACETIC ACID

注意:
- 以下配体不在下拉菜单里,而是任意 CCD 代码:3C4。上传 job JSON 会自动填入 “CCD Code” 栏;只有从零手工搭建输入时才需自己敲
- 以下链在实验结构里并未全部解出(Auxin transporte 88%),比较时只对齐两者共有的残基
- Server 不知道膜平面,跨膜螺旋的相对排布是本题最可能出错的地方

文件:`job_files/ding_6_9M2H.json` · `sequences/ding_6_9M2H.fasta`

#### 丁班 7. 9N93 — 第 6 类 膜蛋白

Human TMEM63A mutant V53M lipid-open state

- ELECTRON MICROSCOPY / 2.95 Å · deposit 2025-02-10 · release 2025-06-11 · 估算 token 807
- 这一类考察什么:带 PDBTM / MemProtMD / mpstruc 跨膜注释。Server 不建模膜平面,跨膜螺旋排布与构象态最易错;本类优先选训练窗口内查不到同源体的
- 评估重点:单链:pLDDT + Cα RMSD / TM-score

输入:
- 蛋白质链:807 aa — CSC1-like protein 1
- 　└ 全长序列 O94886(构建体 807 aa,实验结构里只解出 650 个残基)
- 　└ 训练截止前同源体:无同源体

注意:
- 以下链在实验结构里并未全部解出(CSC1-like protei 80%),比较时只对齐两者共有的残基
- Server 不知道膜平面,跨膜螺旋的相对排布是本题最可能出错的地方

文件:`job_files/ding_7_9N93.json` · `sequences/ding_7_9N93.fasta`

#### 丁班 8. 9IA3 — 第 7 类 抗原抗体复合物

Bc8.108 Fab bound to preS2 peptide

- X-RAY DIFFRACTION / 1.11 Å · deposit 2025-02-07 · release 2025-11-26 · 估算 token 634
- 这一类考察什么:抗体 CDR 环构象 + 表位定位。抗体骨架本身高度保守(实测 ≥80%),无法规避,真正未知的是表位识别;建议跑多个 seed 按 ipTM 排序
- 评估重点:多链:ipTM 与 chain_pair_iptm 看界面,pLDDT 看各链自身

输入:
- 蛋白质链:191 aa — Large S protein
- 　└ 全长序列 B2Y6K4(构建体 23 aa,实验结构里只解出 23 个残基)
- 蛋白质链:228 aa — Fab Bc8.108 heavy chain
- 　└ 用沉积序列:合成肽 / 人工设计蛋白,无天然全长序列
- 蛋白质链:215 aa — Fab Bc8.108 light chain
- 　└ 用沉积序列:合成肽 / 人工设计蛋白,无天然全长序列

注意:
- 去垢剂/结晶助剂/不支持的重原子已剔除,不要输入:PEG×1、EDO×2
- 建议跑 3–5 个不同 seed,按 ranking_score / chain_pair_iptm 选最优模型

文件:`job_files/ding_8_9IA3.json` · `sequences/ding_8_9IA3.fasta`

#### 丁班 9. 9NW4 — 第 7 类 抗原抗体复合物

Structure of CISV1 antibody bound to PvCSP repeat peptide

- X-RAY DIFFRACTION / 1.82 Å · deposit 2025-03-21 · release 2026-03-25 · 估算 token 827
- 这一类考察什么:抗体 CDR 环构象 + 表位定位。抗体骨架本身高度保守(实测 ≥80%),无法规避,真正未知的是表位识别;建议跑多个 seed 按 ipTM 排序
- 评估重点:多链:ipTM 与 chain_pair_iptm 看界面,pLDDT 看各链自身

输入:
- 蛋白质链:229 aa — CISV1 Fab Heavy Chain
- 　└ 用沉积序列:合成肽 / 人工设计蛋白,无天然全长序列
- 蛋白质链:220 aa — CISV1 Fab Light Chain
- 　└ 用沉积序列:合成肽 / 人工设计蛋白,无天然全长序列
- 蛋白质链:378 aa — PvCSPvk210 peptide from Circumsporozoite protein
- 　└ 全长序列 P08677(构建体 18 aa,实验结构里只解出 10 个残基)

注意:
- 以下链在实验结构里并未全部解出(CISV1 Fab Heavy  92%、PvCSPvk210 pepti 55%),比较时只对齐两者共有的残基
- 建议跑 3–5 个不同 seed,按 ranking_score / chain_pair_iptm 选最优模型

文件:`job_files/ding_9_9NW4.json` · `sequences/ding_9_9NW4.fasta`

#### 丁班 10. 9Q1L — 第 7 类 抗原抗体复合物

Crystal structure of the walnut allergen Jug r 2 bound to the human-derived Fab 6D12

- X-RAY DIFFRACTION / 1.56 Å · deposit 2025-08-14 · release 2026-01-21 · 估算 token 1038
- 这一类考察什么:抗体 CDR 环构象 + 表位定位。抗体骨架本身高度保守(实测 ≥80%),无法规避,真正未知的是表位识别;建议跑多个 seed 按 ipTM 排序
- 评估重点:多链:ipTM 与 chain_pair_iptm 看界面,pLDDT 看各链自身

输入:
- 蛋白质链:214 aa — 6D12 Fab light chain
- 　└ 用沉积序列:合成肽 / 人工设计蛋白,无天然全长序列
- 蛋白质链:230 aa — 6D12 Fab Heavy Chain
- 　└ 用沉积序列:合成肽 / 人工设计蛋白,无天然全长序列
- 蛋白质链:593 aa — Vicilin Jug r 2.0102 hairpinin alpha 4
- 　└ 全长序列 Q9SEW4(构建体 47 aa,实验结构里只解出 39 个残基)
- 离子:CL ×1

注意:
- 去垢剂/结晶助剂/不支持的重原子已剔除,不要输入:EDO×4
- 以下链在实验结构里并未全部解出(Vicilin Jug r 2. 83%),比较时只对齐两者共有的残基
- 建议跑 3–5 个不同 seed,按 ranking_score / chain_pair_iptm 选最优模型

文件:`job_files/ding_10_9Q1L.json` · `sequences/ding_10_9Q1L.fasta`

## 目录

```
README.md          本文件
targets.json       41 题的结构化数据
job_files/         41 个可直接 Upload JSON 导入的作业文件
sequences/         41 个 FASTA
all_41_jobs.json   一次性导入全部 41 题
pipeline/          筛选流程的全部脚本(可复现)
evidence/          筛选依据:同源筛查、界面实测、装配拷贝数等原始结果
```

## 数据来源

结构数据来自 [RCSB PDB](https://www.rcsb.org)(CC0)。
AlphaFold Server 的限制依据其 [FAQ](https://alphafoldserver.com/faq) 与
[Release Updates](https://alphafoldserver.com/release-updates);
任意 CCD 配体输入自 2026-08-19 起开放。
模型见 Abramson et al., *Nature* 630:493–500 (2024)。