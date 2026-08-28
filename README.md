# PrimerFold · DNA 引物二级结构本地分析 GUI

PrimerFold 是一个完全在本机运行的 DNA 引物二级结构分析界面：输入引物序列，
一键获得发卡结构、自二聚体与交叉二聚体的热力学预测、彩色结构图与逐位对照，
每一项指标都附带解释信息。序列不上传网络。

后端使用 ViennaRNA 2.7.2 的 `RNAfold`、`RNAduplex` 和 `RNAplot`，固定启用
DNA 参数（Mathews 2004）并禁止 T→U 自动转换。

## 功能特性

- 🔬 **三种分析**：发卡结构（RNAfold）、自二聚体与交叉二聚体（RNAduplex，仅允许链间配对）。
- 🧾 **引物属性**：Tm（SantaLucia 1998 最近邻法，含单价盐校正）、GC 含量、分子量、
  3′GC clamp 与最长同聚串，纯序列本地计算。
- 🔥 **互作热图**：N×N 矩阵总览所有引物两两二聚体强度，颜色越红结合越强，
  对角线为自二聚体，适合快速定位多重 PCR 中需要关注的组合。
- 📄 **HTML 报告导出**：一键生成自包含报告（含指标表、热图、结构图、逐条解读与
  逐位对照），可存档、打印或分享；CSV 导出同步保留。
- 🎨 **美化结构图**：碱基按类型着色，实心圆＝配对/空心圆＝未配对，G≡C 与 A=T 配对短线区分，
  紫色虚线圈标出 3′端最后 5 nt，图内含 5′/3′ 标签与完整图例，可导出 SVG。
- 📖 **每块数据都有解释**：列头 ⓘ 提示、指标说明对照表、逐条结果的自动“结果解读”。
- 🔤 **序列/结构逐位对照**：配对括号按嵌套深度着色、链间配对洋红高亮、3′端下划线。
- 🖥️ **本地运行**：仅监听 127.0.0.1，不上传任何序列；不使用通用阈值做“合格/不合格”判定。

## 安装

需要 Python ≥ 3.10 与 [ViennaRNA](https://www.tbi.univie.ac.at/RNA/) 2.7+
命令行工具（`RNAfold`、`RNAduplex`、`RNAplot` 需在 PATH 中）：

```bash
# 方式一：conda 安装 ViennaRNA
conda create -n rnafold -c bioconda -c conda-forge viennarna=2.7 python=3.11
conda activate rnafold

# 安装 PrimerFold
git clone https://github.com/X-1546520470/rnafold.git
cd rnafold
pip install .
```

macOS 快捷方式：双击 `launch.command`（默认使用
`/opt/anaconda3/envs/rnafold` 环境中的 Python，可在脚本开头修改路径）。

## 启动

```bash
primerfold-gui            # pip 安装后的命令行入口
# 或
python -m primerfold      # 等价入口
python -m primerfold --port 8501 --headless   # 自定义端口 / 不自动开浏览器
```

浏览器会自动打开 `http://127.0.0.1:8501`。停止服务时在终端按 `Control-C`。

## 输入格式

支持以下任意一种格式，序列方向均为 5′→3′：

```text
>Primer_F
ATGACCATGATTACGCCAAG
>Primer_R
GCGCGCTTTTTGCGCGC
```

```text
Primer_F,ATGACCATGATTACGCCAAG
Primer_R,GCGCGCTTTTTGCGCGC
```

```text
ATGACCATGATTACGCCAAG
GCGCGCTTTTTGCGCGC
```

仅接受 A、C、G、T。一次最多 30 条引物，每条 4–500 nt。输入框下方会实时显示
解析到的引物名称与长度，便于在上传前发现格式问题。

## 结果展示与解释

界面中每个数据块都配有对应的解释信息：

- 汇总指标（引物数、分析任务、温度、单价盐）悬停可见解释。
- 结果表格的每个列头带 ⓘ 提示；“3′末 5 nt 配对”以进度条显示（0–5）。
- “指标说明”展开面板给出所有指标的含义与解读建议。
- 每条结果附带自动生成的“结果解读”段落，用描述性语言说明结构组成、
  MFE 强度与 3′末端状态；工具不做“合格/不合格”判定。
- 序列与结构逐位对照视图：配对括号按嵌套深度着色，链间配对统一为洋红，
  3′端最后 5 nt 带紫色下划线，每 10 nt 一组并标注起始位置。

## 结构图

结构图沿用 RNAplot（naview 布局）计算的坐标，仅重新着色与标注，不改变
预测结果：

- 碱基按类型着色：A 绿、C 蓝、G 橙、T 红；实心圆＝参与配对，空心圆＝未配对。
- 配对短线按氢键数着色：G≡C 青色（3 个氢键）、A=T 灰色（2 个氢键）。
- 紫色虚线圈标出每条链 3′端最后 5 nt，与表格中的末端指标一一对应。
- 二聚体在图例中标注两条链的名称与长度；图内含 5′/3′ 端标签与完整图例。

## 输出解释

- `MFE` 是给定模型条件下最低自由能结构的自由能，单位为 kcal/mol。
- 匹配的 `(` 和 `)` 表示配对，`.` 表示未配对，`&` 分隔两条链。
- 发卡结果的“3′末 5 nt 配对”统计链内配对；二聚体结果只统计链间配对。
- 二聚体由 `RNAduplex` 计算，只允许链间碱基配对，避免把单条引物自身发卡误报为二聚体。
- 温度和单价盐浓度可设置；模型不显式表示 Mg²⁺、dNTP、引物浓度或完整 PCR 缓冲液。

## 开发与测试

```bash
pip install -e .
python -m unittest discover -s tests   # 需要本机可用 ViennaRNA
```

代码结构：

- `primerfold/core.py` — 输入解析、ViennaRNA 调用与指标计算。
- `primerfold/qc.py` — 基础引物属性（Tm、GC%、分子量、GC clamp、同聚串）。
- `primerfold/plot.py` — 解析 RNAplot SVG 几何并重绘为美化结构图。
- `primerfold/present.py` — 结果解读、对照 HTML、互作热图与 HTML 报告。
- `primerfold/app.py` — Streamlit 界面。
- `primerfold/__main__.py` — `python -m primerfold` 启动入口。

## 许可证

[MIT](LICENSE)
