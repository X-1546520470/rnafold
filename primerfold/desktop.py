#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""desktop.py — PrimerFold 引物二级结构分析 (Tkinter 桌面版)

运行:  python -m primerfold.desktop   (或双击 run_desktop.command)
依赖:  Python 标准库 (tkinter) + ViennaRNA 命令行工具 (RNAfold/RNAduplex/
       RNAcofold/RNAplot); 结构图预览使用 macOS 自带的 qlmanage。

布局对齐 T7 盘 Analysis Tools 桌面工具风格:
  ① 输入引物 (FASTA 文本框, 实时识别) → ② 左侧参数 / 右侧多页签结果
"""

from __future__ import annotations

import atexit
import datetime
import shutil
import os
import queue
import subprocess
import sys
import tempfile
import threading
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from primerfold import (
    AnalysisJob,
    FoldExecutionError,
    FoldResult,
    Primer,
    PrimerInputError,
    __version__,
    build_html_report,
    build_jobs,
    cofold_energies,
    fold_dimer,
    fold_hairpin,
    parse_primers,
    qc_rows,
    render_styled_structure_svg,
    svg_pixel_dimensions,
)
from primerfold.equilibrium import interpret_equilibrium, species_equilibrium
from primerfold.present import interpret_result

KIND_TABS = (
    ("hairpin", " 发卡结构 "),
    ("self_dimer", " 自二聚体 "),
    ("cross_dimer", " 交叉二聚体 "),
)
KIND_LABELS = {
    "hairpin": "发卡结构",
    "self_dimer": "自二聚体",
    "cross_dimer": "交叉二聚体",
}

DEFAULT_INPUT = """>Primer_F
ATGACCATGATTACGCCAAG
>Primer_R
GCGCGCTTTTTGCGCGC
"""

MAX_CONCENTRATION_PAIRS = 40

HELP_TEXT = f"""\
【使用方法】
1. 在"① 输入引物"粘贴引物序列(支持 FASTA、名称+序列、CSV/TSV 或每行一条
   纯序列; 左下角会实时显示识别结果)。
2. 在左侧"② 分析条件"设置温度、单价盐与要分析的类型, 点击"▶ 开始分析"。
3. 在右侧页签查看结果: 引物属性 / 互作热图 / 三类结构明细 / 浓度平衡。

【输入格式】
• FASTA:  >Primer_F 换行 ATGC… (推荐, 可自带引物名称)
• CSV:    Primer_F,ATGC…
• 纯序列: 每行一条, 自动命名 Primer_1、Primer_2…
• 仅接受 A/C/G/T; 一次最多 30 条引物, 每条 4–500 nt。

【分析条件说明】
• 温度 (°C): ViennaRNA 热力学参数重标定温度; 37 °C 便于与常见 ΔG 报告比较。
• 单价盐 (M): Na⁺ 等效浓度, 50 mM 填 0.05; 不等同于 Mg²⁺。
• 发卡结构: 单条引物内部互补形成的茎环 (RNAfold)。
• 自二聚体: 引物与自身拷贝的链间配对 (RNAduplex, 仅允许链间配对)。
• 交叉二聚体: 不同引物之间的链间配对; 占用双方 3′末端时最危险。
• 浓度平衡: 基于 RNAcofold 集合自由能, 在给定引物浓度下计算游离单链与
  AB/AA/BB 二聚体的平衡比例 (五物种统计力学自洽求解)。

【指标解读】
• MFE (kcal/mol): 模型预测的最稳定结构的自由能, 越负结构越稳定;
  仅适合同一引物在不同条件间相对比较, 不是实验测量值。
• 总碱基对: 配对碱基总数 (发卡=链内, 二聚体=链间)。
• 3′末 5 nt 配对 / 3′末端: 3′端最后 5 nt 中参与配对的数目与最末位状态。
  3′端是 DNA 聚合酶延伸起点, 被占用可能降低延伸效率 (描述性信息)。
• Tm (°C): SantaLucia 1998 最近邻法, 按当前单价盐校正, 总链浓度 250 nM
  (非自互补假设); 不含 Mg²⁺/dNTP。
• GC (%): 常见建议 40–60%; 3′GC clamp: 3′端最后 5 nt 中 G/C 个数,
  常见建议 1–3; 最长同聚串 ≥4 (尤其 G/C) 可能增加滑移风险。
• A 被占用 (%): (AB + 2·AA)/[A]₀ — 该引物分子处于任意二聚体中的比例。

【结构图图例】
• 实心圆 = 参与配对的碱基; 空心圆 = 未配对 (A 绿 / C 蓝 / G 橙 / T 红)。
• 青色短线 = G≡C (3 氢键); 灰色短线 = A=T (2 氢键)。
• 紫色虚线圈 = 3′端最后 5 nt; 图内有 5′/3′ 端标签与完整图例。

【导出】
• "保存完整报告 (HTML)": 自包含报告 (属性表/热图/结构图/逐条解读)。
• 各结构页签可单独保存当前美化结构图 (SVG)。

【模型与边界】
• ViennaRNA 2.7.2, Mathews 2004 DNA 参数 (--paramFile=DNA), 保留 T (--noconv)。
• MFE 模型不显式表示 Mg²⁺、dNTP、引物浓度或完整 PCR 缓冲液。
• 本工具只报告描述性指标, 不使用通用阈值判定引物"合格/不合格";
  所有结果均为热力学模型预测, 请结合实际实验体系判断。

PrimerFold v{__version__} · 本地运行, 序列不上传网络
"""


class App:
    def __init__(self, root: tk.Tk):
        self.root = root
        root.title(f"PrimerFold 引物二级结构分析 v{__version__}")
        root.geometry("1240x860")
        root.minsize(1060, 720)

        self.en_hairpin = tk.BooleanVar(value=True)
        self.en_self = tk.BooleanVar(value=True)
        self.en_cross = tk.BooleanVar(value=True)
        self.en_conc = tk.BooleanVar(value=False)

        self._results: list[FoldResult] = []
        self._primers: list[Primer] = []
        self._salt_m: float = 0.05
        self._temperature_c: float = 37.0
        self._kind_results: dict[str, list[FoldResult]] = {}
        self._svg_cache: dict[str, str] = {}
        self._photo = None  # 防 PhotoImage 被回收
        self._fig_dir = tempfile.mkdtemp(prefix="primerfold-figs-")
        atexit.register(self._cleanup_figs)

        self._queue: queue.Queue = queue.Queue()
        self._worker: threading.Thread | None = None

        self._build_ui()
        self.input_txt.insert("1.0", DEFAULT_INPUT)
        self._refresh_detect()
        root.protocol("WM_DELETE_WINDOW", self._on_close)

    # ------------------------------------------------------------------
    # 界面构建
    # ------------------------------------------------------------------
    def _build_ui(self):
        style = ttk.Style(self.root)
        try:
            style.theme_use("aqua")
        except tk.TclError:
            pass

        outer = ttk.Frame(self.root, padding=8)
        outer.pack(fill="both", expand=True)

        inp = ttk.LabelFrame(outer, text=" ① 输入引物 (支持 FASTA / CSV / 纯序列) ", padding=6)
        inp.pack(fill="both", expand=False)

        bar = ttk.Frame(inp)
        bar.pack(fill="x")
        ttk.Button(bar, text="载入 FASTA 文件", command=self.load_fasta).pack(side="left")
        ttk.Button(bar, text="载入示例", command=self.load_example).pack(side="left", padx=6)
        ttk.Button(bar, text="清空", command=self.clear_input).pack(side="left")
        self.detect_var = tk.StringVar(value="待输入…")
        ttk.Label(bar, textvariable=self.detect_var, foreground="#555555").pack(side="right")

        wrap = ttk.Frame(inp)
        wrap.pack(fill="both", expand=True)
        self.input_txt = tk.Text(wrap, height=7, wrap="none", undo=True,
                                 font=("Menlo", 12), bg="#fbfbfd")
        ys = ttk.Scrollbar(wrap, orient="vertical", command=self.input_txt.yview)
        self.input_txt.configure(yscrollcommand=ys.set)
        self.input_txt.pack(side="left", fill="both", expand=True)
        ys.pack(side="right", fill="y")
        self._detect_job = None
        self.input_txt.bind("<<Modified>>", self._on_input_change)

        paned = ttk.PanedWindow(outer, orient="horizontal")
        paned.pack(fill="both", expand=True, pady=(8, 0))

        left = ttk.Frame(paned)
        paned.add(left, weight=0)
        self._build_left(left)

        right = ttk.Frame(paned)
        paned.add(right, weight=1)
        self._build_right(right)

    def _build_left(self, parent):
        panel = ttk.LabelFrame(parent, text=" ② 分析条件 ", padding=8)
        panel.pack(fill="both", expand=True)

        grid = ttk.Frame(panel)
        grid.pack(fill="x")
        ttk.Label(grid, text="温度 (°C)").grid(row=0, column=0, sticky="w", pady=2)
        self.temp_var = tk.DoubleVar(value=37.0)
        ttk.Spinbox(grid, from_=0.0, to=100.0, increment=0.5, width=7,
                    textvariable=self.temp_var).grid(row=0, column=1, sticky="w", padx=6)
        ttk.Label(grid, text="单价盐 (M)").grid(row=1, column=0, sticky="w", pady=2)
        self.salt_var = tk.StringVar(value="0.05")
        ttk.Entry(grid, width=8, textvariable=self.salt_var).grid(row=1, column=1, sticky="w", padx=6)
        ttk.Label(grid, text="50 mM → 0.05", foreground="#777777").grid(
            row=1, column=2, sticky="w")

        ttk.Separator(panel).pack(fill="x", pady=6)
        ttk.Label(panel, text="分析类型:").pack(anchor="w")
        ttk.Checkbutton(panel, text="发卡结构", variable=self.en_hairpin).pack(anchor="w")
        ttk.Checkbutton(panel, text="自二聚体", variable=self.en_self).pack(anchor="w")
        ttk.Checkbutton(panel, text="交叉二聚体", variable=self.en_cross).pack(anchor="w")

        ttk.Separator(panel).pack(fill="x", pady=6)
        ttk.Checkbutton(panel, text="启用浓度平衡分析", variable=self.en_conc).pack(anchor="w")
        conc_row = ttk.Frame(panel)
        conc_row.pack(anchor="w", pady=(2, 0))
        ttk.Label(conc_row, text="引物浓度 (nM)").pack(side="left")
        self.conc_var = tk.IntVar(value=500)
        ttk.Spinbox(conc_row, from_=1, to=100000, increment=50, width=7,
                    textvariable=self.conc_var).pack(side="left", padx=6)

        ttk.Separator(panel).pack(fill="x", pady=6)
        self.run_btn = ttk.Button(panel, text="▶ 开始分析", command=self.run_analysis)
        self.run_btn.pack(fill="x", ipady=4)
        ttk.Button(panel, text="保存完整报告 (HTML)", command=self.save_report).pack(
            fill="x", pady=(6, 0))

        self.status_var = tk.StringVar(value="就绪")
        self.status_lbl = ttk.Label(panel, textvariable=self.status_var,
                                    foreground="#0a7d32", wraplength=360, justify="left")
        self.status_lbl.pack(fill="x", pady=(8, 0))

    def _build_right(self, parent):
        nb = ttk.Notebook(parent)
        nb.pack(fill="both", expand=True)
        self.nb = nb

        # 引物属性
        f_qc = ttk.Frame(nb, padding=4)
        nb.add(f_qc, text=" 引物属性 ")
        self.qc_txt = self._make_text(f_qc)

        # 互作热图
        f_hm = ttk.Frame(nb, padding=4)
        nb.add(f_hm, text=" 互作热图 ")
        self.heat_txt = self._make_text(f_hm)

        # 三个结构页签
        self.kind_combo: dict[str, ttk.Combobox] = {}
        self.kind_detail: dict[str, tk.Text] = {}
        self.kind_fig: dict[str, ttk.Label] = {}
        self.kind_save_btn: dict[str, ttk.Button] = {}
        for kind, title in KIND_TABS:
            f = ttk.Frame(nb, padding=4)
            nb.add(f, text=title)
            bar = ttk.Frame(f)
            bar.pack(fill="x")
            ttk.Label(bar, text="选择结果:").pack(side="left")
            combo = ttk.Combobox(bar, width=42, state="readonly")
            combo.pack(side="left", padx=6)
            combo.bind("<<ComboboxSelected>>",
                       lambda _e, k=kind: self._show_kind_result(k))
            self.kind_combo[kind] = combo
            btn = ttk.Button(bar, text="保存结构图 (SVG)",
                             command=lambda k=kind: self.save_structure_svg(k))
            btn.pack(side="right")
            btn.configure(state="disabled")
            self.kind_save_btn[kind] = btn
            wrap = ttk.Frame(f)
            wrap.pack(fill="both", expand=True)
            detail = self._make_text(wrap, height_ratio=False)
            detail.pack(fill="both", expand=True)
            fig = ttk.Label(wrap, anchor="center", text="(暂无结果)")
            fig.pack(fill="x", pady=(4, 0))
            self.kind_detail[kind] = detail
            self.kind_fig[kind] = fig

        # 浓度平衡
        f_conc = ttk.Frame(nb, padding=4)
        nb.add(f_conc, text=" 浓度平衡 ")
        self.conc_txt = self._make_text(f_conc)

        # 使用说明
        f_help = ttk.Frame(nb, padding=4)
        nb.add(f_help, text=" 使用说明 ")
        self.help_txt = self._make_text(f_help)
        self._set_text(self.help_txt, HELP_TEXT)

    def _make_text(self, parent, height_ratio: bool = True) -> tk.Text:
        txt = tk.Text(parent, wrap="none", font=("Menlo", 12),
                      bg="#ffffff", relief="flat", state="disabled")
        ys = ttk.Scrollbar(parent, orient="vertical", command=txt.yview)
        xs = ttk.Scrollbar(parent, orient="horizontal", command=txt.xview)
        txt.configure(yscrollcommand=ys.set, xscrollcommand=xs.set)
        if height_ratio:
            txt.grid(row=0, column=0, sticky="nsew")
            ys.grid(row=0, column=1, sticky="ns")
            xs.grid(row=1, column=0, sticky="we")
            parent.rowconfigure(0, weight=1)
            parent.columnconfigure(0, weight=1)
        else:
            txt.pack(fill="both", expand=True)
        return txt

    def _set_text(self, widget: tk.Text, content: str):
        widget.configure(state="normal")
        widget.delete("1.0", "end")
        widget.insert("1.0", content)
        widget.configure(state="disabled")

    # ------------------------------------------------------------------
    # 输入处理
    # ------------------------------------------------------------------
    def load_fasta(self):
        path = filedialog.askopenfilename(
            title="选择序列文件",
            filetypes=[("序列文件", "*.fa *.fasta *.fna *.txt *.seq *.csv *.tsv"),
                       ("所有文件", "*.*")])
        if not path:
            return
        try:
            with open(path, encoding="utf-8-sig", errors="replace") as f:
                self.input_txt.delete("1.0", "end")
                self.input_txt.insert("1.0", f.read())
        except OSError as e:
            messagebox.showerror("读取失败", str(e))

    def load_example(self):
        self.input_txt.delete("1.0", "end")
        self.input_txt.insert("1.0", DEFAULT_INPUT)

    def clear_input(self):
        self.input_txt.delete("1.0", "end")
        self.detect_var.set("待输入…")

    def _on_input_change(self, event=None):
        self.input_txt.edit_modified(False)
        if self._detect_job is not None:
            self.root.after_cancel(self._detect_job)
        self._detect_job = self.root.after(250, self._refresh_detect)

    def _refresh_detect(self):
        self._detect_job = None
        text = self.input_txt.get("1.0", "end")
        if not text.strip():
            self.detect_var.set("待输入…")
            return
        try:
            primers = parse_primers(text)
        except PrimerInputError as exc:
            self.detect_var.set(f"⚠ {exc}")
            return
        lengths = [len(p.sequence) for p in primers]
        self.detect_var.set(
            f"已识别: {len(primers)} 条引物 (长度 {min(lengths)}–{max(lengths)} nt)"
        )

    # ------------------------------------------------------------------
    # 主流程
    # ------------------------------------------------------------------
    def _read_conditions(self) -> tuple[float, float] | None:
        temperature = float(self.temp_var.get())
        try:
            salt = float(self.salt_var.get())
        except ValueError:
            messagebox.showwarning("提示", "单价盐浓度需为数字, 例如 0.05")
            return None
        if not 0.0 <= temperature <= 100.0:
            messagebox.showwarning("提示", "温度必须在 0–100 °C 之间。")
            return None
        if not 0.0001 <= salt <= 2.0:
            messagebox.showwarning("提示", "单价盐浓度必须在 0.0001–2.0 M 之间。")
            return None
        return temperature, salt

    def run_analysis(self):
        text = self.input_txt.get("1.0", "end")
        if self._worker is not None and self._worker.is_alive():
            messagebox.showinfo("提示", "分析正在进行中, 请稍候。")
            return
        try:
            primers = parse_primers(text)
        except PrimerInputError as exc:
            messagebox.showwarning("输入有误", str(exc))
            return
        conditions = self._read_conditions()
        if conditions is None:
            return
        temperature, salt = conditions
        if not (self.en_hairpin.get() or self.en_self.get() or self.en_cross.get()):
            messagebox.showwarning("提示", "请至少勾选一种分析类型。")
            return
        try:
            jobs = build_jobs(
                primers,
                hairpin=self.en_hairpin.get(),
                self_dimer=self.en_self.get(),
                cross_dimer=self.en_cross.get(),
            )
        except PrimerInputError as exc:
            messagebox.showwarning("提示", str(exc))
            return

        self._primers = primers
        self._temperature_c = temperature
        self._salt_m = salt
        self._svg_cache.clear()
        self._results = []
        self._kind_results = {}
        self.run_btn.configure(state="disabled")
        self.status_var.set(f"准备运行 {len(jobs)} 个任务…")

        conc_enabled = self.en_conc.get()
        conc_nm = int(self.conc_var.get())
        self._worker = threading.Thread(
            target=self._analysis_worker,
            args=(jobs, primers, temperature, salt, conc_enabled, conc_nm),
            daemon=True,
        )
        self._worker.start()
        self.root.after(120, self._poll_worker)

    def _analysis_worker(self, jobs, primers, temperature, salt, conc_enabled, conc_nm):
        try:
            results: list[FoldResult] = []
            for index, job in enumerate(jobs, start=1):
                self._queue.put(("progress", f"正在分析 ({index}/{len(jobs)}): {job.label}"))
                result = self._execute_job(job, temperature, salt)
                results.append(result)
            self._queue.put(("results", (primers, results, temperature, salt)))
            if conc_enabled:
                self._queue.put(("progress", "正在计算浓度平衡 (RNAcofold)…"))
                conc = self._compute_concentration(results, temperature, salt, conc_nm)
                self._queue.put(("conc", conc))
            self._queue.put(("done", f"✔ 完成: {len(results)} 个任务"))
        except (FoldExecutionError, PrimerInputError) as exc:
            self._queue.put(("error", str(exc)))

    def _execute_job(self, job: AnalysisJob, temperature: float, salt: float) -> FoldResult:
        if job.kind == "hairpin":
            return fold_hairpin(job.primer_a, temperature_c=temperature, salt_m=salt)
        if job.primer_b is None:
            raise FoldExecutionError("二聚体任务缺少第二条引物。")
        return fold_dimer(
            job.primer_a, job.primer_b, kind=job.kind,
            temperature_c=temperature, salt_m=salt,
        )

    def _compute_concentration(self, results, temperature, salt, conc_nm):
        dimers = [r for r in results if r.kind != "hairpin"]
        seen: set[tuple[str, str]] = set()
        unique = []
        for result in sorted(dimers, key=lambda item: item.mfe_kcal_mol):
            key = tuple(sorted((result.primer_a.name, result.primer_b.name)))
            if key in seen:
                continue
            seen.add(key)
            unique.append(result)
        self_pairs = [r for r in unique if r.primer_a.name == r.primer_b.name]
        cross_pairs = [r for r in unique if r.primer_a.name != r.primer_b.name]
        chosen = self_pairs + cross_pairs[:MAX_CONCENTRATION_PAIRS]
        conc_m = conc_nm * 1e-9
        rows = []
        strongest = None
        for result in chosen:
            energies = cofold_energies(
                result.primer_a, result.primer_b,
                temperature_c=temperature, salt_m=salt,
            )
            eq = species_equilibrium(
                energies, temperature_c=temperature, conc_a=conc_m, conc_b=conc_m,
            )
            kind = "自二聚体" if result.primer_a.name == result.primer_b.name else "交叉二聚体"
            rows.append((result, kind, eq))
            dimerized = eq.ab + eq.aa + eq.bb
            label = f"{result.primer_a.name} × {result.primer_b.name}"
            if strongest is None or dimerized > strongest[2]:
                strongest = (label, eq, dimerized)
        return rows, strongest, conc_nm

    def _poll_worker(self):
        try:
            while True:
                kind, payload = self._queue.get_nowait()
                if kind == "progress":
                    self.status_var.set(payload)
                elif kind == "results":
                    _primers, results, _t, _s = payload
                    self._results = results
                    self._populate_result_tabs(results)
                elif kind == "conc":
                    self._show_concentration(payload)
                elif kind == "done":
                    self.status_var.set(payload)
                    self.run_btn.configure(state="normal")
                    self._worker = None
                    return
                elif kind == "error":
                    self.status_var.set("✘ 分析失败")
                    self.run_btn.configure(state="normal")
                    self._worker = None
                    messagebox.showerror("分析失败", payload)
                    return
        except queue.Empty:
            pass
        self.root.after(120, self._poll_worker)

    # ------------------------------------------------------------------
    # 结果页填充
    # ------------------------------------------------------------------
    def _populate_result_tabs(self, results: list[FoldResult]):
        self._fill_qc_tab()
        self._fill_heatmap_tab(results)
        for kind, _ in KIND_TABS:
            group = [r for r in results if r.kind == kind]
            self._kind_results[kind] = group
            combo = self.kind_combo[kind]
            combo.configure(values=[r.label for r in group])
            if group:
                combo.current(0)
                self._show_kind_result(kind)
            else:
                self._set_text(self.kind_detail[kind], "(该类型未勾选或无结果)")
                self.kind_fig[kind].configure(image="", text="(暂无结果)")
                self._photo = None
                self.kind_save_btn[kind].configure(state="disabled")

    def _fill_qc_tab(self):
        rows = qc_rows(self._primers, na_m=self._salt_m)
        headers = ["引物", "长度", "GC(%)", "Tm(°C)", "分子量(Da)", "3′GCclamp", "同聚串"]
        widths = [max(14, max(len(str(r["引物"])) for r in rows)) if rows else 14, 6, 7, 8, 12, 10, 8]
        lines = ["  ".join(h.center(w) for h, w in zip(headers, widths))]
        lines.append("  ".join("-" * w for w in widths))
        for row in rows:
            cells = [
                str(row["引物"]), str(row["长度 (nt)"]), f"{row['GC (%)']}",
                f"{row['Tm (°C)']}", f"{row['分子量 (Da)']}",
                str(row["3′GC clamp"]), str(row["最长同聚串"]),
            ]
            lines.append("  ".join(c.center(w) for c, w in zip(cells, widths)))
        lines += [
            "",
            "说明: Tm = SantaLucia 1998 最近邻法 (按当前单价盐校正, Ct=250 nM,",
            "非自互补假设), 不含 Mg²⁺/dNTP; 参考范围 GC 40–60%、clamp 1–3、",
            "同聚串 ≥4 提示滑移风险。均为描述性指标, 非合格判定。",
        ]
        self._set_text(self.qc_txt, "\n".join(lines))

    def _fill_heatmap_tab(self, results: list[FoldResult]):
        dimers = [r for r in results if r.kind != "hairpin"]
        if not dimers:
            self._set_text(self.heat_txt, "(未勾选二聚体分析, 无互作热图)")
            return
        labels: list[str] = []
        for r in dimers:
            for p in (r.primer_a, r.primer_b):
                if p.name not in labels:
                    labels.append(p.name)
        lookup = {}
        for r in dimers:
            lookup[tuple(sorted((r.primer_a.name, r.primer_b.name)))] = r
        col_w = max(10, max(len(n) for n in labels) + 2)
        lines = [" " * col_w + "".join(n.rjust(col_w) for n in labels)]
        for row_name in labels:
            cells = []
            for col_name in labels:
                result = lookup.get(tuple(sorted((row_name, col_name))))
                if result is None:
                    cells.append("—".rjust(col_w))
                elif row_name == col_name:
                    cells.append(f"[{result.mfe_kcal_mol:.1f}]".rjust(col_w))
                else:
                    cells.append(f"{result.mfe_kcal_mol:.1f}".rjust(col_w))
            lines.append(row_name.ljust(col_w) + "".join(cells))
        lines += [
            "",
            "图例: [ ] 对角线 = 自二聚体; 数值 = 链间 MFE (kcal/mol),",
            "越负表示结合预测越强; — = 该组合未分析。引物较多时建议",
            "优先关注最负的交叉组合是否占用双方 3′末端。",
        ]
        self._set_text(self.heat_txt, "\n".join(lines))

    def _show_kind_result(self, kind: str):
        group = self._kind_results.get(kind, [])
        if not group:
            return
        label = self.kind_combo[kind].get()
        result = next((r for r in group if r.label == label), group[0])

        lines = [
            f"MFE = {result.mfe_kcal_mol:.2f} kcal/mol   |   "
            f"总碱基对 {result.base_pairs}   |   分子间 {result.intermolecular_pairs}",
            f"3′末 5 nt 配对: A {result.a_last5_paired}/5"
            + (f" , B {result.b_last5_paired}/5" if result.b_last5_paired is not None else "")
            + f"   |   3′末端: A {'配对' if result.a_terminal_paired else '游离'}"
            + (f" / B {'配对' if result.b_terminal_paired else '游离'}"
               if result.b_terminal_paired is not None else ""),
            "",
            "结果解读 (描述性, 非合格判定):",
        ]
        lines += [f"  • {b}" for b in interpret_result(result)]
        lines += ["", f"序列: {result.folded_sequence}", f"结构: {result.structure}", ""]
        self._set_text(self.kind_detail[kind], "\n".join(lines))

        self._render_figure(kind, result)

    def _render_figure(self, kind: str, result: FoldResult):
        fig_label = self.kind_fig[kind]
        if result.base_pairs == 0:
            fig_label.configure(image="", text="(该结构无碱基对, 结构图省略)")
            self._photo = None
            self.kind_save_btn[kind].configure(state="disabled")
            return
        svg = self._svg_cache.get(result.label)
        if svg is None:
            try:
                labels = [result.primer_a.name]
                if result.primer_b:
                    labels.append(result.primer_b.name)
                svg = render_styled_structure_svg(
                    result.folded_sequence, result.structure, strand_labels=labels,
                )
            except FoldExecutionError as exc:
                fig_label.configure(image="", text=f"⚠ 结构图生成失败: {exc}")
                self._photo = None
                return
            self._svg_cache[result.label] = svg
        png_path = self._svg_to_png(svg)
        if png_path is None:
            fig_label.configure(
                image="", text="⚠ 结构图预览不可用 (可用上方按钮导出 SVG)"
            )
            self._photo = None
            self.kind_save_btn[kind].configure(state="normal")
            return
        photo = tk.PhotoImage(file=png_path)
        self._photo = photo
        fig_label.configure(image=photo, text="")
        self.kind_save_btn[kind].configure(state="normal")

    def _svg_to_png(self, svg: str) -> str | None:
        qlmanage = shutil.which("qlmanage")
        if qlmanage is None:
            return None
        width, height = svg_pixel_dimensions(svg)
        size = int(min(1400, max(700, max(width, height))))
        svg_path = os.path.join(self._fig_dir, "figure.svg")
        with open(svg_path, "w", encoding="utf-8") as f:
            f.write(svg)
        try:
            subprocess.run(
                [qlmanage, "-t", "-s", str(size), "-o", self._fig_dir, svg_path],
                capture_output=True, timeout=30, check=False,
            )
        except (subprocess.TimeoutExpired, OSError):
            return None
        png_path = svg_path + ".png"
        return png_path if os.path.isfile(png_path) else None

    def _show_concentration(self, payload):
        rows, strongest, conc_nm = payload
        if not rows:
            self._set_text(self.conc_txt, "(无二聚体结果, 无法计算浓度平衡)")
            return
        headers = ["组合", "类型", "AB(nM)", "AA(nM)", "BB(nM)", "A游离(%)", "B游离(%)", "A占用(%)"]
        widths = [max(20, max(len(f"{r.primer_a.name} × {r.primer_b.name}")
                              for r, _k, _e in rows)) + 2, 10, 9, 9, 9, 9, 9, 9]
        lines = [f"引物浓度平衡 (每条引物初始 {conc_nm} nM, RNAcofold 集合自由能)", ""]
        lines.append("  ".join(h.rjust(w) if i else h.ljust(w)
                               for i, (h, w) in enumerate(zip(headers, widths))))
        lines.append("  ".join("-" * w for w in widths))
        for result, kind, eq in rows:
            cells = [
                f"{result.primer_a.name} × {result.primer_b.name}", kind,
                f"{eq.ab * 1e9:.2f}", f"{eq.aa * 1e9:.2f}", f"{eq.bb * 1e9:.2f}",
                f"{eq.free_a / eq.conc_a * 100:.1f}", f"{eq.free_b / eq.conc_b * 100:.1f}",
                f"{eq.a_occupied_fraction * 100:.1f}",
            ]
            lines.append("  ".join(c.rjust(w) if i else c.ljust(w)
                                   for i, (c, w) in enumerate(zip(cells, widths))))
        lines += [
            "",
            "列说明: AB/AA/BB = 平衡时各二聚体物种浓度; A占用(%) = (AB+2·AA)/[A]₀,",
            "即可被聚合酶利用的 A 引物损失比例 (描述性指标)。",
        ]
        if strongest is not None and strongest[1] is not None:
            label, eq, _ = strongest
            lines += ["", f"── 浓度解读: {label} 的二聚体化程度最高 ──"]
            lines += [f"  • {b}" for b in interpret_equilibrium(eq, conc_nm=conc_nm)]
        self._set_text(self.conc_txt, "\n".join(lines))
        self.nb.select(self.nb.tabs().index(str(self.nb.tabs()[-2])))

    # ------------------------------------------------------------------
    # 导出
    # ------------------------------------------------------------------
    def save_structure_svg(self, kind: str):
        group = self._kind_results.get(kind, [])
        label = self.kind_combo[kind].get()
        result = next((r for r in group if r.label == label), None)
        if result is None:
            return
        svg = self._svg_cache.get(result.label)
        if svg is None:
            messagebox.showinfo("提示", "当前结果暂无结构图。")
            return
        path = filedialog.asksaveasfilename(
            title="保存结构图", defaultextension=".svg",
            initialfile=f"{result.primer_a.name}_{result.kind}.svg",
            filetypes=[("SVG", "*.svg"), ("所有文件", "*.*")])
        if not path:
            return
        try:
            with open(path, "w", encoding="utf-8") as f:
                f.write(svg)
            self.status_var.set(f"已保存: {path}")
        except OSError as e:
            messagebox.showerror("保存失败", str(e))

    def save_report(self):
        if not self._results:
            messagebox.showinfo("提示", "请先运行一次分析。")
            return
        figures: dict[str, str] = {}
        for result in self._results:
            if result.base_pairs == 0:
                continue
            svg = self._svg_cache.get(result.label)
            if svg is None:
                labels = [result.primer_a.name]
                if result.primer_b:
                    labels.append(result.primer_b.name)
                try:
                    svg = render_styled_structure_svg(
                        result.folded_sequence, result.structure, strand_labels=labels,
                    )
                except FoldExecutionError:
                    continue
                self._svg_cache[result.label] = svg
            figures[result.label] = svg
        report = build_html_report(
            primers=self._primers,
            results=self._results,
            figures=figures,
            temperature_c=self._temperature_c,
            salt_m=self._salt_m,
            tool_version=__version__,
        )
        stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        path = filedialog.asksaveasfilename(
            title="保存 HTML 报告", defaultextension=".html",
            initialfile=f"primerfold_report_{stamp}.html",
            filetypes=[("HTML", "*.html"), ("所有文件", "*.*")])
        if not path:
            return
        try:
            with open(path, "w", encoding="utf-8") as f:
                f.write(report)
            self.status_var.set(f"已保存报告: {path}")
        except OSError as e:
            messagebox.showerror("保存失败", str(e))

    # ------------------------------------------------------------------
    def _cleanup_figs(self):
        shutil.rmtree(self._fig_dir, ignore_errors=True)

    def _on_close(self):
        self._cleanup_figs()
        self.root.destroy()


def main() -> None:
    root = tk.Tk()
    App(root)
    root.mainloop()


if __name__ == "__main__":
    main()
