"""Local interactive GUI for DNA-primer secondary-structure analysis."""

from __future__ import annotations

import html as html_module

import pandas as pd
import streamlit as st

from primerfold import (
    AnalysisJob,
    FoldExecutionError,
    FoldResult,
    PrimerInputError,
    alignment_html,
    alignment_legend_html,
    build_jobs,
    fold_dimer,
    fold_hairpin,
    interpret_result,
    parse_primers,
    render_styled_structure_svg,
    render_structure_svg,
    svg_pixel_dimensions,
)


DEFAULT_INPUT = """>Primer_F
ATGACCATGATTACGCCAAG
>Primer_R
GCGCGCTTTTTGCGCGC
"""

KIND_LABELS = {
    "hairpin": "发卡结构",
    "self_dimer": "自二聚体",
    "cross_dimer": "交叉二聚体",
}
KIND_ICONS = {
    "hairpin": "🪝",
    "self_dimer": "🪞",
    "cross_dimer": "🔗",
}
KIND_DESCRIPTIONS = {
    "hairpin": "单条引物内部互补配对形成的茎环结构；3′端参与发卡会直接占用延伸起点。",
    "self_dimer": "一条引物与自身另一拷贝之间的链间配对，高引物浓度下更易发生。",
    "cross_dimer": "不同引物之间的链间配对；占用双方 3′末端时最易引发引物二聚体。",
}

APP_CSS = """
<style>
  .block-container {padding-top: 1.6rem; padding-bottom: 3rem; max-width: 1500px;}
  h1 {letter-spacing: -0.02em;}
  /* Metric cards */
  [data-testid="stMetric"] {
    background: #FFFFFF;
    border: 1px solid #E3EAE6;
    border-radius: 14px;
    padding: 14px 18px;
    box-shadow: 0 1px 2px rgba(23, 32, 38, 0.05);
  }
  [data-testid="stMetricLabel"] p {font-size: .86rem; color: #5C6B73; font-weight: 600;}
  [data-testid="stMetricValue"] {font-variant-numeric: tabular-nums;}
  /* Tabs */
  [data-baseweb="tab-list"] {gap: 2px; border-bottom: 2px solid #E3EAE6;}
  [data-baseweb="tab"] {padding: 10px 16px; font-weight: 600; color: #5C6B73;}
  [data-baseweb="tab"][aria-selected="true"] {color: #087F5B;}
  [data-baseweb="tab-highlight"] {background-color: #087F5B; height: 3px;}
  /* Buttons */
  .stButton > button, [data-testid="stBaseButton"] {border-radius: 10px;}
  /* Sidebar */
  [data-testid="stSidebar"] {
    background: linear-gradient(180deg, #F3F8F5 0%, #FBFDFC 100%);
  }
  [data-testid="stSidebar"] hr {margin: .8rem 0;}
  /* Data frame */
  [data-testid="stDataFrame"] {border-radius: 12px; overflow: hidden;}
  /* Expander */
  [data-testid="stExpander"] details {border-radius: 12px;}
  /* Custom widgets */
  .pf-kicker {
    color: #087F5B; font-weight: 700; letter-spacing: .1em;
    text-transform: uppercase; font-size: .78rem; margin-bottom: .1rem;
  }
  .pf-hero-sub {color: #52606D; font-size: 1.02rem; margin-top: -.5rem;}
  .pf-chip-row {display: flex; flex-wrap: wrap; gap: 6px; margin-top: 2px;}
  .pf-chip {
    background: #E9F5EF; color: #086952; border: 1px solid #BFE0D0;
    border-radius: 999px; padding: 2px 10px; font-size: .8rem; font-weight: 600;
  }
  .pf-chip-neutral {
    background: #F1F4F6; color: #52606D; border: 1px solid #DBE3E8;
    border-radius: 999px; padding: 2px 10px; font-size: .78rem;
  }
  .pf-chip-warn {
    background: #FDF0E7; color: #B2400F; border: 1px solid #F5CBA8;
    border-radius: 999px; padding: 2px 10px; font-size: .8rem;
  }
  .pf-interpret {
    background: #F4F9F6; border: 1px solid #DCEAE2; border-left: 4px solid #087F5B;
    border-radius: 10px; padding: 12px 16px; margin: 4px 0 8px 0;
  }
  .pf-interpret-title {font-weight: 700; color: #086952; margin-bottom: 4px;}
  .pf-interpret ul {margin: 0 0 0 1.1rem; padding: 0;}
  .pf-interpret li {margin: 3px 0; color: #33414B; font-size: .92rem;}
  .pf-section-h {
    font-size: .78rem; font-weight: 700; color: #087F5B;
    letter-spacing: .06em; margin: 1rem 0 .2rem 0;
  }
  .pf-caption {color: #6B7A85; font-size: .8rem; margin-top: -.4rem;}
</style>
"""

GLOSSARY_ROWS = [
    (
        "MFE（kcal/mol）",
        "模型预测的最稳定结构的最小自由能。",
        "越负表示结构越稳定；仅适合同一引物在不同条件间的相对比较，不是实验测量值。",
    ),
    (
        "总碱基对",
        "结构中参与配对的碱基总数。",
        "发卡结果全部为链内配对；二聚体结果全部为链间配对。",
    ),
    (
        "分子间碱基对",
        "两条链之间形成的配对数。",
        "发卡恒为 0；二聚体中数值越大表示链间结合越强。",
    ),
    (
        "3′末 5 nt 配对",
        "3′端最后 5 个核苷酸中参与配对的数目。",
        "3′端是 DNA 聚合酶延伸的起点；末端被占用可能降低延伸效率。发卡统计链内配对，二聚体只统计链间配对。",
    ),
    (
        "3′末端配对",
        "3′最末位碱基是否处于配对状态。",
        "“配对”提示延伸起点被占用；“游离”表示末端可正常延伸。",
    ),
    (
        "点括号结构",
        "结构的文本表示。",
        "成对的 ( ) 表示配对，. 表示未配对，& 分隔两条链（各自均为 5′→3′）。",
    ),
    (
        "温度 / 单价盐",
        "热力学参数的重标定条件。",
        "仅含单价盐（Na⁺ 等效）校正；模型不显式表示 Mg²⁺、dNTP 与引物浓度。",
    ),
]


st.set_page_config(
    page_title="PrimerFold · DNA 引物二级结构",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="auto",
)

st.markdown(APP_CSS, unsafe_allow_html=True)


@st.cache_data(show_spinner=False)
def cached_styled_svg(sequence: str, structure: str, strand_labels: tuple[str, ...]) -> str:
    return render_styled_structure_svg(sequence, structure, strand_labels=list(strand_labels))


@st.cache_data(show_spinner=False)
def cached_raw_svg(sequence: str, structure: str) -> str:
    return render_structure_svg(sequence, structure)


def execute_job(job: AnalysisJob, *, temperature_c: float, salt_m: float) -> FoldResult:
    if job.kind == "hairpin":
        return fold_hairpin(
            job.primer_a,
            temperature_c=temperature_c,
            salt_m=salt_m,
        )
    if job.primer_b is None:
        raise FoldExecutionError("二聚体任务缺少第二条引物。")
    return fold_dimer(
        job.primer_a,
        job.primer_b,
        kind=job.kind,
        temperature_c=temperature_c,
        salt_m=salt_m,
    )


def result_to_row(
    result: FoldResult,
    *,
    temperature_c: float,
    salt_m: float,
) -> dict[str, object]:
    primer_b = result.primer_b
    return {
        "对象": result.label,
        "分析类型": KIND_LABELS[result.kind],
        "引物 A": result.primer_a.name,
        "序列 A (5′→3′)": result.primer_a.sequence,
        "长度 A (nt)": len(result.primer_a.sequence),
        "引物 B": primer_b.name if primer_b else "",
        "序列 B (5′→3′)": primer_b.sequence if primer_b else "",
        "长度 B (nt)": len(primer_b.sequence) if primer_b else "",
        "MFE (kcal/mol)": result.mfe_kcal_mol,
        "点括号结构": result.structure,
        "总碱基对": result.base_pairs,
        "分子间碱基对": result.intermolecular_pairs,
        "A 末5nt 配对": result.a_last5_paired,
        "B 末5nt 配对": result.b_last5_paired if primer_b else "",
        "A 3′末端": "配对" if result.a_terminal_paired else "游离",
        "B 3′末端": (
            "配对" if result.b_terminal_paired else "游离"
        ) if primer_b else "",
        "温度 (°C)": temperature_c,
        "单价盐 (M)": salt_m,
        "ViennaRNA 提示": " | ".join(result.warnings),
    }


def make_dataframe(
    results: list[FoldResult],
    *,
    temperature_c: float,
    salt_m: float,
) -> pd.DataFrame:
    return pd.DataFrame(
        [
            result_to_row(
                result,
                temperature_c=temperature_c,
                salt_m=salt_m,
            )
            for result in results
        ]
    )


COLUMN_HELP = {
    "对象": "分析对象：发卡与自二聚体为单条引物，交叉二聚体为两条引物的组合。",
    "MFE (kcal/mol)": "最低自由能。越负表示模型认为该结构越稳定；仅用于相对比较，非实验值。",
    "总碱基对": "结构中配对碱基总数：发卡为链内配对，二聚体为链间配对。",
    "A 末5nt 配对": "引物 A 的 3′端最后 5 nt 中参与配对的数目（发卡计链内，二聚体只计链间）。",
    "B 末5nt 配对": "引物 B 的 3′端最后 5 nt 中参与链间配对的数目。",
    "A 3′末端": "引物 A 的 3′最末位碱基是否处于配对状态。“配对”提示延伸起点被占用。",
    "B 3′末端": "引物 B 的 3′最末位碱基是否处于配对状态。",
    "点括号结构": "成对 ( ) = 配对，. = 未配对，& = 两条链分隔符。",
}

TABLE_COLUMNS = {
    "hairpin": [
        "对象",
        "MFE (kcal/mol)",
        "总碱基对",
        "A 末5nt 配对",
        "A 3′末端",
        "点括号结构",
    ],
    "self_dimer": [
        "对象",
        "MFE (kcal/mol)",
        "总碱基对",
        "A 末5nt 配对",
        "B 末5nt 配对",
        "A 3′末端",
        "B 3′末端",
        "点括号结构",
    ],
}
TABLE_COLUMNS["cross_dimer"] = TABLE_COLUMNS["self_dimer"]


def _column_config(name: str):
    if name == "MFE (kcal/mol)":
        return st.column_config.NumberColumn(format="%.2f", help=COLUMN_HELP[name])
    if name in ("A 末5nt 配对", "B 末5nt 配对"):
        return st.column_config.ProgressColumn(
            format="%.0f / 5",
            min_value=0,
            max_value=5,
            help=COLUMN_HELP[name],
        )
    if name == "点括号结构":
        return st.column_config.TextColumn(width="large", help=COLUMN_HELP[name])
    return st.column_config.Column(help=COLUMN_HELP[name])


def render_table(results: list[FoldResult], kind: str, temperature_c: float, salt_m: float) -> None:
    dataframe = make_dataframe(
        results,
        temperature_c=temperature_c,
        salt_m=salt_m,
    )
    visible = dataframe[TABLE_COLUMNS[kind]]
    st.dataframe(
        visible,
        hide_index=True,
        width="stretch",
        column_config={name: _column_config(name) for name in TABLE_COLUMNS[kind]},
    )
    st.caption(
        "💡 将鼠标悬停在列名旁的 ⓘ 图标上可查看每个指标的详细解释；"
        "进度条表示 3′端最后 5 nt 中的配对数（0–5）。"
    )


def render_metric_with_caption(
    column,
    label: str,
    value: str,
    *,
    help_text: str,
    caption: str,
) -> None:
    column.metric(label, value, help=help_text)
    column.markdown(f'<div class="pf-caption">{caption}</div>', unsafe_allow_html=True)


def render_interpretation(selected: FoldResult) -> None:
    bullets = interpret_result(selected)
    items = "".join(
        f"<li>{html_module.escape(bullet)}</li>" for bullet in bullets
    )
    st.markdown(
        '<div class="pf-interpret"><div class="pf-interpret-title">📖 结果解读'
        "（描述性说明，非合格判定）</div>"
        f"<ul>{items}</ul></div>",
        unsafe_allow_html=True,
    )


def render_alignment(selected: FoldResult) -> None:
    st.caption(
        "序列与结构逐位对照（每 10 nt 一组，组上方数字为该组起始位置；"
        "两条链各自均按 5′→3′ 方向显示）"
    )
    st.markdown(
        alignment_html(selected.folded_sequence, selected.structure),
        unsafe_allow_html=True,
    )
    st.markdown(alignment_legend_html(), unsafe_allow_html=True)


def render_figure(selected: FoldResult, key_prefix: str) -> None:
    labels = [selected.primer_a.name]
    if selected.primer_b:
        labels.append(selected.primer_b.name)

    styled = True
    try:
        svg = cached_styled_svg(selected.folded_sequence, selected.structure, tuple(labels))
    except FoldExecutionError:
        styled = False
        svg = cached_raw_svg(selected.folded_sequence, selected.structure)

    if not styled:
        st.caption("⚠️ 自定义美化渲染失败，已回退到 RNAplot 原始结构图。")
    else:
        st.caption(
            "结构图说明：实心圆＝参与配对的碱基，空心圆＝未配对，青色短线＝G≡C、"
            "灰色短线＝A=T，紫色虚线圈＝3′端最后 5 nt。"
        )

    _, svg_height = svg_pixel_dimensions(svg)
    embed_height = int(min(760, max(380, svg_height)))
    svg_start = svg.find("<svg")
    svg_fragment = svg[svg_start:] if svg_start >= 0 else svg
    st.iframe(
        """
        <div style="display:flex;justify-content:center;align-items:flex-start;"
             aria-label="二级结构图">
        """
        + svg_fragment
        + "</div>",
        height=embed_height,
        tab_index=-1,
    )
    st.download_button(
        "下载当前结构图（SVG）",
        data=svg.encode("utf-8"),
        file_name=f"{selected.primer_a.name}_{selected.kind}.svg",
        mime="image/svg+xml",
        key=f"{key_prefix}_svg_download",
        on_click="ignore",
    )


def render_structure_viewer(results: list[FoldResult], key_prefix: str) -> None:
    if not results:
        return
    st.markdown("#### 结构详情")
    chosen_label = st.selectbox(
        "选择结果",
        options=[result.label for result in results],
        key=f"{key_prefix}_structure_choice",
    )
    selected = next(result for result in results if result.label == chosen_label)

    metric_columns = st.columns(4)
    render_metric_with_caption(
        metric_columns[0],
        "MFE (kcal/mol)",
        f"{selected.mfe_kcal_mol:.2f}",
        help_text="最低自由能：模型预测的最稳定结构的能量，越负越稳定。",
        caption="自由能，越负越稳定；模型预测值，非实验测量。",
    )
    render_metric_with_caption(
        metric_columns[1],
        "总碱基对",
        str(selected.base_pairs),
        help_text="结构中参与配对的碱基总数。",
        caption="发卡＝链内配对数；二聚体＝链间配对数。",
    )
    render_metric_with_caption(
        metric_columns[2],
        "分子间碱基对",
        str(selected.intermolecular_pairs),
        help_text="两条链之间形成的碱基对数。",
        caption="发卡恒为 0；二聚体中越大链间结合越强。",
    )
    if selected.kind == "hairpin":
        terminal_value = f"{selected.a_last5_paired} / 5"
        terminal_help = "该引物 3′端最后 5 nt 中参与发卡配对的核苷酸数。"
        terminal_caption = "3′端被占用会直接影响聚合酶延伸。"
    else:
        terminal_value = f"{selected.a_last5_paired} / {selected.b_last5_paired}"
        terminal_help = "A 链 / B 链各自 3′端最后 5 nt 中参与链间配对的核苷酸数。"
        terminal_caption = "两端 3′末端同时被占用最易引发引物二聚体。"
    render_metric_with_caption(
        metric_columns[3],
        "3′末 5 nt 配对",
        terminal_value,
        help_text=terminal_help,
        caption=terminal_caption,
    )

    render_interpretation(selected)

    info_columns = st.columns(2)
    with info_columns[0]:
        st.caption("序列（两条链以 & 分隔，均 5′→3′）")
        st.code(selected.folded_sequence, language=None, wrap_lines=True)
    with info_columns[1]:
        st.caption("点括号结构（( )＝配对，.＝未配对，&＝链分隔）")
        st.code(selected.structure, language=None, wrap_lines=True)

    st.markdown("##### 结构图与逐位对照")
    figure_column, alignment_column = st.columns([1, 1])
    with figure_column:
        render_figure(selected, key_prefix)
    with alignment_column:
        render_alignment(selected)


def render_glossary() -> None:
    with st.expander("指标说明：每个数字代表什么、如何解读"):
        rows = "".join(
            f"<tr><td style='font-weight:600;white-space:nowrap;'>{html_module.escape(metric)}</td>"
            f"<td>{html_module.escape(meaning)}</td>"
            f"<td>{html_module.escape(reading)}</td></tr>"
            for metric, meaning, reading in GLOSSARY_ROWS
        )
        st.markdown(
            "<div style='font-size:.9rem;'>"
            "<table style='width:100%;border-collapse:collapse;'>"
            "<thead><tr>"
            "<th style='text-align:left;padding:6px 10px;border-bottom:2px solid #DCEAE2;'>指标</th>"
            "<th style='text-align:left;padding:6px 10px;border-bottom:2px solid #DCEAE2;'>含义</th>"
            "<th style='text-align:left;padding:6px 10px;border-bottom:2px solid #DCEAE2;'>如何解读</th>"
            "</tr></thead>"
            f"<tbody>{rows}</tbody></table></div>",
            unsafe_allow_html=True,
        )
        st.caption(
            "本工具只报告描述性指标，不使用通用阈值判定引物“合格/不合格”；"
            "请结合具体实验体系判断。"
        )


def render_results(payload: dict[str, object]) -> None:
    primers = payload["primers"]
    results = payload["results"]
    temperature_c = float(payload["temperature_c"])
    salt_m = float(payload["salt_m"])
    assert isinstance(primers, list)
    assert isinstance(results, list)

    st.divider()
    st.subheader("分析结果")
    summary_columns = st.columns(4)
    summary_columns[0].metric(
        "引物数",
        len(primers),
        help="本次解析并参与分析的引物条数。",
    )
    summary_columns[1].metric(
        "分析任务",
        len(results),
        help="依据勾选的分析类型生成的折叠任务总数：发卡 N + 自二聚体 N + 交叉二聚体 C(N,2)。",
    )
    summary_columns[2].metric(
        "温度",
        f"{temperature_c:g} °C",
        help="ViennaRNA 热力学参数重标定温度，可在左侧边栏修改。",
    )
    summary_columns[3].metric(
        "单价盐",
        f"{salt_m * 1000:g} mM",
        help="Na⁺ 等效单价盐浓度；不包含 Mg²⁺ 与 dNTP。",
    )

    st.info(
        "**如何阅读这些结果**：MFE 越负表示该条件下预测结构越稳定；"
        "对 PCR 影响最大的是 **3′末端被配对占用**（见“3′末 5 nt 配对”与“3′末端”列）。"
        "所有数字都是热力学模型推断，不是实验测量。点击下方页签查看各类结构的"
        "数据表与逐条解读。"
    )
    render_glossary()

    grouped = {
        kind: [result for result in results if result.kind == kind]
        for kind in ("hairpin", "self_dimer", "cross_dimer")
    }
    nonempty_kinds = [kind for kind, values in grouped.items() if values]
    tabs = st.tabs(
        [
            f"{KIND_ICONS[kind]} {KIND_LABELS[kind]}（{len(grouped[kind])}）"
            for kind in nonempty_kinds
        ]
    )
    for tab, kind in zip(tabs, nonempty_kinds, strict=True):
        with tab:
            st.caption(KIND_DESCRIPTIONS[kind])
            render_table(grouped[kind], kind, temperature_c, salt_m)
            render_structure_viewer(grouped[kind], kind)

    all_rows = make_dataframe(
        results,
        temperature_c=temperature_c,
        salt_m=salt_m,
    )
    csv_bytes = all_rows.to_csv(index=False).encode("utf-8-sig")
    st.download_button(
        "下载全部结果（CSV，含全部列与解释性数值）",
        data=csv_bytes,
        file_name="primerfold_results.csv",
        mime="text/csv",
        type="primary",
        on_click="ignore",
    )

    unique_warnings = list(
        dict.fromkeys(
            warning
            for result in results
            for warning in result.warnings
        )
    )
    if unique_warnings:
        with st.expander("ViennaRNA 运行提示（stderr 输出，通常无需处理）"):
            st.code("\n".join(unique_warnings), language=None)


def render_primer_preview(source_text: str) -> None:
    """Live parse preview shown as chips under the input area."""

    if not source_text.strip():
        return
    try:
        primers = parse_primers(source_text)
    except PrimerInputError as exc:
        st.markdown(
            f'<span class="pf-chip-warn">⚠️ {html_module.escape(str(exc))}</span>',
            unsafe_allow_html=True,
        )
        return
    chips = [
        f'<span class="pf-chip">{html_module.escape(primer.name)} · '
        f"{len(primer.sequence)} nt</span>"
        for primer in primers[:8]
    ]
    if len(primers) > 8:
        chips.append(f'<span class="pf-chip-neutral">… 共 {len(primers)} 条</span>')
    else:
        chips.append(f'<span class="pf-chip-neutral">共 {len(primers)} 条，解析正常</span>')
    st.markdown(f'<div class="pf-chip-row">{"".join(chips)}</div>', unsafe_allow_html=True)


with st.sidebar:
    st.header("🧬 PrimerFold")
    st.caption("DNA 引物发卡与二聚体分析 · 全程本地运行，序列不出机。")

    st.markdown('<div class="pf-section-h">① 热力学条件</div>', unsafe_allow_html=True)
    temperature_c = st.number_input(
        "温度（°C）",
        min_value=0.0,
        max_value=100.0,
        value=37.0,
        step=0.5,
        key="pf_temperature",
        help="ViennaRNA 用该温度重标定热力学参数；37 °C 便于与常见 ΔG 报告比较。",
    )
    st.markdown('<div class="pf-caption">退火或延伸温度通常取 37–60 °C。</div>', unsafe_allow_html=True)
    salt_m = st.number_input(
        "单价盐浓度（M）",
        min_value=0.0001,
        max_value=2.0,
        value=0.05,
        step=0.01,
        format="%.4f",
        key="pf_salt",
        help="Na⁺ 等效单价盐浓度，例如 50 mM 输入 0.05。该参数不等同于 Mg²⁺ 浓度。",
    )
    st.markdown('<div class="pf-caption">50 mM → 填 0.05；模型不含 Mg²⁺ 与 dNTP 校正。</div>', unsafe_allow_html=True)

    st.markdown('<div class="pf-section-h">② 分析类型</div>', unsafe_allow_html=True)
    analyze_hairpin = st.checkbox(
        "发卡结构",
        value=True,
        help="单条引物内部互补配对形成的茎环结构。",
    )
    analyze_self_dimer = st.checkbox(
        "自二聚体",
        value=True,
        help="一条引物与自身另一拷贝之间的链间配对。",
    )
    analyze_cross_dimer = st.checkbox(
        "引物间交叉二聚体",
        value=True,
        help="不同引物之间的链间配对；占用双方 3′末端时最危险。",
    )

    with st.expander("模型与边界"):
        st.markdown(
            "- ViennaRNA 2.7.2\n"
            "- Mathews 2004 DNA 参数（`--paramFile=DNA`）\n"
            "- 保留 T，不转换为 U（`--noconv`）\n"
            "- 发卡使用 RNAfold；二聚体使用仅允许链间配对的 RNAduplex\n"
            "- MFE 模型不显式表示 Mg²⁺、dNTP 或引物浓度\n"
            "- GUI 只报告描述性指标，不用通用阈值判定引物合格/不合格"
        )
    st.caption(
        "结构图布局由 RNAplot 计算，界面仅重新着色与标注，不改变化学预测结果。"
    )

st.markdown('<div class="pf-kicker">ViennaRNA local GUI</div>', unsafe_allow_html=True)
st.title("PrimerFold")
st.markdown(
    '<p class="pf-hero-sub">DNA 引物发卡与二聚体分析 · 每项指标均附解释 · 本地运行，不上传序列</p>',
    unsafe_allow_html=True,
)
st.markdown(
    '<div class="pf-chip-row">'
    '<span class="pf-chip">🖥️ 本地分析 · 序列不出机</span>'
    '<span class="pf-chip">⚙️ ViennaRNA 2.7.2 · DNA 参数</span>'
    '<span class="pf-chip">📖 描述性指标 · 不做合格判定</span>'
    "</div>",
    unsafe_allow_html=True,
)

with st.container(border=True):
    input_mode = st.radio(
        "输入方式",
        ("粘贴序列", "上传文件"),
        horizontal=True,
        help="支持 FASTA、CSV/TSV、名称+序列或每行一条纯序列。",
    )
    source_text = ""
    if input_mode == "粘贴序列":
        source_text = st.text_area(
            "引物序列",
            value=DEFAULT_INPUT,
            height=210,
            help="支持 FASTA（>名称 换行 序列）、`名称,序列` 或每行一条纯序列；方向均为 5′→3′。",
        )
    else:
        uploaded = st.file_uploader(
            "上传 FASTA / CSV / TSV / TXT",
            type=("fa", "fasta", "fna", "csv", "tsv", "txt"),
            help="文件需为 UTF-8 编码文本，最大 5 MB。",
        )
        if uploaded is not None:
            try:
                source_text = uploaded.getvalue().decode("utf-8-sig")
            except UnicodeDecodeError:
                st.error("文件不是 UTF-8 文本，请转换编码后重试。")

    render_primer_preview(source_text)

    with st.expander("支持的输入格式示例"):
        st.markdown(
            "**FASTA**（推荐，可自带引物名称）\n"
            "```text\n>Primer_F\nATGACCATGATTACGCCAAG\n>Primer_R\nGCGCGCTTTTTGCGCGC\n```\n"
            "**CSV / 名称+序列**\n"
            "```text\nPrimer_F,ATGACCATGATTACGCCAAG\n```\n"
            "**纯序列**（每行一条，自动命名 Primer_1、Primer_2…）\n"
            "```text\nATGACCATGATTACGCCAAG\n```"
        )
        st.caption("仅接受 A/C/G/T；一次最多 30 条引物，每条 4–500 nt。")

    analyze_clicked = st.button(
        "开始分析",
        type="primary",
        width="stretch",
    )

if analyze_clicked:
    st.session_state.pop("primerfold_payload", None)
    try:
        primers = parse_primers(source_text)
        jobs = build_jobs(
            primers,
            hairpin=analyze_hairpin,
            self_dimer=analyze_self_dimer,
            cross_dimer=analyze_cross_dimer,
        )
    except PrimerInputError as exc:
        st.error(str(exc))
    else:
        progress = st.progress(0.0, text=f"准备运行 {len(jobs)} 个分析任务…")
        results: list[FoldResult] = []
        execution_error: Exception | None = None
        for index, job in enumerate(jobs, start=1):
            progress.progress(
                (index - 1) / len(jobs),
                text=f"正在分析：{job.label}（{index}/{len(jobs)}）",
            )
            try:
                results.append(
                    execute_job(
                        job,
                        temperature_c=float(temperature_c),
                        salt_m=float(salt_m),
                    )
                )
            except (FoldExecutionError, PrimerInputError) as exc:
                execution_error = exc
                break

        if execution_error is not None:
            progress.empty()
            st.error(f"分析中止：{execution_error}")
        else:
            progress.progress(1.0, text="分析完成")
            st.session_state["primerfold_payload"] = {
                "primers": primers,
                "results": results,
                "temperature_c": float(temperature_c),
                "salt_m": float(salt_m),
            }

payload = st.session_state.get("primerfold_payload")
if isinstance(payload, dict):
    render_results(payload)

st.divider()
st.caption(
    "点括号结构：成对的 ( ) 表示碱基配对，. 表示未配对，& 表示两条链的分隔；"
    "结构图为热力学模型预测的最稳定（MFE）结构。"
)
