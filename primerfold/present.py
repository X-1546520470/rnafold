"""Pure-presentation helpers shared by the GUI: interpretation text, a
colour-coded sequence/structure alignment, an interaction heatmap and a
standalone HTML report.

Kept free of Streamlit so the logic stays unit-testable.
"""

from __future__ import annotations

import html as html_module
from datetime import datetime
from xml.sax.saxutils import escape

from .core import FoldResult, Primer, base_pair_coordinates
from .plot import BASE_COLORS
from .qc import qc_rows

PAIR_COLOR_INTER = "#D6336C"
PAIR_COLOR_DEPTH = ("#087F5B", "#0B7285", "#6741D9", "#C2255C", "#E8590C")
UNPAIRED_COLOR = "#9AA9B7"
SEQ_MUTED_COLOR = "#C1CCD4"
TERMINAL_UNDERLINE = "#7048E8"
SEPARATOR_COLOR = "#7A8B98"
MONO_STACK = "ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace"


def strand_starts(sequence: str) -> list[int]:
    """Global start index of every strand inside the &-joined sequence."""

    starts: list[int] = []
    offset = 0
    for part in sequence.split("&"):
        starts.append(offset)
        offset += len(part) + 1
    return starts


def terminal_positions(sequence: str, *, last_n: int = 5) -> set[int]:
    """Global positions of the final ``last_n`` nucleotides of every strand."""

    positions: set[int] = set()
    for start, part in zip(strand_starts(sequence), sequence.split("&"), strict=True):
        positions.update(range(start + max(0, len(part) - last_n), start + len(part)))
    return positions


def classified_pairs(
    sequence: str, structure: str
) -> list[tuple[int, int, str]]:
    """Return (global pos A, global pos B, colour) for every base pair.

    Inter-strand pairs share one magenta colour; intra-strand pairs are
    coloured by nesting depth so individual stems stay distinguishable.
    """

    starts = strand_starts(sequence)
    colours: list[tuple[int, int, str]] = []
    stack: list[tuple[int, int]] = []
    position = 0
    for character in structure:
        if character == "&":
            position += 1
            continue
        if character == "(":
            strand = max(index for index, start in enumerate(starts) if start <= position)
            stack.append((position, strand))
        elif character == ")":
            if not stack:  # pragma: no cover - validated by core first
                raise ValueError("不平衡的点括号结构")
            open_position, open_strand = stack.pop()
            strand = max(index for index, start in enumerate(starts) if start <= position)
            if open_strand != strand:
                colour = PAIR_COLOR_INTER
            else:
                colour = PAIR_COLOR_DEPTH[len(stack) % len(PAIR_COLOR_DEPTH)]
            colours.append((open_position, position, colour))
        position += 1
    return colours


def structure_char_colors(sequence: str, structure: str) -> dict[int, str]:
    """Colour for each '(' / ')' character position in the structure string."""

    colors: dict[int, str] = {}
    for first, second, colour in classified_pairs(sequence, structure):
        colors[first] = colour
        colors[second] = colour
    return colors


def interpret_result(result: FoldResult, *, last_n: int = 5) -> list[str]:
    """Human-readable, strictly descriptive interpretation of one result."""

    bullets: list[str] = []
    is_dimer = result.kind != "hairpin"

    if not is_dimer:
        if result.base_pairs == 0:
            bullets.append(
                "未预测到链内碱基对：该引物在此条件下自身折叠不稳定。"
            )
        else:
            loop_unpaired = _hairpin_loop_unpaired(result.structure)
            bullets.append(
                f"发卡（茎环）结构：共 {result.base_pairs} 个链内碱基对，"
                f"茎环内部有 {loop_unpaired} 个未配对碱基。"
            )
        bullets.append(
            f"MFE = {result.mfe_kcal_mol:.2f} kcal/mol，"
            f"属{_mfe_descriptor(result.mfe_kcal_mol)}水平"
            "（描述性说法，非合格判定；同一引物换条件才有比较意义）。"
        )
        bullets.append(
            f"3′末端：最后 {last_n} nt 中 {result.a_last5_paired} 个参与链内配对，"
            f"3′最末位碱基{'处于配对状态' if result.a_terminal_paired else '保持游离'}。"
        )
    else:
        name_a = result.primer_a.name
        name_b = result.primer_b.name if result.primer_b else "B"
        if result.intermolecular_pairs == 0:
            bullets.append(
                "未检出链间碱基对：模型预测两条链在该条件下不会形成稳定二聚体"
                "（MFE 记为 0 kcal/mol）。"
            )
            bullets.append("两条链的 3′末端均保持游离。")
        else:
            bullets.append(
                f"链间结合：{name_a} 与 {name_b} 通过 {result.intermolecular_pairs} "
                f"个碱基对结合，{_binding_region_text(result)}。"
            )
            bullets.append(
                f"MFE = {result.mfe_kcal_mol:.2f} kcal/mol，"
                f"结合强度属{_mfe_descriptor(result.mfe_kcal_mol)}水平"
                "（描述性说法，非合格判定）。"
            )
            bullets.append(
                f"3′末端（只统计链间配对）：{name_a} 最后 {last_n} nt 中 "
                f"{result.a_last5_paired} 个被占用，"
                f"{name_b} 最后 {last_n} nt 中 "
                f"{result.b_last5_paired} 个被占用；"
                f"末端状态 {'/'.join(_terminal_state(result))}。"
            )
    bullets.append(
        "3′末端是 DNA 聚合酶延伸的起点：末端若被配对占用，延伸可能受影响。"
        "以上均为热力学模型预测，请结合 Mg²⁺、dNTP、引物浓度等实际条件判断。"
    )
    return bullets


def _mfe_descriptor(mfe: float) -> str:
    strength = abs(mfe)
    if strength >= 8.0:
        return "很强"
    if strength >= 5.0:
        return "较强"
    if strength >= 2.5:
        return "中等"
    return "较弱"


def _hairpin_loop_unpaired(structure: str) -> int:
    open_index = structure.find("(")
    close_index = structure.rfind(")")
    if open_index < 0 or close_index <= open_index:
        return 0
    return structure.count(".", open_index + 1, close_index)


def _binding_region_text(result: FoldResult) -> str:
    starts = strand_starts(result.folded_sequence)
    a_positions: list[int] = []
    b_positions: list[int] = []
    for first, second, _ in classified_pairs(
        result.folded_sequence, result.structure
    ):
        a_positions.append(first - starts[0] + 1)
        b_positions.append(second - starts[1] + 1)
    a_positions.sort()
    b_positions.sort()

    def describe(positions: list[int]) -> str:
        if not positions:
            return "无配对区域"
        if len(positions) == 1:
            return f"第 {positions[0]} 位"
        if positions == list(range(positions[0], positions[0] + len(positions))):
            return f"第 {positions[0]}–{positions[-1]} 位"
        return f"共 {len(positions)} 个配对位点（第 {positions[0]}–{positions[-1]} 位）"

    return f"涉及 A 链{describe(a_positions)}与 B 链{describe(b_positions)}"


def _terminal_state(result: FoldResult) -> list[str]:
    states = []
    for terminal in (result.a_terminal_paired, result.b_terminal_paired):
        if terminal is None:
            continue
        states.append("配对" if terminal else "游离")
    return states or ["游离"]


def alignment_html(
    sequence: str,
    structure: str,
    *,
    group: int = 10,
    row_groups: int = 5,
) -> str:
    """Colour-coded sequence/structure alignment as self-contained HTML."""

    char_colors = structure_char_colors(sequence, structure)
    terminal = terminal_positions(sequence)
    paired_positions = set(char_colors)

    row_width = group * row_groups
    rows_html: list[str] = []
    for row_start in range(0, len(sequence), row_width):
        groups_html: list[str] = []
        for group_start in range(row_start, min(row_start + row_width, len(sequence)), group):
            group_end = min(group_start + group, len(sequence))
            seq_spans: list[str] = []
            struct_spans: list[str] = []
            for position in range(group_start, group_end):
                seq_char = sequence[position]
                struct_char = structure[position]
                if seq_char == "&":
                    seq_spans.append(
                        f'<span style="color:{SEPARATOR_COLOR};font-weight:700;'
                        f'padding:0 6px;">&amp;</span>'
                    )
                    struct_spans.append(
                        f'<span style="color:{SEPARATOR_COLOR};font-weight:700;'
                        f'padding:0 6px;">&amp;</span>'
                    )
                    continue
                is_terminal = position in terminal
                underline = (
                    f"border-bottom:2px solid {TERMINAL_UNDERLINE};"
                    if is_terminal
                    else ""
                )
                if position in paired_positions:
                    seq_color = BASE_COLORS.get(seq_char, "#334155")
                else:
                    seq_color = SEQ_MUTED_COLOR
                seq_spans.append(
                    f'<span style="color:{seq_color};{underline}">{escape(seq_char)}</span>'
                )
                if struct_char in "()":
                    struct_color = char_colors.get(position, "#334155")
                    struct_spans.append(
                        f'<span style="color:{struct_color};{underline}">'
                        f"{escape(struct_char)}</span>"
                    )
                else:
                    struct_spans.append(
                        f'<span style="color:{UNPAIRED_COLOR};{underline}">.</span>'
                    )
            label = str(group_start + 1)
            groups_html.append(
                '<div style="display:inline-block;margin-right:14px;'
                'vertical-align:top;">'
                f'<div style="font-size:10px;line-height:1.3;color:#8FA3AD;'
                f'font-family:{MONO_STACK};">{label}</div>'
                f'<pre style="margin:0;font-family:{MONO_STACK};font-size:14px;'
                f'line-height:1.5;">{"".join(seq_spans)}</pre>'
                f'<pre style="margin:0;font-family:{MONO_STACK};font-size:14px;'
                f'line-height:1.5;">{"".join(struct_spans)}</pre>'
                "</div>"
            )
        rows_html.append(f'<div style="margin-bottom:10px;">{"".join(groups_html)}</div>')

    container = (
        '<div style="overflow-x:auto;padding:4px 2px;">' + "".join(rows_html) + "</div>"
    )
    return container


def alignment_legend_html() -> str:
    """Legend explaining the colour code used by :func:`alignment_html`."""

    def dot(color: str, label: str) -> str:
        return (
            f'<span style="display:inline-block;width:10px;height:10px;'
            f'border-radius:50%;background:{color};margin:0 4px 0 12px;'
            f'vertical-align:middle;"></span>{label}'
        )

    def underline_sample(label: str) -> str:
        return (
            f'<span style="text-decoration:underline;text-decoration-color:'
            f'{TERMINAL_UNDERLINE};text-decoration-thickness:2px;'
            f'text-underline-offset:3px;margin-left:12px;">{label}</span>'
        )

    parts = [
        '<div style="font-size:12.5px;color:#5C6B73;margin-top:6px;">图例：'
        + dot(BASE_COLORS["A"], "A")
        + dot(BASE_COLORS["C"], "C")
        + dot(BASE_COLORS["G"], "G")
        + dot(BASE_COLORS["T"], "T")
        + dot(PAIR_COLOR_INTER, "链间配对（结构行括号）")
        + dot(PAIR_COLOR_DEPTH[0], "链内配对·按嵌套深度着色")
        + dot(SEQ_MUTED_COLOR, "未配对")
        + underline_sample("紫色下划线＝3′端最后 5 nt")
        + "</div>"
    ]
    return "".join(parts)


# ---------------------------------------------------------------------------
# Interaction heatmap
# ---------------------------------------------------------------------------

HEATMAP_MAX_ABS_MFE = 8.0
KIND_LABELS_LOCAL = {
    "hairpin": "发卡结构",
    "self_dimer": "自二聚体",
    "cross_dimer": "交叉二聚体",
}


def _heatmap_color(mfe: float) -> str:
    intensity = max(0.0, min(1.0, -mfe / HEATMAP_MAX_ABS_MFE))
    red = round(255 - 54 * intensity)
    green = round(255 - 213 * intensity)
    blue = round(255 - 213 * intensity)
    return f"rgb({red},{green},{blue})"


def heatmap_html(results: list[FoldResult]) -> str:
    """N×N dimer-interaction heatmap (self dimers on the diagonal)."""

    dimers = [result for result in results if result.kind != "hairpin"]
    if not dimers:
        return ""

    labels: list[str] = []
    for result in dimers:
        for primer in (result.primer_a, result.primer_b):
            if primer.name not in labels:
                labels.append(primer.name)

    lookup: dict[tuple[str, str], FoldResult] = {}
    for result in dimers:
        key = tuple(sorted((result.primer_a.name, result.primer_b.name)))
        lookup[key] = result

    header_cells = "".join(
        f'<th style="padding:6px 8px;border:1px solid #E3EAE6;background:#F4F9F6;'
        f'font-size:12px;">{html_module.escape(label)}</th>'
        for label in labels
    )
    rows: list[str] = []
    for row_label in labels:
        cells = [
            f'<th style="padding:6px 8px;border:1px solid #E3EAE6;background:#F4F9F6;'
            f'font-size:12px;text-align:left;">{html_module.escape(row_label)}</th>'
        ]
        for column_label in labels:
            if row_label == column_label:
                key = (row_label, row_label)
            else:
                key = tuple(sorted((row_label, column_label)))
            result = lookup.get(key)
            if result is None:
                cells.append(
                    '<td style="padding:6px 8px;border:1px solid #E3EAE6;'
                    'text-align:center;color:#B3BEC7;font-size:12px;">—</td>'
                )
                continue
            background = _heatmap_color(result.mfe_kcal_mol)
            text_color = "#FFFFFF" if -result.mfe_kcal_mol > 4.0 else "#33414B"
            diagonal = (
                'box-shadow:inset 0 0 0 2px #C77D00;'
                if row_label == column_label
                else ""
            )
            cells.append(
                f'<td title="{html_module.escape(result.label)} · '
                f'MFE = {result.mfe_kcal_mol:.2f} kcal/mol" '
                f'style="padding:6px 8px;border:1px solid #E3EAE6;text-align:center;'
                f'font-size:12px;font-variant-numeric:tabular-nums;color:{text_color};'
                f'background:{background};{diagonal}">{result.mfe_kcal_mol:.1f}</td>'
            )
        rows.append(
            f'<tr>{"".join(cells)}</tr>'
        )

    legend_swatches = "".join(
        f'<span style="display:inline-block;width:22px;height:12px;'
        f'background:{_heatmap_color(-value)};margin:0 3px;vertical-align:middle;'
        f'border:1px solid #E3EAE6;"></span>{value:+.0f}'
        for value in (0, -2, -4, -6, -8)
    )
    return (
        '<div style="overflow-x:auto;">'
        f'<table style="border-collapse:collapse;min-width:60%;">'
        f'<thead><tr><th></th>{header_cells}</tr></thead>'
        f'<tbody>{"".join(rows)}</tbody></table></div>'
        f'<div style="font-size:12px;color:#5C6B73;margin-top:8px;">'
        f'MFE 色阶（kcal/mol）：{legend_swatches}　'
        "对角线（橙框）＝自二聚体；—＝未分析该组合。数值越负颜色越红，"
        "表示链间结合预测越强。</div>"
    )


# ---------------------------------------------------------------------------
# Standalone HTML report
# ---------------------------------------------------------------------------

_REPORT_CSS = """
  body {font-family: -apple-system, 'Segoe UI', 'Helvetica Neue', Arial, sans-serif;
        color:#172026; margin:0; background:#FBFDFC;}
  .wrap {max-width: 1000px; margin: 0 auto; padding: 32px 24px 48px;}
  h1 {font-size: 26px; letter-spacing:-.02em;}
  h2 {font-size: 19px; color:#087F5B; border-left:4px solid #087F5B;
      padding-left:10px; margin-top:40px;}
  .meta {color:#5C6B73; font-size:13px; margin: 4px 0 18px;}
  .chips span {display:inline-block; background:#E9F5EF; color:#086952;
      border:1px solid #BFE0D0; border-radius:999px; padding:2px 10px;
      font-size:12px; margin-right:6px;}
  table.data {border-collapse:collapse; width:100%; font-size:13px;}
  table.data th {background:#F4F9F6; text-align:left; padding:7px 9px;
      border:1px solid #E3EAE6;}
  table.data td {padding:7px 9px; border:1px solid #E3EAE6;
      font-variant-numeric:tabular-nums;}
  table.data td.mono {font-family:ui-monospace, Menlo, Consolas, monospace;}
  .card {border:1px solid #E3EAE6; border-radius:12px; padding:14px 16px;
      margin:14px 0; background:#FFFFFF; break-inside: avoid;}
  .card h3 {margin:0 0 6px; font-size:15px;}
  .metrics {color:#5C6B73; font-size:13px; margin-bottom:8px;}
  .interpret {background:#F4F9F6; border:1px solid #DCEAE2;
      border-left:4px solid #087F5B; border-radius:8px; padding:10px 14px;
      font-size:13px;}
  .interpret ul {margin:0 0 0 18px; padding:0;}
  .figure {text-align:center; margin:10px 0; break-inside: avoid;}
  .figure svg {max-width:100%; height:auto; border:1px solid #E3EAE6;
      border-radius:12px;}
  .note {font-size:12px; color:#5C6B73;}
  footer {margin-top:48px; padding-top:14px; border-top:1px solid #E3EAE6;
      font-size:12px; color:#5C6B73;}
  @media print { body {background:#FFFFFF;} .card, .figure {page-break-inside:avoid;} }
"""


def _report_result_table(results: list[FoldResult]) -> str:
    has_b = any(result.primer_b for result in results)
    headers = ["对象", "MFE (kcal/mol)", "总碱基对", "A 末5nt", "A 3′末端"]
    if has_b:
        headers += ["B 末5nt", "B 3′末端"]
    headers.append("点括号结构")

    head = "".join(f"<th>{html_module.escape(name)}</th>" for name in headers)
    rows = []
    for result in results:
        cells = [
            html_module.escape(result.label),
            f"{result.mfe_kcal_mol:.2f}",
            str(result.base_pairs),
            f"{result.a_last5_paired} / 5",
            "配对" if result.a_terminal_paired else "游离",
        ]
        if has_b:
            cells += [
                f"{result.b_last5_paired} / 5",
                "配对" if result.b_terminal_paired else "游离",
            ]
        cells.append(f'<td class="mono">{html_module.escape(result.structure)}</td>')
        row_cells = "".join(
            cell if cell.startswith("<td") else f"<td>{cell}</td>" for cell in cells
        )
        rows.append(f"<tr>{row_cells}</tr>")
    return (
        f'<table class="data"><thead><tr>{head}</tr></thead>'
        f'<tbody>{"".join(rows)}</tbody></table>'
    )


def build_html_report(
    *,
    primers: list[Primer],
    results: list[FoldResult],
    figures: dict[str, str],
    temperature_c: float,
    salt_m: float,
    na_m_for_tm: float | None = None,
    generated_at: datetime | None = None,
    tool_version: str = "1.0.0",
) -> str:
    """Assemble a standalone, self-contained HTML report."""

    generated_at = generated_at or datetime.now()
    na_for_tm = salt_m if na_m_for_tm is None else na_m_for_tm
    dimers_exist = any(result.kind != "hairpin" for result in results)

    parts: list[str] = [
        "<!DOCTYPE html><html lang='zh-CN'><head><meta charset='utf-8'>",
        "<meta name='viewport' content='width=device-width, initial-scale=1'>",
        "<title>PrimerFold 分析报告</title>",
        f"<style>{_REPORT_CSS}</style></head><body><div class='wrap'>",
        "<h1>🧬 PrimerFold 分析报告</h1>",
        f"<div class='meta'>生成时间：{generated_at:%Y-%m-%d %H:%M} · "
        f"PrimerFold v{html_module.escape(tool_version)} · ViennaRNA 本地分析</div>",
        "<div class='chips'>"
        f"<span>引物数 {len(primers)}</span>"
        f"<span>分析任务 {len(results)}</span>"
        f"<span>温度 {temperature_c:g} °C</span>"
        f"<span>单价盐 {salt_m * 1000:g} mM</span></div>",
        "<h2>引物属性</h2>",
        _report_qc_table(primers, na_for_tm),
        '<p class="note">Tm 采用 SantaLucia 1998 最近邻法，按单价盐浓度校正，'
        "总链浓度 250 nM（非自互补假设）；不含 Mg²⁺ 与 dNTP 校正。"
        "所有数值均为描述性参考，非合格判定。</p>",
    ]

    if dimers_exist:
        heatmap = heatmap_html(results)
        if heatmap:
            parts += ["<h2>互作热图</h2>", heatmap]

    grouped = {
        kind: [result for result in results if result.kind == kind]
        for kind in ("hairpin", "self_dimer", "cross_dimer")
    }
    for kind, group in grouped.items():
        if not group:
            continue
        parts += [
            f"<h2>{KIND_LABELS_LOCAL[kind]}（{len(group)}）</h2>",
            _report_result_table(group),
        ]
        for result in group:
            bullets = "".join(
                f"<li>{html_module.escape(bullet)}</li>"
                for bullet in interpret_result(result)
            )
            parts.append(
                f"<div class='card'><h3>{html_module.escape(result.label)}</h3>"
                f"<div class='metrics'>MFE {result.mfe_kcal_mol:.2f} kcal/mol · "
                f"碱基对 {result.base_pairs} · 分子间 {result.intermolecular_pairs} · "
                f"3′末5nt {result.a_last5_paired}"
                + (
                    f" / {result.b_last5_paired}"
                    if result.b_last5_paired is not None
                    else ""
                )
                + f"</div><div class='interpret'><ul>{bullets}</ul></div>"
            )
            figure = figures.get(result.label)
            if figure:
                figure_inner = figure[figure.find("<svg"):] if "<svg" in figure else figure
                parts.append(f"<div class='figure'>{figure_inner}</div>")
            parts.append(
                "<div class='note'>序列与结构逐位对照：</div>"
                + alignment_html(result.folded_sequence, result.structure)
                + alignment_legend_html()
                + "</div>"
            )

    parts += [
        "<footer>点括号结构：成对的 ( ) 表示碱基配对，. 表示未配对，& 表示两条链分隔。"
        "所有结果均为热力学模型预测（MFE 结构），不是实验测量；"
        "请结合 Mg²⁺、dNTP、引物浓度等实际 PCR 条件判断。"
        "由 PrimerFold 本地生成，序列未离开本机。</footer>",
        "</div></body></html>",
    ]
    return "".join(parts)


def _report_qc_table(primers: list[Primer], na_m: float) -> str:
    rows = qc_rows(primers, na_m=na_m)
    for primer, row in zip(primers, rows, strict=True):
        row_with_sequence = {"序列 (5′→3′)": primer.sequence, **row}
        row.clear()
        row.update(row_with_sequence)
    headers = list(rows[0].keys()) if rows else []
    head = "".join(
        f"<th>{html_module.escape(str(name))}</th>" for name in headers
    )
    body = "".join(
        "<tr>"
        + "".join(
            f"<td>{html_module.escape(str(row[name]))}</td>" for name in headers
        )
        + "</tr>"
        for row in rows
    )
    if not head:
        return "<p class='note'>无引物。</p>"
    return f'<table class="data"><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>'
