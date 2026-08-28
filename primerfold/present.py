"""Pure-presentation helpers shared by the GUI: interpretation text and a
colour-coded sequence/structure alignment rendered as inline-styled HTML.

Kept free of Streamlit so the logic stays unit-testable.
"""

from __future__ import annotations

from xml.sax.saxutils import escape

from .core import FoldResult, base_pair_coordinates
from .plot import BASE_COLORS

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
