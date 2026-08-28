"""Styled secondary-structure figures rebuilt from RNAplot geometry.

RNAplot handles the difficult part (the 2D coordinate layout). This module
parses its SVG output and re-renders the same geometry with the GUI's visual
language: base-coloured discs, pair rungs coloured by hydrogen-bond count,
highlighted 3' ends and an inline legend. Everything is drawn with inline
presentation attributes, so the SVG carries no global CSS when embedded or
downloaded.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
import re
from typing import Sequence
from xml.sax.saxutils import escape

from .core import (
    FoldExecutionError,
    base_pair_coordinates,
    render_structure_svg,
)

BASE_COLORS = {
    "A": "#2F9E44",
    "C": "#1971C2",
    "G": "#E8890C",
    "T": "#E03131",
    "U": "#E03131",
}
GC_RUNG_COLOR = "#1098AD"
AT_RUNG_COLOR = "#AAB4BE"
BACKBONE_COLOR = "#9AA9B7"
HALO_COLOR = "#7048E8"
LABEL_COLOR = "#5C6B73"
LEGEND_TEXT_COLOR = "#4A5568"
CARD_BORDER_COLOR = "#DCE7E1"
FONT_STACK = "-apple-system, 'Segoe UI', 'Helvetica Neue', Arial, sans-serif"

DISC_RADIUS = 5.6
HALO_RADIUS = 7.4
LETTER_FONT_SIZE = 6.2
BACKBONE_WIDTH = 2.4
CANVAS_PADDING = 17.0

_NUCLEOTIDE_RE = re.compile(
    r'<text class="nucleotide"\s+x="([-\d.]+)"\s+y="([-\d.]+)">([^<]+)</text>'
)
_BACKBONE_RE = re.compile(
    r'<polyline class="backbone"[^>]*?points="(.*?)"', re.DOTALL
)
_PAIR_RE = re.compile(r'<line class="basepairs"\s+id="(\d+),(\d+)"')
_DIMENSIONS_RE = re.compile(r'width="(\d+)"\s+height="(\d+)"')


@dataclass(frozen=True, slots=True)
class StructureGeometry:
    """Nucleotide coordinates and pairs extracted from an RNAplot SVG."""

    letters: tuple[str, ...]
    xs: tuple[float, ...]
    ys: tuple[float, ...]
    pairs: tuple[tuple[int, int], ...]
    strand_lengths: tuple[int, ...]

    @property
    def strand_starts(self) -> tuple[int, ...]:
        starts: list[int] = []
        offset = 0
        for length in self.strand_lengths:
            starts.append(offset)
            offset += length + 1
        return tuple(starts)


def _strand_lengths(sequence: str) -> tuple[int, ...]:
    return tuple(len(part) for part in sequence.split("&"))


def _expected_pairs(sequence: str, structure: str) -> set[tuple[int, int]]:
    starts = []
    offset = 0
    for length in _strand_lengths(sequence):
        starts.append(offset)
        offset += length + 1

    expected: set[tuple[int, int]] = set()
    for (strand_a, pos_a), (strand_b, pos_b) in base_pair_coordinates(structure):
        expected.add((starts[strand_a] + pos_a, starts[strand_b] + pos_b))
    return expected


def parse_geometry(svg: str, sequence: str, structure: str) -> StructureGeometry:
    """Extract nucleotide coordinates and base pairs from an RNAplot SVG."""

    texts = _NUCLEOTIDE_RE.findall(svg)
    if len(texts) != len(sequence):
        raise FoldExecutionError(
            f"RNAplot 输出包含 {len(texts)} 个碱基标注，与序列长度 {len(sequence)} 不符。"
        )
    letters: list[str] = []
    xs: list[float] = []
    ys: list[float] = []
    for index, (x_text, y_text, letter) in enumerate(texts):
        if letter.upper() != sequence[index].upper():
            raise FoldExecutionError(
                f"RNAplot 碱基标注与序列不一致：位置 {index + 1} 应为 "
                f"{sequence[index]}，实际为 {letter}。"
            )
        letters.append(letter)
        xs.append(float(x_text))
        ys.append(float(y_text))

    backbone = _BACKBONE_RE.search(svg)
    if backbone is None:
        raise FoldExecutionError("RNAplot 输出中缺少骨架折线。")
    points = [
        (float(piece.split(",")[0]), float(piece.split(",")[1]))
        for piece in backbone.group(1).split()
        if piece
    ]
    if len(points) != len(sequence):
        raise FoldExecutionError(
            f"RNAplot 骨架包含 {len(points)} 个坐标点，与序列长度 {len(sequence)} 不符。"
        )

    pairs: set[tuple[int, int]] = set()
    for first, second in _PAIR_RE.findall(svg):
        pair = (int(first) - 1, int(second) - 1)
        if not (0 <= pair[0] < len(sequence) and 0 <= pair[1] < len(sequence)):
            raise FoldExecutionError("RNAplot 碱基对编号越界。")
        pairs.add(pair)
    if pairs != _expected_pairs(sequence, structure):
        raise FoldExecutionError("RNAplot 碱基对与点括号结构不一致。")

    return StructureGeometry(
        letters=tuple(letters),
        xs=tuple(xs),
        ys=tuple(ys),
        pairs=tuple(sorted(pairs)),
        strand_lengths=_strand_lengths(sequence),
    )


def svg_pixel_dimensions(svg: str) -> tuple[int, int]:
    match = _DIMENSIONS_RE.search(svg)
    if match is None:
        raise FoldExecutionError("SVG 缺少画布尺寸。")
    return int(match.group(1)), int(match.group(2))


def _text_width(text: str, font_size: float) -> float:
    width = 0.0
    for character in text:
        width += font_size if ord(character) > 0x2E80 else font_size * 0.58
    return width


def render_styled_structure_svg(
    sequence: str,
    structure: str,
    *,
    strand_labels: Sequence[str] | None = None,
    highlight_last_n: int = 5,
    target_size_px: int = 780,
) -> str:
    """Render sequence/structure as a styled SVG using RNAplot's layout."""

    raw = render_structure_svg(sequence, structure)
    geometry = parse_geometry(raw, sequence, structure)
    return build_styled_svg(
        geometry,
        strand_labels=strand_labels,
        highlight_last_n=highlight_last_n,
        target_size_px=target_size_px,
    )


def build_styled_svg(
    geometry: StructureGeometry,
    *,
    strand_labels: Sequence[str] | None = None,
    highlight_last_n: int = 5,
    target_size_px: int = 780,
) -> str:
    base_positions = [
        index
        for index, letter in enumerate(geometry.letters)
        if letter != "&"
    ]
    if not base_positions:
        raise FoldExecutionError("结构中没有可绘制的碱基。")

    min_x = min(geometry.xs[index] for index in base_positions)
    max_x = max(geometry.xs[index] for index in base_positions)
    min_y = min(geometry.ys[index] for index in base_positions)
    max_y = max(geometry.ys[index] for index in base_positions)

    def px(index: int) -> tuple[float, float]:
        return geometry.xs[index] - min_x + CANVAS_PADDING, geometry.ys[index] - min_y + CANVAS_PADDING

    main_width = max_x - min_x + 2 * CANVAS_PADDING
    main_height = max_y - min_y + 2 * CANVAS_PADDING

    strand_members: list[list[int]] = []
    cursor = 0
    for length in geometry.strand_lengths:
        strand_members.append(list(range(cursor, cursor + length)))
        cursor += length + 1

    pair_colors: dict[int, str] = {}
    for first, second in geometry.pairs:
        letter_pair = {geometry.letters[first], geometry.letters[second]}
        if letter_pair <= {"G", "C"}:
            color = GC_RUNG_COLOR
        elif letter_pair <= {"A", "T", "U"}:
            color = AT_RUNG_COLOR
        else:
            color = BACKBONE_COLOR
        pair_colors[first] = color
        pair_colors[second] = color

    parts: list[str] = []

    legend_items = _legend_items(geometry, strand_labels)
    legend_fragment, legend_height = _render_legend(legend_items, main_width)
    total_width = main_width
    total_height = main_height + legend_height
    scale = max(2.0, min(4.4, target_size_px / max(total_width, total_height)))

    parts.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{total_width * scale:.0f}" '
        f'height="{total_height * scale:.0f}" viewBox="0 0 {total_width:.0f} {total_height:.0f}" '
        'role="img" aria-label="DNA 二级结构图">'
    )
    parts.append(
        f'<rect x="0.5" y="0.5" width="{total_width - 1:.1f}" height="{total_height - 1:.1f}" '
        f'rx="14" fill="#FFFFFF" stroke="{CARD_BORDER_COLOR}"/>'
    )

    # 3' end halos sit underneath everything else as soft highlights.
    if highlight_last_n > 0:
        for members in strand_members:
            if not members:
                continue
            for index in members[-highlight_last_n:]:
                x, y = px(index)
                parts.append(
                    f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{HALO_RADIUS}" '
                    f'fill="{HALO_COLOR}" fill-opacity="0.08" '
                    f'stroke="{HALO_COLOR}" stroke-width="1.4" stroke-dasharray="2.4 1.7"/>'
                )

    # Pair rungs, then backbone, then discs, then letters.
    for first, second in geometry.pairs:
        x1, y1 = px(first)
        x2, y2 = px(second)
        color = pair_colors[first]
        width = 2.2 if color == GC_RUNG_COLOR else 1.7
        parts.append(
            f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" '
            f'stroke="{color}" stroke-width="{width}" stroke-linecap="round"/>'
        )

    for members in strand_members:
        if len(members) < 2:
            continue
        points = " ".join(f"{px(index)[0]:.2f},{px(index)[1]:.2f}" for index in members)
        parts.append(
            f'<polyline points="{points}" fill="none" stroke="{BACKBONE_COLOR}" '
            f'stroke-width="{BACKBONE_WIDTH}" stroke-linecap="round" stroke-linejoin="round"/>'
        )

    paired = {index for pair in geometry.pairs for index in pair}
    for index in base_positions:
        letter = geometry.letters[index]
        x, y = px(index)
        color = BASE_COLORS.get(letter, BACKBONE_COLOR)
        if index in paired:
            parts.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{DISC_RADIUS}" fill="{color}"/>')
            text_fill = "#FFFFFF"
        else:
            parts.append(
                f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{DISC_RADIUS}" '
                f'fill="#FFFFFF" stroke="{color}" stroke-width="1.4"/>'
            )
            text_fill = color
        parts.append(
            f'<text class="pf-letter" x="{x:.2f}" y="{y + 2.1:.2f}" text-anchor="middle" '
            f'font-family="{FONT_STACK}" font-size="{LETTER_FONT_SIZE}" '
            f'font-weight="600" fill="{text_fill}">{escape(letter)}</text>'
        )

    # 5'/3' end labels, pushed outward from the structure centroid.
    centroid_x = sum(geometry.xs[index] for index in base_positions) / len(base_positions)
    centroid_y = sum(geometry.ys[index] for index in base_positions) / len(base_positions)
    for members in strand_members:
        if not members:
            continue
        for index, label in ((members[0], "5′"), (members[-1], "3′")):
            dx = geometry.xs[index] - centroid_x
            dy = geometry.ys[index] - centroid_y
            norm = math.hypot(dx, dy)
            if norm < 1e-6:
                dx, dy = 0.0, -1.0
            else:
                dx, dy = dx / norm, dy / norm
            offset = HALO_RADIUS + 4.6
            x = geometry.xs[index] - min_x + CANVAS_PADDING + dx * offset
            y = geometry.ys[index] - min_y + CANVAS_PADDING + dy * offset
            parts.append(
                f'<text x="{x:.2f}" y="{y + 1.9:.2f}" text-anchor="middle" '
                f'font-family="{FONT_STACK}" font-size="5.4" font-weight="700" '
                f'fill="{LABEL_COLOR}">{label}</text>'
            )

    parts.append(f'<g transform="translate(0, {main_height:.1f})">{legend_fragment}</g>')
    parts.append("</svg>")
    return "".join(parts)


def _legend_items(
    geometry: StructureGeometry,
    strand_labels: Sequence[str] | None,
) -> list[tuple[str, str]]:
    items: list[tuple[str, str]] = []
    if len(geometry.strand_lengths) > 1:
        labels = list(strand_labels or [])
        names = []
        for index, length in enumerate(geometry.strand_lengths):
            name = labels[index] if index < len(labels) else f"链 {index + 1}"
            names.append(f"链 {index + 1}：{name}（{length} nt）")
        items.append(("", "　".join(names)))
    items.append(
        (
            f'<circle cx="4" cy="0" r="3.2" fill="{LABEL_COLOR}"/>',
            "实心圆＝参与配对的碱基",
        )
    )
    items.append(
        (
            f'<circle cx="4" cy="0" r="3.2" fill="#FFFFFF" stroke="{LABEL_COLOR}" stroke-width="1.1"/>',
            "空心圆＝未配对碱基",
        )
    )
    items.append(
        (
            f'<line x1="0" y1="0" x2="8" y2="0" stroke="{GC_RUNG_COLOR}" '
            'stroke-width="2" stroke-linecap="round"/>',
            "G≡C（3 个氢键）",
        )
    )
    items.append(
        (
            f'<line x1="0" y1="0" x2="8" y2="0" stroke="{AT_RUNG_COLOR}" '
            'stroke-width="2" stroke-linecap="round"/>',
            "A＝T（2 个氢键）",
        )
    )
    items.append(
        (
            f'<circle cx="4" cy="0" r="3.2" fill="none" stroke="{HALO_COLOR}" '
            'stroke-width="1.3" stroke-dasharray="1.8 1.2"/>',
            "紫色虚线圈＝3′端最后 5 nt",
        )
    )
    return items


def _render_legend(
    items: list[tuple[str, str]],
    available_width: float,
) -> tuple[str, float]:
    font_size = 5.8
    row_height = 9.6
    icon_width = 9.0
    gap = 13.0
    margin = CANVAS_PADDING

    fragments: list[str] = []
    cursor_x = margin
    row = 0
    for icon, text in items:
        item_width = (
            (icon_width + 2.6 if icon else 0.0)
            + _text_width(text, font_size)
        )
        if cursor_x > margin and cursor_x + item_width > available_width - margin:
            row += 1
            cursor_x = margin
        y = row * row_height + row_height / 2
        if icon:
            fragments.append(
                f'<g transform="translate({cursor_x:.1f}, {y:.1f})">{icon}</g>'
            )
            cursor_x += icon_width + 2.6
        fragments.append(
            f'<text x="{cursor_x:.1f}" y="{y + 2.0:.1f}" font-family="{FONT_STACK}" '
            f'font-size="{font_size}" fill="{LEGEND_TEXT_COLOR}">{escape(text)}</text>'
        )
        cursor_x += item_width + gap
    return "".join(fragments), (row + 1) * row_height + 6.0
