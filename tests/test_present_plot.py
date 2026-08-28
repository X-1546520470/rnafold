from __future__ import annotations

import unittest

from primerfold.core import FoldResult, Primer
from primerfold.plot import (
    parse_geometry,
    render_structure_svg,
    render_styled_structure_svg,
    svg_pixel_dimensions,
)
from primerfold.present import (
    alignment_html,
    classified_pairs,
    interpret_result,
    terminal_positions,
)

HAIRPIN_SEQUENCE = "GCGCGCTTTTTGCGCGC"
HAIRPIN_STRUCTURE = "((((((.....))))))"


def _hairpin_result() -> FoldResult:
    from primerfold.core import fold_hairpin

    return fold_hairpin(
        Primer("Primer_F", HAIRPIN_SEQUENCE),
        temperature_c=37.0,
        salt_m=0.05,
    )


class StyledSvgTests(unittest.TestCase):
    def test_styled_hairpin_contains_expected_features(self) -> None:
        svg = render_styled_structure_svg(
            HAIRPIN_SEQUENCE, HAIRPIN_STRUCTURE, strand_labels=["Primer_F"]
        )
        self.assertIn("<svg", svg)
        self.assertIn('class="pf-letter"', svg)
        self.assertIn("5′", svg)
        self.assertIn("3′", svg)
        self.assertIn("实心圆", svg)
        self.assertIn("3′端最后 5 nt", svg)
        # Paired G discs use the amber base colour; loop T discs stay hollow.
        self.assertIn("#E8890C", svg)
        # No top-level <style> block: all styling must stay inline.
        self.assertNotIn("<style", svg)
        self.assertNotIn("<script", svg)

    def test_styled_duplex_lists_strands_and_skips_ampersand(self) -> None:
        sequence = f"{HAIRPIN_SEQUENCE}&{HAIRPIN_SEQUENCE}"
        structure = "((((((.....((((((&)))))).....))))))"
        svg = render_styled_structure_svg(
            sequence, structure, strand_labels=["Primer_F", "Primer_R"]
        )
        self.assertIn("Primer_F", svg)
        self.assertIn("Primer_R", svg)
        # One letter disc per nucleotide; the & separator gets no disc.
        self.assertEqual(svg.count('class="pf-letter"'), len(sequence) - 1)
        width, height = svg_pixel_dimensions(svg)
        self.assertGreater(width, 100)
        self.assertGreater(height, 100)

    def test_parse_geometry_matches_structure_pairs(self) -> None:
        raw = render_structure_svg(HAIRPIN_SEQUENCE, HAIRPIN_STRUCTURE)
        geometry = parse_geometry(raw, HAIRPIN_SEQUENCE, HAIRPIN_STRUCTURE)
        self.assertEqual(len(geometry.letters), len(HAIRPIN_SEQUENCE))
        self.assertEqual(len(geometry.pairs), 6)

    def test_parse_geometry_rejects_mismatched_sequence(self) -> None:
        raw = render_structure_svg(HAIRPIN_SEQUENCE, HAIRPIN_STRUCTURE)
        from primerfold.core import FoldExecutionError

        with self.assertRaises(FoldExecutionError):
            parse_geometry(raw, "A" * len(HAIRPIN_SEQUENCE), HAIRPIN_STRUCTURE)


class PresentTests(unittest.TestCase):
    def test_terminal_positions_covers_each_strand_tail(self) -> None:
        sequence = "AAAACCCCC&GGGGTTTTT"
        positions = terminal_positions(sequence, last_n=5)
        # Each strand is 9 nt: its final 5 positions are 4-8 and 14-18 globally.
        self.assertEqual(positions, {4, 5, 6, 7, 8, 14, 15, 16, 17, 18})

    def test_classified_pairs_colours_inter_strand(self) -> None:
        sequence = "GCGC&CGCG"
        structure = "((((&))))"
        pairs = classified_pairs(sequence, structure)
        self.assertEqual(len(pairs), 4)
        for _, _, colour in pairs:
            self.assertEqual(colour, "#D6336C")

    def test_alignment_html_renders_rows_and_legend_entities(self) -> None:
        html = alignment_html(
            f"{HAIRPIN_SEQUENCE}&{HAIRPIN_SEQUENCE}",
            "((((((.....((((((&)))))).....))))))",
        )
        self.assertIn("&amp;", html)
        self.assertNotIn("pf-letter", html)  # alignment uses plain spans
        self.assertIn('border-bottom:2px solid', html)
        # Position ruler starts at 1.
        self.assertIn(">1<", html)

    def test_interpret_result_hairpin_mentions_pairs_and_mfe(self) -> None:
        result = _hairpin_result()
        bullets = interpret_result(result)
        joined = "\n".join(bullets)
        self.assertIn("发卡", joined)
        self.assertIn(f"{result.base_pairs} 个链内碱基对", joined)
        self.assertIn("3′末端", joined)
        self.assertIn("非合格判定", joined)

    def test_interpret_result_unbound_dimer(self) -> None:
        from primerfold.core import fold_dimer

        result = fold_dimer(
            Primer("A", "AAAAAAAAAA"),
            Primer("B", "AAAAAAAAAA"),
            kind="cross_dimer",
            temperature_c=37.0,
            salt_m=0.05,
        )
        bullets = interpret_result(result)
        joined = "\n".join(bullets)
        self.assertIn("未检出链间碱基对", joined)
        self.assertIn("保持游离", joined)


if __name__ == "__main__":
    unittest.main()
