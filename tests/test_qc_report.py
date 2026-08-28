from __future__ import annotations

import unittest

from primerfold.core import Primer, parse_primers, fold_dimer, fold_hairpin
from primerfold.present import build_html_report, heatmap_html
from primerfold.qc import (
    analyze_primer_qc,
    gc_clamp_3p,
    gc_percent,
    max_homopolymer,
    molecular_weight,
    qc_rows,
    tm_nearest_neighbor,
)

PRIMER_F = "ATGACCATGATTACGCCAAG"
PRIMER_R = "GCGCGCTTTTTGCGCGC"


class QcMetricTests(unittest.TestCase):
    def test_gc_percent(self) -> None:
        self.assertAlmostEqual(gc_percent(PRIMER_F), 45.0)
        self.assertAlmostEqual(gc_percent("GGCC"), 100.0)
        self.assertAlmostEqual(gc_percent("ATAT"), 0.0)

    def test_molecular_weight_known_value(self) -> None:
        # A=313.21, C=289.18, G=329.21, T=304.2; minus 61.96 offset.
        self.assertAlmostEqual(molecular_weight("ACGT"), 1173.84, places=2)

    def test_tm_in_expected_range_and_salt_dependent(self) -> None:
        tm_low_salt = tm_nearest_neighbor(PRIMER_F, na_m=0.005)
        tm_default = tm_nearest_neighbor(PRIMER_F, na_m=0.05)
        tm_high_salt = tm_nearest_neighbor(PRIMER_F, na_m=1.0)
        # Monovalent-only NN Tm for a 45% GC 20-mer sits near 50-60 C.
        self.assertTrue(45.0 < tm_default < 62.0, tm_default)
        self.assertLess(tm_low_salt, tm_default)
        self.assertLess(tm_default, tm_high_salt)

    def test_gc_clamp_and_homopolymer(self) -> None:
        self.assertEqual(gc_clamp_3p("ACGTG"), 3)
        self.assertEqual(gc_clamp_3p("ACGTA"), 2)
        self.assertEqual(gc_clamp_3p("AAAAA"), 0)
        self.assertEqual(max_homopolymer("AAAGGGGT"), (4, "G"))
        self.assertEqual(max_homopolymer("ACACAC"), (1, "A"))

    def test_analyze_primer_qc_row_shape(self) -> None:
        rows = qc_rows([Primer("F", PRIMER_F)], na_m=0.05)
        self.assertEqual(len(rows), 1)
        row = rows[0]
        self.assertEqual(row["引物"], "F")
        self.assertEqual(row["长度 (nt)"], 20)
        self.assertEqual(row["GC (%)"], 45.0)
        self.assertEqual(row["3′GC clamp"], 3)
        qc = analyze_primer_qc(Primer("F", PRIMER_F), na_m=0.05)
        self.assertAlmostEqual(qc.tm_c, tm_nearest_neighbor(PRIMER_F, na_m=0.05))


class HeatmapTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        primers = parse_primers(f">F\n{PRIMER_F}\n>R\n{PRIMER_R}\n")
        cls.results = [
            fold_hairpin(primers[0], temperature_c=37.0, salt_m=0.05),
            fold_hairpin(primers[1], temperature_c=37.0, salt_m=0.05),
            fold_dimer(
                primers[0], primers[0], kind="self_dimer",
                temperature_c=37.0, salt_m=0.05,
            ),
            fold_dimer(
                primers[0], primers[1], kind="cross_dimer",
                temperature_c=37.0, salt_m=0.05,
            ),
        ]

    def test_heatmap_contains_matrix_cells_and_legend(self) -> None:
        html = heatmap_html(self.results)
        self.assertIn("<table", html)
        self.assertIn(">F</th>", html)
        self.assertIn(">R</th>", html)
        self.assertIn("对角线", html)
        self.assertIn("MFE 色阶", html)

    def test_heatmap_empty_without_dimers(self) -> None:
        primers = parse_primers(f">F\n{PRIMER_F}\n")
        hairpin_only = [fold_hairpin(primers[0], temperature_c=37.0, salt_m=0.05)]
        self.assertEqual(heatmap_html(hairpin_only), "")


class HtmlReportTests(unittest.TestCase):
    def test_report_is_standalone_and_complete(self) -> None:
        primers = parse_primers(f">F\n{PRIMER_F}\n>R\n{PRIMER_R}\n")
        hairpin = fold_hairpin(primers[0], temperature_c=37.0, salt_m=0.05)
        cross = fold_dimer(
            primers[0], primers[1], kind="cross_dimer",
            temperature_c=37.0, salt_m=0.05,
        )
        figures = {
            hairpin.label: "<svg xmlns='x'><g>fig</g></svg>",
        }
        report = build_html_report(
            primers=primers,
            results=[hairpin, cross],
            figures=figures,
            temperature_c=37.0,
            salt_m=0.05,
            generated_at=__import__("datetime").datetime(2026, 1, 2, 3, 4),
            tool_version="1.0.0",
        ).encode("utf-8").decode("utf-8")
        self.assertTrue(report.startswith("<!DOCTYPE html>"))
        self.assertIn("PrimerFold 分析报告", report)
        self.assertIn("引物属性", report)
        self.assertIn("互作热图", report)
        self.assertIn("ATGACCATGATTACGCCAAG", report)
        self.assertIn("链间结合", report)
        self.assertIn("<svg", report)
        self.assertIn("序列与结构逐位对照", report)
        self.assertIn("2026-01-02 03:04", report)


if __name__ == "__main__":
    unittest.main()
