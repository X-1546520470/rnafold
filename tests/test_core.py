from __future__ import annotations

import unittest

from primerfold.core import (
    Primer,
    PrimerInputError,
    build_jobs,
    fold_dimer,
    fold_hairpin,
    parse_primers,
    render_structure_svg,
)


class PrimerParsingTests(unittest.TestCase):
    def test_parse_fasta(self) -> None:
        primers = parse_primers(
            ">Forward primer\nATGACCATGATTACGCCAAG\n>Reverse\nGCGCGCTTTTTGCGCGC\n"
        )
        self.assertEqual([primer.name for primer in primers], ["Forward", "Reverse"])
        self.assertEqual(primers[0].sequence, "ATGACCATGATTACGCCAAG")

    def test_parse_csv_with_header(self) -> None:
        primers = parse_primers(
            "name,sequence\nForward,ATGACCATGATTACGCCAAG\nReverse,GCGCGCTTTTTGCGCGC"
        )
        self.assertEqual(len(primers), 2)
        self.assertEqual(primers[1].name, "Reverse")

    def test_parse_sequence_only(self) -> None:
        primers = parse_primers("ATGACCATGATTACGCCAAG\nGCGCGCTTTTTGCGCGC")
        self.assertEqual([primer.name for primer in primers], ["Primer_1", "Primer_2"])

    def test_reject_invalid_or_duplicate(self) -> None:
        with self.assertRaises(PrimerInputError):
            parse_primers(">P1\nACGTNACGT")
        with self.assertRaises(PrimerInputError):
            parse_primers(">P1\nACGTACGT\n>P1\nTGCATGCA")

    def test_build_all_jobs(self) -> None:
        primers = [Primer("P1", "ACGTACGT"), Primer("P2", "TGCATGCA")]
        jobs = build_jobs(
            primers,
            hairpin=True,
            self_dimer=True,
            cross_dimer=True,
        )
        self.assertEqual(len(jobs), 5)


class ViennaIntegrationTests(unittest.TestCase):
    def test_hairpin_and_svg(self) -> None:
        result = fold_hairpin(
            Primer("hairpin", "GCGCGCTTTTTGCGCGC"),
            temperature_c=37.0,
            salt_m=1.021,
        )
        self.assertEqual(result.structure, "((((((.....))))))")
        self.assertLess(result.mfe_kcal_mol, -8.0)
        self.assertEqual(result.base_pairs, 6)
        self.assertEqual(result.a_last5_paired, 5)
        self.assertTrue(result.a_terminal_paired)

        svg = render_structure_svg(result.folded_sequence, result.structure)
        self.assertIn("<svg", svg)
        self.assertIn("class=\"basepairs\"", svg)
        self.assertNotIn("<script", svg)

    def test_cross_dimer_metrics(self) -> None:
        result = fold_dimer(
            Primer("A", "GCGCGCTTTTTGCGCGC"),
            Primer("B", "GCGCGCTTTTTGCGCGC"),
            kind="cross_dimer",
            temperature_c=37.0,
            salt_m=1.021,
        )
        self.assertIn("&", result.structure)
        self.assertGreater(result.intermolecular_pairs, 0)
        self.assertEqual(result.base_pairs, result.intermolecular_pairs)
        self.assertIsNotNone(result.b_last5_paired)

    def test_no_duplex_sentinel_becomes_unbound_state(self) -> None:
        result = fold_dimer(
            Primer("A", "AAAAAAAAAA"),
            Primer("B", "AAAAAAAAAA"),
            kind="cross_dimer",
            temperature_c=37.0,
            salt_m=0.05,
        )
        self.assertEqual(result.mfe_kcal_mol, 0.0)
        self.assertEqual(result.intermolecular_pairs, 0)
        self.assertEqual(result.structure, "..........&..........")


if __name__ == "__main__":
    unittest.main()
