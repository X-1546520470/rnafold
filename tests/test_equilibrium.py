from __future__ import annotations

import math
import unittest

from primerfold.core import CofoldEnergies, Primer, cofold_energies, parse_primers
from primerfold.equilibrium import (
    equilibrium_constants,
    interpret_equilibrium,
    species_equilibrium,
)

PRIMER_R = "GCGCGCTTTTTGCGCGC"
PRIMER_F = "ATGACCATGATTACGCCAAG"


class EquilibriumSolverTests(unittest.TestCase):
    def test_equilibrium_constants_from_association_energies(self) -> None:
        energies = CofoldEnergies(ab=-5.32, aa=-2.29, bb=-13.31, a=-0.086, b=-6.74)
        k_ab, k_aa, k_bb = equilibrium_constants(energies, temperature_c=37.0)
        rt = 1.9872041e-3 * (37.0 + 273.15)
        self.assertAlmostEqual(k_ab, math.exp(5.32 / rt), places=8)
        self.assertAlmostEqual(k_aa, math.exp(2.29 / rt), places=8)
        self.assertAlmostEqual(k_bb, math.exp(13.31 / rt), places=8)

    def test_weak_association_stays_mostly_free(self) -> None:
        # Association energy ~0 -> K ~ 1/M -> almost nothing dimerised.
        energies = CofoldEnergies(ab=-0.05, aa=-0.05, bb=-0.05, a=-0.01, b=-0.01)
        equilibrium = species_equilibrium(
            energies, temperature_c=37.0, conc_a=5e-7, conc_b=5e-7
        )
        self.assertGreater(equilibrium.free_a / equilibrium.conc_a, 0.99)
        self.assertLess(equilibrium.ab / equilibrium.conc_a, 0.001)

    def test_mass_balance_and_equilibrium_conditions(self) -> None:
        energies = CofoldEnergies(ab=-12.0, aa=-6.0, bb=-9.0, a=-1.0, b=-2.0)
        conc_a, conc_b = 8e-7, 5e-7
        eq = species_equilibrium(
            energies, temperature_c=37.0, conc_a=conc_a, conc_b=conc_b
        )
        # Mass balance for each strand (AB consumes one of each, AA/BB two).
        self.assertAlmostEqual(eq.free_a + eq.ab + 2 * eq.aa, conc_a, places=15)
        self.assertAlmostEqual(eq.free_b + eq.ab + 2 * eq.bb, conc_b, places=15)
        # Equilibrium conditions AB = K_AB [A][B] etc.
        k_ab, k_aa, k_bb = equilibrium_constants(energies, temperature_c=37.0)
        self.assertAlmostEqual(eq.ab, k_ab * eq.free_a * eq.free_b, places=18)
        self.assertAlmostEqual(eq.aa, k_aa * eq.free_a**2, places=18)
        self.assertAlmostEqual(eq.bb, k_bb * eq.free_b**2, places=18)

    def test_dimerisation_increases_with_concentration(self) -> None:
        energies = CofoldEnergies(ab=-10.0, aa=-5.0, bb=-5.0, a=-0.5, b=-0.5)
        fractions = []
        for conc in (1e-7, 1e-6, 1e-5):
            eq = species_equilibrium(
                energies, temperature_c=37.0, conc_a=conc, conc_b=conc
            )
            fractions.append((eq.ab + eq.aa + eq.bb) / conc)
        self.assertLess(fractions[0], fractions[1])
        self.assertLess(fractions[1], fractions[2])

    def test_occupied_fraction_and_interpretation(self) -> None:
        energies = CofoldEnergies(ab=-13.31, aa=-13.31, bb=-13.31, a=-6.74, b=-6.74)
        eq = species_equilibrium(
            energies, temperature_c=37.0, conc_a=5e-7, conc_b=5e-7
        )
        self.assertAlmostEqual(eq.a_occupied_fraction, (eq.ab + 2 * eq.aa) / 5e-7)
        bullets = interpret_equilibrium(eq, conc_nm=500)
        joined = "\n".join(bullets)
        self.assertIn("500 nM", joined)
        self.assertIn("AB", joined)
        self.assertIn("不含 Mg²⁺", joined)


class CofoldEnergiesTests(unittest.TestCase):
    def test_parse_real_rnacofold_energies(self) -> None:
        primers = parse_primers(f">F\n{PRIMER_F}\n>R\n{PRIMER_R}\n")
        energies = cofold_energies(
            primers[0], primers[1], temperature_c=37.0, salt_m=0.05
        )
        # The printed AB column is the association free energy; an ensemble
        # association energy can never be below the MFE of the duplex (-6.3).
        self.assertLess(energies.ab, 0.0)
        self.assertGreater(energies.ab, -6.4)
        self.assertLess(energies.bb, -5.0)

    def test_identical_sequences_have_equal_ab_and_aa(self) -> None:
        primer = Primer("R", PRIMER_R)
        energies = cofold_energies(primer, primer, temperature_c=37.0, salt_m=0.05)
        self.assertAlmostEqual(energies.ab, energies.aa, places=6)
        self.assertAlmostEqual(energies.aa, energies.bb, places=6)

    def test_strong_self_dimer_highly_dimerised_at_5uM(self) -> None:
        primer = Primer("R", PRIMER_R)
        energies = cofold_energies(primer, primer, temperature_c=37.0, salt_m=0.05)
        eq = species_equilibrium(
            energies, temperature_c=37.0, conc_a=5e-6, conc_b=5e-6
        )
        self.assertGreater(eq.a_occupied_fraction, 0.9)


if __name__ == "__main__":
    unittest.main()
