"""Two-strand equilibrium concentrations from ensemble free energies.

Five-species statistical-mechanics calculation, equivalent to RNAcofold's
built-in concentration routine: equilibrium constants are derived from the
ensemble free energies (AB, AA, BB, A, B) and the mass balance is solved by
monotone bisection.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

from .core import CofoldEnergies

GAS_CONSTANT_KCAL = 1.9872041e-3  # kcal/(mol·K)


@dataclass(frozen=True, slots=True)
class EquilibriumResult:
    conc_a: float  # initial strand concentration (M)
    conc_b: float
    free_a: float  # equilibrium concentration of unpaired A (M)
    free_b: float
    ab: float
    aa: float
    bb: float

    @property
    def a_occupied_fraction(self) -> float:
        """Fraction of A molecules inside any dimer (AB or AA)."""

        return (self.ab + 2.0 * self.aa) / self.conc_a if self.conc_a else 0.0

    @property
    def b_occupied_fraction(self) -> float:
        return (self.ab + 2.0 * self.bb) / self.conc_b if self.conc_b else 0.0


def equilibrium_constants(
    energies: CofoldEnergies,
    *,
    temperature_c: float,
) -> tuple[float, float, float]:
    """Association constants (1/M) for AB, AA and BB from ensemble energies.

    RNAcofold's ``Free Energies`` block reports the association free energies
    ``F_AB − F_A − F_B`` and ``F_AA − 2·F_A`` / ``F_BB − 2·F_B`` in the first
    three columns (the last two columns are the absolute monomer ensemble
    energies), so ``K = exp(−ΔG°/RT)`` applies to them directly.
    """

    rt = GAS_CONSTANT_KCAL * (temperature_c + 273.15)
    k_ab = math.exp(-energies.ab / rt)
    k_aa = math.exp(-energies.aa / rt)
    k_bb = math.exp(-energies.bb / rt)
    return k_ab, k_aa, k_bb


def species_equilibrium(
    energies: CofoldEnergies,
    *,
    temperature_c: float,
    conc_a: float,
    conc_b: float,
) -> EquilibriumResult:
    """Solve [A], [B], [AB], [AA], [BB] for the given initial concentrations."""

    if conc_a <= 0 or conc_b <= 0:
        raise ValueError("初始浓度必须为正数。")
    k_ab, k_aa, k_bb = equilibrium_constants(energies, temperature_c=temperature_c)

    def free_a_given_b(free_b: float) -> float:
        # 2·k_aa·x² + (1 + k_ab·y)·x − conc_a = 0
        factor = 1.0 + k_ab * free_b
        if k_aa <= 0.0:
            return conc_a / factor
        return (-factor + math.sqrt(factor * factor + 8.0 * k_aa * conc_a)) / (4.0 * k_aa)

    def residual(free_b: float) -> float:
        free_a = free_a_given_b(free_b)
        return free_b + k_ab * free_a * free_b + 2.0 * k_bb * free_b * free_b - conc_b

    low, high = 0.0, conc_b  # residual(0) <= 0, residual(conc_b) >= 0, monotone
    for _ in range(200):
        middle = 0.5 * (low + high)
        if residual(middle) < 0.0:
            low = middle
        else:
            high = middle
        if high - low <= max(conc_b * 1e-14, 1e-24):
            break
    free_b = 0.5 * (low + high)
    free_a = free_a_given_b(free_b)
    ab = k_ab * free_a * free_b
    aa = k_aa * free_a * free_a
    bb = k_bb * free_b * free_b
    return EquilibriumResult(
        conc_a=conc_a,
        conc_b=conc_b,
        free_a=free_a,
        free_b=free_b,
        ab=ab,
        aa=aa,
        bb=bb,
    )


def interpret_equilibrium(
    result: EquilibriumResult,
    *,
    conc_nm: float,
) -> list[str]:
    """Descriptive interpretation of one equilibrium calculation."""

    n_m = 1e9
    return [
        f"在 {conc_nm:g} nM 起始浓度下：AB 异源二聚体 ≈ {result.ab * n_m:.1f} nM，"
        f"A 链自二聚体 AA ≈ {result.aa * n_m:.1f} nM，"
        f"B 链自二聚体 BB ≈ {result.bb * n_m:.1f} nM。",
        f"平衡时游离单链：A {result.free_a / result.conc_a * 100:.1f}%、"
        f"B {result.free_b / result.conc_b * 100:.1f}%；"
        f"A 链有 {result.a_occupied_fraction * 100:.1f}% 处于二聚体状态，"
        f"B 链有 {result.b_occupied_fraction * 100:.1f}%。",
        "该计算基于 RNAcofold 集合自由能（允许链内与链间配对），"
        "与 RNAduplex 仅链间模型数值略有差异；理想溶液假设，不含 Mg²⁺ 与 dNTP。",
    ]
