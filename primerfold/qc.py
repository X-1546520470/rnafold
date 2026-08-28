"""Sequence-level primer QC metrics computed locally (no ViennaRNA needed).

Tm uses the SantaLucia 1998 unified nearest-neighbour parameters with the
salt correction recommended therein (ΔS += 0.368·(N−1)·ln[Na⁺]). The duplex
is assumed non-self-complementary (Ct/4). All values are descriptive
reference points, not pass/fail judgements.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

from .core import Primer

# Nearest-neighbour enthalpy (kcal/mol) and entropy (cal/(mol·K)), SantaLucia 1998.
_NN_DH = {
    "AA": -7.9, "TT": -7.9,
    "AT": -7.2,
    "TA": -7.2,
    "CA": -8.5, "TG": -8.5,
    "GT": -8.4, "AC": -8.4,
    "CT": -7.8, "AG": -7.8,
    "GA": -8.2, "TC": -8.2,
    "CG": -10.6,
    "GC": -9.8,
    "GG": -8.0, "CC": -8.0,
}
_NN_DS = {
    "AA": -22.2, "TT": -22.2,
    "AT": -20.4,
    "TA": -21.3,
    "CA": -22.7, "TG": -22.7,
    "GT": -22.4, "AC": -22.4,
    "CT": -21.0, "AG": -21.0,
    "GA": -22.2, "TC": -22.2,
    "CG": -27.2,
    "GC": -24.4,
    "GG": -19.9, "CC": -19.9,
}
_INIT_DH = 0.2        # kcal/mol
_INIT_DS = -5.7       # cal/(mol·K)
_TERMINAL_AT_DH = 2.2  # per terminal A·T pair
_TERMINAL_AT_DS = 6.9
_SALT_DS_COEFF = 0.368  # per stacking step, multiplied by ln[Na+]
_GAS_CONSTANT = 1.9872041  # cal/(mol·K)

_DNA_MW = {"A": 313.21, "T": 304.2, "C": 289.18, "G": 329.21}
_DNA_MW_OFFSET = 61.96  # dephosphorylated ssDNA correction


def tm_nearest_neighbor(
    sequence: str,
    *,
    na_m: float = 0.05,
    ct_m: float = 2.5e-7,
) -> float:
    """Melting temperature in °C (SantaLucia 1998 NN + Mon Na⁺ correction)."""

    seq = sequence.upper()
    if len(seq) < 2:
        raise ValueError("Tm 计算至少需要 2 nt 序列。")
    for base in seq:
        if base not in "ACGT":
            raise ValueError(f"Tm 计算不支持碱基：{base}")

    dh = _INIT_DH
    ds = _INIT_DS
    for first, second in zip(seq, seq[1:]):
        dh += _NN_DH[first + second]
        ds += _NN_DS[first + second]

    terminal_at = sum(1 for base in (seq[0], seq[-1]) if base in "AT")
    dh += _TERMINAL_AT_DH * terminal_at
    ds += _TERMINAL_AT_DS * terminal_at

    ds += _SALT_DS_COEFF * (len(seq) - 1) * math.log(na_m)
    tm_k = (1000.0 * dh) / (ds + _GAS_CONSTANT * math.log(ct_m / 4.0))
    return tm_k - 273.15


def gc_percent(sequence: str) -> float:
    seq = sequence.upper()
    return 100.0 * sum(1 for base in seq if base in "GC") / len(seq)


def molecular_weight(sequence: str) -> float:
    """Single-stranded DNA molecular weight in Da (dephosphorylated)."""

    return sum(_DNA_MW[base] for base in sequence.upper()) - _DNA_MW_OFFSET


def gc_clamp_3p(sequence: str, *, last_n: int = 5) -> int:
    """Number of G/C bases in the 3'-terminal ``last_n`` nucleotides."""

    tail = sequence.upper()[-last_n:]
    return sum(1 for base in tail if base in "GC")


def max_homopolymer(sequence: str) -> tuple[int, str]:
    """Longest run of one repeated base as (length, base)."""

    best_length = 0
    best_base = ""
    current_length = 0
    current_base = ""
    for base in sequence.upper():
        if base == current_base:
            current_length += 1
        else:
            current_base = base
            current_length = 1
        if current_length > best_length:
            best_length = current_length
            best_base = current_base
    return best_length, best_base


@dataclass(frozen=True, slots=True)
class PrimerQc:
    name: str
    sequence: str
    length: int
    gc_percent: float
    tm_c: float
    molecular_weight: float
    gc_clamp_3p: int
    max_homopolymer_length: int
    max_homopolymer_base: str


def analyze_primer_qc(
    primer: Primer,
    *,
    na_m: float = 0.05,
    ct_m: float = 2.5e-7,
) -> PrimerQc:
    run_length, run_base = max_homopolymer(primer.sequence)
    return PrimerQc(
        name=primer.name,
        sequence=primer.sequence,
        length=len(primer.sequence),
        gc_percent=gc_percent(primer.sequence),
        tm_c=tm_nearest_neighbor(primer.sequence, na_m=na_m, ct_m=ct_m),
        molecular_weight=molecular_weight(primer.sequence),
        gc_clamp_3p=gc_clamp_3p(primer.sequence),
        max_homopolymer_length=run_length,
        max_homopolymer_base=run_base,
    )


def qc_rows(primers: list[Primer], *, na_m: float = 0.05) -> list[dict[str, object]]:
    return [
        {
            "引物": qc.name,
            "长度 (nt)": qc.length,
            "GC (%)": round(qc.gc_percent, 1),
            "Tm (°C)": round(qc.tm_c, 1),
            "分子量 (Da)": round(qc.molecular_weight, 1),
            "3′GC clamp": qc.gc_clamp_3p,
            "最长同聚串": f"{qc.max_homopolymer_base}×{qc.max_homopolymer_length}",
        }
        for qc in (analyze_primer_qc(primer, na_m=na_m) for primer in primers)
    ]
