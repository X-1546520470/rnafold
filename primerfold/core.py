"""Validated DNA-primer analysis using ViennaRNA command-line tools."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
import os
import re
import shutil
import subprocess
import sys
import tempfile
from typing import Literal, Sequence


AnalysisKind = Literal["hairpin", "self_dimer", "cross_dimer"]

MAX_PRIMERS = 30
MAX_SEQUENCE_LENGTH = 500
MIN_SEQUENCE_LENGTH = 4
FOLD_TIMEOUT_SECONDS = 30


class PrimerInputError(ValueError):
    """Raised when user-supplied primer text cannot be validated."""


class FoldExecutionError(RuntimeError):
    """Raised when a ViennaRNA executable fails or returns malformed output."""


@dataclass(frozen=True, slots=True)
class Primer:
    name: str
    sequence: str


@dataclass(frozen=True, slots=True)
class AnalysisJob:
    kind: AnalysisKind
    primer_a: Primer
    primer_b: Primer | None = None

    @property
    def label(self) -> str:
        if self.kind == "hairpin":
            return f"{self.primer_a.name} · 发卡"
        if self.kind == "self_dimer":
            return f"{self.primer_a.name} · 自二聚体"
        if self.primer_b is None:
            raise ValueError("cross-dimer job is missing primer_b")
        return f"{self.primer_a.name} × {self.primer_b.name} · 交叉二聚体"


@dataclass(frozen=True, slots=True)
class FoldResult:
    kind: AnalysisKind
    label: str
    primer_a: Primer
    primer_b: Primer | None
    folded_sequence: str
    structure: str
    mfe_kcal_mol: float
    base_pairs: int
    intermolecular_pairs: int
    a_last5_paired: int
    b_last5_paired: int | None
    a_terminal_paired: bool
    b_terminal_paired: bool | None
    warnings: tuple[str, ...]


_STRUCTURE_LINE = re.compile(
    r"^([().,&]+)\s+\(\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+))\s*\)\s*$"
)
_DUPLEX_LINE = re.compile(
    r"^([().]+)&([().]+)\s+"
    r"(\d+)\s*,\s*(\d+)\s*:\s*(\d+)\s*,\s*(\d+)\s+"
    r"\(\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+))\s*\)\s*$"
)
_VALID_SEQUENCE = re.compile(r"^[ACGT]+$")
_HEADER_NAMES = {"name", "primer", "id", "名称", "引物"}
_HEADER_SEQUENCES = {"sequence", "seq", "序列"}


def _clean_sequence(raw: str) -> str:
    sequence = re.sub(r"[\s\-]", "", raw).upper()
    if not sequence:
        raise PrimerInputError("检测到空序列。")
    if not _VALID_SEQUENCE.fullmatch(sequence):
        invalid = "".join(sorted(set(sequence) - set("ACGT")))
        raise PrimerInputError(
            f"序列仅允许 A/C/G/T；检测到无效字符：{invalid or '未知'}。"
        )
    if len(sequence) < MIN_SEQUENCE_LENGTH:
        raise PrimerInputError(
            f"序列长度至少为 {MIN_SEQUENCE_LENGTH} nt；当前为 {len(sequence)} nt。"
        )
    if len(sequence) > MAX_SEQUENCE_LENGTH:
        raise PrimerInputError(
            f"单条序列最长支持 {MAX_SEQUENCE_LENGTH} nt；当前为 {len(sequence)} nt。"
        )
    return sequence


def _validate_primers(records: Sequence[tuple[str, str]]) -> list[Primer]:
    if not records:
        raise PrimerInputError("没有检测到引物序列。")
    if len(records) > MAX_PRIMERS:
        raise PrimerInputError(
            f"一次最多分析 {MAX_PRIMERS} 条引物；当前检测到 {len(records)} 条。"
        )

    primers: list[Primer] = []
    seen: set[str] = set()
    for index, (raw_name, raw_sequence) in enumerate(records, start=1):
        name = raw_name.strip() or f"Primer_{index}"
        if name in seen:
            raise PrimerInputError(f"引物名称重复：{name}。请为每条引物使用唯一名称。")
        seen.add(name)
        try:
            sequence = _clean_sequence(raw_sequence)
        except PrimerInputError as exc:
            raise PrimerInputError(f"{name}：{exc}") from exc
        primers.append(Primer(name=name, sequence=sequence))
    return primers


def _parse_fasta(lines: Sequence[str]) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    name: str | None = None
    chunks: list[str] = []

    for line in lines:
        if not line or line.startswith(";"):
            continue
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(chunks)))
            name = line[1:].strip().split(maxsplit=1)[0] if line[1:].strip() else ""
            chunks = []
            continue
        if name is None:
            raise PrimerInputError("FASTA 序列前缺少以 > 开头的名称行。")
        chunks.append(line)

    if name is not None:
        records.append((name, "".join(chunks)))
    return records


def _parse_delimited(lines: Sequence[str]) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    auto_index = 1

    for line in lines:
        if not line or line.startswith("#"):
            continue

        if "\t" in line:
            fields = [field.strip() for field in line.split("\t") if field.strip()]
        elif "," in line:
            fields = [field.strip() for field in line.split(",") if field.strip()]
        else:
            fields = line.split()

        if len(fields) >= 2:
            if (
                fields[0].lower() in _HEADER_NAMES
                and fields[-1].lower() in _HEADER_SEQUENCES
            ):
                continue
            records.append((fields[0], fields[-1]))
        elif len(fields) == 1:
            records.append((f"Primer_{auto_index}", fields[0]))
            auto_index += 1
        else:
            continue
    return records


def parse_primers(text: str) -> list[Primer]:
    """Parse FASTA, CSV/TSV, name+sequence lines, or sequence-only lines."""

    normalized = text.replace("\ufeff", "").replace("\r\n", "\n").replace("\r", "\n")
    lines = [line.strip() for line in normalized.split("\n")]
    nonempty = [line for line in lines if line]
    if not nonempty:
        raise PrimerInputError("请输入或上传至少一条引物序列。")

    if any(line.startswith(">") for line in nonempty):
        records = _parse_fasta(lines)
    else:
        records = _parse_delimited(lines)
    return _validate_primers(records)


def build_jobs(
    primers: Sequence[Primer],
    *,
    hairpin: bool,
    self_dimer: bool,
    cross_dimer: bool,
) -> list[AnalysisJob]:
    jobs: list[AnalysisJob] = []
    if hairpin:
        jobs.extend(AnalysisJob("hairpin", primer) for primer in primers)
    if self_dimer:
        jobs.extend(AnalysisJob("self_dimer", primer, primer) for primer in primers)
    if cross_dimer:
        jobs.extend(
            AnalysisJob("cross_dimer", primer_a, primer_b)
            for primer_a, primer_b in combinations(primers, 2)
        )
    if not jobs:
        raise PrimerInputError("请至少选择一种分析类型。")
    return jobs


def _find_binary(name: str) -> str:
    environment_candidate = Path(sys.executable).resolve().parent / name
    if environment_candidate.is_file() and os.access(environment_candidate, os.X_OK):
        return str(environment_candidate)
    discovered = shutil.which(name)
    if discovered:
        return discovered
    raise FoldExecutionError(
        f"找不到 {name}。请确认已激活 /opt/anaconda3/envs/rnafold 环境。"
    )


def _run_vienna(
    name: str,
    input_text: str,
    temperature_c: float,
    salt_m: float,
    *,
    no_ps: bool,
) -> tuple[str, tuple[str, ...]]:
    if not 0.0 <= temperature_c <= 100.0:
        raise PrimerInputError("温度必须在 0–100 °C 之间。")
    if not 0.0001 <= salt_m <= 2.0:
        raise PrimerInputError("单价盐浓度必须在 0.0001–2.0 M 之间。")

    command = [
        _find_binary(name),
        "--noconv",
        "--paramFile=DNA",
        f"--temp={temperature_c:.3f}",
        f"--salt={salt_m:.6f}",
    ]
    if no_ps:
        command.insert(1, "--noPS")
    environment = os.environ.copy()
    environment["LC_ALL"] = "C"
    environment["PYTHONNOUSERSITE"] = "1"
    try:
        completed = subprocess.run(
            command,
            input=input_text,
            text=True,
            capture_output=True,
            timeout=FOLD_TIMEOUT_SECONDS,
            env=environment,
            check=False,
        )
    except subprocess.TimeoutExpired as exc:
        raise FoldExecutionError(f"{name} 运行超过 {FOLD_TIMEOUT_SECONDS} 秒。") from exc
    except OSError as exc:
        raise FoldExecutionError(f"无法启动 {name}：{exc}") from exc

    if completed.returncode != 0:
        detail = completed.stderr.strip() or completed.stdout.strip() or "无错误详情"
        raise FoldExecutionError(f"{name} 运行失败：{detail}")

    warnings = tuple(
        line.strip()
        for line in completed.stderr.splitlines()
        if line.strip()
    )
    return completed.stdout, warnings


def _parse_fold_output(output: str) -> tuple[str, float]:
    for line in reversed(output.splitlines()):
        match = _STRUCTURE_LINE.match(line.strip())
        if match:
            return match.group(1), float(match.group(2))
    raise FoldExecutionError(f"无法解析 ViennaRNA 输出：\n{output.strip()}")


def _parse_duplex_output(
    output: str,
    primer_a: Primer,
    primer_b: Primer,
) -> tuple[str, float]:
    for line in reversed(output.splitlines()):
        match = _DUPLEX_LINE.match(line.strip())
        if not match:
            continue
        a_segment, b_segment = match.group(1), match.group(2)
        a_start, a_end = int(match.group(3)), int(match.group(4))
        b_start, b_end = int(match.group(5)), int(match.group(6))
        mfe = float(match.group(7))

        # RNAduplex uses 100000 kcal/mol and a 0,0 range as a sentinel for
        # "no intermolecular base pair found". Represent the unbound state as
        # an all-unpaired structure with zero interaction energy in the GUI.
        if mfe >= 99999 or min(a_start, a_end, b_start, b_end) == 0:
            return (
                "." * len(primer_a.sequence)
                + "&"
                + "." * len(primer_b.sequence),
                0.0,
            )

        if len(a_segment) != a_end - a_start + 1:
            raise FoldExecutionError("RNAduplex 返回的 A 链范围与结构长度不一致。")
        if len(b_segment) != b_end - b_start + 1:
            raise FoldExecutionError("RNAduplex 返回的 B 链范围与结构长度不一致。")
        if not (1 <= a_start <= a_end <= len(primer_a.sequence)):
            raise FoldExecutionError("RNAduplex 返回的 A 链范围越界。")
        if not (1 <= b_start <= b_end <= len(primer_b.sequence)):
            raise FoldExecutionError("RNAduplex 返回的 B 链范围越界。")

        a_structure = (
            "." * (a_start - 1)
            + a_segment
            + "." * (len(primer_a.sequence) - a_end)
        )
        b_structure = (
            "." * (b_start - 1)
            + b_segment
            + "." * (len(primer_b.sequence) - b_end)
        )
        return f"{a_structure}&{b_structure}", mfe
    raise FoldExecutionError(f"无法解析 RNAduplex 输出：\n{output.strip()}")


def base_pair_coordinates(structure: str) -> list[tuple[tuple[int, int], tuple[int, int]]]:
    """Decode a dot-bracket structure into ((strand, pos), (strand, pos)) pairs."""

    pairs: list[tuple[tuple[int, int], tuple[int, int]]] = []
    stack: list[tuple[int, int]] = []
    strand = 0
    position = 0

    for character in structure:
        if character == "&":
            strand += 1
            position = 0
            continue
        if character == "(":
            stack.append((strand, position))
        elif character == ")":
            if not stack:
                raise FoldExecutionError("ViennaRNA 返回了不平衡的点括号结构。")
            pairs.append((stack.pop(), (strand, position)))
        elif character not in ".,":
            raise FoldExecutionError(f"ViennaRNA 返回了未知结构符号：{character}")
        position += 1

    if stack:
        raise FoldExecutionError("ViennaRNA 返回了不平衡的点括号结构。")
    return pairs


def _strand_metrics(
    sequence: str,
    pairs: Sequence[tuple[tuple[int, int], tuple[int, int]]],
    strand_index: int,
    *,
    intermolecular_only: bool,
) -> tuple[int, bool]:
    lengths = [len(part) for part in sequence.split("&")]
    if strand_index >= len(lengths):
        raise FoldExecutionError("结构中的链数与序列不一致。")

    paired_positions: set[int] = set()
    for first, second in pairs:
        if intermolecular_only and first[0] == second[0]:
            continue
        if first[0] == strand_index:
            paired_positions.add(first[1])
        if second[0] == strand_index:
            paired_positions.add(second[1])

    length = lengths[strand_index]
    last_five = set(range(max(0, length - 5), length))
    return len(paired_positions & last_five), (length - 1) in paired_positions


def _make_result(
    job: AnalysisJob,
    folded_sequence: str,
    structure: str,
    mfe: float,
    warnings: tuple[str, ...],
) -> FoldResult:
    pairs = base_pair_coordinates(structure)
    intermolecular = sum(first[0] != second[0] for first, second in pairs)
    is_dimer = job.kind != "hairpin"
    a_last5, a_terminal = _strand_metrics(
        folded_sequence,
        pairs,
        0,
        intermolecular_only=is_dimer,
    )

    b_last5: int | None = None
    b_terminal: bool | None = None
    if is_dimer:
        b_last5, b_terminal = _strand_metrics(
            folded_sequence,
            pairs,
            1,
            intermolecular_only=True,
        )

    return FoldResult(
        kind=job.kind,
        label=job.label,
        primer_a=job.primer_a,
        primer_b=job.primer_b,
        folded_sequence=folded_sequence,
        structure=structure,
        mfe_kcal_mol=mfe,
        base_pairs=len(pairs),
        intermolecular_pairs=intermolecular,
        a_last5_paired=a_last5,
        b_last5_paired=b_last5,
        a_terminal_paired=a_terminal,
        b_terminal_paired=b_terminal,
        warnings=warnings,
    )


def fold_hairpin(primer: Primer, *, temperature_c: float, salt_m: float) -> FoldResult:
    job = AnalysisJob("hairpin", primer)
    output, warnings = _run_vienna(
        "RNAfold",
        f">primerfold\n{primer.sequence}\n",
        temperature_c,
        salt_m,
        no_ps=True,
    )
    structure, mfe = _parse_fold_output(output)
    return _make_result(job, primer.sequence, structure, mfe, warnings)


def fold_dimer(
    primer_a: Primer,
    primer_b: Primer,
    *,
    kind: Literal["self_dimer", "cross_dimer"],
    temperature_c: float,
    salt_m: float,
) -> FoldResult:
    job = AnalysisJob(kind, primer_a, primer_b)
    folded_sequence = f"{primer_a.sequence}&{primer_b.sequence}"
    output, warnings = _run_vienna(
        "RNAduplex",
        f"{primer_a.sequence}\n{primer_b.sequence}\n",
        temperature_c,
        salt_m,
        no_ps=False,
    )
    structure, mfe = _parse_duplex_output(output, primer_a, primer_b)
    return _make_result(job, folded_sequence, structure, mfe, warnings)


def render_structure_svg(sequence: str, structure: str) -> str:
    """Render a ViennaRNA structure to an SVG string with RNAplot."""

    binary = _find_binary("RNAplot")
    environment = os.environ.copy()
    environment["LC_ALL"] = "C"
    environment["PYTHONNOUSERSITE"] = "1"
    input_text = f">primerfold_structure\n{sequence}\n{structure}\n"

    with tempfile.TemporaryDirectory(prefix="primerfold-svg-") as temporary:
        try:
            completed = subprocess.run(
                [binary, "--output-format=svg", "--jobs=1"],
                input=input_text,
                text=True,
                capture_output=True,
                timeout=FOLD_TIMEOUT_SECONDS,
                env=environment,
                cwd=temporary,
                check=False,
            )
        except (subprocess.TimeoutExpired, OSError) as exc:
            raise FoldExecutionError(f"RNAplot 无法生成结构图：{exc}") from exc

        if completed.returncode != 0:
            detail = completed.stderr.strip() or completed.stdout.strip() or "无错误详情"
            raise FoldExecutionError(f"RNAplot 生成结构图失败：{detail}")

        svg_path = Path(temporary) / "primerfold_structure_ss.svg"
        if not svg_path.is_file():
            raise FoldExecutionError("RNAplot 未生成预期的 SVG 文件。")
        svg = svg_path.read_text(encoding="utf-8")

    # RNAplot embeds a click handler that is unnecessary inside the local GUI.
    svg = re.sub(r"\s*<script\b.*?</script>\s*", "\n", svg, flags=re.DOTALL)
    svg = re.sub(r'\s+onclick="[^"]*"', "", svg)
    return svg
