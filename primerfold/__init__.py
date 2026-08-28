"""PrimerFold GUI core package."""

__version__ = "1.0.0"

from .core import (
    AnalysisJob,
    FoldExecutionError,
    FoldResult,
    Primer,
    PrimerInputError,
    base_pair_coordinates,
    build_jobs,
    fold_hairpin,
    fold_dimer,
    parse_primers,
    render_structure_svg,
)
from .plot import (
    parse_geometry,
    render_styled_structure_svg,
    svg_pixel_dimensions,
)
from .present import (
    alignment_html,
    alignment_legend_html,
    interpret_result,
)

__all__ = [
    "AnalysisJob",
    "FoldExecutionError",
    "FoldResult",
    "Primer",
    "PrimerInputError",
    "alignment_html",
    "alignment_legend_html",
    "base_pair_coordinates",
    "build_jobs",
    "fold_hairpin",
    "fold_dimer",
    "interpret_result",
    "parse_geometry",
    "parse_primers",
    "render_styled_structure_svg",
    "render_structure_svg",
    "svg_pixel_dimensions",
]
