"""PrimerFold GUI core package."""

__version__ = "1.0.0"

from .core import (
    AnalysisJob,
    CofoldEnergies,
    FoldExecutionError,
    FoldResult,
    Primer,
    PrimerInputError,
    base_pair_coordinates,
    build_jobs,
    cofold_energies,
    fold_hairpin,
    fold_dimer,
    parse_primers,
    render_structure_svg,
)
from .equilibrium import (
    EquilibriumResult,
    equilibrium_constants,
    interpret_equilibrium,
    species_equilibrium,
)
from .plot import (
    parse_geometry,
    render_styled_structure_svg,
    svg_pixel_dimensions,
)
from .present import (
    alignment_html,
    alignment_legend_html,
    build_html_report,
    heatmap_html,
    interpret_result,
)
from .qc import (
    PrimerQc,
    analyze_primer_qc,
    qc_rows,
)

__all__ = [
    "AnalysisJob",
    "CofoldEnergies",
    "EquilibriumResult",
    "FoldExecutionError",
    "FoldResult",
    "Primer",
    "PrimerQc",
    "PrimerInputError",
    "alignment_html",
    "alignment_legend_html",
    "analyze_primer_qc",
    "base_pair_coordinates",
    "build_html_report",
    "build_jobs",
    "cofold_energies",
    "equilibrium_constants",
    "fold_hairpin",
    "fold_dimer",
    "heatmap_html",
    "interpret_equilibrium",
    "interpret_result",
    "parse_geometry",
    "parse_primers",
    "qc_rows",
    "render_styled_structure_svg",
    "render_structure_svg",
    "species_equilibrium",
    "svg_pixel_dimensions",
]
