"""sim_config.py — loader for config.yaml.

Central place that reads config.yaml and exposes its contents as typed,
dot-accessible dataclasses (Config.physics.compass_engine.damping, etc.).

Every script in this repo that previously hardcoded a physical constant,
numerical default, or run/output setting now sources that default from here
instead. Values are unchanged; this module only removes duplication of the
literals across files.

Usage:
    from sim_config import load_config
    CFG = load_config()
    parser.add_argument("--damping", type=float,
                         default=CFG.physics.compass_engine.damping)
"""

from __future__ import annotations

from dataclasses import dataclass, fields, is_dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, get_type_hints

import yaml

_DEFAULT_CONFIG_PATH = Path(__file__).resolve().parent / "config.yaml"


# =============================================================================
# physics
# =============================================================================


@dataclass
class Constants:
    mu0_over_4pi: float


@dataclass
class ForcPhysics:
    n_curves: int
    t_sat: float
    t_ramp_down: float
    t_ramp_up: float


@dataclass
class DemagPhysics:
    freq: float
    cycles: int
    t_relax_after: float


@dataclass
class CompassEnginePhysics:
    steel_density: float
    steel_ms_saturation: float
    center_distance: float
    needle_len: float
    needle_width: float
    needle_thickness: float
    needle_frac_legacy: float
    use_legacy_size_from_R: int
    damping: float
    damping_noise: float
    pivot_radius: float
    pivot_thickness: float
    pivot_density: float
    field_mode: str
    B_max_factor: float
    phi_ext_deg: float
    field_freq: float
    field_delay: float
    hyst_spacing: str
    hyst_log_k: float
    t_sim: float
    noise: float
    forc: ForcPhysics
    demag: DemagPhysics


@dataclass
class CalcDefaultsTempPhysics:
    density: float
    thickness: float
    needle_len: float
    needle_width: float
    pivot_radius: float
    pivot_thickness: float
    pivot_density: float
    Ms: float


@dataclass
class DampingSweepPhysics:
    R_default: float
    needle_frac: float
    needle_thickness: float
    steel_density: float
    pivot_radius: float
    pivot_thickness: float
    pivot_density: float
    B_max_factor: float
    noise: float
    cutoff_default: float


@dataclass
class SweepDampingHysteresisPhysics:
    R_default: float
    needle_frac: float
    needle_thickness: float
    steel_density: float
    pivot_radius: float
    pivot_thickness: float
    needle_width_to_length_ratio: float
    B_max_factor: float
    noise_init: float


@dataclass
class NeedleRenderPhysics:
    colors: Dict[str, str]
    pivot_radius_frac: float
    pivot_inner_radius_frac: float


@dataclass
class CompassGenerateImagesPhysics:
    r_nn_fallback: float
    needle_len_to_r_nn_fallback_ratio: float
    needle_width_to_length_fallback_ratio: float


@dataclass
class PhysicsConfig:
    constants: Constants
    compass_engine: CompassEnginePhysics
    calc_defaults_temp: CalcDefaultsTempPhysics
    damping_sweep: DampingSweepPhysics
    sweep_damping_hysteresis: SweepDampingHysteresisPhysics
    compass_generate_images: CompassGenerateImagesPhysics
    needle_render: NeedleRenderPhysics


# =============================================================================
# numerics
# =============================================================================


@dataclass
class CompassEngineGrid:
    geometry: str
    N: int
    M: int
    cutoff_shells: float
    n_images: int
    tensor_mem_limit_gb: float
    float32: bool
    pbc: bool


@dataclass
class CompassEngineTime:
    dt_factor: float
    log_every: int


@dataclass
class CompassEngineTolerances:
    flip_angle_deg: float
    flip_band_deg: float
    flip_dwell_T0: float
    flip_settle_frac: float
    dt_guard_alpha: float
    domain_tol_deg: float


@dataclass
class CompassEngineNumerics:
    grid: CompassEngineGrid
    time: CompassEngineTime
    tolerances: CompassEngineTolerances


@dataclass
class DampingSweepQuickTest:
    n_seeds: int
    n_dampings: int
    t_sim_periods: float
    grid_n: int
    forc_n_curves: int


@dataclass
class DampingSweepNumerics:
    grid_n: int
    n_seeds: int
    n_dampings: int
    Q_min: float
    Q_max: float
    t_sim_periods: float
    dt_factor: float
    domain_tol: float
    forc_n_curves: int
    quick_test: DampingSweepQuickTest


@dataclass
class SweepDampingHysteresisNumerics:
    N_grid: int
    M_grid: int
    n_damping: int
    Q_min: float
    Q_max: float
    cutoff: float
    dt_factor: float
    t_sim_in_T0: float


@dataclass
class DampingSweepAnalysisNumerics:
    mad_threshold: float


@dataclass
class AvalancheProcessorNumerics:
    t_link_T0: str
    r_link_rnn: float
    min_tail: int


@dataclass
class RenderingNumerics:
    dpi_default: int
    figsize_default: float


@dataclass
class NumericsConfig:
    compass_engine: CompassEngineNumerics
    damping_sweep: DampingSweepNumerics
    sweep_damping_hysteresis: SweepDampingHysteresisNumerics
    damping_sweep_analysis: DampingSweepAnalysisNumerics
    avalanche_processor: AvalancheProcessorNumerics
    rendering: RenderingNumerics


# =============================================================================
# run
# =============================================================================


@dataclass
class CompassEngineRun:
    out_dir: str
    tag: Optional[str]
    seed: Optional[int]
    use_gpu: bool
    progress: bool
    verbose: bool
    make_plot: bool
    event_log: bool
    dt_guard_substep: bool
    png_dpi: int
    png_transparent: bool
    png_with_axes: bool
    png_no_panel_titles: bool


@dataclass
class DampingSweepRun:
    geometries: str
    pbc: int
    use_gpu: int
    damping_noise: float
    skip_relax_stage: int
    skip_forc_stage: int
    quick_test: int
    n_workers: int
    resume: int


@dataclass
class SweepDampingHysteresisRun:
    seeds: List[int]
    geometries: List[str]
    gpu: int
    only_index: Optional[int]


@dataclass
class AvalancheProcessorRun:
    out_dir: str
    channels: str
    group_by: str
    make_plot: bool


@dataclass
class CompassGenerateImagesRun:
    run_dir: str
    recursive: bool
    transparent: bool
    with_axes: bool
    no_panel_titles: bool
    verbose: bool


@dataclass
class PlotDefaultGeometriesCleanRun:
    input_dir: str
    input_file: Optional[str]
    output_file: Optional[str]
    output_dir: Optional[str]
    single: bool
    tag_suffix: str
    transparent: bool


@dataclass
class DampingSweepAnalysisRun:
    geometry_colors: Dict[str, str]
    geometry_markers: Dict[str, str]


@dataclass
class RunConfig:
    compass_engine: CompassEngineRun
    damping_sweep: DampingSweepRun
    sweep_damping_hysteresis: SweepDampingHysteresisRun
    avalanche_processor: AvalancheProcessorRun
    compass_generate_images: CompassGenerateImagesRun
    plot_default_geometries_clean: PlotDefaultGeometriesCleanRun
    damping_sweep_analysis: DampingSweepAnalysisRun


# =============================================================================
# top-level config + loader
# =============================================================================


@dataclass
class Config:
    physics: PhysicsConfig
    numerics: NumericsConfig
    run: RunConfig


def _build(cls: Any, data: Dict[str, Any]) -> Any:
    """Recursively construct a (possibly nested) dataclass from a dict.

    Every dataclass field must be present in `data`; this is intentional so
    that a truncated or misspelled config.yaml fails loudly at startup
    rather than silently falling back to some other default.

    Uses get_type_hints() rather than raw Field.type: this module uses
    `from __future__ import annotations`, so Field.type is only the
    unevaluated string ("PhysicsConfig"), not the class itself.
    """
    hints = get_type_hints(cls)
    kwargs = {}
    for f in fields(cls):
        if f.name not in data:
            raise KeyError(f"config.yaml is missing '{f.name}' required by {cls.__name__}")
        value = data[f.name]
        field_type = hints[f.name]
        if is_dataclass(field_type):
            kwargs[f.name] = _build(field_type, value)
        else:
            kwargs[f.name] = value
    return cls(**kwargs)


def load_config(path: Optional[str | Path] = None) -> Config:
    """Load config.yaml (or the given path) into a Config dataclass tree."""
    cfg_path = Path(path) if path is not None else _DEFAULT_CONFIG_PATH
    with open(cfg_path, "r") as fh:
        raw = yaml.safe_load(fh)
    return _build(Config, raw)
