"""Shared definition of the reference test cases.

Both generate_reference.py (offline: builds fixtures + reference outputs)
and test_reference_runs.py (pytest: re-runs + compares) import CASES from
here, so the two never drift apart.

Each case is a small, fast compassV2_2.py run: a tiny lattice, few logged
steps, a fixed seed, and an explicit tag. The `overrides` dict is applied on
top of a full copy of the repo's config.yaml to produce that case's own
variant config file (see build_variant_config in generate_reference.py /
conftest.py) -- this is deliberately a *config.yaml* variant, not a set of
CLI flags, so the tests exercise the same physics/numerics/run.compass_engine
path real users go through.
"""

CASES = [
    {
        "name": "square_static",
        "description": "3x3 square lattice, static field, plain relaxation.",
        "overrides": {
            "numerics": {
                "compass_engine": {
                    "grid": {"geometry": "square", "N": 3, "M": 3},
                    "time": {"dt_factor": 0.12, "log_every": 1},
                }
            },
            "physics": {
                "compass_engine": {
                    "field_mode": "static",
                    "t_sim": 0.02,
                }
            },
            "run": {
                "compass_engine": {
                    "seed": 12345,
                    "tag": "square_static_ref",
                    "png_dpi": 60,
                }
            },
        },
    },
    {
        "name": "triangular_hysteresis",
        "description": "3x3 triangular lattice, log-spaced hysteresis cycle.",
        "overrides": {
            "numerics": {
                "compass_engine": {
                    "grid": {"geometry": "triangular", "N": 3, "M": 3},
                    "time": {"dt_factor": 0.15, "log_every": 1},
                }
            },
            "physics": {
                "compass_engine": {
                    "field_mode": "hysteresis",
                    "t_sim": 0.02,
                    "hyst_spacing": "log",
                }
            },
            "run": {
                "compass_engine": {
                    "seed": 23456,
                    "tag": "triangular_hysteresis_ref",
                    "png_dpi": 60,
                }
            },
        },
    },
    {
        "name": "honeycomb_demag",
        "description": "3x3 honeycomb lattice, one rotational-demag cycle.",
        "overrides": {
            "numerics": {
                "compass_engine": {
                    "grid": {"geometry": "honeycomb", "N": 3, "M": 3},
                    "time": {"dt_factor": 0.3, "log_every": 1},
                }
            },
            "physics": {
                "compass_engine": {
                    "field_mode": "demag_rot",
                    "demag": {"freq": 20.0, "cycles": 1, "t_relax_after": 0.005},
                }
            },
            "run": {
                "compass_engine": {
                    "seed": 34567,
                    "tag": "honeycomb_demag_ref",
                    "png_dpi": 60,
                }
            },
        },
    },
]
