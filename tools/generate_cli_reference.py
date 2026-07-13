#!/usr/bin/env python3
"""tools/generate_cli_reference.py

Introspects the canonical engine's build_parser() and emits a LaTeX
longtable body (rows only, grouped by argparse argument group) describing
every CLI flag: name, default value (as it actually is in the running
code), and its argparse `help=` string.

Purpose: docs/project_report.tex's Appendix A (the CLI reference table)
has drifted from the code before (it was missing the six V2.1 flip/dt_guard
flags -- see docs/AUDIT.md). This script lets anyone re-check -- or
regenerate -- that table directly from build_parser(), so flag names and
defaults can be verified mechanically rather than by manual cross-reading.

Known limitation (see the report's own note next to Appendix A): roughly
a third of the flags in compassV2_2.py's build_parser() have no `help=`
string in the source at all. For those, this script emits a placeholder
rather than inventing prose. Appendix A's DESCRIPTIONS for those flags are
therefore still hand-authored and can in
principle drift from actual behavior; only flag NAMES and DEFAULTS are
guaranteed correct by construction here. Closing that gap fully would mean
adding `help=` strings to every remaining argparse call in compassV2_2.py
(a code change, out of scope for a documentation-generation tool).

Usage:
    python3 tools/generate_cli_reference.py            # print to stdout
    python3 tools/generate_cli_reference.py --check    # exit 1 if any flag
                                                        # is missing help=
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

import compassV2_2  # noqa: E402  (path insert must happen first)

_NO_HELP_PLACEHOLDER = r"\textit{(no CLI help string in source; see code)}"

# Skip argparse's own auto-added group.
_SKIP_GROUPS = {"options", "positional arguments", "optional arguments"}


def _format_default(default) -> str:
    if default is None:
        return r"\code{None}"
    if isinstance(default, bool):
        return r"\code{True}" if default else r"\code{False}"
    if isinstance(default, float):
        # Match the report's style: plain decimal for "nice" numbers,
        # scientific notation only when that is how the value would
        # otherwise render (e.g. 5e-08).
        text = repr(default)
        return f"\\code{{{text}}}"
    return f"\\code{{{default}}}"


def _escape_latex(text: str) -> str:
    return (
        text.replace("\\", r"\textbackslash{}")
        .replace("_", r"\_")
        .replace("%", r"\%")
        .replace("#", r"\#")
    )


def iter_flags(parser: argparse.ArgumentParser):
    """Yield (group_title, flag, default, help_text, choices) in
    build_parser()'s own declared order."""
    for group in parser._action_groups:
        if group.title in _SKIP_GROUPS or not group._group_actions:
            continue
        for action in group._group_actions:
            flag = max(action.option_strings, key=len)  # prefer --long_form
            yield group.title, flag, action.default, action.help, action.choices


def build_longtable_rows(parser: argparse.ArgumentParser) -> str:
    lines = []
    last_group = None
    for group_title, flag, default, help_text, choices in iter_flags(parser):
        if group_title != last_group:
            if last_group is not None:
                lines.append("")
            lines.append(f"% --- {group_title} ---")
            last_group = group_title

        default_str = _format_default(default)
        if help_text:
            desc = _escape_latex(help_text)
        else:
            desc = _NO_HELP_PLACEHOLDER
        if choices:
            choice_str = ", ".join(f"\\code{{{c}}}" for c in choices)
            desc = f"{desc} (choices: {choice_str})" if help_text else f"Choices: {choice_str}."
        lines.append(f"\\code{{{flag}}} & {desc} Default: {default_str}.\\\\")
    return "\n".join(lines)


def check_missing_help(parser: argparse.ArgumentParser) -> list[str]:
    return [flag for _, flag, _, help_text, _ in iter_flags(parser) if not help_text]


def main() -> int:
    parser = compassV2_2.build_parser()

    if "--check" in sys.argv:
        missing = check_missing_help(parser)
        if missing:
            print(f"{len(missing)} flag(s) with no help= string in build_parser():")
            for flag in missing:
                print(f"  {flag}")
            return 1
        print("All flags have a help= string.")
        return 0

    print(build_longtable_rows(parser))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
