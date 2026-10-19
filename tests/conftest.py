import sys
from pathlib import Path

# Defensive: make `import cases` / `import helpers` work regardless of how
# pytest was invoked (rootdir, -p no:cacheprovider, etc.) or which import
# mode is configured.
_TESTS_DIR = Path(__file__).resolve().parent
if str(_TESTS_DIR) not in sys.path:
    sys.path.insert(0, str(_TESTS_DIR))


# Slow physics checks (tens of seconds each) are skipped unless --runslow is
# given, so the default suite stays fast:  python3 -m pytest -q tests --runslow
import pytest  # noqa: E402


def pytest_addoption(parser):
    parser.addoption("--runslow", action="store_true", default=False, help="also run tests marked slow")


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: slow physics check, run with --runslow")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--runslow"):
        return
    skip = pytest.mark.skip(reason="slow; use --runslow")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip)

