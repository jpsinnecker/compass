import sys
from pathlib import Path

# Defensive: make `import cases` / `import helpers` work regardless of how
# pytest was invoked (rootdir, -p no:cacheprovider, etc.) or which import
# mode is configured.
_TESTS_DIR = Path(__file__).resolve().parent
if str(_TESTS_DIR) not in sys.path:
    sys.path.insert(0, str(_TESTS_DIR))
