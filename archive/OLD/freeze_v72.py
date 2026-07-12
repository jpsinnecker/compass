#!/usr/bin/env python3
"""
freeze_v72.py — Creates frozen snapshot compass_simV72.py from compass.py
"""
import shutil
import os

ws = os.path.dirname(os.path.abspath(__file__))
src = os.path.join(ws, "..", "compass.py")
dst = os.path.join(ws, "compass_simV72.py")
shutil.copy2(src, dst)
print(f"✓ Created: {dst}")
