#!/usr/bin/env python3
"""
freeze_v74.py — Creates frozen snapshot compass_simV74.py from compass.py
"""
import shutil
import os

ws = os.path.dirname(os.path.abspath(__file__))
src = os.path.join(ws, "..", "compass.py")
dst = os.path.join(ws, "compass_simV74.py")
shutil.copy2(src, dst)
print(f"✓ Created: {dst}")
