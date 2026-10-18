# Reconstructed versions V55, V56 and V57 (05/07/2026)

Source: Antigravity conversation `3a8ccbb7` (and `551fff93`), in
`archive/transcripts/antigravity_2026-07/conversations/`. Step numbers are line
indices of the `*_transcript_full.jsonl` files.

| File | How it was rebuilt | Check |
|---|---|---|
| `compass_simV55_original.py` | `archive/OLD/compass_simV55.py` with edit 70 of `3a8ccbb7` (05/07 15:43) reversed | Byte-identical to the full file the agent wrote when it created V55 (`551fff93`, line 862, 15:25:46) |
| `compass_simV56_original.py` | `archive/OLD/compass_simV56.py` with edits 357 (16:26) and 30 (15:41) of `3a8ccbb7` reversed | Byte-identical to the full file written at creation (`551fff93`, line 901, 15:32:35) |
| `compass_simV57.py` | `archive/OLD/compass_simV58.py` with edit 433 of `3a8ccbb7` (19:52, the only edit between V57 and V58) reversed | Independent path agrees: archived V56 + edit 357 differs from it only in the version strings and two comment lines, as expected from edits 357/363 |

Why the archived files are not the originals: on 05/07 the agent edited the frozen
snapshots in place. It applied the "dial radius +30 %, centred arrow" change to
`compass_simV55.py` and `compass_simV56.py` when the user asked to "return to V55
and disconsider V56", and the bold-arrow change of V57 to `compass_simV56.py`.
`compass_simV57.py` was never written: the sandbox blocked `cp` into the Dropbox
folder, the agent wrote `freeze_v57.py` to its own scratch folder, and the session
ended with an error before running it. From V58 on, every snapshot was made with a
`freeze_vNN.py` helper (`archive/OLD/freeze_v58.py` … `freeze_v77.py`).

Times are local (UTC−3). The archived files in `archive/OLD/` are kept unchanged.
