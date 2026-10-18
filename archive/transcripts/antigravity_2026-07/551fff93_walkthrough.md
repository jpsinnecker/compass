# Walkthrough — Version 56 Updates

This walkthrough documents the transition from Version 55 to Version 56 of the simulator, and the report update to V22.

## Changes Made

### 1. HUD Dial Scaling and Minimum Length Calibration (V56)
- Modified [compass.py](file:///Users/jps/Library/CloudStorage/Dropbox/simula%20agulhas/compass.py) and [compass_simV56.py](file:///Users/jps/Library/CloudStorage/Dropbox/simula%20agulhas/compass_simV56.py):
  - Defined a minimum arrow length ($0.35 R_{\mathrm{dial}}$) for active fields so that the shaft remains longer than the arrowhead at small fields, reaching full radius ($1.0 R_{\mathrm{dial}}$) at maximum field.

### 2. Report Update to V22
- Created [relatorio_projeto_V22_EN.tex](file:///Users/jps/Library/CloudStorage/Dropbox/simula%20agulhas/relatorio_projeto_V22_EN.tex).
- Fixed undefined LaTeX macro warnings and commands (`\figref` $\to$ `\ref` and `\lastfoot` $\to$ `\endlastfoot`).
- Compiled the report to PDF inside the scratch folder successfully.

### 3. PDF Copy Helper
- Updated [copy_pdf.py](file:///Users/jps/Library/CloudStorage/Dropbox/simula%20agulhas/copy_pdf.py) to point to the new V22 PDF file.
