# VASPtool_VCDD-smooth

This script allows users to extract the smoothed local charge difference around specific atoms in a system, based on VASP CHGCAR difference files. It uses Gaussian-based masking to smoothly highlight the charge density changes near selected atoms, providing cleaner isosurface visualization for analysis in VESTA or other tools.

---

## Features

- Applies a **Gaussian-shaped spatial mask** around user-selected atoms
- Accepts VASP `CHGCAR`-format charge difference files as input
- Allows interactive input for:
  - Target atom indices (0-based)
  - Smoothing factor `sigma` (Å)
- Outputs a modified `CHGCAR_diff_smooth.vasp` file with masked charge density
- Improves clarity of charge difference isosurfaces and minimizes visual noise

---

## Requirements

- Python 3.6+
- [pymatgen](https://pymatgen.org/)
- numpy

---

## Usage

1. Place your VASP difference file (e.g., `CHGCAR_diff.vasp`) in the same directory as the script.
2. Run the script:
   ```bash
   python vcdsmooth.py
