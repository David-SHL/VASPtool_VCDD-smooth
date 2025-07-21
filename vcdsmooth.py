#!/usr/bin/env python3
import numpy as np
from pymatgen.io.vasp import Chgcar

# ---------------------
# Interactive user input
# ---------------------
filename = input("Enter the CHGCAR difference filename (e.g., CHGCAR_diff.vasp): ").strip()
atom_str = input("Enter target atom indices (starting from 0, comma-separated, e.g., 4,5,8): ").strip()
sigma_str = input("Enter sigma value (in Å, default = 2.0): ").strip()

# Parse inputs
target_indices = [int(i) for i in atom_str.split(',') if i.strip().isdigit()]
sigma = float(sigma_str) if sigma_str else 2.0
output_filename = "CHGCAR_diff_smooth.vasp"

# ---------------------
# Load CHGCAR and structure
# ---------------------
print(f"Loading file: {filename}")
chgcar = Chgcar.from_file(filename)
chg_data = chgcar.data['total']
structure = chgcar.structure
grid = np.array(chg_data.shape)
lattice = structure.lattice

# Validate atom indices
n_atoms = len(structure)
print(f"Number of atoms in structure: {n_atoms}")
for idx in target_indices:
    if idx >= n_atoms:
        raise IndexError(f"Atom index {idx} is out of bounds (max index = {n_atoms - 1})")

# ---------------------
# Create Gaussian weight mask
# ---------------------
print("Generating Gaussian weight mask...")
total_weights = np.zeros(grid)
coords = np.indices(grid).reshape(3, -1).T
real_coords = np.dot(coords / grid, lattice.matrix)  # shape: (N, 3)

for idx in target_indices:
    center = structure.cart_coords[idx]
    distances = np.linalg.norm(real_coords - center, axis=1)
    weights = np.exp(- (distances / sigma) ** 2)
    weights = weights.reshape(grid)
    total_weights += weights

# Normalize weights
max_weight = np.max(total_weights)
if max_weight > 0:
    total_weights /= max_weight
else:
    raise ValueError("All weights are zero. Check atom indices or sigma value.")

# ---------------------
# Apply weights and save smoothed CHGCAR
# ---------------------
chg_data_smooth = chg_data * total_weights
chgcar.data['total'] = chg_data_smooth
chgcar.write_file(output_filename)
print(f"✅ Smoothed charge difference written to: {output_filename}")
