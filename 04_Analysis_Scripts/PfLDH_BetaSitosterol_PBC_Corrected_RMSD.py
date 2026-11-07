from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.lib.distances import minimize_vectors

# Original simulation files
root = Path.home() / "Downloads" / "PfLDH_Unrestrained_NPT_100ps_Backup"
topology = root / "PfLDH_NAD_BetaSitosterol_Solvated.prmtop"

trajectories = [
    str(root / "Unrestrained_NPT" / f"NPT_segment_{i:02d}.dcd")
    for i in range(1, 21)
]

# Master project results folder
output = (
    Path.home() / "Documents" /
    "PfLDH_BetaSitosterol_MD_2026" /
    "05_Analysis" / "Analysis_Results"
)
output.mkdir(parents=True, exist_ok=True)

u = mda.Universe(str(topology), trajectories)

protein = u.select_atoms("protein and name CA")
ligand = u.select_atoms("resname MOL and not name H*")

assert len(u.trajectory) == 100
assert len(protein) == 315
assert len(ligand) == 30

# Reference frame
u.trajectory[0]

reference_protein = protein.positions.astype(float).copy()
reference_center = reference_protein.mean(axis=0)

def correct_ligand(protein_positions, ligand_positions, box):
    """Place ligand atoms together and near the protein using PBC."""
    center = protein_positions.mean(axis=0)

    # Reconstruct ligand around its first heavy atom
    anchor = ligand_positions[0].copy()
    internal_vectors = ligand_positions - anchor
    internal_vectors = minimize_vectors(internal_vectors, box)

    whole_ligand = anchor + internal_vectors

    # Place ligand image nearest to protein center
    ligand_center = whole_ligand.mean(axis=0)
    displacement = ligand_center - center
    corrected_displacement = minimize_vectors(
        displacement.reshape(1, 3), box
    )[0]

    return whole_ligand + (corrected_displacement - displacement)

reference_ligand = correct_ligand(
    protein.positions.astype(float),
    ligand.positions.astype(float),
    u.trajectory.ts.dimensions
)

rmsd_values = []

for ts in u.trajectory:
    mobile_protein = protein.positions.astype(float)
    mobile_ligand = correct_ligand(
        mobile_protein,
        ligand.positions.astype(float),
        ts.dimensions
    )

    mobile_center = mobile_protein.mean(axis=0)

    rotation, _ = align.rotation_matrix(
        mobile_protein - mobile_center,
        reference_protein - reference_center
    )

    aligned_ligand = (
        (mobile_ligand - mobile_center) @ rotation.T
        + reference_center
    )

    displacement = aligned_ligand - reference_ligand

    rmsd_values.append(
        np.sqrt(np.mean(np.sum(displacement ** 2, axis=1)))
    )

results = pd.DataFrame({
    "Time_ps": np.arange(1, 101),
    "BetaSitosterol_PBC_Corrected_RMSD_A": rmsd_values
})

results.to_csv(
    output / "PfLDH_BetaSitosterol_PBC_Corrected_RMSD.csv",
    index=False
)

plt.figure(figsize=(8, 5))
plt.plot(results["Time_ps"],
         results["BetaSitosterol_PBC_Corrected_RMSD_A"])
plt.xlabel("Unrestrained NPT time (ps)")
plt.ylabel("Beta-sitosterol RMSD (Angstrom)")
plt.title("PBC-corrected, protein-aligned beta-sitosterol RMSD")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(
    output / "PfLDH_BetaSitosterol_PBC_Corrected_RMSD.png",
    dpi=300
)
plt.close()

print("PBC-CORRECTED RMSD: PASS")
print("Mean RMSD (A):", round(np.mean(rmsd_values), 4))
print("Maximum RMSD (A):", round(np.max(rmsd_values), 4))
print("Results saved:", output)
