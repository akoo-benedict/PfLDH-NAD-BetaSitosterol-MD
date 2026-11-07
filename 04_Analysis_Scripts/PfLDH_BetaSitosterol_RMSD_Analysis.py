from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import align

# Original simulation files
root = Path.home() / "Downloads" / "PfLDH_Unrestrained_NPT_100ps_Backup"

topology = root / "PfLDH_NAD_BetaSitosterol_Solvated.prmtop"

trajectories = [
    str(root / "Unrestrained_NPT" / f"NPT_segment_{i:02d}.dcd")
    for i in range(1, 21)
]

# Central project results folder
output = (
    Path.home() / "Documents" /
    "PfLDH_BetaSitosterol_MD_2026" /
    "05_Analysis" / "Analysis_Results"
)
output.mkdir(parents=True, exist_ok=True)

# Load trajectory
u = mda.Universe(str(topology), trajectories)

protein = u.select_atoms("protein and name CA")
ligand = u.select_atoms("resname MOL and not name H*")

assert len(u.trajectory) == 100
assert len(protein) == 315
assert len(ligand) > 0

print("Protein CA atoms:", len(protein))
print("Ligand heavy atoms:", len(ligand))

# Reference coordinates from first trajectory frame
u.trajectory[0]
reference_protein = protein.positions.copy()
reference_ligand = ligand.positions.copy()

# Calculate ligand RMSD after protein alignment
rmsd_values = []

for ts in u.trajectory:
    mobile_protein = protein.positions.copy()
    mobile_ligand = ligand.positions.copy()

    mobile_center = mobile_protein.mean(axis=0)
    reference_center = reference_protein.mean(axis=0)

    rotation, _ = align.rotation_matrix(
        mobile_protein - mobile_center,
        reference_protein - reference_center
    )

    aligned_ligand = (
        (mobile_ligand - mobile_center) @ rotation.T
        + reference_center
    )

    displacement = aligned_ligand - reference_ligand

    ligand_rmsd = np.sqrt(
        np.mean(np.sum(displacement ** 2, axis=1))
    )

    rmsd_values.append(ligand_rmsd)

# Save numerical results
results = pd.DataFrame({
    "Time_ps": np.arange(1, len(rmsd_values) + 1),
    "BetaSitosterol_RMSD_A": rmsd_values
})

results.to_csv(
    output / "PfLDH_BetaSitosterol_RMSD.csv",
    index=False
)

# Save figure
plt.figure(figsize=(8, 5))
plt.plot(
    results["Time_ps"],
    results["BetaSitosterol_RMSD_A"]
)

plt.xlabel("Unrestrained NPT time (ps)")
plt.ylabel("Beta-sitosterol RMSD (Angstrom)")
plt.title("Protein-aligned beta-sitosterol RMSD")
plt.grid(alpha=0.3)
plt.tight_layout()

plt.savefig(
    output / "PfLDH_BetaSitosterol_RMSD.png",
    dpi=300
)
plt.close()

print("BETA-SITOSTEROL RMSD: PASS")
print("Mean RMSD (A):", round(np.mean(rmsd_values), 4))
print("Maximum RMSD (A):", round(np.max(rmsd_values), 4))
print("Results saved:", output)
