from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array

# Original simulation files
root = Path.home() / "Downloads" / "PfLDH_Unrestrained_NPT_100ps_Backup"

topology = root / "PfLDH_NAD_BetaSitosterol_Solvated.prmtop"

trajectories = [
    str(root / "Unrestrained_NPT" / f"NPT_segment_{i:02d}.dcd")
    for i in range(1, 21)
]

# Central analysis results folder
output = (
    Path.home() / "Documents" /
    "PfLDH_BetaSitosterol_MD_2026" /
    "05_Analysis" / "Analysis_Results"
)
output.mkdir(parents=True, exist_ok=True)

# Load the trajectory
u = mda.Universe(str(topology), trajectories)

protein = u.select_atoms("protein and not name H*")
ligand = u.select_atoms("resname MOL and not name H*")

assert len(u.trajectory) == 100
assert len(protein) > 0
assert len(ligand) == 30

cutoff = 4.0  # Angstrom

frame_results = []
residue_contacts = {}

# Analyze all trajectory frames
for frame_number, ts in enumerate(u.trajectory, start=1):

    distances = distance_array(
        protein.positions,
        ligand.positions,
        box=ts.dimensions
    )

    contacting_atoms = np.any(distances <= cutoff, axis=1)
    contacting_residues = protein[contacting_atoms].residues

    frame_results.append({
        "Time_ps": frame_number,
        "Contacting_Residues": len(contacting_residues),
        "Atom_Pair_Contacts": int(np.sum(distances <= cutoff))
    })

    for residue in contacting_residues:
        key = f"{residue.resname}{residue.resid}"
        residue_contacts[key] = residue_contacts.get(key, 0) + 1

# Save frame-by-frame results
df = pd.DataFrame(frame_results)

df.to_csv(
    output / "PfLDH_BetaSitosterol_Contacts.csv",
    index=False
)

# Save residue contact frequencies
frequency = pd.DataFrame([
    {
        "Residue": residue,
        "Frames_in_Contact": count,
        "Occupancy_Percent": 100 * count / len(u.trajectory)
    }
    for residue, count in residue_contacts.items()
])

frequency = frequency.sort_values(
    "Occupancy_Percent", ascending=False
)

frequency.to_csv(
    output / "PfLDH_BetaSitosterol_Residue_Occupancy.csv",
    index=False
)

# Generate contact plot
plt.figure(figsize=(8, 5))
plt.plot(df["Time_ps"], df["Contacting_Residues"])

plt.xlabel("Unrestrained NPT time (ps)")
plt.ylabel("PfLDH residues within 4 Angstrom")
plt.title("PfLDH-beta-sitosterol contacts during 100 ps NPT")
plt.grid(alpha=0.3)
plt.tight_layout()

plt.savefig(
    output / "PfLDH_BetaSitosterol_Contacts.png",
    dpi=300
)
plt.close()

print("CONTACT ANALYSIS: PASS")
print("Mean contacting residues:",
      round(df["Contacting_Residues"].mean(), 2))
print("Minimum contacting residues:",
      df["Contacting_Residues"].min())
print("Maximum contacting residues:",
      df["Contacting_Residues"].max())
print("Top contacting residues:")
print(frequency.head(10).to_string(index=False))
print("Results saved:", output)
