from pathlib import Path
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.lib.distances import (
    distance_array,
    minimize_vectors
)

root = Path.home() / "Downloads" / "PfLDH_Unrestrained_NPT_100ps_Backup"

topology = root / "PfLDH_NAD_BetaSitosterol_Solvated.prmtop"

trajectories = [
    str(root / "Unrestrained_NPT" / f"NPT_segment_{i:02d}.dcd")
    for i in range(1, 21)
]

output = (
    Path.home() / "Documents" /
    "PfLDH_BetaSitosterol_MD_2026" /
    "05_Analysis" / "Analysis_Results"
)
output.mkdir(parents=True, exist_ok=True)

u = mda.Universe(str(topology), trajectories)

protein = u.select_atoms("protein and not name H*")
ligand = u.select_atoms("resname MOL and not name H*")

assert len(u.trajectory) == 100
assert len(ligand) == 30

records = []

for frame, ts in enumerate(u.trajectory, start=1):
    box = ts.dimensions

    if box is None or not np.all(np.isfinite(box)):
        raise ValueError(f"Invalid box at frame {frame}")

    # Reconstruct ligand around its first heavy atom
    anchor = ligand.positions[0].astype(float)
    vectors = ligand.positions.astype(float) - anchor
    vectors = minimize_vectors(vectors, box)
    reconstructed = anchor + vectors

    # Measure molecular continuity
    raw_span = np.max(
        np.linalg.norm(
            ligand.positions - ligand.positions[0],
            axis=1
        )
    )

    reconstructed_span = np.max(
        np.linalg.norm(reconstructed - reconstructed[0], axis=1)
    )

    # Compare protein-ligand minimum distances
    pbc_distance = distance_array(
        protein.positions,
        ligand.positions,
        box=box
    ).min()

    reconstructed_distance = distance_array(
        protein.positions,
        reconstructed
    ).min()

    records.append({
        "Time_ps": frame,
        "Raw_Ligand_Span_A": raw_span,
        "Reconstructed_Ligand_Span_A": reconstructed_span,
        "Protein_Ligand_MinDist_PBC_A": pbc_distance,
        "Protein_Ligand_MinDist_Reconstructed_A":
            reconstructed_distance
    })

df = pd.DataFrame(records)

df.to_csv(
    output / "PfLDH_PBC_Reconstruction_Validation.csv",
    index=False
)

print("PBC RECONSTRUCTION CHECK: PASS")
print("Frames:", len(df))

print(
    "Maximum raw ligand span (A):",
    round(df["Raw_Ligand_Span_A"].max(), 4)
)

print(
    "Maximum reconstructed ligand span (A):",
    round(df["Reconstructed_Ligand_Span_A"].max(), 4)
)

print(
    "Mean reconstructed ligand span (A):",
    round(df["Reconstructed_Ligand_Span_A"].mean(), 4)
)

print(
    "Maximum PBC protein-ligand minimum distance (A):",
    round(df["Protein_Ligand_MinDist_PBC_A"].max(), 4)
)

print(
    "Results saved:", output
)
