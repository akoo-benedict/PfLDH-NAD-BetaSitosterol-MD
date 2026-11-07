from pathlib import Path
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array

# Original simulation files
root = Path.home() / "Downloads" / "PfLDH_Unrestrained_NPT_100ps_Backup"

topology = root / "PfLDH_NAD_BetaSitosterol_Solvated.prmtop"

trajectories = [
    str(root / "Unrestrained_NPT" / f"NPT_segment_{i:02d}.dcd")
    for i in range(1, 21)
]

# Central results folder
output = (
    Path.home() / "Documents" /
    "PfLDH_BetaSitosterol_MD_2026" /
    "05_Analysis" / "Analysis_Results"
)
output.mkdir(parents=True, exist_ok=True)

u = mda.Universe(str(topology), trajectories)

protein = u.select_atoms("protein and not name H*")
ligand = u.select_atoms("resname MOL and not name H*")
nad = u.select_atoms("resname NDP and not name H*")

assert len(u.trajectory) == 100
assert len(ligand) == 30
assert len(nad) == 44

records = []

for frame, ts in enumerate(u.trajectory, start=1):

    box = np.asarray(ts.dimensions, dtype=float)

    valid_box = (
        box.shape == (6,)
        and np.all(np.isfinite(box))
        and np.all(box[:3] > 0)
    )

    if not valid_box:
        raise ValueError(
            f"Invalid periodic box at frame {frame}: {box}"
        )

    # Distances with and without periodic boundaries
    ligand_pbc = distance_array(
        protein.positions, ligand.positions, box=box
    ).min()

    ligand_raw = distance_array(
        protein.positions, ligand.positions
    ).min()

    nad_pbc = distance_array(
        protein.positions, nad.positions, box=box
    ).min()

    nad_raw = distance_array(
        protein.positions, nad.positions
    ).min()

    records.append({
        "Time_ps": frame,
        "Box_X_A": box[0],
        "Box_Y_A": box[1],
        "Box_Z_A": box[2],
        "Ligand_MinDist_PBC_A": ligand_pbc,
        "Ligand_MinDist_Raw_A": ligand_raw,
        "NAD_MinDist_PBC_A": nad_pbc,
        "NAD_MinDist_Raw_A": nad_raw
    })

df = pd.DataFrame(records)

df.to_csv(
    output / "PfLDH_PBC_Quality_Control.csv",
    index=False
)

print("PBC QUALITY CONTROL: PASS")
print("Frames analyzed:", len(df))

for molecule in ["Ligand", "NAD"]:
    difference = abs(
        df[f"{molecule}_MinDist_Raw_A"]
        - df[f"{molecule}_MinDist_PBC_A"]
    )

    print(
        molecule,
        "frames with distance difference > 0.5 A:",
        int((difference > 0.5).sum())
    )

    print(
        molecule,
        "maximum distance difference (A):",
        round(difference.max(), 4)
    )

print("Results saved:", output)
