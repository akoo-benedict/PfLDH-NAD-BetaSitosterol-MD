from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.lib.mdamath import triclinic_vectors

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

u = mda.Universe(str(topology), trajectories)

assert len(u.trajectory) == 100
assert len(u.atoms) == 64071

# Total system mass from the topology, in daltons
total_mass_da = float(u.atoms.masses.sum())

assert np.isfinite(total_mass_da) and total_mass_da > 0

records = []

for frame, ts in enumerate(u.trajectory, start=1):
    box = np.asarray(ts.dimensions, dtype=float)

    if (
        box.shape != (6,)
        or not np.all(np.isfinite(box))
        or np.any(box[:3] <= 0)
    ):
        raise ValueError(f"Invalid box at frame {frame}: {box}")

    # Supports orthorhombic and triclinic boxes
    vectors = triclinic_vectors(box).astype(float)
    volume_a3 = abs(np.linalg.det(vectors))

    if volume_a3 <= 0:
        raise ValueError(f"Invalid volume at frame {frame}")

    # 1 dalton = 1.66053906660e-24 grams
    # 1 cubic angstrom = 1e-24 cubic centimetres
    density_g_cm3 = total_mass_da * 1.66053906660 / volume_a3

    records.append({
        "Time_ps": frame,
        "Volume_A3": volume_a3,
        "Density_g_cm3": density_g_cm3
    })

df = pd.DataFrame(records)

# Save numerical data
df.to_csv(
    output / "PfLDH_NPT_Volume_Density.csv",
    index=False
)

# Volume figure
plt.figure(figsize=(8, 5))
plt.plot(df["Time_ps"], df["Volume_A3"])
plt.xlabel("Unrestrained NPT time (ps)")
plt.ylabel("Simulation box volume (cubic Angstrom)")
plt.title("PfLDH system volume during 100 ps NPT")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(output / "PfLDH_NPT_Volume.png", dpi=300)
plt.close()

# Density figure
plt.figure(figsize=(8, 5))
plt.plot(df["Time_ps"], df["Density_g_cm3"])
plt.xlabel("Unrestrained NPT time (ps)")
plt.ylabel("System density (g/cm3)")
plt.title("PfLDH system density during 100 ps NPT")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(output / "PfLDH_NPT_Density.png", dpi=300)
plt.close()

print("NPT VOLUME AND DENSITY: PASS")
print("Total system mass (Da):", round(total_mass_da, 2))
print("Mean volume (A3):", round(df["Volume_A3"].mean(), 2))
print("Mean density (g/cm3):",
      round(df["Density_g_cm3"].mean(), 4))
print("Minimum density (g/cm3):",
      round(df["Density_g_cm3"].min(), 4))
print("Maximum density (g/cm3):",
      round(df["Density_g_cm3"].max(), 4))
print("Results saved:", output)
