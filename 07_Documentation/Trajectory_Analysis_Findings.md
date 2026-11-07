# PfLDH-NAD+-Beta-Sitosterol Trajectory Analysis Findings

## Analysis scope
These findings describe a 100 ps unrestrained NPT molecular dynamics trajectory containing 100 recorded frames.

## Protein and cofactor RMSD
- Protein C-alpha RMSD: mean 0.7888 Angstrom; maximum 1.0153 Angstrom.
- NAD+ protein-aligned RMSD: mean 0.9219 Angstrom; maximum 1.5606 Angstrom.

## Beta-sitosterol RMSD
- Original protein-aligned RMSD: mean 4.6171 Angstrom; maximum 7.1300 Angstrom.
- PBC-corrected RMSD: mean 3.0501 Angstrom; maximum 4.3524 Angstrom.

The PBC-corrected results remain provisional pending independent validation of the reconstruction procedure.

## Protein-ligand contacts
The mean number of protein residues within 4 Angstrom of beta-sitosterol was 6.07, with a range of 3 to 9.

Residue contact occupancy:
- Ile227: 94%
- Tyr160: 91%
- Ala237: 80%
- Tyr235: 76%
- Val236: 75%
- Pro234: 55%
- Ser156: 48%
- His231: 42%
- Pro170: 16%
- Tyr161: 15%

## Periodic-boundary quality control
- Beta-sitosterol: 100 of 100 frames showed a PBC versus non-PBC minimum-distance difference greater than 0.5 Angstrom.
- Maximum difference: 43.7662 Angstrom.
- NAD+: 0 of 100 frames exceeded the same threshold.
- Maximum reconstructed ligand span: 14.7245 Angstrom.
- Mean reconstructed ligand span: 12.7158 Angstrom.
- Maximum PBC protein-ligand minimum distance: 3.817 Angstrom.

## System volume and density
- Mean simulation box volume: 642671.88 cubic Angstrom.
- Mean density: 1.0135 g/cm^3.
- Density range: 0.9875 to 1.0181 g/cm^3.

## Interpretation and limitations
These findings represent short-timescale structural observations, not proof of equilibrium or long-term stability.

Contact occupancy measures the frequency of geometric proximity within the analyzed trajectory and does not independently establish binding affinity or enzyme inhibition.

The corrected beta-sitosterol RMSD values are provisional because bonded geometry and periodic-image selection have not been independently validated.

The numerical CSV files and original PNG figures are preserved in the corresponding GitHub repository folders.
