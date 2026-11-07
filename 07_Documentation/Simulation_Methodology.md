# PfLDH-NAD+-Beta-Sitosterol Molecular Dynamics Methodology

## Simulation system
The solvated PfLDH-NAD+-beta-sitosterol system contained 64,071 atoms, including 19,653 water molecules, 54 sodium ions, and 53 chloride ions. The protein contained 315 C-alpha atoms. NAD+ was designated NDP and beta-sitosterol was designated MOL.

## Simulation conditions
Molecular dynamics was performed using OpenMM 8.6 with an NVIDIA RTX A6000 GPU, CUDA mixed precision, particle mesh Ewald electrostatics, a 1 nm nonbonded cutoff, hydrogen-bond constraints, and rigid water molecules.

A LangevinMiddleIntegrator was used at 300 K, with a friction coefficient of 1 ps^-1 and a timestep of 1 fs. Pressure was controlled at 1 bar using a Monte Carlo barostat with a frequency of 25 steps. The recorded random seed was 20261106.

## Production trajectory
The unrestrained NPT simulation started from a saved state following 112 ps of restrained equilibration. The recorded production segment covered 100 ps, comprising 100,000 integration steps organized into 20 segments of 5 ps each.

Trajectory frames were saved every 1,000 steps, corresponding to 1 ps intervals and 100 recorded frames. Twenty XML checkpoint states were retained.

## Trajectory analysis
The analyses included protein C-alpha RMSD, NAD+ RMSD, beta-sitosterol RMSD, periodic-boundary-condition quality control, protein-ligand contacts, simulation box volume, and system density.

## Limitations
The 100 ps trajectory provides short-timescale structural observations and cannot independently establish equilibration, long-term complex stability, binding affinity, or enzyme inhibition.

The PBC-corrected beta-sitosterol RMSD is provisional because the reconstruction procedure has not been independently validated for bonded geometry and image selection.

The original system preparation and parameterization details should be consulted before making claims about force fields, water models, or cofactor parameters.
