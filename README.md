# Mathematical Optimisation of Spatial Selectivity in Deep Brain Stimulation under Anisotropic Brain Tissue Conductivity

## Project Overview

This project investigates how anisotropic brain tissue conductivity affects the spatial selectivity of Deep Brain Stimulation (DBS), and how stimulation parameters can be optimised to maximise activation of a desired target region while limiting stimulation of surrounding tissue.

The model compares isotropic and anisotropic conductivity assumptions using both two-dimensional finite-difference modelling and three-dimensional finite-element modelling.

The code also performs an optimisation over electrode contact and stimulation voltage to determine DBS settings that provide high target coverage while minimising stimulation of an avoidance region.

---

## Model

The DBS electrode is based on the geometry of a Medtronic 3389 lead with four contacts.

The model includes:

- Four electrode contacts
- An encapsulation layer surrounding the DBS lead
- Isotropic and anisotropic brain tissue conductivity
- Different anisotropy ratios and fibre orientations
- An electric-field activation threshold of 0.2 V/mm

The Volume of Tissue Activated (VTA) is defined as the region where the magnitude of the electric field satisfies

|E| >= 0.2 V/mm.

---

## Code Structure

The complete project is executed using:

`run_all.m`

This runs the four main parts of the project in sequence.

### Part 1 – 2D finite-difference model

`part1_2d5.m`

This section develops and verifies a two-dimensional finite-difference model of the electric potential surrounding the DBS electrode.

It performs:

- a manufactured-solution convergence test
- comparison with an analytical isotropic solution
- VTA area convergence analysis
- Richardson extrapolation
- trapezium-rule convergence analysis
- investigation of truncation and round-off error
- condition-number analysis
- domain-size sensitivity testing
- comparison between isotropic and anisotropic conductivity

The resulting electric fields and VTA shapes are used to demonstrate how tissue anisotropy changes the spatial distribution of stimulation.

### Part 2 – 3D finite-element verification

`part2_3d_verify14.m`

This section verifies the three-dimensional finite-element solver.

Verification includes:

- manufactured solutions for isotropic and anisotropic conductivity
- comparison with an analytical coaxial-cylinder solution
- mesh-convergence analysis
- Richardson extrapolation
- domain-size sensitivity analysis

These tests are used to determine whether the numerical 3D model produces stable and convergent results.

### Part 3 – Isotropic versus anisotropic tissue

`part3_3d_sweep13.m`

This section investigates how anisotropy changes the Volume of Tissue Activated.

For a selected DBS contact, the model varies:

- anisotropy ratio κ = 1, 2, 4.4, 7 and 10
- fibre orientation θ = 0°, 45° and 90°

For each combination, the code calculates:

- VTA volume
- VTA spread in the x, y and z directions
- sphericity
- Dice coefficient relative to the isotropic VTA
- relative electric-field difference

The code also determines the anisotropy levels at which selected changes in VTA volume and spatial overlap occur.

### Part 4 – Optimisation of DBS settings

`part4_optimise15.m`

The final section performs an exhaustive optimisation over DBS electrode contact and stimulation voltage.

The stimulation voltage is varied from 0.5 V to 5 V in increments of 0.05 V for each of the four electrode contacts.

A target region representing the desired stimulation region and an avoidance region representing tissue that should not be stimulated are defined.

For every possible setting, the code calculates:

- target coverage
- spill into the avoidance region
- total VTA volume
- proportion of the VTA outside the target

The optimisation objective is

J = coverage - λ × spill

where λ controls the penalty assigned to stimulation of the avoidance region.

The code determines the optimal electrode contact and stimulation voltage for several isotropic and anisotropic tissue models.

It also investigates model mismatch by applying settings optimised under an isotropic conductivity assumption to anisotropic tissue.

---

## Running the Project

The project is configured to run directly through GitHub Actions using GNU Octave.

To execute the complete project:

1. Open the **Actions** tab of this repository.
2. Select **Execute Project**.
3. Click **Run workflow**.
4. Open the workflow run to view the execution output.

The workflow executes:

`run_all.m`

which automatically runs all four parts of the project.

No modification of file paths or source code is required.

---

## Output

During execution, the code creates the `Finalresults` directory.

The following numerical result files are generated:

- `results2d.mat`
- `verify3d.mat`
- `sweep3d.mat`
- `opt3d.mat`

The following text summaries are generated:

- `summary2d.txt`
- `summary3d_verify.txt`
- `summary3d_sweep.txt`
- `summary3d_opt.txt`

The project also generates six figures:

- `fig1_geometry.png` – DBS lead, target region and avoidance-region geometry
- `fig2_verification.png` – numerical convergence and solver verification
- `fig3_2d_fields.png` – 2D isotropic and anisotropic electric-field comparison
- `fig4_crosssections.png` – 3D VTA cross-sections for different fibre orientations
- `fig5_sensitivity.png` – sensitivity of the VTA to anisotropy ratio and fibre orientation
- `fig6_optimisation.png` – target coverage, spill and optimisation results

When the project is executed through GitHub Actions, the generated `Finalresults` directory is uploaded as a workflow artifact so that the results can be downloaded after execution.

---

## Supporting Functions

The remaining MATLAB files contain the numerical and modelling functions used by the main scripts.

These include functions for:

- defining DBS model parameters
- constructing the 3D mesh
- constructing conductivity tensors
- assembling the finite-element stiffness matrix
- solving the finite-element system
- calculating electric fields
- calculating VTA area and volume
- calculating VTA shape metrics
- numerical bisection

---

## Software

The original code was written in MATLAB syntax.

The GitHub workflow executes the code using GNU Octave to allow the full project to be tested directly on GitHub without requiring manual setup by the user.
