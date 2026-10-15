# Biomathematics Project

This repository contains the code used for the numerical modelling project.

## Running the project

The complete project is executed using:

run_all.m

The repository has been configured to execute automatically using GitHub Actions and GNU Octave.

To execute the project on GitHub:

1. Open the **Actions** tab.
2. Select **Execute Project**.
3. Click **Run workflow**.
4. Open the completed workflow to view the output.

The workflow automatically installs the required execution environment and runs:

octave --no-gui --quiet run_all.m

No modification of file paths or source code is required.

## Project structure

- `run_all.m` — runs the complete project
- `part1_2d5.m` — Part 1
- `part2_3d_verify14.m` — Part 2
- `part3_3d_sweep13.m` — Part 3
- `part4_optimise15.m` — Part 4
- remaining `.m` files contain supporting functions

## Results

Generated figures and numerical results are saved in the `Finalresults` directory.
