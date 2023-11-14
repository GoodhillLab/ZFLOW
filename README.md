# Zebrafish Fluid-dynamics Locomotion Optimization Workspace (ZFLOW)

C++ code for running CFD simulation of larval zebrafish, and MATLAB code for preprocessing and postprocessing.
See the companion preprint [**Behavioral adaptation to changing energy constraints via altered frequency of movement selection**](https://www.biorxiv.org/content/10.1101/2023.11.08.566262v1) by Thomas Darveniza, Shuyu I Zhu, Zac Pujic, Biao Sun, Matthew Levendosky, Robert Wong, Ramesh Agarwal, Michael H McCullough, Geoffrey J Goodhill for more information.

## Overview

To facilitate use of this CFD software by others, we demonstrate how to run a randomly-selected 9 dpf larval zebrafish swim bout through the CFD simulation pipeline.

## Preprocessing
Preprocess the swim bout according to steps described in the paper.
Following successful preprocessing, the swim bout should be in a similar form to /simulation/zebrafish_2D_coords.txt. (To better understand this form, run the script preprocess/preSimVis.m)

## Simulation
1. Download and install [IBAMR](https://ibamr.github.io) on your local computer. Alternatively, [docker containers](https://hub.docker.com/repository/docker/d0ckaaa/ibamr) have also been made available by other IBAMR users.
2. Navigate to `simulation` folder.
3. Adjust the `Makefile` for your specific system.
4. Make and run the executable (run-time of this example on an Intel Cascade Lake core is approximately 1-2 hours).

## Postprocessing
1. Following successful simulation, there should be two new folders: `ZFISH3d`, and `examplespline.txt`.
	+ `ZFISH3d` contains various IBAMR outputs, including power expenditure
	+ `examplespline.txt` contains tracking points on the 3D fish during simulation, which is used during validation
2. Navigate to `postprocess/postprocess.m`, and run the post-processing functions

## Requirements

This code requires MATLAB (tested on R2021b --- R2023b) for the pre/postprocessing, `make` or `cmake` and [IBAMR](https://ibamr.github.io) library for the main CFD simulation.

## Contact
Thomas Darveniza darveniza.t@wustl.edu
