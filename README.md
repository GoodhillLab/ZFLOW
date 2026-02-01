# Zebrafish Fluid-dynamics Locomotion Optimization Workspace (ZFLOW)

C++ code for running Computational fluid dynamics (CFD) simulation of larval zebrafish, and MATLAB code for preprocessing and postprocessing.
See the companion

- Publication [Larval zebrafish minimize energy consumption during hunting via adaptive movement selection](https://doi.org/10.1073/pnas.2513853123) by Darveniza, T., Wong, R., Zhu, S., Pujic, Z., Sun, B., Levendosky, M., Agarwal, R., McCullough, M. H., & Goodhill, G. (2026). Proceedings of the National Academy of Sciences of the United States of America.
- Preprint [Behavioral adaptation to changing energy constraints via altered frequency of movement selection](https://www.biorxiv.org/content/10.1101/2023.11.08.566262v1) by by Darveniza, T., Zhu, S., Pujic, Z., Sun, B., Levendosky, M., Wong, R., Agarwal, R., McCullough, M. H., & Goodhill, G. for more information.

## Overview

To facilitate use of this CFD software by others, we demonstrate how to run a randomly-selected 9 dpf larval zebrafish swim bout through the CFD simulation pipeline.

## Preprocessing
Preprocess the swim bout according to steps described in the paper.
Following successful preprocessing, the swim bout should be in a similar form to [`simulation/zebrafish_2D_coords.txt`](simulation/zebrafish_2D_coords.txt). To better understand this form, run the script [`preprocess/preSimVis.m`](preprocess/preSimVis.m).

## Simulation

See [simulation/README.md](simulation/README.md)

## Postprocessing
1. Following successful simulation, there should be two new folders: `ZFISH3d`, and `examplespline.txt`.
	+ `ZFISH3d` contains various IBAMR outputs, including power expenditure
	+ `examplespline.txt` contains tracking points on the 3D fish during simulation, which is used during validation
2. Navigate to [`postprocess/postprocess.m`](postprocess/postprocess.m), and run the post-processing functions

## Requirements

This code requires MATLAB (tested on R2021b&ndash;R2023b) for the pre/postprocessing, `make` or `cmake` and [IBAMR](https://ibamr.github.io) library (tested on 0.{8,9,13,14,15}.0, currently 0.15.0 is used) for the main CFD simulation.

## Contact

- Robert Wong robert.wong@wustl.edu
- Thomas Darveniza darveniza.t@wustl.edu
