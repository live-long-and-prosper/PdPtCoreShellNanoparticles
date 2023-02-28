# PdPtCoreShellNanoparticles
This repository contains the code and data for the following paper:

**Probing the atomically diffuse interfaces in core-shell nanoparticles in three dimensions**

Zezhou Li<sup>1</sup>, Zhiheng Xie<sup>1</sup>, Yao Zhang<sup>1</sup>, Xilong Mu<sup>1</sup>, Jisheng Xie<sup>1</sup>,  Haijing Yin<sup>1</sup>, Ya-wen Zhang<sup>1</sup>, Colin Ophus<sup>2</sup>, Jihan Zhou<sup>1†</sup>

*<sup>1</sup>Beijing National Laboratory for Molecular Sciences, College of Chemistry and Molecular Engineering, Peking University, Beijing, 100871, China.*

*<sup>2</sup>National Center for Electron Microscopy, Molecular Foundry, Lawrence Berkeley National Laboratory, Berkeley, CA 94720, USA.*

*†Correspondence and requests for materials should be addressed to J. Z. (email: jhzhou@pku.edu.cn).*

# Repositary Contents

### 1. Experiment Data

Folder: [Measured_data](./1_Measured_data)

This folder contains denoised and aligned ADF-STEM projections and corresponding finalized tilt angles for three Pd@Pt core-shell nanoparticles. Three particles are named PB (pentagonal bipyramid shaped), EPB (elongated pentagonal bipyramid shaped) and TO (truncated octahedron shaped), respectively.

### 2. Reconstructed 3D Volume

Folder: [Final_reconstruction_volume](./2_Final_reconstruction_volume)

This folder contains 3D tomographic reconstruction volumes of three particles. For the source code of RESIRE algorithm used in these reconstructions, please see the [source code](https://github.com/AET-MetallicGlass/Supplementary-Data-Codes/tree/master/2_RESIRE_package) of Yao Yang's paper on github.

### 3. Atom Tracing and Classification

Folder: [Tracing_and_classification](./3_Tracing_and_classification)

This folder contains the source code to trace and classify atoms in the 3D volume.

### 4. Experimental Atomic Models

Folder: [Final_coordinates](./4_Final_coordinates)

This folder contains the final coordinates of three nanoparticles.

### 5. Analysis of core-shell interface and others

Folder: [Analysis_of_interface](./5_Analysis_of_interface)

This folder contains the codes to analyse the pair distribution function (PDF), the core-shell interface, the local coordination structure (PTM), the displacement and strain map of three nanoparticles.
