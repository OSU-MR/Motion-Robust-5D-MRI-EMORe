# EMORe: Motion-Robust 5D MRI Reconstruction via Expectation-Maximization–Guided Binning Correction and Outlier Rejection

[![DOI](https://img.shields.io/badge/DOI-10.1016/j.jocmr.2024.101509-blue)](https://doi.org/10.1016/j.jocmr.2024.101509)
[![arXiv](https://img.shields.io/badge/arXiv-24XX.XXXXX-b31b1b.svg)](https://arxiv.org/abs/24XX.XXXXX) Official MATLAB implementation for the paper **EMORe: Motion-Robust 5D MRI Reconstruction via Expectation-Maximization–Guided Binning Correction and Outlier Rejection**.

**Authors**: Syed M. Arshad, Lee C. Potter, Chong Chen, Preethi C. Subramanian, Yingmin Liu, Christopher Crabtree, Xuan Lei, Juliet Varghese, Matthew S. Tong, and Rizwan Ahmad

**Contact**: [arshad.32@osu.edu](mailto:arshad.32@osu.edu)

---

### Table of Contents
1.  [Abstract](#abstract)
2.  [Citation](#citation)
3.  [Repository Structure](#repository-structure)
4.  [System Requirements](#system-requirements)
5.  [Installation and Setup](#installation-and-setup)
6.  [How to Run](#how-to-run)
7.  [Configuration Parameters](#configuration-parameters)
8.  [Expected Output](#expected-output)
9.  [License](#license)

### Abstract

We propose EMORe, an adaptive reconstruction method designed to enhance motion robustness in free-running, free-breathing self-gated 5D cardiac magnetic resonance imaging (MRI). Traditional self-gating-based motion binning for 5D MRI often results in residual motion artifacts due to inaccuracies in cardiac and respiratory signal extraction and sporadic bulk motion, compromising clinical utility. EMORe addresses these issues by integrating adaptive inter-bin correction and explicit outlier rejection within an expectation-maximization (EM) framework, whereby the E-step and M-step are executed alternately until convergence. In the E-step, probabilistic (soft) bin assignments are refined by correcting misassignment of valid data and rejecting motion-corrupted data to a dedicated outlier bin. In the M-step, the image estimate is improved using the refined soft bin assignments. Validation in a simulated 5D MRXCAT phantom and in vivo volunteers confirmed EMORe's robustness, significantly enhancing image quality and reducing motion artifacts compared to standard compressed sensing. Although EMORe incurs a modest increase in computational complexity, its adaptability and robust handling of bulk motion artifacts significantly enhance the clinical applicability and diagnostic confidence of 5D cardiac MRI.

### Citation
If you use this code or the EMORe framework in your research, please cite our paper:

```bibtex
@article{Arshad2025_EMORe,
  author    = {Arshad, Syed M. and Potter, Lee C. and Chen, Chong and Subramanian, Preethi C. and Liu, Yingmin and Crabtree, Christopher and Lei, Xuan and Varghese, Juliet and Tong, Matthew S. and Ahmad, Rizwan},
  title     = {{EMORe}: Motion-Robust {5D} {MRI} Reconstruction via Expectation-Maximization–Guided Binning Correction and Outlier Rejection},
  journal   = {},
  year      = {2025},
  volume    = {},
  number    = {},
  pages     = {},
  doi       = {},
  publisher = {}
}
Repository Structure
The repository is organized as follows:

Github Repo For EMORe/
│
├── Datasets/
│   └── (Place your raw .mat data files here)
│
├── Recon Functions/
│   ├── runEMOReRecon.m
│   ├── runCSRecon.m
│   └── (Other helper functions for SG, binning, ADMM, etc.)
│
├── Recons/
│   └── (Output directory, created automatically)
│
└── mainRecon.m
mainRecon.m: The main script to execute the entire reconstruction pipeline.
Datasets/: Folder to store input raw data files (.mat).
Recon Functions/: Contains the core reconstruction algorithms and all necessary helper functions.
Recons/: The default output directory where all results will be saved.
System Requirements
Software:
MATLAB (tested on R2022b or later)
Parallel Computing Toolbox™ (for GPU acceleration)
Hardware:
An NVIDIA GPU with CUDA support is strongly recommended.
Installation and Setup
Clone the repository:
Bash

git clone [https://github.com/your-username/EMORe-Recon.git](https://github.com/your-username/EMORe-Recon.git)
Prepare Data: Place your raw .mat k-space data files into the Datasets/ folder.
How to Run
Open MATLAB.
Navigate to the repository's root directory.
Run the main script:
Matlab

mainRecon
A dialog box will appear. Select the raw data file from the Datasets/ folder.
If you have multiple GPUs, a prompt will appear to select one.
The script will execute the full pipeline, displaying progress in the command window.
Configuration Parameters
All key reconstruction parameters can be adjusted in the Reconstruction Parameters section of mainRecon.m.

Matlab

% Comparison flag
opt.cs         = 1;                   % 1 = run CS recon; 0 = skip CS
% ADMM tuning
opt.lam        = 5e-4 .* [1 1 1 5 3];  % [λ_sx, λ_sy, λ_sz, λ_c, λ_r]
opt.gStp       = 4e-2;                % gradient step size (γ)
opt.mu         = 5e-1;                % Lagrange multiplier (μ)
% EMORe priors
opt.out_factor = 0.05;                % α_o (outlier prior)
opt.sg_factor  = 0.85;                % α_g (self‐gating prior)
% ... and others
Expected Output
A unique, timestamped directory will be created inside Recons/<dataset_name>/, containing:

params/: .mat files with final images and a .txt log file.
gif/: GIFs of the reconstructed image series.
outliers/: Plots showing signal traces and outlier rejection.
