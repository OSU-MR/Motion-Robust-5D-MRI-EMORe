# EMORe: Motion-Robust 5D MRI Reconstruction via Expectation-Maximization–Guided Binning Correction and Outlier Rejection

[![DOI](https://img.shields.io/badge/DOI-10.1016/j.jocmr.2024.101509-blue)](https://doi.org/10.1016/j.jocmr.2024.101509)
[![arXiv](https://img.shields.io/badge/arXiv-24XX.XXXXX-b31b1b.svg)](https://arxiv.org/abs/24XX.XXXXX)
[![Zenodo](https://img.shields.io/badge/DOI-10.5281/zenodo.15769665-blue)](https://zenodo.org/records/15769665)

Official MATLAB implementation for the paper **EMORe: Motion-Robust 5D MRI Reconstruction via Expectation-Maximization–Guided Binning Correction and Outlier Rejection**.

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
If you use this code or the EMORe framework in your research, please cite published peer-reviewed abstract (presented at SCMR'25):

```bibtex
@article{arshad2025emore,
  title={EMORe: motion-robust XD-CMR reconstruction using expectation-maximization (EM) algorithm},
  author={Arshad, Syed Murtaza and Potter, Lee C and Lei, Xuan and Ahmad, Rizwan},
  journal={Journal of Cardiovascular Magnetic Resonance},
  volume={27},
  year={2025},
  publisher={Elsevier}
}
```

### Repository Structure
The repository is organized as follows:

```
EMORe/
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
```

-   **`mainRecon.m`**: The main script to execute the entire reconstruction pipeline.
-   **`Datasets/`**: Folder to store input raw data files (`.mat`).
-   **`Recon Functions/`**: Contains the core reconstruction algorithms and all necessary helper functions.
-   **`Recons/`**: The default output directory where all results will be saved.

### System Requirements

#### Software:
-   MATLAB (tested on R2022b or later)
-   Parallel Computing Toolbox™ (for GPU acceleration)

#### Hardware:
-   An NVIDIA GPU with CUDA support is strongly recommended.

### Installation and Setup

1.  **Clone the repository:**
    ```bash
    git clone [https://github.com/your-username/EMORe-Recon.git](https://github.com/your-username/EMORe-Recon.git)
    ```
2.  **Download Sample Data:** Download the sample *in vivo* dataset from Zenodo:
    - **Link:** [https://zenodo.org/records/15769665](https://zenodo.org/records/15769665)
    - Download the `volunteer1.mat` file.

3.  **Prepare Data:** Place the downloaded `volunteer1.mat` file (or your own `.mat` files) into the `Datasets/` folder.

### How to Run

1.  Open MATLAB.
2.  Navigate to the repository's root directory.
3.  Run the main script from the MATLAB command window:
    ```matlab
    mainRecon
    ```
4.  A dialog box will appear. Select the `dataset_invivo.mat` file from the `Datasets/` folder.
5.  If you have multiple GPUs, a prompt will appear to select one.
6.  The script will execute the full pipeline, displaying progress in the command window.

### Configuration Parameters
All key reconstruction parameters can be adjusted in the `Reconstruction Parameters` section of `mainRecon.m`.

```matlab
% Comparison flag
opt.cs         = 1;                    % 1 = run CS recon; 0 = skip CS

% ADMM tuning
opt.lam        = 5e-4 .* [1 1 1 5 3];  % [λ_sx, λ_sy, λ_sz, λ_c, λ_r]
opt.gStp       = 4e-2;                 % gradient step size (γ)
opt.mu         = 5e-1;                 % Lagrange multiplier (μ)

% EMORe priors
opt.out_factor = 0.05;                  % α_o (outlier prior)
opt.sg_factor  = 0.85;                  % α_g (self‐gating prior)

% ... and others
```

### Expected Output
A unique, timestamped directory will be created inside `Recons/volunteer1/`, containing:

-   **`params/`**: `.mat` files with final images and a `.txt` log file.
-   **`gif/`**: GIFs of the reconstructed image series.
-   **`outliers/`**: Plots showing signal traces and outlier rejection.
