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
  journal   = {IEEE Transactions on Medical Imaging},
  year      = {2025},
  volume    = {},
  number    = {},
  pages     = {},
  doi       = {},
  publisher = {IEEE}
}
