LCM3 – NMCC Spatial Analysis Dataset

This repository provides the dataset and code used for spatial clustering analysis.

Contents
data.csv
Dataset of 94 spatially validated WENN sites
nmcc_ripleys_k_analysis.py
Python script for spatial clustering analysis
Reproducibility (Simple Guide)

You can reproduce the analysis in a few simple steps:

Step 1 — Download files

Download the following files from this repository:

data.csv
nmcc_ripleys_k_analysis.py
Step 2 — Place files

Put both files in the same location (for example, your Desktop)

Step 3 — Run

Open Terminal and run:

python3 nmcc_ripleys_k_analysis.py
Step 4 — Done

The script will print the results directly in the terminal.

Requirements
Python 3
pandas
numpy
scipy

If needed, install packages using:

pip install pandas numpy scipy
Data Description

The dataset (data.csv) contains 94 archaeological sites with:

Geographic coordinates
Classification into two groups (K-A and K-B)

Only valid and complete records were included in the analysis.

Method Overview

The script performs:

Distance calculation between sites
Spatial clustering analysis
Monte Carlo simulation (999 runs) for statistical evaluation

All parameters are fixed to ensure reproducible results.

Notes
Results correspond to those reported in the study
Minor numerical differences may occur due to system precision
Overall patterns and statistical outcomes remain consistent
Reference

Detailed computational workflow is described in
Appendix 2 of the study
