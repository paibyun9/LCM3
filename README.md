# LCM3 – NMCC Spatial Analysis Dataset

## Overview
This repository provides the dataset and code used to reproduce the spatial clustering analysis presented in the study.

## Repository Structure

* data.csv  
  Dataset of 94 WENN archaeological sites

* nmcc_ripleys_k_analysis.py  
  Python script used for spatial analysis

## Quick Start (Reproducibility)

1. Download files  
   Download both files:

   * data.csv

   * nmcc_ripleys_k_analysis.py

2. Place files together  
   Put both files in the same folder  
   (e.g., Desktop)

3. Run the script  
   Open Terminal and run:

   ```bash
   python3 nmcc_ripleys_k_analysis.py

Output
Results will be printed directly in the terminal.

Requirements

Python 3
pandas
numpy
scipy

Install if needed:
Bashpip install pandas numpy scipy
Data Description
The dataset contains 94 spatially validated archaeological sites with:

Geographic coordinates
Classification into two groups (K-A and K-B)

Only complete and valid records were included.
What the Script Does
The script performs:

Distance calculation between sites
Spatial clustering analysis
Monte Carlo simulation (999 iterations)

All parameters are fixed for reproducibility.
Reproducibility Notes

Running the script reproduces the reported results
Minor numerical differences may occur due to system precision
Overall statistical patterns remain unchanged

Additional Information

