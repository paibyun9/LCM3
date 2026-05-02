# LCM3 – NMCC Spatial Analysis Dataset

## Overview
This repository provides the dataset and code used to reproduce the spatial clustering analysis presented in the study.

## Repository Structure

- data.csv  
  Dataset of 94 WENN archaeological sites  

- nmcc_ripleys_k_analysis.py  
  Python script used for spatial analysis  

## Quick Start (Reproducibility)

### 1. Download files

Download both files:

- data.csv  
- nmcc_ripleys_k_analysis.py  

### 2. Place files together

Put both files in the same folder (for example, Desktop)

### 3. Run the script

Open Terminal and run:

```bash
python3 nmcc_ripleys_k_analysis.py
4. Output

Results will be printed directly in the terminal.

Requirements
Python 3
pandas
numpy
scipy

Install if needed:

pip install pandas numpy scipy
Data Description

The dataset contains 94 archaeological sites with:

Geographic coordinates
Classification into two groups (K-A and K-B)

Only valid records were included.

Method Summary

The analysis script performs:

Distance calculation
Spatial clustering analysis
Monte Carlo simulation (999 iterations)

All parameters are fixed for reproducibility.

Notes
Running the script reproduces the reported results
Minor numerical differences may occur depending on system precision
Overall patterns remain unchanged
Implementation

## The full computational workflow is implemented in:

nmcc_ripleys_k_analysis.py

