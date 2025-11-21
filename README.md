# VNA Validation of Saline Dielectric Properties

This repository contains MATLAB scripts used to compare measured dielectric properties of saline to known reference values and evaluate the overall measurement error. The idea is to validate VNA-based measurements of permittivity and conductivity across a range of frequencies.

The code reads experimental data from CSV files, computes error metrics, and generates visual comparisons between measured and reference properties.

---

## What this project does

- Reads dielectric measurement data from CSV files
- Analyzes the real and imaginary parts of permittivity
- Converts permittivity to conductivity
- Compares measured values to reference values
- Calculates percent error across the frequency range
- Generates plots for visual comparison

This provides a simple way to validate measurements and assess system or probe accuracy.

---

## Files in this repository

| File | Description |
|------|-------------|
| [`PlotSaline.m`](https://github.com/EMT-Lab/VNA-Saline-Validation/blob/main/PlotSaline.m) | Automatically finds multiple saline measurement files, averages the data, computes errors, and generates comparison plots. |
| [`VNA_Validation_Saline.m`](https://github.com/EMT-Lab/VNA-Saline-Validation/blob/main/VNA_Validation_Saline.m) | Analyzes a single saline dataset, computes error metrics, and produces comparison plots. |
| `saline_test2.csv` | Example measurement dataset used for testing and validation. |

---

## How to use

1. Do later
