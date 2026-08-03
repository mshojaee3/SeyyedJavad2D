# Linear Elastic Ensemble Analysis of 2D Heterogeneous Representative Unit Cells (RUCs)

## Overview

This project provides an automated MATLAB–Abaqus workflow for performing ensemble analyses of two-dimensional heterogeneous representative unit cells (RUCs) containing circular inclusions. The workflow generates multiple Abaqus simulations by translating the inclusion within the unit cell, computes the mechanical response for a library of loading conditions, and statistically averages the resulting fields.

## Repository Structure

```
Main2D_LE.py
Run_LE_Ensemble.m
Plot_LE_Ensemble.m
```

## Workflow

```
Run_LE_Ensemble.m
        │
        ├── Updates simulation parameters
        ▼
Main2D_LE.py
        │
        ├── Builds geometry
        ├── Assigns materials
        ├── Meshes the model
        ├── Applies boundary conditions
        ├── Runs Abaqus
        └── Exports CSV results
        ▼
Raw nodal and summary CSV files
        │
        ├── Ensemble averaging
        ▼
LoadCase_ensemble.csv
        │
        ▼
Plot_LE_Ensemble.m
```

# Files

## Main2D_LE.py

Performs one Abaqus analysis for a single inclusion translation.

- Builds geometry and partitions inclusions.
- Assigns matrix and inclusion materials.
- Generates the finite element mesh.
- Applies displacement boundary conditions.
- Runs a linear static analysis.
- Exports nodal displacement, stress, strain, reaction forces, strain energy, and volume-averaged quantities.

## Run_LE_Ensemble.m

Automates the complete simulation campaign.

- Updates the Abaqus Python script.
- Executes Abaqus for every translation and load case.
- Reads all CSV outputs.
- Computes ensemble means and standard deviations.
- Computes deformation gradients and gradient-derived strains.
- Writes one ensemble CSV for each load case.

## Plot_LE_Ensemble.m

Visualizes ensemble-averaged results.

Produces figures of:

- Displacement
- Ensemble fluctuation
- Deformation gradient
- Stress
- Strain
- Von Mises stress
- Centre-line profiles

## Geometry

The total specimen dimensions are fixed:

```
Lx = Lx_tot / NX
Ly = Ly_tot / NY
```

The inclusion radius and mesh size scale with the unit-cell dimensions:

```
R = Rfrac * min(Lx, Ly)
mesh_size = meshFrac * min(Lx, Ly)
```

## Outputs

- Raw nodal CSV files
- Summary CSV files
- Ensemble CSV files
- Ensemble summary table

## Requirements

- Abaqus/CAE
- MATLAB

## Typical Workflow

1. Configure parameters in `Run_LE_Ensemble.m`.
2. Execute `Run_LE_Ensemble.m`.
3. Visualize results using `Plot_LE_Ensemble.m`.
