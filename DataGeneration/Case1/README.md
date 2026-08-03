# Abaqus 2D RVE Simulation Script

## Overview

This Python script automates the generation, analysis, and post-processing of a two-dimensional Representative Volume Element (RVE) in Abaqus. The model consists of a square matrix containing periodically shifted circular inclusions. The script supports multiple loading conditions and exports nodal displacement, reaction force, and strain energy data for further analysis.

---

# Features

- Automatic generation of a 2D RVE geometry
- Circular inclusions with arbitrary translations (`zeta1`, `zeta2`)
- Multiple repeated unit cells (`NX × NY`)
- Automatic material and section assignment
- Automatic meshing
- Multiple mechanical loading cases
- Automatic Abaqus job submission
- Export of:
  - nodal coordinates
  - nodal displacements
  - total strain energy (ALLSE)
  - boundary reaction forces

---

# Model Description

The model consists of

- Square matrix domain
- Circular inclusions
- Structured periodic arrangement
- Optional inclusion translation within the unit cell

The overall specimen dimensions are

```
Width  = NX × Lx
Height = NY × Ly
```

where

- `Lx` = unit cell width
- `Ly` = unit cell height

---

# User Parameters

The most important parameters are located at the beginning of the script.

## Load Case

```python
loadCaseName = 'TensionBiaxial'
```

Available options

- TensionUniaxial
- TensionBiaxial
- CompressiveUniaxial
- CompressiveBiaxial
- SimpleShear
- PureShear
- PureBending
- ParabolicBending

---

## Inclusion Translation

```python
zeta1 = 0.15
zeta2 = -0.15
```

These parameters shift every inclusion by

- `zeta1` in x-direction
- `zeta2` in y-direction

Periodic image circles are automatically generated so that inclusions crossing specimen boundaries remain continuous.

---

## Unit Cell Geometry

```python
Lx = 1.0
Ly = 1.0
R = 0.30 * Lx
```

where

- `Lx` = unit cell width
- `Ly` = unit cell height
- `R` = inclusion radius

---

## Number of Cells

```python
NX = 4
NY = 4
```

Total model size

```
Width  = NX × Lx
Height = NY × Ly
```

---

## Mesh

```python
mesh_size = 0.05 * Lx
```

The model is meshed using

- free triangular mesh
- quadratic elements

Element types

- CPE6 (plane strain)
- CPS6 (plane stress)

depending on

```python
STATE_2D
```

---

# Material Models

The script supports three material models selected using

```python
materialID
```

## Linear Elastic

```
HEl1
NEl1
```

Matrix

- Young's modulus = `E_m`
- Poisson ratio = `nu_m`

Inclusion

- Young's modulus = `E_i`
- Poisson ratio = `nu_i`

A linear static analysis is performed with

```
nlgeom = OFF
```

---

## Neo-Hookean Hyperelastic

```
HNe1
NNe1
```

Material constants

- C10
- D1

Finite deformation analysis

```
nlgeom = ON
```

---

## Mooney-Rivlin Hyperelastic

```
HMo1
NMo1
```

Material constants

- C10
- C01
- D1

Finite deformation analysis

```
nlgeom = ON
```

---

# Geometry Generation

The script automatically

1. Creates the rectangular specimen.
2. Generates every inclusion.
3. Adds periodic-image circles.
4. Partitions the matrix.
5. Identifies inclusion faces using centroid testing.
6. Assigns sections.

---

# Loading Cases

The following loading cases are implemented.

## TensionUniaxial

Vertical displacement applied to the top edge.

Bottom edge fixed in the vertical direction.

---

## TensionBiaxial

Horizontal and vertical tensile displacement.

---

## CompressiveUniaxial

Vertical compression.

---

## CompressiveBiaxial

Horizontal and vertical compression.

---

## SimpleShear

Bottom fixed.

Horizontal displacement applied to top.

---

## PureShear

Simultaneous horizontal and vertical shear.

---

## PureBending

Linear displacement field using an Abaqus ExpressionField.

---

## ParabolicBending

Parabolic displacement field using an ExpressionField.

---

# Boundary Conditions

Boundary nodes are automatically detected using

```python
TOL = 1e-6
```

Node sets created automatically

- LEFT
- RIGHT
- TOP
- BOTTOM
- CORNER00 (when required)

---

# Output Requests

The following Abaqus field outputs are requested

```
S
E
LE
U
RF
COORD
```

These correspond to

- Stress
- Engineering strain
- Logarithmic strain
- Displacement
- Reaction force
- Coordinates

---

# Job Submission

The script automatically

- creates the job
- submits the analysis
- waits for completion

Job names are generated from

```
LoadCase
zeta1
zeta2
```

Example

```
TensionBiaxial_zx_0p150_zy_m0p150
```

---

# Exported Results

After completion the script exports several CSV files.

## 1. Nodal Displacements

Directory

```
OUTPUT_CSV/<LoadCase>/
```

Each frame contains

- Node label
- Original coordinates
- U1
- U2

---

## 2. Total Strain Energy

Directory

```
OUTPUT_ENERGY/<LoadCase>/
```

Contains

- Step time
- ALLSE

---

## 3. Boundary Reaction Forces

Also stored in

```
OUTPUT_ENERGY/<LoadCase>/
```

For every frame the script exports

- Left boundary reaction
- Right boundary reaction
- Bottom boundary reaction
- Top boundary reaction

Both x and y components are summed.

---

# Directory Structure

```
Project/

│
├── AbaqusScript.py
├── OUTPUT_CSV/
│     └── LoadCase/
│            *.csv
│
└── OUTPUT_ENERGY/
      └── LoadCase/
             *_ALLSE.csv
             *_RF_sides.csv
```

---

# Typical Workflow

1. Select the loading case.
2. Specify material model.
3. Define inclusion translation.
4. Set mesh size.
5. Run the Abaqus script.
6. Wait for analysis completion.
7. Post-process exported CSV files.

---

# Notes

- Inclusion faces are identified using centroid classification.
- Periodic image circles ensure inclusions intersecting specimen boundaries are represented correctly.
- The script is intended for automated parameter studies involving different loading conditions and inclusion translations.
- Linear elastic analyses use small-strain kinematics (`nlgeom=OFF`), while hyperelastic analyses use finite deformation (`nlgeom=ON`).
