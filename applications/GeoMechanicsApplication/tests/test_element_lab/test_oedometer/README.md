# Oedometer Test

This test simulates an oedometer test with a prescribed stepwise top loading process on a linear elastic soil model. It replicates a laboratory test used to determine soil properties, such as consolidation characteristics.
In the laboratory, this is performed on a cylindrical volume of soil, where an increasing pressure is applied on the top of the cylinder. In the model test, the soil sample is emulated by two axisymmetric differential order elements.

A schematic overview of the model is displayed in the figure below. Here the nodes are displayed as black dots. The two elements are numbered in red and are separated by the red line. The symmetry axis is displayed as a red dashed line.

![MeshStructure](drawing.svg)

## Setup

The test is performed under the following conditions:

- **Constraints**:
    - The bottom nodes (1, 2, 3) are fixed in the X and Y directions.
    - The side nodes and the symmetry axis (3, 6, 9 and 1, 4, 7) are fixed in the X direction.
- **Material**:
    - The material is described by the linear elastic model with the following parameters:
        - Poisson's ratio = 0.25
        - Young's modulus = 1e+07 Pa
- **Loading Conditions**:
    - A top load is applied in 4 equal steps:
        - t = 0.01 to 0.25: load = -2.5e+05 Pa
        - t = 0.26 to 0.50: load = -5e+05 Pa
        - t = 0.51 to 0.75: load = -7.5e+05 Pa
        - t = 0.76 to 1.0: load = -1e+06 Pa

## Assertions
For this test, the stresses at t = 1.0 are verified.
Additionally, the displacement in the Y-direction of the top nodes is verified at four time steps: t = 0.25, 0.5, 0.75, and 1.0.