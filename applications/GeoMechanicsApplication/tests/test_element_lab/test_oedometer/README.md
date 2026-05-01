# Oedometer test

This test is an oedometer test with a prescribed top load on a linear elastic model. It mimics a lab test, where soil properties such as the consolidation are determined.
In the lab this is performed on a cylindric volume of soil, where an increasing pressure is applied on the top of the cylinder. In the model test, the soil sample is emulated by two 2 differential order elements.

A schematic overview of the model is displayed in the figure below. Here the nodes are displayed as black dots. The 2 elements are numbered in red and are seperated by the red line.

![MeshStructure](drawing.svg)

## Setup

The test is performed with the following conditions:

- Constraints:
    - The bottom nodes (1, 2, 3) are fixed in the X and Y direction.
    - The side nodes (1, 4, 7 and 3, 6, 9) are fixed in the X direction.
- Material:
    - The material is described by the linear elastic model with the following parameters:
        - Poisson ratio = 0.0,
        - Young's modulus = 10000 $kN/m^2$
- Conditions:
    - A top load with a value of -1000 $kN/m^2$ is applied.

## Assertions
For this test the normal stresses at $t$ = 1 are asserted.
On top of this, the displacement of the top nodes is asserted at 4 time steps: $t$ = 0.25, 0.5, 0.75 and 1.