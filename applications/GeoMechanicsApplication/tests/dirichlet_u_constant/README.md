# Prescribed displacement test with.

This test is a displacements test of 1 UPwSmallStrainElement2D4N element.
It checks whether displacement Dirichlet conditions work with output of incremental and total displacements.
The mesh is a single element, which is displayed in the figure below:

![MeshStructure](MeshStructure.svg)

## Setup
The test is performed in two stages, with the following common conditions for both stages:

- Constraints:
    - Displacements on the bottom are fixed in all directions.
    - Displacements on the left side are fixed in the X direction.
    - Displacements on the top are prescribed in the Y direction. These appear at the start time and remain constant.
- Material:
    - The material is elastic according to the GeoLinearElasticPlaneStrain2DLaw.
- Load:
    - Apply a top vertical displacement of -0.1 [m]

The result is a uniform strain and stress field. Total, incremental and stage displacement are equal in the first step. In the second step total and stage displacement remain constant, incremental displacement is zero.

## Assertions
The calculated displacements, incremental displacements, total displacements, strains and effective stresses from the Kratos Geomechanics calculations for steps 1 and 2 reflect the solution described above.

