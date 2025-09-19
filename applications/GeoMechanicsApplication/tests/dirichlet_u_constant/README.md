# Prescribed displacement test with incremental and total displacement output.

This test is a displacements test of 1 UPwSmallStrainElement2D4N element.
It checks whether output of incremental and total displacements works correctly in combination with non-zero Dirichlet conditions.
The mesh is a single element as shown in the figure below:

![MeshStructure](MeshStructure.svg)

## Setup
The test is performed in a single stage of 2 steps, with the following conditions:

- Constraints:
    - Displacements on the bottom are fixed in Y direction.
    - Displacements on the left side are fixed in the X direction.
    - Displacements on the top are prescribed in the Y direction. These appear at the start time and remain constant.
- Material:
    - The material behaves elastically following the GeoLinearElasticPlaneStrain2DLaw.
- Load:
    - A vertical displacement of -0.1 [m] is instantly applied at the top

The result is a uniform strain and stress field.
- In the first step total, incremental and stage displacement are equal.
- In the second step, total and stage displacement remain constant, whilst the incremental displacement is zero.

## Assertions
The calculated displacements, incremental displacements, total displacements, strains and effective stresses from the Kratos Geomechanics calculations for steps 1 and 2 reflect the expected behavior described above.

