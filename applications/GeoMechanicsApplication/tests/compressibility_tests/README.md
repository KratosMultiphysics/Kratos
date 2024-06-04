# Water pressure test on a hexahedron with 8 or 20 nodes

This is a family of single-element regression test, using one of the following elements:
- UPwSmallStrainElement3D8N
- UPwSmallStrainFICElement3D8N
- UPwUpdatedLagrangianElement3D8N
- UPwUpdatedLagrangianFICElement3D8N
- SmallStrainUPwDiffOrderElement3D20N
- UpdatedLagrangianUPwDiffOrderElement3D20N

It checks the water pressure after a single stage of calculation on these single element cubes with a characteristic length of 5, with gravity loading. They are supposed to fail when additions to the right hand side are changed in these elements.

## Setup
-   Constraints:
    -   The X, Y and Z displacements in the bottom nodes are fixed.
    - The water pressure in the top nodes is kept at 0.0.
-   Material:
    -   The material is described using a linear elastic material with a `LinearElastic3DLaw`, a Young's modulus
        of 10000 [kPa] and a Poisson ratio of 0.2.
-   Loads:
    -   Using the `ApplyVectorConstraintTableProcess`, a negative VOLUME_ACCELERATION of 9.81 [m/s^2] is applied in the z-direction to the entire cube.
    
## Assertions
The following two assertions are made in this test:
1. The water pressure in the corner nodes is regression tested for the bottom and top nodes. The top nodes should be 0.0, while the bottom node results are all around -50.0 [kPa].