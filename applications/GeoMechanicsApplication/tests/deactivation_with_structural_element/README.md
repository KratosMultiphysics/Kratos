# Deactivation workflow test with a structural element

This test consists of a single square UPwSmallStrainElement2D4N element with a LinearTrussElement2D2N from the StructuralMechanicsApplication directly attached to the right nodes of the square (no interfaces). It checks the deactivation/activation workflow of the structural element, combined with the `"input_type" : "rest"` option in the second stage.

## Setup

The test is performed in two stages, with the following common conditions for both stages:

- Constraints:
    - Displacements in the X-direction are all fixed.
    - Displacements in the X- and Y-direction of the bottom nodes are fixed
    - The water pressures are fixed to a hydrostatic gradient.

- Material:
    - The soil material is elastic according to the GeoLinearElasticPlaneStrain2DLaw.
    - The truss material is described by the TrussConstitutiveLaw.

The following staged analysis is performed:
- Stage 1: the truss element is deactivated, meaning the water-pressure induced displacement is the same for both top nodes.
- Stage 2: the truss element is activated and a normal load is applied to the top of the cube. The characteristics of the truss are chosen in such a way that the truss doubles the 'stiffness' of the right side of the cube, yielding half the displacement of the left side.
 

## Assertions

In the first stage, the truss is deactivated, so the expected y displacements are the same for nodes 3 and 4. In the second stage, the truss is activated and since it's only connected to node 4, the displacement of node 4 is half that of node 3 (the stiffness is doubled due to the truss). Therefore, the expected Y displacements for node 3 (left) are 0.05 and -0.1, and for node 4 (right) 0.05 and -0.05 in stage 1 and 2 respectively.
