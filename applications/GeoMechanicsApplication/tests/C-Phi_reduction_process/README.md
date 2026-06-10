# $c-\phi$ reduction process

The $c-\phi$ reduction process is tested here. This test allows one to perform a real 2D computation. 
However, the end time is decreased and tolerances are increased to perform just a few steps. 
This is done to have enough data to verify the implemented process and to not over-spend CPU resources. 

## Setup
This test consists of 553 6-noded elements (using the UPwSmallStrainElement2D6N class). These elements are shown in the figure below:

![MeshStructure](mesh.svg)

The test is performed with the following conditions:

- Constraints:
    - The bottom nodes are fixed in both the X and Y direction,
    - The nodes on left and right vertical boundaries are fixed in the X direction. 
- Conditions:
    - A gravitational force is applied in Y direction.

The computation is done in two stages. The first stage is a settlement due to the gravitation and it uses K0 procedure.
The second stage is the $c-\phi$ reduction process and it is done only for two time steps. The figure below shows the calculated deformation.

![Deformation](deformation.svg)

Different material properties are used at the first and second stages. 

The material at the first stage is described using:

-  Young’s modulus: 14.0 $MPa$,
-  Poisson’s ratio: 0.3,
-  Unit weight: 20 $kN/m^3$,
-  Friction angle: 35.0 $\degree$.

From the unit weight the "DENSITY_SOLID" in the material parameters is calculated. Given the porosity of 0.0, the density is calculated using the following formula:

```math
\rho = \frac{\text{Unit weight}}{\text{g}} = \frac{20000}{9.81} = 2.0387 \times 10^3 \ \text{kg/m}^3.
```
The friction angle is used to calculate the "K0_VALUE_XX" and "K0_VALUE_ZZ" using the formula for $K_0^{nc}$ as documented in the [$K_0$ procedure process](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/GeoMechanicsApplication/custom_processes/README.md#K_0-procedure-process).

For the second stage, two distinct implementations of the Mohr–Coulomb model are employed:
- one provided by the MohrCoulomb64 library,
- one implemented internally in Kratos.

The material parameters are as follows:
- Young’s modulus: 30.0 $MPa$,
- Poisson’s ratio: 0.2,
- Cohesion: 1000.0 $Pa$,
- Friction angle: 30.0 $\degree$,
- Tension cut-off: 1000.0 $Pa$,
- Dilatancy angle: 0.0 $\degree$.

The Kratos model uses specific input keywords for these parameters, where the MohrCoulomb64 library uses an input in the form of an array of UMAT parameters. The meaning of the UMAT parameters in the array can be found [here](https://deltares.github.io/Kratos-GeoMechanicsApplication-Documentation/input_parameters/material_parameters/#mohr-coulomb-model-mohrcoulomb64dll).

## Assertions

The test asserts movement in the horizontal direction at three nodes, which are chosen because their movement is larger than the movement for most of the other nodes. 

