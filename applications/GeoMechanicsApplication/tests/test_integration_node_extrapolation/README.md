# Integration Point to Node Extrapolation Tests

This folder contains a series of tests to test the process that extrapolates integration point values to nodes. The functionality is found in `GeoIntegrationValuesExtrapolationToNodesProcess`. The meshes used are simple rectangular domains, built up using either a single 4- or 8-noded quadrilateral, or four 3- or 6-noded triangles.

### 3-, 4-, 6-noded test cases
For these testcases the SteadyStatePwElement is used, in 2D3N, 2D6N or 2D4N configurations.

-   Constraints:
    -   The X, Y and Z displacements are fixed in the entire domain (the focus is on water pressure)
    - The water pressure at the top and the base is fixed using a `ApplyScalarConstraintTableProcess` with a "Hydrostatic" fluid pressure type.
-   Material:
    -   The material is described using a linear elastic material with a `GeoLinearElasticPlaneStrain2DLaw`, a Young's modulus
        of 10000 [kPa] and a Poisson ratio of 0.2 and can be found in the [Material Parameters](common/MaterialParameters.json).
-   Loads:
    -  Self-weight is induced by adding a -10.0 $\mathrm{[m/s^2]}$ VOLUME_ACCELERATION to the domain (using the `ApplyVectorConstraintTableProcess`)

The `GeoIntegrationValuesExtrapolationToNodesProcess` is configured to extrapolate HYDRAULIC HEAD.

### 8 noded test case
For one of the tests, a more extensive set of variables was tested for the extrapolation. 
-   Constraints:
    - The X, Y and Z displacements at the bottom are fixed.
    - The X and Z displacements are fixed in the entire domain.
    - The water pressure at the top is fixed to follow a linear profile from 0 to -5000 [kPa] in the x direction to promote a varying FLUID_FLUX_VECTOR.
    - The water pressure at the bottom is fixed to follow a linear profile from -40000 to -30000 in the x direction to promote a varying FLUID_FLUX_VECTOR.
-   Material: the material is the same as the other cases (see [Material Parameters](common/MaterialParameters.json))
-   Loads:
    -  Self-weight is induced by adding a -10.0 $\mathrm{[m/s^2]}$ VOLUME_ACCELERATION to the domain (using the `ApplyVectorConstraintTableProcess`)

The `GeoIntegrationValuesExtrapolationToNodesProcess` is configured to extrapolate the HYDRAULIC HEAD, CAUCHY_STRESS_TENSOR and FLUID_FLUX_VECTOR.

## Assertions
These tests are in essence regression tests and check that the extrapolated values stay the same as previous versions.


