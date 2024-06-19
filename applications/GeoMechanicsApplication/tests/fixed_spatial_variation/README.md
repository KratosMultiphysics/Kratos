# Prescribed spatially varying fixed water pressure

This test is a water pressure tests of a single UPwSmallStrainElement2D8N element.
It checks whether fixed water pressure conditions stay fixed when using a spatially varying fixed water pressure. In this case, it means we use the `ApplyScalarConstraintTableProcess` with the following tables (x, water_pressure pairs) for the top and bottom sides of the element respectively:
```
Begin Table 1 X WATER_PRESSURE
0.0000000000 0.0000000000
1.0000000000 -5000
End Table

Begin Table 2 X WATER_PRESSURE
0.0000000000 -40000
1.0000000000 -30000
End Table
```

## Setup
- Constraints:
    - Displacements on the bottom are fixed in all directions.
    - Displacements on the sides are fixed in the X direction.
    - Water pressures on the top are varying from 0 to -5000, and on the bottom from -40000 to -30000. These are constant in time.

- Material:
    - The material is elastic according to the GeoLinearElasticPlaneStrain2DLaw.

A VOLUME_ACCELERATION is also set, to make sure that if the water pressure is not fixed, the water pressure would get different values.

## Assertions
The assertion is very simple: the water pressure should be constant in time and space for the top and bottom nodes, and equal to the prescribed values.
