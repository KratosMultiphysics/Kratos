# Grouted anchor test

This test models a 2D plane-strain soil body with a sheetpile wall that is held by a grouted anchor. The anchor consists of a free (ungrouted) part, which is a truss attached to the sheetpile, and a grouted part, which is a truss connected to the surrounding soil by a line interface element. It checks the interaction between soil, a beam, truss elements and a grout-soil interface in a single static stage.

## Setup

The test consists of a single stage (`1_Initial_stage`) with a static, drained analysis (`U_Pw`, Newton-Raphson, one time step of 1.0 s). A schematic depiction of the model is shown below.

![Mesh.svg](Mesh.svg)

- Geometry and elements (see `Sheetpile_with_anchor.mdpa`):
    - `Soil` (yellow in schematic): `SmallStrainUPwDiffOrderElement2D6N`.
    - `Sheetpile` (blue in the schematic): `LinearTimoshenkoCurvedBeamElement2D3N`.
    - `Anchor_grouted` (brown in the schematic): `LinearTrussElement2D3N`.
    - `Anchor_free` (green in the schematic): `LinearTrussElement2D2N`.
    - `Grout_soil_interface` (between the grouted anchor and the edges of the soil elements): `Geo_UPwLineInterfacePlaneStrainDiffOrderElement3Plus3N`.

- Constraints:
    - Displacements in the X- and Y-direction of the bottom boundary (`Bottom_Fixed`) are fixed.
    - Displacements in the X-direction of the left and right edges (`Left_Edge`, `Right_Edge`) are fixed.
    - The water pressure in the soil is fixed to 0.

- Loads:
    - Gravity is applied as a `VOLUME_ACCELERATION` of -9.81 m/s² in the Y-direction.

- Materials (see `MaterialParameters.json`):
    - Soil: `GeoLinearElasticPlaneStrain2DLaw` (Young's modulus 2.0e7 Pa, Poisson's ratio 0.3) with a constant pore-water pressure field.
    - Sheetpile: `TimoshenkoBeamElasticConstitutiveLaw` (Young's modulus 2.1e11 Pa, cross area 1.5e-2 m², I33 3.0e-4 m⁴).
    - Grouted anchor part: `TrussConstitutiveLaw` (Young's modulus 2.1e11 Pa, cross area 4.938e-2 m²).
    - Free anchor part: `TrussConstitutiveLaw` (Young's modulus 2.1e11 Pa, cross area 4.938e-3 m²).
    - Grout-soil interface: `GeoIncrementalLinearElasticInterfaceLaw` (normal stiffness 7.69e8 Pa/m, shear stiffness 7.69e7 Pa/m).

## Assertions

The test (`test_grouted_anchor.py`) runs the stage with the sequential orchestrator and reads the `DISPLACEMENT` results from `Sheetpile_with_anchor.post.res`. At time 1.0, the displacement vectors of nodes 1 to 5 are compared with regression values:

| Node | Displacement X [m] | Displacement Y [m] |
|------|--------------------|--------------------|
| 1    | -0.00075417        | -0.104971          |
| 2    | -2.74607e-05       | -0.105343          |
| 3    | 0.000989452        | -0.105352          |
| 4    | -0.00286768        | -0.109109          |
| 5    | 7.98339e-05        | -0.10327           |

These five nodes are also displayed in the schematic. The Z-displacement is 0 for all nodes. These values are regression results rather than analytical solutions.
