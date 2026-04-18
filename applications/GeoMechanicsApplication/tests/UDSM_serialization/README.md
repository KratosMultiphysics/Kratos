# UDSM serialization

Compact 2D column benchmarks for verifying that a UDSM-based constitutive law serializes and restarts correctly.

## Test outline
 - Run a full three-stage simulation and save checkpoint dumps after each stage.
 - Run a restart simulation that loads the stage 2 checkpoint and executes stage 3 only.
 - Compare the final nodal `TOTAL_DISPLACEMENT` and `WATER_PRESSURE` fields from the two workflows: they should match within numerical tolerance.

## Mesh & elements
 - Diff-order elements `SmallStrainUPwDiffOrderElement2D6N` (4 triangles with midside nodes).
 - Setup illustration: ![Setup](setup.svg)


## Geometry and boundary conditions
 - Domain: $1\ \mathrm{m}$ (width) by $2\ \mathrm{m}$ (height).
 - Bottom: fixed in both directions.
 - Left and right faces: fixed in horizontal direction.
 - Stage 1: gravity loading and $K_0$ initialization.
 - Stages 2 & 3: UDSM constitutive law; top boundary uses a prescribed vertical displacement (table to $-0.001\ \mathrm{m}$) and top water pressure is prescribed to $10000\ \mathrm{Pa}$.

## Material models

The entries below list the physical symbol, a short description, and the corresponding JSON key used in the case material files.

### Stage 1 (initialization)

 - Constitutive law: `GeoLinearElasticPlaneStrain2DLaw`
 - $E = 1.0e9\ \mathrm{[Pa]}$, 
 - $\nu = 0.2$
 - $\rho_s = 2000\ \mathrm{[kg/m^3]}$, $\rho_w = 1000\ \mathrm{[kg/m^3]}$
 - $n = 0.3$, $K_s = 1.0e12\ \mathrm{[Pa]}$, $K_w = 2.2e9\ \mathrm{[Pa]}$
 - $k_{xx}=k_{yy}=4.5e-10\ \mathrm{[m^2]}$, 
 - $\mu=1.0e-3\ \mathrm{[Pa\cdot s]}$
 - $\mathrm{K0\_MAIN\_DIRECTION}$: index of K0 main direction (1)
 - $\mathrm{K0\_VALUE\_XX},\;\mathrm{K0\_VALUE\_YY},\;\mathrm{K0\_VALUE\_ZZ}$: K0 values used by the K0 procedure (e.g. 0.5, 1.0, 0.5)
 - $k_{xy}$ (in-plane cross permeability): 0.0
 - $S_{\mathrm{sat}}$ (saturated saturation): 1.0

### Stages 2 and 3 (UDSM)

 - Constitutive law: `SmallStrainUDSM2DPlaneStrainLaw`
 - Same elastic/hydraulic keys as stage 1 (`YOUNG_MODULUS`, `POISSON_RATIO`, `POROSITY`, `PERMEABILITY_*`, `DYNAMIC_VISCOSITY`, ...)
 - UMAT parameters (Props): [
                        4.998729486706428e-02,
                        0.19999260891644746,
                        1.999926089164475e-02,
                        8.64e+04,
                        1.5]

## Assertions
 - The test compares final `TOTAL_DISPLACEMENT` and `WATER_PRESSURE` fields from the full run and the restart run for the third stage. Differences should be within numerical tolerance.
