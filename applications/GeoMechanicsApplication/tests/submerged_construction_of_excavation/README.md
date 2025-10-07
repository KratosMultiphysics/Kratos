# Submerged construction of an excavation

This document describes the changes that have been made to the Kratos model compared to the corresponding Plaxis tutorial.


## Material properties

GeoMechanicsApplication expects a slightly different set of material properties compared to Plaxis.

In Plaxis, the user needs to provide the unsaturated unit weight $`\gamma_{\mathrm{unsat}}`$ and saturated unit weight $`\gamma_{\mathrm{sat}}`$ of the clay layer as well as the sand layer.  GeoMechanicsApplication, however, expects to receive solid densities and porosity values.  The following formulas show how to calculate the unsaturated unit weight as well as the saturated unit weight:

```math
\gamma_{\mathrm{unsat}} = (1 - n) \cdot \rho_{\mathrm{g}} \cdot g + S_{\mathrm{res}} \cdot n \cdot \rho_{\mathrm{w}} \cdot g \\
\gamma_{\mathrm{sat}}   = (1 - n) \cdot \rho_{\mathrm{g}} \cdot g + S_{\mathrm{sat}} \cdot n \cdot \rho_{\mathrm{w}} \cdot g
```

given porosity $`n`$, grain density $`\rho_{\mathrm{g}}`$, residual saturation $`S_{\mathrm{res}}`$, "saturated" saturation $`S_{\mathrm{sat}}`$, water density $`\gamma_{\mathrm{w}}`$, and gravity acceleration $`g`$.

For the initial stage (which uses a $`K_0`$ procedure), we need to assign linear elastic material properties which are defined by a Young's modulus $`E`$ and a Poisson's ratio $`\nu`$.  Since the Plaxis tutorial specifies three different stiffness values (i.e., a secant stiffness in standard drained triaxial test $`E_{50}^{\mathrm{ref}}`$, a tangent stiffness for primary oedometer loading $`E_{\mathrm{oed}}^{\mathrm{ref}}`$, and an unloading / reloading stiffness $`E_{\mathrm{ur}}^{\mathrm{ref}}`$), we have opted for the maximum value of these (being the unloading / reloading stiffness $`E_{\mathrm{ur}}^{\mathrm{ref}}`$).

The following table lists the material properties of the soil layers that have been adopted by the Kratos model.

| Property                                  | Kratos input parameter | Clay | Sand    | Unit                           |
|-------------------------------------------|------------------------|------|---------|--------------------------------|
| Grain density $`\rho_{\mathrm{g}}`$       | `DENSITY_SOLID` | 2048.66 | 2496.33 | $`\mathrm{kg} / \mathrm{m}^3`$ |
| Water density $`\rho_{\mathrm{w}}`$       | `DENSITY_WATER` | 1000.0 | 1000.0 | $`\mathrm{kg} / \mathrm{m}^3`$ |
| Porosity $`n`$                            | `POROSITY` | 0.203874 | 0.305810 | $`[-]`$                        |
| Retention law type                        | `RETENTION_LAW` | `SaturatedBelowPhreaticLevelLaw` | `SaturatedLaw` | N/A                            |
| Residual saturation $`S_{\mathrm{res}}`$  | `RESIDUAL_SATURATION` | $`1 \cdot 10^{-10}`$ | $`1 \cdot 10^{-10}`$ | $`[-]`$                        |
| Saturated saturation $`S_{\mathrm{sat}}`$ | `SATURATED_SATURATION` | 1.0 | 1.0 | $`[-]`$                        |
| Young's modulus $`E`$ | `YOUNG_MODULUS` | $`12 \cdot 10^3`$ | $`120 \cdot 10^3`$ | $`\mathrm{kN} / \mathrm{m}^2`$ |
| Poisson's ratio $`\nu`$ | `POISSON_RATIO` | 0.15 | 0.20 | $`[-]`$ |


## Staged analysis

This section lists the stages that are taken into account:

1. **Initial stage:**: in this stage, the complete soil domain is active, and the diaphragm wall as well the interfaces at both sides of it are inactive.  Where the soil will later be separated by the diaphragm wall, we will apply master-slave constraints to ensure continuity of the displacement field in this early stage of analysis.  The only load that is being applied is self-weight.  At the end of the stage, a $`K_0`$ procedure is performed to initialize the horizontal stress field.  Note that the $`K_0`$ procedure requires the use of linear elastic materials for all soil parts.
2. **Null step:** in this stage, we will change the material models from linear elastic to hardening soil (or Mohr-Coulomb if we don't have a working hardening soil model at hand).  By changing the material models, we will trigger a stiffness redistribution, and hence a stress redistribution.
3. **Diaphragm wall installation and applying the external load:** in this stage, we activate the diaphragm wall and the interfaces that are attached to its left and right sides.  At the same time, we also need to deactivate the master-slave constraints, since we no longer require the displacement field to be continuous across the entire domain.  In addition, a surface load is applied to a part of the top of the soil on the left hand side.
4. **First excavation stage:** In this stage, the top part of the clay to the right of the diaphragm wall is excavated.  In other words, the corresponding model part is deactivated.  Also, the interface elements that connect the now excavated clay layer to the diaphragm wall need to be deactivated. 
5. **Installation of a strut:** In this stage, we activate the strut.
6. **Second excavation stage:** In this stage, the now top-most part of the clay to the right of the diaphragm wall is excavated.  In other words, the corresponding model part is deactivated.  Also, the interface elements that connect the now excavated clay layer to the diaphragm wall need to be deactivated.  In addition, we need to apply a normal contact stress to the "naked" part of the diaphragm wall as well as the bottom of the excavation pit.
7. **Third excavation stage:** In this final stage, the now top-most part of the clay to the right of the diaphragm wall is excavated.  In other words, the corresponding model part is deactivated.  Also, the interface elements that connect the now excavated clay layer to the diaphragm wall need to be deactivated.  In addition, we need to apply a normal contact stress to the "naked" part of the diaphragm wall as well as the bottom of the excavation pit.