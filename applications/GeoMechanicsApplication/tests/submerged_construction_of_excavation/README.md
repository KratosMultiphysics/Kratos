# Submerged construction of an excavation

This document describes the changes that have been made to the Kratos model compared to the corresponding Plaxis tutorial.


## Material properties

GeoMechanicsApplication expects a slightly different set of material properties compared to Plaxis.

In Plaxis, the user needs to provide the unsaturated unit weight $`\gamma_{\mathrm{unsat}}`$ and saturated unit weight $`\gamma_{\mathrm{sat}}`$ of the clay layer as well as the sand layer.  GeoMechanicsApplication, however, expects to receive solid densities and porosity values.  The following formulas show how to calculate the unsaturated unit weight as well as the saturated unit weight:

```math
\gamma_{\mathrm{unsat}} = (1 - n) \cdot \rho_{\mathrm{g}} \cdot g + s_{\mathrm{res}} \cdot n \cdot \rho_{\mathrm{w}} \cdot g \\
\gamma_{\mathrm{sat}}   = (1 - n) \cdot \rho_{\mathrm{g}} \cdot g + s_{\mathrm{sat}} \cdot n \cdot \rho_{\mathrm{w}} \cdot g
```

given porosity $`n`$, grain density $`\rho_{\mathrm{g}}`$, residual saturation $`s_{\mathrm{res}}`$, "saturated" saturation $`s_{\mathrm{sat}}`$, water density $`\gamma_{\mathrm{w}}`$, and gravity acceleration $`g`$.  The following table lists the material properties of the soil layers that have been adopted by the Kratos model.

| Property                                  | Kratos input parameter | Clay | Sand    | Unit                           |
|-------------------------------------------|------------------------|------|---------|--------------------------------|
| Grain density $`\rho_{\mathrm{g}}`$       | `DENSITY_SOLID` | 2048.66 | 2496.33 | $`\mathrm{kg} / \mathrm{m}^3`$ |
| Water density $`\rho_{\mathrm{w}}`$       | `DENSITY_WATER` | 1000.0 | 1000.0 | $`\mathrm{kg} / \mathrm{m}^3`$ |
| Porosity $`n`$                            | `POROSITY` | 0.203874 | 0.305810 | $`[-]`$ |
| Retention law type | `RETENTION_LAW` | `SaturatedBelowPhreaticLevelLaw` | `SaturatedLaw` | N/A |
| Residual saturation $`s_{\mathrm{res}}`$  | `RESIDUAL_SATURATION` | $`1 \cdot 10^{-10}`$ | $`1 \cdot 10^{-10}`$ | $`[-]`$ |
| Saturated saturation $`s_{\mathrm{sat}}`$ | `SATURATED_SATURATION` | 1.0 | 1.0 | $`[-]`$ |
