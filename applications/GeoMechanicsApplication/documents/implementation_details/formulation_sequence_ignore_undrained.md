# IGNORE_UNDRAINED formulation implementation sequence

This note documents how the `IGNORE_UNDRAINED` property is implemented in the GeoMechanicsApplication and how each U–Pw formulation uses it during checks, initialization, and assembly.

The `IGNORE_UNDRAINED` flag is used to disable the hydraulic flow part of the coupled hydro-mechanical formulation.
When enabled, the pressure-flow equations are skipped while the mechanical response and displacement–pressure coupling remain active.

## Summary table

| Formulation | Allowed | Property must be present | Missing when `true` skips validation of | Unconditionally read despite flag | Skipped LHS terms | Skipped RHS terms |
|---|---|---|---|---|---|---|
| `UPwSmallStrainElement` | Yes | Yes | `BULK_MODULUS_FLUID`, `DYNAMIC_VISCOSITY`, permeability | `DYNAMIC_VISCOSITY`, permeability (in `InitializeProperties`) | PP permeability, PP compressibility, PU coupling | Compressibility flow, permeability flow, fluid body flow, pressure-side coupling |
| `UPwUpdatedLagrangianElement` | Yes | Yes | (inherits from above) | (inherits from above) | (inherits from above) | (inherits from above) |
| `SmallStrainUPwDiffOrderElement` | Yes | Yes | Permeability | `DYNAMIC_VISCOSITY`, permeability (in `InitializeProperties`) | PP permeability, PP compressibility, PU coupling | Compressibility flow, permeability flow, fluid body flow, pressure-side coupling; pressure transfer to intermediate nodes |
| `UpdatedLagrangianUPwDiffOrderElement` | Yes | Yes | (inherits from above) | (inherits from above) | (inherits from above) | (inherits from above) |
| `UPwSmallStrainInterfaceElement` (legacy) | Yes | Yes | `TRANSVERSAL_PERMEABILITY`, `BULK_MODULUS_FLUID`, `DYNAMIC_VISCOSITY` | `DYNAMIC_VISCOSITY` (in `InitializeElementVariables`), `BULK_MODULUS_FLUID` (in per-GP `BiotModulusInverse`) | Coupling matrix, compressibility matrix, permeability matrix | Compressibility flow, permeability flow, fluid body flow, pressure-side coupling |
| `UPwInterfaceElement` (modular) | Yes | No (defaults to `false`) | — | — | `PUCoupling`, `Permeability` | `PUCoupling`, `Permeability`, `FluidBodyFlow` |
| `UPwSmallStrainFICElement` | **No** | Yes | — | — | — | — |
| `UPwUpdatedLagrangianFICElement` | **No** | Yes | — | — | — | — |

## 1. Variable definition and registration

1. Declared as application variable in:
   - `geo_mechanics_application_variables.h`
2. Defined in:
   - `geo_mechanics_application_variables.cpp`
3. Registered in C++ and Python bindings in:
   - `geo_mechanics_application.cpp`
   - `custom_python/geo_mechanics_python_application.cpp`
4. Set by material input, for example in `MaterialParameters*.json` as:

```json
"IGNORE_UNDRAINED": true
```

## 2. Common implementation pattern

Across U–Pw formulations the use of `IGNORE_UNDRAINED` follows a consistent sequence.

1. **Check stage**

   `Check(...)` verifies that the property exists.  
   When `IGNORE_UNDRAINED == false`, additional hydraulic properties are validated (e.g. permeability, fluid bulk modulus, viscosity).

2. **Initialization**

   During element initialization the property value is copied into an element variable: `IgnoreUndrained`

3. **Assembly**

During LHS/RHS assembly:

- If `IgnoreUndrained == true`, pressure-flow equation terms are skipped.
- If `IgnoreUndrained == false`, the full coupled hydro-mechanical formulation is assembled.

4. **Utility-level behavior**

Some utilities adjust compressibility behavior.  
When `IGNORE_UNDRAINED == true`, the inverse Biot modulus is computed using a very small fluid bulk modulus (`TINY`).

## 3. Formulation-by-formulation sequence

### 3.1 UPwSmallStrainElement (equal-order) and UPwUpdatedLagrangianElement

Main files:
- `custom_elements/U_Pw_small_strain_element.cpp`
- `custom_elements/U_Pw_updated_lagrangian_element.cpp`

Sequence:

1. `Check(...)`
   - requires `IGNORE_UNDRAINED`
   - if drained (`false`), requires `BULK_MODULUS_FLUID`, `DYNAMIC_VISCOSITY`, and permeability properties
2. `InitializeProperties(...)`
   - sets `rVariables.IgnoreUndrained`
3. `CalculateAndAddLHS(...)`
   - always: stiffness (`UU`) + U-P coupling (`UP`)
   - only when drained: permeability (`PP`), compressibility (`PP`), and P-U coupling (`PU`)
4. `CalculateAndAddRHS(...)`
   - always: stiffness force, mixture body force, coupling force to displacement
   - only when drained: compressibility flow, permeability flow, fluid body flow, coupling flow to pressure
5. Updated Lagrangian variant
   - calls base small-strain assembly first
   - then optionally adds geometric stiffness

### 3.2 SmallStrainUPwDiffOrderElement and UpdatedLagrangianUPwDiffOrderElement

Main files:
- `custom_elements/small_strain_U_Pw_diff_order_element.cpp`
- `custom_elements/updated_lagrangian_U_Pw_diff_order_element.cpp`

Sequence:

1. `Check(...)`
   - requires `IGNORE_UNDRAINED`
   - if drained (`false`), checks permeability properties
2. `InitializeProperties(...)`
   - sets `rVariables.IgnoreUndrained`
3. `CalculateAndAddLHS(...)`
   - always: stiffness (`UU`) + U-P coupling (`UP`)
   - only when drained: permeability (`PP`), compressibility (`PP`), and P-U coupling (`PU`)
4. Internal and external force construction
   - always: mechanical and U-side coupling terms
   - only when drained: pressure-flow terms (compressibility, permeability, fluid body flow)
5. `FinalizeSolutionStep(...)`
   - pressure transfer to intermediate nodes is done only when drained
6. Updated Lagrangian diff-order variant
   - reuses small-strain diff-order behavior
   - then optionally adds geometric stiffness

### 3.3 UPwSmallStrainInterfaceElement (legacy interface implementation)

Main file:
- `custom_elements/U_Pw_small_strain_interface_element.cpp`

Sequence:

1. `Check(...)`
   - requires `IGNORE_UNDRAINED`
   - if drained (`false`), checks `TRANSVERSAL_PERMEABILITY`, `BULK_MODULUS_FLUID`, `DYNAMIC_VISCOSITY`
2. `InitializeElementVariables(...)`
   - sets `rVariables.IgnoreUndrained`
3. `CalculateAndAddLHS(...)`
   - always: interface stiffness (`UU`)
   - only when drained: coupling matrix, compressibility matrix, permeability matrix
4. `CalculateAndAddRHS(...)`
   - always: stiffness force, mixture body force, coupling force to displacement
   - only when drained: compressibility flow, permeability flow, fluid body flow, coupling flow to pressure

### 3.4 UPwInterfaceElement (modular interface implementation)

Main file:
- `custom_elements/U_Pw_interface_element.cpp`

Sequence:

1. `GetIgnoreUndrained(...)`
   - reads property when available
   - defaults to `false` when property is missing
2. LHS contribution loop
   - always includes contributions configured in `mContributions`
   - when `IGNORE_UNDRAINED == true`, explicitly skips:
     - `PUCoupling`
     - `Permeability`
3. RHS contribution loop
   - when `IGNORE_UNDRAINED == true`, explicitly skips:
     - `PUCoupling`
     - `Permeability`
     - `FluidBodyFlow`

### 3.5 FIC formulations

Main file:
- `custom_elements/U_Pw_small_strain_FIC_element.cpp`

Sequence:

1. `Check(...)` throws when `IGNORE_UNDRAINED == true`
2. Message: use non-FIC elements for `IGNORE_UNDRAINED`
3. Updated Lagrangian FIC inherits from small-strain FIC, so this restriction also applies there

## 4. Utility behavior tied to IGNORE_UNDRAINED

Main utility files:
- `custom_utilities/transport_equation_utilities.cpp`
- `custom_elements/contribution_calculators/compressibility_calculator.hpp`

Behavior:

1. Inverse Biot modulus is computed with a very small fluid bulk modulus (`TINY`) when `IGNORE_UNDRAINED == true`
2. This makes inverse Biot modulus very large, which is consistent with suppressing normal drained compressibility behavior
3. Unit tests verify this behavior in:
   - `tests/cpp_tests/custom_utilities/test_transport_equation_utilities.cpp`

## 5. Typical staged analysis workflow

Example:
- `tests/dsettlement/fully_saturated_column_uniform_load/README.md`

Observed usage pattern:

1. stage with `IGNORE_UNDRAINED = true` to keep hydrostatic pressure unchanged
2. later stage with `IGNORE_UNDRAINED = false` to enable consolidation and pressure evolution

## 6. Current implementation caveats

Several formulations conditionally validate flow properties in `Check(...)` but still read them unconditionally during element initialization or at every integration point.

### UPwSmallStrainElement and SmallStrainUPwDiffOrderElement (and their Updated Lagrangian variants)

`Check(...)` validates `DYNAMIC_VISCOSITY` and all permeability properties only when `IGNORE_UNDRAINED == false`.
However, `InitializeProperties(...)` always executes:

```cpp
rVariables.DynamicViscosityInverse = 1.0 / r_properties[DYNAMIC_VISCOSITY];
GeoElementUtilities::FillPermeabilityMatrix(rVariables.PermeabilityMatrix, r_properties);
```

Both reads happen regardless of `IGNORE_UNDRAINED`.
If those properties are missing, the code will throw a property-not-found exception at run time, not at the `Check(...)` stage.

Note: `BULK_MODULUS_FLUID` is safe to omit for these two formulations because they delegate Biot modulus computation to `GeoTransportEquationUtilities::CalculateInverseBiotModulus`, which explicitly substitutes `TINY` for the fluid bulk modulus when `IGNORE_UNDRAINED == true`.

### UPwSmallStrainInterfaceElement (legacy interface)

`Check(...)` validates `DYNAMIC_VISCOSITY`, `BULK_MODULUS_FLUID`, and `TRANSVERSAL_PERMEABILITY` only when `IGNORE_UNDRAINED == false`.
However:

- `InitializeElementVariables(...)` always reads:

  ```cpp
  rVariables.DynamicViscosityInverse = 1.0 / rProperties[DYNAMIC_VISCOSITY];
  ```

- A per-integration-point helper always computes `BiotModulusInverse` directly, without delegating to the utility or checking the flag:

  ```cpp
  rVariables.BiotModulusInverse =
      (rVariables.BiotCoefficient - r_properties[POROSITY]) / r_properties[BULK_MODULUS_SOLID] +
      r_properties[POROSITY] / r_properties[BULK_MODULUS_FLUID];
  ```

Both `DYNAMIC_VISCOSITY` and `BULK_MODULUS_FLUID` must therefore be present in the material properties even when `IGNORE_UNDRAINED == true`, or the element will fail at run time.

### Practical recommendation

Keep the full hydraulic property set (`DYNAMIC_VISCOSITY`, `BULK_MODULUS_FLUID`, all permeability entries) in materials even when setting `IGNORE_UNDRAINED = true`.
