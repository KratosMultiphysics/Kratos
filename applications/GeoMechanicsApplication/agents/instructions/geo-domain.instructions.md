---
description: "Use when changing GeoMechanics elements, conditions, constitutive or retention laws, processes or strategies. Collects the physical conventions documented in the custom_* README files: DOF ordering, sign conventions, units, stress states, process naming and documentation."
applyTo: ["applications/GeoMechanicsApplication/custom_elements/**", "applications/GeoMechanicsApplication/custom_conditions/**", "applications/GeoMechanicsApplication/custom_constitutive/**", "applications/GeoMechanicsApplication/custom_retention/**", "applications/GeoMechanicsApplication/custom_processes/**", "applications/GeoMechanicsApplication/custom_strategies/**", "applications/GeoMechanicsApplication/custom_geometries/**"]
---

# Domain conventions: GeoMechanicsApplication

The rules below come from the README files in the `custom_*` folders. The README is the source;
read it before you change the code it describes, and update it in the same change when the
documented behaviour changes.

## Degrees of freedom

U-Pw elements order their DOFs as all displacement DOFs first, grouped per node, then all water
pressure DOFs (`custom_elements/README.md`). For a 6-noded diff-order triangle:
`u_x1, u_y1, ..., u_x6, u_y6, p_1, p_2, p_3`.

- Diff-order elements interpolate the water pressure one order lower than the displacement;
  non-diff-order elements use the same order for both.
- Interface elements use relative displacements as strain measure. In their B-matrix the
  normal component comes before the tangential component(s) (`custom_elements/Policies.md`).
- Assembly helpers and test helpers (`AssertLHSMatrixBlocksAreNear`, ...) assume this order.

## Sign conventions and units

- Mechanics: tension positive for stress and strain (`PORE_PRESSURE_SIGN_FACTOR = 1.0` in
  `geo_mechanics_application_constants.h`, "mechanical sign convention").
- Water pressure: tension positive (`custom_retention/README.md`). `p > 0` is suction (above
  the phreatic level, unsaturated); `p <= 0` is compressive pore water pressure (saturated).
- In axisymmetric elements the radial coordinate is `x` (`custom_elements/Policies.md`).
- The application does not convert units. All values in `ProjectParameters.json` and
  `MaterialParameters.json` must use one consistent unit system (most tests use SI; some use
  kPa, see `tests/Simple_Dike_Gravity_Loading/line_body_all_stage_new_units_kPa`). State the
  unit of every new parameter in the README of its folder.

## Stress states and integration

- Strain vector layouts: 3D `[xx, yy, zz, xy, yz, xz]`; plane strain `[xx, yy, 0, xy]`;
  axisymmetric `[xx, yy, u_r/r, xy]` (`custom_elements/Policies.md`). Use the matching
  `StressStatePolicy` (`ThreeDimensionalStressState`, `PlaneStrainStressState`,
  `AxisymmetricStressState`, `Line2DInterfaceStressState`, `SurfaceInterfaceStressState`).
- Integration coefficients: `w_i * det(J)` by default; the axisymmetric modifier multiplies by
  `2 * pi * r`; the line modifier multiplies by `CROSS_AREA`.

## Governing equations

The U-Pw system and the names of its matrices (M mass, D damping, K stiffness, Q coupling,
C compressibility, H permeability) are defined in `custom_elements/README.md` and
`custom_utilities/README.md`. Use these names in code, comments and tests.

## Material parameters and laws

- Recent material parameters use the `GEO_` prefix (`GEO_YOUNGS_MODULUS_FORMULATION`,
  `GEO_PRESSURE_REFERENCE`, `GEO_FRICTION_ANGLE`, ...; see `custom_constitutive/README.md`).
  Use it for new parameters.
- Retention law parameters and their admissible ranges are listed in
  `custom_retention/README.md` (for example `VAN_GENUCHTEN_GL` may be negative). Keep the
  `Check()` bounds in the code equal to the documented ranges.
- UDSM and UMAT libraries are not part of the repository (`README.md`).

## Processes

- C++ process names end with `Process`. Prefix the name with `Geo` when a core Kratos process
  with the same name exists (for example `GeoApplyConstantScalarValueProcess` next to the core
  `ApplyConstantScalarValueProcess`).
- Each documented process has a section in `custom_processes/README.md` with its purpose, its
  JSON parameters and an example. Add a section for every new process.
- A process that must also work in the C++ settlement workflow is registered in
  `custom_workflows/dgeosettlement.cpp` with `mProcessFactory->AddCreator(...)`. That route is
  tested by the `*CppRoute` Python tests.

## Strategies and schemes

Time integration schemes derive from `GeoMechanicsTimeIntegrationScheme`
(`custom_strategies/schemes/README.md`).

- Derive new schemes from the appropriate existing scheme and override only the operations
  required by the new formulation. Reuse inherited time-factor, derivative-update and
  system-contribution implementations where applicable, and add a README section for the new
  scheme.
- In load stepping (`LoadSteppingScheme`), time represents load progression rather than
  physical transient evolution. Do not combine the load-stepping scheme with time-dependent
  processes.

