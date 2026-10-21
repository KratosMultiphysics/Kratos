---
description: "Use when writing or reviewing C++ code in applications/GeoMechanicsApplication. Covers the GeoMechanics file header, clang-format, C++20 usage, Geo aliases and constants, checks, and where GeoMechanics practice differs from the general Kratos C++ conventions."
applyTo: ["applications/GeoMechanicsApplication/**/*.cpp", "applications/GeoMechanicsApplication/**/*.h", "applications/GeoMechanicsApplication/**/*.hpp"]
---

# C++ conventions: GeoMechanicsApplication

The general Kratos rules in `/agents/instructions/cpp-conventions.instructions.md` apply
(naming with `m`/`r`/`p` prefixes, `PascalCase` classes and methods, `snake_case` locals,
`KRATOS_ERROR_IF` for checks, no `std::cout`). This file adds the GeoMechanics rules and lists
where this application differs.

## File header

Every new `.h`, `.hpp` and `.cpp` file (sources and tests) starts with this block. Put the real
author names on the `Main authors` lines, one per line.

```cpp
// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    First Author
//                   Second Author
//
```

Then `#pragma once` for headers. Two legacy headers still use `#if !defined` guards; do not copy that.

## Formatting

- Format every changed C++ file with `clang-format` using
  `applications/GeoMechanicsApplication/.clang-format`. Main settings: 4-space indent, no tabs,
  column limit 100, `BreakBeforeBraces: Linux` (function and class braces on a new line,
  control-statement braces on the same line), `PointerAlignment: Left`, consecutive assignments
  and declarations aligned, `SeparateDefinitionBlocks: Always`, no namespace indentation.
- Reformat only the code you touch. Whole-file reformat commits are listed in
  `.git-blame-ignore-revs`; add the hash there if you make one.
- Use `// clang-format off` / `// clang-format on` only for tables that must stay aligned
  (for example variable registration lists).

## Language level

The targets are built with `CXX_STANDARD 20` (see `CMakeLists.txt`). C++20 features are used
and welcome: `std::ranges` algorithms, `requires` clauses, `using enum`, designated initializers.
The `.clang-format` file says `Standard: c++17`; that only affects formatting.

## Naming and files

- File names are `snake_case` and match the main class. New headers in this application use
  `.hpp` or `.h`; follow the neighbouring files of the folder.
- Exported classes use `KRATOS_API(GEO_MECHANICS_APPLICATION)`. Test helper classes in
  `test_setup_utilities/` use `KRATOS_API(KRATOS_GEO_TEST_UTILS)`.
- Application-specific variables and material parameters added for new features use the
  `GEO_` prefix (for example `GEO_FRICTION_ANGLE`, `GEO_COHESION`). Do not rename existing
  variables.
- Utilities are classes with `static` member functions (`CheckUtilities`, `StressStrainUtilities`,
  `GeoElementUtilities`, ...). Add a function to the existing utility class that matches its
  topic before creating a new class.

## Shared types and constants

- Use the aliases in `geo_aliases.h` (namespace `Kratos::Geo`), for example
  `Geo::ConstVariableRefs`, `Geo::IntegrationPointVectorType`, `Geo::PropertiesGetter`,
  `Geo::BMatricesGetter`. Add a new alias there when it is shared by more than one class.
- Use the constants in `geo_mechanics_application_constants.h` (`N_DIM_2D`, `N_DIM_3D`,
  `VOIGT_SIZE_3D`, `VOIGT_SIZE_2D_PLANE_STRAIN`, `STRESS_TENSOR_SIZE_2D`, `TINY`, `LARGE`, ...)
  and its enums (`IsDiffOrderElement`, `PlasticityStatus`, `DrainageType`) instead of literal
  numbers.

## Input checks

- Element, condition and law `Check()` functions use `CheckProperties` and `CheckUtilities`
  from `custom_utilities/check_utilities.hpp`, as in `Pw_element.cpp`:

  ```cpp
  const CheckProperties check_properties(r_properties, "material properties at element",
                                         this->Id(), CheckProperties::Bounds::AllInclusive);
  check_properties.Check(DENSITY_WATER);
  check_properties.Check(POROSITY, max_value_porosity);
  ```
- Error messages name the variable, the value, the allowed range and the element or property Id,
  and end with `std::endl`.

## Design patterns in use

- Element behaviour that varies with the stress state or the integration volume is a policy
  object (`StressStatePolicy`, `IntegrationCoefficientModifier`), described in
  `custom_elements/Policies.md`. Add a new policy implementation instead of `if` branches on
  element type.
- New element LHS/RHS terms go into a contribution calculator in
  `custom_elements/contribution_calculators/` (see its `README.md`): a class with an
  `InputProvider` of getter functions, returning `std::optional<LHSMatrixType>` and a bounded
  `RHSVectorType`.
- Interface elements use the interface geometries in `custom_geometries/` (see its `README.md`
  for node numbering).

## Where GeoMechanics practice differs from the general Kratos rules

- `KRATOS_TRY` / `KRATOS_CATCH` are used in about 70 of roughly 390 source files. They are not
  required in every method. Use them in the long Kratos interface methods of elements,
  conditions and processes when the surrounding class already does.
- Pointer creation in `Create()` and `Clone()` of legacy elements and conditions still uses
  `new` inside `Kratos::make_intrusive`/`Element::Pointer(new ...)`. New code uses
  `Kratos::make_intrusive`, `Kratos::make_shared` or `std::make_unique`.
- Prefer `std::ranges` algorithms (`std::ranges::transform`, `find_if`, `adjacent_find`, ...)
  over raw index loops when they make the intent clearer. Mark pure query functions
  `[[nodiscard]]`.




