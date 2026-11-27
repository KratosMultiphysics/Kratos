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

Every new `.h`, `.hpp` and `.cpp` file (sources and tests) starts with this block. Put the result
of `git config user.name` on the `Main authors` line.

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
//  Main authors:    <result of `git config user.name` goes here>
//
```

Then `#pragma once` for headers. Two legacy headers still use `#if !defined` guards; do not copy that.

## Formatting

Format changed C++ code as described in `agents/skills/geo-format/SKILL.md`. The settings are
in `.clang-format`.

## Comments

- Comment only to explain non-obvious intent, constraints or reasoning; do not restate the code.
- Do not annotate obvious standard-library includes, such as
  `#include <vector> // Required for std::vector`. Use `#include <vector>` without a comment.
  Add an include comment only when it explains a non-obvious dependency or requirement.

## Language level

The targets are built with `CXX_STANDARD 20` (see `CMakeLists.txt`). C++20 features are used
and welcome: `std::ranges` algorithms, `requires` clauses, `using enum`, designated initializers.
The `.clang-format` file says `Standard: c++17`; that only affects formatting.

## String literals

- Prefer `std::string_literals` (`"text"s`) when explicitly constructing a `std::string`
  improves clarity or follows nearby code. Do not flag implicit conversion from a string
  literal as a defect solely for style, and do not add the suffix where a string object is
  unnecessary.

## Move semantics

- Use `std::move` when intentionally transferring a value or ownership to a destination that
  supports moving, and the source's previous value is no longer needed. A typical case is
  moving a by-value string, container or smart-pointer parameter into a member or base class.
- Check the destination's parameter type: passing `std::move(value)` to a `const T&` parameter
  does not move the value. Moving a const object usually copies; do not add `std::move` when
  it has no effect.
- Do not move from caller-owned objects unless the interface explicitly permits consumption.
  After a move, do not rely on the source's previous value; only use operations allowed by
  that type's moved-from contract.
- Return local values and by-value parameters directly (`return result;`), without
  `std::move`, so copy elision or implicit move can apply. Do not wrap temporaries in
  `std::move`.
- Include `<utility>` when using `std::move`. In reviews, suggest a move only when the transfer
  is safe and avoids an identifiable copy or is required for move-only ownership, not as a
  blanket style rule.

## Naming and files

- File names are `snake_case` and match the main class. Use the `.h` extension when the header
  contains the interface only (so no function bodies). When it does contain function bodies
  (which may be the case for member template functions), use the `.hpp` extension, to signal
  that there is code in the file.
- Data members that are references or pointers combine the prefixes: `mr` for a reference
  (`mrModelPart`) and `mp` for a pointer (`mpElement`), as in the Kratos Core.
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

- In element, condition and law `Check()` functions, use `CheckProperties` and `CheckUtilities`
  from `custom_utilities/check_utilities.hpp` when validating properties, as in `Pw_element.cpp`:

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



