# Python bindings

pybind11 module `KratosContactStructuralMechanicsApplication` (imported by `KratosMultiphysics/ContactStructuralMechanicsApplication/__init__.py`, which imports the `StructuralMechanicsApplication` first).

| File | Exports |
|---|---|
| `contact_structural_mechanics_python_application.cpp` | The application class and the variables (all application variables except `CONSTRAINT_POINTER`, `PARENT_ELEMENT` and `FRICTIONAL_LAW`; `TANGENT_FACTOR` re-registered; `ACTIVE_CHECK_FACTOR` registered twice), the enum `NormalDerivativesComputation`. |
| `add_custom_strategies_to_python.cpp` | `ResidualBasedNewtonRaphsonContactStrategy`, `LineSearchContactStrategy`, `ResidualBasedNewtonRaphsonMPCContactStrategy`; the 17 convergence criteria; `ContactResidualBasedBlockBuilderAndSolver`, `ContactResidualBasedEliminationBuilderAndSolver`, `ContactResidualBasedEliminationBuilderAndSolverWithConstraints`. |
| `add_custom_processes_to_python.cpp` | `ALMFastInit`, `MasterSlaveProcess`, `ComputeDynamicFactorProcess`, `ALMVariablesCalculationProcess`, `ContactSPRErrorProcess2D/3D`, `ContactSearchProcess` (wrapper), `MPCContactSearchProcess` (wrapper), `NormalCheckProcess`, `FindIntersectedGeometricalObjectsWithOBBContactSearchProcess`, `AssignParentElementConditionsProcess`, and the templated `SimpleContactSearchProcess<…>`, `AdvancedContactSearchProcess<…>`, `MPCContactSearchProcess<…>`, `NormalGapProcess<…>` for `2D2N`, `3D3N`, `3D4N`, `3D3N4N`, `3D4N3N`. |
| `add_custom_utilities_to_python.cpp` | `ProcessFactoryUtility`, `ContactUtilities`, `InterfacePreprocessCondition`, submodules `ActiveSetUtilities` and `SelfContactUtilities`. |
| `add_custom_linear_solvers_to_python.cpp` | `MixedULMLinearSolver`. |
| `add_custom_frictional_laws_to_python.cpp` | `FrictionalLaw`, `FrictionalLaw<suffix>`, `CoulombFrictionalLaw<suffix>`, `TrescaFrictionalLaw<suffix>` (suffixes `2D2N`, `3D3N`, `3D4N`, `3D3N4N`, `3D4N3N` and their `NV` versions). |
| `process_factory_utility.h/.cpp` | `ProcessFactoryUtility`: holds a Python list of processes so that the C++ strategies (adaptive time stepping) can call their `Execute*` methods. |

Full documentation: [Architecture](https://kratosmultiphysics.github.io/Kratos/pages/Applications/Contact_Structural_Mechanics_Application/Implementation/Architecture.html) ([source](../../../docs/pages/Applications/Contact_Structural_Mechanics_Application/Implementation/Architecture.md)).
