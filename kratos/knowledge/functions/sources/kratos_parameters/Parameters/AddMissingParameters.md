---
type: C++ Method
title: AddMissingParameters
resource: sources/kratos_parameters.cpp#L1408-L1423
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/IsSubParameter
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/pointer_vector_map_iterator/PointerVectorMapIterator/key
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/future/linear_solvers/cg_solver/CGSolver/GetDefaultParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/linear_solvers/direct_solver/DirectSolver/GetDefaultParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/linear_solvers/iterative_solver/IterativeSolver/GetDefaultParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/linear_solvers/skyline_lu_factorization_solver/SkylineLUFactorizationSolver/GetDefaultParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/python_scripts/distributed_import_model_part_utility/DistributedImportModelPartUtility/ImportModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/python_scripts/distributed_restart_utility/DistributedRestartUtility/_GetDefaultParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/modelers/delete_model_parts_modeler/DeleteModelPartsModeler/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/modelers/import_mdpa_modeler/ImportMDPAModeler/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_kratos_parameters/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_kratos_parameters/TestParameters/test_add_missing_parameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/single_import_model_part/SingleImportModelPart/Import
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Parameters::AddMissingParameters(const Parameters& rDefaultParameters)`

# Calls

- [IsSubParameter](../../../../functions/sources/kratos_parameters/Parameters/IsSubParameter.md)
- [key](../../../../functions/containers/pointer_vector_map_iterator/PointerVectorMapIterator/key.md)

# Called by

- [GetDefaultParameters](../../../../functions/future/linear_solvers/cg_solver/CGSolver/GetDefaultParameters.md)
- [GetDefaultParameters](../../../../functions/future/linear_solvers/direct_solver/DirectSolver/GetDefaultParameters.md)
- [GetDefaultParameters](../../../../functions/future/linear_solvers/iterative_solver/IterativeSolver/GetDefaultParameters.md)
- [GetDefaultParameters](../../../../functions/future/linear_solvers/skyline_lu_factorization_solver/SkylineLUFactorizationSolver/GetDefaultParameters.md)
- [ImportModelPart](../../../../functions/mpi/python_scripts/distributed_import_model_part_utility/DistributedImportModelPartUtility/ImportModelPart.md)
- [_GetDefaultParameters](../../../../functions/mpi/python_scripts/distributed_restart_utility/DistributedRestartUtility/_GetDefaultParameters.md)
- [__init__](../../../../functions/python_scripts/modelers/delete_model_parts_modeler/DeleteModelPartsModeler/__init__.md)
- [__init__](../../../../functions/python_scripts/modelers/import_mdpa_modeler/ImportMDPAModeler/__init__.md)
- [ValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_kratos_parameters/KRATOS_TEST_CASE_IN_SUITE.md)
- [test_add_missing_parameters](../../../../functions/tests/test_kratos_parameters/TestParameters/test_add_missing_parameters.md)
- [Import](../../../../functions/utilities/single_import_model_part/SingleImportModelPart/Import.md)