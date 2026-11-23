---
type: C++ Method
title: AddInt
resource: sources/kratos_parameters.cpp#L922-L930
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/SetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddValue
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/ImplicitScheme
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_output/CadJsonOutput/GetParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_output/CadJsonOutput/GetBrepSurfaceParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_output/CadJsonOutput/GetBoundaryLoopParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/linear_solvers/amgcl_solver_impl/AMGCLSolver<TSparse,TDense>/ApplySettings
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/voxel_mesh_generator_modeler/VoxelMeshGeneratorModeler/ApplyColoring
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/ValidateTheDefaultParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_kratos_parameters/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_kratos_parameters/TestParameters/test_add_methods
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/profiler/Profiler/Write
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Parameters::AddInt( const std::string& rEntry, const int Value )`

# Calls

- [SetInt](../../../../functions/sources/kratos_parameters/Parameters/SetInt.md)
- [AddValue](../../../../functions/sources/kratos_parameters/Parameters/AddValue.md)

# Called by

- [ImplicitScheme](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/ImplicitScheme.md)
- [GetParameters](../../../../functions/input_output/cad_json_output/CadJsonOutput/GetParameters.md)
- [GetBrepSurfaceParameters](../../../../functions/input_output/cad_json_output/CadJsonOutput/GetBrepSurfaceParameters.md)
- [GetBoundaryLoopParameters](../../../../functions/input_output/cad_json_output/CadJsonOutput/GetBoundaryLoopParameters.md)
- [ApplySettings](../../../../functions/linear_solvers/amgcl_solver_impl/AMGCLSolver<TSparse,TDense>/ApplySettings.md)
- [ApplyColoring](../../../../functions/modeler/voxel_mesh_generator_modeler/VoxelMeshGeneratorModeler/ApplyColoring.md)
- [ValidateTheDefaultParameters](../../../../functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/ValidateTheDefaultParameters.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_kratos_parameters/KRATOS_TEST_CASE_IN_SUITE.md)
- [test_add_methods](../../../../functions/tests/test_kratos_parameters/TestParameters/test_add_methods.md)
- [Write](../../../../functions/utilities/profiler/Profiler/Write.md)