---
type: C++ Method
title: AddString
resource: sources/kratos_parameters.cpp#L948-L956
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/SetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddValue
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/cad_json_output/CadJsonOutput/GetBrepSurfaceParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/linear_solvers/amgcl_solver_impl/AMGCLSolver<TSparse,TDense>/AMGCLSolver
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/linear_solvers/amgcl_solver_impl/AMGCLSolver<TSparse,TDense>/ApplySettings
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/copy_properties_modeler/CopyPropertiesModeler/CopyPropertiesModeler
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/serial_model_part_combinator_modeler/SerialModelPartCombinatorModeler/SetupModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_hdf5_postprocess/GiDHDF5OutputProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_output_process/GiDOutputProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/__PrepareOutputSettings
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/project/Project/__init__
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

`void Parameters::AddString( const std::string& rEntry, const std::string& rValue )`

# Calls

- [SetString](../../../../functions/sources/kratos_parameters/Parameters/SetString.md)
- [AddValue](../../../../functions/sources/kratos_parameters/Parameters/AddValue.md)

# Called by

- [GetBrepSurfaceParameters](../../../../functions/input_output/cad_json_output/CadJsonOutput/GetBrepSurfaceParameters.md)
- [AMGCLSolver](../../../../functions/linear_solvers/amgcl_solver_impl/AMGCLSolver<TSparse,TDense>/AMGCLSolver.md)
- [ApplySettings](../../../../functions/linear_solvers/amgcl_solver_impl/AMGCLSolver<TSparse,TDense>/ApplySettings.md)
- [CopyPropertiesModeler](../../../../functions/modeler/copy_properties_modeler/CopyPropertiesModeler/CopyPropertiesModeler.md)
- [SetupModelPart](../../../../functions/modeler/serial_model_part_combinator_modeler/SerialModelPartCombinatorModeler/SetupModelPart.md)
- [__init__](../../../../functions/python_scripts/gid_hdf5_postprocess/GiDHDF5OutputProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/gid_output_process/GiDOutputProcess/__init__.md)
- [__PrepareOutputSettings](../../../../functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/__PrepareOutputSettings.md)
- [__init__](../../../../functions/python_scripts/project/Project/__init__.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_kratos_parameters/KRATOS_TEST_CASE_IN_SUITE.md)
- [test_add_methods](../../../../functions/tests/test_kratos_parameters/TestParameters/test_add_methods.md)
- [Write](../../../../functions/utilities/profiler/Profiler/Write.md)