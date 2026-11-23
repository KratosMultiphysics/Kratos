---
type: C++ Method
title: AddBool
resource: sources/kratos_parameters.cpp#L935-L943
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/SetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddValue
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/cad_json_output/CadJsonOutput/GetBrepSurfaceParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_output/CadJsonOutput/GetBoundaryLoopParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_mpi_processes/TestMPIProcesses/testDiscontinuousDistanceProcessCutOnEdge2D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_mpi_processes/TestMPIProcesses/testDiscontinuousDistanceProcessCutThroughNode2D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/processes/test_calculate_discontinuous_distance_to_skin_process/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_kratos_parameters/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_kratos_parameters/TestParameters/test_add_methods
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Parameters::AddBool( const std::string& rEntry, const bool Value )`

# Calls

- [SetBool](../../../../functions/sources/kratos_parameters/Parameters/SetBool.md)
- [AddValue](../../../../functions/sources/kratos_parameters/Parameters/AddValue.md)

# Called by

- [GetBrepSurfaceParameters](../../../../functions/input_output/cad_json_output/CadJsonOutput/GetBrepSurfaceParameters.md)
- [GetBoundaryLoopParameters](../../../../functions/input_output/cad_json_output/CadJsonOutput/GetBoundaryLoopParameters.md)
- [testDiscontinuousDistanceProcessCutOnEdge2D](../../../../functions/mpi/tests/test_mpi_processes/TestMPIProcesses/testDiscontinuousDistanceProcessCutOnEdge2D.md)
- [testDiscontinuousDistanceProcessCutThroughNode2D](../../../../functions/mpi/tests/test_mpi_processes/TestMPIProcesses/testDiscontinuousDistanceProcessCutThroughNode2D.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/processes/test_calculate_discontinuous_distance_to_skin_process/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_kratos_parameters/KRATOS_TEST_CASE_IN_SUITE.md)
- [test_add_methods](../../../../functions/tests/test_kratos_parameters/TestParameters/test_add_methods.md)