---
type: C++ Function
title: IsDefined
resource: containers/flags.h#L325-L328
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/includes/node/IsActive
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/ensight_output/EnSightOutput/WriteNodalFlagToFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/ensight_output/EnSightOutput/WriteGeometricalFlagToFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteFlagContainerVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/cpp_tests/sources/test_flag_synchronization/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/cpp_tests/sources/test_mpi_communicator/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/apply_constant_scalarvalue_process/ApplyConstantScalarValueProcess/ApplyConstantScalarValueProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/apply_constant_vectorvalue_process/ApplyConstantVectorValueProcess/ApplyConstantVectorValueProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/geometrical_object/GeometricalObject/IsActive
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/master_slave_constraint/MasterSlaveConstraint/IsActive
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/flags_tensor_adaptor/FlagsTensorAdaptor/CollectData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/containers/test_flags/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/integration/test_integration_info/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_flags/TestFlags/testFlagSetReset
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_flags/TestFlags/testFlagFlip
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_flags/TestFlags/testFlagAndEqual
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_flags/TestFlags/testFlagOrEqual
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_tensor_adaptors/TestTensorAdaptors/__TestFlagsTensorAdaptor
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_variable_utils/TestVariableUtils/test_set_flag
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/TestIODebug
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/mortar_utilities/InvertNormalForFlag
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool IsDefined(Flags const & rOther) const`

# Called by

- [IsActive](../../../functions/includes/node/IsActive.md)
- [WriteNodalFlagToFile](../../../functions/input_output/ensight_output/EnSightOutput/WriteNodalFlagToFile.md)
- [WriteGeometricalFlagToFile](../../../functions/input_output/ensight_output/EnSightOutput/WriteGeometricalFlagToFile.md)
- [WriteFlagContainerVariable](../../../functions/input_output/vtk_output/VtkOutput/WriteFlagContainerVariable.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/mpi/tests/cpp_tests/sources/test_flag_synchronization/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/mpi/tests/cpp_tests/sources/test_mpi_communicator/KRATOS_TEST_CASE_IN_SUITE.md)
- [ApplyConstantScalarValueProcess](../../../functions/processes/apply_constant_scalarvalue_process/ApplyConstantScalarValueProcess/ApplyConstantScalarValueProcess.md)
- [ApplyConstantVectorValueProcess](../../../functions/processes/apply_constant_vectorvalue_process/ApplyConstantVectorValueProcess/ApplyConstantVectorValueProcess.md)
- [IsActive](../../../functions/sources/geometrical_object/GeometricalObject/IsActive.md)
- [IsActive](../../../functions/sources/master_slave_constraint/MasterSlaveConstraint/IsActive.md)
- [CollectData](../../../functions/tensor_adaptors/flags_tensor_adaptor/FlagsTensorAdaptor/CollectData.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/containers/test_flags/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/integration/test_integration_info/KRATOS_TEST_CASE_IN_SUITE.md)
- [testFlagSetReset](../../../functions/tests/test_flags/TestFlags/testFlagSetReset.md)
- [testFlagFlip](../../../functions/tests/test_flags/TestFlags/testFlagFlip.md)
- [testFlagAndEqual](../../../functions/tests/test_flags/TestFlags/testFlagAndEqual.md)
- [testFlagOrEqual](../../../functions/tests/test_flags/TestFlags/testFlagOrEqual.md)
- [__TestFlagsTensorAdaptor](../../../functions/tests/test_tensor_adaptors/TestTensorAdaptors/__TestFlagsTensorAdaptor.md)
- [test_set_flag](../../../functions/tests/test_variable_utils/TestVariableUtils/test_set_flag.md)
- [TestIODebug](../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/TestIODebug.md)
- [InvertNormalForFlag](../../../functions/utilities/mortar_utilities/InvertNormalForFlag.md)