---
type: C++ Module
title: vtu_output
resource: input_output/vtu_output.cpp#L1-L1656
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/algorithm
    resolved_by: tree-sitter
    confidence: exact
  - target: external/filesystem
    resolved_by: tree-sitter
    confidence: exact
  - target: external/numeric
    resolved_by: tree-sitter
    confidence: exact
  - target: external/type-traits
    resolved_by: tree-sitter
    confidence: exact
  - target: external/cstdint
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-data-communicator-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/input-output-base-64-encoded-output-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/input-output-vtk-definitions-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/tensor-adaptors-flags-tensor-adaptor-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/tensor-adaptors-historical-variable-tensor-adaptor-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/tensor-adaptors-node-position-tensor-adaptor-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/tensor-adaptors-variable-tensor-adaptor-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-container-io-utils-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-data-type-traits-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-parallel-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-reduction-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-string-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-xml-utilities-xml-appended-data-element-wrapper-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-xml-utilities-xml-elements-array-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-xml-utilities-xml-in-place-data-element-wrapper-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-xml-utilities-xml-utils-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/vtu-output-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [IsBigEndian](../../functions/input_output/vtu_output/IsBigEndian.md)
- [GetBaseName](../../functions/input_output/vtu_output/GetBaseName.md)
- [GetFileName](../../functions/input_output/vtu_output/GetFileName.md)
- [GetEndianness](../../functions/input_output/vtu_output/GetEndianness.md)
- [CheckDataArrayName](../../functions/input_output/vtu_output/CheckDataArrayName.md)
- [CheckDataArrayName](../../functions/input_output/vtu_output/CheckDataArrayName-2.md)
- [CheckDataArrayName](../../functions/input_output/vtu_output/CheckDataArrayName-3.md)
- [for](../../functions/input_output/vtu_output/for.md)
- [GetEntityName](../../functions/input_output/vtu_output/GetEntityName.md)
- [GetWritingNDData](../../functions/input_output/vtu_output/GetWritingNDData.md)
- [GetGeometryTypes](../../functions/input_output/vtu_output/GetGeometryTypes.md)
- [GetOffsets](../../functions/input_output/vtu_output/GetOffsets.md)
- [GetConnectivities](../../functions/input_output/vtu_output/GetConnectivities.md)
- [AddConnectivityData](../../functions/input_output/vtu_output/AddConnectivityData.md)
- [AddFieldsUsingTensorAdaptorImpl](../../functions/input_output/vtu_output/AddFieldsUsingTensorAdaptorImpl.md)
- [AddFieldsUsingTensorAdaptor](../../functions/input_output/vtu_output/AddFieldsUsingTensorAdaptor.md)
- [for](../../functions/input_output/vtu_output/for-2.md)
- [AddFields](../../functions/input_output/vtu_output/AddFields.md)
- [GetNodesContainer](../../functions/input_output/vtu_output/GetNodesContainer.md)
- [AddUnstructuredGridData](../../functions/input_output/vtu_output/AddUnstructuredGridData.md)
- [PrintLocationData](../../functions/input_output/vtu_output/PrintLocationData.md)
- [FindUnstructuredGridData](../../functions/input_output/vtu_output/FindUnstructuredGridData.md)
- [VtuOutput](../../functions/input_output/vtu_output/VtuOutput/VtuOutput.md)
- [GetOutputContainerList](../../functions/input_output/vtu_output/VtuOutput/GetOutputContainerList.md)
- [AddFlag](../../functions/input_output/vtu_output/VtuOutput/AddFlag.md)
- [AddVariable](../../functions/input_output/vtu_output/VtuOutput/AddVariable.md)
- [AddIntegrationPointVariable](../../functions/input_output/vtu_output/VtuOutput/AddIntegrationPointVariable.md)
- [AddTensorAdaptor](../../functions/input_output/vtu_output/VtuOutput/AddTensorAdaptor.md)
- [ReplaceTensorAdaptor](../../functions/input_output/vtu_output/VtuOutput/ReplaceTensorAdaptor.md)
- [EmplaceTensorAdaptor](../../functions/input_output/vtu_output/VtuOutput/EmplaceTensorAdaptor.md)
- [WriteUnstructuredGridData](../../functions/input_output/vtu_output/VtuOutput/WriteUnstructuredGridData.md)
- [WriteIntegrationPointData](../../functions/input_output/vtu_output/VtuOutput/WriteIntegrationPointData.md)
- [WriteData](../../functions/input_output/vtu_output/VtuOutput/WriteData.md)
- [PrintOutput](../../functions/input_output/vtu_output/VtuOutput/PrintOutput.md)
- [PrintOutput](../../functions/input_output/vtu_output/VtuOutput/PrintOutput-2.md)
- [Info](../../functions/input_output/vtu_output/VtuOutput/Info.md)
- [PrintInfo](../../functions/input_output/vtu_output/VtuOutput/PrintInfo.md)
- [PrintData](../../functions/input_output/vtu_output/VtuOutput/PrintData.md)

# Imports

- `<algorithm>`
- `<filesystem>`
- `<numeric>`
- `<type_traits>`
- `<cstdint>`
- `"includes/data_communicator.h"`
- `"input_output/base_64_encoded_output.h"`
- `"input_output/vtk_definitions.h"`
- `"tensor_adaptors/flags_tensor_adaptor.h"`
- `"tensor_adaptors/historical_variable_tensor_adaptor.h"`
- `"tensor_adaptors/node_position_tensor_adaptor.h"`
- `"tensor_adaptors/variable_tensor_adaptor.h"`
- `"utilities/container_io_utils.h"`
- `"utilities/data_type_traits.h"`
- `"utilities/parallel_utilities.h"`
- `"utilities/reduction_utilities.h"`
- `"utilities/string_utilities.h"`
- `"utilities/xml_utilities/xml_appended_data_element_wrapper.h"`
- `"utilities/xml_utilities/xml_elements_array.h"`
- `"utilities/xml_utilities/xml_in_place_data_element_wrapper.h"`
- `"utilities/xml_utilities/xml_utils.h"`
- `"vtu_output.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)