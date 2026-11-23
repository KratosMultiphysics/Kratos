---
type: C++ Function
title: FullName
resource: includes/model_part.h#L1830-L1837
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/model_part/IsSubModelPart
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/containers/edge_based_data_structure/EdgeBasedDataStructure/CalculateEdgeDataValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/AddUnstructuredGridData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/VtuOutput/VtuOutput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/VtuOutput/WriteUnstructuredGridData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/VtuOutput/WriteIntegrationPointData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/VtuOutput/PrintOutput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/VtuOutput/Info
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/VtuOutput/PrintData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/combine_model_part_modeler/CombineModelPartModeler/CheckOriginModelPartsAndAssignRoot
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/create_entities_from_geometries_modeler/CreateEntitiesFromGeometriesModeler/RemoveModelPartEntities
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/operation/find_contacts_in_skin_model_part/FindContactsInSkinModelPart/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/future/tests/test_mpi_vtu_output/check_pvtu
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/future/tests/test_mpi_vtu_output/check
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/utilities/distributed_model_part_initializer/GetSubModelPartHierarchy
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/compute_nodal_gradient_process/ComputeNodalGradientProcess<TOutputHistorical>/ComputeElementalContributionsAndVolume
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/find_global_nodal_entity_neighbours_process/FindGlobalNodalEntityNeighboursProcess<ModelPart/ConditionsContainerType>/FindGlobalNodalEntityNeighboursProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/find_global_nodal_entity_neighbours_process/FindGlobalNodalEntityNeighboursProcess<ModelPart/ElementsContainerType>/FindGlobalNodalEntityNeighboursProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/find_global_nodal_entity_neighbours_process/FindGlobalNodalEntityNeighboursProcess<TContainerType>/FindGlobalNodalEntityNeighboursProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_hdf5_postprocess/GiDHDF5OutputProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_output_process/GiDOutputProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/modelers/delete_model_parts_modeler/DeleteModelPartsModeler/__GetSubModelPartsNames
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/vtu_output_process/VtuOutputProcess/PrintOutput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/AddEntitiesFromIds
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/AddGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/RemoveSubModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/ErrorNonExistingSubModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_model_part/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtu_output/TestVtuOutput/test_PointVariableAddition
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtu_output/TestVtuOutput/GetUnstructuredGridList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtu_output/TestVtuOutput/GetUnstructuredGridName
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtu_output/TestVtuOutput/CheckGaussVtuFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/dof_utilities/block_build_dof_array_utility/BlockBuildDofArrayUtility/SetUpDofArray
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/dof_utilities/dof_array_utilities/DofArrayUtilities/SetUpDofArray
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/dof_utilities/elimination_build_dof_array_utility/EliminationBuildDofArrayUtility/SetUpDofArray
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/file_name_data_collector/FileNameDataCollector/FileNameDataCollector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/FillOutputSubModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/CheckValidityOfModelPartsForOperations
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/model_part_operation_utilities/FillSubModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/rve_periodicity_utility/RVEPeriodicityUtility/AssignPeriodicity
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/CalculateGradientBasedConformingExtensionBasis
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/CalculateMeshlessBasedConformingExtensionBasis
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/CalculateMeshlessBasedNonConformingExtensionBasis
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/variable_utils/CopyModelPartFlaggedNodalHistoricalVarToHistoricalVar
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/variable_utils/CopyModelPartFlaggedNodalNonHistoricalVarToNonHistoricalVar
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/variable_utils/CopyModelPartFlaggedElementVar
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/variable_utils/CopyModelPartFlaggedConditionVar
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::string FullName() const`

# Calls

- [IsSubModelPart](../../../functions/includes/model_part/IsSubModelPart.md)

# Called by

- [CalculateEdgeDataValues](../../../functions/containers/edge_based_data_structure/EdgeBasedDataStructure/CalculateEdgeDataValues.md)
- [AddUnstructuredGridData](../../../functions/input_output/vtu_output/AddUnstructuredGridData.md)
- [VtuOutput](../../../functions/input_output/vtu_output/VtuOutput/VtuOutput.md)
- [WriteUnstructuredGridData](../../../functions/input_output/vtu_output/VtuOutput/WriteUnstructuredGridData.md)
- [WriteIntegrationPointData](../../../functions/input_output/vtu_output/VtuOutput/WriteIntegrationPointData.md)
- [PrintOutput](../../../functions/input_output/vtu_output/VtuOutput/PrintOutput.md)
- [Info](../../../functions/input_output/vtu_output/VtuOutput/Info.md)
- [PrintData](../../../functions/input_output/vtu_output/VtuOutput/PrintData.md)
- [CheckOriginModelPartsAndAssignRoot](../../../functions/modeler/combine_model_part_modeler/CombineModelPartModeler/CheckOriginModelPartsAndAssignRoot.md)
- [RemoveModelPartEntities](../../../functions/modeler/create_entities_from_geometries_modeler/CreateEntitiesFromGeometriesModeler/RemoveModelPartEntities.md)
- [Execute](../../../functions/modeler/operation/find_contacts_in_skin_model_part/FindContactsInSkinModelPart/Execute.md)
- [check_pvtu](../../../functions/mpi/future/tests/test_mpi_vtu_output/check_pvtu.md)
- [check](../../../functions/mpi/future/tests/test_mpi_vtu_output/check.md)
- [GetSubModelPartHierarchy](../../../functions/mpi/utilities/distributed_model_part_initializer/GetSubModelPartHierarchy.md)
- [ComputeElementalContributionsAndVolume](../../../functions/processes/compute_nodal_gradient_process/ComputeNodalGradientProcess<TOutputHistorical>/ComputeElementalContributionsAndVolume.md)
- [FindGlobalNodalEntityNeighboursProcess](../../../functions/processes/find_global_nodal_entity_neighbours_process/FindGlobalNodalEntityNeighboursProcess<ModelPart/ConditionsContainerType>/FindGlobalNodalEntityNeighboursProcess.md)
- [FindGlobalNodalEntityNeighboursProcess](../../../functions/processes/find_global_nodal_entity_neighbours_process/FindGlobalNodalEntityNeighboursProcess<ModelPart/ElementsContainerType>/FindGlobalNodalEntityNeighboursProcess.md)
- [FindGlobalNodalEntityNeighboursProcess](../../../functions/processes/find_global_nodal_entity_neighbours_process/FindGlobalNodalEntityNeighboursProcess<TContainerType>/FindGlobalNodalEntityNeighboursProcess.md)
- [__init__](../../../functions/python_scripts/gid_hdf5_postprocess/GiDHDF5OutputProcess/__init__.md)
- [__init__](../../../functions/python_scripts/gid_output_process/GiDOutputProcess/__init__.md)
- [__GetSubModelPartsNames](../../../functions/python_scripts/modelers/delete_model_parts_modeler/DeleteModelPartsModeler/__GetSubModelPartsNames.md)
- [PrintOutput](../../../functions/python_scripts/vtu_output_process/VtuOutputProcess/PrintOutput.md)
- [AddEntitiesFromIds](../../../functions/sources/model_part/AddEntitiesFromIds.md)
- [AddGeometry](../../../functions/sources/model_part/ModelPart/AddGeometry.md)
- [RemoveSubModelPart](../../../functions/sources/model_part/ModelPart/RemoveSubModelPart.md)
- [ErrorNonExistingSubModelPart](../../../functions/sources/model_part/ModelPart/ErrorNonExistingSubModelPart.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/sources/test_model_part/KRATOS_TEST_CASE_IN_SUITE.md)
- [test_PointVariableAddition](../../../functions/tests/test_vtu_output/TestVtuOutput/test_PointVariableAddition.md)
- [GetUnstructuredGridList](../../../functions/tests/test_vtu_output/TestVtuOutput/GetUnstructuredGridList.md)
- [GetUnstructuredGridName](../../../functions/tests/test_vtu_output/TestVtuOutput/GetUnstructuredGridName.md)
- [CheckGaussVtuFile](../../../functions/tests/test_vtu_output/TestVtuOutput/CheckGaussVtuFile.md)
- [SetUpDofArray](../../../functions/utilities/dof_utilities/block_build_dof_array_utility/BlockBuildDofArrayUtility/SetUpDofArray.md)
- [SetUpDofArray](../../../functions/utilities/dof_utilities/dof_array_utilities/DofArrayUtilities/SetUpDofArray.md)
- [SetUpDofArray](../../../functions/utilities/dof_utilities/elimination_build_dof_array_utility/EliminationBuildDofArrayUtility/SetUpDofArray.md)
- [FileNameDataCollector](../../../functions/utilities/file_name_data_collector/FileNameDataCollector/FileNameDataCollector.md)
- [FillOutputSubModelPart](../../../functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/FillOutputSubModelPart.md)
- [CheckValidityOfModelPartsForOperations](../../../functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/CheckValidityOfModelPartsForOperations.md)
- [FillSubModelPart](../../../functions/utilities/model_part_operation_utilities/FillSubModelPart.md)
- [AssignPeriodicity](../../../functions/utilities/rve_periodicity_utility/RVEPeriodicityUtility/AssignPeriodicity.md)
- [CalculateGradientBasedConformingExtensionBasis](../../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/CalculateGradientBasedConformingExtensionBasis.md)
- [CalculateMeshlessBasedConformingExtensionBasis](../../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/CalculateMeshlessBasedConformingExtensionBasis.md)
- [CalculateMeshlessBasedNonConformingExtensionBasis](../../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/CalculateMeshlessBasedNonConformingExtensionBasis.md)
- [CopyModelPartFlaggedNodalHistoricalVarToHistoricalVar](../../../functions/utilities/variable_utils/CopyModelPartFlaggedNodalHistoricalVarToHistoricalVar.md)
- [CopyModelPartFlaggedNodalNonHistoricalVarToNonHistoricalVar](../../../functions/utilities/variable_utils/CopyModelPartFlaggedNodalNonHistoricalVarToNonHistoricalVar.md)
- [CopyModelPartFlaggedElementVar](../../../functions/utilities/variable_utils/CopyModelPartFlaggedElementVar.md)
- [CopyModelPartFlaggedConditionVar](../../../functions/utilities/variable_utils/CopyModelPartFlaggedConditionVar.md)