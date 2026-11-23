---
type: Python Method
title: GetModel
resource: python_scripts/project.py#L36-L39
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/input_output/vtk_output/VtkOutput/WriteModelPartWithoutNodesToFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/copy_properties_modeler/CopyPropertiesModeler/CopyPropertiesModeler
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/Clear
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/GenerateIntersectedEdgesElementsModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/SetObtainedEmbeddedNodalValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/SetLinearStrategy
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/check_same_modelpart_using_skin_distance_process/CheckSameModelPartUsingSkinDistanceProcess<TDim>/Execute
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
  - target: functions/processes/levelset_convection_process/LevelSetConvectionProcess/LevelSetConvectionProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/VariationalDistanceCalculationProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_model_part_to_python/AddModelPartToPython
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_hdf5_postprocess/GiDHDF5OutputProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_output_process/GiDOutputProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_output_process/GiDOutputProcess/ExecuteInitialize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/orchestrators/orchestrator/Orchestrator/CreateStage
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/orchestrators/orchestrator/Orchestrator/__CreateListOfModelers
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/orchestrators/orchestrator/Orchestrator/__CreateListOfOperations
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/containers/test_model/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_combine_model_part_modeler/TestCombineModelPartModeler/test_CombineModelParts
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_obj_io/debug_vtk
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_processes/SolutionLoopPointOutputProcesses
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/TestIODebug
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def GetModel(self) -> KratosMultiphysics.Model:`

# Called by

- [WriteModelPartWithoutNodesToFile](../../../../functions/input_output/vtk_output/VtkOutput/WriteModelPartWithoutNodesToFile.md)
- [CopyPropertiesModeler](../../../../functions/modeler/copy_properties_modeler/CopyPropertiesModeler/CopyPropertiesModeler.md)
- [Clear](../../../../functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/Clear.md)
- [GenerateIntersectedEdgesElementsModelPart](../../../../functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/GenerateIntersectedEdgesElementsModelPart.md)
- [SetObtainedEmbeddedNodalValues](../../../../functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/SetObtainedEmbeddedNodalValues.md)
- [SetLinearStrategy](../../../../functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/SetLinearStrategy.md)
- [Execute](../../../../functions/processes/check_same_modelpart_using_skin_distance_process/CheckSameModelPartUsingSkinDistanceProcess<TDim>/Execute.md)
- [FindGlobalNodalEntityNeighboursProcess](../../../../functions/processes/find_global_nodal_entity_neighbours_process/FindGlobalNodalEntityNeighboursProcess<ModelPart/ConditionsContainerType>/FindGlobalNodalEntityNeighboursProcess.md)
- [FindGlobalNodalEntityNeighboursProcess](../../../../functions/processes/find_global_nodal_entity_neighbours_process/FindGlobalNodalEntityNeighboursProcess<ModelPart/ElementsContainerType>/FindGlobalNodalEntityNeighboursProcess.md)
- [FindGlobalNodalEntityNeighboursProcess](../../../../functions/processes/find_global_nodal_entity_neighbours_process/FindGlobalNodalEntityNeighboursProcess<TContainerType>/FindGlobalNodalEntityNeighboursProcess.md)
- [LevelSetConvectionProcess](../../../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/LevelSetConvectionProcess.md)
- [VariationalDistanceCalculationProcess](../../../../functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/VariationalDistanceCalculationProcess.md)
- [AddModelPartToPython](../../../../functions/python/add_model_part_to_python/AddModelPartToPython.md)
- [__init__](../../../../functions/python_scripts/gid_hdf5_postprocess/GiDHDF5OutputProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/gid_output_process/GiDOutputProcess/__init__.md)
- [ExecuteInitialize](../../../../functions/python_scripts/gid_output_process/GiDOutputProcess/ExecuteInitialize.md)
- [CreateStage](../../../../functions/python_scripts/orchestrators/orchestrator/Orchestrator/CreateStage.md)
- [__CreateListOfModelers](../../../../functions/python_scripts/orchestrators/orchestrator/Orchestrator/__CreateListOfModelers.md)
- [__CreateListOfOperations](../../../../functions/python_scripts/orchestrators/orchestrator/Orchestrator/__CreateListOfOperations.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/containers/test_model/KRATOS_TEST_CASE_IN_SUITE.md)
- [test_CombineModelParts](../../../../functions/tests/test_combine_model_part_modeler/TestCombineModelPartModeler/test_CombineModelParts.md)
- [debug_vtk](../../../../functions/tests/test_obj_io/debug_vtk.md)
- [SolutionLoopPointOutputProcesses](../../../../functions/tests/test_processes/SolutionLoopPointOutputProcesses.md)
- [TestIODebug](../../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/TestIODebug.md)