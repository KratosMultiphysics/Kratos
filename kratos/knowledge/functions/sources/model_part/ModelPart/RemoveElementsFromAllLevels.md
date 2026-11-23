---
type: C++ Method
title: RemoveElementsFromAllLevels
resource: sources/model_part.cpp#L1144-L1148
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part/ModelPart/RemoveElements
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/combine_model_part_modeler/CombineModelPartModeler/ResetModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/connectivity_preserve_modeler/ConnectivityPreserveModeler/ResetModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/create_entities_from_geometries_modeler/CreateEntitiesFromGeometriesModeler/RemoveModelPartEntities
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/entity_erase_process/EntitiesEraseProcess<Element>/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/ExecuteExplicitMapping
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/ExecuteImplicitMapping
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_model_part_to_python/ModelPartRemoveElementsFromAllLevels
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/processes/test_skin_detection_process/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/Clear
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPart::RemoveElementsFromAllLevels(Flags IdentifierFlag)`

# Calls

- [RemoveElements](../../../../functions/sources/model_part/ModelPart/RemoveElements.md)

# Called by

- [ResetModelPart](../../../../functions/modeler/combine_model_part_modeler/CombineModelPartModeler/ResetModelPart.md)
- [ResetModelPart](../../../../functions/modeler/connectivity_preserve_modeler/ConnectivityPreserveModeler/ResetModelPart.md)
- [RemoveModelPartEntities](../../../../functions/modeler/create_entities_from_geometries_modeler/CreateEntitiesFromGeometriesModeler/RemoveModelPartEntities.md)
- [Execute](../../../../functions/processes/entity_erase_process/EntitiesEraseProcess<Element>/Execute.md)
- [ExecuteExplicitMapping](../../../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/ExecuteExplicitMapping.md)
- [ExecuteImplicitMapping](../../../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/ExecuteImplicitMapping.md)
- [ModelPartRemoveElementsFromAllLevels](../../../../functions/python/add_model_part_to_python/ModelPartRemoveElementsFromAllLevels.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/processes/test_skin_detection_process/KRATOS_TEST_CASE_IN_SUITE.md)
- [Clear](../../../../functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/Clear.md)