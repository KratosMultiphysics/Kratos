---
type: C++ Method
title: RecursivelyHasProperties
resource: sources/model_part.cpp#L650-L665
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/model_part/IsSubModelPart
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/obj_io/ObjIO/ReadModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/stl_io/StlIO/ReadModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/entity_generation/voxel_mesher_entity_generation/VoxelMesherEntityGeneration/CreateAndGetProperty
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/skin_detection_process/SkinDetectionProcess<TDim>/CreateConditions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/sub_model_part_skin_detection_process/SubModelPartSkinDetectionProcess<TDim>/CreateConditions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_model_part_to_python/AddModelPartToPython
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_model_part/TestModelPart/test_model_part_properties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/read_materials_utility/ReadMaterialsUtility/GetPropertyBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/read_materials_utility/ReadMaterialsUtility/CreateSubProperties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/read_materials_utility/ReadMaterialsUtility/AssignPropertyBlock
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool ModelPart::RecursivelyHasProperties( IndexType PropertiesId, IndexType MeshIndex ) const`

# Calls

- [IsSubModelPart](../../../../functions/includes/model_part/IsSubModelPart.md)

# Called by

- [ReadModelPart](../../../../functions/input_output/obj_io/ObjIO/ReadModelPart.md)
- [ReadModelPart](../../../../functions/input_output/stl_io/StlIO/ReadModelPart.md)
- [CreateAndGetProperty](../../../../functions/modeler/entity_generation/voxel_mesher_entity_generation/VoxelMesherEntityGeneration/CreateAndGetProperty.md)
- [CreateConditions](../../../../functions/processes/skin_detection_process/SkinDetectionProcess<TDim>/CreateConditions.md)
- [CreateConditions](../../../../functions/processes/sub_model_part_skin_detection_process/SubModelPartSkinDetectionProcess<TDim>/CreateConditions.md)
- [AddModelPartToPython](../../../../functions/python/add_model_part_to_python/AddModelPartToPython.md)
- [test_model_part_properties](../../../../functions/tests/test_model_part/TestModelPart/test_model_part_properties.md)
- [GetPropertyBlock](../../../../functions/utilities/read_materials_utility/ReadMaterialsUtility/GetPropertyBlock.md)
- [CreateSubProperties](../../../../functions/utilities/read_materials_utility/ReadMaterialsUtility/CreateSubProperties.md)
- [AssignPropertyBlock](../../../../functions/utilities/read_materials_utility/ReadMaterialsUtility/AssignPropertyBlock.md)