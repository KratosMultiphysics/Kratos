---
type: C++ Method
title: GenerateConditions
resource: modeler/structured_mesh_refinement_modeler.h#L341-L366
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/structured_mesh_refinement_modeler/StructuredMeshRefinementModeler/GenerateMesh
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void GenerateConditions(ModelPart& rThisModelPart, Condition& rThisCondition, NodesVectorType& VolumeNodesVector, NodesVectorType& SurfaceNodesVector, vector<int>& number_of_divisions, SizeType StartConditionId)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [GenerateMesh](../../../../functions/modeler/structured_mesh_refinement_modeler/StructuredMeshRefinementModeler/GenerateMesh.md)