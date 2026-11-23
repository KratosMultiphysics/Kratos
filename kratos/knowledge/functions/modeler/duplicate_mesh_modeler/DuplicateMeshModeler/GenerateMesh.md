---
type: C++ Method
title: GenerateMesh
resource: modeler/duplicate_mesh_modeler.cpp#L24-L97
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/timer/Timer/Start
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/timer/Timer/Stop
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/duplicate_mesh_modeler/DuplicateMeshModeler/GenerateSubModelParts
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void DuplicateMeshModeler::GenerateMesh( ModelPart& rThisModelPart, Element const& rReferenceElement, Condition const& rReferenceCondition)`

# Calls

- [Start](../../../../functions/utilities/timer/Timer/Start.md)
- [Stop](../../../../functions/utilities/timer/Timer/Stop.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [GenerateSubModelParts](../../../../functions/modeler/duplicate_mesh_modeler/DuplicateMeshModeler/GenerateSubModelParts.md)