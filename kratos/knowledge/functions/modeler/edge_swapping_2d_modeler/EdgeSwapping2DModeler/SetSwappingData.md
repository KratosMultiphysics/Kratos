---
type: C++ Method
title: SetSwappingData
resource: modeler/edge_swapping_2d_modeler.h#L629-L655
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/modeler/edge_swapping_2d_modeler/EdgeSwapping2DModeler/FindNodalNeighbours
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/edge_swapping_2d_modeler/EdgeSwapping2DModeler/FindNeighbourElement
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/edge_swapping_2d_modeler/EdgeSwapping2DModeler/IsInElementSphere
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/edge_swapping_2d_modeler/EdgeSwapping2DModeler/Remesh
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SetSwappingData(ModelPart& rThisModelPart)`

# Calls

- [FindNodalNeighbours](../../../../functions/modeler/edge_swapping_2d_modeler/EdgeSwapping2DModeler/FindNodalNeighbours.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [FindNeighbourElement](../../../../functions/modeler/edge_swapping_2d_modeler/EdgeSwapping2DModeler/FindNeighbourElement.md)
- [IsInElementSphere](../../../../functions/modeler/edge_swapping_2d_modeler/EdgeSwapping2DModeler/IsInElementSphere.md)

# Called by

- [Remesh](../../../../functions/modeler/edge_swapping_2d_modeler/EdgeSwapping2DModeler/Remesh.md)