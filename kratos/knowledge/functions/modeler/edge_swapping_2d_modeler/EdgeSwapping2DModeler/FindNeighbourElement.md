---
type: C++ Method
title: FindNeighbourElement
resource: modeler/edge_swapping_2d_modeler.h#L658-L692
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/edge_swapping_2d_modeler/EdgeSwapping2DModeler/SetSwappingData
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void FindNeighbourElement(unsigned int NodeId1, unsigned int NodeId2, int ElementId, ModelPart::ElementsContainerType::ContainerType const& ElementsArray, SwappingData& rSwappingData, int EdgeIndex)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [SetSwappingData](../../../../functions/modeler/edge_swapping_2d_modeler/EdgeSwapping2DModeler/SetSwappingData.md)