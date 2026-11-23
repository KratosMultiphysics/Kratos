---
type: C++ Function
title: GetNodesContainer
resource: input_output/vtu_output.cpp#L564-L578
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/vtu_output/AddUnstructuredGridData
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`ModelPart::NodesContainerType::Pointer GetNodesContainer(TEntityContainerType& rContainer)`

# Calls

- [GetGeometry](../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [AddUnstructuredGridData](../../../functions/input_output/vtu_output/AddUnstructuredGridData.md)