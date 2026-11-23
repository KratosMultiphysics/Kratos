---
type: C++ Function
title: GenerateNode
resource: modeler/structured_mesh_refinement_modeler.h#L240-L303
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/variables_list_data_value_container/TotalDataSize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/dof/Dof/FixDof
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void GenerateNode(ComponentType& rThisComponent, NodesVectorType& VolumeNodesVector, NodesVectorType& SurfaceNodesVector, Element::GeometryType::CoordinatesArrayType& rLocalCoordinates)`

# Calls

- [GetGeometry](../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [TotalDataSize](../../../functions/containers/variables_list_data_value_container/TotalDataSize.md)
- [GetVariable](../../../functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable.md)
- [FixDof](../../../functions/includes/dof/Dof/FixDof.md)