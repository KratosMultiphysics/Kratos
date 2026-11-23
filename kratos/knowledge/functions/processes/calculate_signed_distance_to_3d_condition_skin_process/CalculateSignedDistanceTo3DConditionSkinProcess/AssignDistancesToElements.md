---
type: C++ Method
title: AssignDistancesToElements
resource: processes/calculate_signed_distance_to_3d_condition_skin_process.h#L1137-L1158
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalcNodalDistanceTo3DSkin
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AssignDistancesToElements(ModelPart::ElementsContainerType::iterator& i_fluid_element, array_1d<double,4> ElementalDistances)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [CalcNodalDistanceTo3DSkin](../../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalcNodalDistanceTo3DSkin.md)