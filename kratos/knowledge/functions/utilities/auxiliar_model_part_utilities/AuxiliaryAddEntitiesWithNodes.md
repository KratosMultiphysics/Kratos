---
type: C++ Function
title: AuxiliaryAddEntitiesWithNodes
resource: utilities/auxiliar_model_part_utilities.h#L878-L913
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/model_part/IsSubModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AuxiliaryAddEntitiesWithNodes( TEntitiesContainer& rEntitiesContainer, const std::vector<IndexType>& rEntitiesIds )`

# Calls

- [GetGeometry](../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [IsSubModelPart](../../../functions/includes/model_part/IsSubModelPart.md)