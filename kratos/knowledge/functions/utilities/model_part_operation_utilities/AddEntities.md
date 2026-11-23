---
type: C++ Function
title: AddEntities
resource: utilities/model_part_operation_utilities.h#L164-L196
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void AddEntities( std::vector<typename TContainerType::value_type*>& rEntitiesOutput, TContainerType& rMainEntityContainer, TUnaryFunc&& rIsValidEntity)`

# Calls

- [GetGeometry](../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)