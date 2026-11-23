---
type: C++ Method
title: FindNeighbourEntities
resource: utilities/model_part_operation_utilities.cpp#L99-L120
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/AddNeighbours
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartOperationUtilities::FindNeighbourEntities( std::vector<typename TContainerType::value_type*>& rOutputEntities, const Flags& rNodeSelectionFlag, TContainerType& rMainEntities)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [AddNeighbours](../../../../functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/AddNeighbours.md)