---
type: C++ Method
title: RetrieveConditionsNeighbourConditionsIds
resource: utilities/auxiliar_model_part_utilities.cpp#L535-L575
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_find_conditions_neighbours_process/TestFindConditionsNeighboursProcess/test_sphere
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::unordered_map<IndexType, std::vector<IndexType>> AuxiliarModelPartUtilities::RetrieveConditionsNeighbourConditionsIds()`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [test_sphere](../../../../functions/tests/test_find_conditions_neighbours_process/TestFindConditionsNeighboursProcess/test_sphere.md)