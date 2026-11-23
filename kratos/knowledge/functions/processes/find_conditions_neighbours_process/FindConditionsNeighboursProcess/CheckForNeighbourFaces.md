---
type: C++ Method
title: CheckForNeighbourFaces
resource: processes/find_conditions_neighbours_process.cpp#L145-L165
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/find_conditions_neighbours_process/FindConditionsNeighboursProcess/Execute
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`Condition::WeakPointer FindConditionsNeighboursProcess::CheckForNeighbourFaces( const unsigned int Id1, const unsigned int Id2, GlobalPointersVector<Condition>& rNeighbourFace, const unsigned int Face )`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [Execute](../../../../functions/processes/find_conditions_neighbours_process/FindConditionsNeighboursProcess/Execute.md)