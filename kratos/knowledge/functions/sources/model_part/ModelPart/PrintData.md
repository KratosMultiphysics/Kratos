---
type: C++ Method
title: PrintData
resource: sources/model_part.cpp#L2245-L2281
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/model_part/IsSubModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/model_part/NumberOfTables
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/model_part/NumberOfSubModelParts
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/communicator/Communicator/TotalProcesses
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPart::PrintData(std::ostream& rOStream) const`

# Calls

- [IsSubModelPart](../../../../functions/includes/model_part/IsSubModelPart.md)
- [NumberOfTables](../../../../functions/includes/model_part/NumberOfTables.md)
- [NumberOfSubModelParts](../../../../functions/includes/model_part/NumberOfSubModelParts.md)
- [TotalProcesses](../../../../functions/sources/communicator/Communicator/TotalProcesses.md)