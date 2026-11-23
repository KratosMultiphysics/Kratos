---
type: C++ Method
title: PrintData
resource: input_output/vtu_output.cpp#L1591-L1622
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/vtu_output/PrintLocationData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/model_part/FullName
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/variable_indexing_python/visit
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void VtuOutput::PrintData(std::ostream& rOStream) const`

# Calls

- [PrintLocationData](../../../../functions/input_output/vtu_output/PrintLocationData.md)
- [FullName](../../../../functions/includes/model_part/FullName.md)
- [visit](../../../../functions/python/variable_indexing_python/visit.md)