---
type: C++ Function
title: as_integer
resource: input_output/unv_output.h#L109-L113
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/input_output/unv_output/UnvOutput/WriteNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/unv_output/UnvOutput/WriteElements
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/unv_output/WriteNodalResultRecords
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`auto as_integer(Enumeration const value) -> typename std::underlying_type<Enumeration>::type`

# Called by

- [WriteNodes](../../../functions/input_output/unv_output/UnvOutput/WriteNodes.md)
- [WriteElements](../../../functions/input_output/unv_output/UnvOutput/WriteElements.md)
- [WriteNodalResultRecords](../../../functions/input_output/unv_output/WriteNodalResultRecords.md)