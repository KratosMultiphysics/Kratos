---
type: C++ Method
title: GetCurrentControlValue
resource: controllers/output_controller.cpp#L120-L126
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python/variable_indexing_python/visit
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/unv_output_process/UnvOutputProcess/PrintOutput
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::variant<int, double> OutputController::GetCurrentControlValue() const`

# Calls

- [visit](../../../../functions/python/variable_indexing_python/visit.md)

# Called by

- [PrintOutput](../../../../functions/python_scripts/unv_output_process/UnvOutputProcess/PrintOutput.md)