---
type: C++ Method
title: GetOutputContainerList
resource: input_output/vtu_output.cpp#L792-L821
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python/variable_indexing_python/visit
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_vtu_output/TestVtuOutput/test_GetOutputContainerList
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::vector<VtuOutput::SupportedContainerPointerType> VtuOutput::GetOutputContainerList() const`

# Calls

- [visit](../../../../functions/python/variable_indexing_python/visit.md)

# Called by

- [test_GetOutputContainerList](../../../../functions/tests/test_vtu_output/TestVtuOutput/test_GetOutputContainerList.md)