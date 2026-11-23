---
type: C++ Method
title: GetNextPositionSkippingWhiteSpaces
resource: sources/code_location.cpp#L148-L154
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/sources/code_location/CodeLocation/ReduceTemplateArgumentsToFirstN
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::size_t CodeLocation::GetNextPositionSkippingWhiteSpaces(std::string const& ThisString, std::size_t Position)`

# Called by

- [ReduceTemplateArgumentsToFirstN](../../../../functions/sources/code_location/CodeLocation/ReduceTemplateArgumentsToFirstN.md)