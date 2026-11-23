---
type: C++ Method
title: InternalApplyValue
resource: processes/apply_constant_vectorvalue_process.cpp#L165-L192
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/node/HasDofFor
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/Fix
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ApplyConstantVectorValueProcess::InternalApplyValue( const Variable<double>& rVariable, const bool ToBeFixed, const double Value )`

# Calls

- [HasDofFor](../../../../functions/includes/node/HasDofFor.md)
- [Fix](../../../../functions/includes/node/Fix.md)