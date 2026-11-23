---
type: C++ Function
title: VariableData
resource: containers/variable_data.h#L66-L72
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/variable_data/GenerateKey
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`constexpr VariableData(const VariableData& rOtherVariable) : mName(rOtherVariable.mName), mKey(rOtherVariable.mKey), mSize(rOtherVariable.mSize), mpSourceVariable(rOtherVariable.mpSourceVariable) , mIsComponent(rOtherVariable.mIsComponent)`

# Calls

- [GenerateKey](../../../functions/containers/variable_data/GenerateKey.md)