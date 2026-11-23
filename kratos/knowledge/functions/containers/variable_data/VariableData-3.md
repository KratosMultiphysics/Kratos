---
type: C++ Function
title: VariableData
resource: containers/variable_data.h#L350-L360
generated:
  by: okf-rs/0.7.0
---

# Signature

`constexpr VariableData( const std::string_view& NewName, const std::size_t NewSize, const VariableData* pSourceVariable, const char ComponentIndex) : mName(NewName), mKey(GenerateKey(pSourceVariable->mName, NewSize, true, ComponentIndex)), mSize(NewSize), mpSourceVariable(pSourceVariable), mIsComponent(true)`