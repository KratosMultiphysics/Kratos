---
type: C++ Function
title: VariableData
resource: containers/variable_data.h#L339-L347
generated:
  by: okf-rs/0.7.0
---

# Signature

`constexpr VariableData( const std::string_view& NewName, const std::size_t NewSize) : mName(NewName), mKey(GenerateKey(mName, NewSize, false, 0)), mSize(NewSize), mpSourceVariable(this), mIsComponent(false)`