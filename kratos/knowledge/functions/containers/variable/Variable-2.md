---
type: C++ Function
title: Variable
resource: containers/variable.h#L111-L118
generated:
  by: okf-rs/0.7.0
---

# Signature

`constexpr explicit Variable( const std::string_view& NewName, TSourceVariableType* pSourceVariable, char ComponentIndex, const VariableType* pTimeDerivativeVariable) : VariableData(NewName, sizeof(TDataType), pSourceVariable, ComponentIndex), mpTimeDerivativeVariable(pTimeDerivativeVariable)`