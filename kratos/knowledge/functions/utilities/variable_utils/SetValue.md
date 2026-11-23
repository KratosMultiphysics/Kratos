---
type: C++ Function
title: SetValue
resource: utilities/variable_utils.h#L1637-L1664
generated:
  by: okf-rs/0.7.0
---

# Signature

`static void SetValue(TEntity& rEntity, const Variable<TValue>& rVariable, std::conditional_t<std::is_arithmetic_v<TValue>, TValue, /*pass scalar types by value*/ const TValue&> /*pass non-scalar types by reference*/ Value)`