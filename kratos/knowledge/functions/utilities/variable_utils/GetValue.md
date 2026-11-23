---
type: C++ Function
title: GetValue
resource: utilities/variable_utils.h#L1575-L1601
generated:
  by: okf-rs/0.7.0
---

# Signature

`static std::conditional_t<std::is_arithmetic_v<TValue>, TValue, // <== return by value if scalar type const TValue&> // <== return by reference if non-scalar type GetValue(const TEntity& rEntity, const Variable<TValue>& rVariable)`