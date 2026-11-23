---
type: C++ Function
title: SetValue
resource: utilities/proxies.h#L126-L132
generated:
  by: okf-rs/0.7.0
---

# Signature

`void SetValue(const Variable<TValue>& rVariable, std::conditional_t<std::is_integral_v<TValue> || std::is_floating_point_v<TValue>, TValue, /*pass scalar types by value*/ const TValue&> /*pass non-scalar types by reference*/ Value)`