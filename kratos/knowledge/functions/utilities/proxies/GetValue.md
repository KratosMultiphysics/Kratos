---
type: C++ Function
title: GetValue
resource: utilities/proxies.h#L109-L115
generated:
  by: okf-rs/0.7.0
---

# Signature

`std::conditional_t<std::is_integral_v<TValue> || std::is_floating_point_v<TValue>, TValue, // <== return by value if scalar type const TValue&> // <== return by reference in non-scalar type GetValue(const Variable<TValue>& rVariable) const`