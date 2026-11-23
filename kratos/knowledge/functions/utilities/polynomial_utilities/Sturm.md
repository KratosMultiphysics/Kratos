---
type: C++ Function
title: Sturm
resource: utilities/polynomial_utilities.cpp#L49-L56
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/polynomial_utilities/Divide
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/polynomial_utilities/IsolateRoots
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`PolynomialType Sturm(const PolynomialType& rS1, const PolynomialType& rS2)`

# Calls

- [Divide](../../../functions/utilities/polynomial_utilities/Divide.md)

# Called by

- [IsolateRoots](../../../functions/utilities/polynomial_utilities/IsolateRoots.md)