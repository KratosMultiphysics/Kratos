---
type: C++ Method
title: ConstructorHelper
resource: constraints/slip_constraint.cpp#L59-L99
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/constraints/slip_constraint/SlipConstraint/SlipConstraint
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SlipConstraint::ConstructorHelper( DofPointerVectorType& rAllDofs, array_1d<double,3>& rNormalVector )`

# Calls

- [norm_2](../../../../functions/includes/amatrix_interface/norm_2.md)

# Called by

- [SlipConstraint](../../../../functions/constraints/slip_constraint/SlipConstraint/SlipConstraint.md)