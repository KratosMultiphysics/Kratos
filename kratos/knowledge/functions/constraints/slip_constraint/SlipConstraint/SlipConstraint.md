---
type: C++ Method
title: SlipConstraint
resource: constraints/slip_constraint.cpp#L25-L39
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/constraints/slip_constraint/SlipConstraint/ConstructorHelper
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`SlipConstraint::SlipConstraint( IndexType Id, Dof<double>* pDofX, Dof<double>* pDofY, array_1d<double,3> NormalVector //note that a copy is made ) : BaseType(Id)`

# Calls

- [ConstructorHelper](../../../../functions/constraints/slip_constraint/SlipConstraint/ConstructorHelper.md)