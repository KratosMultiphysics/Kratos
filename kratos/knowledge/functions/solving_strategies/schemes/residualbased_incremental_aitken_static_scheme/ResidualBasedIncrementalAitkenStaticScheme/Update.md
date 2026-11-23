---
type: C++ Method
title: Update
resource: solving_strategies/schemes/residualbased_incremental_aitken_static_scheme.h#L156-L204
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/dof/Dof/IsFree
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Update(ModelPart &r_model_part, DofsArrayType &rDofSet, TSystemMatrixType &A, TSystemVectorType &Dx, TSystemVectorType &b) override`

# Calls

- [IsFree](../../../../../functions/includes/dof/Dof/IsFree.md)