---
type: C++ Function
title: Dof
resource: includes/dof.h#L119-L133
generated:
  by: okf-rs/0.7.0
---

# Signature

`Dof(NodalData* pThisNodalData, const TVariableType& rThisVariable) : mIsFixed(false), mVariableType(DofTrait<TDataType, TVariableType>::Id), mReactionType(DofTrait<TDataType, Variable<TDataType> >::Id), mEquationId(IndexType()), mEffectiveEquationId(IndexType()), mpNodalData(pThisNodalData)`