---
type: C++ Function
title: Dof
resource: includes/dof.h#L157-L176
generated:
  by: okf-rs/0.7.0
---

# Signature

`Dof(NodalData* pThisNodalData, const TVariableType& rThisVariable, const TReactionType& rThisReaction) : mIsFixed(false), mVariableType(DofTrait<TDataType, TVariableType>::Id), mReactionType(DofTrait<TDataType, TReactionType>::Id), mEquationId(IndexType()), mEffectiveEquationId(IndexType()), mpNodalData(pThisNodalData)`