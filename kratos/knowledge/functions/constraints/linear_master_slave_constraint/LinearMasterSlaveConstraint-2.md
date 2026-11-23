---
type: C++ Function
title: LinearMasterSlaveConstraint
resource: constraints/linear_master_slave_constraint.h#L118-L130
generated:
  by: okf-rs/0.7.0
---

# Signature

`LinearMasterSlaveConstraint( IndexType Id, DofPointerVectorType& rMasterDofsVector, DofPointerVectorType& rSlaveDofsVector, const MatrixType& rRelationMatrix, const VectorType& rConstantVector ) : BaseType(Id), mSlaveDofsVector(rSlaveDofsVector), mMasterDofsVector(rMasterDofsVector), mRelationMatrix(rRelationMatrix), mConstantVector(rConstantVector)`