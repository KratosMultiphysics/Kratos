---
type: C++ Method
title: AllocateConstraints
resource: solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler.cpp#L152-L302
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/Matrix/fill
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/csr_matrix/CsrMatrix/nnz
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void MasterSlaveConstraintAssembler<TSparse,TDense>::AllocateConstraints(PointerVectorSet<MasterSlaveConstraint,IndexedObject>::const_iterator itConstraintBegin, PointerVectorSet<MasterSlaveConstraint,IndexedObject>::const_iterator itConstraintEnd, const ProcessInfo& rProcessInfo, typename Base::DofSet::const_iterator itDofBegin, typename Base::DofSet::const_iterator itDofEnd)`

# Calls

- [fill](../../../../../../functions/includes/amatrix_interface/Matrix/fill.md)
- [nnz](../../../../../../functions/containers/csr_matrix/CsrMatrix/nnz.md)