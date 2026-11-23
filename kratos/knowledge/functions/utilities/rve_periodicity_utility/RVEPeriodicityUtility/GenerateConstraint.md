---
type: C++ Method
title: GenerateConstraint
resource: utilities/rve_periodicity_utility.cpp#L173-L192
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/rve_periodicity_utility/RVEPeriodicityUtility/Finalize
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`MasterSlaveConstraint::Pointer RVEPeriodicityUtility::GenerateConstraint( IndexType& rConstraintId, const DoubleVariableType& rVar, NodeType::Pointer pSlaveNode, const std::vector<IndexType>& rMasterIds, const Matrix& rRelationMatrix, const Vector& rTranslationVector)`

# Called by

- [Finalize](../../../../functions/utilities/rve_periodicity_utility/RVEPeriodicityUtility/Finalize.md)