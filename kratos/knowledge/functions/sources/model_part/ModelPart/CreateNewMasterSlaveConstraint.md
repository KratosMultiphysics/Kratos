---
type: C++ Method
title: CreateNewMasterSlaveConstraint
resource: sources/model_part.cpp#L1189-L1237
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/model_part/IsSubModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/HasDofFor
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`ModelPart::MasterSlaveConstraintType::Pointer ModelPart::CreateNewMasterSlaveConstraint(const std::string& ConstraintName, IndexType Id, ModelPart::DofsVectorType& rMasterDofsVector, ModelPart::DofsVectorType& rSlaveDofsVector, const ModelPart::MatrixType& RelationMatrix, const ModelPart::VectorType& ConstantVector, IndexType ThisIndex)`

# Calls

- [IsSubModelPart](../../../../functions/includes/model_part/IsSubModelPart.md)
- [HasDofFor](../../../../functions/includes/node/HasDofFor.md)