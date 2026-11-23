---
type: C++ Method
title: CreateNewMasterSlaveConstraint
resource: sources/model_part.cpp#L1239-L1301
visibility: private
generated:
  by: okf-rs/0.7.0
---

# Signature

`ModelPart::MasterSlaveConstraintType::Pointer ModelPart::CreateNewMasterSlaveConstraint(const std::string& ConstraintName, ModelPart::IndexType Id, ModelPart::NodeType& rMasterNode, const ModelPart::DoubleVariableType& rMasterVariable, ModelPart::NodeType& rSlaveNode, const ModelPart::DoubleVariableType& rSlaveVariable, const double Weight, const double Constant, IndexType ThisIndex)`