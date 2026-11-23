---
type: C++ Method
title: CalculateEmbeddedNodalVariableFromSkinProcess
resource: processes/calculate_embedded_nodal_variable_from_skin_process.h#L212-L264
generated:
  by: okf-rs/0.7.0
---

# Signature

`CalculateEmbeddedNodalVariableFromSkinProcess( ModelPart &rBaseModelPart, ModelPart &rSkinModelPart, Parameters LinearSolverSettings, const Variable<TVarType> &rSkinVariable, const Variable<TVarType> &rEmbeddedNodalVariable, const double GradientPenaltyCoefficient = 0.0, const unsigned int BufferPosition = 0, const std::string& AuxPartName = "IntersectedElementsModelPart", const std::size_t EchoLevel = 0) : Process() , mEchoLevel(EchoLevel) , mBufferPosition(BufferPosition), mAuxModelPartName(AuxPartName), mGradientPenaltyCoefficient(GradientPenaltyCoefficient), mrBaseModelPart(rBaseModelPart), mrSkinModelPart(rSkinModelPart), mrSkinVariable(rSkinVariable), mrEmbeddedNodalVariable(rEmbeddedNodalVariable)`