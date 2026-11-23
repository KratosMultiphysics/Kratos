---
type: C++ Method
title: SimpleMortarMapperProcess
resource: processes/simple_mortar_mapper_process.cpp#L80-L109
visibility: private
generated:
  by: okf-rs/0.7.0
---

# Signature

`SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>::SimpleMortarMapperProcess( ModelPart& rOriginModelPart, ModelPart& rDestinationModelPart, const TVarType& rThisVariable, Parameters ThisParameters, LinearSolverType::Pointer pThisLinearSolver ): mOriginModelPart(rOriginModelPart), mDestinationModelPart(rDestinationModelPart), mpOriginVariable(&rThisVariable), mpDestinationVariable(&rThisVariable), mThisParameters(ThisParameters), mpThisLinearSolver(pThisLinearSolver)`