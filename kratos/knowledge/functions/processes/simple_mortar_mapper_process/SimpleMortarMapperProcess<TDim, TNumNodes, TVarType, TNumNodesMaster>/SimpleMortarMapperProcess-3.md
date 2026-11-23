---
type: C++ Method
title: SimpleMortarMapperProcess
resource: processes/simple_mortar_mapper_process.cpp#L115-L145
visibility: private
generated:
  by: okf-rs/0.7.0
---

# Signature

`SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>::SimpleMortarMapperProcess( ModelPart& rOriginModelPart, ModelPart& rDestinationModelPart, TVarType& rOriginVariable, TVarType& rDestinationVariable, Parameters ThisParameters, LinearSolverType::Pointer pThisLinearSolver ): mOriginModelPart(rOriginModelPart), mDestinationModelPart(rDestinationModelPart), mpOriginVariable(&rOriginVariable), mpDestinationVariable(&rDestinationVariable), mThisParameters(ThisParameters), mpThisLinearSolver(pThisLinearSolver)`