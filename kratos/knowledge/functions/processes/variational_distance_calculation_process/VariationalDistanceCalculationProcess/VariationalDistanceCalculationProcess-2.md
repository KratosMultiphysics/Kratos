---
type: C++ Method
title: VariationalDistanceCalculationProcess
resource: processes/variational_distance_calculation_process.h#L157-L187
generated:
  by: okf-rs/0.7.0
---

# Signature

`VariationalDistanceCalculationProcess( Model& rModel, typename TLinearSolver::Pointer pLinearSolver, BuilderSolverPointerType pBuilderAndSolver, Parameters ThisParameters) : mDistancePartIsInitialized(false) , mrModel(rModel) , mrBaseModelPart(rModel.GetModelPart(ThisParameters["model_part_name"].GetString()))`