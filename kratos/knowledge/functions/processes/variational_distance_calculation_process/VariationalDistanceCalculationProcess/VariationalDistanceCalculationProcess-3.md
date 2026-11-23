---
type: C++ Method
title: VariationalDistanceCalculationProcess
resource: processes/variational_distance_calculation_process.h#L189-L221
generated:
  by: okf-rs/0.7.0
---

# Signature

`VariationalDistanceCalculationProcess( ModelPart& rBaseModelPart, typename TLinearSolver::Pointer pLinearSolver, unsigned int MaxIterations = 10, Flags Options = CALCULATE_EXACT_DISTANCES_TO_PLANE.AsFalse(), std::string AuxPartName = "RedistanceCalculationPart", double Coefficient1 = 0.01, double Coefficient2 = 0.1) : mDistancePartIsInitialized(false), mMaxIterations(MaxIterations), mrModel( rBaseModelPart.GetModel() ), mrBaseModelPart (rBaseModelPart), mOptions( Options ), mAuxModelPartName( AuxPartName ), mCoefficient1(Coefficient1), mCoefficient2(Coefficient2)`