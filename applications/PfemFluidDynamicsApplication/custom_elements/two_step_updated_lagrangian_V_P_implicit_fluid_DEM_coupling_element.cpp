//
//   Project Name:        KratosFluidDynamicsApplication $
//   Last modified by:    $Author:               AFranci $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 0.0 $
//
//   Implementation of the Gauss-Seidel two step Updated Lagrangian Velocity-Pressure element
//     ( There is a ScalingConstant to multiply the mass balance equation for a number because i read it somewhere)
//

// System includes

// External includes

// Project includes
#include "custom_elements/two_step_updated_lagrangian_V_P_implicit_fluid_DEM_coupling_element.h"
#include "includes/cfd_variables.h"
#include <math.h>

namespace Kratos
{

/*
   * public TwoStepUpdatedLagrangianVPImplicitFluidDEMcouplingElement<TDim> functions
   */

template <unsigned int TDim>
Element::Pointer TwoStepUpdatedLagrangianVPImplicitFluidDEMcouplingElement<TDim>::Clone(IndexType NewId, NodesArrayType const &rThisNodes) const
{

  TwoStepUpdatedLagrangianVPImplicitFluidDEMcouplingElement NewElement(NewId, this->GetGeometry().Create(rThisNodes), this->pGetProperties());

  NewElement.SetData(this->GetData());
  NewElement.SetFlags(this->GetFlags());

  return Element::Pointer(new TwoStepUpdatedLagrangianVPImplicitFluidDEMcouplingElement(NewElement));
}

template <unsigned int TDim>
void TwoStepUpdatedLagrangianVPImplicitFluidDEMcouplingElement<TDim>::Initialize()
{
  KRATOS_TRY;
  KRATOS_CATCH("");
}

template <unsigned int TDim>
void TwoStepUpdatedLagrangianVPImplicitFluidDEMcouplingElement<TDim>::InitializeSolutionStep(ProcessInfo &rCurrentProcessInfo)
{
}

template <unsigned int TDim>
void TwoStepUpdatedLagrangianVPImplicitFluidDEMcouplingElement<TDim>::InitializeNonLinearIteration(ProcessInfo &rCurrentProcessInfo)
{
  KRATOS_TRY;
  KRATOS_CATCH("");
}

template <unsigned int TDim>
void TwoStepUpdatedLagrangianVPImplicitFluidDEMcouplingElement<TDim>::ComputeMaterialParameters(double &Density,
                                                                                                double &DeviatoricCoeff,
                                                                                                double &VolumetricCoeff,
                                                                                                ProcessInfo &currentProcessInfo,
                                                                                                ElementalVariables &rElementalVariables)
{
  double timeStep = currentProcessInfo[DELTA_TIME];

  Density = this->GetProperties()[DENSITY];
  double FluidBulkModulus = this->GetProperties()[BULK_MODULUS];
  double FluidYieldShear = this->GetProperties()[YIELD_SHEAR];
  double staticFrictionCoefficient = this->GetProperties()[STATIC_FRICTION];

  if (FluidBulkModulus == 0)
  {
    FluidBulkModulus = 1000000000.0;
  }
  VolumetricCoeff = FluidBulkModulus * timeStep;

  if (FluidYieldShear != 0)
  {
    // std::cout<<"For a Newtonian fluid I should not enter here"<<std::endl;
    DeviatoricCoeff = this->ComputeNonLinearViscosity(rElementalVariables.EquivalentStrainRate);
  }
  else if (staticFrictionCoefficient != 0)
  {
    DeviatoricCoeff = this->ComputePapanastasiouMuIrheologyViscosity(rElementalVariables);
  }
  else
  {
    // std::cout<<"For a Newtonian fluid I should  enter here"<<std::endl;
    DeviatoricCoeff = this->GetProperties()[DYNAMIC_VISCOSITY];
  }

  // this->ComputeMaterialParametersGranularGas(rElementalVariables,VolumetricCoeff,DeviatoricCoeff);
  // std::cout<<"Density "<<Density<<std::endl;
  // std::cout<<"FluidBulkModulus "<<FluidBulkModulus<<std::endl;
  // std::cout<<"staticFrictionCoefficient "<<staticFrictionCoefficient<<std::endl;
  // std::cout<<"DeviatoricCoeff "<<DeviatoricCoeff<<std::endl;

  this->mMaterialDeviatoricCoefficient = DeviatoricCoeff;
  this->mMaterialVolumetricCoefficient = VolumetricCoeff;
  this->mMaterialDensity = Density;

  // const SizeType NumNodes = this->GetGeometry().PointsNumber();
  // for (SizeType i = 0; i < NumNodes; ++i)
  //   {
  // 	this->GetGeometry()[i].FastGetSolutionStepValue(ADAPTIVE_EXPONENT)=VolumetricCoeff;
  // 	this->GetGeometry()[i].FastGetSolutionStepValue(ALPHA_PARAMETER)=DeviatoricCoeff;
  // 	this->GetGeometry()[i].FastGetSolutionStepValue(FLOW_INDEX)=rElementalVariables.EquivalentStrainRate;
  //   }
}

template <unsigned int TDim>
int TwoStepUpdatedLagrangianVPImplicitFluidDEMcouplingElement<TDim>::Check(const ProcessInfo &rCurrentProcessInfo)
{
  KRATOS_TRY;

  // Base class checks for positive Jacobian and Id > 0
  int ierr = Element::Check(rCurrentProcessInfo);
  if (ierr != 0)
    return ierr;

  // Check that all required variables have been registered
  if (VELOCITY.Key() == 0)
    KRATOS_THROW_ERROR(std::invalid_argument, "VELOCITY Key is 0. Check that the application was correctly registered.", "");
  if (ACCELERATION.Key() == 0)
    KRATOS_THROW_ERROR(std::invalid_argument, "ACCELERATION Key is 0. Check that the application was correctly registered.", "");
  if (PRESSURE.Key() == 0)
    KRATOS_THROW_ERROR(std::invalid_argument, "PRESSURE Key is 0. Check that the application was correctly registered.", "");
  if (BODY_FORCE.Key() == 0)
    KRATOS_THROW_ERROR(std::invalid_argument, "BODY_FORCE Key is 0. Check that the application was correctly registered.", "");
  if (DENSITY.Key() == 0)
    KRATOS_THROW_ERROR(std::invalid_argument, "DENSITY Key is 0. Check that the application was correctly registered.", "");
  if (DYNAMIC_VISCOSITY.Key() == 0)
    KRATOS_THROW_ERROR(std::invalid_argument, "DYNAMIC_VISCOSITY Key is 0. Check that the application was correctly registered.", "");
  if (DELTA_TIME.Key() == 0)
    KRATOS_THROW_ERROR(std::invalid_argument, "DELTA_TIME Key is 0. Check that the application was correctly registered.", "");

  // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
  for (unsigned int i = 0; i < this->GetGeometry().size(); ++i)
  {
    if (this->GetGeometry()[i].SolutionStepsDataHas(VELOCITY) == false)
      KRATOS_THROW_ERROR(std::invalid_argument, "missing VELOCITY variable on solution step data for node ", this->GetGeometry()[i].Id());
    if (this->GetGeometry()[i].SolutionStepsDataHas(PRESSURE) == false)
      KRATOS_THROW_ERROR(std::invalid_argument, "missing PRESSURE variable on solution step data for node ", this->GetGeometry()[i].Id());
    if (this->GetGeometry()[i].SolutionStepsDataHas(BODY_FORCE) == false)
      KRATOS_THROW_ERROR(std::invalid_argument, "missing BODY_FORCE variable on solution step data for node ", this->GetGeometry()[i].Id());
    if (this->GetGeometry()[i].SolutionStepsDataHas(DENSITY) == false)
      KRATOS_THROW_ERROR(std::invalid_argument, "missing DENSITY variable on solution step data for node ", this->GetGeometry()[i].Id());
    if (this->GetGeometry()[i].SolutionStepsDataHas(DYNAMIC_VISCOSITY) == false)
      KRATOS_THROW_ERROR(std::invalid_argument, "missing DYNAMIC_VISCOSITY variable on solution step data for node ", this->GetGeometry()[i].Id());
    if (this->GetGeometry()[i].HasDofFor(VELOCITY_X) == false ||
        this->GetGeometry()[i].HasDofFor(VELOCITY_Y) == false ||
        this->GetGeometry()[i].HasDofFor(VELOCITY_Z) == false)
      KRATOS_THROW_ERROR(std::invalid_argument, "missing VELOCITY component degree of freedom on node ", this->GetGeometry()[i].Id());
    if (this->GetGeometry()[i].HasDofFor(PRESSURE) == false)
      KRATOS_THROW_ERROR(std::invalid_argument, "missing PRESSURE component degree of freedom on node ", this->GetGeometry()[i].Id());
  }

  // If this is a 2D problem, check that nodes are in XY plane
  if (this->GetGeometry().WorkingSpaceDimension() == 2)
  {
    for (unsigned int i = 0; i < this->GetGeometry().size(); ++i)
    {
      if (this->GetGeometry()[i].Z() != 0.0)
        KRATOS_THROW_ERROR(std::invalid_argument, "Node with non-zero Z coordinate found. Id: ", this->GetGeometry()[i].Id());
    }
  }

  return ierr;

  KRATOS_CATCH("");
}

template <>
void TwoStepUpdatedLagrangianVPImplicitFluidDEMcouplingElement<2>::ComputeBoundRHSVectorComplete(VectorType &BoundRHSVector,
                                                                                                 const double TimeStep,
                                                                                                 const double BoundRHSCoeffAcc,
                                                                                                 const double BoundRHSCoeffDev,
                                                                                                 const VectorType SpatialDefRate)
{
  GeometryType &rGeom = this->GetGeometry();

  if (rGeom[0].Is(FREE_SURFACE) && rGeom[1].Is(FREE_SURFACE))
  {
    array_1d<double, 3> AccA(3, 0.0);
    array_1d<double, 3> AccB(3, 0.0);
    array_1d<double, 3> MeanAcc(3, 0.0);
    array_1d<double, 3> NormalVector(3, 0.0);
    const double factor = 0.5 / TimeStep;
    const double one_third = 1.0 / 3.0;

    this->GetOutwardsUnitNormalForTwoPoints(NormalVector, 0, 1, 2);

    double SpatialDefRateNormalProjection = this->CalcNormalProjectionDefRate(SpatialDefRate, NormalVector);

    noalias(AccA) = factor * (rGeom[0].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[0].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION, 1);
    noalias(AccB) = factor * (rGeom[1].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[1].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION, 1);
    noalias(MeanAcc) = 0.5 * AccA + 0.5 * AccB;

    const double accelerationsNormalProjection = MeanAcc[0] * NormalVector[0] + MeanAcc[1] * NormalVector[1];

    if (rGeom[0].IsNot(INLET)) //to change into moving wall!!!!!
      BoundRHSVector[0] += one_third * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

    if (rGeom[1].IsNot(INLET))
      BoundRHSVector[1] += one_third * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);
  }

  if (rGeom[0].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE))
  {

    array_1d<double, 3> AccA(3, 0.0);
    array_1d<double, 3> AccB(3, 0.0);
    array_1d<double, 3> MeanAcc(3, 0.0);
    array_1d<double, 3> NormalVector(3, 0.0);
    const double factor = 0.5 / TimeStep;
    const double one_third = 1.0 / 3.0;

    this->GetOutwardsUnitNormalForTwoPoints(NormalVector, 0, 2, 1);

    double SpatialDefRateNormalProjection = this->CalcNormalProjectionDefRate(SpatialDefRate, NormalVector);

    noalias(AccA) = factor * (rGeom[0].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[0].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION, 1);
    noalias(AccB) = factor * (rGeom[2].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[2].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION, 1);
    noalias(MeanAcc) = 0.5 * AccA + 0.5 * AccB;

    const double accelerationsNormalProjection = MeanAcc[0] * NormalVector[0] + MeanAcc[1] * NormalVector[1];

    if (rGeom[0].IsNot(INLET)) //to change into moving wall!!!!!
      BoundRHSVector[0] += one_third * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

    if (rGeom[2].IsNot(INLET))
      BoundRHSVector[2] += one_third * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);
  }

  if (rGeom[1].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE))
  {

    array_1d<double, 3> AccA(3, 0.0);
    array_1d<double, 3> AccB(3, 0.0);
    array_1d<double, 3> MeanAcc(3, 0.0);
    array_1d<double, 3> NormalVector(3, 0.0);
    const double factor = 0.5 / TimeStep;
    const double one_third = 1.0 / 3.0;

    this->GetOutwardsUnitNormalForTwoPoints(NormalVector, 1, 2, 0);

    double SpatialDefRateNormalProjection = this->CalcNormalProjectionDefRate(SpatialDefRate, NormalVector);

    noalias(AccA) = factor * (rGeom[1].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[1].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION, 1);
    noalias(AccB) = factor * (rGeom[2].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[2].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION, 1);
    noalias(MeanAcc) = 0.5 * AccA + 0.5 * AccB;

    const double accelerationsNormalProjection = MeanAcc[0] * NormalVector[0] + MeanAcc[1] * NormalVector[1];

    if (rGeom[1].IsNot(INLET))
      BoundRHSVector[1] += one_third * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

    if (rGeom[2].IsNot(INLET))
      BoundRHSVector[2] += one_third * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);
  }
}

template <>
void TwoStepUpdatedLagrangianVPImplicitFluidDEMcouplingElement<3>::ComputeBoundRHSVectorComplete(VectorType &BoundRHSVector,
                                                                                                 const double TimeStep,
                                                                                                 const double BoundRHSCoeffAcc,
                                                                                                 const double BoundRHSCoeffDev,
                                                                                                 const VectorType SpatialDefRate)
{
  GeometryType &rGeom = this->GetGeometry();

  if (rGeom[0].Is(FREE_SURFACE) && rGeom[1].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE))
  {

    array_1d<double, 3> AccA(3, 0.0);
    array_1d<double, 3> AccB(3, 0.0);
    array_1d<double, 3> AccC(3, 0.0);
    array_1d<double, 3> MeanAcc(3, 0.0);
    array_1d<double, 3> NormalVector(3, 0.0);
    const double factor = 0.5 / TimeStep;
    const double one_third = 1.0 / 3.0;

    this->GetOutwardsUnitNormalForThreePoints(NormalVector, 0, 1, 2, 3);

    double SpatialDefRateNormalProjection = this->CalcNormalProjectionDefRate(SpatialDefRate, NormalVector);

    noalias(AccA) = factor * (rGeom[0].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[0].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION, 1);
    noalias(AccB) = factor * (rGeom[1].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[1].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION, 1);
    noalias(AccC) = factor * (rGeom[2].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[2].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION, 1);

    noalias(MeanAcc) = one_third * AccA + one_third * AccB + one_third * AccC;

    const double accelerationsNormalProjection = MeanAcc[0] * NormalVector[0] + MeanAcc[1] * NormalVector[1] + MeanAcc[2] * NormalVector[2];

    if (rGeom[0].IsNot(INLET))
      BoundRHSVector[0] += 0.25 * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

    if (rGeom[1].IsNot(INLET))
      BoundRHSVector[1] += 0.25 * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

    if (rGeom[2].IsNot(INLET))
      BoundRHSVector[2] += 0.25 * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);
  }

  if (rGeom[0].Is(FREE_SURFACE) && rGeom[1].Is(FREE_SURFACE) && rGeom[3].Is(FREE_SURFACE))
  {

    array_1d<double, 3> AccA(3, 0.0);
    array_1d<double, 3> AccB(3, 0.0);
    array_1d<double, 3> AccC(3, 0.0);
    array_1d<double, 3> MeanAcc(3, 0.0);
    array_1d<double, 3> NormalVector(3, 0.0);
    const double factor = 0.5 / TimeStep;
    const double one_third = 1.0 / 3.0;

    this->GetOutwardsUnitNormalForThreePoints(NormalVector, 0, 1, 3, 2);

    double SpatialDefRateNormalProjection = this->CalcNormalProjectionDefRate(SpatialDefRate, NormalVector);

    noalias(AccA) = factor * (rGeom[0].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[0].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION, 1);
    noalias(AccB) = factor * (rGeom[1].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[1].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION, 1);
    noalias(AccC) = factor * (rGeom[3].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[3].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[3].FastGetSolutionStepValue(ACCELERATION, 1);

    noalias(MeanAcc) = one_third * AccA + one_third * AccB + one_third * AccC;

    const double accelerationsNormalProjection = MeanAcc[0] * NormalVector[0] + MeanAcc[1] * NormalVector[1] + MeanAcc[2] * NormalVector[2];

    if (rGeom[0].IsNot(INLET))
      BoundRHSVector[0] += 0.25 * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

    if (rGeom[1].IsNot(INLET))
      BoundRHSVector[1] += 0.25 * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

    if (rGeom[3].IsNot(INLET))
      BoundRHSVector[3] += 0.25 * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);
  }

  if (rGeom[0].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE) && rGeom[3].Is(FREE_SURFACE))
  {

    array_1d<double, 3> AccA(3, 0.0);
    array_1d<double, 3> AccB(3, 0.0);
    array_1d<double, 3> AccC(3, 0.0);
    array_1d<double, 3> MeanAcc(3, 0.0);
    array_1d<double, 3> NormalVector(3, 0.0);
    const double factor = 0.5 / TimeStep;
    const double one_third = 1.0 / 3.0;

    this->GetOutwardsUnitNormalForThreePoints(NormalVector, 0, 2, 3, 1);

    double SpatialDefRateNormalProjection = this->CalcNormalProjectionDefRate(SpatialDefRate, NormalVector);

    noalias(AccA) = factor * (rGeom[0].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[0].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION, 1);
    noalias(AccB) = factor * (rGeom[2].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[2].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION, 1);
    noalias(AccC) = factor * (rGeom[3].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[3].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[3].FastGetSolutionStepValue(ACCELERATION, 1);

    noalias(MeanAcc) = one_third * AccA + one_third * AccB + one_third * AccC;

    const double accelerationsNormalProjection = MeanAcc[0] * NormalVector[0] + MeanAcc[1] * NormalVector[1] + MeanAcc[2] * NormalVector[2];

    if (rGeom[0].IsNot(INLET))
      BoundRHSVector[0] += 0.25 * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

    if (rGeom[2].IsNot(INLET))
      BoundRHSVector[2] += 0.25 * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

    if (rGeom[3].IsNot(INLET))
      BoundRHSVector[3] += 0.25 * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);
  }

  if (rGeom[1].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE) && rGeom[3].Is(FREE_SURFACE))
  {

    array_1d<double, 3> AccA(3, 0.0);
    array_1d<double, 3> AccB(3, 0.0);
    array_1d<double, 3> AccC(3, 0.0);
    array_1d<double, 3> MeanAcc(3, 0.0);
    array_1d<double, 3> NormalVector(3, 0.0);
    const double factor = 0.5 / TimeStep;
    const double one_third = 1.0 / 3.0;

    this->GetOutwardsUnitNormalForThreePoints(NormalVector, 1, 2, 3, 0);

    double SpatialDefRateNormalProjection = this->CalcNormalProjectionDefRate(SpatialDefRate, NormalVector);

    noalias(AccA) = factor * (rGeom[1].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[1].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION, 1);
    noalias(AccB) = factor * (rGeom[2].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[2].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION, 1);
    noalias(AccC) = factor * (rGeom[3].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[3].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[3].FastGetSolutionStepValue(ACCELERATION, 1);

    noalias(MeanAcc) = one_third * AccA + one_third * AccB + one_third * AccC;

    const double accelerationsNormalProjection = MeanAcc[0] * NormalVector[0] + MeanAcc[1] * NormalVector[1] + MeanAcc[2] * NormalVector[2];

    if (rGeom[1].IsNot(INLET))
      BoundRHSVector[1] += 0.25 * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

    if (rGeom[2].IsNot(INLET))
      BoundRHSVector[2] += 0.25 * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

    if (rGeom[3].IsNot(INLET))
      BoundRHSVector[3] += 0.25 * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);
  }
}

template <>
void TwoStepUpdatedLagrangianVPImplicitFluidDEMcouplingElement<2>::ComputeBoundRHSVector(VectorType &BoundRHSVector,
                                                                                         const ShapeFunctionsType &rN,
                                                                                         const double TimeStep,
                                                                                         const double BoundRHSCoeffAcc,
                                                                                         const double BoundRHSCoeffDev)
{
  GeometryType &rGeom = this->GetGeometry();
  //const SizeType NumNodes = rGeom.PointsNumber();
  array_1d<double, 3> AccA(3, 0.0);
  array_1d<double, 3> AccB(3, 0.0);

  const double factor = 0.5 / TimeStep;
  const double one_third = 1.0 / 3.0;

  if (rGeom[0].Is(FREE_SURFACE) && rGeom[1].Is(FREE_SURFACE))
  {
    noalias(AccA) = factor * (rGeom[0].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[0].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION, 1);
    noalias(AccB) = factor * (rGeom[1].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[1].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION, 1);
    // noalias(AccA)=rGeom[0].FastGetSolutionStepValue(ACCELERATION,0);
    // noalias(AccB)=rGeom[1].FastGetSolutionStepValue(ACCELERATION,0);
    const array_1d<double, 3> &NormalA = rGeom[0].FastGetSolutionStepValue(NORMAL);
    const array_1d<double, 3> &NormalB = rGeom[1].FastGetSolutionStepValue(NORMAL);
    if (rGeom[0].IsNot(INLET)) //to change into moving wall!!!!!
      BoundRHSVector[0] += one_third * (BoundRHSCoeffAcc * (AccA[0] * NormalA[0] + AccA[1] * NormalA[1]) + BoundRHSCoeffDev);
    if (rGeom[1].IsNot(INLET))
      BoundRHSVector[1] += one_third * (BoundRHSCoeffAcc * (AccB[0] * NormalB[0] + AccB[1] * NormalB[1]) + BoundRHSCoeffDev);
  }
  if (rGeom[0].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE))
  {
    noalias(AccA) = factor * (rGeom[0].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[0].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION, 1);
    noalias(AccB) = factor * (rGeom[2].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[2].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION, 1);
    const array_1d<double, 3> &NormalA = rGeom[0].FastGetSolutionStepValue(NORMAL);
    const array_1d<double, 3> &NormalB = rGeom[2].FastGetSolutionStepValue(NORMAL);
    if (rGeom[0].IsNot(INLET))
      BoundRHSVector[0] += one_third * (BoundRHSCoeffAcc * (AccA[0] * NormalA[0] + AccA[1] * NormalA[1]) + BoundRHSCoeffDev);
    if (rGeom[2].IsNot(INLET))
      BoundRHSVector[2] += one_third * (BoundRHSCoeffAcc * (AccB[0] * NormalB[0] + AccB[1] * NormalB[1]) + BoundRHSCoeffDev);
  }
  if (rGeom[1].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE))
  {
    noalias(AccA) = factor * (rGeom[1].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[1].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION, 1);
    noalias(AccB) = factor * (rGeom[2].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[2].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION, 1);
    const array_1d<double, 3> &NormalA = rGeom[1].FastGetSolutionStepValue(NORMAL);
    const array_1d<double, 3> &NormalB = rGeom[2].FastGetSolutionStepValue(NORMAL);
    if (rGeom[1].IsNot(INLET))
      BoundRHSVector[1] += one_third * (BoundRHSCoeffAcc * (AccA[0] * NormalA[0] + AccA[1] * NormalA[1]) + BoundRHSCoeffDev);
    if (rGeom[2].IsNot(INLET))
      BoundRHSVector[2] += one_third * (BoundRHSCoeffAcc * (AccB[0] * NormalB[0] + AccB[1] * NormalB[1]) + BoundRHSCoeffDev);
  }
}

template <>
void TwoStepUpdatedLagrangianVPImplicitFluidDEMcouplingElement<3>::ComputeBoundRHSVector(VectorType &BoundRHSVector,
                                                                                         const ShapeFunctionsType &rN,
                                                                                         const double TimeStep,
                                                                                         const double BoundRHSCoeffAcc,
                                                                                         const double BoundRHSCoeffDev)
{
  GeometryType &rGeom = this->GetGeometry();
  //const SizeType NumNodes = rGeom.PointsNumber();
  array_1d<double, 3> AccA(3, 0.0);
  array_1d<double, 3> AccB(3, 0.0);
  array_1d<double, 3> AccC(3, 0.0);

  const double factor = 0.5 / TimeStep;

  if (rGeom[0].Is(FREE_SURFACE) && rGeom[1].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE))
  {
    noalias(AccA) = factor * (rGeom[0].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[0].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION, 1);
    noalias(AccB) = factor * (rGeom[1].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[1].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION, 1);
    noalias(AccC) = factor * (rGeom[2].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[2].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION, 1);
    const array_1d<double, 3> &NormalA = rGeom[0].FastGetSolutionStepValue(NORMAL);
    const array_1d<double, 3> &NormalB = rGeom[1].FastGetSolutionStepValue(NORMAL);
    const array_1d<double, 3> &NormalC = rGeom[2].FastGetSolutionStepValue(NORMAL);
    if (rGeom[0].IsNot(INLET))
      BoundRHSVector[0] += 0.25 * (BoundRHSCoeffAcc * (AccA[0] * NormalA[0] + AccA[1] * NormalA[1] + AccA[2] * NormalA[2]) +
                                   BoundRHSCoeffDev);
    if (rGeom[1].IsNot(INLET))
      BoundRHSVector[1] += 0.25 * (BoundRHSCoeffAcc * (AccB[0] * NormalB[0] + AccB[1] * NormalB[1] + AccB[2] * NormalB[2]) +
                                   BoundRHSCoeffDev);
    if (rGeom[2].IsNot(INLET))
      BoundRHSVector[2] += 0.25 * (BoundRHSCoeffAcc * (AccC[0] * NormalC[0] + AccC[1] * NormalC[1] + AccC[2] * NormalC[2]) +
                                   BoundRHSCoeffDev);
  }
  if (rGeom[0].Is(FREE_SURFACE) && rGeom[1].Is(FREE_SURFACE) && rGeom[3].Is(FREE_SURFACE))
  {
    noalias(AccA) = factor * (rGeom[0].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[0].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION, 1);
    noalias(AccB) = factor * (rGeom[1].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[1].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION, 1);
    noalias(AccC) = factor * (rGeom[3].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[3].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[3].FastGetSolutionStepValue(ACCELERATION, 1);
    const array_1d<double, 3> &NormalA = rGeom[0].FastGetSolutionStepValue(NORMAL);
    const array_1d<double, 3> &NormalB = rGeom[1].FastGetSolutionStepValue(NORMAL);
    const array_1d<double, 3> &NormalC = rGeom[3].FastGetSolutionStepValue(NORMAL);
    if (rGeom[0].IsNot(INLET))
      BoundRHSVector[0] += 0.25 * (BoundRHSCoeffAcc * (AccA[0] * NormalA[0] + AccA[1] * NormalA[1] + AccA[2] * NormalA[2]) +
                                   BoundRHSCoeffDev);
    if (rGeom[1].IsNot(INLET))
      BoundRHSVector[1] += 0.25 * (BoundRHSCoeffAcc * (AccB[0] * NormalB[0] + AccB[1] * NormalB[1] + AccB[2] * NormalB[2]) +
                                   BoundRHSCoeffDev);
    if (rGeom[3].IsNot(INLET))
      BoundRHSVector[3] += 0.25 * (BoundRHSCoeffAcc * (AccC[0] * NormalC[0] + AccC[1] * NormalC[1] + AccC[2] * NormalC[2]) +
                                   BoundRHSCoeffDev);
  }
  if (rGeom[0].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE) && rGeom[3].Is(FREE_SURFACE))
  {
    noalias(AccA) = factor * (rGeom[0].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[0].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION, 1);
    noalias(AccB) = factor * (rGeom[2].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[2].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION, 1);
    noalias(AccC) = factor * (rGeom[3].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[3].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[3].FastGetSolutionStepValue(ACCELERATION, 1);
    const array_1d<double, 3> &NormalA = rGeom[0].FastGetSolutionStepValue(NORMAL);
    const array_1d<double, 3> &NormalB = rGeom[2].FastGetSolutionStepValue(NORMAL);
    const array_1d<double, 3> &NormalC = rGeom[3].FastGetSolutionStepValue(NORMAL);
    if (rGeom[0].IsNot(INLET))
      BoundRHSVector[0] += 0.25 * (BoundRHSCoeffAcc * (AccA[0] * NormalA[0] + AccA[1] * NormalA[1] + AccA[2] * NormalA[2]) +
                                   BoundRHSCoeffDev);
    if (rGeom[2].IsNot(INLET))
      BoundRHSVector[2] += 0.25 * (BoundRHSCoeffAcc * (AccB[0] * NormalB[0] + AccB[1] * NormalB[1] + AccB[2] * NormalB[2]) +
                                   BoundRHSCoeffDev);
    if (rGeom[3].IsNot(INLET))
      BoundRHSVector[3] += 0.25 * (BoundRHSCoeffAcc * (AccC[0] * NormalC[0] + AccC[1] * NormalC[1] + AccC[2] * NormalC[2]) +
                                   BoundRHSCoeffDev);
  }
  if (rGeom[1].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE) && rGeom[3].Is(FREE_SURFACE))
  {
    noalias(AccA) = factor * (rGeom[1].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[1].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION, 1);
    noalias(AccB) = factor * (rGeom[2].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[2].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION, 1);
    noalias(AccC) = factor * (rGeom[3].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[3].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[3].FastGetSolutionStepValue(ACCELERATION, 1);
    const array_1d<double, 3> &NormalA = rGeom[1].FastGetSolutionStepValue(NORMAL);
    const array_1d<double, 3> &NormalB = rGeom[2].FastGetSolutionStepValue(NORMAL);
    const array_1d<double, 3> &NormalC = rGeom[3].FastGetSolutionStepValue(NORMAL);
    if (rGeom[1].IsNot(INLET))
      BoundRHSVector[1] += 0.25 * (BoundRHSCoeffAcc * (AccA[0] * NormalA[0] + AccA[1] * NormalA[1] + AccA[2] * NormalA[2]) +
                                   BoundRHSCoeffDev);
    if (rGeom[2].IsNot(INLET))
      BoundRHSVector[2] += 0.25 * (BoundRHSCoeffAcc * (AccB[0] * NormalB[0] + AccB[1] * NormalB[1] + AccB[2] * NormalB[2]) +
                                   BoundRHSCoeffDev);
    if (rGeom[3].IsNot(INLET))
      BoundRHSVector[3] += 0.25 * (BoundRHSCoeffAcc * (AccC[0] * NormalC[0] + AccC[1] * NormalC[1] + AccC[2] * NormalC[2]) +
                                   BoundRHSCoeffDev);
  }
}

template <unsigned int TDim>
void TwoStepUpdatedLagrangianVPImplicitFluidDEMcouplingElement<TDim>::CalculateTauFIC(double &Tau,
                                                                                      double ElemSize,
                                                                                      const double Density,
                                                                                      const double Viscosity,
                                                                                      const ProcessInfo &rCurrentProcessInfo)
{
  double DeltaTime = rCurrentProcessInfo.GetValue(DELTA_TIME);
  if (rCurrentProcessInfo.GetValue(DELTA_TIME) < rCurrentProcessInfo.GetValue(PREVIOUS_DELTA_TIME))
  {
    DeltaTime = 0.5 * rCurrentProcessInfo.GetValue(DELTA_TIME) + 0.5 * rCurrentProcessInfo.GetValue(PREVIOUS_DELTA_TIME);
  }

  double MeanVelocity = 0;
  this->CalcMeanVelocity(MeanVelocity, 0);

  // Tau = 1.0 / (2.0 * Density *(0.5 * MeanVelocity / ElemSize + 0.5/DeltaTime) +  8.0 * Viscosity / (ElemSize * ElemSize) );
  Tau = (ElemSize * ElemSize * DeltaTime) / (Density * MeanVelocity * DeltaTime * ElemSize + Density * ElemSize * ElemSize + 8.0 * Viscosity * DeltaTime);
  // if(Tau<0.0000001){
  //   Tau=0.0000001;
  // }
  // if(Tau>0.0001){
  //   Tau=0.0001;
  // }

  if (MeanVelocity == 0)
  {
    Tau = 0;
  }
}

template <unsigned int TDim>
void TwoStepUpdatedLagrangianVPImplicitFluidDEMcouplingElement<TDim>::CalculateLocalContinuityEqForPressure(MatrixType &rLeftHandSideMatrix, VectorType &rRightHandSideVector, ProcessInfo &rCurrentProcessInfo)
{

  GeometryType &rGeom = this->GetGeometry();
  const unsigned int NumNodes = rGeom.PointsNumber();

  // Check sizes and initialize
  if (rLeftHandSideMatrix.size1() != NumNodes)
    rLeftHandSideMatrix.resize(NumNodes, NumNodes, false);

  rLeftHandSideMatrix = ZeroMatrix(NumNodes, NumNodes);

  if (rRightHandSideVector.size() != NumNodes)
    rRightHandSideVector.resize(NumNodes);

  rRightHandSideVector = ZeroVector(NumNodes);

  // Shape functions and integration points
  ShapeFunctionDerivativesArrayType DN_DX;
  Matrix NContainer;
  VectorType GaussWeights;
  this->CalculateGeometryData(DN_DX, NContainer, GaussWeights);
  const unsigned int NumGauss = GaussWeights.size();

  double TimeStep = rCurrentProcessInfo[DELTA_TIME];
  double theta = this->GetThetaContinuity();
  double ElemSize = this->ElementSize();

  ElementalVariables rElementalVariables;
  this->InitializeElementalVariables(rElementalVariables);

  double maxViscousValueForStabilization = 0.1;
  double Density = this->mMaterialDensity;
  double VolumetricCoeff = this->mMaterialVolumetricCoefficient;
  double DeviatoricCoeff = this->mMaterialDeviatoricCoefficient;
  double FluidFraction = 0.0;
  double FluidFractionRate = 0.0;

  if (DeviatoricCoeff > maxViscousValueForStabilization)
  {
    DeviatoricCoeff = maxViscousValueForStabilization;
  }

  double Tau = 0;
  this->CalculateTauFIC(Tau, ElemSize, Density, DeviatoricCoeff, rCurrentProcessInfo);

  double totalVolume = 0;
  bool computeElement = false;
  // Loop on integration points
  for (unsigned int g = 0; g < NumGauss; g++)
  {
    const double GaussWeight = GaussWeights[g];
    totalVolume += GaussWeight;
    const ShapeFunctionsType &N = row(NContainer, g);
    const ShapeFunctionDerivativesType &rDN_DX = DN_DX[g];
    computeElement = this->CalcCompleteStrainRate(rElementalVariables, rCurrentProcessInfo, rDN_DX, theta);
    bool wallElement = false;
    if (computeElement == true)
    {

      this->EvaluateInPoint(FluidFraction, FLUID_FRACTION, N);
      this->EvaluateInPoint(FluidFractionRate, FLUID_FRACTION_RATE, N);

      for (SizeType i = 0; i < NumNodes; ++i)
      {
        if (rGeom[i].Is(RIGID))
        {
          wallElement = true;
          break;
        }
      }
      if (wallElement == true)
      {
        FluidFractionRate = this->GetProperties()[FLUID_FRACTION_RATE];
      }

      if (std::abs(FluidFraction) < 1.0e-12)
      {
        FluidFraction = 1.0;
        FluidFractionRate = 0.0;
      }
      FluidFractionRate = 0.0;

      // double BulkCoeff =GaussWeight/(VolumetricCoeff);
      // this->ComputeBulkMatrix(BulkVelMatrix,N,BulkCoeff);
      // double BulkStabCoeff=BulkCoeff*Tau*Density/TimeStep;
      // this->ComputeBulkMatrix(BulkAccMatrix,N,BulkStabCoeff);

      double BoundLHSCoeff = Tau * 4.0 * GaussWeight / (ElemSize * ElemSize);
      // if(TDim==3){
      //   BoundLHSCoeff=Tau*2*GaussWeight/(0.81649658*ElemSize*ElemSize);
      // }

      this->ComputeBoundLHSMatrix(rLeftHandSideMatrix, N, BoundLHSCoeff);

      double BoundRHSCoeffAcc = Tau * Density * 2 * GaussWeight / ElemSize;
      double BoundRHSCoeffDev = Tau * 8.0 * DeviatoricCoeff * GaussWeight / (ElemSize * ElemSize);
      // double NProjSpatialDefRate=this->CalcNormalProjectionDefRate(rElementalVariables.SpatialDefRate);
      // double BoundRHSCoeffDev=Tau*8.0*NProjSpatialDefRate*DeviatoricCoeff*GaussWeight/(ElemSize*ElemSize);

      // this->ComputeBoundRHSVector(rRightHandSideVector,N,TimeStep,BoundRHSCoeffAcc,BoundRHSCoeffDev);
      this->ComputeBoundRHSVectorComplete(rRightHandSideVector, TimeStep, BoundRHSCoeffAcc, BoundRHSCoeffDev, rElementalVariables.SpatialDefRate);

      double StabLaplacianWeight = Tau * GaussWeight;
      this->ComputeStabLaplacianMatrix(rLeftHandSideMatrix, rDN_DX, StabLaplacianWeight);

      for (SizeType i = 0; i < NumNodes; ++i)
      {
        // RHS contribution
        // Velocity divergence
        rRightHandSideVector[i] += GaussWeight * N[i] * rElementalVariables.VolumetricDefRate;

        rRightHandSideVector[i] += GaussWeight * N[i] * FluidFractionRate / FluidFraction;

        this->AddStabilizationNodalTermsRHS(rRightHandSideVector, Tau, Density, GaussWeight, rDN_DX, i);
      }
    }
  }

  if (computeElement == true)
  {

    VectorType PressureValues = ZeroVector(NumNodes);
    VectorType PressureValuesForRHS = ZeroVector(NumNodes);
    this->GetPressureValues(PressureValuesForRHS, 0);
    //the LHS matrix up to now just contains the laplacian term and the bound term
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, PressureValuesForRHS);

    this->GetPressureValues(PressureValues, 1);
    noalias(PressureValuesForRHS) += -PressureValues;
    MatrixType BulkMatrix = ZeroMatrix(NumNodes, NumNodes);
    MatrixType BulkMatrixConsistent = ZeroMatrix(NumNodes, NumNodes);
    double lumpedBulkCoeff = totalVolume / (VolumetricCoeff);
    double lumpedBulkStabCoeff = lumpedBulkCoeff * Tau * Density / TimeStep;

    this->ComputeBulkMatrixLump(BulkMatrix, lumpedBulkCoeff);
    // this->ComputeBulkMatrixConsistent(BulkMatrixConsistent,lumpedBulkCoeff);
    noalias(rLeftHandSideMatrix) += BulkMatrix;
    // noalias(rLeftHandSideMatrix)+=BulkMatrixConsistent;
    noalias(rRightHandSideVector) -= prod(BulkMatrix, PressureValuesForRHS);
    // noalias(rRightHandSideVector) -=prod(BulkMatrixConsistent,PressureValuesForRHS);

    this->GetPressureVelocityValues(PressureValues, 0);
    noalias(PressureValuesForRHS) += -PressureValues * TimeStep;
    noalias(BulkMatrix) = ZeroMatrix(NumNodes, NumNodes);
    this->ComputeBulkMatrixLump(BulkMatrix, lumpedBulkStabCoeff);
    // this->ComputeBulkMatrixConsistent(BulkMatrixConsistent,lumpedBulkStabCoeff);
    noalias(rLeftHandSideMatrix) += BulkMatrix;
    // noalias(rLeftHandSideMatrix)+=BulkMatrixConsistent;
    noalias(rRightHandSideVector) -= prod(BulkMatrix, PressureValuesForRHS);
    // noalias(rRightHandSideVector) -=prod(BulkMatrixConsistent,PressureValuesForRHS);
  }
  else
  {
    double lumpedBulkCoeff = totalVolume * Tau * Density / (TimeStep * VolumetricCoeff);
    MatrixType BulkVelMatrixLump = ZeroMatrix(NumNodes, NumNodes);
    this->ComputeBulkMatrixLump(BulkVelMatrixLump, lumpedBulkCoeff);
    noalias(rLeftHandSideMatrix) += BulkVelMatrixLump;
    VectorType PressureValues = ZeroVector(NumNodes);
    VectorType PressureValuesForRHS = ZeroVector(NumNodes);
    this->GetPressureValues(PressureValuesForRHS, 0);
    this->GetPressureValues(PressureValues, 1);
    noalias(PressureValuesForRHS) += -PressureValues;
    noalias(rRightHandSideVector) -= prod(BulkVelMatrixLump, PressureValuesForRHS);
  }
}

template class TwoStepUpdatedLagrangianVPImplicitFluidDEMcouplingElement<2>;
template class TwoStepUpdatedLagrangianVPImplicitFluidDEMcouplingElement<3>;

} // namespace Kratos
