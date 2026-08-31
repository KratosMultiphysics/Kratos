//
//   Project Name:        KratosFluidDynamicsApplication $
//   Last modified by:    $Author:               AFranci $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_elements/two_step_updated_lagrangian_V_P_implicit_fluid_FIC_element.h"
#include "includes/cfd_variables.h"
#include <cmath>

namespace Kratos
{

  /*
   * public TwoStepUpdatedLagrangianVPImplicitFluidFicElement<TDim> functions
   */

  template <unsigned int TDim>
  Element::Pointer TwoStepUpdatedLagrangianVPImplicitFluidFicElement<TDim>::Clone(IndexType NewId, NodesArrayType const &rThisNodes) const
  {

    TwoStepUpdatedLagrangianVPImplicitFluidFicElement NewElement(NewId, this->GetGeometry().Create(rThisNodes), this->pGetProperties());

    NewElement.SetData(this->GetData());
    NewElement.SetFlags(this->GetFlags());

    return Element::Pointer(new TwoStepUpdatedLagrangianVPImplicitFluidFicElement(NewElement));
  }

  template <unsigned int TDim>
  int TwoStepUpdatedLagrangianVPImplicitFluidFicElement<TDim>::Check(const ProcessInfo &rCurrentProcessInfo) const
  {
    KRATOS_TRY;

    // Base class checks for positive Jacobian and Id > 0
    int ierr = Element::Check(rCurrentProcessInfo);
    if (ierr != 0)
      return ierr;

    // Check that all required variables have been registered
    if (VELOCITY.Key() == 0)
      KRATOS_THROW_ERROR(std::invalid_argument,
                         "VELOCITY Key is 0. Check that the application was correctly registered.", "");
    if (ACCELERATION.Key() == 0)
      KRATOS_THROW_ERROR(std::invalid_argument,
                         "ACCELERATION Key is 0. Check that the application was correctly registered.", "");
    if (PRESSURE.Key() == 0)
      KRATOS_THROW_ERROR(std::invalid_argument,
                         "PRESSURE Key is 0. Check that the application was correctly registered.", "");
    if (BODY_FORCE.Key() == 0)
      KRATOS_THROW_ERROR(std::invalid_argument,
                         "BODY_FORCE Key is 0. Check that the application was correctly registered.", "");
    if (DENSITY.Key() == 0)
      KRATOS_THROW_ERROR(std::invalid_argument,
                         "DENSITY Key is 0. Check that the application was correctly registered.", "");
    if (DYNAMIC_VISCOSITY.Key() == 0)
      KRATOS_THROW_ERROR(std::invalid_argument,
                         "DYNAMIC_VISCOSITY Key is 0. Check that the application was correctly registered.", "");
    if (DELTA_TIME.Key() == 0)
      KRATOS_THROW_ERROR(std::invalid_argument,
                         "DELTA_TIME Key is 0. Check that the application was correctly registered.", "");

    const GeometryType &r_geom = this->GetGeometry();
    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for (unsigned int i = 0; i < r_geom.size(); ++i)
    {
      if (!r_geom[i].SolutionStepsDataHas(VELOCITY))
        KRATOS_THROW_ERROR(std::invalid_argument, "missing VELOCITY variable on solution step data for node ", r_geom[i].Id());
      if (!r_geom[i].SolutionStepsDataHas(PRESSURE))
        KRATOS_THROW_ERROR(std::invalid_argument, "missing PRESSURE variable on solution step data for node ", r_geom[i].Id());
      if (!r_geom[i].SolutionStepsDataHas(BODY_FORCE))
        KRATOS_THROW_ERROR(std::invalid_argument, "missing BODY_FORCE variable on solution step data for node ", r_geom[i].Id());
      if (!r_geom[i].SolutionStepsDataHas(DENSITY))
        KRATOS_THROW_ERROR(std::invalid_argument, "missing DENSITY variable on solution step data for node ", r_geom[i].Id());
      if (!r_geom[i].SolutionStepsDataHas(DYNAMIC_VISCOSITY))
        KRATOS_THROW_ERROR(std::invalid_argument, "missing DYNAMIC_VISCOSITY variable on solution step data for node ", r_geom[i].Id());
      if (!r_geom[i].HasDofFor(VELOCITY_X) ||
          !r_geom[i].HasDofFor(VELOCITY_Y) ||
          !r_geom[i].HasDofFor(VELOCITY_Z))
        KRATOS_THROW_ERROR(std::invalid_argument, "missing VELOCITY component degree of freedom on node ", r_geom[i].Id());
      if (!r_geom[i].HasDofFor(PRESSURE))
        KRATOS_THROW_ERROR(std::invalid_argument, "missing PRESSURE component degree of freedom on node ", r_geom[i].Id());
    }

    // If this is a 2D problem, check that nodes are in XY plane
    if (r_geom.WorkingSpaceDimension() == 2)
    {
      for (unsigned int i = 0; i < this->GetGeometry().size(); ++i)
      {
        if (r_geom[i].Z() != 0.0)
          KRATOS_THROW_ERROR(std::invalid_argument, "Node with non-zero Z coordinate found. Id: ", r_geom[i].Id());
      }
    }

    // Consitutive law checks
    const auto &r_properties = this->GetProperties();
    const auto &r_geometry = this->GetGeometry();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();

    // WARNING THIS MUST BE REMOVED ASAP
    const_cast<TwoStepUpdatedLagrangianVPImplicitFluidFicElement<TDim> *>(this)->mpConstitutiveLaw = const_cast<TwoStepUpdatedLagrangianVPImplicitFluidFicElement<TDim> *>(this)->GetProperties().GetValue(CONSTITUTIVE_LAW);
    // mpConstitutiveLaw = this->GetProperties().GetValue(CONSTITUTIVE_LAW);

    // Verify that the constitutive law exists
    KRATOS_ERROR_IF_NOT(r_properties.Has(CONSTITUTIVE_LAW))
        << "Constitutive law not provided for property " << r_properties.Id() << std::endl;

    // Verify that the constitutive law has the correct dimension
    const SizeType strain_size = r_properties.GetValue(CONSTITUTIVE_LAW)->GetStrainSize();

    if (dimension == 2)
    {
      KRATOS_ERROR_IF(strain_size < 3 || strain_size > 4)
          << "Wrong constitutive law used. This is a 2D element! expected strain size is 3 or 4 "
             "(el id = ) "
          << this->Id() << std::endl;
    }
    else
    {
      KRATOS_ERROR_IF_NOT(strain_size == 6) << "Wrong constitutive law used. This is a 3D element! "
                                               "expected strain size is 6 (el id = ) "
                                            << this->Id() << std::endl;
    }

    // Check constitutive law
    return this->mpConstitutiveLaw->Check(r_properties, r_geometry, rCurrentProcessInfo);

    return ierr;

    KRATOS_CATCH("");
  }

  template <>
  void TwoStepUpdatedLagrangianVPImplicitFluidFicElement<2>::ComputeBoundLHSMatrix(Matrix &BoundLHSMatrix,
                                                                                   const ShapeFunctionsType &rN,
                                                                                   const double Weight)
  {
    const GeometryType &rGeom = this->GetGeometry();
    constexpr double coeff = 1.0 / 3.0;

    if (rGeom[0].Is(FREE_SURFACE) && rGeom[1].Is(FREE_SURFACE))
    {
      if (rGeom[0].IsNot(INLET))
        BoundLHSMatrix(0, 0) += Weight * coeff;
      if (rGeom[1].IsNot(INLET))
        BoundLHSMatrix(1, 1) += Weight * coeff;
    }
    if (rGeom[0].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE))
    {
      if (rGeom[0].IsNot(INLET))
        BoundLHSMatrix(0, 0) += Weight * coeff;
      if (rGeom[2].IsNot(INLET))
        BoundLHSMatrix(2, 2) += Weight * coeff;
    }
    if (rGeom[1].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE))
    {
      if (rGeom[1].IsNot(INLET))
        BoundLHSMatrix(1, 1) += Weight * coeff;
      if (rGeom[2].IsNot(INLET))
        BoundLHSMatrix(2, 2) += Weight * coeff;
    }
    // }
  }

  template <>
  void TwoStepUpdatedLagrangianVPImplicitFluidFicElement<3>::ComputeBoundLHSMatrix(Matrix &BoundLHSMatrix,
                                                                                   const ShapeFunctionsType &rN,
                                                                                   const double Weight)
  {
    const GeometryType &rGeom = this->GetGeometry();
    constexpr double coeff = 0.25;

    if (rGeom[0].Is(FREE_SURFACE) && rGeom[1].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE))
    {
      if (rGeom[0].IsNot(INLET))
        BoundLHSMatrix(0, 0) += Weight * coeff;
      if (rGeom[1].IsNot(INLET))
        BoundLHSMatrix(1, 1) += Weight * coeff;
      if (rGeom[2].IsNot(INLET))
        BoundLHSMatrix(2, 2) += Weight * coeff;
    }
    if (rGeom[0].Is(FREE_SURFACE) && rGeom[1].Is(FREE_SURFACE) && rGeom[3].Is(FREE_SURFACE))
    {
      if (rGeom[0].IsNot(INLET))
        BoundLHSMatrix(0, 0) += Weight * coeff;
      if (rGeom[1].IsNot(INLET))
        BoundLHSMatrix(1, 1) += Weight * coeff;
      if (rGeom[3].IsNot(INLET))
        BoundLHSMatrix(3, 3) += Weight * coeff;
    }
    if (rGeom[0].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE) && rGeom[3].Is(FREE_SURFACE))
    {
      if (rGeom[0].IsNot(INLET))
        BoundLHSMatrix(0, 0) += Weight * coeff;
      if (rGeom[2].IsNot(INLET))
        BoundLHSMatrix(2, 2) += Weight * coeff;
      if (rGeom[3].IsNot(INLET))
        BoundLHSMatrix(3, 3) += Weight * coeff;
    }
    if (rGeom[1].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE) && rGeom[3].Is(FREE_SURFACE))
    {
      if (rGeom[1].IsNot(INLET))
        BoundLHSMatrix(1, 1) += Weight * coeff;
      if (rGeom[2].IsNot(INLET))
        BoundLHSMatrix(2, 2) += Weight * coeff;
      if (rGeom[3].IsNot(INLET))
        BoundLHSMatrix(3, 3) += Weight * coeff;
    }
  }

  template <>
  void TwoStepUpdatedLagrangianVPImplicitFluidFicElement<2>::ComputeBoundRHSVectorComplete(VectorType &BoundRHSVector,
                                                                                           const double TimeStep,
                                                                                           const double BoundRHSCoeffAcc,
                                                                                           const double BoundRHSCoeffDev,
                                                                                           const VectorType SpatialDefRate)
  {
    const GeometryType &rGeom = this->GetGeometry();
    constexpr double coeff = 1.0 / 3.0;
    const double timeFactor = 0.5 / TimeStep;

    if (rGeom[0].Is(FREE_SURFACE) && rGeom[1].Is(FREE_SURFACE))
    {
      array_1d<double, 3> AccA(3, 0.0);
      array_1d<double, 3> AccB(3, 0.0);
      array_1d<double, 3> MeanAcc(3, 0.0);
      array_1d<double, 3> NormalVector(3, 0.0);

      this->GetOutwardsUnitNormalForTwoPoints(NormalVector, 0, 1, 2);

      double SpatialDefRateNormalProjection = this->CalcNormalProjectionDefRate(SpatialDefRate, NormalVector);

      noalias(AccA) = timeFactor * (rGeom[0].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[0].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION, 1);
      noalias(AccB) = timeFactor * (rGeom[1].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[1].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION, 1);
      noalias(MeanAcc) = 0.5 * (AccA + AccB);

      const double accelerationsNormalProjection = MeanAcc[0] * NormalVector[0] + MeanAcc[1] * NormalVector[1];

      if (rGeom[0].IsNot(INLET)) // to change into moving wall!!!!!
        BoundRHSVector[0] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

      if (rGeom[1].IsNot(INLET))
        BoundRHSVector[1] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);
    }

    if (rGeom[0].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE))
    {

      array_1d<double, 3> AccA(3, 0.0);
      array_1d<double, 3> AccB(3, 0.0);
      array_1d<double, 3> MeanAcc(3, 0.0);
      array_1d<double, 3> NormalVector(3, 0.0);

      this->GetOutwardsUnitNormalForTwoPoints(NormalVector, 0, 2, 1);

      double SpatialDefRateNormalProjection = this->CalcNormalProjectionDefRate(SpatialDefRate, NormalVector);

      noalias(AccA) = timeFactor * (rGeom[0].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[0].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION, 1);
      noalias(AccB) = timeFactor * (rGeom[2].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[2].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION, 1);
      noalias(MeanAcc) = 0.5 * (AccA + AccB);

      const double accelerationsNormalProjection = MeanAcc[0] * NormalVector[0] + MeanAcc[1] * NormalVector[1];

      if (rGeom[0].IsNot(INLET)) // to change into moving wall!!!!!
        BoundRHSVector[0] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

      if (rGeom[2].IsNot(INLET))
        BoundRHSVector[2] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);
    }

    if (rGeom[1].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE))
    {

      array_1d<double, 3> AccA(3, 0.0);
      array_1d<double, 3> AccB(3, 0.0);
      array_1d<double, 3> MeanAcc(3, 0.0);
      array_1d<double, 3> NormalVector(3, 0.0);

      this->GetOutwardsUnitNormalForTwoPoints(NormalVector, 1, 2, 0);

      double SpatialDefRateNormalProjection = this->CalcNormalProjectionDefRate(SpatialDefRate, NormalVector);

      noalias(AccA) = timeFactor * (rGeom[1].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[1].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION, 1);
      noalias(AccB) = timeFactor * (rGeom[2].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[2].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION, 1);
      noalias(MeanAcc) = 0.5 * (AccA + AccB);

      const double accelerationsNormalProjection = MeanAcc[0] * NormalVector[0] + MeanAcc[1] * NormalVector[1];

      if (rGeom[1].IsNot(INLET))
        BoundRHSVector[1] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

      if (rGeom[2].IsNot(INLET))
        BoundRHSVector[2] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);
    }
  }

  template <>
  void TwoStepUpdatedLagrangianVPImplicitFluidFicElement<3>::ComputeBoundRHSVectorComplete(VectorType &BoundRHSVector,
                                                                                           const double TimeStep,
                                                                                           const double BoundRHSCoeffAcc,
                                                                                           const double BoundRHSCoeffDev,
                                                                                           const VectorType SpatialDefRate)
  {
    const GeometryType &rGeom = this->GetGeometry();
    constexpr double coeff = 0.25;
    const double timeFactor = 0.5 / TimeStep;
    constexpr double one_third = 1.0 / 3.0;

    if (rGeom[0].Is(FREE_SURFACE) && rGeom[1].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE))
    {

      array_1d<double, 3> AccA(3, 0.0);
      array_1d<double, 3> AccB(3, 0.0);
      array_1d<double, 3> AccC(3, 0.0);
      array_1d<double, 3> MeanAcc(3, 0.0);
      array_1d<double, 3> NormalVector(3, 0.0);

      this->GetOutwardsUnitNormalForThreePoints(NormalVector, 0, 1, 2, 3);

      double SpatialDefRateNormalProjection = this->CalcNormalProjectionDefRate(SpatialDefRate, NormalVector);

      noalias(AccA) = timeFactor * (rGeom[0].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[0].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION, 1);
      noalias(AccB) = timeFactor * (rGeom[1].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[1].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION, 1);
      noalias(AccC) = timeFactor * (rGeom[2].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[2].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION, 1);

      noalias(MeanAcc) = (AccA + AccB + AccC) * one_third;

      const double accelerationsNormalProjection = MeanAcc[0] * NormalVector[0] + MeanAcc[1] * NormalVector[1] + MeanAcc[2] * NormalVector[2];

      if (rGeom[0].IsNot(INLET))
        BoundRHSVector[0] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

      if (rGeom[1].IsNot(INLET))
        BoundRHSVector[1] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

      if (rGeom[2].IsNot(INLET))
        BoundRHSVector[2] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);
    }

    if (rGeom[0].Is(FREE_SURFACE) && rGeom[1].Is(FREE_SURFACE) && rGeom[3].Is(FREE_SURFACE))
    {

      array_1d<double, 3> AccA(3, 0.0);
      array_1d<double, 3> AccB(3, 0.0);
      array_1d<double, 3> AccC(3, 0.0);
      array_1d<double, 3> MeanAcc(3, 0.0);
      array_1d<double, 3> NormalVector(3, 0.0);
      this->GetOutwardsUnitNormalForThreePoints(NormalVector, 0, 1, 3, 2);

      double SpatialDefRateNormalProjection = this->CalcNormalProjectionDefRate(SpatialDefRate, NormalVector);

      noalias(AccA) = timeFactor * (rGeom[0].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[0].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION, 1);
      noalias(AccB) = timeFactor * (rGeom[1].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[1].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION, 1);
      noalias(AccC) = timeFactor * (rGeom[3].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[3].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[3].FastGetSolutionStepValue(ACCELERATION, 1);

      noalias(MeanAcc) = (AccA + AccB + AccC) * one_third;

      const double accelerationsNormalProjection = MeanAcc[0] * NormalVector[0] + MeanAcc[1] * NormalVector[1] + MeanAcc[2] * NormalVector[2];

      if (rGeom[0].IsNot(INLET))
        BoundRHSVector[0] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

      if (rGeom[1].IsNot(INLET))
        BoundRHSVector[1] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

      if (rGeom[3].IsNot(INLET))
        BoundRHSVector[3] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);
    }

    if (rGeom[0].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE) && rGeom[3].Is(FREE_SURFACE))
    {

      array_1d<double, 3> AccA(3, 0.0);
      array_1d<double, 3> AccB(3, 0.0);
      array_1d<double, 3> AccC(3, 0.0);
      array_1d<double, 3> MeanAcc(3, 0.0);
      array_1d<double, 3> NormalVector(3, 0.0);

      this->GetOutwardsUnitNormalForThreePoints(NormalVector, 0, 2, 3, 1);

      double SpatialDefRateNormalProjection = this->CalcNormalProjectionDefRate(SpatialDefRate, NormalVector);

      noalias(AccA) = timeFactor * (rGeom[0].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[0].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION, 1);
      noalias(AccB) = timeFactor * (rGeom[2].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[2].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION, 1);
      noalias(AccC) = timeFactor * (rGeom[3].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[3].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[3].FastGetSolutionStepValue(ACCELERATION, 1);

      noalias(MeanAcc) = (AccA + AccB + AccC) * one_third;

      const double accelerationsNormalProjection = MeanAcc[0] * NormalVector[0] + MeanAcc[1] * NormalVector[1] + MeanAcc[2] * NormalVector[2];

      if (rGeom[0].IsNot(INLET))
        BoundRHSVector[0] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

      if (rGeom[2].IsNot(INLET))
        BoundRHSVector[2] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

      if (rGeom[3].IsNot(INLET))
        BoundRHSVector[3] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);
    }

    if (rGeom[1].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE) && rGeom[3].Is(FREE_SURFACE))
    {

      array_1d<double, 3> AccA(3, 0.0);
      array_1d<double, 3> AccB(3, 0.0);
      array_1d<double, 3> AccC(3, 0.0);
      array_1d<double, 3> MeanAcc(3, 0.0);
      array_1d<double, 3> NormalVector(3, 0.0);

      this->GetOutwardsUnitNormalForThreePoints(NormalVector, 1, 2, 3, 0);

      double SpatialDefRateNormalProjection = this->CalcNormalProjectionDefRate(SpatialDefRate, NormalVector);

      noalias(AccA) = timeFactor * (rGeom[1].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[1].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION, 1);
      noalias(AccB) = timeFactor * (rGeom[2].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[2].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION, 1);
      noalias(AccC) = timeFactor * (rGeom[3].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[3].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[3].FastGetSolutionStepValue(ACCELERATION, 1);

      noalias(MeanAcc) = (AccA + AccB + AccC) * one_third;

      const double accelerationsNormalProjection = MeanAcc[0] * NormalVector[0] + MeanAcc[1] * NormalVector[1] + MeanAcc[2] * NormalVector[2];

      if (rGeom[1].IsNot(INLET))
        BoundRHSVector[1] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

      if (rGeom[2].IsNot(INLET))
        BoundRHSVector[2] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

      if (rGeom[3].IsNot(INLET))
        BoundRHSVector[3] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);
    }
  }

  template <>
  void TwoStepUpdatedLagrangianVPImplicitFluidFicElement<2>::ComputeBoundRHSVector(VectorType &BoundRHSVector,
                                                                                   const ShapeFunctionsType &rN,
                                                                                   const double TimeStep,
                                                                                   const double BoundRHSCoeffAcc,
                                                                                   const double BoundRHSCoeffDev)
  {
    const GeometryType &rGeom = this->GetGeometry();
    array_1d<double, 3> AccA(3, 0.0);
    array_1d<double, 3> AccB(3, 0.0);
    const double factor = 0.5 / TimeStep;
    constexpr double one_third = 1.0 / 3.0;

    if (rGeom[0].Is(FREE_SURFACE) && rGeom[1].Is(FREE_SURFACE))
    {
      noalias(AccA) = factor * (rGeom[0].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[0].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION, 1);
      noalias(AccB) = factor * (rGeom[1].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[1].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION, 1);
      const array_1d<double, 3> &NormalA = rGeom[0].FastGetSolutionStepValue(NORMAL);
      const array_1d<double, 3> &NormalB = rGeom[1].FastGetSolutionStepValue(NORMAL);
      if (rGeom[0].IsNot(INLET)) // to change into moving wall!!!!!
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
  void TwoStepUpdatedLagrangianVPImplicitFluidFicElement<3>::ComputeBoundRHSVector(VectorType &BoundRHSVector,
                                                                                   const ShapeFunctionsType &rN,
                                                                                   const double TimeStep,
                                                                                   const double BoundRHSCoeffAcc,
                                                                                   const double BoundRHSCoeffDev)
  {
    const GeometryType &rGeom = this->GetGeometry();
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
  void TwoStepUpdatedLagrangianVPImplicitFluidFicElement<TDim>::CalculateTauFIC(double &Tau,
                                                                                const double ElemSize,
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
    this->CalcMeanVelocityNorm(MeanVelocity, 0);

    Tau = (ElemSize * ElemSize * DeltaTime) / (Density * MeanVelocity * DeltaTime * ElemSize + Density * ElemSize * ElemSize + 8.0 * Viscosity * DeltaTime);

    constexpr double tolerance = 1.0e-13;
    if (MeanVelocity < tolerance)
    {
      Tau = 0;
    }
  }

  template <unsigned int TDim>
  void TwoStepUpdatedLagrangianVPImplicitFluidFicElement<TDim>::AddStabilizationMatrixLHS(MatrixType &rLeftHandSideMatrix,
                                                                                          Matrix &BulkAccMatrix,
                                                                                          const ShapeFunctionsType &rN,
                                                                                          const double Weight)
  {
    const SizeType NumNodes = this->GetGeometry().PointsNumber();

    if (BulkAccMatrix.size1() != NumNodes)
      BulkAccMatrix.resize(NumNodes, NumNodes, false);

    noalias(BulkAccMatrix) = ZeroMatrix(NumNodes, NumNodes);
    for (SizeType i = 0; i < NumNodes; ++i)
    {
      // LHS contribution
      for (SizeType j = 0; j < NumNodes; ++j)
      {
        BulkAccMatrix(i, j) += Weight * rN[i] * rN[j];
      }
    }
    rLeftHandSideMatrix += BulkAccMatrix;
  }

  template <unsigned int TDim>
  void TwoStepUpdatedLagrangianVPImplicitFluidFicElement<TDim>::ComputeStabLaplacianMatrix(MatrixType &StabLaplacianMatrix,
                                                                                           const ShapeFunctionDerivativesType &rDN_DX,
                                                                                           const double Weight)

  {
    // LHS contribution
    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    for (SizeType i = 0; i < NumNodes; ++i)
    {
      for (SizeType j = 0; j < NumNodes; ++j)
      {
        for (SizeType d = 0; d < TDim; ++d)
        {
          StabLaplacianMatrix(i, j) += Weight * rDN_DX(i, d) * rDN_DX(j, d);
        }
      }
    }
  }

  template <unsigned int TDim>
  void TwoStepUpdatedLagrangianVPImplicitFluidFicElement<TDim>::AddStabilizationNodalTermsLHS(MatrixType &rLeftHandSideMatrix,
                                                                                              const double Tau,
                                                                                              const double Weight,
                                                                                              const ShapeFunctionDerivativesType &rDN_DX,
                                                                                              const SizeType i)
  {
    // LHS contribution
    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    for (SizeType j = 0; j < NumNodes; ++j)
    {
      double Lij = 0.0;
      for (SizeType d = 0; d < TDim; ++d)
      {
        Lij += rDN_DX(i, d) * rDN_DX(j, d);
      }
      Lij *= Tau;

      rLeftHandSideMatrix(i, j) += Weight * Lij;
    }
  }

  template <unsigned int TDim>
  void TwoStepUpdatedLagrangianVPImplicitFluidFicElement<TDim>::AddStabilizationNodalTermsRHS(VectorType &rRightHandSideVector,
                                                                                              const double Tau,
                                                                                              const double Density,
                                                                                              const double Weight,
                                                                                              const ShapeFunctionDerivativesType &rDN_DX,
                                                                                              const SizeType i)
  {

    double RHSi = 0;
    if (this->GetGeometry()[i].SolutionStepsDataHas(VOLUME_ACCELERATION))
    { // it must be checked once at the beginning only
      const array_1d<double, 3> &VolumeAcceleration = this->GetGeometry()[i].FastGetSolutionStepValue(VOLUME_ACCELERATION);

      for (SizeType d = 0; d < TDim; ++d)
      {
        RHSi += -rDN_DX(i, d) * Tau * (Density * VolumeAcceleration[d]);
      }
    }
    rRightHandSideVector[i] += Weight * RHSi;
  }

  template <unsigned int TDim>
  void TwoStepUpdatedLagrangianVPImplicitFluidFicElement<TDim>::CalculateLocalContinuityEqForPressure(MatrixType &rLeftHandSideMatrix,
                                                                                                      VectorType &rRightHandSideVector,
                                                                                                      const ProcessInfo &rCurrentProcessInfo)
  {

    const GeometryType &rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();

    // Check sizes and initialize
    if (rLeftHandSideMatrix.size1() != NumNodes)
      rLeftHandSideMatrix.resize(NumNodes, NumNodes, false);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(NumNodes, NumNodes);
    MatrixType LaplacianMatrix = ZeroMatrix(NumNodes, NumNodes);

    if (rRightHandSideVector.size() != NumNodes)
      rRightHandSideVector.resize(NumNodes, false);

    noalias(rRightHandSideVector) = ZeroVector(NumNodes);

    // Shape functions and integration points
    ShapeFunctionDerivativesArrayType DN_DX;
    Matrix NContainer;
    VectorType GaussWeights;
    this->CalculateGeometryData(DN_DX, NContainer, GaussWeights);
    const unsigned int NumGauss = GaussWeights.size();

    const double TimeStep = rCurrentProcessInfo[DELTA_TIME];
    const double theta = this->GetThetaContinuity();
    const double ElemSize = this->ElementSize();

    ElementalVariables rElementalVariables;
    this->InitializeElementalVariables(rElementalVariables);

    double maxViscousValueForStabilization = 0.1;
    const double Density = this->mMaterialDensity;
    double VolumetricCoeff = this->mMaterialVolumetricCoefficient;
    double DeviatoricCoeff = this->mMaterialDeviatoricCoefficient;

    if (DeviatoricCoeff > maxViscousValueForStabilization)
    {
      DeviatoricCoeff = maxViscousValueForStabilization;
    }

    double Tau = 0;
    this->CalculateTauFIC(Tau, ElemSize, Density, DeviatoricCoeff, rCurrentProcessInfo);

    double totalVolume = 0;
    bool computeElement = false;
    // Loop on integration points
    for (unsigned int g = 0; g < NumGauss; ++g)
    {
      const double GaussWeight = GaussWeights[g];
      totalVolume += GaussWeight;
      const ShapeFunctionsType &N = row(NContainer, g);
      const ShapeFunctionDerivativesType &rDN_DX = DN_DX[g];
      computeElement = this->CalcCompleteStrainRate(rElementalVariables, rCurrentProcessInfo, rDN_DX, theta);

      if (computeElement && this->IsNot(BLOCKED) && this->IsNot(ISOLATED))
      {

        double BoundLHSCoeff = Tau * 4.0 * GaussWeight / (ElemSize * ElemSize);

        this->ComputeBoundLHSMatrix(rLeftHandSideMatrix, N, BoundLHSCoeff);

        double BoundRHSCoeffAcc = Tau * Density * 2 * GaussWeight / ElemSize;
        double BoundRHSCoeffDev = Tau * 8.0 * DeviatoricCoeff * GaussWeight / (ElemSize * ElemSize);
        this->ComputeBoundRHSVectorComplete(rRightHandSideVector, TimeStep, BoundRHSCoeffAcc, BoundRHSCoeffDev, rElementalVariables.SpatialDefRate);

        double StabLaplacianWeight = Tau * GaussWeight;
        this->ComputeStabLaplacianMatrix(LaplacianMatrix, rDN_DX, StabLaplacianWeight);

        array_1d<double, TDim> OldPressureGradient = ZeroVector(TDim);
        this->EvaluateGradientInPoint(OldPressureGradient, PRESSURE, rDN_DX);

        for (SizeType i = 0; i < NumNodes; ++i)
        {
          // RHS contribution
          // Velocity divergence
          rRightHandSideVector[i] += GaussWeight * N[i] * rElementalVariables.VolumetricDefRate;
          this->AddStabilizationNodalTermsRHS(rRightHandSideVector, Tau, Density, GaussWeight, rDN_DX, i);
          double laplacianRHSi = 0;
          for (SizeType d = 0; d < TDim; ++d)
          {
            laplacianRHSi += StabLaplacianWeight * rDN_DX(i, d) * OldPressureGradient[d];
          }
          rRightHandSideVector[i] += -laplacianRHSi;
        }
      }
    }

    if (computeElement && this->IsNot(BLOCKED) && this->IsNot(ISOLATED))
    {

      VectorType PressureValues = ZeroVector(NumNodes);
      VectorType PressureValuesForRHS = ZeroVector(NumNodes);
      this->GetPressureValues(PressureValuesForRHS, 0);
      // the LHS matrix up to now just contains the laplacian term and the bound term
      noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, PressureValuesForRHS);
      rLeftHandSideMatrix += LaplacianMatrix;

      this->GetPressureValues(PressureValues, 1);
      noalias(PressureValuesForRHS) += -PressureValues;
      MatrixType BulkMatrix = ZeroMatrix(NumNodes, NumNodes);
      double lumpedBulkCoeff = totalVolume / (VolumetricCoeff);
      double lumpedBulkStabCoeff = lumpedBulkCoeff * Tau * Density / TimeStep;

      this->ComputeBulkMatrixLump(BulkMatrix, lumpedBulkCoeff);
      noalias(rLeftHandSideMatrix) += BulkMatrix;
      noalias(rRightHandSideVector) -= prod(BulkMatrix, PressureValuesForRHS);

      this->GetPressureVelocityValues(PressureValues, 0);
      noalias(PressureValuesForRHS) += -PressureValues * TimeStep;
      noalias(BulkMatrix) = ZeroMatrix(NumNodes, NumNodes);
      this->ComputeBulkMatrixLump(BulkMatrix, lumpedBulkStabCoeff);
      noalias(rLeftHandSideMatrix) += BulkMatrix;
      noalias(rRightHandSideVector) -= prod(BulkMatrix, PressureValuesForRHS);
    }
    else if (this->IsNot(BLOCKED) && this->IsNot(ISOLATED))
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
    else if (this->Is(BLOCKED) && this->IsNot(ISOLATED))
    {
      VectorType PressureValues = ZeroVector(NumNodes);
      VectorType PressureValuesForRHS = ZeroVector(NumNodes);
      this->GetPressureValues(PressureValuesForRHS, 0);
      // the LHS matrix up to now is void

      this->GetPressureValues(PressureValues, 1);
      noalias(PressureValuesForRHS) += -PressureValues;
      MatrixType BulkMatrix = ZeroMatrix(NumNodes, NumNodes);
      double lumpedBulkCoeff = totalVolume / (VolumetricCoeff);

      this->ComputeBulkMatrixLump(BulkMatrix, lumpedBulkCoeff);
      noalias(rLeftHandSideMatrix) += BulkMatrix;
      noalias(rRightHandSideVector) -= prod(BulkMatrix, PressureValuesForRHS);
    }
    else if (this->Is(ISOLATED))
    {
      MatrixType BulkMatrix = ZeroMatrix(NumNodes, NumNodes);
      double lumpedBulkCoeff = totalVolume / (VolumetricCoeff);

      this->ComputeBulkMatrixLump(BulkMatrix, lumpedBulkCoeff);
      noalias(rLeftHandSideMatrix) += BulkMatrix;
    }
  }

  template <unsigned int TDim>
  void TwoStepUpdatedLagrangianVPImplicitFluidFicElement<TDim>::GetPressureAccelerationValues(Vector &rValues,
                                                                                              const int Step)
  {
    const GeometryType &rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();

    if (rValues.size() != NumNodes)
      rValues.resize(NumNodes, false);

    for (SizeType i = 0; i < NumNodes; ++i)
    {
      rValues[i] = rGeom[i].FastGetSolutionStepValue(PRESSURE_ACCELERATION, Step);
    }
  }

  template class TwoStepUpdatedLagrangianVPImplicitFluidFicElement<2>;
  template class TwoStepUpdatedLagrangianVPImplicitFluidFicElement<3>;

} // namespace Kratos
