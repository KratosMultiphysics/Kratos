//
//   Project Name:        KratosFluidDynamicsApplication $
//   Last modified by:    $Author:               AFranci $
//   Date:                $Date:                May 2021 $
//   Revision:            $Revision:                 0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_elements/two_step_updated_lagrangian_V_P_implicit_fluid_PSPG_element.h"
#include "includes/cfd_variables.h"
#include <cmath>

namespace Kratos
{

  /*
   * public TwoStepUpdatedLagrangianVPImplicitFluidPspgElement<TDim> functions
   */

  template <unsigned int TDim>
  Element::Pointer TwoStepUpdatedLagrangianVPImplicitFluidPspgElement<TDim>::Clone(IndexType NewId, NodesArrayType const &rThisNodes) const
  {

    TwoStepUpdatedLagrangianVPImplicitFluidPspgElement NewElement(NewId, this->GetGeometry().Create(rThisNodes), this->pGetProperties());

    NewElement.SetData(this->GetData());
    NewElement.SetFlags(this->GetFlags());

    return Element::Pointer(new TwoStepUpdatedLagrangianVPImplicitFluidPspgElement(NewElement));
  }

  template <unsigned int TDim>
  int TwoStepUpdatedLagrangianVPImplicitFluidPspgElement<TDim>::Check(const ProcessInfo &rCurrentProcessInfo) const
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

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for (unsigned int i = 0; i < this->GetGeometry().size(); ++i)
    {
      if (this->GetGeometry()[i].SolutionStepsDataHas(VELOCITY) == false)
        KRATOS_THROW_ERROR(std::invalid_argument, "missing VELOCITY variable on solution step data for node ",
                           this->GetGeometry()[i].Id());
      if (this->GetGeometry()[i].SolutionStepsDataHas(PRESSURE) == false)
        KRATOS_THROW_ERROR(std::invalid_argument, "missing PRESSURE variable on solution step data for node ",
                           this->GetGeometry()[i].Id());
      if (this->GetGeometry()[i].SolutionStepsDataHas(BODY_FORCE) == false)
        KRATOS_THROW_ERROR(std::invalid_argument, "missing BODY_FORCE variable on solution step data for node ",
                           this->GetGeometry()[i].Id());
      if (this->GetGeometry()[i].SolutionStepsDataHas(DENSITY) == false)
        KRATOS_THROW_ERROR(std::invalid_argument, "missing DENSITY variable on solution step data for node ",
                           this->GetGeometry()[i].Id());
      if (this->GetGeometry()[i].SolutionStepsDataHas(DYNAMIC_VISCOSITY) == false)
        KRATOS_THROW_ERROR(std::invalid_argument,
                           "missing DYNAMIC_VISCOSITY variable on solution step data for node ",
                           this->GetGeometry()[i].Id());
      if (this->GetGeometry()[i].HasDofFor(VELOCITY_X) == false ||
          this->GetGeometry()[i].HasDofFor(VELOCITY_Y) == false ||
          this->GetGeometry()[i].HasDofFor(VELOCITY_Z) == false)
        KRATOS_THROW_ERROR(std::invalid_argument, "missing VELOCITY component degree of freedom on node ",
                           this->GetGeometry()[i].Id());
      if (this->GetGeometry()[i].HasDofFor(PRESSURE) == false)
        KRATOS_THROW_ERROR(std::invalid_argument, "missing PRESSURE component degree of freedom on node ",
                           this->GetGeometry()[i].Id());
    }

    // If this is a 2D problem, check that nodes are in XY plane
    if (this->GetGeometry().WorkingSpaceDimension() == 2)
    {
      for (unsigned int i = 0; i < this->GetGeometry().size(); ++i)
      {
        if (this->GetGeometry()[i].Z() != 0.0)
          KRATOS_THROW_ERROR(std::invalid_argument,
                             "Node with non-zero Z coordinate found. Id: ", this->GetGeometry()[i].Id());
      }
    }

    // Consitutive law checks
    const auto &r_properties = this->GetProperties();
    const auto &r_geometry = this->GetGeometry();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();

    //WARNING THIS MUST BE REMOVED ASAP
    const_cast<TwoStepUpdatedLagrangianVPImplicitFluidPspgElement<TDim> *>(this)->mpConstitutiveLaw = const_cast<TwoStepUpdatedLagrangianVPImplicitFluidPspgElement<TDim> *>(this)->GetProperties().GetValue(CONSTITUTIVE_LAW);
    //mpConstitutiveLaw = this->GetProperties().GetValue(CONSTITUTIVE_LAW);

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
    return mpConstitutiveLaw->Check(r_properties, r_geometry, rCurrentProcessInfo);

    return ierr;

    KRATOS_CATCH("");
  }

  template <unsigned int TDim>
  void TwoStepUpdatedLagrangianVPImplicitFluidPspgElement<TDim>::CalculateLocalContinuityEqForPressure(MatrixType &rLeftHandSideMatrix,
                                                                                                       VectorType &rRightHandSideVector,
                                                                                                       const ProcessInfo &rCurrentProcessInfo)
  {

    GeometryType &rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();

    // Check sizes and initialize
    if (rLeftHandSideMatrix.size1() != NumNodes)
      rLeftHandSideMatrix.resize(NumNodes, NumNodes, false);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(NumNodes, NumNodes);

    if (rRightHandSideVector.size() != NumNodes)
      rRightHandSideVector.resize(NumNodes, false);

    noalias(rRightHandSideVector) = ZeroVector(NumNodes);

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

    double Density = this->mMaterialDensity;
    double VolumetricCoeff = this->mMaterialVolumetricCoefficient;
    double DeviatoricCoeff = this->mMaterialDeviatoricCoefficient;

    double Tau = 0;
    this->CalculateTauPSPG(Tau, ElemSize, Density, DeviatoricCoeff, rCurrentProcessInfo);

    double totalVolume = 0;
    bool computeElement = false;

    MatrixType DynamicStabilizationMatrix = ZeroMatrix(NumNodes, NumNodes);

    // Loop on integration points
    for (unsigned int g = 0; g < NumGauss; ++g)
    {
      const double GaussWeight = GaussWeights[g];
      totalVolume += GaussWeight;
      const ShapeFunctionsType &N = row(NContainer, g);
      const ShapeFunctionDerivativesType &rDN_DX = DN_DX[g];
      computeElement = this->CalcCompleteStrainRate(rElementalVariables, rCurrentProcessInfo, rDN_DX, theta);
      computeElement = true;
      if (computeElement == true && this->IsNot(BLOCKED) && this->IsNot(ISOLATED))
      {

        double StabilizedWeight = Tau * GaussWeight;
        this->ComputeStabLaplacianMatrix(rLeftHandSideMatrix, rDN_DX, StabilizedWeight);

        array_1d<double, TDim> OldPressureGradient = ZeroVector(TDim);
        this->EvaluateGradientInPoint(OldPressureGradient, PRESSURE, rDN_DX);

        for (SizeType i = 0; i < NumNodes; ++i)
        {
          // RHS contribution
          // Velocity divergence
          rRightHandSideVector[i] += GaussWeight * N[i] * rElementalVariables.VolumetricDefRate;

          // this->AddStabilizationNodalTermsRHS(rRightHandSideVector, Tau, Density, GaussWeight, rDN_DX, i);

          this->AddPspgDynamicPartStabilization(rRightHandSideVector, Tau, Density, GaussWeight, TimeStep, rDN_DX, N, i);

          double laplacianRHSi = 0;
          double bodyForceStabilizedRHSi = 0;
          array_1d<double, 3> &VolumeAcceleration = this->GetGeometry()[i].FastGetSolutionStepValue(VOLUME_ACCELERATION);
          for (SizeType d = 0; d < TDim; ++d)
          {
            laplacianRHSi += StabilizedWeight * rDN_DX(i, d) * OldPressureGradient[d];

            bodyForceStabilizedRHSi += StabilizedWeight * rDN_DX(i, d) * (Density * VolumeAcceleration[d]);
          }
          rRightHandSideVector[i] += -laplacianRHSi - bodyForceStabilizedRHSi;

          // for (SizeType j = 0; j < NumNodes; ++j)
          // {
          //   double Tij = 0.0;
          //   for (SizeType d = 0; d < TDim; ++d)
          //   {
          //     Tij += rDN_DX(i, d) * N[j];
          //   }
          //   DynamicStabilizationMatrix(i, j) += StabilizedWeight * Density * Tij;
          // }
        }
      }
    }

    if (computeElement == true && this->IsNot(BLOCKED) && this->IsNot(ISOLATED))
    {

      VectorType PressureValues = ZeroVector(NumNodes);
      VectorType PressureValuesForRHS = ZeroVector(NumNodes);
      this->GetPressureValues(PressureValuesForRHS, 0);

      VectorType AccelerationValues = ZeroVector(NumNodes);
      this->GetAccelerationValues(AccelerationValues, 0);

      // noalias(rRightHandSideVector) += prod(DynamicStabilizationMatrix, AccelerationValues);

      // the LHS matrix up to now just contains the laplacian term and the bound term
      //  noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, PressureValuesForRHS);

      this->GetPressureValues(PressureValues, 1);
      noalias(PressureValuesForRHS) += -PressureValues;
      MatrixType BulkMatrix = ZeroMatrix(NumNodes, NumNodes);
      double lumpedBulkCoeff = totalVolume / (VolumetricCoeff);

      this->ComputeBulkMatrixLump(BulkMatrix, lumpedBulkCoeff);
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
  void TwoStepUpdatedLagrangianVPImplicitFluidPspgElement<TDim>::CalculateTauPSPG(double &Tau,
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
    this->CalcMeanVelocityNorm(MeanVelocity, 0);

    // // Tau Fic
    // Tau = (ElemSize * ElemSize * DeltaTime) / (Density * MeanVelocity * DeltaTime * ElemSize + Density * ElemSize * ElemSize + 8.0 * Viscosity * DeltaTime);

    // // Tau Tezduyar first proposal
    // double timeScale = 0.5 * (ElemSize / (2 * MeanVelocity) + DeltaTime);
    // double Reynolds = MeanVelocity * ElemSize * Density / (2.0 * Viscosity);
    // double ReynoldsFactor = 1.0;
    // if (Reynolds < 3)
    // {
    //   ReynoldsFactor = (ElemSize * ReynoldsFactor / (2 * MeanVelocity)) / 3.0;
    // }
    // Tau = timeScale * ReynoldsFactor / Density;

    // // Tau Tezduyar second proposal
    // double timeStepTerm = 2.0 / DeltaTime;
    // double timeScaleTerm = 2.0 * MeanVelocity / ElemSize;
    // double viscousTerm = 3.0 * (4 * Viscosity / (Density * pow(ElemSize, 2)));
    // // std::cout << "timeStepTerm= "<<timeStepTerm<< " timeScaleTerm= "<<timeScaleTerm<< " viscousTerm= "<<viscousTerm << std::endl;
    // double alternativeTau = pow(timeStepTerm, 2) + pow(timeScaleTerm, 2) + pow(viscousTerm, 2);
    // Tau = 1.0 / (std::sqrt(alternativeTau) * Density);

    // // Simplest Tau
    Tau = DeltaTime / Density;

    // std::cout << "Tau= " << Tau << std::endl;

    const double tolerance = 1.0e-13;
    if (MeanVelocity < tolerance)
    {
      Tau = 0;
    }
  }

  template <unsigned int TDim>
  void TwoStepUpdatedLagrangianVPImplicitFluidPspgElement<TDim>::ComputeStabLaplacianMatrix(MatrixType &StabLaplacianMatrix,
                                                                                            const ShapeFunctionDerivativesType &rDN_DX,
                                                                                            const double Weight)

  {
    // LHS contribution
    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    for (SizeType i = 0; i < NumNodes; ++i)
    {
      for (SizeType j = 0; j < NumNodes; ++j)
      {
        double Lij = 0.0;
        for (SizeType d = 0; d < TDim; ++d)
        {
          Lij += rDN_DX(i, d) * rDN_DX(j, d);
        }
        StabLaplacianMatrix(i, j) += Weight * Lij;
      }
    }
  }

  template <unsigned int TDim>
  void TwoStepUpdatedLagrangianVPImplicitFluidPspgElement<TDim>::AddStabilizationNodalTermsLHS(MatrixType &rLeftHandSideMatrix,
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
  void TwoStepUpdatedLagrangianVPImplicitFluidPspgElement<TDim>::AddStabilizationNodalTermsRHS(VectorType &rRightHandSideVector,
                                                                                               const double Tau,
                                                                                               const double Density,
                                                                                               const double Weight,
                                                                                               const ShapeFunctionDerivativesType &rDN_DX,
                                                                                               const SizeType i)
  {

    double RHSi = 0;
    if (this->GetGeometry()[i].SolutionStepsDataHas(VOLUME_ACCELERATION))
    { // it must be checked once at the begining only
      array_1d<double, 3> &VolumeAcceleration = this->GetGeometry()[i].FastGetSolutionStepValue(VOLUME_ACCELERATION);

      for (SizeType d = 0; d < TDim; ++d)
      {
        RHSi += -rDN_DX(i, d) * Tau * (Density * VolumeAcceleration[d]);
      }
    }
    rRightHandSideVector[i] += Weight * RHSi;
  }

  template <>
  void TwoStepUpdatedLagrangianVPImplicitFluidPspgElement<2>::AddPspgDynamicPartStabilization(VectorType &rRightHandSideVector,
                                                                                              const double Tau,
                                                                                              const double Density,
                                                                                              const double Weight,
                                                                                              const double TimeStep,
                                                                                              const ShapeFunctionDerivativesType &rDN_DX,
                                                                                              const ShapeFunctionsType &rN,
                                                                                              const SizeType i)
  {

    double RHSi = 0;
    // LHS contribution
    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    for (SizeType j = 0; j < NumNodes; ++j)
    {
      // double acc_X = (this->GetGeometry()[j].FastGetSolutionStepValue(VELOCITY_X,0) - this->GetGeometry()[j].FastGetSolutionStepValue(VELOCITY_X,1)) / TimeStep;
      // double acc_Y = (this->GetGeometry()[j].FastGetSolutionStepValue(VELOCITY_Y,0) - this->GetGeometry()[j].FastGetSolutionStepValue(VELOCITY_Y,1)) / TimeStep;

      RHSi += rDN_DX(i, 0) * rN[j] * this->GetGeometry()[j].FastGetSolutionStepValue(ACCELERATION_X, 0) + rDN_DX(i, 1) * rN[j] * this->GetGeometry()[j].FastGetSolutionStepValue(ACCELERATION_Y, 0);
      // RHSi += Tau * Density * (rDN_DX(i, 0) + rDN_DX(i, 1)) * rN[j] * (acc_X + acc_Y) ;
      //  RHSi +=  rN[i] *(rDN_DX(j, 0)*this->GetGeometry()[j].FastGetSolutionStepValue(VELOCITY_X,0) + rDN_DX(j, 1)*this->GetGeometry()[j].FastGetSolutionStepValue(VELOCITY_Y,0));
    }
    rRightHandSideVector[i] += Weight * Tau * Density * RHSi;
  }

  template <>
  void TwoStepUpdatedLagrangianVPImplicitFluidPspgElement<3>::AddPspgDynamicPartStabilization(VectorType &rRightHandSideVector,
                                                                                              const double Tau,
                                                                                              const double Density,
                                                                                              const double Weight,
                                                                                              const double TimeStep,
                                                                                              const ShapeFunctionDerivativesType &rDN_DX,
                                                                                              const ShapeFunctionsType &rN,
                                                                                              const SizeType i)
  {

    double RHSi = 0;

    // LHS contribution
    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    for (SizeType j = 0; j < NumNodes; ++j)
    {
      double termX = rDN_DX(i, 0) * rN[j] * this->GetGeometry()[j].FastGetSolutionStepValue(ACCELERATION_X, 0);
      double termY = rDN_DX(i, 1) * rN[j] * this->GetGeometry()[j].FastGetSolutionStepValue(ACCELERATION_Y, 0);
      double termZ = rDN_DX(i, 2) * rN[j] * this->GetGeometry()[j].FastGetSolutionStepValue(ACCELERATION_Z, 0);
      RHSi += Tau * Density * (termX + termY + termZ);
    }
    rRightHandSideVector[i] += Weight * RHSi;
  }

  template class TwoStepUpdatedLagrangianVPImplicitFluidPspgElement<2>;
  template class TwoStepUpdatedLagrangianVPImplicitFluidPspgElement<3>;

} // namespace Kratos
