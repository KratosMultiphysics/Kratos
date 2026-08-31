//
//   Project Name:        KratosFluidDynamicsApplication $
//   Last modified by:    $Author:               AFranci $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_elements/two_step_updated_lagrangian_V_P_implicit_element.h"
#include "includes/cfd_variables.h"

namespace Kratos
{

  template <unsigned int TDim>
  Element::Pointer TwoStepUpdatedLagrangianVPImplicitElement<TDim>::Clone(IndexType NewId, NodesArrayType const &rThisNodes) const
  {
    KRATOS_TRY;

    TwoStepUpdatedLagrangianVPImplicitElement NewElement(NewId, this->GetGeometry().Create(rThisNodes), this->pGetProperties());
    return Element::Pointer(new TwoStepUpdatedLagrangianVPImplicitElement(NewElement));

    KRATOS_CATCH("");
  }

  template <unsigned int TDim>
  void TwoStepUpdatedLagrangianVPImplicitElement<TDim>::CalculateLocalSystem(MatrixType &rLeftHandSideMatrix,
                                                                             VectorType &rRightHandSideVector,
                                                                             const ProcessInfo &rCurrentProcessInfo)
  {
    KRATOS_TRY;

    switch (rCurrentProcessInfo[FRACTIONAL_STEP])
    {
    case 1:
    {
      this->CalculateLocalMomentumEquations(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
      break;
    }
    case 5:
    {
      this->CalculateLocalContinuityEqForPressure(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
      break;
    }

    default:
    {
      KRATOS_THROW_ERROR(std::logic_error, "Unexpected value for TWO_STEP_UPDATED_LAGRANGIAN_V_P_ELEMENT index: ", rCurrentProcessInfo[FRACTIONAL_STEP]);
    }
    }

    KRATOS_CATCH("");
  }

  template <unsigned int TDim>
  void TwoStepUpdatedLagrangianVPImplicitElement<TDim>::CalculateLocalMomentumEquations(MatrixType &rLeftHandSideMatrix,
                                                                                        VectorType &rRightHandSideVector,
                                                                                        const ProcessInfo &rCurrentProcessInfo)
  {
    KRATOS_TRY;

    const GeometryType &rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = TDim * NumNodes;

    MatrixType MassMatrix = ZeroMatrix(LocalSize, LocalSize);
    MatrixType StiffnessMatrix = ZeroMatrix(LocalSize, LocalSize);

    // Check sizes and initialize
    if (rLeftHandSideMatrix.size1() != LocalSize)
      rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);

    if (rRightHandSideVector.size() != LocalSize)
      rRightHandSideVector.resize(LocalSize, false);

    noalias(rRightHandSideVector) = ZeroVector(LocalSize);

    // Shape functions and integration points
    ShapeFunctionDerivativesArrayType DN_DX;
    Matrix NContainer;
    VectorType GaussWeights;
    this->CalculateGeometryData(DN_DX, NContainer, GaussWeights);
    const SizeType NumGauss = GaussWeights.size();
    const double TimeStep = rCurrentProcessInfo[DELTA_TIME];

    const double theta = this->GetThetaMomentum();

    ElementalVariables rElementalVariables;
    this->InitializeElementalVariables(rElementalVariables);

    double totalVolume = 0.0;
    double MeanValueMass = 0.0;
    double Density = 0.0;
    double DeviatoricCoeff = 0.0;
    double VolumetricCoeff = 0.0;
    bool computeElement = false;
    // Loop on integration points
    for (SizeType g = 0; g < NumGauss; ++g)
    {
      const double GaussWeight = GaussWeights[g];
      totalVolume += GaussWeight;
      const ShapeFunctionsType &N = row(NContainer, g);
      const ShapeFunctionDerivativesType &rDN_DX = DN_DX[g];

      double Pressure = 0.0;
      double OldPressure = 0.0;

      this->EvaluateInPoint(Pressure, PRESSURE, N, 0);

      this->EvaluateInPoint(OldPressure, PRESSURE, N, 1);

      rElementalVariables.MeanPressure = OldPressure * (1 - theta) + Pressure * theta;

      computeElement = this->CalcMechanicsUpdated(rElementalVariables, rCurrentProcessInfo, rDN_DX);

      this->CalcElasticPlasticCauchySplitted(rElementalVariables, g, N, rCurrentProcessInfo, Density, DeviatoricCoeff, VolumetricCoeff);

      if (computeElement && this->IsNot(BLOCKED) && this->IsNot(ISOLATED))
      {
        this->AddExternalForces(rRightHandSideVector, Density, N, GaussWeight);

        this->AddInternalForces(rRightHandSideVector, rDN_DX, rElementalVariables, GaussWeight);

        this->ComputeCompleteTangentTerm(rElementalVariables, StiffnessMatrix, rDN_DX, DeviatoricCoeff, VolumetricCoeff, theta, GaussWeight);
      }
      // else
      // {
      //   std::cout << "in the momentum BLOCKED ELEMENT: " << this->Id() << std::endl;
      //   mpConstitutiveLaw = this->GetProperties().GetValue(CONSTITUTIVE_LAW);
      //   auto constitutive_law_values = ConstitutiveLaw::Parameters(this->GetGeometry(), this->GetProperties(), rCurrentProcessInfo);
      //   Density = mpConstitutiveLaw->CalculateValue(constitutive_law_values, DENSITY, Density);
      // }
    }

    const double lumpedDynamicWeight = totalVolume * Density;
    // FIXME: This has to be computed with the corresponding shape functions...
    this->ComputeLumpedMassMatrix(MassMatrix, lumpedDynamicWeight, MeanValueMass);
    if (computeElement && this->IsNot(BLOCKED) && this->IsNot(ISOLATED))
    {
      double BulkReductionCoefficient = 1.0;
      double MeanValueStiffness = 0.0;
      this->ComputeBulkReductionCoefficient(MassMatrix, StiffnessMatrix, MeanValueStiffness, BulkReductionCoefficient, TimeStep);
      if (BulkReductionCoefficient != 1.0)
      {
        // VolumetricCoeff*=BulkReductionCoefficient;
        VolumetricCoeff *= MeanValueMass * 2.0 / (TimeStep * MeanValueStiffness);
        noalias(StiffnessMatrix) = ZeroMatrix(LocalSize, LocalSize);

        for (SizeType g = 0; g < NumGauss; ++g)
        {
          const double GaussWeight = GaussWeights[g];
          const ShapeFunctionDerivativesType &rDN_DX = DN_DX[g];
          this->ComputeCompleteTangentTerm(rElementalVariables, StiffnessMatrix, rDN_DX, DeviatoricCoeff, VolumetricCoeff, theta, GaussWeight);
        }
      }
    }
    // Add residual of previous iteration to RHS
    VectorType VelocityValues = ZeroVector(LocalSize);
    VectorType AccelerationValues = ZeroVector(LocalSize);

    // //1st order
    // this->GetVelocityValues(VelocityValues,0);
    // AccelerationValues = VelocityValues/TimeStep;
    // this->GetAccelerationValues(LastAccValues,0);
    // this->GetVelocityValues(VelocityValues,1);
    // AccelerationValues += -VelocityValues/TimeStep;
    // noalias( rRightHandSideVector ) += -prod(MassMatrix,AccelerationValues);
    // noalias( rLeftHandSideMatrix ) +=  MassMatrix/TimeStep;

    // 2nd order
    this->GetAccelerationValues(AccelerationValues, 0);
    this->GetVelocityValues(VelocityValues, 0);
    noalias(AccelerationValues) += -2.0 * VelocityValues / TimeStep;
    this->GetVelocityValues(VelocityValues, 1);
    noalias(AccelerationValues) += 2.0 * VelocityValues / TimeStep; // these are negative accelerations
    noalias(rRightHandSideVector) += prod(MassMatrix, AccelerationValues);
    noalias(rLeftHandSideMatrix) += StiffnessMatrix + MassMatrix * 2 / TimeStep;

    // // Add residual of previous iteration to RHS
    // VectorType VelocityValues = ZeroVector(LocalSize);
    // VectorType UpdatedAccelerations = ZeroVector(LocalSize);
    // VectorType LastAccValues = ZeroVector(LocalSize);

    // // //1st order
    // // this->GetVelocityValues(VelocityValues,0);
    // // UpdatedAccelerations = VelocityValues/TimeStep;
    // // this->GetAccelerationValues(LastAccValues,0);
    // // this->GetVelocityValues(VelocityValues,1);
    // // UpdatedAccelerations += -VelocityValues/TimeStep;
    // // // UpdatedAccelerations =LastAccValues;
    // // noalias( rRightHandSideVector ) += -prod(MassMatrix,UpdatedAccelerations);
    // // noalias( rLeftHandSideMatrix ) +=  MassMatrix/TimeStep;

    // //2nd order
    // this->GetVelocityValues(VelocityValues,0);
    // UpdatedAccelerations = 2.0*VelocityValues/TimeStep;
    // this->GetAccelerationValues(LastAccValues,0);
    // this->GetVelocityValues(VelocityValues,1);
    // UpdatedAccelerations += -2.0*VelocityValues/TimeStep - LastAccValues;
    // noalias( rRightHandSideVector ) += -prod(MassMatrix,UpdatedAccelerations);
    // noalias( rLeftHandSideMatrix ) +=  StiffnessMatrix;
    // noalias( rLeftHandSideMatrix ) +=  MassMatrix*2/TimeStep;

    KRATOS_CATCH("");
  }

  template <unsigned int TDim>
  void TwoStepUpdatedLagrangianVPImplicitElement<TDim>::CalculateOnIntegrationPoints(const Variable<bool> &rVariable,
                                                                                     std::vector<bool> &rOutput,
                                                                                     const ProcessInfo &rCurrentProcessInfo)
  {
    if (rVariable == YIELDED)
    {
      rOutput[0] = this->GetValue(YIELDED);
    }
  }

  template <unsigned int TDim>
  void TwoStepUpdatedLagrangianVPImplicitElement<TDim>::CalculateOnIntegrationPoints(const Variable<double> &rVariable,
                                                                                     std::vector<double> &rOutput,
                                                                                     const ProcessInfo &rCurrentProcessInfo)
  {
    if (rVariable == EQ_STRAIN_RATE)
    {
      rOutput[0] = this->GetValue(EQ_STRAIN_RATE);
    }
    if (rVariable == MECHANICAL_DISSIPATION)
    {
      rOutput[0] = this->GetValue(MECHANICAL_DISSIPATION);
    }
  }

  template <unsigned int TDim>
  void TwoStepUpdatedLagrangianVPImplicitElement<TDim>::CalculateOnIntegrationPoints(const Variable<Vector> &rVariable,
                                                                                     std::vector<Vector> &rOutput,
                                                                                     const ProcessInfo &rCurrentProcessInfo)
  {
    if (rVariable == CAUCHY_STRESS_VECTOR)
    {
      rOutput[0] = this->GetValue(CAUCHY_STRESS_VECTOR);
    }
  }

  template <>
  void TwoStepUpdatedLagrangianVPImplicitElement<2>::ComputeCompleteTangentTerm(ElementalVariables &rElementalVariables,
                                                                                MatrixType &rDampingMatrix,
                                                                                const ShapeFunctionDerivativesType &rDN_DX,
                                                                                const double secondLame,
                                                                                const double bulkModulus,
                                                                                const double theta,
                                                                                const double Weight)
  {
    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    constexpr double FourThirds = 4.0 / 3.0;
    constexpr double nTwoThirds = -2.0 / 3.0;

    SizeType FirstRow = 0;
    SizeType FirstCol = 0;

    for (SizeType j = 0; j < NumNodes; ++j)
    {
      for (SizeType i = 0; i < NumNodes; ++i)
      {
        const double lagDNXi = rDN_DX(i, 0) * rElementalVariables.InvFgrad(0, 0) + rDN_DX(i, 1) * rElementalVariables.InvFgrad(1, 0);
        const double lagDNYi = rDN_DX(i, 0) * rElementalVariables.InvFgrad(0, 1) + rDN_DX(i, 1) * rElementalVariables.InvFgrad(1, 1);
        const double lagDNXj = rDN_DX(j, 0) * rElementalVariables.InvFgrad(0, 0) + rDN_DX(j, 1) * rElementalVariables.InvFgrad(1, 0);
        const double lagDNYj = rDN_DX(j, 0) * rElementalVariables.InvFgrad(0, 1) + rDN_DX(j, 1) * rElementalVariables.InvFgrad(1, 1);

        // First Row
        rDampingMatrix(FirstRow, FirstCol) += Weight * ((FourThirds * secondLame + bulkModulus) * lagDNXi * lagDNXj + lagDNYi * lagDNYj * secondLame) * theta;
        rDampingMatrix(FirstRow, FirstCol + 1) += Weight * ((nTwoThirds * secondLame + bulkModulus) * lagDNXi * lagDNYj + lagDNYi * lagDNXj * secondLame) * theta;

        // Second Row
        rDampingMatrix(FirstRow + 1, FirstCol) += Weight * ((nTwoThirds * secondLame + bulkModulus) * lagDNYi * lagDNXj + lagDNXi * lagDNYj * secondLame) * theta;
        rDampingMatrix(FirstRow + 1, FirstCol + 1) += Weight * ((FourThirds * secondLame + bulkModulus) * lagDNYi * lagDNYj + lagDNXi * lagDNXj * secondLame) * theta;

        // Update Counter
        FirstRow += 2;
      }
      FirstRow = 0;
      FirstCol += 2;
    }
  }

  template <>
  void TwoStepUpdatedLagrangianVPImplicitElement<3>::ComputeCompleteTangentTerm(ElementalVariables &rElementalVariables,
                                                                                MatrixType &rDampingMatrix,
                                                                                const ShapeFunctionDerivativesType &rDN_DX,
                                                                                const double secondLame,
                                                                                const double bulkModulus,
                                                                                const double theta,
                                                                                const double Weight)
  {

    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    constexpr double FourThirds = 4.0 / 3.0;
    constexpr double nTwoThirds = -2.0 / 3.0;

    SizeType FirstRow = 0;
    SizeType FirstCol = 0;

    for (SizeType j = 0; j < NumNodes; ++j)
    {
      for (SizeType i = 0; i < NumNodes; ++i)
      {
        const double lagDNXi = rDN_DX(i, 0) * rElementalVariables.InvFgrad(0, 0) + rDN_DX(i, 1) * rElementalVariables.InvFgrad(1, 0) + rDN_DX(i, 2) * rElementalVariables.InvFgrad(2, 0);
        const double lagDNYi = rDN_DX(i, 0) * rElementalVariables.InvFgrad(0, 1) + rDN_DX(i, 1) * rElementalVariables.InvFgrad(1, 1) + rDN_DX(i, 2) * rElementalVariables.InvFgrad(2, 1);
        const double lagDNZi = rDN_DX(i, 0) * rElementalVariables.InvFgrad(0, 2) + rDN_DX(i, 1) * rElementalVariables.InvFgrad(1, 2) + rDN_DX(i, 2) * rElementalVariables.InvFgrad(2, 2);
        const double lagDNXj = rDN_DX(j, 0) * rElementalVariables.InvFgrad(0, 0) + rDN_DX(j, 1) * rElementalVariables.InvFgrad(1, 0) + rDN_DX(j, 2) * rElementalVariables.InvFgrad(2, 0);
        const double lagDNYj = rDN_DX(j, 0) * rElementalVariables.InvFgrad(0, 1) + rDN_DX(j, 1) * rElementalVariables.InvFgrad(1, 1) + rDN_DX(j, 2) * rElementalVariables.InvFgrad(2, 1);
        const double lagDNZj = rDN_DX(j, 0) * rElementalVariables.InvFgrad(0, 2) + rDN_DX(j, 1) * rElementalVariables.InvFgrad(1, 2) + rDN_DX(j, 2) * rElementalVariables.InvFgrad(2, 2);

        // First Row
        rDampingMatrix(FirstRow, FirstCol) += Weight * ((FourThirds * secondLame + bulkModulus) * lagDNXi * lagDNXj + (lagDNYi * lagDNYj + lagDNZi * lagDNZj) * secondLame) * theta;
        rDampingMatrix(FirstRow, FirstCol + 1) += Weight * ((nTwoThirds * secondLame + bulkModulus) * lagDNXi * lagDNYj + lagDNYi * lagDNXj * secondLame) * theta;
        rDampingMatrix(FirstRow, FirstCol + 2) += Weight * ((nTwoThirds * secondLame + bulkModulus) * lagDNXi * lagDNZj + lagDNZi * lagDNXj * secondLame) * theta;

        // Second Row
        rDampingMatrix(FirstRow + 1, FirstCol) += Weight * ((nTwoThirds * secondLame + bulkModulus) * lagDNYi * lagDNXj + lagDNXi * lagDNYj * secondLame) * theta;
        rDampingMatrix(FirstRow + 1, FirstCol + 1) += Weight * ((FourThirds * secondLame + bulkModulus) * lagDNYi * lagDNYj + (lagDNXi * lagDNXj + lagDNZi * lagDNZj) * secondLame) * theta;
        rDampingMatrix(FirstRow + 1, FirstCol + 2) += Weight * ((nTwoThirds * secondLame + bulkModulus) * lagDNYi * lagDNZj + lagDNZi * lagDNYj * secondLame) * theta;

        // Third Row
        rDampingMatrix(FirstRow + 2, FirstCol) += Weight * ((nTwoThirds * secondLame + bulkModulus) * lagDNZi * lagDNXj + lagDNXi * lagDNZj * secondLame) * theta;
        rDampingMatrix(FirstRow + 2, FirstCol + 1) += Weight * ((nTwoThirds * secondLame + bulkModulus) * lagDNZi * lagDNYj + lagDNYi * lagDNZj * secondLame) * theta;
        rDampingMatrix(FirstRow + 2, FirstCol + 2) += Weight * ((FourThirds * secondLame + bulkModulus) * lagDNZi * lagDNZj + (lagDNXi * lagDNXj + lagDNYi * lagDNYj) * secondLame) * theta;

        // Update Counter
        FirstRow += 3;
      }
      FirstRow = 0;
      FirstCol += 3;
    }
  }

  template <unsigned int TDim>
  int TwoStepUpdatedLagrangianVPImplicitElement<TDim>::Check(const ProcessInfo &rCurrentProcessInfo) const
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

    return ierr;

    KRATOS_CATCH("");
  }

  /*
   * Template class definition (this should allow us to compile the desired template instantiations)
   */

  template class TwoStepUpdatedLagrangianVPImplicitElement<2>;
  template class TwoStepUpdatedLagrangianVPImplicitElement<3>;

} // namespace Kratos
