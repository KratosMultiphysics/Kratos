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
#include "custom_elements/two_step_updated_lagrangian_V_P_implicit_nodally_integrated_solid_element.h"
#include "includes/cfd_variables.h"

namespace Kratos
{

  /*
   * public TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedSolidElement<TDim> functions
   */

  template <unsigned int TDim>
  Element::Pointer TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedSolidElement<TDim>::Clone(IndexType NewId, NodesArrayType const &rThisNodes) const
  {
    // return Element::Pointer( BaseType::Clone(NewId,rThisNodes) );
    TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedSolidElement NewElement(NewId, this->GetGeometry().Create(rThisNodes), this->pGetProperties());

    if (NewElement.mCurrentTotalCauchyStress.size() != this->mCurrentTotalCauchyStress.size())
      NewElement.mCurrentTotalCauchyStress.resize(this->mCurrentTotalCauchyStress.size());

    for (unsigned int i = 0; i < this->mCurrentTotalCauchyStress.size(); i++)
    {
      NewElement.mCurrentTotalCauchyStress[i] = this->mCurrentTotalCauchyStress[i];
    }

    if (NewElement.mCurrentDeviatoricCauchyStress.size() != this->mCurrentDeviatoricCauchyStress.size())
      NewElement.mCurrentDeviatoricCauchyStress.resize(this->mCurrentDeviatoricCauchyStress.size());

    for (unsigned int i = 0; i < this->mCurrentDeviatoricCauchyStress.size(); i++)
    {
      NewElement.mCurrentDeviatoricCauchyStress[i] = this->mCurrentDeviatoricCauchyStress[i];
    }

    if (NewElement.mUpdatedTotalCauchyStress.size() != this->mUpdatedTotalCauchyStress.size())
      NewElement.mUpdatedTotalCauchyStress.resize(this->mUpdatedTotalCauchyStress.size());

    for (unsigned int i = 0; i < this->mUpdatedTotalCauchyStress.size(); i++)
    {
      NewElement.mUpdatedTotalCauchyStress[i] = this->mUpdatedTotalCauchyStress[i];
    }

    if (NewElement.mUpdatedDeviatoricCauchyStress.size() != this->mUpdatedDeviatoricCauchyStress.size())
      NewElement.mUpdatedDeviatoricCauchyStress.resize(this->mUpdatedDeviatoricCauchyStress.size());

    for (unsigned int i = 0; i < this->mUpdatedDeviatoricCauchyStress.size(); i++)
    {
      NewElement.mUpdatedDeviatoricCauchyStress[i] = this->mUpdatedDeviatoricCauchyStress[i];
    }

    NewElement.SetData(this->GetData());
    NewElement.SetFlags(this->GetFlags());

    return Element::Pointer(new TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedSolidElement(NewElement));
  }

  template <unsigned int TDim>
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedSolidElement<TDim>::Initialize()
  {
    KRATOS_TRY;

    // LargeDisplacementElement::Initialize();
    const GeometryType &rGeom = this->GetGeometry();
    SizeType integration_points_number = rGeom.IntegrationPointsNumber(GeometryData::GI_GAUSS_1);
    // SizeType integration_points_number = rGeom.IntegrationPointsNumber(GeometryData::GI_GAUSS_4);
    if (this->mCurrentTotalCauchyStress.size() != integration_points_number)
      this->mCurrentTotalCauchyStress.resize(integration_points_number);

    if (this->mCurrentDeviatoricCauchyStress.size() != integration_points_number)
      this->mCurrentDeviatoricCauchyStress.resize(integration_points_number);

    if (this->mUpdatedTotalCauchyStress.size() != integration_points_number)
      this->mUpdatedTotalCauchyStress.resize(integration_points_number);

    if (this->mUpdatedDeviatoricCauchyStress.size() != integration_points_number)
      this->mUpdatedDeviatoricCauchyStress.resize(integration_points_number);

    unsigned int voigtsize = 3;
    if (TDim == 3)
    {
      voigtsize = 6;
    }
    for (unsigned int PointNumber = 0; PointNumber < integration_points_number; PointNumber++)
    {
      this->mCurrentTotalCauchyStress[PointNumber] = ZeroVector(voigtsize);
      this->mCurrentDeviatoricCauchyStress[PointNumber] = ZeroVector(voigtsize);
      this->mUpdatedTotalCauchyStress[PointNumber] = ZeroVector(voigtsize);
      this->mUpdatedDeviatoricCauchyStress[PointNumber] = ZeroVector(voigtsize);
    }

    KRATOS_CATCH("");
  }

  template <unsigned int TDim>
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedSolidElement<TDim>::UpdateCauchyStress(unsigned int g,
                                                                                                 ProcessInfo &rCurrentProcessInfo)
  {
    // double theta = this->GetThetaContinuity();
    double theta = 1.0;
    ElementalVariables rElementalVariables;
    this->InitializeElementalVariables(rElementalVariables);
    ShapeFunctionDerivativesArrayType DN_DX;
    Matrix NContainer;
    VectorType GaussWeights;
    this->CalculateGeometryData(DN_DX, NContainer, GaussWeights);
    const ShapeFunctionDerivativesType &rDN_DX = DN_DX[g];
    // bool computeElement=this->CalcStrainRate(rElementalVariables,rCurrentProcessInfo,rDN_DX,theta);
    bool computeElement = this->CalcCompleteStrainRate(rElementalVariables, rCurrentProcessInfo, rDN_DX, theta);
    const double TimeStep = rCurrentProcessInfo[DELTA_TIME];

    if (computeElement == true)
    {
      double Density = 0;
      double DeviatoricCoeff = 0;
      double VolumetricCoeff = 0;
      CalcElasticPlasticCauchySplitted(rElementalVariables, TimeStep, g, rCurrentProcessInfo, Density,
                                       DeviatoricCoeff, VolumetricCoeff);
    }

    this->mCurrentTotalCauchyStress[g] = this->mUpdatedTotalCauchyStress[g];
    this->mCurrentDeviatoricCauchyStress[g] = this->mUpdatedDeviatoricCauchyStress[g];
  }

  template <unsigned int TDim>
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedSolidElement<TDim>::InitializeSolutionStep(ProcessInfo &rCurrentProcessInfo)
  {
    KRATOS_TRY;
    const GeometryType &rGeom = this->GetGeometry();
    SizeType integration_points_number = rGeom.IntegrationPointsNumber(GeometryData::GI_GAUSS_1);
    // SizeType integration_points_number = rGeom.IntegrationPointsNumber(GeometryData::GI_GAUSS_4);

    for (unsigned int PointNumber = 0; PointNumber < integration_points_number; PointNumber++)
    {
      this->UpdateCauchyStress(PointNumber, rCurrentProcessInfo);
    }

    KRATOS_CATCH("");
  }

  template <unsigned int TDim>
  int TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedSolidElement<TDim>::Check(const ProcessInfo &rCurrentProcessInfo)
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

  template <unsigned int TDim>
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedSolidElement<TDim>::InitializeElementalVariables(ElementalVariables &rElementalVariables)
  {
    unsigned int voigtsize = 3;
    if (TDim == 3)
    {
      voigtsize = 6;
    }
    rElementalVariables.voigtsize = voigtsize;
    rElementalVariables.ConstitutiveMatrix = ZeroMatrix(voigtsize, voigtsize);
    rElementalVariables.DetFgrad = 1;
    rElementalVariables.DetFgradVel = 1;
    rElementalVariables.DeviatoricInvariant = 1;
    rElementalVariables.EquivalentStrainRate = 1;
    rElementalVariables.VolumetricDefRate = 1;
    rElementalVariables.SpatialDefRate.resize(voigtsize);
    rElementalVariables.MDGreenLagrangeMaterial.resize(voigtsize);
    rElementalVariables.Fgrad.resize(TDim, TDim);
    rElementalVariables.InvFgrad.resize(TDim, TDim);
    rElementalVariables.FgradVel.resize(TDim, TDim);
    rElementalVariables.InvFgradVel.resize(TDim, TDim);
    rElementalVariables.SpatialVelocityGrad.resize(TDim, TDim);

    rElementalVariables.MeanPressure = 0;
    rElementalVariables.CurrentTotalCauchyStress.resize(voigtsize);
    rElementalVariables.UpdatedTotalCauchyStress.resize(voigtsize);
    rElementalVariables.CurrentDeviatoricCauchyStress.resize(voigtsize);
    rElementalVariables.UpdatedDeviatoricCauchyStress.resize(voigtsize);
  }

  template <unsigned int TDim>
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedSolidElement<TDim>::CalculateLocalMomentumEquations(
      MatrixType &rLeftHandSideMatrix,
      VectorType &rRightHandSideVector,
      ProcessInfo &rCurrentProcessInfo)
  {
    KRATOS_TRY;

    GeometryType &rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();
    const unsigned int LocalSize = TDim * NumNodes;

    MatrixType MassMatrix = ZeroMatrix(LocalSize, LocalSize);
    MatrixType StiffnessMatrix = ZeroMatrix(LocalSize, LocalSize);

    // Check sizes and initialize
    if (rLeftHandSideMatrix.size1() != LocalSize)
      rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);

    rLeftHandSideMatrix = ZeroMatrix(LocalSize, LocalSize);

    if (rRightHandSideVector.size() != LocalSize)
      rRightHandSideVector.resize(LocalSize);

    rRightHandSideVector = ZeroVector(LocalSize);

    // Shape functions and integration points
    ShapeFunctionDerivativesArrayType DN_DX;
    Matrix NContainer;
    VectorType GaussWeights;
    this->CalculateGeometryData(DN_DX, NContainer, GaussWeights);
    const unsigned int NumGauss = GaussWeights.size();
    const double TimeStep = rCurrentProcessInfo[DELTA_TIME];

    double theta = 1.0;
    ElementalVariables rElementalVariables;
    this->InitializeElementalVariables(rElementalVariables);

    double totalVolume = 0;
    double Density = 0.0;
    double DeviatoricCoeff = 0;
    double VolumetricCoeff = 0;

    // Loop on integration points
    for (unsigned int g = 0; g < NumGauss; ++g)
    {
      const double GaussWeight = GaussWeights[g];
      totalVolume += GaussWeight;
      const ShapeFunctionsType &N = row(NContainer, g);
      const ShapeFunctionDerivativesType &rDN_DX = DN_DX[g];

      // double Pressure = 0;
      // double OldPressure = 0;

      // this->EvaluateInPoint(Pressure, PRESSURE, N, 0);

      // this->EvaluateInPoint(OldPressure, PRESSURE, N, 1);

      // rElementalVariables.MeanPressure = OldPressure * (1 - theta) + Pressure * theta;

      bool computeElement = this->CalcCompleteStrainRate(rElementalVariables, rCurrentProcessInfo, rDN_DX, theta);

      CalcElasticPlasticCauchySplitted(rElementalVariables, TimeStep, g, rCurrentProcessInfo, Density, DeviatoricCoeff, VolumetricCoeff);

      if (computeElement == true)
      {
        // this->AddExternalForces(rRightHandSideVector, Density, N, GaussWeight);

        double hybridCoeff = 0.5; // half nodal - half elemental
        // hybridCoeff = 1.0; // fully elemental
        double reducedElementalWeight = GaussWeight * hybridCoeff;

        this->AddInternalForces(rRightHandSideVector, rDN_DX, rElementalVariables, reducedElementalWeight);

        this->ComputeCompleteTangentTerm(rElementalVariables, StiffnessMatrix, rDN_DX, DeviatoricCoeff, VolumetricCoeff, theta, reducedElementalWeight);
      }
    }
    noalias(rLeftHandSideMatrix) += StiffnessMatrix;

    KRATOS_CATCH("");
  }


  template <>
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedSolidElement<2>::CalcElasticPlasticCauchySplitted(
      ElementalVariables &rElementalVariables, double TimeStep, unsigned int g, const ProcessInfo &rCurrentProcessInfo,
      double &Density, double &DeviatoricCoeff, double &VolumetricCoeff)
  {

    mpConstitutiveLaw = this->GetProperties().GetValue(CONSTITUTIVE_LAW);
    auto constitutive_law_values =
        ConstitutiveLaw::Parameters(this->GetGeometry(), this->GetProperties(), rCurrentProcessInfo);

    Flags &constitutive_law_options = constitutive_law_values.GetOptions();
    constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

    rElementalVariables.CurrentTotalCauchyStress = this->mCurrentTotalCauchyStress[g];
    rElementalVariables.CurrentDeviatoricCauchyStress = this->mCurrentDeviatoricCauchyStress[g];

    const Vector &r_shape_functions = row((this->GetGeometry()).ShapeFunctionsValues(), g);
    constitutive_law_values.SetShapeFunctionsValues(r_shape_functions);
    constitutive_law_values.SetStrainVector(rElementalVariables.SpatialDefRate);
    constitutive_law_values.SetStressVector(rElementalVariables.CurrentDeviatoricCauchyStress);

    mpConstitutiveLaw->CalculateMaterialResponseCauchy(constitutive_law_values);

    Density = mpConstitutiveLaw->CalculateValue(constitutive_law_values, DENSITY, Density);

    double poisson_ratio = mpConstitutiveLaw->CalculateValue(constitutive_law_values, POISSON_RATIO, poisson_ratio);
    double young_modulus = mpConstitutiveLaw->CalculateValue(constitutive_law_values, YOUNG_MODULUS, young_modulus);
    const double time_step = rCurrentProcessInfo[DELTA_TIME];
    DeviatoricCoeff = time_step * young_modulus / (2.0 * (1 + poisson_ratio));
    VolumetricCoeff =
        time_step * poisson_ratio * young_modulus / ((1.0 + poisson_ratio) * (1.0 - 2.0 * poisson_ratio)) +
        2.0 / 3.0 * DeviatoricCoeff;

    const double current_first_lame = time_step * poisson_ratio * young_modulus / ((1.0 + poisson_ratio) * (1.0 - 2.0 * poisson_ratio));

    this->mMaterialDeviatoricCoefficient = DeviatoricCoeff;
    this->mMaterialVolumetricCoefficient = VolumetricCoeff;
    this->mMaterialDensity = Density;

    rElementalVariables.UpdatedDeviatoricCauchyStress[0] = rElementalVariables.CurrentDeviatoricCauchyStress[0];
    rElementalVariables.UpdatedDeviatoricCauchyStress[1] = rElementalVariables.CurrentDeviatoricCauchyStress[1];
    rElementalVariables.UpdatedDeviatoricCauchyStress[2] = rElementalVariables.CurrentDeviatoricCauchyStress[2];

    rElementalVariables.UpdatedTotalCauchyStress[0] = +current_first_lame * rElementalVariables.VolumetricDefRate +
                                                      2.0 * DeviatoricCoeff * rElementalVariables.SpatialDefRate[0] +
                                                      rElementalVariables.CurrentTotalCauchyStress[0];
    rElementalVariables.UpdatedTotalCauchyStress[1] = current_first_lame * rElementalVariables.VolumetricDefRate +
                                                      2.0 * DeviatoricCoeff * rElementalVariables.SpatialDefRate[1] +
                                                      rElementalVariables.CurrentTotalCauchyStress[1];
    rElementalVariables.UpdatedTotalCauchyStress[2] =
        2.0 * DeviatoricCoeff * rElementalVariables.SpatialDefRate[2] + rElementalVariables.CurrentTotalCauchyStress[2];

    this->mUpdatedTotalCauchyStress[g] = rElementalVariables.UpdatedTotalCauchyStress;
    this->mUpdatedDeviatoricCauchyStress[g] = rElementalVariables.UpdatedDeviatoricCauchyStress;
  }

  template <>
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedSolidElement<3>::CalcElasticPlasticCauchySplitted(
      ElementalVariables &rElementalVariables, double TimeStep, unsigned int g, const ProcessInfo &rCurrentProcessInfo,
      double &Density, double &DeviatoricCoeff, double &VolumetricCoeff)
  {

    mpConstitutiveLaw = this->GetProperties().GetValue(CONSTITUTIVE_LAW);
    auto constitutive_law_values =
        ConstitutiveLaw::Parameters(this->GetGeometry(), this->GetProperties(), rCurrentProcessInfo);

    Flags &constitutive_law_options = constitutive_law_values.GetOptions();
    constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

    rElementalVariables.CurrentTotalCauchyStress = this->mCurrentTotalCauchyStress[g];
    rElementalVariables.CurrentDeviatoricCauchyStress = this->mCurrentDeviatoricCauchyStress[g];

    const Vector &r_shape_functions = row((this->GetGeometry()).ShapeFunctionsValues(), g);
    constitutive_law_values.SetShapeFunctionsValues(r_shape_functions);
    constitutive_law_values.SetStrainVector(rElementalVariables.SpatialDefRate);
    constitutive_law_values.SetStressVector(rElementalVariables.CurrentDeviatoricCauchyStress);

    mpConstitutiveLaw->CalculateMaterialResponseCauchy(constitutive_law_values);

    Density = mpConstitutiveLaw->CalculateValue(constitutive_law_values, DENSITY, Density);

    double poisson_ratio = mpConstitutiveLaw->CalculateValue(constitutive_law_values, POISSON_RATIO, poisson_ratio);
    double young_modulus = mpConstitutiveLaw->CalculateValue(constitutive_law_values, YOUNG_MODULUS, young_modulus);
    const double time_step = rCurrentProcessInfo[DELTA_TIME];
    DeviatoricCoeff = time_step * young_modulus / (2.0 * (1 + poisson_ratio));
    VolumetricCoeff =
        time_step * poisson_ratio * young_modulus / ((1.0 + poisson_ratio) * (1.0 - 2.0 * poisson_ratio)) +
        2.0 / 3.0 * DeviatoricCoeff;

    const double current_first_lame = VolumetricCoeff - 2.0 / 3.0 * DeviatoricCoeff;

    this->mMaterialDeviatoricCoefficient = DeviatoricCoeff;
    this->mMaterialVolumetricCoefficient = VolumetricCoeff;
    this->mMaterialDensity = Density;

    rElementalVariables.UpdatedDeviatoricCauchyStress[0] = rElementalVariables.CurrentDeviatoricCauchyStress[0];
    rElementalVariables.UpdatedDeviatoricCauchyStress[1] = rElementalVariables.CurrentDeviatoricCauchyStress[1];
    rElementalVariables.UpdatedDeviatoricCauchyStress[2] = rElementalVariables.CurrentDeviatoricCauchyStress[2];
    rElementalVariables.UpdatedDeviatoricCauchyStress[3] = rElementalVariables.CurrentDeviatoricCauchyStress[3];
    rElementalVariables.UpdatedDeviatoricCauchyStress[4] = rElementalVariables.CurrentDeviatoricCauchyStress[4];
    rElementalVariables.UpdatedDeviatoricCauchyStress[5] = rElementalVariables.CurrentDeviatoricCauchyStress[5];

    rElementalVariables.UpdatedTotalCauchyStress[0] = +current_first_lame * rElementalVariables.VolumetricDefRate +
                                                      2.0 * DeviatoricCoeff * rElementalVariables.SpatialDefRate[0] +
                                                      rElementalVariables.CurrentTotalCauchyStress[0];
    rElementalVariables.UpdatedTotalCauchyStress[1] = current_first_lame * rElementalVariables.VolumetricDefRate +
                                                      2.0 * DeviatoricCoeff * rElementalVariables.SpatialDefRate[1] +
                                                      rElementalVariables.CurrentTotalCauchyStress[1];
    rElementalVariables.UpdatedTotalCauchyStress[2] = current_first_lame * rElementalVariables.VolumetricDefRate +
                                                      2.0 * DeviatoricCoeff * rElementalVariables.SpatialDefRate[2] +
                                                      rElementalVariables.CurrentTotalCauchyStress[2];
    rElementalVariables.UpdatedTotalCauchyStress[3] =
        2.0 * DeviatoricCoeff * rElementalVariables.SpatialDefRate[3] + rElementalVariables.CurrentTotalCauchyStress[3];
    rElementalVariables.UpdatedTotalCauchyStress[4] =
        2.0 * DeviatoricCoeff * rElementalVariables.SpatialDefRate[4] + rElementalVariables.CurrentTotalCauchyStress[4];
    rElementalVariables.UpdatedTotalCauchyStress[5] =
        2.0 * DeviatoricCoeff * rElementalVariables.SpatialDefRate[5] + rElementalVariables.CurrentTotalCauchyStress[5];

    this->mUpdatedTotalCauchyStress[g] = rElementalVariables.UpdatedTotalCauchyStress;
    this->mUpdatedDeviatoricCauchyStress[g] = rElementalVariables.UpdatedDeviatoricCauchyStress;
  }

  template class TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedSolidElement<2>;
  template class TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedSolidElement<3>;

} // namespace Kratos
