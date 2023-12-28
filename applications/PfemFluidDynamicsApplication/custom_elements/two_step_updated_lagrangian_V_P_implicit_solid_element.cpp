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
#include "custom_elements/two_step_updated_lagrangian_V_P_implicit_solid_element.h"
#include "includes/cfd_variables.h"

namespace Kratos
{

  /*
   * public TwoStepUpdatedLagrangianVPImplicitSolidElement<TDim> functions
   */

  template <unsigned int TDim>
  Element::Pointer TwoStepUpdatedLagrangianVPImplicitSolidElement<TDim>::Clone(IndexType NewId, NodesArrayType const &rThisNodes) const
  {
    // return Element::Pointer( BaseType::Clone(NewId,rThisNodes) );
    TwoStepUpdatedLagrangianVPImplicitSolidElement NewElement(NewId, this->GetGeometry().Create(rThisNodes), this->pGetProperties());

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

    return Element::Pointer(new TwoStepUpdatedLagrangianVPImplicitSolidElement(NewElement));
  }

  template <unsigned int TDim>
  void TwoStepUpdatedLagrangianVPImplicitSolidElement<TDim>::Initialize(const ProcessInfo &rCurrentProcessInfo)
  {
    KRATOS_TRY;

    const GeometryType &rGeom = this->GetGeometry();
    SizeType integration_points_number = rGeom.IntegrationPointsNumber(GeometryData::IntegrationMethod::GI_GAUSS_1);

    if (this->mCurrentTotalCauchyStress.size() != integration_points_number)
      this->mCurrentTotalCauchyStress.resize(integration_points_number);

    if (this->mCurrentDeviatoricCauchyStress.size() != integration_points_number)
      this->mCurrentDeviatoricCauchyStress.resize(integration_points_number);

    if (this->mUpdatedTotalCauchyStress.size() != integration_points_number)
      this->mUpdatedTotalCauchyStress.resize(integration_points_number);

    if (this->mUpdatedDeviatoricCauchyStress.size() != integration_points_number)
      this->mUpdatedDeviatoricCauchyStress.resize(integration_points_number);

    unsigned int voigtsize = 3;
    if constexpr (TDim == 3)
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

    // If we are restarting, the constitutive law will be already defined
    if (mpConstitutiveLaw == nullptr)
    {
      const Properties &r_properties = this->GetProperties();
      KRATOS_ERROR_IF_NOT(r_properties.Has(CONSTITUTIVE_LAW))
          << "In initialization of Element " << this->Info() << ": No CONSTITUTIVE_LAW defined for property "
          << r_properties.Id() << "." << std::endl;
      mpConstitutiveLaw = r_properties[CONSTITUTIVE_LAW]->Clone();
    }

    KRATOS_CATCH("");
  }

    template <unsigned int TDim>
    void TwoStepUpdatedLagrangianVPImplicitSolidElement<TDim>::InitializeSolutionStep(const ProcessInfo &rCurrentProcessInfo)
    {
        KRATOS_TRY;

        // TODO: Temporary solution until the mesher calls the Initialize() after the elements creation
        if (!mpConstitutiveLaw) {
            this->Initialize(rCurrentProcessInfo);
        }

        KRATOS_CATCH("");
    }

  template <unsigned int TDim>
  void TwoStepUpdatedLagrangianVPImplicitSolidElement<TDim>::InitializeNonLinearIteration(const ProcessInfo &rCurrentProcessInfo)
  {
    KRATOS_TRY;

    const GeometryType &rGeom = this->GetGeometry();
    SizeType integration_points_number = rGeom.IntegrationPointsNumber(GeometryData::IntegrationMethod::GI_GAUSS_1);
    // SizeType integration_points_number = rGeom.IntegrationPointsNumber(GeometryData::IntegrationMethod::GI_GAUSS_4);

    for (unsigned int PointNumber = 0; PointNumber < integration_points_number; PointNumber++)
    {
      this->UpdateCauchyStress(PointNumber, rCurrentProcessInfo);
    }

    KRATOS_CATCH("");
  }

  template <unsigned int TDim>
  int TwoStepUpdatedLagrangianVPImplicitSolidElement<TDim>::Check(const ProcessInfo &rCurrentProcessInfo) const
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

    // WARNING THIS MUST BE REMOVED ASAP
    const_cast<TwoStepUpdatedLagrangianVPImplicitSolidElement<TDim> *>(this)->mpConstitutiveLaw = const_cast<TwoStepUpdatedLagrangianVPImplicitSolidElement<TDim> *>(this)->GetProperties().GetValue(CONSTITUTIVE_LAW);
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
    return mpConstitutiveLaw->Check(r_properties, r_geometry, rCurrentProcessInfo);

    return ierr;

    KRATOS_CATCH("");
  }

  template <unsigned int TDim>
  void TwoStepUpdatedLagrangianVPImplicitSolidElement<TDim>::ComputeBulkMatrixForPressureVelLump(Matrix &BulkVelMatrix,
                                                                                                 const double Weight)
  {
    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    double coeff = 1.0 + TDim;
    // coeff=6.0;
    if (TDim == 2 && NumNodes == 6)
    {
      double Mij = Weight / 57.0;
      double consistent = 1.0;
      for (SizeType i = 0; i < NumNodes; ++i)
      {
        if (i < 3)
        {
          consistent = coeff;
        }
        else
        {
          consistent = 16.0;
        }

        BulkVelMatrix(i, i) += Mij * consistent;
      }
    }
    else
    {
      for (SizeType i = 0; i < NumNodes; ++i)
      {
        // LHS contribution
        double Mij = Weight / coeff;
        BulkVelMatrix(i, i) += Mij;
      }
      if (NumNodes > 4)
      {
        std::cout << "ComputeBulkMatrixForPressureVelLump 3D not yet implemented!" << std::endl;
      }
    }
  }

  template <unsigned int TDim>
  void TwoStepUpdatedLagrangianVPImplicitSolidElement<TDim>::ComputeBulkMatrixForPressureVel(Matrix &BulkVelMatrix,
                                                                                             const ShapeFunctionsType &rN,
                                                                                             const double Weight)
  {

    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    for (SizeType i = 0; i < NumNodes; ++i)
    {
      for (SizeType j = 0; j < NumNodes; ++j)
      {
        // LHS contribution
        double Mij = Weight * rN[i] * rN[j];
        BulkVelMatrix(i, j) += Mij;
      }
    }
  }

  template <unsigned int TDim>
  void TwoStepUpdatedLagrangianVPImplicitSolidElement<TDim>::InitializeElementalVariables(ElementalVariables &rElementalVariables)
  {
    unsigned int voigtsize = 3;
    if constexpr (TDim == 3)
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
    rElementalVariables.Fgrad.resize(TDim, TDim, false);
    rElementalVariables.InvFgrad.resize(TDim, TDim, false);
    rElementalVariables.FgradVel.resize(TDim, TDim, false);
    rElementalVariables.InvFgradVel.resize(TDim, TDim, false);
    rElementalVariables.SpatialVelocityGrad.resize(TDim, TDim, false);

    rElementalVariables.MeanPressure = 0;
    rElementalVariables.CurrentTotalCauchyStress.resize(voigtsize);
    rElementalVariables.UpdatedTotalCauchyStress.resize(voigtsize);
    rElementalVariables.CurrentDeviatoricCauchyStress.resize(voigtsize);
    rElementalVariables.UpdatedDeviatoricCauchyStress.resize(voigtsize);
  }

  template <>
  void TwoStepUpdatedLagrangianVPImplicitSolidElement<2>::UpdateStressTensor(ElementalVariables &rElementalVariables)
  {
    KRATOS_TRY;

    rElementalVariables.UpdatedDeviatoricCauchyStress[0] = rElementalVariables.CurrentDeviatoricCauchyStress[0];
    rElementalVariables.UpdatedDeviatoricCauchyStress[1] = rElementalVariables.CurrentDeviatoricCauchyStress[1];
    rElementalVariables.UpdatedDeviatoricCauchyStress[2] = rElementalVariables.CurrentDeviatoricCauchyStress[2];

    rElementalVariables.UpdatedTotalCauchyStress[0] =
        rElementalVariables.CurrentDeviatoricCauchyStress[0] + rElementalVariables.MeanPressure;
    rElementalVariables.UpdatedTotalCauchyStress[1] =
        rElementalVariables.CurrentDeviatoricCauchyStress[1] + rElementalVariables.MeanPressure;
    rElementalVariables.UpdatedTotalCauchyStress[2] = rElementalVariables.CurrentDeviatoricCauchyStress[2];

    this->SetValue(CAUCHY_STRESS_VECTOR, rElementalVariables.UpdatedTotalCauchyStress);

    KRATOS_CATCH("");
  }

  template <>
  void TwoStepUpdatedLagrangianVPImplicitSolidElement<3>::UpdateStressTensor(ElementalVariables &rElementalVariables)
  {
    KRATOS_TRY;

    rElementalVariables.UpdatedDeviatoricCauchyStress[0] = rElementalVariables.CurrentDeviatoricCauchyStress[0];
    rElementalVariables.UpdatedDeviatoricCauchyStress[1] = rElementalVariables.CurrentDeviatoricCauchyStress[1];
    rElementalVariables.UpdatedDeviatoricCauchyStress[2] = rElementalVariables.CurrentDeviatoricCauchyStress[2];
    rElementalVariables.UpdatedDeviatoricCauchyStress[3] = rElementalVariables.CurrentDeviatoricCauchyStress[3];
    rElementalVariables.UpdatedDeviatoricCauchyStress[4] = rElementalVariables.CurrentDeviatoricCauchyStress[4];
    rElementalVariables.UpdatedDeviatoricCauchyStress[5] = rElementalVariables.CurrentDeviatoricCauchyStress[5];

    rElementalVariables.UpdatedTotalCauchyStress[0] =
        rElementalVariables.CurrentDeviatoricCauchyStress[0] + rElementalVariables.MeanPressure;
    rElementalVariables.UpdatedTotalCauchyStress[1] =
        rElementalVariables.CurrentDeviatoricCauchyStress[1] + rElementalVariables.MeanPressure;
    rElementalVariables.UpdatedTotalCauchyStress[2] =
        rElementalVariables.CurrentDeviatoricCauchyStress[2] + rElementalVariables.MeanPressure;
    rElementalVariables.UpdatedTotalCauchyStress[3] = rElementalVariables.CurrentDeviatoricCauchyStress[3];
    rElementalVariables.UpdatedTotalCauchyStress[4] = rElementalVariables.CurrentDeviatoricCauchyStress[4];
    rElementalVariables.UpdatedTotalCauchyStress[5] = rElementalVariables.CurrentDeviatoricCauchyStress[5];

    this->SetValue(CAUCHY_STRESS_VECTOR, rElementalVariables.UpdatedTotalCauchyStress);

    KRATOS_CATCH("");
  }
  template <>
  void TwoStepUpdatedLagrangianVPImplicitSolidElement<2>::CalcElasticPlasticCauchySplitted(
        ElementalVariables &rElementalVariables,
        const unsigned int g,
        const Vector &rN,
        const ProcessInfo &rCurrentProcessInfo,
        double &Density,
        double &DeviatoricCoeff,
        double &VolumetricCoeff)
  {

    rElementalVariables.CurrentTotalCauchyStress = this->mCurrentTotalCauchyStress[g];
    rElementalVariables.CurrentDeviatoricCauchyStress = this->mCurrentDeviatoricCauchyStress[g];

    // mpConstitutiveLaw = this->GetProperties().GetValue(CONSTITUTIVE_LAW);
    auto constitutive_law_values =
        ConstitutiveLaw::Parameters(this->GetGeometry(), this->GetProperties(), rCurrentProcessInfo);

    Flags &constitutive_law_options = constitutive_law_values.GetOptions();
    constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

    constitutive_law_values.SetShapeFunctionsValues(rN);
    constitutive_law_values.SetStrainVector(rElementalVariables.SpatialDefRate);
    constitutive_law_values.SetStressVector(rElementalVariables.CurrentDeviatoricCauchyStress);

    mpConstitutiveLaw->CalculateMaterialResponseCauchy(constitutive_law_values);

    this->UpdateStressTensor(rElementalVariables);

    this->mUpdatedTotalCauchyStress[g] = rElementalVariables.UpdatedTotalCauchyStress;
    this->mUpdatedDeviatoricCauchyStress[g] = rElementalVariables.UpdatedDeviatoricCauchyStress;

    Density = mpConstitutiveLaw->CalculateValue(constitutive_law_values, DENSITY, Density);

    double poisson_ratio = mpConstitutiveLaw->CalculateValue(constitutive_law_values, POISSON_RATIO, poisson_ratio);
    double young_modulus = mpConstitutiveLaw->CalculateValue(constitutive_law_values, YOUNG_MODULUS, young_modulus);
    const double time_step = rCurrentProcessInfo[DELTA_TIME];
    DeviatoricCoeff = time_step * young_modulus / (2.0 * (1.0 + poisson_ratio));
    VolumetricCoeff =
        time_step * poisson_ratio * young_modulus / ((1.0 + poisson_ratio) * (1.0 - 2.0 * poisson_ratio)) +
        2.0 / 3.0 * DeviatoricCoeff;

    this->mMaterialDeviatoricCoefficient = DeviatoricCoeff;
    this->mMaterialVolumetricCoefficient = VolumetricCoeff;
    this->mMaterialDensity = Density;

    this->ComputeMechanicalDissipation(rElementalVariables);
  }

  template <>
  void TwoStepUpdatedLagrangianVPImplicitSolidElement<3>::CalcElasticPlasticCauchySplitted(
        ElementalVariables &rElementalVariables,
        const unsigned int g,
        const Vector &rN,
        const ProcessInfo &rCurrentProcessInfo,
        double &Density,
        double &DeviatoricCoeff,
        double &VolumetricCoeff)
  {

    rElementalVariables.CurrentTotalCauchyStress = this->mCurrentTotalCauchyStress[g];
    rElementalVariables.CurrentDeviatoricCauchyStress = this->mCurrentDeviatoricCauchyStress[g];

    // mpConstitutiveLaw = this->GetProperties().GetValue(CONSTITUTIVE_LAW);
    auto constitutive_law_values =
        ConstitutiveLaw::Parameters(this->GetGeometry(), this->GetProperties(), rCurrentProcessInfo);

    Flags &constitutive_law_options = constitutive_law_values.GetOptions();
    constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

    constitutive_law_values.SetShapeFunctionsValues(rN);
    constitutive_law_values.SetStrainVector(rElementalVariables.SpatialDefRate);
    constitutive_law_values.SetStressVector(rElementalVariables.CurrentDeviatoricCauchyStress);

    mpConstitutiveLaw->CalculateMaterialResponseCauchy(constitutive_law_values);

    this->UpdateStressTensor(rElementalVariables);

    this->mUpdatedTotalCauchyStress[g] = rElementalVariables.UpdatedTotalCauchyStress;
    this->mUpdatedDeviatoricCauchyStress[g] = rElementalVariables.UpdatedDeviatoricCauchyStress;

    Density = mpConstitutiveLaw->CalculateValue(constitutive_law_values, DENSITY, Density);

    double poisson_ratio = mpConstitutiveLaw->CalculateValue(constitutive_law_values, POISSON_RATIO, poisson_ratio);
    double young_modulus = mpConstitutiveLaw->CalculateValue(constitutive_law_values, YOUNG_MODULUS, young_modulus);
    const double time_step = rCurrentProcessInfo[DELTA_TIME];
    DeviatoricCoeff = time_step * young_modulus / (2.0 * (1.0 + poisson_ratio));
    VolumetricCoeff =
        time_step * poisson_ratio * young_modulus / ((1.0 + poisson_ratio) * (1.0 - 2.0 * poisson_ratio)) +
        2.0 / 3.0 * DeviatoricCoeff;

    this->mMaterialDeviatoricCoefficient = DeviatoricCoeff;
    this->mMaterialVolumetricCoefficient = VolumetricCoeff;
    this->mMaterialDensity = Density;

    this->ComputeMechanicalDissipation(rElementalVariables);
  }

  template <unsigned int TDim>
  void TwoStepUpdatedLagrangianVPImplicitSolidElement<TDim>::UpdateCauchyStress(unsigned int g,
                                                                                const ProcessInfo &rCurrentProcessInfo)
  {
    double theta = this->GetThetaContinuity();
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
      const auto& r_N = row(NContainer, g);
      this->CalcElasticPlasticCauchySplitted(rElementalVariables, g, r_N, rCurrentProcessInfo, Density, DeviatoricCoeff, VolumetricCoeff);
    }

    this->mCurrentTotalCauchyStress[g] = this->mUpdatedTotalCauchyStress[g];
    this->mCurrentDeviatoricCauchyStress[g] = this->mUpdatedDeviatoricCauchyStress[g];
  }

  template <unsigned int TDim>
  void TwoStepUpdatedLagrangianVPImplicitSolidElement<TDim>::CalculateLocalContinuityEqForPressure(MatrixType &rLeftHandSideMatrix,
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

    // MatrixType BulkVelMatrix = ZeroMatrix(NumNodes,NumNodes);

    // Shape functions and integration points
    ShapeFunctionDerivativesArrayType DN_DX;
    Matrix NContainer;
    VectorType GaussWeights;
    this->CalculateGeometryData(DN_DX, NContainer, GaussWeights);
    const unsigned int NumGauss = GaussWeights.size();

    // const double TimeStep=rCurrentProcessInfo[DELTA_TIME];
    double theta = this->GetThetaContinuity();

    ElementalVariables rElementalVariables;
    this->InitializeElementalVariables(rElementalVariables);

    // double Density  = mMaterialDensity;
    // double DeviatoricCoeff = mMaterialVolumetricCoefficient;
    double VolumetricCoeff = this->mMaterialVolumetricCoefficient;

    double totalVolume = 0;

    // Loop on integration points
    for (unsigned int g = 0; g < NumGauss; ++g)
    {
      const double GaussWeight = GaussWeights[g];
      totalVolume += GaussWeight;
      const ShapeFunctionsType &N = row(NContainer, g);
      const ShapeFunctionDerivativesType &rDN_DX = DN_DX[g];
      // bool computeElement=this->CalcStrainRate(rElementalVariables,rCurrentProcessInfo,rDN_DX,theta);
      bool computeElement = this->CalcCompleteStrainRate(rElementalVariables, rCurrentProcessInfo, rDN_DX, theta);
      if (computeElement == true)
      {
        // double BulkCoeff =GaussWeight/(VolumetricCoeff);
        // this->ComputeBulkMatrixForPressureVel(BulkVelMatrix,N,BulkCoeff);

        for (SizeType i = 0; i < NumNodes; ++i)
        {
          // RHS contribution
          // Velocity divergence
          double RHSi = N[i] * rElementalVariables.VolumetricDefRate;
          rRightHandSideVector[i] += GaussWeight * RHSi;
        }
      }
    }

    MatrixType BulkVelMatrixLump = ZeroMatrix(NumNodes, NumNodes);
    double lumpedBulkCoeff = totalVolume / (VolumetricCoeff);
    this->ComputeBulkMatrixForPressureVelLump(BulkVelMatrixLump, lumpedBulkCoeff);

    rLeftHandSideMatrix += BulkVelMatrixLump;
    // rLeftHandSideMatrix+=BulkVelMatrix;

    VectorType UpdatedPressure = ZeroVector(NumNodes);
    VectorType CurrentPressure = ZeroVector(NumNodes);
    ;

    this->GetPressureValues(UpdatedPressure, 0);
    this->GetPressureValues(CurrentPressure, 1);

    VectorType DeltaPressure = UpdatedPressure - CurrentPressure;

    rRightHandSideVector -= prod(BulkVelMatrixLump, DeltaPressure);
    // rRightHandSideVector -= prod(BulkVelMatrix,DeltaPressure);
  }

  template class TwoStepUpdatedLagrangianVPImplicitSolidElement<2>;
  template class TwoStepUpdatedLagrangianVPImplicitSolidElement<3>;

} // namespace Kratos
