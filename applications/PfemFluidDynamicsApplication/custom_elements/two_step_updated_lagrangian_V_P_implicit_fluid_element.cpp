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
#include "custom_elements/two_step_updated_lagrangian_V_P_implicit_fluid_element.h"
#include "includes/cfd_variables.h"
#include <cmath>

namespace Kratos
{

  /*
   * public TwoStepUpdatedLagrangianVPImplicitFluidElement<TDim> functions
   */

  template <unsigned int TDim>
  Element::Pointer TwoStepUpdatedLagrangianVPImplicitFluidElement<TDim>::Clone(IndexType NewId, NodesArrayType const &rThisNodes) const
  {

    TwoStepUpdatedLagrangianVPImplicitFluidElement NewElement(NewId, this->GetGeometry().Create(rThisNodes), this->pGetProperties());

    NewElement.SetData(this->GetData());
    NewElement.SetFlags(this->GetFlags());

    return Element::Pointer(new TwoStepUpdatedLagrangianVPImplicitFluidElement(NewElement));
  }

  template <unsigned int TDim>
  void TwoStepUpdatedLagrangianVPImplicitFluidElement<TDim>::Initialize(const ProcessInfo &rCurrentProcessInfo)
  {
    KRATOS_TRY;

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
    void TwoStepUpdatedLagrangianVPImplicitFluidElement<TDim>::InitializeSolutionStep(const ProcessInfo &rCurrentProcessInfo)
    {
        KRATOS_TRY;

        // TODO: Temporary solution until the mesher calls the Initialize() after the elements creation
        if (!mpConstitutiveLaw) {
            this->Initialize(rCurrentProcessInfo);
        }

        KRATOS_CATCH("");
    }

  template <unsigned int TDim>
  int TwoStepUpdatedLagrangianVPImplicitFluidElement<TDim>::Check(const ProcessInfo &rCurrentProcessInfo) const
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
    const_cast<TwoStepUpdatedLagrangianVPImplicitFluidElement<TDim> *>(this)->mpConstitutiveLaw = const_cast<TwoStepUpdatedLagrangianVPImplicitFluidElement<TDim> *>(this)->GetProperties().GetValue(CONSTITUTIVE_LAW);
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

  template <>
  void TwoStepUpdatedLagrangianVPImplicitFluidElement<2>::ComputeMeanValueMaterialTangentMatrix(ElementalVariables &rElementalVariables, double &MeanValue, const ShapeFunctionDerivativesType &rDN_DX, const double secondLame, double &bulkModulus, const double Weight, double &MeanValueMass, const double TimeStep)
  {
    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    const double FourThirds = 4.0 / 3.0;
    const double nTwoThirds = -2.0 / 3.0;

    double theta = 0.5;
    double Count = 0;
    for (SizeType j = 0; j < NumNodes; ++j)
    {
      for (SizeType i = 0; i < NumNodes; ++i)
      {
        const double lagDNXi = rDN_DX(i, 0) * rElementalVariables.InvFgrad(0, 0) + rDN_DX(i, 1) * rElementalVariables.InvFgrad(1, 0);
        const double lagDNYi = rDN_DX(i, 0) * rElementalVariables.InvFgrad(0, 1) + rDN_DX(i, 1) * rElementalVariables.InvFgrad(1, 1);
        const double lagDNXj = rDN_DX(j, 0) * rElementalVariables.InvFgrad(0, 0) + rDN_DX(j, 1) * rElementalVariables.InvFgrad(1, 0);
        const double lagDNYj = rDN_DX(j, 0) * rElementalVariables.InvFgrad(0, 1) + rDN_DX(j, 1) * rElementalVariables.InvFgrad(1, 1);
        // lagDNXi=rDN_DX(i,0);
        // lagDNYi=rDN_DX(i,1);
        // lagDNXj=rDN_DX(j,0);
        // lagDNYj=rDN_DX(j,1);

        // First Row
        MeanValue += fabs(Weight * ((FourThirds * secondLame + bulkModulus) * lagDNXi * lagDNXj + lagDNYi * lagDNYj * secondLame) * theta);
        MeanValue += fabs(Weight * ((nTwoThirds * secondLame + bulkModulus) * lagDNXi * lagDNYj + lagDNYi * lagDNXj * secondLame) * theta);

        // Second Row
        MeanValue += fabs(Weight * ((nTwoThirds * secondLame + bulkModulus) * lagDNYi * lagDNXj + lagDNXi * lagDNYj * secondLame) * theta);
        MeanValue += fabs(Weight * ((FourThirds * secondLame + bulkModulus) * lagDNYi * lagDNYj + lagDNXi * lagDNXj * secondLame) * theta);

        Count += 4.0;
      }
    }

    MeanValue *= 1.0 / Count;

    if (MeanValueMass != 0 && MeanValue != 0)
    {
      bulkModulus *= MeanValueMass * 2.0 / TimeStep / MeanValue;
    }
    else
    {
      std::cout << " DANGEROUS ELEMENT!!!!!!!!!!!!!!!!!!!!!" << std::endl;
      std::cout << " MeanValueMass=" << MeanValueMass;
      std::cout << "\t MeanValueMaterial= " << MeanValue;
      std::cout << "\t VolumetricCoeff= " << bulkModulus << std::endl;
      bulkModulus *= TimeStep;
    }
  }

  template <>
  void TwoStepUpdatedLagrangianVPImplicitFluidElement<3>::ComputeMeanValueMaterialTangentMatrix(ElementalVariables &rElementalVariables, double &MeanValue, const ShapeFunctionDerivativesType &rDN_DX, const double secondLame, double &bulkModulus, const double Weight, double &MeanValueMass, const double TimeStep)
  {

    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    const double FourThirds = 4.0 / 3.0;
    const double nTwoThirds = -2.0 / 3.0;

    double theta = 0.5;
    double Count = 0;
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
        // lagDNXi=rDN_DX(i,0);
        // lagDNYi=rDN_DX(i,1);
        // lagDNZi=rDN_DX(i,2);
        // lagDNXj=rDN_DX(j,0);
        // lagDNYj=rDN_DX(j,1);
        // lagDNZj=rDN_DX(j,2);

        // First Row
        MeanValue += fabs(Weight * ((FourThirds * secondLame + bulkModulus) * lagDNXi * lagDNXj + (lagDNYi * lagDNYj + lagDNZi * lagDNZj) * secondLame) * theta);
        MeanValue += fabs(Weight * ((nTwoThirds * secondLame + bulkModulus) * lagDNXi * lagDNYj + lagDNYi * lagDNXj * secondLame) * theta);
        MeanValue += fabs(Weight * ((nTwoThirds * secondLame + bulkModulus) * lagDNXi * lagDNZj + lagDNZi * lagDNXj * secondLame) * theta);

        // Second Row
        MeanValue += fabs(Weight * ((nTwoThirds * secondLame + bulkModulus) * lagDNYi * lagDNXj + lagDNXi * lagDNYj * secondLame) * theta);
        MeanValue += fabs(Weight * ((FourThirds * secondLame + bulkModulus) * lagDNYi * lagDNYj + (lagDNXi * lagDNXj + lagDNZi * lagDNZj) * secondLame) * theta);
        MeanValue += fabs(Weight * ((nTwoThirds * secondLame + bulkModulus) * lagDNYi * lagDNZj + lagDNZi * lagDNYj * secondLame) * theta);

        // Third Row
        MeanValue += fabs(Weight * ((nTwoThirds * secondLame + bulkModulus) * lagDNZi * lagDNXj + lagDNXi * lagDNZj * secondLame) * theta);
        MeanValue += fabs(Weight * ((nTwoThirds * secondLame + bulkModulus) * lagDNZi * lagDNYj + lagDNYi * lagDNZj * secondLame) * theta);
        MeanValue += fabs(Weight * ((FourThirds * secondLame + bulkModulus) * lagDNZi * lagDNZj + (lagDNXi * lagDNXj + lagDNYi * lagDNYj) * secondLame) * theta);
        Count += 9.0;
      }
    }
    MeanValue *= 1.0 / Count;

    if (MeanValueMass != 0 && MeanValue != 0)
    {
      bulkModulus *= MeanValueMass * 2.0 / TimeStep / MeanValue;
    }
    else
    {
      std::cout << " DANGEROUS ELEMENT!!!!!!!!!!!!!!!!!!!!!" << std::endl;
      std::cout << " MeanValueMass=" << MeanValueMass;
      std::cout << "\t MeanValueMaterial= " << MeanValue;
      std::cout << "\t VolumetricCoeff= " << bulkModulus << std::endl;
      bulkModulus *= TimeStep;
    }
  }

  template <>
  void TwoStepUpdatedLagrangianVPImplicitFluidElement<2>::ComputeBulkReductionCoefficient(MatrixType MassMatrix,
                                                                                          MatrixType StiffnessMatrix,
                                                                                          double &meanValueStiff,
                                                                                          double &bulkCoefficient,
                                                                                          double timeStep)
  {
    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    IndexType FirstRow = 0;
    IndexType FirstCol = 0;
    double meanValueMass = 0;
    double countStiff = 0;
    double countMass = 0;
    for (SizeType j = 0; j < NumNodes; ++j)
    {
      for (SizeType i = 0; i < NumNodes; ++i)
      {
        meanValueStiff += fabs(StiffnessMatrix(FirstRow, FirstCol));
        meanValueStiff += fabs(StiffnessMatrix(FirstRow, FirstCol + 1));
        meanValueStiff += fabs(StiffnessMatrix(FirstRow + 1, FirstCol));
        meanValueStiff += fabs(StiffnessMatrix(FirstRow + 1, FirstCol + 1));
        // Update Counter
        countStiff += 4.0;

        meanValueMass += fabs(MassMatrix(FirstRow, FirstCol));
        meanValueMass += fabs(MassMatrix(FirstRow + 1, FirstCol + 1));
        // Update Counter
        countMass += 2.0;

        FirstRow += 2;
      }
      FirstRow = 0;
      FirstCol += 2;
    }

    meanValueStiff *= 1.0 / countStiff;
    meanValueMass *= 1.0 / countMass;

    if (meanValueMass != 0 && meanValueStiff != 0)
    {
      bulkCoefficient = meanValueMass * 4 / (timeStep * meanValueStiff);
    }
    else
    {
      std::cout << " DANGEROUS ELEMENT!!!!!!!!!!!!!!!!!!!!!" << std::endl;
      std::cout << " coordinates " << this->GetGeometry()[0].X() << " " << this->GetGeometry()[0].Y() << std::endl;
      std::cout << " MeanValueMass=" << meanValueMass;
      std::cout << "\t MeanValueMaterial= " << meanValueStiff;
      bulkCoefficient = timeStep;
    }
  }

  template <>
  void TwoStepUpdatedLagrangianVPImplicitFluidElement<3>::ComputeBulkReductionCoefficient(MatrixType MassMatrix,
                                                                                          MatrixType StiffnessMatrix,
                                                                                          double &meanValueStiff,
                                                                                          double &bulkCoefficient,
                                                                                          double timeStep)
  {
    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    IndexType FirstRow = 0;
    IndexType FirstCol = 0;
    double meanValueMass = 0;
    double countStiff = 0;
    double countMass = 0;
    for (SizeType j = 0; j < NumNodes; ++j)
    {
      for (SizeType i = 0; i < NumNodes; ++i)
      {
        meanValueStiff += fabs(StiffnessMatrix(FirstRow, FirstCol));
        meanValueStiff += fabs(StiffnessMatrix(FirstRow, FirstCol + 1));
        meanValueStiff += fabs(StiffnessMatrix(FirstRow, FirstCol + 2));
        meanValueStiff += fabs(StiffnessMatrix(FirstRow + 1, FirstCol));
        meanValueStiff += fabs(StiffnessMatrix(FirstRow + 1, FirstCol + 1));
        meanValueStiff += fabs(StiffnessMatrix(FirstRow + 1, FirstCol + 2));
        meanValueStiff += fabs(StiffnessMatrix(FirstRow + 2, FirstCol));
        meanValueStiff += fabs(StiffnessMatrix(FirstRow + 2, FirstCol + 1));
        meanValueStiff += fabs(StiffnessMatrix(FirstRow + 2, FirstCol + 2));
        countStiff += 9.0;

        meanValueMass += fabs(MassMatrix(FirstRow, FirstCol));
        meanValueMass += fabs(MassMatrix(FirstRow + 1, FirstCol + 1));
        meanValueMass += fabs(MassMatrix(FirstRow + 2, FirstCol + 2));
        countMass += 3.0;

        // Update Counter
        FirstRow += 3;
      }
      FirstRow = 0;
      FirstCol += 3;
    }

    meanValueStiff *= 1.0 / countStiff;
    meanValueMass *= 1.0 / countMass;

    if (meanValueMass != 0 && meanValueStiff != 0)
    {
      bulkCoefficient = meanValueMass * 2.0 / timeStep / meanValueStiff;
    }
    else
    {
      std::cout << " DANGEROUS ELEMENT!!!!!!!!!!!!!!!!!!!!!" << std::endl;
      std::cout << " MeanValueMass=" << meanValueMass;
      std::cout << "\t MeanValueMaterial= " << meanValueStiff;
      bulkCoefficient = timeStep;
    }
  }

  template <unsigned int TDim>
  void TwoStepUpdatedLagrangianVPImplicitFluidElement<TDim>::ComputeBulkMatrixLump(Matrix &BulkMatrix,
                                                                                   const double Weight)
  {
    const SizeType NumNodes = this->GetGeometry().PointsNumber();

    double coeff = 1.0 + TDim;
    if ((NumNodes == 3 && TDim == 2) || (NumNodes == 4 && TDim == 3))
    {
      for (SizeType i = 0; i < NumNodes; ++i)
      {
        // LHS contribution
        double Mij = Weight / coeff;
        BulkMatrix(i, i) += Mij;
      }
    }
    else
    {
      KRATOS_ERROR << "... ComputeBulkMatrixLump TO IMPLEMENT" << std::endl;
    }
  }

  template <unsigned int TDim>
  void TwoStepUpdatedLagrangianVPImplicitFluidElement<TDim>::InitializeElementalVariables(ElementalVariables &rElementalVariables)
  {
    KRATOS_TRY;

    unsigned int voigtsize = 3;
    if constexpr (TDim == 3)
    {
      voigtsize = 6;
    }
    rElementalVariables.voigtsize = voigtsize;

    rElementalVariables.ConstitutiveMatrix = ZeroMatrix(voigtsize, voigtsize);

    rElementalVariables.DetFgrad = 1.0;

    rElementalVariables.DetFgradVel = 1.0;

    rElementalVariables.DeviatoricInvariant = 1.0;

    rElementalVariables.EquivalentStrainRate = 1.0;

    rElementalVariables.VolumetricDefRate = 1.0;

    rElementalVariables.SpatialDefRate = ZeroVector(voigtsize);

    rElementalVariables.MDGreenLagrangeMaterial.resize(voigtsize, false);

    noalias(rElementalVariables.MDGreenLagrangeMaterial) = ZeroVector(voigtsize);

    rElementalVariables.Fgrad = ZeroMatrix(TDim, TDim);

    rElementalVariables.InvFgrad = ZeroMatrix(TDim, TDim);

    rElementalVariables.FgradVel = ZeroMatrix(TDim, TDim);

    rElementalVariables.InvFgradVel = ZeroMatrix(TDim, TDim);

    rElementalVariables.SpatialVelocityGrad = ZeroMatrix(TDim, TDim);

    rElementalVariables.MeanPressure = 0;

    rElementalVariables.CurrentTotalCauchyStress = ZeroVector(voigtsize);

    rElementalVariables.UpdatedTotalCauchyStress = ZeroVector(voigtsize);

    rElementalVariables.CurrentDeviatoricCauchyStress = ZeroVector(voigtsize);

    rElementalVariables.UpdatedDeviatoricCauchyStress = ZeroVector(voigtsize);

    KRATOS_CATCH("");
  }

  template <>
  void TwoStepUpdatedLagrangianVPImplicitFluidElement<2>::UpdateStressTensor(ElementalVariables &rElementalVariables)
  {
    KRATOS_TRY;

    rElementalVariables.UpdatedTotalCauchyStress[0] =
        rElementalVariables.UpdatedDeviatoricCauchyStress[0] + rElementalVariables.MeanPressure;
    rElementalVariables.UpdatedTotalCauchyStress[1] =
        rElementalVariables.UpdatedDeviatoricCauchyStress[1] + rElementalVariables.MeanPressure;
    rElementalVariables.UpdatedTotalCauchyStress[2] = rElementalVariables.UpdatedDeviatoricCauchyStress[2];

    this->SetValue(CAUCHY_STRESS_VECTOR, rElementalVariables.UpdatedTotalCauchyStress);
    KRATOS_CATCH("");
  }

  template <>
  void TwoStepUpdatedLagrangianVPImplicitFluidElement<3>::UpdateStressTensor(ElementalVariables &rElementalVariables)
  {
    KRATOS_TRY;

    rElementalVariables.UpdatedTotalCauchyStress[0] =
        rElementalVariables.UpdatedDeviatoricCauchyStress[0] + rElementalVariables.MeanPressure;
    rElementalVariables.UpdatedTotalCauchyStress[1] =
        rElementalVariables.UpdatedDeviatoricCauchyStress[1] + rElementalVariables.MeanPressure;
    rElementalVariables.UpdatedTotalCauchyStress[2] =
        rElementalVariables.UpdatedDeviatoricCauchyStress[2] + rElementalVariables.MeanPressure;
    rElementalVariables.UpdatedTotalCauchyStress[3] = rElementalVariables.UpdatedDeviatoricCauchyStress[3];
    rElementalVariables.UpdatedTotalCauchyStress[4] = rElementalVariables.UpdatedDeviatoricCauchyStress[4];
    rElementalVariables.UpdatedTotalCauchyStress[5] = rElementalVariables.UpdatedDeviatoricCauchyStress[5];

    this->SetValue(CAUCHY_STRESS_VECTOR, rElementalVariables.UpdatedTotalCauchyStress);

    KRATOS_CATCH("");
  }

  template <>
  void TwoStepUpdatedLagrangianVPImplicitFluidElement<2>::SetYieldedElements(ElementalVariables &rElementalVariables)
  {
    KRATOS_TRY;
    if (this->GetProperties().Has(YIELD_SHEAR) && this->Has(YIELDED))
    {
      const double tolerance = 1e-10;
      if (this->GetProperties()[YIELD_SHEAR] > tolerance)
      {
        const double TauNorm = sqrt(0.5 * rElementalVariables.UpdatedDeviatoricCauchyStress[0] * rElementalVariables.UpdatedDeviatoricCauchyStress[0] +
                              0.5 * rElementalVariables.UpdatedDeviatoricCauchyStress[1] * rElementalVariables.UpdatedDeviatoricCauchyStress[1] +
                              rElementalVariables.UpdatedDeviatoricCauchyStress[2] * rElementalVariables.UpdatedDeviatoricCauchyStress[2]);

        if (TauNorm > this->GetProperties()[YIELD_SHEAR])
        {
          this->SetValue(YIELDED, true);
        }
        else
        {
          this->SetValue(YIELDED, false);
        }
      }
      else
      {
        if (this->Has(YIELDED))
          this->SetValue(YIELDED, false);
      }
    }
    KRATOS_CATCH("");
  }

  template <>
  void TwoStepUpdatedLagrangianVPImplicitFluidElement<3>::SetYieldedElements(ElementalVariables &rElementalVariables)
  {
    KRATOS_TRY;

    if (this->GetProperties().Has(YIELD_SHEAR) && this->Has(YIELDED))
    {
      const double tolerance = 1e-10;
      if (this->GetProperties()[YIELD_SHEAR] > tolerance)
      {

        const double TauNorm = sqrt(2.0 * rElementalVariables.UpdatedDeviatoricCauchyStress[0] * rElementalVariables.UpdatedDeviatoricCauchyStress[0] +
                              2.0 * rElementalVariables.UpdatedDeviatoricCauchyStress[1] * rElementalVariables.UpdatedDeviatoricCauchyStress[1] +
                              2.0 * rElementalVariables.UpdatedDeviatoricCauchyStress[2] * rElementalVariables.UpdatedDeviatoricCauchyStress[2] +
                              4.0 * rElementalVariables.UpdatedDeviatoricCauchyStress[3] * rElementalVariables.UpdatedDeviatoricCauchyStress[3] +
                              4.0 * rElementalVariables.UpdatedDeviatoricCauchyStress[4] * rElementalVariables.UpdatedDeviatoricCauchyStress[4] +
                              4.0 * rElementalVariables.UpdatedDeviatoricCauchyStress[5] * rElementalVariables.UpdatedDeviatoricCauchyStress[5]);

        if (TauNorm > this->GetProperties()[YIELD_SHEAR])
        {
          this->SetValue(YIELDED, true);
        }
        else
        {
          this->SetValue(YIELDED, false);
        }
      }
      else
      {
        if (this->Has(YIELDED))
          this->SetValue(YIELDED, false);
      }
    }

    KRATOS_CATCH("");
  }

  template <>
  void TwoStepUpdatedLagrangianVPImplicitFluidElement<2>::CalcElasticPlasticCauchySplitted(
      ElementalVariables &rElementalVariables,
      const unsigned int g,
      const Vector& rN,
      const ProcessInfo &rCurrentProcessInfo,
      double &Density,
      double &DeviatoricCoeff,
      double &VolumetricCoeff)
  {

    //mpConstitutiveLaw = this->GetProperties().GetValue(CONSTITUTIVE_LAW);
    auto constitutive_law_values =
        ConstitutiveLaw::Parameters(this->GetGeometry(), this->GetProperties(), rCurrentProcessInfo);

    Flags &constitutive_law_options = constitutive_law_values.GetOptions();
    constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    constitutive_law_values.SetShapeFunctionsValues(rN);
    constitutive_law_values.SetStrainVector(rElementalVariables.SpatialDefRate);
    constitutive_law_values.SetStressVector(rElementalVariables.UpdatedDeviatoricCauchyStress);
    constitutive_law_values.SetConstitutiveMatrix(rElementalVariables.ConstitutiveMatrix);

    mpConstitutiveLaw->CalculateMaterialResponseCauchy(constitutive_law_values);

    this->UpdateStressTensor(rElementalVariables);

    this->SetYieldedElements(rElementalVariables); // only for non-newtonian laws

    const double time_step = rCurrentProcessInfo[DELTA_TIME];
    const double bulk_modulus = this->GetProperties()[BULK_MODULUS];
    const int voigt_size = (this->GetGeometry().WorkingSpaceDimension() - 1) * 3;

    DeviatoricCoeff = rElementalVariables.ConstitutiveMatrix(voigt_size - 1, voigt_size - 1);
    VolumetricCoeff = bulk_modulus * time_step;
    Density = mpConstitutiveLaw->CalculateValue(constitutive_law_values, DENSITY, Density);

    this->mMaterialDeviatoricCoefficient = DeviatoricCoeff;
    this->mMaterialVolumetricCoefficient = VolumetricCoeff;
    this->mMaterialDensity = Density;

    this->ComputeMechanicalDissipation(rElementalVariables);
  }

  template <>
  void TwoStepUpdatedLagrangianVPImplicitFluidElement<3>::CalcElasticPlasticCauchySplitted(
      ElementalVariables &rElementalVariables,
      const unsigned int g,
      const Vector& rN,
      const ProcessInfo &rCurrentProcessInfo,
      double &Density,
      double &DeviatoricCoeff,
      double &VolumetricCoeff)
  {

    // mpConstitutiveLaw = this->GetProperties().GetValue(CONSTITUTIVE_LAW);
    auto constitutive_law_values =
        ConstitutiveLaw::Parameters(this->GetGeometry(), this->GetProperties(), rCurrentProcessInfo);

    Flags &constitutive_law_options = constitutive_law_values.GetOptions();
    constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    constitutive_law_values.SetShapeFunctionsValues(rN);
    constitutive_law_values.SetStrainVector(rElementalVariables.SpatialDefRate);
    constitutive_law_values.SetStressVector(rElementalVariables.UpdatedDeviatoricCauchyStress);
    constitutive_law_values.SetConstitutiveMatrix(rElementalVariables.ConstitutiveMatrix);

    mpConstitutiveLaw->CalculateMaterialResponseCauchy(constitutive_law_values);

    this->UpdateStressTensor(rElementalVariables);

    this->SetYieldedElements(rElementalVariables); // only for non-newtonian laws

    const double time_step = rCurrentProcessInfo[DELTA_TIME];
    const double bulk_modulus = this->GetProperties()[BULK_MODULUS];
    const int voigt_size = (this->GetGeometry().WorkingSpaceDimension() - 1) * 3;

    DeviatoricCoeff = rElementalVariables.ConstitutiveMatrix(voigt_size - 1, voigt_size - 1);
    VolumetricCoeff = bulk_modulus * time_step;
    Density = mpConstitutiveLaw->CalculateValue(constitutive_law_values, DENSITY, Density);

    this->mMaterialDeviatoricCoefficient = DeviatoricCoeff;
    this->mMaterialVolumetricCoefficient = VolumetricCoeff;
    this->mMaterialDensity = Density;

    this->ComputeMechanicalDissipation(rElementalVariables);
  }

  template class TwoStepUpdatedLagrangianVPImplicitFluidElement<2>;
  template class TwoStepUpdatedLagrangianVPImplicitFluidElement<3>;

} // namespace Kratos
