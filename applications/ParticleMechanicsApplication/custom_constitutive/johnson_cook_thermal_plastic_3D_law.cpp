//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    JMCarbonell
//					 (adapted to Particle Mechanics by Peter Wilson)
//

// System includes

// External includes

// Project includes

#include "custom_constitutive/johnson_cook_thermal_plastic_3D_law.hpp"
#include "particle_mechanics_application_variables.h"

namespace Kratos
{

  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

    JohnsonCookThermalPlastic3DLaw::JohnsonCookThermalPlastic3DLaw()
  : HyperElastic3DLaw()
  {
    mpHardeningLaw   = ParticleHardeningLaw::Pointer( new JohnsonCookThermalHardeningLaw() );
    mpYieldCriterion = ParticleYieldCriterion::Pointer( new JohnsonCookThermalYieldCriterion(mpHardeningLaw) );
    mpFlowRule       = ParticleFlowRule::Pointer( new JohnsonCookPlasticFlowRule(mpYieldCriterion) );
  }


  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

    JohnsonCookThermalPlastic3DLaw::JohnsonCookThermalPlastic3DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
  {
    mpHardeningLaw    =  pHardeningLaw;
    mpYieldCriterion  = ParticleYieldCriterion::Pointer( new JohnsonCookThermalYieldCriterion(mpHardeningLaw) );
    mpFlowRule        =  pFlowRule;
  }

  //******************************COPY CONSTRUCTOR**************************************
  //************************************************************************************

    JohnsonCookThermalPlastic3DLaw::JohnsonCookThermalPlastic3DLaw(const JohnsonCookThermalPlastic3DLaw& rOther)
  : HyperElastic3DLaw(rOther)
  {

  }

  //********************************CLONE***********************************************
  //************************************************************************************

  ConstitutiveLaw::Pointer JohnsonCookThermalPlastic3DLaw::Clone() const
  {
    return Kratos::make_shared<JohnsonCookThermalPlastic3DLaw>(*this);
  }

  //*******************************DESTRUCTOR*******************************************
  //************************************************************************************

  JohnsonCookThermalPlastic3DLaw::~JohnsonCookThermalPlastic3DLaw()
  {
  }

  //******************************* COMPUTE DOMAIN TEMPERATURE  ************************
  //************************************************************************************


  double & JohnsonCookThermalPlastic3DLaw::CalculateDomainTemperature (const MaterialResponseVariables & rElasticVariables,
										      double & rTemperature)
  {

    //1.-Temperature from nodes
    const GeometryType& DomainGeometry = rElasticVariables.GetElementGeometry();
    const Vector& ShapeFunctionsValues = rElasticVariables.GetShapeFunctionsValues();
    const unsigned int number_of_nodes = DomainGeometry.size();

    rTemperature=0;

    for ( unsigned int j = 0; j < number_of_nodes; j++ )
      {
     	rTemperature += ShapeFunctionsValues[j] * DomainGeometry[j].GetSolutionStepValue(TEMPERATURE);
      }


    return rTemperature;
  }

  //*******************************OPERATIONS FROM BASE CLASS***************************
  //************************************************************************************

  //***********************HAS : DOUBLE - VECTOR - MATRIX*******************************
  //************************************************************************************

  bool JohnsonCookThermalPlastic3DLaw::Has( const Variable<double>& rThisVariable )
  {
    if(rThisVariable == DELTA_PLASTIC_DISSIPATION || rThisVariable == PLASTIC_DISSIPATION ) return true;
    else KRATOS_ERROR << "Variable " << rThisVariable << " not implemented in Johnson Cook 3D material law function Has double.";

    return false;
  }

  bool JohnsonCookThermalPlastic3DLaw::Has(const Variable<Vector>& rThisVariable)
  {
     KRATOS_ERROR << "Variable " << rThisVariable << " not implemented in Johnson Cook 3D material law function Has Vector.";

      return false;
  }

  bool JohnsonCookThermalPlastic3DLaw::Has(const Variable<Matrix>& rThisVariable)
  {
     KRATOS_ERROR << "Variable " << rThisVariable << " not implemented in Johnson Cook 3D material law function Has Matrix.";

      return false;
  }






  void JohnsonCookThermalPlastic3DLaw::InitializeMaterial(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const Vector& rShapeFunctionsValues)
  {
      HyperElastic3DLaw::InitializeMaterial(rMaterialProperties, rElementGeometry, rShapeFunctionsValues);

      mElasticLeftCauchyGreen = identity_matrix<double>(3);

      mpHardeningLaw->SetProperties(rMaterialProperties);

      mpFlowRule->InitializeMaterial(mpYieldCriterion, mpHardeningLaw, rMaterialProperties);
      //mpFlowRule->InitializeMaterial( rMaterialProperties );
  }

  void JohnsonCookThermalPlastic3DLaw::CalculateMaterialResponsePK2(Parameters& rValues)
  {
      // TODO this is hyperelastic atm

      this->CalculateMaterialResponseKirchhoff(rValues);

      //1.- Obtain parameters
      Flags& Options = rValues.GetOptions();

      Vector& StressVector = rValues.GetStressVector();
      Vector& StrainVector = rValues.GetStrainVector();

      const Matrix& DeformationGradientF = rValues.GetDeformationGradientF();
      const double& DeterminantF = rValues.GetDeterminantF();

      Matrix& ConstitutiveMatrix = rValues.GetConstitutiveMatrix();

      //2.-Green-Lagrange Strain:
      if (Options.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN))
      {
          TransformStrains(StrainVector, DeformationGradientF, StrainMeasure_Almansi, StrainMeasure_GreenLagrange);
      }

      //3.-Calculate Total PK2 stress
      if (Options.Is(ConstitutiveLaw::COMPUTE_STRESS))
      {
          TransformStresses(StressVector, DeformationGradientF, DeterminantF, StressMeasure_Kirchhoff, StressMeasure_PK2);
      }

      //4.-Calculate PK2 constitutive tensor
      if (Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR))
      {
          PullBackConstitutiveMatrix(ConstitutiveMatrix, DeformationGradientF);
      }
  }

  void JohnsonCookThermalPlastic3DLaw::CalculateMaterialResponseKirchhoff(Parameters& rValues)
  {
      // TODO this is hyperelastic atm

      //a.-Check if the constitutive parameters are passed correctly to the law calculation
      CheckParameters(rValues);

      //b.- Get Values to compute the constitutive law:
      Flags& Options = rValues.GetOptions();

      const ProcessInfo& CurrentProcessInfo = rValues.GetProcessInfo();
      const Properties& MaterialProperties = rValues.GetMaterialProperties();

      const Matrix& DeformationGradientF = rValues.GetDeformationGradientF();
      const double& DeterminantF = rValues.GetDeterminantF();

      const GeometryType& DomainGeometry = rValues.GetElementGeometry();
      const Vector& ShapeFunctions = rValues.GetShapeFunctionsValues();

      Vector& StrainVector = rValues.GetStrainVector();
      Vector& StressVector = rValues.GetStressVector();
      Matrix& ConstitutiveMatrix = rValues.GetConstitutiveMatrix();

      //-----------------------------//

      //0.- Initialize parameters
      MaterialResponseVariables ElasticVariables;
      ElasticVariables.Identity = identity_matrix<double>(3);

      ElasticVariables.SetElementGeometry(DomainGeometry);
      ElasticVariables.SetShapeFunctionsValues(ShapeFunctions);

      ParticleFlowRule::RadialReturnVariables ReturnMappingVariables;
      ReturnMappingVariables.initialize(); //it has to be called at the start

      // Initialize variables from the process information
      ReturnMappingVariables.DeltaTime = CurrentProcessInfo[DELTA_TIME];

      if (CurrentProcessInfo[IMPLEX] == 1)
          ReturnMappingVariables.Options.Set(ParticleFlowRule::IMPLEX_ACTIVE, true);
      else
          ReturnMappingVariables.Options.Set(ParticleFlowRule::IMPLEX_ACTIVE, false);

      // Initialize Splited Parts: Isochoric and Volumetric stresses and constitutive tensors
      double voigtsize = StressVector.size();
      VectorSplit SplitStressVector;
      MatrixSplit SplitConstitutiveMatrix;

      //1.- Lame constants
      const double& YoungModulus = MaterialProperties[YOUNG_MODULUS];
      const double& PoissonCoefficient = MaterialProperties[POISSON_RATIO];

      ElasticVariables.LameLambda = (YoungModulus * PoissonCoefficient) / ((1 + PoissonCoefficient) * (1 - 2 * PoissonCoefficient));
      ElasticVariables.LameMu = YoungModulus / (2 * (1 + PoissonCoefficient));

      //1.1- Thermal constants
      if (MaterialProperties.Has(THERMAL_EXPANSION_COEFFICIENT))
          ElasticVariables.ThermalExpansionCoefficient = MaterialProperties[THERMAL_EXPANSION_COEFFICIENT];
      else
          ElasticVariables.ThermalExpansionCoefficient = 0;

      if (MaterialProperties.Has(REFERENCE_TEMPERATURE))
          ElasticVariables.ReferenceTemperature = MaterialProperties[REFERENCE_TEMPERATURE];
      else
          ElasticVariables.ReferenceTemperature = 0;


      //2.-Determinant of the Total DeformationGradientF
      ElasticVariables.DeterminantF = DeterminantF;

      //3.-Compute Incremental DeformationGradientF_bar
      double detF = DeterminantF / mDeterminantF0;

      ElasticVariables.J_pow13 = pow(detF, 1.0 / 3.0);

      ElasticVariables.DeformationGradientF = DeformationGradientF;

      ElasticVariables.DeformationGradientF = this->Transform2DTo3D(ElasticVariables.DeformationGradientF);

      ElasticVariables.DeformationGradientF = prod(ElasticVariables.DeformationGradientF, this->mInverseDeformationGradientF0);

      ElasticVariables.DeformationGradientF /= ElasticVariables.J_pow13; //now ElasticVariables.DeformationGradientF is DeformationGradientFbar

      //4.-Left Cauchy-Green tensor b_bar to the new configuration
      ElasticVariables.CauchyGreenMatrix.resize(3, 3, false);
      noalias(ElasticVariables.CauchyGreenMatrix) = prod(mElasticLeftCauchyGreen, trans(ElasticVariables.DeformationGradientF));
      ElasticVariables.CauchyGreenMatrix = prod(ElasticVariables.DeformationGradientF, ElasticVariables.CauchyGreenMatrix);


      //5.-Calculate trace of Left Cauchy-Green tensor b_bar
      ElasticVariables.traceCG = 0;
      for (unsigned int i = 0; i < 3; i++)
      {
          ElasticVariables.traceCG += ElasticVariables.CauchyGreenMatrix(i, i);
      }

      ReturnMappingVariables.LameMu_bar = ElasticVariables.LameMu * (ElasticVariables.traceCG / 3.0);

      //4.-Almansi Strain:
      if (Options.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN))
      {
          // e= 0.5*(1-invbT*invb)
          this->CalculateAlmansiStrain(ElasticVariables.CauchyGreenMatrix, StrainVector);
          // correct b_bar to b
          double J_pow23 = pow(ElasticVariables.DeterminantF, 2.0 / 3.0);
          StrainVector /= (J_pow23 * J_pow23);
      }


      //5.-Calculate Total Kirchhoff stress
      SplitStressVector.Isochoric.resize(voigtsize, false);
      noalias(SplitStressVector.Isochoric) = ZeroVector(voigtsize);
      Matrix IsochoricStressMatrix(3, 3);
      noalias(IsochoricStressMatrix) = ZeroMatrix(3, 3);

      if (Options.Is(ConstitutiveLaw::COMPUTE_STRESS) || Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR))
          this->CalculatePlasticIsochoricStress(ElasticVariables, ReturnMappingVariables, StressMeasure_Kirchhoff, IsochoricStressMatrix, SplitStressVector.Isochoric);

      if (Options.Is(ConstitutiveLaw::COMPUTE_STRESS))
      {

          SplitStressVector.Volumetric.resize(voigtsize, false);
          noalias(SplitStressVector.Volumetric) = ZeroVector(voigtsize);

          ElasticVariables.CauchyGreenMatrix = ElasticVariables.Identity;

          this->CalculateVolumetricStress(ElasticVariables, SplitStressVector.Volumetric);

          //Kirchhoff Stress:
          StressVector = SplitStressVector.Isochoric + SplitStressVector.Volumetric;

          // if( ReturnMappingVariables.Options.Is(FlowRule::PLASTIC_REGION) ){
          //   std::cout<<" StressVector.Isochoric "<<SplitStressVector.Isochoric<<std::endl;
          //   std::cout<<" StressVector.Volumetric "<<SplitStressVector.Volumetric<<std::endl;
          // }

          if (Options.Is(ConstitutiveLaw::ISOCHORIC_TENSOR_ONLY))
          {
              StressVector = SplitStressVector.Isochoric;
          }
          else if (Options.Is(ConstitutiveLaw::VOLUMETRIC_TENSOR_ONLY))
          {
              StressVector = SplitStressVector.Volumetric;
          }

      }


      if (Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR))
      {

          if (ReturnMappingVariables.Options.IsNot(ParticleFlowRule::RETURN_MAPPING_COMPUTED))
          {
              KRATOS_ERROR << " ReturnMappingCall was not performed  ...error in the constitutive calculation..." << std::endl;
          }

          //initialize constitutive tensors
          ConstitutiveMatrix.clear();
          SplitConstitutiveMatrix.Isochoric = ConstitutiveMatrix;
          SplitConstitutiveMatrix.Volumetric = ConstitutiveMatrix;
          SplitConstitutiveMatrix.Plastic = ConstitutiveMatrix;

          ElasticVariables.CauchyGreenMatrix = ElasticVariables.Identity;

          this->CalculateIsochoricConstitutiveMatrix(ElasticVariables, ReturnMappingVariables.TrialIsoStressMatrix, SplitConstitutiveMatrix.Isochoric);

          this->CalculateVolumetricConstitutiveMatrix(ElasticVariables, SplitConstitutiveMatrix.Volumetric);

          if (ReturnMappingVariables.Options.Is(ParticleFlowRule::PLASTIC_REGION))
              this->CalculatePlasticConstitutiveMatrix(ElasticVariables, ReturnMappingVariables, SplitConstitutiveMatrix.Plastic);


          // std::cout<< " Isochoric Constitutive "<<SplitConstitutiveMatrix.Isochoric<<std::endl;
          // std::cout<< " Volumetric Constitutive "<<SplitConstitutiveMatrix.Volumetric<<std::endl;
          //if( ReturnMappingVariables.Options.Is(FlowRule::PLASTIC_REGION) )
            //std::cout<< " Plastic Constitutive   "<<SplitConstitutiveMatrix.Plastic<<std::endl;

          ConstitutiveMatrix = SplitConstitutiveMatrix.Isochoric + SplitConstitutiveMatrix.Volumetric + SplitConstitutiveMatrix.Plastic;

          if (Options.Is(ConstitutiveLaw::ISOCHORIC_TENSOR_ONLY))
          {
              ConstitutiveMatrix = SplitConstitutiveMatrix.Isochoric + SplitConstitutiveMatrix.Plastic;
          }
          else if (Options.Is(ConstitutiveLaw::VOLUMETRIC_TENSOR_ONLY))
          {
              ConstitutiveMatrix = SplitConstitutiveMatrix.Volumetric;
          }
      }



      if (Options.Is(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE))
      {
          mpFlowRule->UpdateInternalVariables(ReturnMappingVariables);

          mElasticLeftCauchyGreen = (IsochoricStressMatrix * (1.0 / ElasticVariables.LameMu));
          mElasticLeftCauchyGreen += (ElasticVariables.traceCG / 3.0) * ElasticVariables.Identity;

      }



      // std::cout<<" StrainVector "<<StrainVector<<std::endl;
      // std::cout<<" StressVector "<<StressVector<<std::endl;
      // std::cout<<" ConstitutiveMatrix "<<ConstitutiveMatrix<<std::endl;

  }

  int JohnsonCookThermalPlastic3DLaw::Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo)
  {
      if (YOUNG_MODULUS.Key() == 0 || rMaterialProperties[YOUNG_MODULUS] <= 0.00)
      {
          KRATOS_ERROR << "YOUNG_MODULUS has Key zero or invalid value " << std::endl;
      }

      const double& nu = rMaterialProperties[POISSON_RATIO];
      const bool check = bool((nu > 0.499 && nu < 0.501) || (nu < -0.999 && nu > -1.01));

      if (POISSON_RATIO.Key() == 0 || check == true)
      {
          KRATOS_ERROR << "POISSON_RATIO has Key zero invalid value " << std::endl;
      }


      if (DENSITY.Key() == 0 || rMaterialProperties[DENSITY] < 0.00)
      {
          KRATOS_ERROR << "DENSITY has Key zero or invalid value " << std::endl;
      }

      return 0;
  }

  double& JohnsonCookThermalPlastic3DLaw::PlasticConstitutiveComponent(double& rCabcd, const MaterialResponseVariables& rElasticVariables, const Matrix& rIsoStressMatrix, const ParticleFlowRule::PlasticFactors& rScalingFactors, const unsigned int& a, const unsigned int& b, const unsigned int& c, const unsigned int& d)
  {
      //Plastic part of the algorithmic moduli

      rCabcd = (1.0 / 3.0) * (rElasticVariables.CauchyGreenMatrix(a, b) * rElasticVariables.CauchyGreenMatrix(c, d));
      rCabcd -= (0.5 * (rElasticVariables.CauchyGreenMatrix(a, c) * rElasticVariables.CauchyGreenMatrix(b, d) + rElasticVariables.CauchyGreenMatrix(a, d) * rElasticVariables.CauchyGreenMatrix(b, c)));

      rCabcd *= rElasticVariables.traceCG * rElasticVariables.LameMu;

      rCabcd += (rElasticVariables.CauchyGreenMatrix(c, d) * rIsoStressMatrix(a, b) + rIsoStressMatrix(c, d) * rElasticVariables.CauchyGreenMatrix(a, b));

      rCabcd *= (-2.0 / 3.0) * ((-1) * rScalingFactors.Beta1);

      rCabcd -= rScalingFactors.Beta3 * 2.0 * (rElasticVariables.LameMu * (rElasticVariables.traceCG / 3.0)) * (rScalingFactors.Normal(a, b) * rScalingFactors.Normal(c, d));

      rCabcd -= rScalingFactors.Beta4 * 2.0 * (rElasticVariables.LameMu * (rElasticVariables.traceCG / 3.0)) * (rScalingFactors.Normal(a, b) * rScalingFactors.Dev_Normal(c, d));

      return rCabcd;
  }

  void JohnsonCookThermalPlastic3DLaw::CalculatePlasticIsochoricStress(MaterialResponseVariables& rElasticVariables, ParticleFlowRule::RadialReturnVariables& rReturnMappingVariables, StressMeasure rStressMeasure, Matrix& rIsoStressMatrix, Vector& rIsoStressVector)
  {
      //note.- rElasticVariables.traceCG is "traceCG_bar"

      if (rStressMeasure == StressMeasure_PK2)
      {

          //rElasticVariables.CauchyGreenMatrix is InverseRightCauchyGreen

          //2.-Isochoric part of the 2nd Piola Kirchhoff Stress Matrix
          rIsoStressMatrix = (rElasticVariables.Identity - (rElasticVariables.traceCG / 3.0) * rElasticVariables.CauchyGreenMatrix);
          rIsoStressMatrix *= rElasticVariables.LameMu;

          //std::cout<<" PK2 "<<std::endl;
      }

      if (rStressMeasure == StressMeasure_Kirchhoff)
      {

          //rElasticVariables.CauchyGreenMatrix is LeftCauchyGreen

          //2.-Isochoric part of the Kirchhoff Stress Matrix
          rIsoStressMatrix = (rElasticVariables.CauchyGreenMatrix - (rElasticVariables.traceCG / 3.0) * rElasticVariables.Identity);
          rIsoStressMatrix *= rElasticVariables.LameMu;

          //std::cout<<" Kirchooff "<<IsoStressMatrix<<std::endl;

      }


      //thermal effects:
      rReturnMappingVariables.Temperature = this->CalculateDomainTemperature(rElasticVariables, rReturnMappingVariables.Temperature);

      rReturnMappingVariables.TrialIsoStressMatrix = rIsoStressMatrix;

      //std::cout<<" TrialIsoStressMatrix "<<rReturnMappingVariables.TrialIsoStressMatrix<<std::endl;

      mpFlowRule->CalculateReturnMapping(rReturnMappingVariables, rIsoStressMatrix);

      rIsoStressVector = MathUtils<double>::StressTensorToVector(rIsoStressMatrix, rIsoStressVector.size());

      //std::cout<<" PLASTICITY "<<rElasticVariables.Plasticity<<" rIsoStressVector "<<rIsoStressVector<<std::endl;
  }

  void JohnsonCookThermalPlastic3DLaw::CalculateGreenLagrangeStrain(const Matrix& rRightCauchyGreen, Vector& rStrainVector)
  {
      // TODO this is currently hyperelastic 2D

      //E= 0.5*(FT*F-1)
      rStrainVector[0] = 0.5 * (rRightCauchyGreen(0, 0) - 1.00);
      rStrainVector[1] = 0.5 * (rRightCauchyGreen(1, 1) - 1.00);
      rStrainVector[2] = rRightCauchyGreen(0, 1);
  }

  void JohnsonCookThermalPlastic3DLaw::CalculateAlmansiStrain(const Matrix& rLeftCauchyGreen, Vector& rStrainVector)
  {
      // TODO this is hyperelastic currently 2D

      // e= 0.5*(1-invbT*invb)
      Matrix InverseLeftCauchyGreen(rLeftCauchyGreen.size1(), rLeftCauchyGreen.size2());
      double det_b = 0;
      MathUtils<double>::InvertMatrix(rLeftCauchyGreen, InverseLeftCauchyGreen, det_b);

      rStrainVector.clear();
      rStrainVector[0] = 0.5 * (1.0 - InverseLeftCauchyGreen(0, 0));
      rStrainVector[1] = 0.5 * (1.0 - InverseLeftCauchyGreen(1, 1));
      rStrainVector[2] = -InverseLeftCauchyGreen(0, 1);
  }

  void JohnsonCookThermalPlastic3DLaw::CalculateIsochoricConstitutiveMatrix(const MaterialResponseVariables& rElasticVariables, const Matrix& rIsoStressMatrix, Matrix& rConstitutiveMatrix)
  {
      // TODO this is 2D

      rConstitutiveMatrix.clear();

      for (unsigned int i = 0; i < 3; i++)
      {
          for (unsigned int j = 0; j < 3; j++)
          {
              rConstitutiveMatrix(i, j) = IsochoricConstitutiveComponent(rConstitutiveMatrix(i, j), rElasticVariables, rIsoStressMatrix,
                  this->msIndexVoigt2D3C[i][0], this->msIndexVoigt2D3C[i][1], this->msIndexVoigt2D3C[j][0], this->msIndexVoigt2D3C[j][1]);
          }

      }
  }

  void JohnsonCookThermalPlastic3DLaw::CalculateVolumetricConstitutiveMatrix(const MaterialResponseVariables& rElasticVariables, Matrix& rConstitutiveMatrix)
  {
      // TODO this is 2D

      rConstitutiveMatrix.clear();

      Vector Factors(3);
      noalias(Factors) = ZeroVector(3);
      Factors = this->CalculateVolumetricPressureFactors(rElasticVariables, Factors);

      for (unsigned int i = 0; i < 3; i++)
      {
          for (unsigned int j = 0; j < 3; j++)
          {
              rConstitutiveMatrix(i, j) = VolumetricConstitutiveComponent(rConstitutiveMatrix(i, j), rElasticVariables, Factors,
                  this->msIndexVoigt2D3C[i][0], this->msIndexVoigt2D3C[i][1], this->msIndexVoigt2D3C[j][0], this->msIndexVoigt2D3C[j][1]);
          }

      }
  }

  void JohnsonCookThermalPlastic3DLaw::CalculatePlasticConstitutiveMatrix(const MaterialResponseVariables& rElasticVariables, ParticleFlowRule::RadialReturnVariables& rReturnMappingVariables, Matrix& rConstitutiveMatrix)
  {
      // TODO this is 2D

      rConstitutiveMatrix.clear();

      Matrix IsoStressMatrix = rReturnMappingVariables.TrialIsoStressMatrix;

      //std::cout<< " TrialStressMatrix 2D "<<IsoStressMatrix<<std::endl;

      ParticleFlowRule::PlasticFactors ScalingFactors;
      mpFlowRule->CalculateScalingFactors(rReturnMappingVariables, ScalingFactors);

      for (unsigned int i = 0; i < 3; i++)
      {
          for (unsigned int j = 0; j < 3; j++)
          {
              rConstitutiveMatrix(i, j) = PlasticConstitutiveComponent(rConstitutiveMatrix(i, j), rElasticVariables, IsoStressMatrix, ScalingFactors,
                  this->msIndexVoigt2D3C[i][0], this->msIndexVoigt2D3C[i][1], this->msIndexVoigt2D3C[j][0], this->msIndexVoigt2D3C[j][1]);
          }

      }
  }

  bool JohnsonCookThermalPlastic3DLaw::CheckParameters(Parameters& rValues)
  {
      return rValues.CheckAllParameters();
  }

  //***********************GET VALUE: DOUBLE - VECTOR - MATRIX**************************
  //************************************************************************************

  void JohnsonCookThermalPlastic3DLaw::GetLawFeatures(Features& rFeatures)
  {
      //Set the type of law
      rFeatures.mOptions.Set(PLANE_STRAIN_LAW);
      rFeatures.mOptions.Set(FINITE_STRAINS);
      rFeatures.mOptions.Set(ISOTROPIC);

      //Set strain measure required by the consitutive law
      rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

      //Set the strain size
      rFeatures.mStrainSize = GetStrainSize();

      //Set the spacedimension
      rFeatures.mSpaceDimension = WorkingSpaceDimension();
  }

  double& JohnsonCookThermalPlastic3DLaw::CalculateValue(Parameters& rParameterValues, const Variable<double>& rThisVariable, double& rValue)
  {
      return (this->GetValue(rThisVariable, rValue));
  }

  double& JohnsonCookThermalPlastic3DLaw::GetValue( const Variable<double>& rThisVariable, double& rValue )
  {
      if (rThisVariable == DETERMINANT_F)
      {
          rValue = mDeterminantF0;
      }
      else if (rThisVariable == PLASTIC_STRAIN)
      {
          const ParticleFlowRule::InternalVariables& InternalVariables = mpFlowRule->GetInternalVariables();
          rValue = InternalVariables.EquivalentPlasticStrain;
      }
      else if (rThisVariable == DELTA_PLASTIC_STRAIN)
      {
          const ParticleFlowRule::InternalVariables& InternalVariables = mpFlowRule->GetInternalVariables();
          rValue = InternalVariables.DeltaPlasticStrain;
      }
      else if (rThisVariable == PLASTIC_DISSIPATION)
      {
          const ParticleFlowRule::ThermalVariables& ThermalVariables = mpFlowRule->GetThermalVariables();
          rValue = ThermalVariables.PlasticDissipation;
      }
      else if (rThisVariable == DELTA_PLASTIC_DISSIPATION)
      {
          const ParticleFlowRule::ThermalVariables& ThermalVariables = mpFlowRule->GetThermalVariables();
          rValue = ThermalVariables.DeltaPlasticDissipation;
      }
      else KRATOS_ERROR << "Variable " << rThisVariable << " not implemented in Johnson Cook 3D material law function GetValue double.";

    return( rValue );
  }

  Vector& JohnsonCookThermalPlastic3DLaw::GetValue(const Variable<Vector>& rThisVariable, Vector& rValue)
  {
      KRATOS_ERROR << "Variable " << rThisVariable << " not implemented in Johnson Cook 3D material law function GetValue Vector.";
      return(rValue);
  }

  Matrix& JohnsonCookThermalPlastic3DLaw::GetValue(const Variable<Matrix>& rThisVariable, Matrix& rValue)
  {
      KRATOS_ERROR << "Variable " << rThisVariable << " not implemented in Johnson Cook 3D material law function GetValue Matrix.";
      return(rValue);
  }

  void JohnsonCookThermalPlastic3DLaw::SetValue(const Variable<double>& rThisVariable, const double& rValue, const ProcessInfo& rCurrentProcessInfo)
  {
      if (rThisVariable == DETERMINANT_F)
      {
          mDeterminantF0 = rValue;
      }
      else if (rThisVariable == PLASTIC_STRAIN)
      {
          mpFlowRule->SetInternalVariables().EquivalentPlasticStrain = rValue;
      }
      else if (rThisVariable == DELTA_PLASTIC_STRAIN)
      {
          mpFlowRule->SetInternalVariables().DeltaPlasticStrain = rValue;
      }
      else KRATOS_ERROR << "Variable " << rThisVariable << " not implemented in Johnson Cook 3D material law function SetValue double.";
  }

  void JohnsonCookThermalPlastic3DLaw::SetValue(const Variable<Vector>& rThisVariable, const Vector& rValue, const ProcessInfo& rCurrentProcessInfo)
  {
      KRATOS_ERROR << "Variable " << rThisVariable << " not implemented in Johnson Cook 3D material law function SetValue Vector.";
  }

  void JohnsonCookThermalPlastic3DLaw::SetValue(const Variable<Matrix>& rThisVariable, const Matrix& rValue, const ProcessInfo& rCurrentProcessInfo)
  {
      KRATOS_ERROR << "Variable " << rThisVariable << " not implemented in Johnson Cook 3D material law function SetValue Matrix.";
  }


} // Namespace Kratos
