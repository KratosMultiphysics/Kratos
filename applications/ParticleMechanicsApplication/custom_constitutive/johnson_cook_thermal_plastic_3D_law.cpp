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


  void JohnsonCookThermalPlastic3DLaw::CalculateMaterialResponseKirchhoff(Parameters& rValues)
  {
      //a.-Check if the constitutive parameters are passed correctly to the law calculation
      CheckParameters(rValues);

      //b.- Get Values to compute the constitutive law:
      Flags& Options = rValues.GetOptions();
      KRATOS_ERROR_IF(Options.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN))
          << "The JohnsonCookThermalPlastic3DLaw cannot accept the deformation graddient F as a strain input."
          << " Please set the strain vector and set USE_ELEMENT_PROVIDED_STRAIN = False in the element.";

      const ProcessInfo& CurrentProcessInfo = rValues.GetProcessInfo();
      const Properties& MaterialProperties = rValues.GetMaterialProperties();

      const GeometryType& DomainGeometry = rValues.GetElementGeometry();
      const Vector& ShapeFunctions = rValues.GetShapeFunctionsValues();

      Vector& StrainVector = rValues.GetStrainVector();
      Vector& StressVector = rValues.GetStressVector();
      Matrix& ConstitutiveMatrix = rValues.GetConstitutiveMatrix();

      // Convert vectors to matrices for easier manipulation
      Matrix stress_old(3, 3);
      Matrix strain_increment(3, 3);
      MakeStrainStressMatrixFromVector((StrainVector - mStrainOld), strain_increment);
      MakeStrainStressMatrixFromVector(StressVector, stress_old);

      const double strain_increment_trace = strain_increment(0, 0) + strain_increment(1, 1) + strain_increment(2, 2);
      const Matrix strain_increment_hydrostatic = strain_increment_trace/3.0 * identity_matrix<double>(3);
      const Matrix strain_increment_deviatoric = strain_increment - strain_increment_hydrostatic;


      // Material moduli
      const double shear_modulus_G = MaterialProperties[YOUNG_MODULUS] / (2.0 + 2.0 * MaterialProperties[POISSON_RATIO]);
      const double bulk_modulus_K = MaterialProperties[YOUNG_MODULUS] / (3.0 - 6.0 * MaterialProperties[POISSON_RATIO]);

      

      // TODO decide if we skip the first timestep
      //if (CurrentProcessInfo.GetSolutionStepIndex() == 0)
      //{
      //
      //}

      const double stress_hydrostatic_old = (stress_old(0, 0) + stress_old(1, 1) + stress_old(2, 2)) / 3.0;
      Matrix stress_deviatoric_old = stress_old - stress_hydrostatic_old * identity_matrix<double>(3);
      const double j2_stress_old = std::sqrt(3.0 / 2.0 * CalculateMatrixDoubleContraction(stress_deviatoric_old));
      Matrix stress_deviatoric_trial = stress_deviatoric_old + 2.0 * shear_modulus_G * strain_increment_deviatoric;
      const double j2_stress_trial = std::sqrt(3.0 / 2.0 * CalculateMatrixDoubleContraction(stress_deviatoric_trial));

      ParticleYieldCriterion::Parameters yield_parameters;
      yield_parameters.SetStressNorm(j2_stress_trial);
      double yield_function;
      mpYieldCriterion->CalculateYieldCondition(yield_function, yield_parameters);

      if (yield_function > 0.0)
      {
          double current_yield = j2_stress_trial - yield_function;

          // Thermal properties
          const double eta = 0.9; // TODO check this
          const double density = MaterialProperties[DENSITY];
          const double specific_heat_Cp = MaterialProperties[SPECIFIC_HEAT];


          // Newton raphson setup
          double gamma = mGammaOld;
          double gamma_min = 0.0;
          double gamma_max = j2_stress_trial / std::sqrt(6.0) / shear_modulus_G;
          bool is_converged = false;
          SizeType iteration_limit = 50;
          IndexType iteration = 1;
          const double tolerance = 1e-9;
          double yield_function, yield_function_gradient, dYield_dGamma, delta_gamma;

          while (!is_converged)
          {
              double predicted_eps = mEquivalentPlasticStrainOld + std::sqrt(2.0 / 3.0 * gamma); // eps = equivalent plastic strain
              double predicted_eps_rate = std::sqrt(2.0 / 3.0 * gamma) / CurrentProcessInfo[DELTA_TIME];

              double predicted_temperature = mTemperatureOld +  eta / 2.0 / density / specific_heat_Cp * // TODO check this, [johnson cool umat pdfp169]
                  std::sqrt(2.0 / 3.0) * (current_yield + j2_stress_old)* gamma;

              //double predicted_temperature = mTemperatureOld +  eta / std::sqrt(6.0) / density / specific_heat_Cp * 
              //    (std::sqrt(2.0 / 3.0) * current_yield + j2_stress_old)* gamma;

              current_yield = mpHardeningLaw->CalculateHardening(current_yield, predicted_eps_rate, predicted_eps, predicted_temperature);
              yield_function = j2_stress_trial - std::sqrt(6.0) * shear_modulus_G * gamma - current_yield;
              if (yield_function < 0.0) gamma_max = gamma;
              else gamma_min = gamma;
              dYield_dGamma = mpFlowRule->CalculatePlasticStrainDerivative();
              dYield_dGamma += mpFlowRule->CalculatePlasticStrainRateDerivative()/ CurrentProcessInfo[DELTA_TIME];
              dYield_dGamma += mpFlowRule->CalculateThermalDerivative()*eta* current_yield/density/specific_heat_Cp;
              dYield_dGamma *= std::sqrt(2.0/3.0);
              yield_function_gradient = -1.0 * std::sqrt(6.0) * shear_modulus_G - dYield_dGamma;
              delta_gamma = -1.0 * yield_function / yield_function_gradient;
              gamma += delta_gamma;

              if (gamma_min > gamma || gamma > gamma_max)
              {
                  delta_gamma = 0.5 * (gamma_max - gamma_min);
                  gamma = gamma_min + delta_gamma;
              }

              if (std::abs(delta_gamma) < tolerance)
              {
                  is_converged = true;
              }

              iteration += 1;
              KRATOS_ERROR_IF(iteration > iteration_limit) << "Johnson Cook iteration limit exceeded";
          }

      }











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


      //3.-Compute Incremental DeformationGradientF_bar


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
          KRATOS_ERROR << "COMPUTE_CONSTITUTIVE_TENSOR not yet implemented in JohnsonCookThermalPlastic3DLaw";
      }



      if (Options.Is(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE))
      {
          mpFlowRule->UpdateInternalVariables(ReturnMappingVariables);

          mElasticLeftCauchyGreen = (IsochoricStressMatrix * (1.0 / ElasticVariables.LameMu));
          mElasticLeftCauchyGreen += (ElasticVariables.traceCG / 3.0) * ElasticVariables.Identity;

      }

      // Update old strain for next timestep
      for (size_t i = 0; i < mStrainOld.size(); ++i) mStrainOld[i] = StrainVector[i];
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

  void JohnsonCookThermalPlastic3DLaw::MakeStrainStressMatrixFromVector(const Vector& rInput, Matrix& rOutput)
  {
      if (rOutput.size1() != 3 || rOutput.size2() != 3)rOutput.resize(3, 3, false);

      // 3D stress arrangement
      rOutput(0, 0) = rInput[0];
      rOutput(1, 1) = rInput[1];
      rOutput(2, 2) = rInput[2];

      rOutput(0, 1) = 0.5*rInput[3]; //xy
      rOutput(1, 2) = 0.5*rInput[4]; //yz
      rOutput(0, 2) = 0.5*rInput[5]; //xz
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
