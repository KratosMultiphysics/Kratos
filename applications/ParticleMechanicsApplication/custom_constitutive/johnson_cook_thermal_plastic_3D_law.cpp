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
  }


  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

    JohnsonCookThermalPlastic3DLaw::JohnsonCookThermalPlastic3DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
  {
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

      mStrainOld = ZeroVector(6);
      mEquivalentPlasticStrainOld = 0.0;
      mPlasticStrainRateOld = 0.0;

      mTemperatureOld = 273.0 + 25.0; // TODO retrieve from element

      mEnergyInternal = 0.0;
      mEnergyDissipated = 0.0;
      mYieldStressOld = CalculateHardenedYieldStress(rMaterialProperties, mEquivalentPlasticStrainOld, mPlasticStrainRateOld, mTemperatureOld);

      mGammaOld = 1e-8;
  }


  void JohnsonCookThermalPlastic3DLaw::CalculateMaterialResponseKirchhoff(Kratos::ConstitutiveLaw::Parameters& rValues)
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
      //const GeometryType& DomainGeometry = rValues.GetElementGeometry();
      //const Vector& ShapeFunctions = rValues.GetShapeFunctionsValues();

      array_1d<double, 6> StrainVector = rValues.GetStrainVector();
      array_1d<double, 6> StressVector = rValues.GetStressVector();
      //Matrix& ConstitutiveMatrix = rValues.GetConstitutiveMatrix();

      // Convert vectors to matrices for easier manipulation
      Matrix stress_old(3, 3);
      Matrix strain_increment(3, 3);
      MakeStrainStressMatrixFromVector((StrainVector - mStrainOld), strain_increment);
      MakeStrainStressMatrixFromVector(StressVector, stress_old);

      // Material moduli
      const double shear_modulus_G = MaterialProperties[YOUNG_MODULUS] / (2.0 + 2.0 * MaterialProperties[POISSON_RATIO]);
      const double bulk_modulus_K = MaterialProperties[YOUNG_MODULUS] / (3.0 - 6.0 * MaterialProperties[POISSON_RATIO]);
      const double density = MaterialProperties[DENSITY];

      // TODO retrieve current temperature
      
      // Calculate strain increments
      const double strain_increment_trace = strain_increment(0, 0) + strain_increment(1, 1) + strain_increment(2, 2);
      const Matrix strain_increment_hydrostatic = strain_increment_trace/3.0 * IdentityMatrix(3);
      const Matrix strain_increment_deviatoric = strain_increment - strain_increment_hydrostatic;

      // Calculate current (old) j2 stress
      const double stress_hydrostatic_old = (stress_old(0, 0) + stress_old(1, 1) + stress_old(2, 2)) / 3.0;
      Matrix stress_deviatoric_old = stress_old - stress_hydrostatic_old * IdentityMatrix(3);
      

      // Calculate trial (predicted) j2 stress
      double stress_hydrostatic_new = stress_hydrostatic_old + bulk_modulus_K * strain_increment_trace; // TODO (1)
      Matrix stress_deviatoric_trial = stress_deviatoric_old + 2.0 * shear_modulus_G * strain_increment_deviatoric;
      const double j2_stress_trial = std::sqrt(3.0 / 2.0 * CalculateMatrixDoubleContraction(stress_deviatoric_trial));

      // Declare deviatoric stress matrix to be used later
      Matrix stress_deviatoric_converged(3, 3);

      // Assume current yield stress is the same as the old (elastic predictor)
      double yield_stress = mYieldStressOld;

      if (j2_stress_trial > yield_stress)
      {
          const double j2_stress_old = std::sqrt(3.0 / 2.0 * CalculateMatrixDoubleContraction(stress_deviatoric_old));

          // Thermal properties
          const double eta = 0.9; // TODO check this
          const double specific_heat_Cp = MaterialProperties[SPECIFIC_HEAT];

          // Newton raphson setup
          double gamma = mGammaOld;
          double gamma_min = 0.0;
          double gamma_max = j2_stress_trial / std::sqrt(6.0) / shear_modulus_G;
          bool is_converged = false;
          const SizeType iteration_limit = 100;
          IndexType iteration = 0;
          const double tolerance = 1e-8;
          double yield_function, yield_function_gradient, dYield_dGamma, 
              delta_gamma, predicted_eps, predicted_eps_rate, 
              predicted_temperature;

          // Initial prediction of quantities
          predicted_eps = mEquivalentPlasticStrainOld + std::sqrt(2.0 / 3.0) * gamma; // eps = equivalent plastic strain
          predicted_eps_rate = std::sqrt(2.0 / 3.0) * gamma / CurrentProcessInfo[DELTA_TIME];
          predicted_temperature = mTemperatureOld + eta / std::sqrt(6.0) / density / specific_heat_Cp * // TODO check this, [johnson cool umat pdfp169]
              (yield_stress + mYieldStressOld) * gamma;

          // Newton Raphson return mapping loop
          while (!is_converged)
          {
              // Calculate predicted yield stress
              yield_stress = CalculateHardenedYieldStress(MaterialProperties, predicted_eps, predicted_eps_rate, predicted_temperature);

              // Compute yield function and derivative
              yield_function = j2_stress_trial - std::sqrt(6.0) * shear_modulus_G * gamma - yield_stress;
              if (yield_function < 0.0) gamma_max = gamma;
              else gamma_min = gamma;
              dYield_dGamma = CalculatePlasticStrainDerivative(MaterialProperties, predicted_eps, predicted_eps_rate, predicted_temperature);
              dYield_dGamma += CalculateThermalDerivative(MaterialProperties, predicted_eps, predicted_eps_rate, predicted_temperature);
              dYield_dGamma += CalculatePlasticStrainRateDerivative(MaterialProperties, predicted_eps, predicted_eps_rate, predicted_temperature);
              dYield_dGamma *= std::sqrt(2.0/3.0);
              yield_function_gradient = -1.0 * std::sqrt(6.0) * shear_modulus_G - dYield_dGamma;

              // Update gamma
              delta_gamma = -1.0 * yield_function / yield_function_gradient;
              gamma += delta_gamma;
              
              // Bisect increment if out of search bounds
              if (gamma_min > gamma || gamma > gamma_max) gamma = gamma_min + 0.5 * (gamma_max - gamma_min);

              if (std::abs(delta_gamma) < tolerance) {
                  is_converged = true;
                  break;
              }

              // Update of quantities
              predicted_eps = mEquivalentPlasticStrainOld + std::sqrt(2.0 / 3.0) * gamma; // eps = equivalent plastic strain
              predicted_eps_rate = std::sqrt(2.0 / 3.0) * gamma / CurrentProcessInfo[DELTA_TIME];
              predicted_temperature = mTemperatureOld + eta / std::sqrt(6.0) / density / specific_heat_Cp * // TODO check this, [johnson cool umat pdfp169]
                  (yield_stress + mYieldStressOld) * gamma;

              iteration += 1;
              KRATOS_ERROR_IF(iteration == iteration_limit) << "Johnson Cook iteration limit exceeded";
          }

          // Correct trial stress
          Matrix flow_direction_normalized = stress_deviatoric_trial / (j2_stress_trial / std::sqrt(3.0 / 2.0));
          stress_deviatoric_converged = stress_deviatoric_trial - 2.0 * shear_modulus_G * gamma * flow_direction_normalized;

          mEnergyDissipated += gamma / std::sqrt(6.0) / density * (mYieldStressOld + yield_stress);
          mEquivalentPlasticStrainOld = predicted_eps;
          mPlasticStrainRateOld = predicted_eps_rate;
          mTemperatureOld = predicted_temperature;
          mYieldStressOld = yield_stress;
      }
      else
      {
          stress_deviatoric_converged = stress_deviatoric_trial;
      }

      
      Matrix stress_converged = stress_deviatoric_converged + stress_hydrostatic_new * IdentityMatrix(3);

      // Store stresses and strains
      MakeStrainStressVectorFromMatrix(stress_converged, StressVector);
      mStrainOld = StrainVector;

      // Udpdate internal energy
      for (size_t i = 0; i < stress_converged.size1(); ++i)
      {
          for (size_t j = 0; j < stress_converged.size2(); ++j)
          {
              mEnergyInternal += 0.5*(stress_converged(i, j) + stress_old(i, j)) * strain_increment(i, j) / density;
          }
      }


      if (Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR))
      {
          KRATOS_ERROR << "COMPUTE_CONSTITUTIVE_TENSOR not yet implemented in JohnsonCookThermalPlastic3DLaw";
      }
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

  bool JohnsonCookThermalPlastic3DLaw::CheckParameters(Parameters& rValues)
  {
      return rValues.CheckAllParameters();
  }

  void JohnsonCookThermalPlastic3DLaw::MakeStrainStressVectorFromMatrix(const Matrix& rInput, array_1d<double, 6>& rOutput)
  {
      if (rOutput.size() != 6) rOutput.resize(6, false);

      // 3D stress arrangement
      rOutput[0] = rInput(0, 0);
      rOutput[1] = rInput(1, 1);
      rOutput[2] = rInput(2, 2);

      rOutput[3] = 2.0* rInput(0, 1); //xy
      rOutput[4] = 2.0* rInput(1, 2); //yz
      rOutput[5] = 2.0* rInput(0, 2); //xz
  }

  void JohnsonCookThermalPlastic3DLaw::MakeStrainStressMatrixFromVector(const array_1d<double, 6>& rInput, Matrix& rOutput)
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

  double JohnsonCookThermalPlastic3DLaw::CalculateHardenedYieldStress(const Properties& MaterialProperties, 
      const double EquivalentPlasticStrain, const double PlasticStrainRate, const double Temperature)
  {
      //Constant Parameters of the -- Johnson and Cook --:
      const double A = MaterialProperties[JC_PARAMETER_A];
      const double B = MaterialProperties[JC_PARAMETER_B];
      const double n = MaterialProperties[JC_PARAMETER_n];

      // Hardening formula is: = (A + B* ep^n) * strain_rate_hardening_factor * thermal_hardening_factor
      double hardened_stress = A+B*std::pow(EquivalentPlasticStrain,n);
      hardened_stress *= CalculateStrainRateHardeningFactor(MaterialProperties, PlasticStrainRate);
      hardened_stress *= CalculateThermalHardeningFactor(MaterialProperties, Temperature);

      return hardened_stress;
  }

  double JohnsonCookThermalPlastic3DLaw::CalculateThermalHardeningFactor(const Properties& MaterialProperties, const double Temperature)
  {
      // Calculate thermal hardening factor
      const double m = MaterialProperties[JC_PARAMETER_m];
      const double ReferenceTemperature = MaterialProperties[REFERENCE_TEMPERATURE];
      const double MeldTemperature = MaterialProperties[MELD_TEMPERATURE];

      double thermal_hardening_factor;
      if (Temperature < ReferenceTemperature)
      {
          thermal_hardening_factor = 1.0;
      }
      else if (Temperature >= MeldTemperature)
      {
          thermal_hardening_factor = 0.0;
      }
      else
      {
          thermal_hardening_factor = 1.0 - std::pow((Temperature - ReferenceTemperature) / (MeldTemperature - ReferenceTemperature), m);
      }

      return thermal_hardening_factor;
  }

  double JohnsonCookThermalPlastic3DLaw::CalculateStrainRateHardeningFactor(const Properties& MaterialProperties, const double PlasticStrainRate)
  {
      // Calculate strain rate hardening factor
      const double ReferenceStrainRate = MaterialProperties[REFERENCE_STRAIN_RATE];
      const double C = MaterialProperties[JC_PARAMETER_C];

      double strain_rate_hardening_factor = 1.0;
      if (PlasticStrainRate > ReferenceStrainRate)
      {
          strain_rate_hardening_factor += C * std::log(PlasticStrainRate / ReferenceStrainRate);
      }

      return strain_rate_hardening_factor;
  }

  double JohnsonCookThermalPlastic3DLaw::CalculateThermalDerivative(const Properties& MaterialProperties, const double EquivalentPlasticStrain, const double PlasticStrainRate, const double Temperature)
  {
      const double A = MaterialProperties[JC_PARAMETER_A];
      const double B = MaterialProperties[JC_PARAMETER_B];
      const double C = MaterialProperties[JC_PARAMETER_C];

      const double n = MaterialProperties[JC_PARAMETER_n];
      const double m = MaterialProperties[JC_PARAMETER_m];

      const double ReferenceTemperature = MaterialProperties[REFERENCE_TEMPERATURE];
      const double MeldTemperature = MaterialProperties[MELD_TEMPERATURE];
      const double ReferenceStrainRate = MaterialProperties[REFERENCE_STRAIN_RATE];

      double thermal_derivative = 0.0;
      if (ReferenceTemperature <= Temperature && Temperature <= MeldTemperature)
      {
          double strain_rate_hardening_factor = CalculateStrainRateHardeningFactor(MaterialProperties,PlasticStrainRate);

          double thermal_hardening_factor = std::pow((Temperature - ReferenceTemperature) / (MeldTemperature - ReferenceTemperature), m);

          double temp = -1.0 * m * (A + B * std::pow(EquivalentPlasticStrain, n)) / (Temperature - ReferenceTemperature);

          thermal_derivative = temp * strain_rate_hardening_factor * thermal_hardening_factor;
      }
      return thermal_derivative;
  }

  double JohnsonCookThermalPlastic3DLaw::CalculatePlasticStrainRateDerivative(const Properties& MaterialProperties, const double EquivalentPlasticStrain, const double PlasticStrainRate, const double Temperature)
  {
      const double A = MaterialProperties[JC_PARAMETER_A];
      const double B = MaterialProperties[JC_PARAMETER_B];
      const double C = MaterialProperties[JC_PARAMETER_C];

      const double n = MaterialProperties[JC_PARAMETER_n];
      const double m = MaterialProperties[JC_PARAMETER_m];

      const double ReferenceTemperature = MaterialProperties[REFERENCE_TEMPERATURE];
      const double MeldTemperature = MaterialProperties[MELD_TEMPERATURE];
      const double ReferenceStrainRate = MaterialProperties[REFERENCE_STRAIN_RATE];

      double plastic_strain_rate_derivative = 0.0;

      if (PlasticStrainRate >= ReferenceStrainRate)
      {
          double thermal_hardening_factor = CalculateThermalHardeningFactor(MaterialProperties,Temperature);

          plastic_strain_rate_derivative = C / PlasticStrainRate * (A + B * std::pow(EquivalentPlasticStrain, n)) * thermal_hardening_factor;
      }

      return plastic_strain_rate_derivative;
  }

  double JohnsonCookThermalPlastic3DLaw::CalculatePlasticStrainDerivative(const Properties& MaterialProperties, const double EquivalentPlasticStrain, const double PlasticStrainRate, const double Temperature)
  {
      const double A = MaterialProperties[JC_PARAMETER_A];
      const double B = MaterialProperties[JC_PARAMETER_B];
      const double C = MaterialProperties[JC_PARAMETER_C];

      const double n = MaterialProperties[JC_PARAMETER_n];
      const double m = MaterialProperties[JC_PARAMETER_m];

      const double ReferenceTemperature = MaterialProperties[REFERENCE_TEMPERATURE];
      const double MeldTemperature = MaterialProperties[MELD_TEMPERATURE];
      const double ReferenceStrainRate = MaterialProperties[REFERENCE_STRAIN_RATE];

      double plastic_strain_derivative = 0.0;

      // Calculate thermal hardening factor
      double thermal_hardening_factor = CalculateThermalHardeningFactor(MaterialProperties,Temperature);

      // Calculate strain rate hardening factor
      double strain_rate_hardening_factor = 1.0;
      if (PlasticStrainRate > ReferenceStrainRate)
      {
          strain_rate_hardening_factor += C * std::log(PlasticStrainRate / ReferenceStrainRate);
      }

      plastic_strain_derivative = n * B * std::pow(EquivalentPlasticStrain, (n - 1.0)) * strain_rate_hardening_factor * thermal_hardening_factor;

      return plastic_strain_derivative;
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
      KRATOS_ERROR << "Variable " << rThisVariable << " not implemented in Johnson Cook 3D material law function GetValue double.";

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
