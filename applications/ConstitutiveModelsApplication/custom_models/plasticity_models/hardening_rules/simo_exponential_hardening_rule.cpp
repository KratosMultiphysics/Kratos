//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_models/plasticity_models/hardening_rules/simo_exponential_hardening_rule.hpp"

namespace Kratos
{

  //*******************************CONSTRUCTOR******************************************
  //************************************************************************************

  SimoExponentialHardeningRule::SimoExponentialHardeningRule()
    :HardeningRule()
  {
    //Combined isotropic-kinematic 0<mTheta<1
    //Pure isotropic hardening mTheta=1;
    //Pure kinematic hardening mTheta=0;
  }


  //*******************************ASSIGMENT OPERATOR***********************************
  //************************************************************************************

  SimoExponentialHardeningRule& SimoExponentialHardeningRule::operator=(SimoExponentialHardeningRule const& rOther)
  {
    HardeningRule::operator=(rOther);
    return *this;
  }

  //*******************************COPY CONSTRUCTOR*************************************
  //************************************************************************************

  SimoExponentialHardeningRule::SimoExponentialHardeningRule(SimoExponentialHardeningRule const& rOther)
    :HardeningRule(rOther)
  {

  }


  //********************************CLONE***********************************************
  //************************************************************************************

  HardeningRule::Pointer SimoExponentialHardeningRule::Clone() const
  {
    return Kratos::make_shared<SimoExponentialHardeningRule>(*this);
  }


  //********************************DESTRUCTOR******************************************
  //************************************************************************************

  SimoExponentialHardeningRule::~SimoExponentialHardeningRule()
  {
  }

  /// Operations.

  //*******************************CALCULATE TOTAL HARDENING****************************
  //************************************************************************************

  double& SimoExponentialHardeningRule::CalculateHardening(const PlasticDataType& rVariables, double& rHardening)
  {
    KRATOS_TRY

    rHardening = this->CalculateAndAddIsotropicHardening(rVariables,rHardening);

    rHardening = this->CalculateAndAddKinematicHardening(rVariables,rHardening);

    return rHardening;

    KRATOS_CATCH(" ")
  }

  //*******************************CALCULATE ISOTROPIC HARDENING************************
  //************************************************************************************

  double& SimoExponentialHardeningRule::CalculateAndAddIsotropicHardening(const PlasticDataType& rVariables, double& rIsotropicHardening)
  {
    KRATOS_TRY

    const ModelDataType& rModelData = rVariables.GetModelData();

    //get values
    const double& rEquivalentPlasticStrain = rVariables.GetInternalVariables()[0];

    //linear hardening properties
    const Properties& rProperties      = rModelData.GetProperties();
    double  YieldStress                = rProperties[YIELD_STRESS];
    double  KinematicHardeningConstant = rProperties[KINEMATIC_HARDENING_MODULUS];

    //exponential saturation properties
    double  K_reference         =  rProperties[REFERENCE_HARDENING_MODULUS];
    double  K_infinity          =  rProperties[INFINITY_HARDENING_MODULUS];
    const double& Delta         =  rProperties[HARDENING_EXPONENT];


    double ThermalFactor        = this->CalculateThermalReferenceEffect(rVariables,ThermalFactor);
    YieldStress                *= ThermalFactor;
    K_reference                *= ThermalFactor;

    ThermalFactor               = this->CalculateThermalCurrentEffect(rVariables,ThermalFactor);
    K_infinity                 *= ThermalFactor;
    KinematicHardeningConstant *= ThermalFactor;


    //Linear Hardening rule: (mTheta = 1)
    rIsotropicHardening += YieldStress + mTheta * KinematicHardeningConstant * rEquivalentPlasticStrain;

    //Exponential Saturation:
    rIsotropicHardening += (K_infinity-K_reference) * (1.0 - exp( (-1.0) * Delta * rEquivalentPlasticStrain ) );

    return rIsotropicHardening;


    KRATOS_CATCH(" ")
  }

  //*******************************CALCULATE KINEMATIC HARDENING************************
  //************************************************************************************

  double& SimoExponentialHardeningRule::CalculateAndAddKinematicHardening(const PlasticDataType& rVariables, double& rKinematicHardening)
  {
    KRATOS_TRY

    const ModelDataType& rModelData = rVariables.GetModelData();

    //linear hardening properties
    double  KinematicHardeningConstant  =  rModelData.GetProperties()[KINEMATIC_HARDENING_MODULUS];

    double ThermalFactor        = this->CalculateThermalCurrentEffect(rVariables,ThermalFactor);
    KinematicHardeningConstant *= ThermalFactor;

    //Linear Hardening rule:
    rKinematicHardening  += (1.0 - mTheta) * KinematicHardeningConstant;

    return rKinematicHardening;

    KRATOS_CATCH(" ")
  }



  //*******************************CALCULATE HARDENING DERIVATIVE***********************
  //************************************************************************************

  double& SimoExponentialHardeningRule::CalculateDeltaHardening(const PlasticDataType& rVariables, double& rDeltaHardening)
  {
    KRATOS_TRY

    rDeltaHardening = this->CalculateAndAddDeltaIsotropicHardening(rVariables, rDeltaHardening);

    rDeltaHardening = this->CalculateAndAddDeltaKinematicHardening(rVariables, rDeltaHardening);

   return rDeltaHardening;

    KRATOS_CATCH(" ")
  }

  //***************************CALCULATE ISOTROPIC HARDENING DERIVATIVE*****************
  //************************************************************************************

  double& SimoExponentialHardeningRule::CalculateAndAddDeltaIsotropicHardening(const PlasticDataType& rVariables, double& rDeltaIsotropicHardening)
  {
    KRATOS_TRY

    const ModelDataType& rModelData = rVariables.GetModelData();

    //get values
    const double& rEquivalentPlasticStrain = rVariables.GetInternalVariables()[0];

    //linear hardening properties
    const Properties& rProperties       =  rModelData.GetProperties();
    double  KinematicHardeningConstant  =  rProperties[KINEMATIC_HARDENING_MODULUS];

    //exponential saturation properties
    double  K_reference           =  rProperties[REFERENCE_HARDENING_MODULUS];
    double  K_infinity            =  rProperties[INFINITY_HARDENING_MODULUS];
    const double& Delta           =  rProperties[HARDENING_EXPONENT];

    double ThermalFactor        = this->CalculateThermalReferenceEffect(rVariables,ThermalFactor);
    K_reference                *= ThermalFactor;

    ThermalFactor               = this->CalculateThermalCurrentEffect(rVariables,ThermalFactor);
    K_infinity                 *= ThermalFactor;
    KinematicHardeningConstant *= ThermalFactor;


    //Linear Hardening rule: (mTheta = 1)
    rDeltaIsotropicHardening += mTheta * KinematicHardeningConstant;

    //Exponential Saturation:
    rDeltaIsotropicHardening += Delta * (K_infinity-K_reference) * ( exp( (-1.0) * Delta * rEquivalentPlasticStrain ) );

    return rDeltaIsotropicHardening;

    KRATOS_CATCH(" ")
  }


  //***************************CALCULATE KINEMATIC HARDENING DERIVATIVE*****************
  //************************************************************************************

  double& SimoExponentialHardeningRule::CalculateAndAddDeltaKinematicHardening(const PlasticDataType& rVariables, double& rDeltaKinematicHardening)
  {
    KRATOS_TRY

    return rDeltaKinematicHardening;

    KRATOS_CATCH(" ")
  }


  //***************************CALCULATE TEMPERATURE EVOLUTION PROPERTIES***************
  //************************************************************************************


  double& SimoExponentialHardeningRule::CalculateThermalReferenceEffect(const PlasticDataType& rVariables, double& rThermalFactor)
  {
    KRATOS_TRY

    rThermalFactor = 1.0;
    return rThermalFactor;

    KRATOS_CATCH(" ")
  }

  //***************************CALCULATE TEMPERATURE EVOLUTION PROPERTIES***************
  //************************************************************************************

  double& SimoExponentialHardeningRule::CalculateThermalCurrentEffect(const PlasticDataType& rVariables, double& rThermalFactor)
  {
    KRATOS_TRY

    rThermalFactor = 1.0;
    return rThermalFactor;

    KRATOS_CATCH(" ")
  }


}  // namespace Kratos.
