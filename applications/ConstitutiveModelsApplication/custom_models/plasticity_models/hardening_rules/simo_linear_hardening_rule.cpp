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
#include "custom_models/plasticity_models/hardening_rules/simo_linear_hardening_rule.hpp"


namespace Kratos
{

  //*******************************CONSTRUCTOR******************************************
  //************************************************************************************

  SimoLinearHardeningRule::SimoLinearHardeningRule()
    :SimoExponentialHardeningRule()
  {
  }


  //*******************************ASSIGMENT OPERATOR***********************************
  //************************************************************************************

  SimoLinearHardeningRule& SimoLinearHardeningRule::operator=(SimoLinearHardeningRule const& rOther)
  {
    SimoExponentialHardeningRule::operator=(rOther);
    return *this;
  }

  //*******************************COPY CONSTRUCTOR*************************************
  //************************************************************************************

  SimoLinearHardeningRule::SimoLinearHardeningRule(SimoLinearHardeningRule const& rOther)
    :SimoExponentialHardeningRule(rOther)
  {

  }


  //********************************CLONE***********************************************
  //************************************************************************************

  HardeningRule::Pointer SimoLinearHardeningRule::Clone() const
  {
    return Kratos::make_shared<SimoLinearHardeningRule>(*this);
  }


  //********************************DESTRUCTOR******************************************
  //************************************************************************************

  SimoLinearHardeningRule::~SimoLinearHardeningRule()
  {
  }

  /// Operations.


  //*******************************CALCULATE ISOTROPIC HARDENING************************
  //************************************************************************************

  double& SimoLinearHardeningRule::CalculateAndAddIsotropicHardening(const PlasticDataType& rVariables, double& rIsotropicHardening)
  {
    KRATOS_TRY

    const ModelDataType& rModelData = rVariables.GetModelData();

    //get values
    const double& rEquivalentPlasticStrain    = rVariables.GetInternalVariables()[0];

    //linear hardening properties
    const Properties& rMaterialProperties     = rModelData.GetMaterialProperties();
    const double& YieldStress                 = rMaterialProperties[YIELD_STRESS];
    const double& KinematicHardeningConstant  = rMaterialProperties[KINEMATIC_HARDENING_MODULUS];


    //Linear Hardening rule: (mTheta = 0)
    rIsotropicHardening  += YieldStress + (1.0 - mTheta) * KinematicHardeningConstant * rEquivalentPlasticStrain;


    return rIsotropicHardening;

    KRATOS_CATCH(" ")

  }


  //*******************************CALCULATE HARDENING DERIVATIVE***********************
  //************************************************************************************

  double& SimoLinearHardeningRule::CalculateDeltaHardening(const PlasticDataType& rVariables, double& rDeltaHardening)
  {
    KRATOS_TRY

    const ModelDataType& rModelData = rVariables.GetModelData();

    //linear hardening properties
    const double& KinematicHardeningConstant  =  rModelData.GetMaterialProperties()[KINEMATIC_HARDENING_MODULUS];

    //Linear Hardening rule: (mTheta = 0)
    rDeltaHardening  = (1.0 - mTheta) * KinematicHardeningConstant;

    return rDeltaHardening;

    KRATOS_CATCH(" ")

  }

  //***************************CALCULATE ISOTROPIC HARDENING DERIVATIVE*****************
  //************************************************************************************

  double& SimoLinearHardeningRule::CalculateAndAddDeltaIsotropicHardening(const PlasticDataType& rVariables, double& rDeltaIsotropicHardening)
  {
    KRATOS_TRY

    const ModelDataType& rModelData = rVariables.GetModelData();

    //linear hardening properties
    const double& KinematicHardeningConstant  =  rModelData.GetMaterialProperties()[KINEMATIC_HARDENING_MODULUS];

    //Linear Hardening rule: (mTheta = 0)
    rDeltaIsotropicHardening  += mTheta * KinematicHardeningConstant;

    return rDeltaIsotropicHardening;

    KRATOS_CATCH(" ")

  }



}  // namespace Kratos.
