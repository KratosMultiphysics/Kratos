//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                  IPouplana $
//   Last modified by:    $Co-Author:             JMCarbonell $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_models/plasticity_models/hardening_rules/exponential_damage_hardening_rule.hpp"

namespace Kratos
{

  //*******************************CONSTRUCTOR******************************************
  //************************************************************************************

  ExponentialDamageHardeningRule::ExponentialDamageHardeningRule()
    :HardeningRule()
  {

  }


  //*******************************ASSIGMENT OPERATOR***********************************
  //************************************************************************************

  ExponentialDamageHardeningRule& ExponentialDamageHardeningRule::operator=(ExponentialDamageHardeningRule const& rOther)
  {
    HardeningRule::operator=(rOther);
    return *this;
  }

  //*******************************COPY CONSTRUCTOR*************************************
  //************************************************************************************

  ExponentialDamageHardeningRule::ExponentialDamageHardeningRule(ExponentialDamageHardeningRule const& rOther)
    :HardeningRule(rOther)
  {

  }


  //********************************CLONE***********************************************
  //************************************************************************************

  HardeningRule::Pointer ExponentialDamageHardeningRule::Clone() const
  {
    return Kratos::make_shared<ExponentialDamageHardeningRule>(*this);
  }


  //********************************DESTRUCTOR******************************************
  //************************************************************************************

  ExponentialDamageHardeningRule::~ExponentialDamageHardeningRule()
  {
  }

  /// Operations.

  //****************************** CALCULATE DAMAGE PARAMETER **************************
  //************************************************************************************

  double& ExponentialDamageHardeningRule::CalculateHardening(const PlasticDataType& rVariables, double& rHardening)
  {
    KRATOS_TRY

    const ModelDataType& rModelData    = rVariables.GetModelData();
    const Properties& rProperties      = rModelData.GetProperties();
    const double& rFractureEnergy      = rProperties[FRACTURE_ENERGY];
    const double& rDamageThreshold     = rProperties[DAMAGE_THRESHOLD];
    const double& rCharacteristicSize  = rModelData.GetCharacteristicSize();
    const double& rStateVariable       = rVariables.GetInternalVariables()[0];


    double A = 1.0/(rFractureEnergy/(rCharacteristicSize*rDamageThreshold*rDamageThreshold)-0.5);

    if(A < 0.0) A = 0.0;

    //Compute Damage variable from the internal historical variable
    rHardening = 1.0-rDamageThreshold/rStateVariable*exp(A*(1.0-rStateVariable/rDamageThreshold));

    if(rHardening < 0.0)
      {
        rHardening = 0.0;
      }
    else if(rHardening > 1.0)
      {
        rHardening = 1.0;
      }

    return rHardening;

    KRATOS_CATCH(" ")

  }


  //***************************** CALCULATE DAMAGE DERIVATIVE **************************
  //************************************************************************************

  double& ExponentialDamageHardeningRule::CalculateDeltaHardening(const PlasticDataType& rVariables, double& rDeltaHardening)
  {
    KRATOS_TRY

    const ModelDataType& rModelData    = rVariables.GetModelData();
    const Properties& rProperties      = rModelData.GetProperties();
    const double& rFractureEnergy      = rProperties[FRACTURE_ENERGY];
    const double& rDamageThreshold     = rProperties[DAMAGE_THRESHOLD];
    const double& rCharacteristicSize  = rModelData.GetCharacteristicSize();
    const double& rStateVariable       = rVariables.GetInternalVariables()[0];

    double A = 1.0/(rFractureEnergy/(rCharacteristicSize*rDamageThreshold*rDamageThreshold)-0.5);

    if(A < 0.0) A = 0.0;

    //Damage derivative with respect to the internal historical variable
    rDeltaHardening = (rDamageThreshold + A*rStateVariable)/(rStateVariable*rStateVariable)*exp(A*(1.0-rStateVariable/rDamageThreshold));

    if(rDeltaHardening < 0.0) rDeltaHardening = 0.0;

    return rDeltaHardening;

    KRATOS_CATCH(" ")

  }




}  // namespace Kratos.
