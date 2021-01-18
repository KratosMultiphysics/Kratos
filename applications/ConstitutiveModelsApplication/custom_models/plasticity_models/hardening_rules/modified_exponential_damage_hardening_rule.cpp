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
#include "custom_models/plasticity_models/hardening_rules/modified_exponential_damage_hardening_rule.hpp"

namespace Kratos
{

  //*******************************CONSTRUCTOR******************************************
  //************************************************************************************

  ModifiedExponentialDamageHardeningRule::ModifiedExponentialDamageHardeningRule()
    :HardeningRule()
  {

  }


  //*******************************ASSIGMENT OPERATOR***********************************
  //************************************************************************************

  ModifiedExponentialDamageHardeningRule& ModifiedExponentialDamageHardeningRule::operator=(ModifiedExponentialDamageHardeningRule const& rOther)
  {
    HardeningRule::operator=(rOther);
    return *this;
  }

  //*******************************COPY CONSTRUCTOR*************************************
  //************************************************************************************

  ModifiedExponentialDamageHardeningRule::ModifiedExponentialDamageHardeningRule(ModifiedExponentialDamageHardeningRule const& rOther)
    :HardeningRule(rOther)
  {

  }


  //********************************CLONE***********************************************
  //************************************************************************************

  HardeningRule::Pointer ModifiedExponentialDamageHardeningRule::Clone() const
  {
    return Kratos::make_shared<ModifiedExponentialDamageHardeningRule>(*this);
  }


  //********************************DESTRUCTOR******************************************
  //************************************************************************************

  ModifiedExponentialDamageHardeningRule::~ModifiedExponentialDamageHardeningRule()
  {
  }

  /// Operations.

  //****************************** CALCULATE DAMAGE PARAMETER **************************
  //************************************************************************************

  double& ModifiedExponentialDamageHardeningRule::CalculateHardening(const PlasticDataType& rVariables, double& rHardening)
  {
    KRATOS_TRY

    const ModelDataType& rModelData  = rVariables.GetModelData();
    const Properties& rProperties    = rModelData.GetProperties();
    const double& rDamageThreshold   = rProperties[DAMAGE_THRESHOLD];
    const double& rResidualStrength  = rProperties[RESIDUAL_STRENGTH];
    const double& rSofteningSlope    = rProperties[SOFTENING_SLOPE];
    const double& rStateVariable     = rVariables.GetInternalVariables()[0];

    //Compute Damage variable from the internal historical variable
    rHardening  = 1.0-rDamageThreshold*(1.0-rResidualStrength)/rStateVariable;
    rHardening -= rResidualStrength*exp(-rSofteningSlope*(rStateVariable-rDamageThreshold));

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

  double& ModifiedExponentialDamageHardeningRule::CalculateDeltaHardening(const PlasticDataType& rVariables, double& rDeltaHardening)
  {
    KRATOS_TRY

    const ModelDataType& rModelData  = rVariables.GetModelData();
    const Properties& rProperties    = rModelData.GetProperties();
    const double& rDamageThreshold   = rProperties[DAMAGE_THRESHOLD];
    const double& rResidualStrength  = rProperties[RESIDUAL_STRENGTH];
    const double& rSofteningSlope    = rProperties[SOFTENING_SLOPE];
    const double& rStateVariable     = rVariables.GetInternalVariables()[0];

    //Damage derivative with respect to the internal historical variable
    rDeltaHardening  = rDamageThreshold*(1.0-rResidualStrength)/(rStateVariable*rStateVariable);
    rDeltaHardening += rResidualStrength*rSofteningSlope*exp(-rSofteningSlope*(rStateVariable-rDamageThreshold));

    if(rDeltaHardening < 0.0) rDeltaHardening = 0.0;

    return rDeltaHardening;

    KRATOS_CATCH(" ")

  }




}  // namespace Kratos.
