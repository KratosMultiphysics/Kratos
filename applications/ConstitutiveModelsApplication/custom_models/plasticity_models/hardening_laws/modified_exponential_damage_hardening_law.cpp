//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                  IPouplana $
//   Last modified by:    $Co-Author:             JMCarbonell $
//   Date:                $Date:                December 2016 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_models/plasticity_models/hardening_laws/modified_exponential_damage_hardening_law.hpp"

namespace Kratos
{

  //*******************************CONSTRUCTOR******************************************
  //************************************************************************************

  ModifiedExponentialDamageHardeningLaw::ModifiedExponentialDamageHardeningLaw()
    :HardeningLaw()
  {
       
  }


  //*******************************ASSIGMENT OPERATOR***********************************
  //************************************************************************************

  ModifiedExponentialDamageHardeningLaw& ModifiedExponentialDamageHardeningLaw::operator=(ModifiedExponentialDamageHardeningLaw const& rOther)
  {
    HardeningLaw::operator=(rOther);
    return *this;
  }

  //*******************************COPY CONSTRUCTOR*************************************
  //************************************************************************************

  ModifiedExponentialDamageHardeningLaw::ModifiedExponentialDamageHardeningLaw(ModifiedExponentialDamageHardeningLaw const& rOther)
    :HardeningLaw(rOther)
  {

  }


  //********************************CLONE***********************************************
  //************************************************************************************

  HardeningLaw::Pointer ModifiedExponentialDamageHardeningLaw::Clone() const
  {
    return ( HardeningLaw::Pointer(new ModifiedExponentialDamageHardeningLaw(*this)) );
  }


  //********************************DESTRUCTOR******************************************
  //************************************************************************************

  ModifiedExponentialDamageHardeningLaw::~ModifiedExponentialDamageHardeningLaw()
  {
  }

  /// Operations.

  //****************************** CALCULATE DAMAGE PARAMETER **************************
  //************************************************************************************

  double& ModifiedExponentialDamageHardeningLaw::CalculateHardening(const PlasticDataType& rVariables, double &rHardening)
  {
    KRATOS_TRY
      
    const ModelData& rModelData           = rVariables.GetModelData();
    const Properties& rMaterialProperties = rModelData.GetMaterialProperties();
    const double& rDamageThreshold        = rMaterialProperties[DAMAGE_THRESHOLD];
    const double& rResidualStrength       = rMaterialProperties[RESIDUAL_STRENGTH];
    const double& rSofteningSlope         = rMaterialProperties[SOFTENING_SLOPE];
    const double& rStateVariable          = rVariables.GetDeltaInternalVariables()[0];
    
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

  double& ModifiedExponentialDamageHardeningLaw::CalculateDeltaHardening(const PlasticDataType& rVariables, double &rDeltaHardening)
  {
    KRATOS_TRY
      
    const ModelData& rModelData           = rVariables.GetModelData();
    const Properties& rMaterialProperties = rModelData.GetMaterialProperties();
    const double& rDamageThreshold        = rMaterialProperties[DAMAGE_THRESHOLD];
    const double& rResidualStrength       = rMaterialProperties[RESIDUAL_STRENGTH];
    const double& rSofteningSlope         = rMaterialProperties[SOFTENING_SLOPE];
    const double& rStateVariable          = rVariables.GetDeltaInternalVariables()[0];
    
    //Damage derivative with respect to the internal historical variable
    rDeltaHardening  = rDamageThreshold*(1.0-rResidualStrength)/(rStateVariable*rStateVariable);
    rDeltaHardening += rResidualStrength*rSofteningSlope*exp(-rSofteningSlope*(rStateVariable-rDamageThreshold));

    if(rDeltaHardening < 0.0) rDeltaHardening = 0.0;
    
    return rDeltaHardening;
    
    KRATOS_CATCH(" ")
    
  }




}  // namespace Kratos.
