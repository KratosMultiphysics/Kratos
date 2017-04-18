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
#include "custom_models/plasticity_models/hardening_laws/exponential_damage_hardening_law.hpp"

namespace Kratos
{

  //*******************************CONSTRUCTOR******************************************
  //************************************************************************************

  ExponentialDamageHardeningLaw::ExponentialDamageHardeningLaw()
    :HardeningLaw()
  {
       
  }


  //*******************************ASSIGMENT OPERATOR***********************************
  //************************************************************************************

  ExponentialDamageHardeningLaw& ExponentialDamageHardeningLaw::operator=(ExponentialDamageHardeningLaw const& rOther)
  {
    HardeningLaw::operator=(rOther);
    return *this;
  }

  //*******************************COPY CONSTRUCTOR*************************************
  //************************************************************************************

  ExponentialDamageHardeningLaw::ExponentialDamageHardeningLaw(ExponentialDamageHardeningLaw const& rOther)
    :HardeningLaw(rOther)
  {

  }


  //********************************CLONE***********************************************
  //************************************************************************************

  HardeningLaw::Pointer ExponentialDamageHardeningLaw::Clone() const
  {
    return ( HardeningLaw::Pointer(new ExponentialDamageHardeningLaw(*this)) );
  }


  //********************************DESTRUCTOR******************************************
  //************************************************************************************

  ExponentialDamageHardeningLaw::~ExponentialDamageHardeningLaw()
  {
  }

  /// Operations.

  //****************************** CALCULATE DAMAGE PARAMETER **************************
  //************************************************************************************

  double& ExponentialDamageHardeningLaw::CalculateHardening(const PlasticDataType& rVariables, double &rHardening)
  {
    KRATOS_TRY
   
    const ModelDataType& rModelData       = rVariables.GetModelData();
    const Properties& rMaterialProperties = rModelData.GetMaterialProperties();
    const double& rFractureEnergy         = rMaterialProperties[FRACTURE_ENERGY];
    const double& rDamageThreshold        = rMaterialProperties[DAMAGE_THRESHOLD];
    const double& rCharacteristicSize     = rModelData.GetCharacteristicSize();
    const double& rStateVariable          = rVariables.GetDeltaInternalVariables()[0];

    
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

  double& ExponentialDamageHardeningLaw::CalculateDeltaHardening(const PlasticDataType& rVariables, double &rDeltaHardening)
  {
    KRATOS_TRY
      
    const ModelDataType& rModelData       = rVariables.GetModelData();
    const Properties& rMaterialProperties = rModelData.GetMaterialProperties();
    const double& rFractureEnergy         = rMaterialProperties[FRACTURE_ENERGY];
    const double& rDamageThreshold        = rMaterialProperties[DAMAGE_THRESHOLD];
    const double& rCharacteristicSize     = rModelData.GetCharacteristicSize();
    const double& rStateVariable          = rVariables.GetDeltaInternalVariables()[0];

    double A = 1.0/(rFractureEnergy/(rCharacteristicSize*rDamageThreshold*rDamageThreshold)-0.5);

    if(A < 0.0) A = 0.0;
    
    //Damage derivative with respect to the internal historical variable
    rDeltaHardening = (rDamageThreshold + A*rStateVariable)/(rStateVariable*rStateVariable)*exp(A*(1.0-rStateVariable/rDamageThreshold));

    if(rDeltaHardening < 0.0) rDeltaHardening = 0.0;
    
    return rDeltaHardening;
    
    KRATOS_CATCH(" ")
	  
  }




}  // namespace Kratos.
