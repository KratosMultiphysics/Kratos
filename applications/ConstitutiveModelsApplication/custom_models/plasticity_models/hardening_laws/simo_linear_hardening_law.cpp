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
#include "custom_models/plasticity_models/hardening_laws/simo_linear_hardening_law.hpp"


namespace Kratos
{

  //*******************************CONSTRUCTOR******************************************
  //************************************************************************************

  SimoLinearHardeningLaw::SimoLinearHardeningLaw()
    :SimoExponentialHardeningLaw()
  {
  }


  //*******************************ASSIGMENT OPERATOR***********************************
  //************************************************************************************

  SimoLinearHardeningLaw& SimoLinearHardeningLaw::operator=(SimoLinearHardeningLaw const& rOther)
  {
    SimoExponentialHardeningLaw::operator=(rOther);
    return *this;
  }

  //*******************************COPY CONSTRUCTOR*************************************
  //************************************************************************************

  SimoLinearHardeningLaw::SimoLinearHardeningLaw(SimoLinearHardeningLaw const& rOther)
    :SimoExponentialHardeningLaw(rOther)
  {

  }


  //********************************CLONE***********************************************
  //************************************************************************************

  HardeningLaw::Pointer SimoLinearHardeningLaw::Clone() const
  {
    return ( HardeningLaw::Pointer(new SimoLinearHardeningLaw(*this)) );
  }


  //********************************DESTRUCTOR******************************************
  //************************************************************************************

  SimoLinearHardeningLaw::~SimoLinearHardeningLaw()
  {
  }

  /// Operations.

  
  //*******************************CALCULATE ISOTROPIC HARDENING************************
  //************************************************************************************

  double& SimoLinearHardeningLaw::CalculateAndAddIsotropicHardening(const PlasticDataType& rVariables, double& rIsotropicHardening)
  {
    KRATOS_TRY

    const ModelDataType& rModelData = rVariables.GetModelData();

    //get values
    const double& rEquivalentPlasticStrain    = rVariables.GetInternalVariables()[0];

    //linear hardening properties
    const Properties& rMaterialProperties     = rModelData.GetMaterialProperties();
    const double& YieldStress                 = rMaterialProperties[YIELD_STRESS];
    const double& KinematicHardeningConstant  = rMaterialProperties[KINEMATIC_HARDENING_MODULUS];
	

    //Linear Hardening law: (mTheta = 0)
    rIsotropicHardening  += YieldStress + (1.0 - mTheta) * KinematicHardeningConstant * rEquivalentPlasticStrain;
	
	
    return rIsotropicHardening;

    KRATOS_CATCH(" ")
    
  }


  //*******************************CALCULATE HARDENING DERIVATIVE***********************
  //************************************************************************************

  double& SimoLinearHardeningLaw::CalculateDeltaHardening(const PlasticDataType& rVariables, double& rDeltaHardening)
  {
    KRATOS_TRY
      
    const ModelDataType& rModelData = rVariables.GetModelData();

    //linear hardening properties
    const double& KinematicHardeningConstant  =  rModelData.GetMaterialProperties()[KINEMATIC_HARDENING_MODULUS];
	
    //Linear Hardening law: (mTheta = 0)
    rDeltaHardening  = (1.0 - mTheta) * KinematicHardeningConstant;
		
    return rDeltaHardening;

    KRATOS_CATCH(" ")
    
  }

  //***************************CALCULATE ISOTROPIC HARDENING DERIVATIVE*****************
  //************************************************************************************

  double& SimoLinearHardeningLaw::CalculateAndAddDeltaIsotropicHardening(const PlasticDataType& rVariables, double& rDeltaIsotropicHardening)
  {
    KRATOS_TRY
      
    const ModelDataType& rModelData = rVariables.GetModelData();

    //linear hardening properties
    const double& KinematicHardeningConstant  =  rModelData.GetMaterialProperties()[KINEMATIC_HARDENING_MODULUS];
	
    //Linear Hardening law: (mTheta = 0)
    rDeltaIsotropicHardening  += mTheta * KinematicHardeningConstant;
	
    return rDeltaIsotropicHardening;

    KRATOS_CATCH(" ")
    
  }



}  // namespace Kratos.
