//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                December 2016 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_models/plasticity_models/hardening_laws/linear_isotropic_kinematic_hardening_law.hpp"


namespace Kratos
{

  //*******************************CONSTRUCTOR******************************************
  //************************************************************************************

  LinearIsotropicKinematicHardeningLaw::LinearIsotropicKinematicHardeningLaw()
    :NonLinearIsotropicKinematicHardeningLaw()
  {
  }


  //*******************************ASSIGMENT OPERATOR***********************************
  //************************************************************************************

  LinearIsotropicKinematicHardeningLaw& LinearIsotropicKinematicHardeningLaw::operator=(LinearIsotropicKinematicHardeningLaw const& rOther)
  {
    NonLinearIsotropicKinematicHardeningLaw::operator=(rOther);
    return *this;
  }

  //*******************************COPY CONSTRUCTOR*************************************
  //************************************************************************************

  LinearIsotropicKinematicHardeningLaw::LinearIsotropicKinematicHardeningLaw(LinearIsotropicKinematicHardeningLaw const& rOther)
    :NonLinearIsotropicKinematicHardeningLaw(rOther)
  {

  }


  //********************************CLONE***********************************************
  //************************************************************************************

  HardeningLaw::Pointer LinearIsotropicKinematicHardeningLaw::Clone() const
  {
    return ( HardeningLaw::Pointer(new LinearIsotropicKinematicHardeningLaw(*this)) );
  }


  //********************************DESTRUCTOR******************************************
  //************************************************************************************

  LinearIsotropicKinematicHardeningLaw::~LinearIsotropicKinematicHardeningLaw()
  {
  }

  /// Operations.

  
  //*******************************CALCULATE ISOTROPIC HARDENING************************
  //************************************************************************************

  double& LinearIsotropicKinematicHardeningLaw::CalculateAndAddIsotropicHardening(const PlasticDataType& rVariables, double &rIsotropicHardening)
  {
    KRATOS_TRY

    const ModelData& rModelData = rVariables.GetModelData();

    //get values
    const double& rEquivalentPlasticStrain = rVariables.GetEquivalentPlasticStrain();

    //linear hardening properties
    const Properties& rMaterialProperties     =  rModelData.GetMaterialProperties();
    const double& YieldStress                 =  rMaterialProperties[YIELD_STRESS];
    const double& KinematicHardeningConstant  =  rMaterialProperties[KINEMATIC_HARDENING_MODULUS];
	

    //Linear Hardening law: (mTheta = 0)
    rIsotropicHardening  += YieldStress + (1.0 - mTheta) * KinematicHardeningConstant * rEquivalentPlasticStrain;
	
	
    return rIsotropicHardening;

    KRATOS_CATCH(" ")
    
  }


  //*******************************CALCULATE HARDENING DERIVATIVE***********************
  //************************************************************************************

  double& LinearIsotropicKinematicHardeningLaw::CalculateDeltaHardening(const PlasticDataType& rVariables, double &rDeltaHardening)
  {
    KRATOS_TRY
      
    const ModelData& rModelData = rVariables.GetModelData();

    //linear hardening properties
    const double& KinematicHardeningConstant  =  rModelData.GetMaterialProperties()[KINEMATIC_HARDENING_MODULUS];
	
    //Linear Hardening law: (mTheta = 0)
    rDeltaHardening  = (1.0 - mTheta) * KinematicHardeningConstant;
		
    return rDeltaHardening;

    KRATOS_CATCH(" ")
    
  }

  //***************************CALCULATE ISOTROPIC HARDENING DERIVATIVE*****************
  //************************************************************************************

  double& LinearIsotropicKinematicHardeningLaw::CalculateAndAddDeltaIsotropicHardening(const PlasticDataType& rVariables, double &rDeltaIsotropicHardening)
  {
    KRATOS_TRY
      
    const ModelData& rModelData = rVariables.GetModelData();

    //linear hardening properties
    const double& KinematicHardeningConstant  =  rModelData.GetMaterialProperties()[KINEMATIC_HARDENING_MODULUS];
	
    //Linear Hardening law: (mTheta = 0)
    rDeltaIsotropicHardening  += mTheta * KinematicHardeningConstant;
	
    return rDeltaIsotropicHardening;

    KRATOS_CATCH(" ")
    
  }



}  // namespace Kratos.
