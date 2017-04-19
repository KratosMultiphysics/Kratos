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
#include "custom_models/plasticity_models/hardening_laws/non_linear_isotropic_kinematic_hardening_law.hpp"

namespace Kratos
{

  //*******************************CONSTRUCTOR******************************************
  //************************************************************************************

  NonLinearIsotropicKinematicHardeningLaw::NonLinearIsotropicKinematicHardeningLaw()
    :HardeningLaw()
  {
    //Combined isotropic-kinematic 0<mTheta<1
    //Pure isotropic hardening mTheta=1;  
    //Pure kinematic hardening mTheta=0;    
  }


  //*******************************ASSIGMENT OPERATOR***********************************
  //************************************************************************************

  NonLinearIsotropicKinematicHardeningLaw& NonLinearIsotropicKinematicHardeningLaw::operator=(NonLinearIsotropicKinematicHardeningLaw const& rOther)
  {
    HardeningLaw::operator=(rOther);
    return *this;
  }

  //*******************************COPY CONSTRUCTOR*************************************
  //************************************************************************************

  NonLinearIsotropicKinematicHardeningLaw::NonLinearIsotropicKinematicHardeningLaw(NonLinearIsotropicKinematicHardeningLaw const& rOther)
    :HardeningLaw(rOther)
  {

  }


  //********************************CLONE***********************************************
  //************************************************************************************

  HardeningLaw::Pointer NonLinearIsotropicKinematicHardeningLaw::Clone() const
  {
    return ( HardeningLaw::Pointer(new NonLinearIsotropicKinematicHardeningLaw(*this)) );
  }


  //********************************DESTRUCTOR******************************************
  //************************************************************************************

  NonLinearIsotropicKinematicHardeningLaw::~NonLinearIsotropicKinematicHardeningLaw()
  {
  }

  /// Operations.

  //*******************************CALCULATE TOTAL HARDENING****************************
  //************************************************************************************

  double& NonLinearIsotropicKinematicHardeningLaw::CalculateHardening(const PlasticDataType& rVariables, double &rHardening)
  {
    KRATOS_TRY

    rHardening = this->CalculateAndAddIsotropicHardening(rVariables,rHardening);

    rHardening = this->CalculateAndAddKinematicHardening(rVariables,rHardening);
	
    return rHardening;
	
    KRATOS_CATCH(" ")
  }
  
  //*******************************CALCULATE ISOTROPIC HARDENING************************
  //************************************************************************************

  double& NonLinearIsotropicKinematicHardeningLaw::CalculateAndAddIsotropicHardening(const PlasticDataType& rVariables, double &rIsotropicHardening)
  {
    KRATOS_TRY

    const ModelDataType& rModelData = rVariables.GetModelData();
      
    //get values   
    const double& rEquivalentPlasticStrain = rVariables.GetInternalVariables()[0];
    const double& rTemperature             = rModelData.GetTemperature();

    //linear hardening properties
    const Properties& rMaterialProperties  = rModelData.GetMaterialProperties();
    double  YieldStress                    = rMaterialProperties[YIELD_STRESS];
    double  KinematicHardeningConstant     = rMaterialProperties[KINEMATIC_HARDENING_MODULUS];
	
    //exponential saturation properties
    double  K_reference         =  rMaterialProperties[REFERENCE_HARDENING_MODULUS];
    double  K_infinity          =  rMaterialProperties[INFINITY_HARDENING_MODULUS];
    const double& Delta         =  rMaterialProperties[HARDENING_EXPONENT];


    double ThermalFactor        = this->CalculateThermalReferenceEffect(rVariables,rTemperature,ThermalFactor);
    YieldStress                *= ThermalFactor;
    K_reference                *= ThermalFactor;

    ThermalFactor               = this->CalculateThermalCurrentEffect(rVariables,rTemperature,ThermalFactor);
    K_infinity                 *= ThermalFactor;
    KinematicHardeningConstant *= ThermalFactor;


    //Linear Hardening law: (mTheta = 1)
    rIsotropicHardening += YieldStress + mTheta * KinematicHardeningConstant * rEquivalentPlasticStrain;
	
    //Exponential Saturation:
    rIsotropicHardening += (K_infinity-K_reference) * (1.0 - exp( (-1.0) * Delta * rEquivalentPlasticStrain ) );
	
    return rIsotropicHardening;	

	
    KRATOS_CATCH(" ")
  }

  //*******************************CALCULATE KINEMATIC HARDENING************************
  //************************************************************************************

  double& NonLinearIsotropicKinematicHardeningLaw::CalculateAndAddKinematicHardening(const PlasticDataType& rVariables, double &rKinematicHardening)
  {
    KRATOS_TRY

    const ModelDataType& rModelData = rVariables.GetModelData();
      
    //get values
    const double& rTemperature          =  rModelData.GetTemperature();

    //linear hardening properties
    double  KinematicHardeningConstant  =  rModelData.GetMaterialProperties()[KINEMATIC_HARDENING_MODULUS];

    double ThermalFactor        = this->CalculateThermalCurrentEffect(rVariables,rTemperature,ThermalFactor);
    KinematicHardeningConstant *= ThermalFactor;

    //Linear Hardening law:
    rKinematicHardening  += (1.0 - mTheta) * KinematicHardeningConstant;
	
    return rKinematicHardening;
 	
    KRATOS_CATCH(" ")
  }



  //*******************************CALCULATE HARDENING DERIVATIVE***********************
  //************************************************************************************

  double& NonLinearIsotropicKinematicHardeningLaw::CalculateDeltaHardening(const PlasticDataType& rVariables, double &rDeltaHardening)
  {
    KRATOS_TRY

    rDeltaHardening = this->CalculateAndAddDeltaIsotropicHardening(rVariables, rDeltaHardening);

    rDeltaHardening = this->CalculateAndAddDeltaKinematicHardening(rVariables, rDeltaHardening);

   return rDeltaHardening;	
	
    KRATOS_CATCH(" ")
  }

  //***************************CALCULATE ISOTROPIC HARDENING DERIVATIVE*****************
  //************************************************************************************

  double& NonLinearIsotropicKinematicHardeningLaw::CalculateAndAddDeltaIsotropicHardening(const PlasticDataType& rVariables, double &rDeltaIsotropicHardening)
  {
    KRATOS_TRY

    const ModelDataType& rModelData = rVariables.GetModelData();

    //get values
    const double& rEquivalentPlasticStrain = rVariables.GetInternalVariables()[0];
    const double& rTemperature             = rModelData.GetTemperature();

    //linear hardening properties
    const Properties& rMaterialProperties  =  rModelData.GetMaterialProperties();
    double  KinematicHardeningConstant     =  rMaterialProperties[KINEMATIC_HARDENING_MODULUS];
	
    //exponential saturation properties
    double  K_reference           =  rMaterialProperties[REFERENCE_HARDENING_MODULUS];
    double  K_infinity            =  rMaterialProperties[INFINITY_HARDENING_MODULUS];
    const double& Delta           =  rMaterialProperties[HARDENING_EXPONENT];

    double ThermalFactor        = this->CalculateThermalReferenceEffect(rVariables,rTemperature,ThermalFactor);
    K_reference                *= ThermalFactor;

    ThermalFactor               = this->CalculateThermalCurrentEffect(rVariables,rTemperature,ThermalFactor);
    K_infinity                 *= ThermalFactor;
    KinematicHardeningConstant *= ThermalFactor;


    //Linear Hardening law: (mTheta = 1)
    rDeltaIsotropicHardening += mTheta * KinematicHardeningConstant;
	
    //Exponential Saturation:
    rDeltaIsotropicHardening += Delta * (K_infinity-K_reference) * ( exp( (-1.0) * Delta * rEquivalentPlasticStrain ) );
	
    return rDeltaIsotropicHardening;	

    KRATOS_CATCH(" ")
  }


  //***************************CALCULATE KINEMATIC HARDENING DERIVATIVE*****************
  //************************************************************************************

  double& NonLinearIsotropicKinematicHardeningLaw::CalculateAndAddDeltaKinematicHardening(const PlasticDataType& rVariables, double &rDeltaKinematicHardening)
  {
    KRATOS_TRY

    return rDeltaKinematicHardening;
	
    KRATOS_CATCH(" ")
  }


  //***************************CALCULATE TEMPERATURE EVOLUTION PROPERTIES***************
  //************************************************************************************


  double& NonLinearIsotropicKinematicHardeningLaw::CalculateThermalReferenceEffect(const PlasticDataType& rVariables, const double &rTemperature, double& rThermalFactor)
  {
    KRATOS_TRY

    rThermalFactor = 1.0;  
    return rThermalFactor;
	
    KRATOS_CATCH(" ")
  }

  //***************************CALCULATE TEMPERATURE EVOLUTION PROPERTIES***************
  //************************************************************************************

  double& NonLinearIsotropicKinematicHardeningLaw::CalculateThermalCurrentEffect(const PlasticDataType& rVariables, const double &rTemperature, double& rThermalFactor)
  {
    KRATOS_TRY

    rThermalFactor = 1.0;  
    return rThermalFactor;
 	
    KRATOS_CATCH(" ")
  }


}  // namespace Kratos.
