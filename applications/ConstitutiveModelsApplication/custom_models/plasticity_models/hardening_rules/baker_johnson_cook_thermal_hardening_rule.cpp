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
#include "custom_models/plasticity_models/hardening_rules/baker_johnson_cook_thermal_hardening_rule.hpp"

namespace Kratos
{

  //*******************************CONSTRUCTOR******************************************
  //************************************************************************************

  BakerJohnsonCookThermalHardeningRule::BakerJohnsonCookThermalHardeningRule()
    :HardeningRule()
  {

  }


  //*******************************ASSIGMENT OPERATOR***********************************
  //************************************************************************************

  BakerJohnsonCookThermalHardeningRule& BakerJohnsonCookThermalHardeningRule::operator=(BakerJohnsonCookThermalHardeningRule const& rOther)
  {
    HardeningRule::operator=(rOther);
    return *this;
  }

  //*******************************COPY CONSTRUCTOR*************************************
  //************************************************************************************

  BakerJohnsonCookThermalHardeningRule::BakerJohnsonCookThermalHardeningRule(BakerJohnsonCookThermalHardeningRule const& rOther)
    :HardeningRule(rOther)
  {

  }


  //********************************CLONE***********************************************
  //************************************************************************************

  HardeningRule::Pointer BakerJohnsonCookThermalHardeningRule::Clone() const
  {
    return Kratos::make_shared<BakerJohnsonCookThermalHardeningRule>(*this);
  }


  //********************************DESTRUCTOR******************************************
  //************************************************************************************

  BakerJohnsonCookThermalHardeningRule::~BakerJohnsonCookThermalHardeningRule()
  {
  }

  /// Operations.

  //*******************************CALCULATE TOTAL HARDENING****************************
  //************************************************************************************

  double& BakerJohnsonCookThermalHardeningRule::CalculateHardening(const PlasticDataType& rVariables, double& rHardening)
  {
    KRATOS_TRY

    const ModelDataType& rModelData = rVariables.GetModelData();

    //get values
    const double& rRateFactor              = rVariables.GetRateFactor();
    const double& rEquivalentPlasticStrain = rVariables.GetInternalVariables()[0];
    const double& rDeltaGamma              = rVariables.GetDeltaInternalVariables()[0];

    const double& rTemperature             = rModelData.GetTemperature();
    const double& rDeltaTime               = rModelData.GetProcessInfo()[DELTA_TIME];

    //Constant Parameters of the -- Johnson and Cook --:
    const Properties& rProperties  = rModelData.GetProperties();
    const double& K = rProperties[JC_PARAMETER_K];
    const double& C = rProperties[JC_PARAMETER_C];

    const double& n = rProperties[JC_PARAMETER_n];
    const double& m = rProperties[JC_PARAMETER_m];

    const double& rReferenceTemperature = rProperties[REFERENCE_TEMPERATURE];
    const double& rMeldTemperature      = rProperties[MELD_TEMPERATURE];
    const double& rPlasticStrainRate    = rProperties[PLASTIC_STRAIN_RATE];

    if(rTemperature - rReferenceTemperature < 0){
      std::cout<<" Initial Temperature conditions not defined properly ("<<rTemperature<<" < "<<rReferenceTemperature<<")"<<std::endl;

    }

    double NormalizedTemperature =  pow( (rTemperature)/(rMeldTemperature), m );

    double Kv = K * exp((-1)*NormalizedTemperature);
    double nv = n * exp((-1)*NormalizedTemperature);


    rHardening   = ( Kv * pow(rEquivalentPlasticStrain, nv) );

    if( rRateFactor != 0 ){

      if( rDeltaGamma == 0 )
	std::cout<<" Something is wrong in the Baker_Johnson_Cook_hardening variables supplied "<<std::endl;

      rHardening  *= (1.0 + rRateFactor * C * std::log( (rDeltaGamma * sqrt(2.0/3.0))/(rPlasticStrainRate * rDeltaTime) ) );

    }

    return rHardening;


    KRATOS_CATCH(" ")
  }


  //*******************************CALCULATE HARDENING DERIVATIVE***********************
  //************************************************************************************

  double& BakerJohnsonCookThermalHardeningRule::CalculateDeltaHardening(const PlasticDataType& rVariables, double& rDeltaHardening)
  {
    KRATOS_TRY

    const ModelDataType& rModelData = rVariables.GetModelData();

    //get values
    const double& rRateFactor              = rVariables.GetRateFactor();
    const double& rEquivalentPlasticStrain = rVariables.GetInternalVariables()[0];
    const double& rDeltaGamma              = rVariables.GetDeltaInternalVariables()[0];

    const double& rTemperature             = rModelData.GetTemperature();
    const double& rDeltaTime               = rModelData.GetProcessInfo()[DELTA_TIME];

    //Constant Parameters of the -- Johnson and Cook --:
    const Properties& rProperties  = rModelData.GetProperties();
    const double& K = rProperties[JC_PARAMETER_K];
    const double& C = rProperties[JC_PARAMETER_C];

    const double& n = rProperties[JC_PARAMETER_n];
    const double& m = rProperties[JC_PARAMETER_m];

    const double& rReferenceTemperature = rProperties[REFERENCE_TEMPERATURE];
    const double& rMeldTemperature      = rProperties[MELD_TEMPERATURE];
    const double& rPlasticStrainRate    = rProperties[PLASTIC_STRAIN_RATE];

    if(rTemperature - rReferenceTemperature < 0){
      std::cout<<" Initial Temperature conditions not defined properly ("<<rTemperature<<" < "<<rReferenceTemperature<<")"<<std::endl;

    }

    double NormalizedTemperature =  pow( (rTemperature)/(rMeldTemperature), m );

    double Kv = K * exp((-1)*NormalizedTemperature);
    double nv = n * exp((-1)*NormalizedTemperature);


    rDeltaHardening  = ( nv * Kv * pow(rEquivalentPlasticStrain, nv-1) );

    if( rRateFactor != 0 ){

      if( rDeltaGamma == 0 )
	std::cout<<" Something is wrong in the Baker_Johnson_Cook_hardening variables supplied "<<std::endl;

      rDeltaHardening *= (1.0 + rRateFactor * C * std::log( (rDeltaGamma * sqrt(2.0/3.0))/(rPlasticStrainRate * rDeltaTime) ) );

      rDeltaHardening += rRateFactor * ( sqrt(3.0/2.0) * ( Kv * pow( rEquivalentPlasticStrain, nv ) ) * C / rDeltaGamma );
    }

    return rDeltaHardening;


    KRATOS_CATCH(" ")
  }


  //***************************CALCULATE HARDENING DERIVATIVE TEMPERATURE***************
  //************************************************************************************

  double& BakerJohnsonCookThermalHardeningRule::CalculateDeltaThermalHardening(const PlasticDataType& rVariables, double& rDeltaThermalHardening)
  {
    KRATOS_TRY

    const ModelDataType& rModelData = rVariables.GetModelData();

    //get values
    const double& rRateFactor              = rVariables.GetRateFactor();
    const double& rEquivalentPlasticStrain = rVariables.GetInternalVariables()[0];
    const double& rDeltaGamma              = rVariables.GetDeltaInternalVariables()[0];

    const double& rTemperature             = rModelData.GetTemperature();
    const double& rDeltaTime               = rModelData.GetProcessInfo()[DELTA_TIME];

    //Constant Parameters of the -- Baker Johnson and Cook --:
    const Properties& rProperties  = rModelData.GetProperties();
    const double& K = rProperties[JC_PARAMETER_K];
    const double& C = rProperties[JC_PARAMETER_C];

    const double& n = rProperties[JC_PARAMETER_n];
    const double& m = rProperties[JC_PARAMETER_m];

    const double& rReferenceTemperature = rProperties[REFERENCE_TEMPERATURE];
    const double& rMeldTemperature      = rProperties[MELD_TEMPERATURE];
    const double& rPlasticStrainRate    = rProperties[PLASTIC_STRAIN_RATE];

    if(rTemperature - rReferenceTemperature < 0){
      std::cout<<" Initial Temperature conditions not defined properly ("<<rTemperature<<" < "<<rReferenceTemperature<<")"<<std::endl;

    }

    double NormalizedTemperature =  pow( (rTemperature)/(rMeldTemperature), m );

    double Kv = K * exp((-1)*NormalizedTemperature);
    double nv = n * exp((-1)*NormalizedTemperature);


    double DeltaNormalizedTemperature = (m / rMeldTemperature) * pow( (rTemperature/rMeldTemperature), m-1 );

    rDeltaThermalHardening  = (1 + nv *std::log( rEquivalentPlasticStrain ) );

    rDeltaThermalHardening *= ( Kv * pow ( rEquivalentPlasticStrain, nv) * DeltaNormalizedTemperature );

    if( rRateFactor != 0 ){

      if( rDeltaGamma == 0 )
	std::cout<<" Something is wrong in the Baker_Johnson_Cook_hardening variables supplied "<<std::endl;

      rDeltaThermalHardening *= ( 1.0 + rRateFactor * C * std::log( (rDeltaGamma * sqrt(2.0/3.0))/(rPlasticStrainRate * rDeltaTime) ) );

    }

    return rDeltaThermalHardening;


    KRATOS_CATCH(" ")
  }




}  // namespace Kratos.
