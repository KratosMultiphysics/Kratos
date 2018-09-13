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
#include "custom_models/plasticity_models/hardening_rules/johnson_cook_thermal_hardening_rule.hpp"

namespace Kratos
{

  //*******************************CONSTRUCTOR******************************************
  //************************************************************************************

  JohnsonCookThermalHardeningRule::JohnsonCookThermalHardeningRule()
    :HardeningRule()
  {

  }


  //*******************************ASSIGMENT OPERATOR***********************************
  //************************************************************************************

  JohnsonCookThermalHardeningRule& JohnsonCookThermalHardeningRule::operator=(JohnsonCookThermalHardeningRule const& rOther)
  {
    HardeningRule::operator=(rOther);
    return *this;
  }

  //*******************************COPY CONSTRUCTOR*************************************
  //************************************************************************************

  JohnsonCookThermalHardeningRule::JohnsonCookThermalHardeningRule(JohnsonCookThermalHardeningRule const& rOther)
    :HardeningRule(rOther)
  {

  }


  //********************************CLONE***********************************************
  //************************************************************************************

  HardeningRule::Pointer JohnsonCookThermalHardeningRule::Clone() const
  {
    return Kratos::make_shared<JohnsonCookThermalHardeningRule>(*this);
  }


  //********************************DESTRUCTOR******************************************
  //************************************************************************************

  JohnsonCookThermalHardeningRule::~JohnsonCookThermalHardeningRule()
  {
  }

  /// Operations.

  //*******************************CALCULATE TOTAL HARDENING****************************
  //************************************************************************************

  double& JohnsonCookThermalHardeningRule::CalculateHardening(const PlasticDataType& rVariables, double& rHardening)
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
    const Properties& rMaterialProperties  = rModelData.GetMaterialProperties();
    const double& A = rMaterialProperties[JC_PARAMETER_A];
    const double& B = rMaterialProperties[JC_PARAMETER_B];
    const double& C = rMaterialProperties[JC_PARAMETER_C];

    const double& n = rMaterialProperties[JC_PARAMETER_n];
    const double& m = rMaterialProperties[JC_PARAMETER_m];

    const double& rReferenceTemperature = rMaterialProperties[REFERENCE_TEMPERATURE];
    const double& rMeldTemperature      = rMaterialProperties[MELD_TEMPERATURE];
    const double& rPlasticStrainRate    = rMaterialProperties[PLASTIC_STRAIN_RATE];

    double DeltaTemperature = rTemperature - rReferenceTemperature;
    if( DeltaTemperature < 0){
      if( DeltaTemperature < -1.0 )
	std::cout<<" Initial Temperature conditions not defined properly ("<<rTemperature<<" < "<<rReferenceTemperature<<") :"<<(rTemperature - rReferenceTemperature)<<std::endl;
      DeltaTemperature = 0;
    }

    double NormalizedTemperature = (1.0 - pow( (DeltaTemperature/(rMeldTemperature - rReferenceTemperature)), m) );

    // if( NormalizedTemperature < 0 )
    //   NormalizedTemperature = 0;

    if( rEquivalentPlasticStrain <= 0 )
      rHardening   =  A * NormalizedTemperature;
    else
      rHardening   = ( A + B * pow(rEquivalentPlasticStrain, n) ) * NormalizedTemperature;

    if( rRateFactor != 0 ){

      if( rDeltaGamma == 0 )
	std::cout<<" H Something is wrong in the Johnson_Cook_hardening variables supplied :: DeltaGamma= "<<rDeltaGamma<<" RateFactor= "<<rRateFactor<<std::endl;

      double RateComparisson = (rDeltaGamma * sqrt(2.0/3.0))/(rPlasticStrainRate * rDeltaTime);

      if( RateComparisson <= 0 )
	RateComparisson = 1e-40;

      rHardening  *= ( 1 + rRateFactor * C * std::log( RateComparisson ) );

    }


    return rHardening;


    KRATOS_CATCH(" ")
  }


  //*******************************CALCULATE HARDENING DERIVATIVE***********************
  //************************************************************************************

  double& JohnsonCookThermalHardeningRule::CalculateDeltaHardening(const PlasticDataType& rVariables, double& rDeltaHardening)
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
    const Properties& rMaterialProperties  = rModelData.GetMaterialProperties();
    const double& A = rMaterialProperties[JC_PARAMETER_A];
    const double& B = rMaterialProperties[JC_PARAMETER_B];
    const double& C = rMaterialProperties[JC_PARAMETER_C];

    const double& n = rMaterialProperties[JC_PARAMETER_n];
    const double& m = rMaterialProperties[JC_PARAMETER_m];

    const double& rReferenceTemperature = rMaterialProperties[REFERENCE_TEMPERATURE];
    const double& rMeldTemperature      = rMaterialProperties[MELD_TEMPERATURE];
    const double& rPlasticStrainRate    = rMaterialProperties[PLASTIC_STRAIN_RATE];


    double DeltaTemperature = rTemperature - rReferenceTemperature;
    if( DeltaTemperature < 0){
      if( DeltaTemperature < -1.0 )
	std::cout<<" Initial Temperature conditions not defined properly ("<<rTemperature<<" < "<<rReferenceTemperature<<") :"<<(rTemperature - rReferenceTemperature)<<std::endl;
      DeltaTemperature = 0;
    }

    double NormalizedTemperature = (1.0 - pow( (DeltaTemperature/(rMeldTemperature - rReferenceTemperature)), m) );

    // if( NormalizedTemperature < 0 )
    //   NormalizedTemperature = 0;

    if( rEquivalentPlasticStrain <= 0 )
      rDeltaHardening = 0;
    else
      rDeltaHardening = ( B * n * pow( rEquivalentPlasticStrain, n-1 ) ) * NormalizedTemperature;


    if( rRateFactor != 0 ){

      if( rDeltaGamma == 0 )
	std::cout<<" DH Something is wrong in the Johnson_Cook_hardening variables supplied :: DeltaGamma= "<<rDeltaGamma<<" RateFactor= "<<rRateFactor<<std::endl;

      double RateComparisson = (rDeltaGamma * sqrt(2.0/3.0))/(rPlasticStrainRate * rDeltaTime);

      if( RateComparisson <= 0 )
	RateComparisson = 1e-40;

      rDeltaHardening *= ( 1.0 + rRateFactor * C * std::log( RateComparisson ) );

      if( rEquivalentPlasticStrain <= 0 )
	rDeltaHardening += rRateFactor * ( sqrt(3.0/2.0) * ( A ) * NormalizedTemperature * C / rDeltaGamma );
      else
	rDeltaHardening += rRateFactor * ( sqrt(3.0/2.0) * ( A + B * pow( rEquivalentPlasticStrain, n ) ) * NormalizedTemperature * C / rDeltaGamma );


    }

    return rDeltaHardening;


    KRATOS_CATCH(" ")
  }


  //***************************CALCULATE HARDENING DERIVATIVE TEMPERATURE***************
  //************************************************************************************

  double& JohnsonCookThermalHardeningRule::CalculateDeltaThermalHardening(const PlasticDataType& rVariables, double& rDeltaThermalHardening)
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
    const Properties& rMaterialProperties  = rModelData.GetMaterialProperties();
    const double& A = rMaterialProperties[JC_PARAMETER_A];
    const double& B = rMaterialProperties[JC_PARAMETER_B];
    const double& C = rMaterialProperties[JC_PARAMETER_C];

    const double& n = rMaterialProperties[JC_PARAMETER_n];
    const double& m = rMaterialProperties[JC_PARAMETER_m];

    const double& rReferenceTemperature = rMaterialProperties[REFERENCE_TEMPERATURE];
    const double& rMeldTemperature      = rMaterialProperties[MELD_TEMPERATURE];
    const double& rPlasticStrainRate    = rMaterialProperties[PLASTIC_STRAIN_RATE];

    double DeltaTemperature = rTemperature - rReferenceTemperature;
    double DeltaNormalizedTemperature = 0;

    if( DeltaTemperature <= 0 ){
      if( DeltaTemperature < -1.0 )
	std::cout<<" Initial Temperature conditions not defined properly ("<<rTemperature<<" < "<<rReferenceTemperature<<") :"<<(rTemperature - rReferenceTemperature)<<std::endl;
      DeltaTemperature = 0;
    }
    else{
      DeltaNormalizedTemperature = ( pow( (DeltaTemperature/(rMeldTemperature - rReferenceTemperature)), m-1)/(rMeldTemperature - rReferenceTemperature) );
    }


    if( rEquivalentPlasticStrain < 0 )
      rDeltaThermalHardening  =  m * ( A ) * DeltaNormalizedTemperature;
    else
      rDeltaThermalHardening  =  m * ( A + B * pow ( rEquivalentPlasticStrain, n) ) * DeltaNormalizedTemperature;


    if( rRateFactor != 0 ){

      if( rDeltaGamma == 0 )
	std::cout<<" DTH Something is wrong in the Johnson_Cook_hardening variables supplied "<<std::endl;

      double RateComparisson = (rDeltaGamma * sqrt(2.0/3.0))/(rPlasticStrainRate * rDeltaTime);

      if( RateComparisson <= 0 )
	RateComparisson =  1e-40;

      rDeltaThermalHardening *=  ( 1.0 + rRateFactor * C * std::log( RateComparisson ) );


    }

    return rDeltaThermalHardening;


    KRATOS_CATCH(" ")
  }




}  // namespace Kratos.
