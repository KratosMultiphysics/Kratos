//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_models/plasticity_models/hardening_rules/gens_nova_hardening_rule.hpp"

namespace Kratos
{

  //*******************************CONSTRUCTOR******************************************
  //************************************************************************************

  GensNovaHardeningRule::GensNovaHardeningRule()
    :HardeningRule()
  {
  }


  //*******************************ASSIGMENT OPERATOR***********************************
  //************************************************************************************

  GensNovaHardeningRule& GensNovaHardeningRule::operator=(GensNovaHardeningRule const& rOther)
  {
    HardeningRule::operator=(rOther);
    return *this;
  }

  //*******************************COPY CONSTRUCTOR*************************************
  //************************************************************************************

  GensNovaHardeningRule::GensNovaHardeningRule(GensNovaHardeningRule const& rOther)
    :HardeningRule(rOther)
  {

  }


  //********************************CLONE***********************************************
  //************************************************************************************

  HardeningRule::Pointer GensNovaHardeningRule::Clone() const
  {
    return ( HardeningRule::Pointer(new GensNovaHardeningRule(*this)) );
  }


  //********************************DESTRUCTOR******************************************
  //************************************************************************************

  GensNovaHardeningRule::~GensNovaHardeningRule()
  {
  }

  /// Operations.

  //*******************************CALCULATE TOTAL HARDENING****************************
  //************************************************************************************

  double& GensNovaHardeningRule::CalculateHardening(const PlasticDataType& rVariables, double& rHardening)
  {
    KRATOS_TRY


    const ModelDataType & rModelData = rVariables.GetModelData();
    const Properties& rMaterialProperties = rModelData.GetProperties();

    // get values
    const double & rVolumetricPlasticDeformation = rVariables.GetInternalVariables()[1];

    // Set constitutive parameters
    const double & rFirstPreconsolidationPressure = rMaterialProperties[PRE_CONSOLIDATION_STRESS];
    const double & rSwellingSlope = rMaterialProperties[SWELLING_SLOPE];
    const double & rOtherSlope = rMaterialProperties[NORMAL_COMPRESSION_SLOPE];


    rHardening = -rFirstPreconsolidationPressure*(std::exp (-rVolumetricPlasticDeformation/(rOtherSlope-rSwellingSlope)) ) ;

    KRATOS_ERROR << " you should not be here " << std::endl;
    return rHardening;
	
    KRATOS_CATCH(" ")
  }
  


  //*******************************CALCULATE HARDENING DERIVATIVE***********************
  //************************************************************************************
  double& GensNovaHardeningRule::CalculateDeltaHardening(const PlasticDataType& rVariables, double& rDeltaHardening)
  {
     KRATOS_ERROR << " not implemented. should not be here " << std::endl;
  }



  //*******************************CALCULATE HARDENING DERIVATIVE***********************
  //************************************************************************************
  double& GensNovaHardeningRule::CalculateDeltaHardening(const PlasticDataType& rVariables, double& rDeltaHardening, const MatrixType & rPlasticPotentialDerivative)
  {
    KRATOS_TRY

      
    const ModelDataType & rModelData = rVariables.GetModelData();
    const Properties& rMaterialProperties = rModelData.GetProperties();
    const MatrixType    & rStressMatrix = rModelData.GetStressMatrix();

    const double & rSwellingSlope = rMaterialProperties[SWELLING_SLOPE];
    const double & rOtherSlope = rMaterialProperties[NORMAL_COMPRESSION_SLOPE];

    const double & rPS     = rVariables.Internal.Variables[3];
    const double & rPT     = rVariables.Internal.Variables[4];
    const double & rPCstar = rVariables.Internal.Variables[5];
    const double & rhos = rMaterialProperties[RHOS];
    const double & rhot = rMaterialProperties[RHOT];
    const double & k = rMaterialProperties[KSIM];

    // Perform the stress translation
    double MeanStressT = 0.0;
    for (unsigned int i = 0; i <3; i++)
       MeanStressT += rStressMatrix(i,i) + rPT;
    MeanStressT /= 3.0;

    double TracePlasticPotDerivative = 0.0;
    for (unsigned int i = 0; i < 3; i++)
       TracePlasticPotDerivative += rPlasticPotentialDerivative(i,i);
    

    MatrixType DevPlasticPotDerivative = rPlasticPotentialDerivative;
    for (unsigned int i = 0; i <3; i++)
       DevPlasticPotDerivative(i,i) -= TracePlasticPotDerivative/3.0;

    double NormPlasticPotDerivative = 0.0;
    for (unsigned int i = 0; i <3; i++) {
       for (unsigned int j = 0; j < 3; j++) {
          NormPlasticPotDerivative += pow( DevPlasticPotDerivative(i,j), 2.0);
       }
    }
    NormPlasticPotDerivative = sqrt( NormPlasticPotDerivative); 

    const double & chis = rMaterialProperties[CHIS];
    const double & chit = rMaterialProperties[CHIT];

    rDeltaHardening = (-MeanStressT) * rhos * (+rPS )* (TracePlasticPotDerivative + chis * sqrt(2.0/3.0) * NormPlasticPotDerivative );
 rDeltaHardening += (-MeanStressT) * ( 1.0 + k) * (rhot) * (-rPT) * ( fabs(TracePlasticPotDerivative) + chit*sqrt(2.0/3.0) * NormPlasticPotDerivative);


    return rDeltaHardening;	

    KRATOS_CATCH(" ")
  }


}  // namespace Kratos.
