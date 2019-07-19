//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                    LHauser $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_models/plasticity_models/hardening_rules/casm_hardening_rule.hpp"
#include "custom_utilities/constitutive_model_utilities.hpp"

namespace Kratos
{

  //*******************************CONSTRUCTOR******************************************
  //************************************************************************************

  CasmHardeningRule::CasmHardeningRule()
    :HardeningRule()
  {
  }


  //*******************************ASSIGMENT OPERATOR***********************************
  //************************************************************************************

  CasmHardeningRule& CasmHardeningRule::operator=(CasmHardeningRule const& rOther)
  {
    HardeningRule::operator=(rOther);
    return *this;
  }

  //*******************************COPY CONSTRUCTOR*************************************
  //************************************************************************************

  CasmHardeningRule::CasmHardeningRule(CasmHardeningRule const& rOther)
    :HardeningRule(rOther)
  {

  }


  //********************************CLONE***********************************************
  //************************************************************************************

  HardeningRule::Pointer CasmHardeningRule::Clone() const
  {
    return ( HardeningRule::Pointer(new CasmHardeningRule(*this)) );
  }


  //********************************DESTRUCTOR******************************************
  //************************************************************************************

  CasmHardeningRule::~CasmHardeningRule()
  {
  }

  /// Operations.

  //*******************************CALCULATE TOTAL HARDENING****************************
  //************************************************************************************

  double& CasmHardeningRule::CalculateHardening(const PlasticDataType& rVariables, double& rHardening)
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
  double& CasmHardeningRule::CalculateDeltaHardening(const PlasticDataType& rVariables, double& rDeltaHardening)
  {
     KRATOS_ERROR << " not implemented. should not be here " << std::endl;
  }



  //*******************************CALCULATE HARDENING DERIVATIVE***********************
  //************************************************************************************
  double& CasmHardeningRule::CalculateDeltaHardening(const PlasticDataType& rVariables, double& rDeltaHardening, const MatrixType & rPlasticPotentialDerivative)
  {
     KRATOS_TRY

     const ModelDataType & rModelData = rVariables.GetModelData();
     const Properties& rMaterialProperties = rModelData.GetProperties();
     // const MatrixType    & rStressMatrix = rModelData.GetStressMatrix();

     const double & rPreconsolidationStress = rVariables.Internal.Variables[5];

     const double& rSpacingR = rMaterialProperties[SPACING_RATIO];
     const double& rOtherSlope   = rMaterialProperties[NORMAL_COMPRESSION_SLOPE];
     const double& rSwellingSlope   = rMaterialProperties[SWELLING_SLOPE];



     Vector PreconDerivativeEpsVol = ZeroVector(6);
     for (unsigned int i = 0; i<3; ++i)
        PreconDerivativeEpsVol(i) = 1.0;
     PreconDerivativeEpsVol *= -rPreconsolidationStress/(rOtherSlope - rSwellingSlope);



     rDeltaHardening = 1.0/( rPreconsolidationStress * std::log( rSpacingR));
     Vector PlasticPot(6);
     PlasticPot = ConstitutiveModelUtilities::StrainTensorToVector( rPlasticPotentialDerivative, PlasticPot);
     rDeltaHardening *= MathUtils<double>::Dot( PreconDerivativeEpsVol, PlasticPot);


     return rDeltaHardening;	

     KRATOS_CATCH(" ")
  }


}  // namespace Kratos.
