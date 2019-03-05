//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                    LHauser $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                February 2019 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_models/plasticity_models/hardening_rules/casm_cemented_hardening_rule.hpp"
#include "custom_utilities/constitutive_model_utilities.hpp"
#include "custom_utilities/stress_invariants_utilities.hpp"
#include "custom_utilities/shape_deviatoric_plane_mcc_utilities.hpp"

namespace Kratos
{

  //*******************************CONSTRUCTOR******************************************
  //************************************************************************************

  CasmCementedHardeningRule::CasmCementedHardeningRule()
    :HardeningRule()
  {
  }


  //*******************************ASSIGMENT OPERATOR***********************************
  //************************************************************************************

  CasmCementedHardeningRule& CasmCementedHardeningRule::operator=(CasmCementedHardeningRule const& rOther)
  {
    HardeningRule::operator=(rOther);
    return *this;
  }

  //*******************************COPY CONSTRUCTOR*************************************
  //************************************************************************************

  CasmCementedHardeningRule::CasmCementedHardeningRule(CasmCementedHardeningRule const& rOther)
    :HardeningRule(rOther)
  {

  }


  //********************************CLONE***********************************************
  //************************************************************************************

  HardeningRule::Pointer CasmCementedHardeningRule::Clone() const
  {
    return ( HardeningRule::Pointer(new CasmCementedHardeningRule(*this)) );
  }


  //********************************DESTRUCTOR******************************************
  //************************************************************************************

  CasmCementedHardeningRule::~CasmCementedHardeningRule()
  {
  }

  /// Operations.

  //*******************************CALCULATE TOTAL HARDENING****************************
  //************************************************************************************

  double& CasmCementedHardeningRule::CalculateHardening(const PlasticDataType& rVariables, double& rHardening)
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
  double& CasmCementedHardeningRule::CalculateDeltaHardening(const PlasticDataType& rVariables, double& rDeltaHardening)
  {
     KRATOS_ERROR << " not implemented. should not be here " << std::endl;
  }



  //*******************************CALCULATE HARDENING DERIVATIVE***********************
  //************************************************************************************
  double& CasmCementedHardeningRule::CalculateDeltaHardening(const PlasticDataType& rVariables, double& rDeltaHardening, const MatrixType & rPlasticPotentialDerivative)
  {
    KRATOS_TRY

    const ModelDataType& rModelData         = rVariables.GetModelData();
    const Properties& rMaterialProperties   = rModelData.GetProperties();
    const MatrixType& rStressMatrix         = rModelData.GetStressMatrix();

    //get constants
    const double& rShearM         = rMaterialProperties[CRITICAL_STATE_LINE];
    const double& rFriction       = rMaterialProperties[FRICTION_ANGLE];
    const double& rSpacingR       = rMaterialProperties[SPACING_RATIO];
    const double& rShapeN         = rMaterialProperties[SHAPE_PARAMETER];
    const double& rOtherSlope     = rMaterialProperties[NORMAL_COMPRESSION_SLOPE];
    const double& rSwellingSlope  = rMaterialProperties[SWELLING_SLOPE];
    const double& rOtherSlope     = rMaterialProperties[NORMAL_COMPRESSION_SLOPE];
    const double& rH1             = rMaterialProperties[DEGRADATION_RATE_COMPRESSION];
    const double& rH2             = rMaterialProperties[DEGRADATION_RATE_SHEAR];
    const double& rOmega          = rMaterialProperties[PLASTIC_DEVIATORIC_STRAIN_HARDENING];
    const double& rAlphaTensile   = rMaterialProperties[ALPHA_TENSILE];

    //get internal variables
    const double& rP0   = rVariables.Internal.Variables[4];
    const double& rB    = rVariables.Internal.Variables[5];
    const double& rPc   = rVariables.Internal.Variables[6];
    const double& rPt   = rVariables.Internal.Variables[7];

    //calculate stress invariants
    double MeanStress, LodeAngle, J2;
    StressInvariantsUtilities::CalculateStressInvariants( rStressMatrix, MeanStress, J2, LodeAngle);

    //calcualte third invariant effect
    double ThirdInvEffect = 1.0;
    ShapeAtDeviatoricPlaneMCCUtility::EvaluateEffectDueToThirdInvariant( ThirdInvEffect, LodeAngle, rFriction);

  //TODO
    //process plastic potential derivatives d_g/d_inv
    const double& FirstDerivativeP  = rPlasticPotentialDerivative(0,0);
    const double& FirstDerivativeJ2 = rPlasticPotentialDerivative(1,1);

    //calculate d_b/d_h * d_h/d_gamma
    double dBdGamma = -rB*( rH1*std::fabs( FirstDerivativeP ) + rH2*std::fabs( FirstDerivativeJ2 ) );

    //calculate d_P0/d_gamma
    double dP0dGamma = ( -FirstDerivativeP + rOmega*FirstDerivativeJ2) * rP0/(rOtherSlope-rSwellingSlope);

    //calcualte hardening modulus H
    rDeltaHardening = (rShapeN*std::pow(J2*std::sqrt(3.0),rShapeN))/( std::pow(rShearM/ThirdInvEffect,rShapeN)*std::pow(-(rPc-rPt),rShapeN+1.0) );
    rDeltaHardening += (rPc-MeanStress)/(std::log(rSpacingR)*(MeanStress+rPt)*(rPc+rPt));
    rDeltaHardening *= -( rAlphaTensile*rB*dP0dGamma + rP0*rAlphaTensile*dBdGamma );
    rDeltaHardening += 1.0/((rPc-rPt)*std::log(rSpacingR)) * ( (1.0+rB)*dP0dGamma + rP0*dBdGamma);
    
    return rDeltaHardening;	

    KRATOS_CATCH(" ")
  }

}  // namespace Kratos.
