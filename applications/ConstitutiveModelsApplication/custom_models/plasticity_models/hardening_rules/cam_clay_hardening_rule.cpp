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
#include "custom_models/plasticity_models/hardening_rules/cam_clay_hardening_rule.hpp"

namespace Kratos
{

  //*******************************CONSTRUCTOR******************************************
  //************************************************************************************

  CamClayHardeningRule::CamClayHardeningRule()
    :HardeningRule()
  {
  }


  //*******************************ASSIGMENT OPERATOR***********************************
  //************************************************************************************

  CamClayHardeningRule& CamClayHardeningRule::operator=(CamClayHardeningRule const& rOther)
  {
    HardeningRule::operator=(rOther);
    return *this;
  }

  //*******************************COPY CONSTRUCTOR*************************************
  //************************************************************************************

  CamClayHardeningRule::CamClayHardeningRule(CamClayHardeningRule const& rOther)
    :HardeningRule(rOther)
  {

  }


  //********************************CLONE***********************************************
  //************************************************************************************

  HardeningRule::Pointer CamClayHardeningRule::Clone() const
  {
    return ( HardeningRule::Pointer(new CamClayHardeningRule(*this)) );
  }


  //********************************DESTRUCTOR******************************************
  //************************************************************************************

  CamClayHardeningRule::~CamClayHardeningRule()
  {
  }

  /// Operations.

  //*******************************CALCULATE TOTAL HARDENING****************************
  //************************************************************************************

  double& CamClayHardeningRule::CalculateHardening(const PlasticDataType& rVariables, double& rHardening)
  {
    KRATOS_TRY

    // get values
      const double & rVolumetricPlasticDeformation = rVariables.GetInternalVariables()[1];

    double FirstPreconsolidationPressure = 70.0;
    double SwellingSlope = 0.016;
    double OtherSlope = 0.1;


    rHardening = -FirstPreconsolidationPressure*(std::exp (-rVolumetricPlasticDeformation/(OtherSlope-SwellingSlope)) ) ;
    return rHardening;
	
    KRATOS_CATCH(" ")
  }
  


  //*******************************CALCULATE HARDENING DERIVATIVE***********************
  //************************************************************************************

  double& CamClayHardeningRule::CalculateDeltaHardening(const PlasticDataType& rVariables, double& rDeltaHardening)
  {
    KRATOS_TRY
      
    double SwellingSlope = 0.016;
    double OtherSlope = 0.1;

    const ModelDataType & rModelData = rVariables.GetModelData();
    const MatrixType    & rStressMatrix = rModelData.GetStressMatrix();

    double MeanStress = 0.0;
    for (unsigned int i = 0; i <3; i++)
       MeanStress += rStressMatrix(i,i)/3.0;

    double PreconsolidationStress = CalculateHardening(rVariables, PreconsolidationStress);


    rDeltaHardening = (2.0*MeanStress-PreconsolidationStress) ;
    rDeltaHardening *= (-MeanStress);
    rDeltaHardening *= PreconsolidationStress/ ( OtherSlope - SwellingSlope);
    return rDeltaHardening;	

    KRATOS_CATCH(" ")
  }


}  // namespace Kratos.
