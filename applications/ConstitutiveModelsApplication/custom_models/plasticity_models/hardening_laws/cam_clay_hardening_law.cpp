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
#include "custom_models/plasticity_models/hardening_laws/cam_clay_hardening_law.hpp"

namespace Kratos
{

  //*******************************CONSTRUCTOR******************************************
  //************************************************************************************

  CamClayHardeningLaw::CamClayHardeningLaw()
    :HardeningLaw()
  {
  }


  //*******************************ASSIGMENT OPERATOR***********************************
  //************************************************************************************

  CamClayHardeningLaw& CamClayHardeningLaw::operator=(CamClayHardeningLaw const& rOther)
  {
    HardeningLaw::operator=(rOther);
    return *this;
  }

  //*******************************COPY CONSTRUCTOR*************************************
  //************************************************************************************

  CamClayHardeningLaw::CamClayHardeningLaw(CamClayHardeningLaw const& rOther)
    :HardeningLaw(rOther)
  {

  }


  //********************************CLONE***********************************************
  //************************************************************************************

  HardeningLaw::Pointer CamClayHardeningLaw::Clone() const
  {
    return ( HardeningLaw::Pointer(new CamClayHardeningLaw(*this)) );
  }


  //********************************DESTRUCTOR******************************************
  //************************************************************************************

  CamClayHardeningLaw::~CamClayHardeningLaw()
  {
  }

  /// Operations.

  //*******************************CALCULATE TOTAL HARDENING****************************
  //************************************************************************************

  double& CamClayHardeningLaw::CalculateHardening(const PlasticDataType& rVariables, double& rHardening)
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

  double& CamClayHardeningLaw::CalculateDeltaHardening(const PlasticDataType& rVariables, double& rDeltaHardening)
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
