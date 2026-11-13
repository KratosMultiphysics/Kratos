//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2015 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/properties.h"
#include "custom_constitutive/custom_hardening_laws/cam_clay_hardening_law.hpp"
#include "pfem_solid_mechanics_application_variables.h"


namespace Kratos
{

//*******************************CONSTRUCTOR******************************************
//************************************************************************************

CamClayHardeningLaw::CamClayHardeningLaw()
	:HardeningLaw()
{
   
}


//******************************ASSIGNMENT OPERATOR***********************************
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


//********************************DESTRUCTOR******************************************
//************************************************************************************

CamClayHardeningLaw::~CamClayHardeningLaw()
{
}

/// Operations.

//*******************************CALCULATE TOTAL HARDENING****************************
//************************************************************************************

double& CamClayHardeningLaw::CalculateHardening(double &rHardening, const double &rAlpha, const double rTemperature)
{
	

    double FirstPreconsolidationPressure = GetProperties()[PRE_CONSOLIDATION_STRESS];
    double SwellingSlope = GetProperties()[SWELLING_SLOPE];
    double OtherSlope = GetProperties()[NORMAL_COMPRESSION_SLOPE];


    rHardening = -FirstPreconsolidationPressure*(std::exp (-rAlpha/(OtherSlope-SwellingSlope)) ) ;
    return rHardening;

}
  


void CamClayHardeningLaw::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, HardeningLaw )

}

void CamClayHardeningLaw::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, HardeningLaw )

}


}  // namespace Kratos.
