//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra
//


// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/properties.h"
#include "custom_constitutive/custom_hardening_laws/cam_clay_hardening_law.hpp"
#include "includes/mat_variables.h"


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


//********************************DESTRUCTOR******************************************
//************************************************************************************

CamClayHardeningLaw::~CamClayHardeningLaw()
{
}

/// Operations.

//*******************************CALCULATE TOTAL HARDENING****************************
//************************************************************************************

/* This function return the updated Preconsolidation Pressure P_c with the following inputs:
    @ rAlpha: Hardening Parameter (Plastic Volumetric Strain)
    @ rOldPreconsolidationPressure: Old Value of Preconsolidation Pressure P_c at the previous time step t_n
*/
double& CamClayHardeningLaw::CalculateHardening(double &rHardening, const double &rAlpha, const double rOldPreconsolidationPressure)
{
    const double swelling_slope = GetProperties()[SWELLING_SLOPE];
    const double other_slope    = GetProperties()[NORMAL_COMPRESSION_SLOPE];

    rHardening = rOldPreconsolidationPressure * (std::exp (- rAlpha / (other_slope-swelling_slope) ) ) ;
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
