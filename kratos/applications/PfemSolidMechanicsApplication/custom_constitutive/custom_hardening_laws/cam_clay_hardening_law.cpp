//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/properties.h"
#include "solid_mechanics_application.h"
#include "custom_constitutive/custom_hardening_laws/cam_clay_hardening_law.hpp"


namespace Kratos
{

//*******************************CONSTRUCTOR******************************************
//************************************************************************************

CamClayKinematicHardeningLaw::CamClayKinematicHardeningLaw()
	:HardeningLaw()
{
   //Combined isotropic-kinematic 0<mTheta<1
   //Pure isotropic hardening mTheta=1;  
   //Pure kinematic hardening mTheta=0; 

   //Hardening law:
   //mTheta = 1; 
   
}


//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

CamClayKinematicHardeningLaw& CamClayKinematicHardeningLaw::operator=(CamClayKinematicHardeningLaw const& rOther)
{
   HardeningLaw::operator=(rOther);
   return *this;
}

//*******************************COPY CONSTRUCTOR*************************************
//************************************************************************************

CamClayKinematicHardeningLaw::CamClayKinematicHardeningLaw(CamClayKinematicHardeningLaw const& rOther)
	:HardeningLaw(rOther)
{

}


//********************************DESTRUCTOR******************************************
//************************************************************************************

CamClayKinematicHardeningLaw::~CamClayKinematicHardeningLaw()
{
}

/// Operations.

//*******************************CALCULATE TOTAL HARDENING****************************
//************************************************************************************

double& CamClayKinematicHardeningLaw::CalculateHardening(double &rHardening, const double &rAlpha, const double rTemperature)
{
	
/*	//linear hardening properties
	const double& YieldStress                 =  GetProperties()[YIELD_STRESS];
	const double& KinematicHardeningConstant  =  GetProperties()[KINEMATIC_HARDENING_MODULUS];
	

	//Linear Hardening law:
	rHardening  = YieldStress + mTheta *  KinematicHardeningConstant * rAlpha;
*/

//std::cout << " Hardening ! Estoy " << std::endl;
    double ReferencePreasure = -80.0;
    double SwellingSlope = 0.0078;
    double OtherSlope = 0.085;
    double AlphaShear = 120.0;
    double Beta = 1.0;
    double VoidRatio = 0.6;


        rHardening = ReferencePreasure*(std::exp (-rAlpha/(OtherSlope-SwellingSlope)) ) ;
	return rHardening;

}
  


void CamClayKinematicHardeningLaw::save( Serializer& rSerializer ) const
{
    //KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, NonLinearIsotropicKinematicHardeningLaw );

}

void CamClayKinematicHardeningLaw::load( Serializer& rSerializer )
{
    //KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, NonLinearIsotropicKinematicHardeningLaw );

}


}  // namespace Kratos.
