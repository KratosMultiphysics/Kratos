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
#include "custom_constitutive/custom_yield_criteria/matsuoka_nakai_yield_criterion.hpp"

#include "pfem_solid_mechanics_application_variables.h"

namespace Kratos
{

//*******************************CONSTRUCTOR******************************************
//************************************************************************************

MatsuokaYieldCriterion::MatsuokaYieldCriterion()
	:YieldCriterion()
{
   
}

//*****************************INITIALIZATION CONSTRUCTOR*****************************
//************************************************************************************

MatsuokaYieldCriterion::MatsuokaYieldCriterion(HardeningLawPointer pHardeningLaw)
	:YieldCriterion(pHardeningLaw)
{
   
}

//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

MatsuokaYieldCriterion& MatsuokaYieldCriterion::operator=(MatsuokaYieldCriterion const& rOther)
{
   YieldCriterion::operator=(rOther);
   return *this;
}

//*******************************COPY CONSTRUCTOR*************************************
//************************************************************************************

MatsuokaYieldCriterion::MatsuokaYieldCriterion(MatsuokaYieldCriterion const& rOther)
	:YieldCriterion(rOther)
{

}


//********************************DESTRUCTOR******************************************
//************************************************************************************

MatsuokaYieldCriterion::~MatsuokaYieldCriterion()
{
}

/// Operations.


//***************************CALCULATE YIELD CONDITION********************************
//************************************************************************************

double& MatsuokaYieldCriterion::CalculateYieldCondition(double & rStateFunction, const Vector& rStressVector, const double& rAlpha)
{	
	double k = 12.0;
	double I1 = rStressVector(0)+rStressVector(1)+rStressVector(2);
	double I2 = rStressVector(0)*rStressVector(1);
	 I2 += rStressVector(1)*rStressVector(2);
	 I2 += rStressVector(2)*rStressVector(0);
	double I3 = rStressVector(0)*rStressVector(1)*rStressVector(2);

	rStateFunction = k*I3 -  I1*I2 ;
	return rStateFunction;
} 

void MatsuokaYieldCriterion::CalculateYieldFunctionDerivative(const Vector& rStressVector, Vector& rFirstDerivative)
{	
	double k = 12.0;
	double I1 = rStressVector(0)+rStressVector(1)+rStressVector(2);
	double I2 = rStressVector(0)*rStressVector(1);
	 I2 += rStressVector(1)*rStressVector(2);
	 I2 += rStressVector(2)*rStressVector(0);
	//double I3 = rStressVector(0)*rStressVector(1)*rStressVector(2);

	Vector I1D = ZeroVector(3);
	I1D(0) = 1.0;
	I1D(1) = 1.0;
	I1D(2) = 1.0;
	Vector I2D = ZeroVector(3);
	I2D(0) = rStressVector(1)+rStressVector(2);
	I2D(1) = rStressVector(2)+rStressVector(0);
	I2D(2) = rStressVector(0)+rStressVector(1);
	Vector I3D = ZeroVector(3);
	I3D(0) = rStressVector(1)*rStressVector(2);
	I3D(1) = rStressVector(2)*rStressVector(0);
	I3D(2) = rStressVector(0)*rStressVector(1);

	rFirstDerivative = k*I3D - I1D*I2 - I1*I2D ;

}







void MatsuokaYieldCriterion::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, YieldCriterion )
}

void MatsuokaYieldCriterion::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, YieldCriterion )
}


}  // namespace Kratos.
