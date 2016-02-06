// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/properties.h"

#include "custom_constitutive/custom_flow_rules/matsuoka_nakai_flow_rule.hpp"

#include "pfem_solid_mechanics_application_variables.h"

namespace Kratos
{



//************ CONSTRUCTOR ***********
MatsuokaNakaiFlowRule::MatsuokaNakaiFlowRule()
   :NonAssociativePlasticFlowRule()
{
}

//*****************************INITIALIZATION CONSTRUCTOR*****************************
//************************************************************************************

MatsuokaNakaiFlowRule::MatsuokaNakaiFlowRule(YieldCriterionPointer pYieldCriterion)
	:NonAssociativePlasticFlowRule(pYieldCriterion)
{
   
}

//********* ASSIGMENT OPERATOR
MatsuokaNakaiFlowRule& MatsuokaNakaiFlowRule::operator=(MatsuokaNakaiFlowRule const& rOther)
{
	NonAssociativePlasticFlowRule::operator=(rOther);
	return *this;

}



//********** COPY CONSTRUCTOR *********
MatsuokaNakaiFlowRule::MatsuokaNakaiFlowRule(MatsuokaNakaiFlowRule const& rOther)
      :NonAssociativePlasticFlowRule(rOther)
{
}




// ********** DESTRUCTOR **************
MatsuokaNakaiFlowRule::~MatsuokaNakaiFlowRule()
{
}


//************* COMPUTE DERIVATIVES OF THE PLASCIT POTENTIAL ********
void MatsuokaNakaiFlowRule::CalculatePlasticPotentialDerivatives(const Vector& rStressVector, Vector& rFirstDerivative, Matrix& rSecondDerivative)
{	
	double k = 11.0;
	// compute Invariants
	InvariantsStructure Inv;
	this->CalculateInvariantsAndDerivatives( rStressVector, Inv);

	// first Derivative
	rFirstDerivative  = k*Inv.I3D;
	rFirstDerivative  += -Inv.I2*Inv.I1D - Inv.I1*Inv.I2D;

	//SecondDerivative
	rSecondDerivative = Inv.I3DD;
	rSecondDerivative*= k;
	rSecondDerivative-= (MathUtils<double>::TensorProduct3(Inv.I1D, Inv.I2D));
	rSecondDerivative-= (MathUtils<double>::TensorProduct3(Inv.I2D, Inv.I1D));
	rSecondDerivative-= Inv.I1*Inv.I2DD; 

}


//************** COMPUTE INV AND DER ************************
void MatsuokaNakaiFlowRule::CalculateInvariantsAndDerivatives(const Vector & rStressVector, InvariantsStructure & Inv)
{
	Inv.I1 = rStressVector(0) + rStressVector(1) + rStressVector(2);
	
	Inv.I2 = rStressVector(1)*rStressVector(2) ;
	Inv.I2+= rStressVector(2)*rStressVector(0) ;
	Inv.I2+= rStressVector(0)*rStressVector(1) ;

	Inv.I3 = rStressVector(0)*rStressVector(1)*rStressVector(2);

	Inv.I1D = ZeroVector(3);
	Inv.I1D(0) = 1.0;
	Inv.I1D(1) = 1.0;
	Inv.I1D(2) = 1.0;
	Inv.I2D = ZeroVector(3);
	Inv.I2D(0) =  rStressVector(1)+rStressVector(2);
	Inv.I2D(1) =  rStressVector(2)+rStressVector(0);
	Inv.I2D(2) =  rStressVector(0)+rStressVector(1);


	Inv.I3D = ZeroVector(3);
	Inv.I3D(0) =  rStressVector(1)*rStressVector(2);
	Inv.I3D(1) =  rStressVector(2)*rStressVector(0);
	Inv.I3D(2) =  rStressVector(0)*rStressVector(1);

	Inv.I2DD = ZeroMatrix(3);
	for (int i = 0; i<3; ++i)
	{
		for (int j = 0; j<3; ++j)
		{
			if (i == j) {
				Inv.I2DD(i,j) = 0.0; }
			else {
				Inv.I2DD(i,j) = 1.0; }
		}
	}
		
	Inv.I3DD = ZeroMatrix(3);
	Inv.I3DD(0,0) = 0.0;
	Inv.I3DD(0,1) = rStressVector(2);
	Inv.I3DD(0,2) = rStressVector(1);
	Inv.I3DD(1,0) = rStressVector(2);
	Inv.I3DD(1,1) = 0.0;
	Inv.I3DD(1,2) = rStressVector(0);
	Inv.I3DD(2,0) = rStressVector(1);
	Inv.I3DD(2,1) = rStressVector(0);
	Inv.I3DD(2,2) = 0.0;


}




} 

