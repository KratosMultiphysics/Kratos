/*
==============================================================================
KratosStructuralApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel 
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/
/* *********************************************************   
*          
*   Last Modified by:    $Author: janosch $
*   Date:                $Date: 2009-01-14 17:14:12 $
*   Revision:            $Revision: 1.3 $
*
* ***********************************************************/


// System includes 
#include <iostream>

// External includes 
#include<cmath>

// Project includes 

#include "includes/define.h"
#include "constitutive_laws/hypoelastic_2d.h"
#include "custom_utilities/sd_math_utils.h"
#include "includes/constitutive_law.h"

#include "utilities/math_utils.h"
#include "includes/variables.h"
#include "includes/process_info.h"
#include "structural_application.h"
#include "includes/properties.h"

namespace Kratos
{


	/**
	 *	TO BE TESTED!!!
	 */
	Hypoelastic2D::Hypoelastic2D() 
	: ConstitutiveLaw()
	{
	}
	/**
	 *	TO BE TESTED!!!
	 */
	Hypoelastic2D::~Hypoelastic2D()
	{
	}
	
	
	bool Hypoelastic2D::Has( const Variable<double>& rThisVariable )
	{
		return false;
	}
	
	bool Hypoelastic2D::Has( const Variable<Vector>& rThisVariable )
	{
		return false;
	}
	
	bool Hypoelastic2D::Has( const Variable<Matrix>& rThisVariable )
	{
		return false;
	}
	
	double& Hypoelastic2D::GetValue( const Variable<double>& rThisVariable, double& rValue )
    {
        return( rValue );
    }
    
    Vector& Hypoelastic2D::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
    {
        return( rValue );
    }
    
    Matrix& Hypoelastic2D::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
    {
        return( rValue );
    }

	void Hypoelastic2D::SetValue( const Variable<double>& rThisVariable, const double& rValue, const ProcessInfo& rCurrentProcessInfo )
	{
	}
	
	void Hypoelastic2D::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue, const ProcessInfo& rCurrentProcessInfo )
	{
	}
	
	void Hypoelastic2D::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue, const ProcessInfo& rCurrentProcessInfo )
	{
	}
	
	/**
	 *	TO BE TESTED!!!
	 */
	void Hypoelastic2D::InitializeMaterial( const Properties& props,
			const GeometryType& geom,
			const Vector& ShapeFunctionsValues )
	{
		mMIU = props[MIU];
		mLAMBDA = props[LAMBDA];
//	KRATOS_WATCH("INSIDE THE CONS LAW")
//		KRATOS_WATCH(mE)
//			KRATOS_WATCH(mNU)
	}

	void Hypoelastic2D::InitializeSolutionStep( const Properties& props,
					const GeometryType& geom, //this is just to give the array of nodes
					const Vector& ShapeFunctionsValues ,
					const ProcessInfo& CurrentProcessInfo)
	{
		mDeltaT = props[DELTA_TIME];
	}
		
	/**
	 *	TO BE TESTED!!!
	 */
	void Hypoelastic2D::CalculateHypoElasticMatrix(Matrix& C, const double base_miu, const double base_lambda, const double deltaT, const double detF)
	{ 

	double miu = base_miu*deltaT/detF;
	double lambda= base_lambda*deltaT/detF;

	double Ehypo = miu*(3*lambda+2*miu)/(lambda+miu);
	double NUhypo = lambda/(2*(lambda+miu));
	double c1 = Ehypo / (1.00 - NUhypo*NUhypo);
	double c2 = c1 * NUhypo;
	double c3 = 0.5*Ehypo / (1 + NUhypo);


		C(0,0) = c1;	C(0,1) = c2;	C(0,2) = 0.0;
		C(1,0) = c2;	C(1,1) = c1;	C(1,2) = 0.0;
		C(2,0) = 0.0;	C(2,1) = 0.0;	C(2,2) = c3;
	}
	
	/**
	 *	TO BE TESTED!!!
	 */
// 	void Hypoelastic2D::CalculateStress(const Vector& StrainVector, Vector& StressVector)
// 	{
// 		double c1 = mE / (1.00 - mNU*mNU);
// 		double c2 = c1 * mNU;
// 		double c3 = 0.5*mE / (1 + mNU);
// 
// 		StressVector[0] = c1*StrainVector[0] + c2 * (StrainVector[1])	;
// 		StressVector[1] = c1*StrainVector[1] + c2 * (StrainVector[0])	;
// 		StressVector[2] = c3*StrainVector[2];
// 
// 	}

	void Hypoelastic2D::CalculatePlaneHypoElasticMatrix(Matrix& C, const double base_miu, const double base_lambda, const double deltaT, const double detF)
	{ 

	double miu = base_miu*deltaT/detF;
	double lambda= base_lambda*deltaT/detF;

// 	double Ehypo = miu*(3*lambda+2*miu)/(lambda+miu);
// 	double NUhypo = lambda/(2*(lambda+miu));
// 	double c1 = Ehypo*(1.00-NUhypo) / ((1.00 + NUhypo)*(1.00 - 2*NUhypo));
// 	double c2 = c1 * NUhypo/(1.00-NUhypo);
// 	double c3 = 0.5*Ehypo / (1 + NUhypo);

	double c1 = 2*miu+lambda;
	double c2 = lambda;
	double c3 = miu;
	
	C(0,0) = c1;	C(0,1) = c2;	C(0,2) = 0.0;
	C(1,0) = c2;	C(1,1) = c1;	C(1,2) = 0.0;
	C(2,0) = 0.0;	C(2,1) = 0.0;	C(2,2) = c3;


	//KRATOS_WATCH("inside D");
	//KRATOS_WATCH(Ehypo);
	//KRATOS_WATCH(NUhypo);
	//KRATOS_WATCH(detF);
	}
	/**
	 *	TO BE REVIEWED!!!
	 */
	void Hypoelastic2D::CalculateStressAndTangentMatrix( Matrix& StressTensor, const Matrix& F,const Matrix& StrainTensor, Matrix& algorithmicTangent)
	{
	double detF = MathUtils<double>::Det2( F );
	CalculatePlaneHypoElasticMatrix( algorithmicTangent, mMIU, mLAMBDA, mDeltaT, detF );


	}




// 	void Hypoelastic2D::CalculateConstitutiveMatrix(const Vector& StrainVector, Matrix& rResult)
// 	{
// 		CalculateElasticMatrix( rResult, mE, mNU );
// 	}

// 	void Hypoelastic2D::PlaneStrainConstitutiveMatrix(const Vector& StrainVector, Matrix& rResult)
// 	{
// 		CalculatePlaneElasticMatrix( rResult, mE, mNU );
// 	}
	
	
	
    //**********************************************************************
	    void Hypoelastic2D::CalculateCauchyStresses(
		Vector& rCauchy_StressVector,
		const Matrix& rF,
		const Vector& rPK2_StressVector,
		const Vector& rGreenLagrangeStrainVector)
    {
		//KRATOS_WATCH("INSIDE SOLID CalculateCauchyStresses");
		Matrix S = MathUtils<double>::StressVectorToTensor( rPK2_StressVector );

		double J = MathUtils<double>::Det2( rF );
		boost::numeric::ublas::bounded_matrix<double,2,2> mstemp;
		boost::numeric::ublas::bounded_matrix<double,2,2> msaux;

		noalias(mstemp) = prod(rF,S);
		noalias(msaux) = prod(mstemp,trans(rF));
		msaux /= J;

		if(rCauchy_StressVector.size() != 3)
			rCauchy_StressVector.resize(3);
		
		rCauchy_StressVector[0] = msaux(0,0);
		rCauchy_StressVector[1] = msaux(1,1);
		rCauchy_StressVector[2] = msaux(0,1);
    }
    
          int Hypoelastic2D::Check(const Properties& props, const GeometryType& geom, const ProcessInfo& CurrentProcessInfo)
       {
	  
            if(DENSITY.Key() == 0 || props[DENSITY]<0.00)
                KRATOS_ERROR(std::invalid_argument,"DENSITY has Key zero or invalid value ","");
	 
	    if(MIU.Key() == 0 || props[MIU]<0.00)
                KRATOS_ERROR(std::invalid_argument,"MIU has Key zero or invalid value ","");
	    
	    if(LAMBDA.Key() == 0 || props[LAMBDA]< 0.00) 
               KRATOS_ERROR(std::invalid_argument,"LAMBDA has Key zero or invalid value ","");
	    
	    return 0;
         }
    
} // Namespace Kratos
