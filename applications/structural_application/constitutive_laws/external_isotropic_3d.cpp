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
*   Date:                $Date: 2009-01-14 17:13:48 $
*   Revision:            $Revision: 1.4 $
*
* ***********************************************************/

// System includes 
#include <iostream>

// External includes 
#include<cmath>

// Project includes 

#include "includes/define.h"
#include "constitutive_laws/external_isotropic_3d.h"

#include "includes/constitutive_law.h"

#include "utilities/math_utils.h"
#include "includes/variables.h"
#include "includes/process_info.h"
// #include "structural_application.h"
#include "includes/properties.h"

namespace Kratos
{
    namespace ExternalExternalIsotropic3DAuxiliaries
    {
        boost::numeric::ublas::bounded_matrix<double,3,3> mstemp;
#ifdef _OPENMP
#pragma omp threadprivate(mstemp)
#endif
        boost::numeric::ublas::bounded_matrix<double,3,3> msaux;
#ifdef _OPENMP
#pragma omp threadprivate(msaux)
#endif
    }
    using namespace ExternalExternalIsotropic3DAuxiliaries;


	/**
	 *	TO BE TESTED!!!
	 */
	ExternalIsotropic3D::ExternalIsotropic3D() 
	: ConstitutiveLaw()
	{
	}
	/**
	 *	TO BE TESTED!!!
	 */
	ExternalIsotropic3D::~ExternalIsotropic3D()
	{
	}
	
	bool ExternalIsotropic3D::Has( const Variable<double>& rThisVariable )
	{
		return false;
	}
	
	bool ExternalIsotropic3D::Has( const Variable<Vector>& rThisVariable )
	{
		if( rThisVariable == INSITU_STRESS )
			return true;
		return false;
	}
	
	bool ExternalIsotropic3D::Has( const Variable<Matrix>& rThisVariable )
	{
		return false;
	}
	
	double& ExternalIsotropic3D::GetValue( const Variable<double>& rThisVariable, double& rValue )
	{
		return( rValue );
	}
	
	Vector& ExternalIsotropic3D::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
	{
		if( rThisVariable == INSITU_STRESS )
            return mInSituStress;
        else if( rThisVariable == MATERIAL_PARAMETERS )
            return mMaterialParameters;
        else
            return rValue;
	}
	
	Matrix& ExternalIsotropic3D::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
	{
		return( rValue );
	}

	void ExternalIsotropic3D::SetValue( const Variable<double>& rThisVariable, const double& rValue, 
								   const ProcessInfo& rCurrentProcessInfo )
	{
	}
	
	void ExternalIsotropic3D::SetValue( const Variable<array_1d<double, 3> >& rThisVariable, 
								   const array_1d<double,3>& rValue, 
								   const ProcessInfo& rCurrentProcessInfo )
	{
	}
	
	void ExternalIsotropic3D::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
								 const ProcessInfo& rCurrentProcessInfo )
	{
		if( rThisVariable == INSITU_STRESS )
		{
			mInSituStress = rValue;
		}
                if( rThisVariable == MATERIAL_PARAMETERS )
                {
                    mMaterialParameters = rValue;
                }
	}
	
	void ExternalIsotropic3D::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
								 const ProcessInfo& rCurrentProcessInfo )
	{
	}
	
	/**
	 *	TO BE TESTED!!!
	 */
	void ExternalIsotropic3D::InitializeMaterial( const Properties& props,
			const GeometryType& geom,
			const Vector& ShapeFunctionsValues )
	{
//             std::cout << "Initilizing Material" << std::endl;
             
            mCurrentStress = ZeroVector(6);
            mCtangent = ZeroMatrix(6,6);
            mInSituStress = ZeroVector(6);
            mMaterialParameters = props[MATERIAL_PARAMETERS];
//             std::cout << "properties: " << props << std::endl;
//             std::cout << "E: " << mMaterialParameters[0] << std::endl;
            CalculateElasticMatrix(mCtangent, mMaterialParameters[0], mMaterialParameters[1]);
	}
	
	void ExternalIsotropic3D::InitializeSolutionStep( const Properties& props,
				     const GeometryType& geom, //this is just to give the array of nodes
				     const Vector& ShapeFunctionsValues ,
				     const ProcessInfo& CurrentProcessInfo)
	{
	}
			
	void ExternalIsotropic3D::FinalizeSolutionStep( const Properties& props,
				   const GeometryType& geom, //this is just to give the array of nodes
				   const Vector& ShapeFunctionsValues ,
				   const ProcessInfo& CurrentProcessInfo)
	{
		if( CurrentProcessInfo[CALCULATE_INSITU_STRESS] )
		{
			mInSituStress -= mCurrentStress;
			//SetValue( INSITU_STRESS, mInSituStress, CurrentProcessInfo );
		}
	}
    
    void ExternalIsotropic3D::UpdateMaterial( const Vector& StrainVector,
                                      const Properties& props,
                                      const GeometryType& geom,
                                      const Vector& ShapeFunctionsValues,
                                      const ProcessInfo& CurrentProcessInfo )
    {
        CalculateElasticMatrix(mCtangent, mMaterialParameters[0], mMaterialParameters[1]);
    }
		
	/**
	 *	TO BE TESTED!!!
	 */
	void ExternalIsotropic3D::CalculateElasticMatrix(Matrix& C, const double E, const double NU)
	{ 
		//setting up material matrix
		double c1 = E / ((1.00+NU)*(1-2*NU));
		double c2 = c1 * (1-NU);
		double c3 = c1 * NU;
		double c4 = c1 * 0.5 * (1 - 2*NU);
		//filling material matrix
		C(0,0) = c2;    C(0,1) = c3;    C(0,2) = c3;    C(0,3) = 0.0;   C(0,4) = 0.0;   C(0,5) = 0.0;
		C(1,0) = c3;    C(1,1) = c2;    C(1,2) = c3;    C(1,3) = 0.0;   C(1,4) = 0.0;   C(1,5) = 0.0;
		C(2,0) = c3;    C(2,1) = c3;    C(2,2) = c2;    C(2,3) = 0.0;   C(2,4) = 0.0;   C(2,5) = 0.0;
		C(3,0) = 0.0;   C(3,1) = 0.0;   C(3,2) = 0.0;   C(3,3) = c4;    C(3,4) = 0.0;   C(3,5) = 0.0;
		C(4,0) = 0.0;   C(4,1) = 0.0;   C(4,2) = 0.0;   C(4,3) = 0.0;   C(4,4) = c4;    C(4,5) = 0.0;
		C(5,0) = 0.0;   C(5,1) = 0.0;   C(5,2) = 0.0;   C(5,3) = 0.0;   C(5,4) = 0.0;   C(5,5) = c4;
		
	}
	
	/**
	 *	TO BE TESTED!!!
	 */
	void ExternalIsotropic3D::CalculateStress(const Vector& StrainVector, Vector& StressVector)
	{
		if( StressVector.size() != 6 )
		{
			StressVector.resize(6);
		}
		noalias(StressVector) = prod(mCtangent,StrainVector);
		mCurrentStress = StressVector;
		noalias(StressVector) -= mInSituStress;
// 		double c1 = mE / ((1.00+mNU)*(1-2*mNU));
// 		double c2 = c1 * (1-mNU);
// 		double c3 = c1 * mNU;
// 		double c4 = c1 * 0.5 * (1 - 2*mNU);
// 
// 		StressVector[0] = c2*StrainVector[0] + c3 * (StrainVector[1] + StrainVector[2])	;
// 		StressVector[1] = c2*StrainVector[1] + c3 * (StrainVector[0] + StrainVector[2])	;
// 		StressVector[2] = c2*StrainVector[2] + c3 * (StrainVector[0] + StrainVector[1])	;
// 		StressVector[3] = c4*StrainVector[3];
// 		StressVector[4] = c4*StrainVector[4];
// 		StressVector[5] = c4*StrainVector[5];
	}
	
	/**
	 *	TO BE REVIEWED!!!
	 */
	void ExternalIsotropic3D::CalculateConstitutiveMatrix(const Vector& StrainVector, Matrix& rResult)
	{
		rResult = mCtangent;
	}
	
    //**********************************************************************
    void ExternalIsotropic3D::CalculateCauchyStresses(
		Vector& rCauchy_StressVector,
		const Matrix& rF,
		const Vector& rPK2_StressVector,
		const Vector& rGreenLagrangeStrainVector)
    {
		Matrix S = MathUtils<double>::StressVectorToTensor( rPK2_StressVector );

        double J = MathUtils<double>::Det3( rF );

		noalias(mstemp) = prod(rF,S);
		noalias(msaux) = prod(mstemp,trans(rF));
		msaux *= J;

		if(rCauchy_StressVector.size() != 6)
			rCauchy_StressVector.resize(6);
		
		rCauchy_StressVector[0] = msaux(0,0);
		rCauchy_StressVector[1] = msaux(1,1);
		rCauchy_StressVector[2] = msaux(2,2);
		rCauchy_StressVector[3] = msaux(1,2);
		rCauchy_StressVector[4] = msaux(1,3);
		rCauchy_StressVector[5] = msaux(2,3);
    }
} // Namespace Kratos
