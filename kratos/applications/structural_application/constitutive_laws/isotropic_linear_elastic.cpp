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
*   Last Modified by:    $Author: hurga $
*   Date:                $Date: 2009-03-05 12:01:22 $
*   Revision:            $Revision: 1.4 $
*
* ***********************************************************/

// System includes 
#include <iostream>

// External includes 
#include<cmath>

// Project includes 

#include "includes/define.h"
#include "structural_application/custom_constitutive_laws/isotropic_linear_elastic.h"

#include "includes/constitutive_law.h"

#include "utilities/math_utils.h"
#include "structural_application/custom_utilities/sd_math_utils.h"
#include "includes/variables.h"
#include "includes/process_info.h"
#include "structural_application/structural_application.h"
#include "includes/properties.h"

namespace Kratos
{
	/**
	 *	TO BE TESTED!!!
	 */
	IsotropicLinearElastic::IsotropicLinearElastic() 
	: ConstitutiveLaw()
	{
	}
	/**
	 *	TO BE TESTED!!!
	 */
	IsotropicLinearElastic::~IsotropicLinearElastic()
	{
	}
	
	bool IsotropicLinearElastic::Has( const Variable<double>& rThisVariable )
	{
		return false;
	}
	
	bool IsotropicLinearElastic::Has( const Variable<Vector>& rThisVariable )
	{
		if( rThisVariable == INSITU_STRESS )
			return true;
		return false;
	}
	
	bool IsotropicLinearElastic::Has( const Variable<Matrix>& rThisVariable )
	{
		return false;
	}
	
	double IsotropicLinearElastic::GetValue( const Variable<double>& rThisVariable )
	{
	}
	
	Vector IsotropicLinearElastic::GetValue( const Variable<Vector>& rThisVariable )
	{
		if( rThisVariable == INSITU_STRESS )
			return mInSituStress;
	}
	
	Matrix IsotropicLinearElastic::GetValue( const Variable<Matrix>& rThisVariable )
	{
	}

	void IsotropicLinearElastic::SetValue( const Variable<double>& rThisVariable, double rValue )
	{
	}
	
	void IsotropicLinearElastic::SetValue( const Variable<Vector>& rThisVariable, Vector rValue )
	{
		if( rThisVariable == INSITU_STRESS )
		{
			mInSituStress = rValue;
		}
	}
	
	void IsotropicLinearElastic::SetValue( const Variable<Matrix>& rThisVariable, Matrix rValue )
	{
	}

	/**
	 *	TO BE TESTED!!!
	 */
	void IsotropicLinearElastic::InitializeMaterial( const Properties& props,
			const GeometryType& geom,
			const Vector& ShapeFunctionsValues )
	{
// 		std::cout << "Initialising IsotropicLinearElastic" << std::endl;
		mCurrentStress = ZeroVector(6);
		SetValue(INSITU_STRESS, ZeroVector(6) );
		mE = props[YOUNG_MODULUS];
		mNU = props[POISSON_RATIO];
		mCtangent = ZeroMatrix(6,6);
		mInSituStress = ZeroVector(6);
		//setting up material matrix
		double c1 = mE / ((1.00+mNU)*(1-2*mNU));
		double c2 = c1 * (1-mNU);
		double c3 = c1 * mNU;
		double c4 = c1 * 0.5 * (1 - 2*mNU);
		//filling material matrix
		mCtangent(0,0) = c2;    mCtangent(0,1) = c3;    mCtangent(0,2) = c3;    mCtangent(0,3) = 0.0;   mCtangent(0,4) = 0.0;   mCtangent(0,5) = 0.0;
		mCtangent(1,0) = c3;    mCtangent(1,1) = c2;    mCtangent(1,2) = c3;    mCtangent(1,3) = 0.0;   mCtangent(1,4) = 0.0;   mCtangent(1,5) = 0.0;
		mCtangent(2,0) = c3;    mCtangent(2,1) = c3;    mCtangent(2,2) = c2;    mCtangent(2,3) = 0.0;   mCtangent(2,4) = 0.0;   mCtangent(2,5) = 0.0;
		mCtangent(3,0) = 0.0;   mCtangent(3,1) = 0.0;   mCtangent(3,2) = 0.0;   mCtangent(3,3) = c4;    mCtangent(3,4) = 0.0;   mCtangent(3,5) = 0.0;
		mCtangent(4,0) = 0.0;   mCtangent(4,1) = 0.0;   mCtangent(4,2) = 0.0;   mCtangent(4,3) = 0.0;   mCtangent(4,4) = c4;    mCtangent(4,5) = 0.0;
		mCtangent(5,0) = 0.0;   mCtangent(5,1) = 0.0;   mCtangent(5,2) = 0.0;   mCtangent(5,3) = 0.0;   mCtangent(5,4) = 0.0;   mCtangent(5,5) = c4;
	}
	
	void IsotropicLinearElastic::InitializeSolutionStep( const Properties& props,
				     const GeometryType& geom, //this is just to give the array of nodes
				     const Vector& ShapeFunctionsValues ,
				     const ProcessInfo& CurrentProcessInfo)
	{
	}
			
	void IsotropicLinearElastic::FinalizeSolutionStep( const Properties& props,
				   const GeometryType& geom, //this is just to give the array of nodes
				   const Vector& ShapeFunctionsValues ,
				   const ProcessInfo& CurrentProcessInfo)
	{
		if( CurrentProcessInfo[CALCULATE_INSITU_STRESS] )
		{
			mInSituStress -= mCurrentStress;
		}
	}
	
	
	/**
	 *	TO BE TESTED!!!
	 */
	void IsotropicLinearElastic::CalculateStress(const Vector& StrainVector, Vector& StressVector)
	{
		if( StressVector.size() != 6 )
		{
            StressVector.resize(6,false);
		}
		noalias(StressVector) = prod(mCtangent,StrainVector);
		mCurrentStress = StressVector;
		noalias(StressVector) -= mInSituStress;
	}
	
	void IsotropicLinearElastic::Calculate( const Variable<Vector >& rVariable, Vector& rResult )
	{
// 		if(rVariable==INSITU_STRESS)
// 		{
// 			for( unsigned int i=0; i<rResult.size(); i++ )
// 			{
// 				rResult[i] -= ((Vector)rVariable)[i];
// 			}	
// 		}
	}
	
	/**
	 *	TO BE REVIEWED!!!
	 */
	void IsotropicLinearElastic::CalculateConstitutiveMatrix(const Vector& StrainVector, Matrix& rResult)
	{
		rResult = mCtangent;
	}
	
	/**
	 *	TO BE TESTED!!!
	 */
	void IsotropicLinearElastic::CalculateStressEig(const Vector& StressVector,
			double& smin, Vector& mineigenvect,
			double& smid, Vector& mideigenvect,
			double& smax, Vector& maxeigenvect)
	{
		CalculateStressEigNonNormalized( StressVector, 
				smin, mineigenvect,
				smid, mideigenvect,
				smax, maxeigenvect);
		SD_MathUtils<double>::Normalize( mineigenvect );
		SD_MathUtils<double>::Normalize( mideigenvect );
		SD_MathUtils<double>::Normalize( maxeigenvect );
	}//CalculateStressEig
	
	/**
	 * TO BE TESTED
	 */ 
	void IsotropicLinearElastic::CalculateStressEigNonNormalized(const Vector& StressVector,
			double& smin, Vector& mineigenvect,
			double& smid, Vector& mideigenvect,
			double& smax, Vector& maxeigenvect)
	{
		Vector Stresses = StressVector;
		Matrix StressMatrix = MathUtils<double>::StressVectorToTensor( Stresses );
		Matrix EigenVectors = ZeroMatrix(3,3);
		double error_tolerance = 1e-10;
		double zero_tolerance = 1e-10;
		//calculating eigenvectors
		SD_MathUtils<double>::EigenVectors( StressMatrix, 
				EigenVectors, 
				error_tolerance, 
				zero_tolerance);
		//storing eigenvalues and eigenvectors into given variables
		smin = StressMatrix(0,0);
		smid = StressMatrix(1,1);
		smax = StressMatrix(2,2);
        mineigenvect.resize(3,false);
        mideigenvect.resize(3,false);
        maxeigenvect.resize(3,false);
		for( int i=0; i<mineigenvect.size(); i++ )
		{
			mineigenvect(i) = EigenVectors(0,i);
			mideigenvect(i) = EigenVectors(1,i);
			maxeigenvect(i) = EigenVectors(2,i);
		}
	}//CalculateStressEigNonNormalized
	
	/**
	 *	TO BE TESTED!!!
	 */
	void IsotropicLinearElastic::PrincipSTRAIN(const Vector& StrainVector,double& eps1, double& eps2, double& eps3)
	{
		Vector Strains = StrainVector;
		Matrix StrainMatrix = SD_MathUtils<double>::StrainVectorToTensor( Strains );
		//preparing matrix of eigenvectors
		Matrix EigenVectors(3,3);
		//setting up tolerances
		double error_tolerance = 1e-10;
		double zero_tolerance = 1e-10;
		//calculating eigenvalues
		SD_MathUtils<double>::EigenVectors( StrainMatrix, 
				EigenVectors, 
				error_tolerance, 
				zero_tolerance);
		eps1 = StrainMatrix(0,0);
		eps2 = StrainMatrix(1,1);
		eps3 = StrainMatrix(2,2);
	}
	
	/**
	 *	TO BE TESTED!!!
	 */
	void IsotropicLinearElastic::PrincipSTRESS(const Vector& StressVector,double& str1, double& str2, double& str3)
	{
		//stress components
		const double& s1 = StressVector[0];
		const double& s2 = StressVector[1];
		const double& s3 = StressVector[2];
		const double& t12 = StressVector[3];
		const double& t23 = StressVector[4];
		const double& t13 = StressVector[5];
		//invariants
		const double I1 = s1+s2+s3;
		const double I2 = s1*s2+s2*s3+s1*s3-t12*t12-t23*t23-t13*t13;
		const double I3 = s1*s2*s3+2*t12*t23*t13-s1*t23*t23-s2*t13*t13-s3*t12*t12;
		//iteration parameters
		int maxIter = 500;
		double tolerance = 1e-12;
		double s = I1;
		//finding first solution of the problem with newton method
		for( int i=0; i<maxIter; i++ )
		{
			s -= (s*s*s-I1*s*s+I2*s-I3)/(3.0*s*s-2.0*I1*s+I2);
			if( fabs(s*s*s-I1*s*s+I2*s-I3) < tolerance ) break;
		}
		//solving principal stresses
		double a = 3*s-I1;
		double b = 3*s*s-2*I1*s+I2;
		double p1 = -a/2+sqrt(fabs((a*a/4)-b));
		double p2 = -a/2-sqrt(fabs((a*a/4)-b));
		double sigma1 = p1+s;
		double sigma2 = p2+s;
		//sorting stresses
		str1 = s;
		if( str1 >= sigma1 )
		{
			str2 = str1;
			str1 = sigma1;
		}
		if( str2 >= sigma2 )
		{
			str3 = str2;
			str2 = sigma2;
			if( str1 >= str2 )
			{
				s = str1;
				str1 = str2;
				str2 = s;
			}
		}
	}//PrincipSTRESS
	
	/**
	 *	:TODO: TO BE TESTED!!!
	 *	Tangent matrix equals linear elastic matrix because of linear elastic
	 *	material definition
	 */
	void IsotropicLinearElastic::CalculateTangentMatrix(const Vector& StrainVector)
	{
		//does nothing but resetting the mEw vector
		noalias(mEw) = ZeroVector(6);
	}
	
	//**********************************************************************
	void IsotropicLinearElastic::UpdateMaterial( const Vector& StrainVector,
			const Properties& props,
			const GeometryType& geom,
			const Vector& ShapeFunctionsValues ,
			const ProcessInfo& CurrentProcessInfo)
	{
		//         does nothing special
// 		mCtangent = CalculateElasticMatrix(mE,mNU);
	}
	
	/**
	 *	:TODO: TO BE TESTED!!!
	 */
	void IsotropicLinearElastic::Calculate(const Variable<Matrix >& rVariable, Matrix& output)
    {
        if(rVariable==AUXILIARY_MATRIX_1)
        {
            output(0,0) = mEw[0];
            output(0,1) = mEw[1];
            output(0,2) = mEw[2];
            output(0,3) = mEw[3]*0.5;
            output(0,4) = mEw[4]*0.5;
            output(0,5) = mEw[5]*0.5;
        }
        KRATOS_ERROR(std::logic_error,"Trying to Calculate an inexisting variable" , "");
    }
    
    //**********************************************************************
    void IsotropicLinearElastic::CalculateCauchyStresses(
		Vector& rCauchy_StressVector,
		const Matrix& rF,
		const Vector& rPK2_StressVector,
		const Vector& rGreenLagrangeStrainVector)
    {
        Matrix F = rF;
        Vector PK2_StressVector = rPK2_StressVector;
        Matrix S = MathUtils<double>::StressVectorToTensor( PK2_StressVector );
        S = SD_MathUtils<double>::Mult( F, S );
        double J = MathUtils<double>::Det3( F );
        F = SD_MathUtils<double>::Transpose( F );
        S = SD_MathUtils<double>::Mult( S, F );
        SD_MathUtils<double>::Mult( S, J );
    }
} // Namespace Kratos
