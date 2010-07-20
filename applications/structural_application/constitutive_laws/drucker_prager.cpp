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
*   Revision:            $Revision: 1.6 $
*
* ***********************************************************/

// System includes 
#include <iostream>

// External includes 
#include<cmath>

// Project includes 

#include "includes/define.h"
#include "constitutive_laws/drucker_prager.h"

#include "includes/constitutive_law.h"

#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"
#include "includes/variables.h"
#include "includes/process_info.h"
#include "structural_application.h"
#include "includes/properties.h"

namespace Kratos
{
	/**
	 *	TO BE TESTED!!!
	 */
	DruckerPrager::DruckerPrager() 
	: ConstitutiveLaw()
	{
	}
	/**
	 *	TO BE TESTED!!!
	 */
	DruckerPrager::~DruckerPrager()
	{
	}
	
	bool DruckerPrager::Has( const Variable<double>& rThisVariable )
	{
		if( rThisVariable == DP_EPSILON )
			return true;
		return false;
	}
	
	bool DruckerPrager::Has( const Variable<Vector>& rThisVariable )
	{
		if( rThisVariable == INSITU_STRESS )
			return true;
		return false;
                if( rThisVariable == INTERNAL_VARIABLES )
                        return true;
	}
	
	bool DruckerPrager::Has( const Variable<Matrix>& rThisVariable )
	{
		return false;
	}
	
// 	template<class TVariableType> bool DruckerPrager::Has(const TVariableType& rThisVariable)
// 	{
// 		if( rThisVariable == DP_EPSILON )
// 			return true;
// 		if( rThisVariable == INSITU_STRESS )
// 			return true;
// 		return false;
// 	}
// 	
// 	template<class TVariableType> typename TVariableType::Type const& DruckerPrager::GetValue( 
// 			const Variable<TVariableType>& rThisVariable )
// 	{
// 		if( rThisVariable == DP_EPSILON )
// 			return mEpsilon;
// 		if( rThisVariable == INSITU_STRESS )
// 		{
// 			return mInSituStress;
// 		}
// 	}
	
	double& DruckerPrager::GetValue( const Variable<double>& rThisVariable, double& rValue )
	{
		if( rThisVariable == DP_EPSILON )
			rValue = mEpsilon;
        return rValue;
	}
	
	Vector& DruckerPrager::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
	{
		if( rThisVariable == INSITU_STRESS )
		{
			rValue = mInSituStress;
		}
        else if( rThisVariable == INTERNAL_VARIABLES )
        {
            rValue = ZeroVector(1);
            rValue[0] = mEpsilon;
        }
        return( rValue );
	}
	
	Matrix& DruckerPrager::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
	{
        return( rValue );
	}
		
	void DruckerPrager::SetValue( const Variable<double>& rThisVariable, const double& rValue, 
									 const ProcessInfo& rCurrentProcessInfo )
	{
		if( rThisVariable == DP_EPSILON )
		{
			mEpsilon = rValue;
		}
	}
	
	void DruckerPrager::SetValue( const Variable<array_1d<double, 3> >& rThisVariable, 
									 const array_1d<double, 3>& rValue, const ProcessInfo& rCurrentProcessInfo 
								   )
	{
	}
	void DruckerPrager::SetValue( const Variable<Vector>& rThisVariable, 
									 const Vector& rValue, const ProcessInfo& rCurrentProcessInfo 
								   )
	{
		if( rThisVariable == INSITU_STRESS )
		{
			mInSituStress = rValue;
		}
	}
	
	void DruckerPrager::SetValue( const Variable<Matrix>& rThisVariable, 
									 const Matrix& rValue, const ProcessInfo& rCurrentProcessInfo 
								   )
	{
	}

	/**
	 *	TO BE TESTED!!!
	 */
	void DruckerPrager::InitializeMaterial( const Properties& props,
			const GeometryType& geom,
			const Vector& ShapeFunctionsValues )
    {
        ResetMaterial( props, geom, ShapeFunctionsValues );
    }
            
    void DruckerPrager::ResetMaterial( const Properties& props, const GeometryType& geom, const Vector& ShapeFunctionsValues )
	{
// 		std::cout << "Initialising DruckerPrager" << std::endl;
		mCurrentStress = ZeroVector(6);
		mPlasticStrain = ZeroVector(6);
		mE = props[YOUNG_MODULUS];
		mNU = props[POISSON_RATIO];
		mALPHA1 = props[DP_ALPHA1];
		mK = props[DP_K];
		mEpsilon = 0.0;
		mDeltaEpsilon = 0.0;
		mInSituStress = ZeroVector(6);
		
		mCelastic = ZeroMatrix(6,6);
		mCtangent = ZeroMatrix(6,6);
		//setting up material matrix
		double c1 = mE / ((1.00+mNU)*(1-2*mNU));
		double c2 = c1 * (1-mNU);
		double c3 = c1 * mNU;
		double c4 = c1 * 0.5 * (1 - 2*mNU);
		//filling material matrix
		mCelastic(0,0) = c2;    mCelastic(0,1) = c3;    mCelastic(0,2) = c3;    mCelastic(0,3) = 0.0;   mCelastic(0,4) = 0.0;   mCelastic(0,5) = 0.0;
		mCelastic(1,0) = c3;    mCelastic(1,1) = c2;    mCelastic(1,2) = c3;    mCelastic(1,3) = 0.0;   mCelastic(1,4) = 0.0;   mCelastic(1,5) = 0.0;
		mCelastic(2,0) = c3;    mCelastic(2,1) = c3;    mCelastic(2,2) = c2;    mCelastic(2,3) = 0.0;   mCelastic(2,4) = 0.0;   mCelastic(2,5) = 0.0;
		mCelastic(3,0) = 0.0;   mCelastic(3,1) = 0.0;   mCelastic(3,2) = 0.0;   mCelastic(3,3) = c4;    mCelastic(3,4) = 0.0;   mCelastic(3,5) = 0.0;
		mCelastic(4,0) = 0.0;   mCelastic(4,1) = 0.0;   mCelastic(4,2) = 0.0;   mCelastic(4,3) = 0.0;   mCelastic(4,4) = c4;    mCelastic(4,5) = 0.0;
		mCelastic(5,0) = 0.0;   mCelastic(5,1) = 0.0;   mCelastic(5,2) = 0.0;   mCelastic(5,3) = 0.0;   mCelastic(5,4) = 0.0;   mCelastic(5,5) = c4;
	}
	
	void DruckerPrager::InitializeSolutionStep( const Properties& props,
			const GeometryType& geom, //this is just to give the array of nodes
			const Vector& ShapeFunctionsValues ,
			const ProcessInfo& CurrentProcessInfo)
	{
	}
	
	void DruckerPrager::FinalizeSolutionStep( const Properties& props,
			const GeometryType& geom, //this is just to give the array of nodes
			const Vector& ShapeFunctionsValues ,
			const ProcessInfo& CurrentProcessInfo)
	{
		if( CurrentProcessInfo[CALCULATE_INSITU_STRESS] )
		{
			mInSituStress -= mCurrentStress;
			SetValue( INSITU_STRESS, mInSituStress, CurrentProcessInfo );
		}
	}
	
	
	/**
	 *	TO BE TESTED!!!
	 */
	void DruckerPrager::CalculateStress(const Vector& StrainVector, Vector& StressVector)
	{
		if( StressVector.size() != 6 )
		{
			StressVector.resize(6);
		}
		StressVector = mCurrentStress;
		noalias(StressVector) -= mInSituStress;
	}
	
    void DruckerPrager::Calculate( const Variable<Vector >& rVariable, Vector& rResult, 
                                   const ProcessInfo& rCurrentProcessInfo )
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
	void DruckerPrager::CalculateConstitutiveMatrix(const Vector& StrainVector, Matrix& rResult)
	{
		rResult = mCtangent;
	}
	
	/**
	 *	TO BE TESTED!!!
	 */
	void DruckerPrager::CalculateStressEig(const Vector& StressVector,
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
	void DruckerPrager::CalculateStressEigNonNormalized(const Vector& StressVector,
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
		mineigenvect.resize(3);
		mideigenvect.resize(3);
		maxeigenvect.resize(3);
		for( unsigned int i=0; i<mineigenvect.size(); i++ )
		{
			mineigenvect(i) = EigenVectors(0,i);
			mideigenvect(i) = EigenVectors(1,i);
			maxeigenvect(i) = EigenVectors(2,i);
		}
	}//CalculateStressEigNonNormalized
	
	/**
	 *	TO BE TESTED!!!
	 */
	void DruckerPrager::PrincipSTRAIN(const Vector& StrainVector,double& eps1, double& eps2, double& eps3)
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
	void DruckerPrager::PrincipSTRESS(const Vector& StressVector,double& str1, double& str2, double& str3)
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
	void DruckerPrager::CalculateTangentMatrix(const Vector& StrainVector)
	{
		//does nothing but resetting the mEw vector
		noalias(mEw) = ZeroVector(6);
	}
	
	//**********************************************************************
	void DruckerPrager::UpdateMaterial( const Vector& StrainVector,
			const Properties& props,
			const GeometryType& geom,
			const Vector& ShapeFunctionsValues ,
			const ProcessInfo& CurrentProcessInfo)
	{
		// setting limits for return-map
		double tolerance = 1.0e-11;
		int maxiter = 15;
		
		// calculate elastic matrix
		mCtangent = mCelastic;
		// calculate trial stresses from current elastic strains and existing plastic strains
		if ( geom.WorkingSpaceDimension() != 3 )
		{
			KRATOS_ERROR(std::logic_error,"This constitutive law is defined for 3D only!" , "");
		}
		Vector trialStressVector = ZeroVector(6);
		noalias(trialStressVector) = prod(mCtangent,(StrainVector-mPlasticStrain));
// 		KRATOS_WATCH( trialStressVector );
		//in-situ stress is calculated with elastic loading only
// 		if( CurrentProcessInfo[CALCULATE_INSITU_STRESS] )
// 		{
// 			std::cout << "INSITU_STRESS: only elastic loading is calculated" << std::endl;
// 			mCurrentStress = trialStressVector;
// 			return;
// 		}
		//subtract insitu-stress from trial stresses
		trialStressVector -= mInSituStress;
		// calculation of deviatoric stresses and invariants
		double I1_tr = trialStressVector[0]+trialStressVector[1]+trialStressVector[2];
		double J2_tr = 0.0;
		for( unsigned int i=0; i<3; i++ )
		{
			J2_tr += (trialStressVector[i]-I1_tr/3.0)*(trialStressVector[i]-I1_tr/3.0);
			J2_tr += 2.0*trialStressVector[i+3]*trialStressVector[i+3];
		}
		J2_tr = sqrt(J2_tr);
		// calculate yield strength for deviatoric stress state
		// not necessary as we leave out the hardening
		
		// calculate trial yield functions
		double f1_tr = mALPHA1*I1_tr + J2_tr - mK;
		double f2_tr = mALPHA1*I1_tr - mK;
		
		// testing if plasticity occurs
		if( f1_tr >= 0.0 )
		{
// 			std::cout << "Yield function1: " << f1_tr << std::endl;
			//moduli
			double modK = mE/(3.0*(1.0-2.0*mNU));
			double modG = mE/(2.0*(1.0+mNU));
			//deviatoric operator
			Matrix iDev = ZeroMatrix(6,6);
			for( int i=0; i<3; i++ )
			{
				for( int j=0; j<3; j++ )
				{
					if( i == j  ) iDev(i,i) = 2.0/3.0;
					else iDev(i,j) = -1.0/3.0;
				}
				iDev(i+3,i+3) = 0.5;
			}
			//consistency parameter gamma
			double gamma = 0.0;
			double gamma1 = J2_tr/2.0/modG;
			double gamma2 = 0.0;
			double residual = 0.0;
			for( int iter = 0; iter < maxiter; iter++ )
			{
				residual = f1_tr - gamma*(9.0*modK*mALPHA1*mALPHA1 + 2.0*modG);
				if( fabs( residual ) < tolerance )
				{
// 					std::cout << "gamma converged within " << iter << "iterations" << std::endl;
					break;
				}
				//update
				
				gamma += residual/(9.0*modK*mALPHA1*mALPHA1+2.0*modG);
			}
			//if two yield surfaces are active
			if( f2_tr >= 0.0 )
			{
// 				std::cout << "Yield function2: " << f2_tr << std::endl;
				residual = 0.0;
				for( int iter = 0; iter < maxiter; iter++ )
				{
					mDeltaEpsilon = gamma1+gamma2/mALPHA1;
					residual = I1_tr-mK/mALPHA1 - gamma1*9.0*modK*mALPHA1 - gamma2*9.0*modK;
					if( fabs( residual ) < tolerance )
					{
// 						std::cout << "gamma2 converged within " << iter << "iterations" << std::endl;
						break;
					}
					gamma2 += residual/(9.0*modK);
				}
				if( gamma < 0.0 ) KRATOS_ERROR(std::logic_error,"Gamma < 0" , "");
				if( (gamma1 < 0.0) || (gamma2 < 0.0) )
				{
					gamma1 = gamma;
					gamma2 = 0.0;
					mDeltaEpsilon = gamma;
				}
			}
			else
			{
				gamma1 = gamma;
				gamma2 = 0;
				mDeltaEpsilon = gamma;
			}
			// plastic flow direction
			Vector plastic_flow_direction = ZeroVector(6);
			Vector auxDirection = ZeroVector(6);
			for( unsigned int i=0; i<3; i++ )
			{
				auxDirection[i] = (trialStressVector[i]-I1_tr/3.0)/J2_tr;
				auxDirection[i+3] = trialStressVector[i+3]/J2_tr;
				plastic_flow_direction[i] = auxDirection[i]+mALPHA1;
				plastic_flow_direction[i+3] = auxDirection[i+3];
			}
			//update incremental plastic strains
// 			Vector incrementalPlasticStrains = ZeroVector(6);
// 			for( unsigned int i=0; i<3; i++ )
// 			{
// 				incrementalPlasticStrains[i] = gamma1*plastic_flow_direction[i]+gamma2;
// 				incrementalPlasticStrains[i+3] = 2.0*gamma1*plastic_flow_direction[i+3];
// 			}
			
			
			// update stresses
			Vector auxStresses = ZeroVector(6);
			for( unsigned int i=0; i<3; i++ )
			{
				for( unsigned int j=0; j<3; j++ )
				{
					auxStresses[i] += mCelastic(i,j)*plastic_flow_direction[j];
				}
				auxStresses[i+3] = 2.0*modG*plastic_flow_direction[i+3];
				mCurrentStress[i] = trialStressVector[i] - 0.9999999*gamma1*auxStresses[i] 
						- 0.9999999*gamma2*3.0*modK;
				mCurrentStress[i+3] = trialStressVector[i+3] - 0.9999999*gamma1*auxStresses[i+3];
			}
			// calculate tangent material matrix
			if( gamma2 == 0.0 )
			{
				// tangent material matrix for one active yield surface
				Vector DGammaDe = ZeroVector(6);
				double edum = 9.0*modK*mALPHA1*mALPHA1 + 2.0*modG;
				for( unsigned int i=0; i<3; i++ )
				{
					DGammaDe[i] = 1.0/(edum)*(3.0*modK*mALPHA1
							+2.0*modG*auxDirection[i] );
					DGammaDe[i+3] = 1.0/edum*2.0*modG*auxDirection[i+3];
				}
				for( unsigned int i=0; i<6; i++ )
				{
					for( unsigned int j=0; j<6; j++ )
					{
						mCtangent(i,j) += (-auxStresses[i]*DGammaDe[j]
								-gamma1*4.0*modG*modG/J2_tr*(iDev(i,j)
								-auxDirection[i]*auxDirection[j] )
								  );
					}
				}
			}
			else //for two active yield surfaces
			{
				Vector DGammaDe1 = ZeroVector(6);
				Vector DGammaDe2 = ZeroVector(6);
				Vector vk = ZeroVector(6);
				double edum = 18.0*modK*modG;
				DGammaDe1 = auxDirection;
				for( int i=0; i<3; i++ )
				{
					DGammaDe2[i] = 6*modG*modK/edum - mALPHA1*auxDirection[i];
					DGammaDe2[i+3] = -mALPHA1*auxDirection[i+3];
					vk[i] = 3.0*modK;
				}
				for( unsigned int i=0; i<6; i++ )
				{
					for( unsigned int j=0; j<6; j++ )
					{
						mCtangent(i,j) += (-auxStresses[i]*DGammaDe1[j]
								-gamma1*4.0*modG*modG/J2_tr
								*(iDev(i,j)-auxDirection[i]*auxDirection[j])
								- vk[i]*DGammaDe2[j]
								  );
					}
				}
			}
			// calculation of elastic strains
			Matrix inverseC = ZeroMatrix(6,6);
			SD_MathUtils<double>::InvertMatrix( mCelastic, inverseC );
			Vector eps_el = ZeroVector(6);
			noalias( eps_el ) = prod(inverseC, mCurrentStress );
			// update strains
// 			StrainVector += incrementalPlasticStrains;
			// plastic strains
			mPlasticStrain = StrainVector - eps_el;
// 			std::cout << "Difference in stresses:" << std::endl;
// 			KRATOS_WATCH( mCurrentStress-trialStressVector );
			
			//update internal variable
			mEpsilon += mDeltaEpsilon;
// 			KRATOS_WATCH( mEpsilon );
		}
		else //elastic loading
		{
// 			std::cout << "ELASTIC LOADING: no plasticity occured" << std::endl;
			mCurrentStress = trialStressVector;
		}
        
        //TESTING: adding insitu stress
        mCurrentStress += mInSituStress;
        
// 		KRATOS_WATCH( mPlasticStrain );
		
	}
	
	/**
	 *	:TODO: TO BE TESTED!!!
	 */
    void DruckerPrager::Calculate(const Variable<Matrix >& rVariable, Matrix& output, 
                                  const ProcessInfo& rCurrentProcessInfo)
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
	void DruckerPrager::CalculateCauchyStresses(
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
