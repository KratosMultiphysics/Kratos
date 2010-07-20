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
*   Date:                $Date: 2009-01-14 17:14:19 $
*   Revision:            $Revision: 1.2 $
*
* ***********************************************************/

// System includes 
#include <iostream>

// External includes 
#include<cmath>

// Project includes 

#include "includes/define.h"
#include "constitutive_laws/von_mises_3d.h"

#include "includes/constitutive_law.h"

#include "utilities/math_utils.h"
#include "includes/variables.h"
#include "includes/process_info.h"
#include "structural_application.h"
#include "includes/properties.h"

namespace Kratos
{
    namespace VonMises3DAuxiliaries
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
    using namespace VonMises3DAuxiliaries;


	/**
	 *	TO BE TESTED!!!
	 */
	VonMises3D::VonMises3D() 
	: ConstitutiveLaw()
	{
	}
	/**
	 *	TO BE TESTED!!!
	 */
	VonMises3D::~VonMises3D()
	{
	}
	
	bool VonMises3D::Has( const Variable<double>& rThisVariable )
	{
		return false;
	}
	
	bool VonMises3D::Has( const Variable<Vector>& rThisVariable )
	{
		if( rThisVariable == INSITU_STRESS )
			return true;
		return false;
	}
	
	bool VonMises3D::Has( const Variable<Matrix>& rThisVariable )
	{
		return false;
	}
	
	double& VonMises3D::GetValue( const Variable<double>& rThisVariable, double& rValue )
	{
	    if( rThisVariable == PLASTICITY_INDICATOR )
            return( mPlasticityIndicator );
        return rValue;
	}
	
	Vector& VonMises3D::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
	{
		if( rThisVariable == INSITU_STRESS )
			return mInSituStress;
        if( rThisVariable == INTERNAL_VARIABLES )
        {
            rValue.resize( 13, false );
            for( int i=0; i<6; i++ )
                rValue[i] = mPlasticStrains[i];
            for( int i=0; i<6; i++ )
                rValue[i+6] = mRho[i];
            rValue[12] = mAlpha;
            return( rValue );
        }
	    return rValue;
	}
	
	Matrix& VonMises3D::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
	{
	    return rValue;
	}

	void VonMises3D::SetValue( const Variable<double>& rThisVariable, const double& rValue, 
								   const ProcessInfo& rCurrentProcessInfo )
	{
	}
	
	void VonMises3D::SetValue( const Variable<array_1d<double, 3> >& rThisVariable, 
								   const array_1d<double,3>& rValue, 
								   const ProcessInfo& rCurrentProcessInfo )
	{
	}
	
	void VonMises3D::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
								 const ProcessInfo& rCurrentProcessInfo )
	{
		if( rThisVariable == INSITU_STRESS )
		{
			mInSituStress = rValue;
		}
	}
	
	void VonMises3D::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
								 const ProcessInfo& rCurrentProcessInfo )
	{
	}
	
	/**
	 *	TO BE TESTED!!!
	 */
	void VonMises3D::InitializeMaterial( const Properties& props,
                                         const GeometryType& geom,
                                         const Vector& ShapeFunctionsValues )
	{
        ResetMaterial( props, geom, ShapeFunctionsValues );
    }
    
    void VonMises3D::ResetMaterial( const Properties& props, const GeometryType& geom, const Vector& ShapeFunctionsValues )
    {
        mPlasticityIndicator = 0.0;
        mCurrentStress = ZeroVector(6);
        mPlasticStrains = ZeroVector(6);
        mCurrentPlasticStrains = ZeroVector(6);
		mCtangent = ZeroMatrix(6,6);
		mInSituStress = ZeroVector(6);
        mH = props[MATERIAL_PARAMETERS][3];
        mBeta = props[MATERIAL_PARAMETERS][4];
        mAlpha = 0.0;
        mCurrentAlpha = 0.0;
        mRho = ZeroVector(6);
        mCurrentRho = ZeroVector(6);
        mIdev = ZeroMatrix(6,6);
        for( int i=0; i<6; i++ )
            mIdev(i,i) = 1.0;
        for( int i=0; i<3; i++ )
            for( int j=0; j<3; j++ )
                mIdev(i,j) -= 1.0/3.0;
        
    }
	
	void VonMises3D::InitializeSolutionStep( const Properties& props,
                                             const GeometryType& geom,
                                             const Vector& ShapeFunctionsValues ,
                                             const ProcessInfo& CurrentProcessInfo )
    {
    }
			
	void VonMises3D::FinalizeSolutionStep( const Properties& props,
                                           const GeometryType& geom,
                                           const Vector& ShapeFunctionsValues,
                                           const ProcessInfo& CurrentProcessInfo )
    {
        mPlasticStrains = mCurrentPlasticStrains;
        mAlpha = mCurrentAlpha;
        mRho = mCurrentRho;
        if( CurrentProcessInfo[CALCULATE_INSITU_STRESS] )
        {
            mInSituStress -= mCurrentStress;
            //SetValue( INSITU_STRESS, mInSituStress, CurrentProcessInfo );
        }
    }
    
    /**
     * MATERIAL_PARAMETERS:
     * [2]: yield stress
     * [3]: hardening parameter
     * [4]: beta
     */
    void VonMises3D::UpdateMaterial( const Vector& StrainVector,
                                     const Properties& props,
                                     const GeometryType& geom,
                                     const Vector& ShapeFunctionsValues,
                                     const ProcessInfo& CurrentProcessInfo )
    {
        //setting up elastic matrix
        mCelastic = ZeroMatrix(6,6);
        //setting up tangent matrix
        mCtangent = ZeroMatrix(6,6);
        //update current plastic strains
        mCurrentPlasticStrains = mPlasticStrains;
        //update current hardening variables
        mCurrentAlpha = mAlpha;
        mCurrentRho = mRho;
        //calculate elastic matrix
        CalculateElasticMatrix( mCelastic, props[YOUNG_MODULUS], props[POISSON_RATIO] );
        //calculate elastic stress
        mCurrentStress = prod( mCelastic, (StrainVector-mCurrentPlasticStrains) );
        
        //setting up material parameters
        double sigma_yield = props[MATERIAL_PARAMETERS][2];
        double G = props[YOUNG_MODULUS]/(2.0+2.0*props[POISSON_RATIO]);
//         double K = props[YOUNG_MODULUS]/(3.0-6.0*props[POISSON_RATIO]);
        
        //calculate deviatoric stress
        Vector SigmaDev = mCurrentStress;
        for( int i=0; i<3; i++ )
            SigmaDev[i] -= (mCurrentStress[0]+mCurrentStress[1]+mCurrentStress[2])/3.0;
        
        //calculate eta
        Vector eta = SigmaDev - mCurrentRho;
        
        //calculate norm of deviatoric stress
        double EtaNorm = 0.0;
        for( int i=0; i<3; i++ )
        {
            EtaNorm += eta[i]*eta[i];
            EtaNorm += 2.0*eta[i+3]*eta[i+3];
        }
        EtaNorm = sqrt(EtaNorm);
        
        //calculate yield function
        double f = EtaNorm - sqrt(2.0/3.0)*(sigma_yield+(mBeta*mH*mCurrentAlpha));
        
        //check for plasticity
        if( f > 0.0 )
        {
            mPlasticityIndicator = 1.0;
            //calculate projection vector
            Vector N = eta/EtaNorm;
            
            //calculate return map factor
            double dgamma = (EtaNorm-sqrt(2.0/3.0)*(sigma_yield+(mBeta*mH*mCurrentAlpha)))/(2.0*G);
            
            //update alpha
            mCurrentAlpha = mCurrentAlpha+sqrt(2.0/3.0)*dgamma;
            
            //update rho
            mCurrentRho = mCurrentRho + sqrt(2.0/3.0)*(1.0-mBeta)*mH*dgamma*N;
            
            //update plastic strains
            mCurrentPlasticStrains = mCurrentPlasticStrains + dgamma*N;
            
            //update stress
            //mCurrentStress = mCurrentStress - 2.0*G*dgamma*N;
            mCurrentStress = prod( mCelastic, (StrainVector-mCurrentPlasticStrains) );
            
            //update tangent matrix
            mCtangent = mCelastic - 2.0*G*((1.0/(1.0+mBeta*mH/3.0/G))-(2.0*G*dgamma/EtaNorm))*outer_prod(N,N);
            mCtangent = mCtangent - (4.0*G*G*dgamma/EtaNorm)*mIdev;
        }
        else
        {
            mCtangent = mCelastic;
        }
    }
    
    /**
     *	returns current stress
     */
    void VonMises3D::CalculateStress(const Vector& StrainVector, Vector& StressVector)
    {
        if( StressVector.size() != 6 )
        {
            StressVector.resize(6);
        }
        noalias(StressVector) = mCurrentStress;
        noalias(StressVector) -= mInSituStress;
        
//         KRATOS_WATCH( StressVector );
    }
    
    /**
     *	returns material tangent matrix
     */
    void VonMises3D::CalculateConstitutiveMatrix(const Vector& StrainVector, Matrix& rResult)
    {
        rResult = mCtangent;
    }
    
    //**********************************************************************
    void VonMises3D::CalculateCauchyStresses( Vector& rCauchy_StressVector,
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
    
    void VonMises3D::CalculateElasticMatrix(Matrix& C, const double E, const double NU)
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
     * MATERIAL_PARAMETERS:
     * [2]: yielding parameter 1 (yield stress)
     * [3]: yielding parameter 2 (hardening stiffness)
     * [4]: hardening parameter
     */
//     void VonMises3D::mate17( const Properties& props, 
//                              const Vector& epsep )
//     {
//         //material parameters
//         int nhistmat = props[INTERNAL_VARIABLES].size();
//         //modulus of compression
//         double kappa = props[YOUNG_MODULUS]/(3.0-6.0*props[POISSON_RATIO]);
//         //shear modulus
//         double xmu = props[YOUNG_MODULUS]/(2.0*(1+props[POISSON_RATIO]));
//         //yielding parameters
//         double ypar1 = props[MATERIAL_PARAMETERS][8];
//         
//         
//         //reset material tangent matrix
//         noalias(mCtangent) = ZeroMatrix(6,6);
//         
//         //calculate elastic tangent matrix
//         for( int i=0; i<3; i++ )
//         {
//             for( int j=0; j<3; j++ )
//             {
//                 mCtangent(i,j) = kappa - 2.0*xmu/3.0;
//             }
//             mCtangent(i,i) = mCtangent(i,i)+2.0*xmu;
//             mCtangent(i+3,i+3) = mCtangent(i+3,i+3)+xmu;
//         }
//         
//         //get plastic strains
//         Vector epsp = ZeroVector(6);
//         noalias(epsp) = mPlasticStrains;
//         
//         //elastic trial strains
//         Vector epset = epsep - epsp;
//         KRATOS_WATCH(epset);
//         //trace of elastic trial strains
//         double epsettr = 0.0;
//         for( int i=0; i<3; i++ )
//             epsettr += epset[i];
//         
//         //deviatoric elastic trial strains
//         Vector epsetdev = ZeroVector(6);
//         for( int i=0; i<3; i++ )
//         {
//             epsetdev[i] = epset[i]-epsettr/3.0;
//             epsetdev[i+3] = epset[i+3];
//         }
//         
//         //deviatoric trial stresses
//         Vector sigtdev = ZeroVector(6);
//         for( int i=0; i<3; i++ )
//         {
//             sigtdev[i] = 2.0*xmu*epsetdev[i];
//             sigtdev[i+3] = xmu*epsetdev[i+3];
//         }
//         
//         //check for yielding
//         double radn = sqrt(2.0/3.0)*ypar1;
//         double xnorm = 0.0;
//         for( int i=0; i<3; i++ )
//         {
//             xnorm += sigtdev[i]*sigtdev[i];
//             xnorm += 2.0*(sigtdev[i+3]*sigtdev[i+3]);
//         }
//         xnorm = sqrt(xnorm);
//         double ftrial = xnorm - radn;
//         KRATOS_WATCH(ftrial);
//         //iterative return mapping
//         int iter = 1;
//         while( ftrial > 1.0e-10 )
//         {
//             std::cout << "#### RETURN MAP: iteration " << iter << std::endl;
//             //determination of lambdap
//             double dlamp = (ftrial-sqrt(2.0/3.0)*ypar1)/(2.0*xmu+2.0/3.0);
//             //rate of plastic strains
//             double fact = dlamp/xnorm;
//             Vector depsp = fact*sigtdev;
//             
//             //determination of N=devS_trial/|devS_trial|
//             Vector en = 1.0/xnorm*sigtdev;
//             
//             //update plastic strains
//             for( int i=0; i<3; i++ )
//             {
//                 epsp[i] += depsp[i];
//                 epsp[i+3] += 2.0*depsp[i+3];
//             }
//             
//             //update internal variables
//             mCurrentPlasticStrains = epsp;
//             
//             //plastic modification for tangent
//             for( int i=0; i<3; i++ )
//             {
//                 for( int j=0; j<3; j++ )
//                 {
//                     mCtangent(i,j) += -2.0*xmu/3.0*2.0*xmu*fact;
//                 }
//                 mCtangent(i,i) += 2.0*xmu*2.0*xmu*fact;
//                 mCtangent(i+3,i+3) += xmu*2.0*xmu*fact;
//             }
//             double x1 = 2.0*xmu*(1.0-fact);
//             for( int i=0; i<6; i++ )
//                 for( int j=0; j<6; j++ )
//                     mCtangent(i,j) += -x1*en[i]*en[j];
//             
//             //update flow rule
//             epset = epsep - epsp;
//             KRATOS_WATCH(epset);
//             epsettr = 0.0;
//             for( int i=0; i<3; i++ )
//                 epsettr += epset[i];
//             for( int i=0; i<3; i++ )
//             {
//                 epsetdev[i] = epset[i]-epsettr/3.0;
//                 epsetdev[i+3] = epset[i+3];
//             }
//             for( int i=0; i<3; i++ )
//             {
//                 sigtdev[i] = 2.0*xmu*epsetdev[i];
//                 sigtdev[i+3] = xmu*epsetdev[i+3];
//             }
//             xnorm = 0.0;
//             for( int i=0; i<3; i++ )
//             {
//                 xnorm += sigtdev[i]*sigtdev[i];
//                 xnorm += 2.0*(sigtdev[i+3]*sigtdev[i+3]);
//             }
//             ftrial = xnorm - radn;
//             KRATOS_WATCH(ftrial);
//             
//             iter++;
//             if( iter > 10 ) break;
//         }//end of return mapping
//         
//         //elastic strains
//         Vector epse = epsep - epsp;
//         
//         //trace of elastic strains
//         double epsetr = 0.0;
//         for( int i=0; i<3; i++ )
//             epsetr += epse[i];
//         
//         //stresses
//         for( int i=0; i<3; i++ )
//         {
//             mCurrentStress[i] = 2.0*xmu*epse[i]+(kappa-2.0/3.0*xmu)*epsetr;
//             mCurrentStress[i+3] = xmu*epse[i+3];
//         }
//     }
} // Namespace Kratos
