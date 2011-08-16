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
*   Date:                $Date: 2008-10-23 12:22:22 $
*   Revision:            $Revision: 1.1 $
*
* ***********************************************************/

// System includes 
#include <iostream>

// External includes 
#include<cmath>

// Project includes 

#include "includes/define.h"
#include "constitutive_laws/continuous_damage_3d.h"
#include "includes/constitutive_law.h"
#include "utilities/math_utils.h"
#include "includes/variables.h"
#include "includes/process_info.h"
#include "nonlinear_materials_application.h"
#include "includes/properties.h"

#define EULER 2.718281828 

namespace Kratos
{
       


    namespace ContinuousDamage3DAuxiliaries
    {
        boost::numeric::ublas::bounded_matrix<double,3,3> mstemp;
        #pragma omp threadprivate(mstemp)
        boost::numeric::ublas::bounded_matrix<double,3,3> msaux;
        #pragma omp threadprivate(msaux)
    }
    using namespace ContinuousDamage3DAuxiliaries;


	/**
	 *	TO BE TESTED!!!
	 */
	ContinuousDamage3D::ContinuousDamage3D() 
	: ConstitutiveLaw< Node<3> >()
	{
	}
	/**
	 *	TO BE TESTED!!!
	 */
	ContinuousDamage3D::~ContinuousDamage3D()
	{
	}
	
	bool ContinuousDamage3D::Has( const Variable<double>& rThisVariable )
	{
		return false;
	}
	
	bool ContinuousDamage3D::Has( const Variable<Vector>& rThisVariable )
	{
		if( rThisVariable == INSITU_STRESS )
			return true;
		return false;
	}
	
	bool ContinuousDamage3D::Has( const Variable<Matrix>& rThisVariable )
	{
		return false;
	}
	
	double ContinuousDamage3D::GetValue( const Variable<double>& rThisVariable )
	{
	    KRATOS_ERROR(std::logic_error, "Vector Variable case not considered" , "");
	}
	
	Vector ContinuousDamage3D::GetValue( const Variable<Vector>& rThisVariable )
	{
		if( rThisVariable == INSITU_STRESS )
			return mInSituStress;
                        
                   if( rThisVariable == INTERNAL_VARIABLES )
                {
                    Vector dummy(5);                
                    //dummy(0): damage 
   		    dummy[0] = mD;                   
                    //dummy(1): stress 1
    		    dummy[1] = mCurrentStress(0);                   
                    //dummy(2): stress 2
      		    dummy[2] = mCurrentStress(1);                  
                    //dummy(3): stress 3  
         	    dummy[3] = mCurrentStress(2);       
                    //dummy(4): process variable
                    dummy[4] = mALPHA;
                    
                    return( dummy );                               
		}  	      
                        
                        
               
               
                
	    KRATOS_ERROR(std::logic_error, "Vector Variable case not considered", "");
	}
	
	Matrix ContinuousDamage3D::GetValue( const Variable<Matrix>& rThisVariable )
	{
	    KRATOS_ERROR(std::logic_error,"Vector Variable case not considered", "");
	}

	void ContinuousDamage3D::SetValue( const Variable<double>& rThisVariable, const double& rValue, 
								   const ProcessInfo& rCurrentProcessInfo )
	{
	}
	
	void ContinuousDamage3D::SetValue( const Variable<array_1d<double, 3> >& rThisVariable, 
								   const array_1d<double,3>& rValue, 
								   const ProcessInfo& rCurrentProcessInfo )
	{
	}
	
	void ContinuousDamage3D::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
								 const ProcessInfo& rCurrentProcessInfo )
	{
		if( rThisVariable == INSITU_STRESS )
		{
			mInSituStress = rValue;
		}
	}
	
	void ContinuousDamage3D::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
								 const ProcessInfo& rCurrentProcessInfo )
	{
	}
	
	/**
	 *	TO BE TESTED!!!
	 */
	void ContinuousDamage3D::InitializeMaterial( const Properties& props,
			const GeometryType& geom,
			const Vector& ShapeFunctionsValues )
	{
                std::cout << "in InitializeMaterial" << std::endl;
		
		mCurrentStress = ZeroVector(6);
		mdevStrainVector=ZeroVector(6);
		mC = ZeroMatrix(6,6);
		mCtangent = ZeroMatrix(6,6);
                
                //compression and shear modulus
                mK = props[MATERIAL_PARAMETERS][0]/(3.0*(1.0-2.0*props[MATERIAL_PARAMETERS][1]));
		mMU = props[MATERIAL_PARAMETERS][0]/(2.0*(1.0+props[MATERIAL_PARAMETERS][1]));

                //initialize process variable alpha, damage parameter d and gradient of damage parameter d
                mDPartial=0.0;
		mALPHA=0.0;
		mD=0.0;

                mInSituStress = ZeroVector(6);
		mStressVector = ZeroVector(6);
		
		CalculateElasticMatrix( mC, mK, mMU );
		
		std::cout << "InitializeMaterial done" << std::endl;

	}
	
	void ContinuousDamage3D::InitializeSolutionStep( const Properties& props,
				     const GeometryType& geom, //this is just to give the array of nodes
				     const Vector& ShapeFunctionsValues ,
				     const ProcessInfo& CurrentProcessInfo)
	{
	}
			
	void ContinuousDamage3D::FinalizeSolutionStep( const Properties& props,
				   const GeometryType& geom, //this is just to give the array of nodes
				   const Vector& ShapeFunctionsValues ,
				   const ProcessInfo& CurrentProcessInfo)
	{
		std::cout<<"this is in FinalizeSolutionStep###################"<<std::endl;
		
		if( CurrentProcessInfo[CALCULATE_INSITU_STRESS] )
		{
			mInSituStress -= mCurrentStress;
			//SetValue( INSITU_STRESS, mInSituStress, CurrentProcessInfo );
		}
		
		//mC=mCtangent;
	}
    
    void ContinuousDamage3D::UpdateMaterial( const Vector& StrainVector,
                                      const Properties& props,
                                      const GeometryType& geom,
                                      const Vector& ShapeFunctionsValues,
                                      const ProcessInfo& CurrentProcessInfo )
        {
	        std::cout << "in UpdateMaterial" << std::endl;
	    
                //calculate trace of strains
		double norm_strains = 0.0;
		for( int i=0; i<3; i++ )
			norm_strains += 1.0/3.0*StrainVector(i);		
		 
		if ( geom.WorkingSpaceDimension() != 3 )
		{
			KRATOS_ERROR(std::logic_error,"This constitutive law is defined for 3D only!" , "");
		}
		
		// calculate deviatoric elastic strains
                Vector devStrainVector = ZeroVector(6);
		
		for ( int i=0; i<3; i++ )
		{
			devStrainVector(i) = StrainVector(i) - norm_strains;    
		}
		for ( int i=3; i<6; i++ )
		{
			devStrainVector(i) =  StrainVector(i);
		}
			
		// calculate deviatoric norm strains sqare : ||dev E||²
		double devStrainVector_norm_sq = (devStrainVector(0)*devStrainVector(0))
		                                +(devStrainVector(1)*devStrainVector(1))
				                +(devStrainVector(2)*devStrainVector(2))
				           +2.0*((devStrainVector(3)*devStrainVector(3))
				                +(devStrainVector(4)*devStrainVector(4))
				                +(devStrainVector(5)*devStrainVector(5)));		
						
		// calculation of deviatoric norm strains : ||dev E||					      
                double devStrainVector_norm = sqrt( devStrainVector_norm_sq );

		// damage parameter		
		mDINF = props[MATERIAL_PARAMETERS][9];
		mDSO = props[MATERIAL_PARAMETERS][10];

	        //calculate damage function
		double phi = (0.5*devStrainVector_norm_sq)-(0.5*(mALPHA*mALPHA));

		bool lambda = false;
		 
                //check for damage criterion
		if( phi > 0.00000001 )
		{
		std::cout << "Yield function > 0:Damage: " << phi << std::endl;
		
		lambda = true;
                //calculate process variable alpha
 	        mALPHA = devStrainVector_norm;

                //calculate gradient of damage funktion D'(alpha)
	        mDPartial = (mDINF/mDSO)*exp(-mALPHA/mDSO);

		}
		else
		{
		std::cout << "Yield function <= 0:No Damage: " << phi << std::endl;
		
		lambda = false;
                
                //calculate process variable alpha
		mALPHA = mALPHA;
		
		}
		
                //calculate damage funktion D(alpha)
		mD = mDINF*(1.0-pow(EULER,(-mALPHA/mDSO)));
		
                //calculate stress vector
		mStressVector[0] = mK*norm_strains+2.0*mMU*(1.0-mD)*devStrainVector(0);
		mStressVector[1] = mK*norm_strains+2.0*mMU*(1.0-mD)*devStrainVector(1);
		mStressVector[2] = mK*norm_strains+2.0*mMU*(1.0-mD)*devStrainVector(2);
		mStressVector[3] = 2.0*mMU*(1.0-mD)*devStrainVector(3);
		mStressVector[4] = 2.0*mMU*(1.0-mD)*devStrainVector(4);
		mStressVector[5] = 2.0*mMU*(1.0-mD)*devStrainVector(5);
                std::cout<<"Stress in Updatematerial"<<std::endl;
		KRATOS_WATCH(mStressVector);
		
                //calculate constitutive damage matrix Cdamage
		mC = ZeroMatrix(6,6);	
		//setting up material matrix
		double c1 = mK-(1.0/3.0)*2.0*mMU*(1.0-mD);	
		double c2 = 0.5*2.0*mMU*(1.0-mD);
		double c3 = 2.0*mMU*(1.0-mD);
		//filling material matrix
		mC(0,0) = c1+c3;    mC(0,1) = c1;      mC(0,2) = c1;       mC(0,3) = 0.0;       mC(0,4) = 0.0;       mC(0,5) = 0.0;
		mC(1,0) = c1;       mC(1,1) = c1+c3;   mC(1,2) = c1;       mC(1,3) = 0.0;       mC(1,4) = 0.0;       mC(1,5) = 0.0;
		mC(2,0) = c1;       mC(2,1) = c1;      mC(2,2) = c1+c3;    mC(2,3) = 0.0;       mC(2,4) = 0.0;       mC(2,5) = 0.0;
		mC(3,0) = 0.0;      mC(3,1) = 0.0;     mC(3,2) = 0.0;      mC(3,3) = c2;        mC(3,4) = 0.0;       mC(3,5) = 0.0;
		mC(4,0) = 0.0;      mC(4,1) = 0.0;     mC(4,2) = 0.0;      mC(4,3) = 0.0;       mC(4,4) = c2;        mC(4,5) = 0.0;
		mC(5,0) = 0.0;      mC(5,1) = 0.0;     mC(5,2) = 0.0;      mC(5,3) = 0.0;       mC(5,4) = 0.0;       mC(5,5) = c2;
		
                //calculate if damage occurs
		if ( lambda == true )
		{
		
		for ( int i=0; i<6; i++ )
		   for ( int j=0; j<6; j++ )
			mC(i,j)=mC(i,j)-((2.0*mMU*mDPartial/devStrainVector_norm)*(devStrainVector(i)*devStrainVector(j)));
		
		}
		KRATOS_WATCH(mC);       
        }
		
	/**
	 *	TO BE TESTED!!!
	 */
	void ContinuousDamage3D::CalculateElasticMatrix(Matrix& C,double K,double MU )
	{ 
		mC = ZeroMatrix(6,6);	
		//setting up material matrix
		
		double c1 = K-(1.0/3.0)*2.0*MU*(1.0-mD);	
		double c2 = 0.5*2.0*MU*(1.0-mD);
		double c3 = 2.0*MU*(1.0-mD);
		//filling material matrix
		C(0,0) = c1+c3;    C(0,1) = c1;      C(0,2) = c1;       C(0,3) = 0.0;       C(0,4) = 0.0;       C(0,5) = 0.0;
		C(1,0) = c1;       C(1,1) = c1+c3;   C(1,2) = c1;       C(1,3) = 0.0;       C(1,4) = 0.0;       C(1,5) = 0.0;
		C(2,0) = c1;       C(2,1) = c1;      C(2,2) = c1+c3;    C(2,3) = 0.0;       C(2,4) = 0.0;       C(2,5) = 0.0;
		C(3,0) = 0.0;      C(3,1) = 0.0;     C(3,2) = 0.0;      C(3,3) = c2;        C(3,4) = 0.0;       C(3,5) = 0.0;
		C(4,0) = 0.0;      C(4,1) = 0.0;     C(4,2) = 0.0;      C(4,3) = 0.0;       C(4,4) = c2;        C(4,5) = 0.0;
		C(5,0) = 0.0;      C(5,1) = 0.0;     C(5,2) = 0.0;      C(5,3) = 0.0;       C(5,4) = 0.0;       C(5,5) = c2;

	}
	
	/**
	 *	TO BE TESTED!!!
	 */
	void ContinuousDamage3D::CalculateStress(const Vector& StrainVector, Vector& StressVector)
	{
	      
		if( StressVector.size() != 6 )
		{
			StressVector.resize(6);
		}

		noalias(StressVector) = mStressVector;

		mCurrentStress = mStressVector;
                std::cout<<"Stress in CalculateStress"<<std::endl;
		KRATOS_WATCH( StressVector );

	}
	
	/**
	 *	TO BE REVIEWED!!!
	 */
	void ContinuousDamage3D::CalculateConstitutiveMatrix(const Vector& StrainVector, Matrix& rResult)
	{
		noalias(rResult) = mC;
	}
	
 
    
        void ContinuousDamage3D::CalculateStressAndTangentMatrix( Vector& StressVector,
                                          const Vector& StrainVector,
                                          Matrix& algorithmicTangent)
    {
    	noalias(StressVector) = mCurrentStress;
	noalias(algorithmicTangent) = mCtangent;
	return;
    }
    
      int ContinuousDamage3D::Check(const Properties& props, const GeometryType& geom, const ProcessInfo& CurrentProcessInfo)
       {
	   if(YOUNG_MODULUS.Key() == 0 || props[YOUNG_MODULUS]<= 0.00)
                KRATOS_ERROR(std::invalid_argument,"YOUNG_MODULUS has Key zero or invalid value ","");

	    const double& nu = props[POISSON_RATIO];
	    const bool check = bool( (nu >0.499 && nu<0.501 ) || (nu < -0.999 && nu > -1.01 ) );
	    if(POISSON_RATIO.Key() == 0 || check==true) // props[POISSON_RATIO] == 1.00 || props[POISSON_RATIO] == -1.00)
                KRATOS_ERROR(std::invalid_argument,"POISSON_RATIO has Key zero invalid value ","");
	  
	    if(DENSITY.Key() == 0 || props[DENSITY]<0.00)
                KRATOS_ERROR(std::invalid_argument,"DENSITY has Key zero or invalid value ","");
	    
	    if(MATERIAL_PARAMETERS.Key() == 0 || props[MATERIAL_PARAMETERS][0]<0.00 || props[MATERIAL_PARAMETERS][1]<0.00 ) 
               KRATOS_ERROR(std::invalid_argument,"MATERIAL_PARAMETERS has Key zero or invalid value ","");
	    
	    return 0;
	    
         }
    
    
} // Namespace Kratos
