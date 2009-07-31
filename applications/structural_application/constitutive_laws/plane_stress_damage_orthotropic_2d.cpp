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
*   Revision:            $Revision: 1.13 $
*
* ***********************************************************/

// System includes 
#include <iostream>

// External includes 
#include<cmath>

// Project includes 

#include "includes/define.h"
#include "constitutive_laws/plane_stress_damage_orthotropic_2d.h"
#include "custom_utilities/sd_math_utils.h"
#include "custom_utilities/maped_space.h"


#include "includes/constitutive_law.h"
#include "utilities/math_utils.h"
#include "includes/variables.h"
#include "includes/process_info.h"
#include "includes/properties.h"


namespace Kratos
{

    namespace Plane_Stress_Damage_Orthotropic_2D_Auxiliaries
    {
        boost::numeric::ublas::bounded_matrix<double,2,2> mstemp;
	#ifdef _OPENMP
	#pragma omp threadprivate(mstemp)
	#endif
        boost::numeric::ublas::bounded_matrix<double,2,2> msaux;
	#ifdef _OPENMP
	#pragma omp threadprivate(msaux)
	#endif
    }
    using namespace Plane_Stress_Damage_Orthotropic_2D_Auxiliaries;


	/**
	 *	TO BE TESTED!!!
	 */
        Plane_Stress_Damage_Orthotropic_2D::Plane_Stress_Damage_Orthotropic_2D()
	{

	}

	Plane_Stress_Damage_Orthotropic_2D::Plane_Stress_Damage_Orthotropic_2D(FluencyCriteriaPointer FluencyCriteria )//)const FluencyCriteriaType& FluencyCriteria)
	: ConstitutiveLaw< Node<3> >()
	{
	   mFluencyCriteria = FluencyCriteria;
          
	}

	/**
	 *	TO BE TESTED!!!
	 */

	Plane_Stress_Damage_Orthotropic_2D::~Plane_Stress_Damage_Orthotropic_2D()
	{
	}
	

	bool Plane_Stress_Damage_Orthotropic_2D::Has( const Variable<double>& rThisVariable )
	{
	    return false;
	}
	
	bool Plane_Stress_Damage_Orthotropic_2D::Has( const Variable<Vector>& rThisVariable )
	{
	   return false;
	}
	
	bool Plane_Stress_Damage_Orthotropic_2D::Has( const Variable<Matrix>& rThisVariable )
	{
	   return false;
	}
	
   
	double Plane_Stress_Damage_Orthotropic_2D::GetValue( const Variable<double>& rThisVariable )
	{
	    KRATOS_ERROR(std::logic_error, "Vector Variable case not considered" , "");
	}
	
	Vector Plane_Stress_Damage_Orthotropic_2D ::GetValue( const Variable<Vector>& rThisVariable )
	{

	   KRATOS_ERROR(std::logic_error, "Vector Variable case not considered" , "");

	}
	
	Matrix Plane_Stress_Damage_Orthotropic_2D ::GetValue( const Variable<Matrix>& rThisVariable )
	{
	  KRATOS_ERROR(std::logic_error, "Vector Variable case not considered" , "");
	}

	void Plane_Stress_Damage_Orthotropic_2D::SetValue( const Variable<double>& rThisVariable, const double rValue, 
								   const ProcessInfo& rCurrentProcessInfo )
	{
	}
	
	void Plane_Stress_Damage_Orthotropic_2D::SetValue( const Variable<array_1d<double, 3> >& rThisVariable, 
								   const array_1d<double,3>& rValue, 
								   const ProcessInfo& rCurrentProcessInfo )
	{
	}
	
	void Plane_Stress_Damage_Orthotropic_2D::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
								 const ProcessInfo& rCurrentProcessInfo )
	{

	}
	
	void Plane_Stress_Damage_Orthotropic_2D::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
								 const ProcessInfo& rCurrentProcessInfo )
	{
	}
	 

	void Plane_Stress_Damage_Orthotropic_2D::InitializeMaterial( const Properties& props,
			const GeometryType& geom,
			const Vector& ShapeFunctionsValues )
	{
                mCurrentStress            = ZeroVector(3);
		mCtangent                 = ZeroMatrix(3,3);
		mInSituStress             = ZeroVector(3);
                mOrtotropic_Elastic_Limit = ZeroVector(3); 
                mOrtotropic_Elastic_Limit =  props(ORTHOTROPIC_ELASTIC_LIMIT);
	        mIsotropic_Elastic_Limit  =  props(ISOTROPIC_ELASTIC_LIMIT);
       
                
                //KRATOS_WATCH(props)
                mFluencyCriteria->InitializeMaterial(props);          
                CalculateElasticMatrix(mCtangent, props[ORTHOTROPIC_YOUNG_MODULUS_2D],props[ORTHOTROPIC_POISSON_RATIO_2D],props[ORTHOTROPIC_ANGLE]);
               
	}
	

	void Plane_Stress_Damage_Orthotropic_2D::InitializeSolutionStep( const Properties& props,
				     const GeometryType& geom, //this is just to give the array of nodes
				     const Vector& ShapeFunctionsValues ,
				     const ProcessInfo& CurrentProcessInfo)
	{
	}
			
	void Plane_Stress_Damage_Orthotropic_2D::FinalizeSolutionStep( const Properties& props,
				   const GeometryType& geom, //this is just to give the array of nodes
				   const Vector& ShapeFunctionsValues ,
				   const ProcessInfo& CurrentProcessInfo)
	{
	}
    
    void Plane_Stress_Damage_Orthotropic_2D::UpdateMaterial( const Vector& StrainVector,
                                      const Properties& props,
                                      const GeometryType& geom,
                                      const Vector& ShapeFunctionsValues,
                                      const ProcessInfo& CurrentProcessInfo )
    {   
    }
		

	void Plane_Stress_Damage_Orthotropic_2D::CalculateElasticMatrix(Matrix& C, const Vector& E, const Vector& NU, const double&  Orthotropic_Angle)
		//setting up material matrix
	    {

		//KRATOS_WATCH(NU)
                //KRATOS_WATCH(E)	

	
		double nu_xy   = NU[0];
		double nu_yx   = NU[1];
		
		double xi = 1.00 - nu_xy*nu_yx;
                
		double dxx =  E[0]/ xi;
		double dxy = (E[0] * nu_yx)/ xi;
	        double dyx = (E[1] * nu_xy)/ xi;	
		double dyy = E[1]  / xi;
                double Gxy = (1+dyx)/E(0) +(1+dxy)/E(1);
                Gxy = 1/Gxy;
                
		if (dxy!=dyx){std::cout<<"Warning: Material is not synmetric"<<std::endl;}

		//filling material matrix
		C(0,0) = dxx;     C(0,1) = dxy;    C(0,2) = 0.00;      
		C(1,0) = dyx;     C(1,1) = dyy;    C(1,2) = 0.00;      
		C(2,0) = 0.00;    C(2,1) = 0.00;   C(2,2) = Gxy;      
		
		Matrix t_matrix(3,3);
                Matrix t_matrix_inv(3,3);
		Plane_Stress_Damage_Orthotropic_2D::CalculateTransformationMatrix(t_matrix,Orthotropic_Angle);
		SD_MathUtils<double>::InvertMatrix(t_matrix, t_matrix_inv);  
                
		Matrix temp(3,3);
                noalias(temp)    = prod(C,trans(t_matrix_inv));
		noalias(C)       = prod(t_matrix_inv,temp);
		
	  }
 	

 	void Plane_Stress_Damage_Orthotropic_2D::CalculateTransformationMatrix(Matrix& T, const double& Orthotropic_Angle)
  	{
  	      double  angle = Orthotropic_Angle*PI/180.00; // angle must be in radians
  	      double  coseno   = cos(angle);
              double  seno   = sin(angle);
  
  	     T(0,0) = coseno*coseno;  T(0,1) = seno*seno;      T(0,2) = 2*coseno*seno;
  	     T(1,0) = seno*seno;      T(1,1) = coseno*coseno;  T(1,2) = -2*coseno*seno;
             T(2,0) = -coseno*seno;   T(2,1) = coseno*seno;    T(2,2) = coseno*coseno-seno*seno;    
 
 	}


	void Plane_Stress_Damage_Orthotropic_2D::CalculateStress( const Vector& StrainVector, Vector& StressVector)
	{
		if( StressVector.size() != 3 )
		{
			StressVector.resize(3,false);
		}

                //double result = 0.00;
		Vector I      = ZeroVector(3);
                Vector J_des  = ZeroVector(3);
                Vector J      = ZeroVector(3);
                Matrix W_Matrix = ZeroMatrix(6,6);
		noalias(StressVector) = prod(mCtangent,StrainVector);
		mCurrentStress = StressVector;
		//noalias(StressVector) -= mInSituStress;

                //Maped_Space<double>::Calculate_Omega_Tensor(mOrtotropic_Elastic_Limit, mIsotropic_Elastic_Limit, W_Matrix);         
		 //mFluencyCriteria->CalculateEquivalentUniaxialStressViaInvariants(StressVector, result);
                 //mFluencyCriteria->CalculateEquivalentUniaxialStressViaPrincipalStress(StressVector, result);
                 //KRATOS_WATCH(result)
	}
	

	void Plane_Stress_Damage_Orthotropic_2D::CalculateConstitutiveMatrix(const Vector& StrainVector, Matrix& rResult)
	{
		rResult = mCtangent;
                
              
	}

	
    //**********************************************************************
    void Plane_Stress_Damage_Orthotropic_2D::CalculateCauchyStresses(
		Vector& rCauchy_StressVector,
		const Matrix& rF,
		const Vector& rPK2_StressVector,
		const Vector& rGreenLagrangeStrainVector)
    {
		Matrix S = MathUtils<double>::StressVectorToTensor( rPK2_StressVector );

		double J = MathUtils<double>::Det2( rF );

		noalias(mstemp) = prod(rF,S);
		noalias(msaux)  = prod(mstemp,trans(rF));
		msaux *= J;

		if(rCauchy_StressVector.size() != 3)
			rCauchy_StressVector.resize(3);
		
		rCauchy_StressVector[0] = msaux(0,0);
		rCauchy_StressVector[1] = msaux(1,1);
		rCauchy_StressVector[2] = msaux(1,2);

    }


//***********************************************************************************************
//***********************************************************************************************	
	 void Plane_Stress_Damage_Orthotropic_2D::CalculateStressAndTangentMatrix(Vector& StressVector,
                    const Vector& StrainVector,
                    Matrix& algorithmicTangent)
   {
			algorithmicTangent = mCtangent;
                        
  }



} // Namespace Kratos
