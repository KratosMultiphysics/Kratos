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

THE  SOFTWARE IS  PROVIDED  "AS  variablesIS", WITHOUT  WARRANTY  OF ANY  KIND,
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
*   Last Modified by:    $Author: nelson $
*   Date:                $Date: 2008-09-03 
\*   Revision:            $Revision: 1.1 $
*
* ***********************************************************/


// System includes 
#include <iostream>

// External includes 
#include<cmath>
#include <limits>



// Project includes 
#include "includes/define.h"
#include "constitutive_laws/Isotropic_Damage_3D.h"
#include "includes/constitutive_law.h"
#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"
#include "custom_utilities/tensor_utils.h"
#include "includes/variables.h"
#include "includes/process_info.h"
#include "includes/properties.h"


namespace Kratos
{
    namespace Isotropic_Damage_3D_Auxiliaries
    {
        boost::numeric::ublas::bounded_matrix<double,3,3> mstemp;
	#ifdef _OPENMP
	#pragma omp threadprivate(mstemp)
	#endif
        boost::numeric::ublas::bounded_matrix<double,3,3> msaux;
	#ifdef _OPENMP
	#pragma omp threadprivate(msaux)
	#endif
        Matrix StressTensor(3,3);
        #ifdef _OPENMP
	#pragma omp threadprivate(StressTensor)
	#endif
	Matrix InvertedMatrix(6,6);
        #ifdef _OPENMP
	#pragma omp threadprivate(InvertedMatrix)
        #endif
        Vector PrincipalStress(3);
	#ifdef _OPENMP
	#pragma omp threadprivate(PrincipalStress)
        #endif
	Vector Aux_Vector(6);
        #ifdef _OPENMP
        #pragma omp threadprivate(Aux_Vector)
        #endif
	Matrix ConstitutiveMatrixAux(6,6);
        #ifdef _OPENMP
        #pragma omp threadprivate(ConstitutiveMatrixAux)
        #endif
        Vector StrainVectorPerturbation(6);
	#ifdef _OPENMP
        #pragma omp threadprivate(StrainVectorPerturbation)
        #endif
	Vector StressVectorPerturbation(6);
	#ifdef _OPENMP
        #pragma omp threadprivate(StressVectorPerturbation)
        #endif
        Vector StressVector_Aux(6);
        #ifdef _OPENMP
        #pragma omp threadprivate(StressVector_Aux)
        #endif
	Vector StrainVector_Aux(6);
	#ifdef _OPENMP
        #pragma omp threadprivate(StrainVector_Aux)
        #endif

    } 


    using namespace Isotropic_Damage_3D_Auxiliaries;

	/**
	 *	TO BE TESTED!!!
	 */
	 Isotropic_Damage_3D::Isotropic_Damage_3D() 
	: ConstitutiveLaw< Node<3> >()
	{
	}
	/**
	 *	TO BE TESTED!!!
	 */
	Isotropic_Damage_3D::~Isotropic_Damage_3D()
	{
	}
	
	
	bool Isotropic_Damage_3D::Has( const Variable<double>& rThisVariable )
	{
		return false;
	}
	
	bool Isotropic_Damage_3D::Has( const Variable<Vector>& rThisVariable )
	{
		return false;
	}
	
	bool Isotropic_Damage_3D::Has( const Variable<Matrix>& rThisVariable )
	{
		return false;
	}
	
	double Isotropic_Damage_3D::GetValue( const Variable<double>& rThisVariable )
	{
		 if( rThisVariable == DAMAGE)			
		{
		    return md;
		}
		else
		{ 
		return 0.00;
	    //KRATOS_ERROR(std::logic_error, "double Variable case not considered", "");
		}
	}

	Vector Isotropic_Damage_3D::GetValue( const Variable<Vector>& rThisVariable )
	{
	    KRATOS_ERROR(std::logic_error, "Vector Variable case not considered", "");
	}
	
	Matrix Isotropic_Damage_3D::GetValue( const Variable<Matrix>& rThisVariable )
	{
	    KRATOS_ERROR(std::logic_error, "Matrix Variable case not considered", "");
	}

    void Isotropic_Damage_3D::SetValue( const Variable<double>& rThisVariable, const double& rValue, 
                                const ProcessInfo& rCurrentProcessInfo )
	{
	}
	
    void Isotropic_Damage_3D::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue, 
                                const ProcessInfo& rCurrentProcessInfo )
	{
	}
	
    void Isotropic_Damage_3D::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue, 
                                const ProcessInfo& rCurrentProcessInfo )
	{
	}
    
    void Isotropic_Damage_3D::Calculate(const Variable<Matrix >& rVariable, Matrix& rResult, 
                   const ProcessInfo& rCurrentProcessInfo)
    {
    }


//***********************************************************************************************
//***********************************************************************************************

void Isotropic_Damage_3D::InitializeMaterial( const Properties& props,
	const GeometryType& geom,
	const Vector& ShapeFunctionsValues )
	
{


//  Nota: Estas Variables no seran necesarias almacenarlas por cada elemento.
// Basta con llamarlas con sus propiedades.
  mFc    = props[FC];
  mFt    = props[FT];
  mEc    = props[CONCRETE_YOUNG_MODULUS_C];
  mEt    = props[CONCRETE_YOUNG_MODULUS_T];
  mNU    = props[POISSON_RATIO];
  mGE    = props[FRACTURE_ENERGY];
  ml     = pow(fabs(geom.Volume()),0.333333333333333333);   // longitud del elemento
  mr_old = mFt/sqrt(mEc);
  ////KRATOS_WATCH(geom.Volume());
 // //KRATOS_WATCH(ml);

}

		

//***********************************************************************************************
//***********************************************************************************************

void Isotropic_Damage_3D::InitializeSolutionStep( const Properties& props,
                    const GeometryType& geom,
                    const Vector& ShapeFunctionsValues ,
                    const ProcessInfo& CurrentProcessInfo)
{
				  ////KRATOS_WATCH(ml);  
}

//***********************************************************************************************
//***********************************************************************************************


void Isotropic_Damage_3D::CalculateNoDamageElasticMatrix(Matrix& C, const double E, const double NU)
	{ 
		C.resize(6,6, false);
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
		////KRATOS_WATCH(C);
		
	}
	

//***********************************************************************************************
//***********************************************************************************************

void Isotropic_Damage_3D::CalculateConstitutiveMatrix(const Vector& StrainVector, Matrix& ConstitutiveMatrix)
	{
		
	   //StressVector_Aux.resize(6,false);
	    Isotropic_Damage_3D::CalculateNoDamageElasticMatrix(ConstitutiveMatrix, mEc ,mNU);
	   //noalias(ConstitutiveMatrix) = (1-md)*ConstitutiveMatrix; 
	   //Isotropic_Damage_3D::CalculateStress(StrainVector,StressVector_Aux);
	   //Isotropic_Damage_3D::CalculateStressAndTangentMatrix(StressVector_Aux,StrainVector,ConstitutiveMatrix);
	   //KRATOS_WATCH(ConstitutiveMatrix)

	}
	

//***********************************************************************************************
//***********************************************************************************************


 void Isotropic_Damage_3D::CalculateDamage(const Matrix& ConstitutiveMatrix, const Vector& StressVector, double& d)
 { 
				
    double Tau 	     = 0.00;
    double A         = 0.00;
    double crit      = 1.0e-9;
    double zero      = 1.0e-9;
    double teta_a    = 0.00;
    double teta_b    = 0.00;
    double teta      = 0.00;
    double r         = 0.00;
    double elev      = 0.00; 
    double norma     = 0.00; 
		  	

    if (Aux_Vector.size() != 6)
	{
	      StressTensor.resize(3,3,false);
	      InvertedMatrix.resize(6,6, false);
	      PrincipalStress.resize(3, false);
	      Aux_Vector.resize(6, false);
	}

	
     double ro	= mFt/sqrt(mEc);       //Suponiendo el Ec=Ect 
     double n   = mFc/mFt;
      A         = 1.00/((mGE*mEc)/(ml*mFt*mFt)-0.5);
      ////KRATOS_WATCH(A);
      
      if (A < 0.00)
       {    
	  A  = 0.00;
	  std::cout<< "Warning: A is less than zero"<<std::endl;
       }

   // evitando terminos nulos en diagonal principal y fuera de ellas			
    StressTensor  = MathUtils<double>::StressVectorToTensor(StressVector);
    SD_MathUtils<double>::InvertMatrix(ConstitutiveMatrix, InvertedMatrix);

   if (fabs(StressTensor(0,0))<1E-10){StressTensor(0,0) = 1E-10; }
   if (fabs(StressTensor(0,1))<1E-10){StressTensor(0,1) = 1E-10; }
   if (fabs(StressTensor(0,2))<1E-10){StressTensor(0,2) = 1E-10; }    
   if (fabs(StressTensor(1,0))<1E-10){StressTensor(1,0) = 1E-10; }
   if (fabs(StressTensor(1,1))<1E-10){StressTensor(1,1) = 1E-10; }
   if (fabs(StressTensor(1,2))<1E-10){StressTensor(1,2) = 1E-10; }
   if (fabs(StressTensor(2,0))<1E-10){StressTensor(2,0) = 1E-10; }
   if (fabs(StressTensor(2,1))<1E-10){StressTensor(2,1) = 1E-10; }
   if (fabs(StressTensor(2,2))<1E-10){StressTensor(2,2) = 1E-10; }


   // Calculando norma de una matriz
   norma = SD_MathUtils<double>::normTensor(StressTensor);
   
   //KRATOS_WATCH(StressVector)
   //KRATOS_WATCH(StressTensor)
   ////KRATOS_WATCH(norma)
    if(norma<1E-5)                            
	{
	  PrincipalStress(0)  = 0.00; 
    	  PrincipalStress(1)  = 0.00;
	  PrincipalStress(2)  = 0.00;
    }

    else
    {
      PrincipalStress  = SD_MathUtils<double>::EigenValues(StressTensor,crit, zero);
    }
    //KRATOS_WATCH(PrincipalStress)
    teta_a =  Tensor_Utils<double>::Mc_aully(PrincipalStress);
    teta_b = norm_1(PrincipalStress);  
			
   // Condicion para no tener divison de 0/0 
    if (teta_a==0.00 && teta_b==0.00)
	{teta = 1.00;}
   else
	{teta = teta_a/teta_b;}
   
			
   noalias(Aux_Vector)     = prod(trans(StressVector),InvertedMatrix);
   Tau                     = inner_prod(StressVector,Aux_Vector);
   Tau                     = sqrt(Tau);			
   Tau                     = (teta + (1.00-teta)/n)*Tau;
  	 								
   r     = std::max(mr_old,Tau);  // rold por ro
   elev  = A*(1.00-r/ro);
   d     = 1.00 - (ro/r)*exp(elev);
   if (d < 0.00) 
   {
   d = fabs(d);
   //std::cout<<"Warning: Damage is less than zero"<<std::endl;
   }

   this-> md      = d;     
   this-> mr_new  = r; 

 }



//***********************************************************************************************
//***********************************************************************************************


void  Isotropic_Damage_3D::FinalizeSolutionStep( const Properties& props,
                                               const GeometryType& geom, 
                                               const Vector& ShapeFunctionsValues ,
                                             const ProcessInfo& CurrentProcessInfo)

			{
				 mr_old = mr_new;
				 //KRATOS_WATCH(mr_old);
			}


//***********************************************************************************************
//***********************************************************************************************

void Isotropic_Damage_3D::CalculateNoDamageStress(const Vector& StrainVector, Vector& StressVector)
	{
 		double c1 = mEc / ((1.00+mNU)*(1-2*mNU));
 		double c2 = c1 * (1-mNU);
 		double c3 = c1 * mNU;
 		double c4 = c1 * 0.5 * (1 - 2*mNU);
 
 		StressVector[0] = c2*StrainVector[0] + c3 * (StrainVector[1] + StrainVector[2])	;
 		StressVector[1] = c2*StrainVector[1] + c3 * (StrainVector[0] + StrainVector[2])	;
 		StressVector[2] = c2*StrainVector[2] + c3 * (StrainVector[0] + StrainVector[1])	;
 		StressVector[3] = c4*StrainVector[3];
 		StressVector[4] = c4*StrainVector[4];
 		StressVector[5] = c4*StrainVector[5];

	}
	

//***********************************************************************************************
//***********************************************************************************************
	

 void Isotropic_Damage_3D::CalculateStress( const Vector& StrainVector, Vector& StressVector)
	{

		ConstitutiveMatrixAux.resize(6,6,false);
                StrainVector_Aux.resize(6,false);
		//mstressVector.resize(6,false);
                noalias(StrainVector_Aux) = StrainVector; 
		double d=0.00;
		for (unsigned int i=0;i<StrainVector.size();i++)
			    {   
			      if  (fabs(StrainVector(i))<1E-10)
				  {
				      StrainVector_Aux(i) = 0.00;
				  }
			     }

		Isotropic_Damage_3D::CalculateNoDamageStress(StrainVector_Aux, StressVector);
		Isotropic_Damage_3D::CalculateNoDamageElasticMatrix(ConstitutiveMatrixAux,mEc,mNU);
		Isotropic_Damage_3D::CalculateDamage(ConstitutiveMatrixAux, StressVector, d);
		noalias(StressVector) = (1.00-d)*StressVector;
		//noalias(mstressVector) = StressVector;											
   }


//***********************************************************************************************
//***********************************************************************************************	

  void Isotropic_Damage_3D::CalculateCauchyStresses(
		Vector& rCauchy_StressVector,
		const Matrix& rF,
		const Vector& rPK2_StressVector,
		const Vector& rGreenLagrangeStrainVector)
    {
		Matrix S = MathUtils<double>::StressVectorToTensor( rPK2_StressVector );

		double J = MathUtils<double>::Det3( rF );

		noalias(mstemp) = prod(rF,S);
		noalias(msaux)  = prod(mstemp,trans(rF));
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

//***********************************************************************************************
//***********************************************************************************************	
	 void Isotropic_Damage_3D::CalculateStressAndTangentMatrix(Vector& StressVector,
                    const Vector& StrainVector,
                    Matrix& algorithmicTangent)
  {
			 // Using perturbation methods
                         long double delta_strain =  0.00;
			 long double factor       =  1E-10;
                         long double max          =  1E-14;
                         double last_damage       =  md;
			 double last_r            =  mr_new;

			 
			 StrainVectorPerturbation.resize(6, false);
			 StressVectorPerturbation.resize(6, false);
			 //StrainVectorPerturbation_aux.resize(6, false);
			 //StressVectorPerturbation_aux.resize(6, false);
			 
			 noalias(StrainVectorPerturbation) = StrainVector;
			 //KRATOS_WATCH(StrainVector)
			 //StrainVectorPerturbation_aux = StrainVector;  			  

			 for (unsigned int i=0;i<StrainVectorPerturbation.size();i++)
			    {
				  
				 if  (fabs(StrainVectorPerturbation(i))<1E-15)
				    {
				     delta_strain = (*std::min_element(StrainVectorPerturbation.begin(),StrainVectorPerturbation.end()))*factor;
				     if (delta_strain==0.00)
					 {
					    delta_strain = factor;
					 }
				    } 
				 else
				    {
				      delta_strain = StrainVectorPerturbation(i)*factor;
				    }
					
				   if (delta_strain < max) {delta_strain=max;}


				    
				     StrainVectorPerturbation(i) += delta_strain;
				     
                                     //StrainVectorPerturbation_aux(i) -= delta_strain;
				    
				     Isotropic_Damage_3D::CalculateStress(StrainVectorPerturbation, StressVectorPerturbation);
                                     //Isotropic_Damage_3D::CalculateStress(StrainVectorPerturbation_aux, StressVectorPerturbation_aux);
				     
                                     // Antiguo procedimiento
			             noalias(StressVectorPerturbation) = StressVectorPerturbation-StressVector;
                                     
				     //noalias(StressVectorPerturbation) = StressVectorPerturbation-StressVectorPerturbation_aux;
				    
                                     noalias(StressVectorPerturbation) = StressVectorPerturbation/delta_strain;
				     //noalias(StressVectorPerturbation) = StressVectorPerturbation/(2.00*delta_strain);

				      	  
			
				    for (unsigned int j = 0; j<StrainVectorPerturbation.size(); j++)
					 {
					    algorithmicTangent(j,i) = StressVectorPerturbation(j); 
					 } 
				    //KRATOS_WATCH(algorithmicTangent)
				    md     = last_damage;
				    mr_new = last_r;
				    StressVectorPerturbation.resize(6, false);
				    StrainVectorPerturbation.resize(6, false);  
				    noalias(StrainVectorPerturbation) = StrainVector;
				    
				  }

				    
				 
		      }


}




