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

// Project includes 
#include "includes/define.h"
#include "constitutive_laws/tension_compression_damage_model.h"
#include "includes/constitutive_law.h"
#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"
#include "custom_utilities/tensor_utils.h"
#include "includes/variables.h"
#include "includes/process_info.h"
#include "includes/properties.h"



namespace Kratos
{
    namespace Tension_Compression_Damage_Model_Auxiliaries
    {

    } 


    using namespace Tension_Compression_Damage_Model_Auxiliaries;

	/**
	 *	TO BE TESTED!!!
	 */
	 Tension_Compression_Damage_Model::Tension_Compression_Damage_Model() 
	: ConstitutiveLaw< Node<3> >()
	{
	}
	/**
	 *	TO BE TESTED!!!
	 */
	Tension_Compression_Damage_Model::~Tension_Compression_Damage_Model()
	{
	}
	
	
	bool Tension_Compression_Damage_Model::Has( const Variable<double>& rThisVariable )
	{
		return false;
	}
	
	bool Tension_Compression_Damage_Model::Has( const Variable<Vector>& rThisVariable )
	{
		return false;
	}
	
	bool Tension_Compression_Damage_Model::Has( const Variable<Matrix>& rThisVariable )
	{
		return false;
	}
	
	double Tension_Compression_Damage_Model::GetValue( const Variable<double>& rThisVariable )
	{
	    if( rThisVariable == DAMAGE)
	    {			
	    //KRATOS_WATCH(md)
	    return 0.00;

	    }
	    else return 0.00;
	    //KRATOS_ERROR(std::logic_error, "double Variable case not considered", "");
	}
	

	Vector Tension_Compression_Damage_Model::GetValue( const Variable<Vector>& rThisVariable )
	{
	    KRATOS_ERROR(std::logic_error, "Vector Variable case not considered", "");
	}
	
	Matrix Tension_Compression_Damage_Model::GetValue( const Variable<Matrix>& rThisVariable )
	{
	    KRATOS_ERROR(std::logic_error, "Matrix Variable case not considered", "");
	}

    void Tension_Compression_Damage_Model::SetValue( const Variable<double>& rThisVariable, const double& rValue, 
                                const ProcessInfo& rCurrentProcessInfo )
	{
	}
	
    void Tension_Compression_Damage_Model::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue, 
                                const ProcessInfo& rCurrentProcessInfo )
	{
	}
	
    void Tension_Compression_Damage_Model::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue, 
                                const ProcessInfo& rCurrentProcessInfo )
	{
	}
    
    void Tension_Compression_Damage_Model::Calculate(const Variable<Matrix >& rVariable, Matrix& rResult, 
                   const ProcessInfo& rCurrentProcessInfo)
    {
    }


//***********************************************************************************************
//***********************************************************************************************

void Tension_Compression_Damage_Model::InitializeMaterial( const Properties& props,
const GeometryType& geom,
const Vector& ShapeFunctionsValues )

      {

      msplit_negative_stress = ZeroVector(6);
      msplit_positive_stress = ZeroVector(6);
      mCurrentStress         = ZeroVector(6);
      mCtangent              = ZeroMatrix(6,6);
      mInSituStress          = ZeroVector(6);
      
      double angle = PI*props[FRICTION_INTERNAL_ANGLE]/180.00;

      malfa      = 2.00*sin(angle)/(sqrt(3)*(3 + sin(angle))); 
      mFc        = props[FC];
      mFt        = props[FT];
      mEc        = props[CONCRETE_YOUNG_MODULUS_C];
      mEt        = props[CONCRETE_YOUNG_MODULUS_T];
      mNU        = props[POISSON_RATIO];
      mGE        = props[FRACTURE_ENERGY];
      ml         = pow(fabs(geom.Volume()),0.333333333333333333); 
      mr_neg_old =  props[FT];
      mr_pos_old = (1.00-malfa)*props[FC];
      

      CalculateNoDamageElasticMatrix(mCtangent, props[YOUNG_MODULUS], props[POISSON_RATIO]);      
      }

		

//***********************************************************************************************
//***********************************************************************************************

void Tension_Compression_Damage_Model::InitializeSolutionStep( const Properties& props,
const GeometryType& geom,
const Vector& ShapeFunctionsValues ,
const ProcessInfo& CurrentProcessInfo)
{
}

//***********************************************************************************************
//***********************************************************************************************


void Tension_Compression_Damage_Model::CalculateNoDamageElasticMatrix(Matrix& C, const double E, const double NU)
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


//***********************************************************************************************
//***********************************************************************************************


void Tension_Compression_Damage_Model::CalculateConstitutiveMatrix(const Vector& StrainVector, Matrix& ConstitutiveMatrix)
      {
      
	  noalias(ConstitutiveMatrix)= mCtangent;
      }


//***********************************************************************************************
//***********************************************************************************************


 void Tension_Compression_Damage_Model::CalculateDamage(const Vector& StressVector)
 { 


    double B_pos = 0.00;
    double B_neg = 0.00;

    Vector scalar_function = ZeroVector(2);
    Vector temp  = ZeroVector(6);
    Vector I     = ZeroVector(3);
    Vector J     = ZeroVector(3);
    Vector J_des = ZeroVector(3);
    
    Matrix invertedmatrix = ZeroMatrix(6,6); 
    
   
    SD_MathUtils<double>::InvertMatrix(mCtangent, invertedmatrix);
    temp = prod(invertedmatrix,msplit_positive_stress);
    scalar_function(0) = sqrt(mEc*(inner_prod(temp,StressVector))); 
    
    
    Tensor_Utils<double>::TensorialInvariants(StressVector,I,J,J_des);
    scalar_function(1) = malfa*I(0) + sqrt(3.00*J_des(1)); 
    
    
    mr_pos_new = std::max(mr_pos_old,scalar_function(0));
    mr_neg_new = std::max(mr_neg_old,scalar_function(1));
    
    // Calculando variables de dano negativas y positivas
   // dano a traccion
   B_pos      =  1.00/((mGE*mEc)/(ml*mFt*mFt)-0.5);
   B_pos      = B_pos*(1.00-mr_pos_new/mFt);
   md_pos = 1.00 - (mFt/mr_pos_new)*exp(B_pos);
   if (md_pos < 0.00) 
   {
     md_pos = fabs(md_pos);
     //std::cout<<"Warning: Damage is less than zero"<<std::endl;
   }
    
    //  No se especifica en articulo y libro matematico de B_neg.
    //  Esta sujeto a ser encontrado por datos experimentales.

    
    
 }



//***********************************************************************************************
//***********************************************************************************************


void  Tension_Compression_Damage_Model::FinalizeSolutionStep( const Properties& props,
                                               const GeometryType& geom, 
                                               const Vector& ShapeFunctionsValues ,
                                             const ProcessInfo& CurrentProcessInfo)

			{
		          
			}


//***********************************************************************************************
//***********************************************************************************************

	void Tension_Compression_Damage_Model::CalculateNoDamageStress(const Vector& StrainVector, Vector& StressVector)
	{

		double zero_tolerance  = 1E-10;
		unsigned int iter      = 100;
		Matrix EigenVectors;
		IdentityMatrix I(6,6);
		Second_Order_Tensor EigenVector(3);
		Vector EigenValues;
		Matrix Stress_Tensor;
		
		EigenValues.resize(3,false);
		Stress_Tensor.resize(3,3,false);
		EigenVectors.resize(3,3,false);
		msplit_positive_stress.resize(6,false);
		msplit_negative_stress.resize(6,false);
		EigenVector[0].resize(3,false);
		EigenVector[1].resize(3,false);
		EigenVector[2].resize(3,false);
		
		if( StressVector.size() != 6 )
		{
			StressVector.resize(6);
		}
		
		
		noalias(StressVector) = prod(mCtangent,StrainVector);
		mCurrentStress = StressVector;
            	
		Stress_Tensor = MathUtils<double>::StressVectorToTensor(StressVector);
		Verify_Matrix(Stress_Tensor);
		SD_MathUtils<double>::EigenVectors(Stress_Tensor, EigenVectors, EigenValues, zero_tolerance, iter); 
		EigenVector[0](0) = EigenVectors(0,0); EigenVector[0](1) = EigenVectors(0,1); EigenVector[0](2) = EigenVectors(0,2);
		EigenVector[1](0) = EigenVectors(1,0); EigenVector[1](1) = EigenVectors(1,1); EigenVector[1](2) = EigenVectors(1,2);
		EigenVector[2](0) = EigenVectors(2,0); EigenVector[2](1) = EigenVectors(2,1); EigenVector[2](2) = EigenVectors(2,2);
		
    
		//double brakets = 0.00;
		Matrix temp(3,3);
		noalias(temp) = ZeroMatrix(3,3);
                Second_Order_Tensor p(3);
		p[0].resize(3);      
		p[1].resize(3);
		p[2].resize(3); 

		Fourth_Order_Tensor P_Tensor;
                Matrix P_Matrix = ZeroMatrix(6,6);
                Matrix temp_a = ZeroMatrix(6,6); 
                Matrix temp_b = ZeroMatrix(3,3);    

		// Forma simplicada para expresar ecuacion 3a de "An energy release rate-based plastic-damage model for concrete".
		// Trabajo todo en un tensor de cuart orden. 
//  		for (unsigned int i=0;i<3;i++)
//  		    {
//  		        brakets = EigenValues(i);
//  			if (brakets < 0.00) {brakets = 0.00;}
//  			noalias(temp) = brakets*outer_prod(EigenVector[i],EigenVector[i]) + temp; 
//  		    }	
 	                
		Compute_pij_Tensor(EigenVector,EigenVector);
		double heaviside = 1.00;
		for(unsigned int i = 0; i<3; i++)
		  { 
		     heaviside = 1.00; 
		     if (EigenValues(i) < 0.00) {heaviside = 0.00;}	  
		     (Tensor_Utils<double>::Prod_Second_Order_Tensor(mp(i,i),mp(i,i),P_Tensor));
                     SD_MathUtils<double>::TensorToMatrix(P_Tensor, temp_a);
                     noalias(P_Matrix) =  P_Matrix + heaviside*temp_a;
                     temp_a = ZeroMatrix(6,6);     
		   } 
		   
		   //Traccion es siempre positiva 
                   SD_MathUtils<double>::TensorToMatrix(P_Tensor, temp_a);
		   noalias(msplit_positive_stress) = prod(P_Matrix,StressVector);
                   noalias(msplit_negative_stress) = StressVector - msplit_positive_stress;  
	}
	

//***********************************************************************************************
//***********************************************************************************************
	

 void Tension_Compression_Damage_Model::CalculateStress( const Vector& StrainVector, Vector& StressVector)
            
	{
	CalculateNoDamageStress(StrainVector, StressVector);	
	CalculateDamage(StressVector);
        }


//***********************************************************************************************
//***********************************************************************************************


 void Tension_Compression_Damage_Model::Compute_Principal_Stress(const Vector& StressVector, Vector& Result)
            
	{
	    
	   double tau_max =  sqrt(((StressVector(0) -StressVector(1))/2.0)*((StressVector(0) -StressVector(1))/2.0) + StressVector(2)*StressVector(2)); 
	   double ten_med =  0.500*(StressVector(0)+ StressVector(1));
	   Result(0) = ten_med + tau_max;  
	   Result(1) = ten_med - tau_max; 									
        }


//***********************************************************************************************
//***********************************************************************************************

  void Tension_Compression_Damage_Model::CalculateCauchyStresses(
		Vector& rCauchy_StressVector,
		const Matrix& rF,
		const Vector& rPK2_StressVector,
		const Vector& rGreenLagrangeStrainVector)
    {		
    
    }

//***********************************************************************************************
//***********************************************************************************************	

void Tension_Compression_Damage_Model::CalculateStressAndTangentMatrix(Vector& StressVector,
                    const Vector& StrainVector,
                    Matrix& algorithmicTangent)
   {
              noalias(algorithmicTangent)= mCtangent;
   }
   
   //***********************************************************************************************
//***********************************************************************************************

 void Tension_Compression_Damage_Model::Verify_Matrix(Matrix& Result)
 
 {
   for (unsigned int i=0; i<Result.size1(); i++)
   {
    for (unsigned int j=0; j<Result.size2();j++)
       {
         if (fabs(Result(i,j))<1E-15){Result(i,j) = 1E-15;} 
       }
   }
    
 }
 
 void Tension_Compression_Damage_Model::Compute_pij_Tensor(const Second_Order_Tensor& A, const Second_Order_Tensor& B)
 
 {
    Matrix temp(3,3);
    Second_Order_Tensor Result(3);
    Result[0].resize(3,false);
    Result[1].resize(3,false);
    Result[2].resize(3,false);
     
    mp.resize(3,3);
 	 for (unsigned int i=0; i<3; i++){
 	      for (unsigned int j=0; j<3; j++){
 		      mp(i,j).resize(3);
		      for (unsigned int k=0; k<3; k++){
 	                mp(i,j)[k].resize(3,false);}
 		    }
	      }
       
    for (unsigned int i=0; i<3; i++){
       for (unsigned int j=0; j<3; j++){
	  //temp = 0.50*(outer_prod(A(i),B(j)) +outer_prod(A(j),B(i))); 
          temp = outer_prod(A(i),B(j));  
	  Result[0](0) = temp(0,0); Result[0](1) = temp(0,1); Result[0](2) = temp(0,2);
	  Result[1](0) = temp(1,0); Result[1](1) = temp(1,1); Result[1](2) = temp(1,2);
	  Result[2](0) = temp(2,0); Result[2](1) = temp(2,1); Result[2](2) = temp(2,2);
          mp(i,j) = Result;}
        }

}
}