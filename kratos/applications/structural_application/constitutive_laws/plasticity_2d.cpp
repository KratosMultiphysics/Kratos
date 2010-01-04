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
#include "constitutive_laws/plasticity_2d.h"
#include "includes/constitutive_law.h"
#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"
#include "custom_utilities/tensor_utils.h"
#include "includes/variables.h"
#include "includes/process_info.h"
#include "includes/properties.h"


namespace Kratos
{
    namespace Plasticity2D_Auxiliaries
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


    using namespace Plasticity2D_Auxiliaries;

	/**
	 *	TO BE TESTED!!!
	 */

         Plasticity2D::Plasticity2D() 
	: ConstitutiveLaw< Node<3> >()
	
	{
	  KRATOS_ERROR(std::logic_error,"Calling the empty constructor.","");
	}

	 Plasticity2D::Plasticity2D(FluencyCriteriaPointer FluencyCriteria, SofteningHardeningCriteriaPointer SofteningBehavior, PropertiesPointer Property) 
	: ConstitutiveLaw< Node<3> >()
	{
	      mpFluencyCriteria   = FluencyCriteria;
              mpSofteningBehavior = SofteningBehavior;
              mpProperties        = Property;

             // mFluencyCriteria = GetProperties()[CONSTITUTIVE_LAW]->Clone();
	      //KRATOS_WATCH(*mProperties)
	}
	/**
	 *	TO BE TESTED!!!
	 */
	Plasticity2D::~Plasticity2D()
	{
	}
	
	
	bool Plasticity2D::Has( const Variable<double>& rThisVariable )
	{
		return false;
	}
	
	bool Plasticity2D::Has( const Variable<Vector>& rThisVariable )
	{
		return false;
	}
	
	bool Plasticity2D::Has( const Variable<Matrix>& rThisVariable )
	{
		return false;
	}
	
	double Plasticity2D::GetValue( const Variable<double>& rThisVariable )
	{
	  if( rThisVariable == DAMAGE)
	  {			
	  return  mlamda;
	  }
	  else return 0.00;
	}   
	

	Vector Plasticity2D::GetValue( const Variable<Vector>& rThisVariable )
	{
//             if(rThisVariable==PLASTIC_STRAIN_TENSOR)
//             {
//                return  mplastic_strain;
//             }
//             else
            { 
	    KRATOS_ERROR(std::logic_error, "Vector Variable case not considered", "");
	    }
        }
	
	Matrix Plasticity2D::GetValue( const Variable<Matrix>& rThisVariable )
	{
	    KRATOS_ERROR(std::logic_error, "Matrix Variable case not considered", "");
	}

    void Plasticity2D::SetValue( const Variable<double>& rThisVariable, const double& rValue, 
                                const ProcessInfo& rCurrentProcessInfo )
	{
	}
	
    void Plasticity2D::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue, 
                                const ProcessInfo& rCurrentProcessInfo )
	{
	}
	
    void Plasticity2D::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue, 
                                const ProcessInfo& rCurrentProcessInfo )
	{
	}
    
    void Plasticity2D::Calculate(const Variable<Matrix >& rVariable, Matrix& rResult, 
                   const ProcessInfo& rCurrentProcessInfo)
    {
      for (unsigned int ii = 0; ii<3; ii++)
           rResult(0,ii) = mplastic_strain(ii); 
    }


//***********************************************************************************************
//***********************************************************************************************

	void Plasticity2D::InitializeMaterial( const Properties& props,
	const GeometryType& geom,
	const Vector& ShapeFunctionsValues )

{
	mE     = (*mpProperties)[YOUNG_MODULUS];
	mNU    = (*mpProperties)[POISSON_RATIO];
	mDE    = (*mpProperties)[DENSITY];
        
	mpFluencyCriteria->InitializeMaterial(*mpProperties);


        noalias(mplastic_strain)         = ZeroVector(4);
        noalias(mcurrent_plastic_strain) = ZeroVector(4);  
        noalias(mDerivate_Fluency)       = ZeroVector(4);     
}


//***********************************************************************************************
//***********************************************************************************************

void Plasticity2D::InitializeSolutionStep( const Properties& props,
                    const GeometryType& geom,
                    const Vector& ShapeFunctionsValues ,
                    const ProcessInfo& CurrentProcessInfo)
{
    mlamda                  = 0.00;   
}

//***********************************************************************************************
//***********************************************************************************************


void Plasticity2D::CalculateElasticMatrix(Matrix& C, const double E, const double NU)
	{ 

                double c1 = 0.00;
		double c2 = 0.00;
		double c3 = 0.00;

		c1 = E*(1.00-NU)/((1.00 + NU)*(1.00-2.00*NU));
		c2 = c1*NU/(1.00 - NU);
		c3 = 0.5*E / (1 + NU);
		C(0,0) = c1;  C(0,1) = c2;	C(0,2) = 0.0;  C(0,3)  = c2;
		C(1,0) = c2;  C(1,1) = c1;	C(1,2) = 0.0;  C(1,3)  = c2;
		C(2,0) = 0.0; C(2,1) = 0.0;	C(2,2) = c3;   C(2,3)  = 0.00;
		C(3,0) = c2;  C(3,1) = c2;	C(3,2) = 0.00; C(3,3)  = c1;
		
	}
	

//***********************************************************************************************
//***********************************************************************************************

void Plasticity2D::CalculateConstitutiveMatrix(const Vector& StrainVector, Matrix& ConstitutiveMatrix)
      {
	   /*   Matrix Aux(4,4);
	      CalculateElasticMatrix(Aux,mE,mNU); 
	      for(unsigned int i = 0; i<3; i ++){
	      for(unsigned int j= 0;j<3; j++)
	      {
	      ConstitutiveMatrix(i,j) = Aux(i,j);  
		    }
	    }
     */
       Vector Aux(3);
       CalculateStressAndTangentMatrix(Aux, StrainVector , ConstitutiveMatrix);
      }



//***********************************************************************************************
//***********************************************************************************************


void  Plasticity2D::FinalizeSolutionStep( const Properties& props,
                                               const GeometryType& geom, 
                                               const Vector& ShapeFunctionsValues ,
                                             const ProcessInfo& CurrentProcessInfo)

			{
                              mefective_plastic_strain =  mcurrent_efective_plastic_strain;
                              noalias(mplastic_strain) =  mcurrent_plastic_strain;
                              Vector Variables(1);
                              Variables(0) = mefective_plastic_strain;
                              mpFluencyCriteria->UpdateVariables(Variables);
			}


//***********************************************************************************************
//***********************************************************************************************

void Plasticity2D::CalculateElasticStress(const Vector& StrainVector, array_1d<double,4>& StressVector)
	{
		
		/*
		double c1 = mE / (1.00 - mNU*mNU);
		double c2 = c1 * mNU;
		double c3 = 0.5*mE / (1.00 + mNU);
	   		  
		StressVector[0] = c1*StrainVector[0] + c2 * (StrainVector[1])	;
		StressVector[1] = c1*StrainVector[1] + c2 * (StrainVector[0])	;
		StressVector[2] = c3*StrainVector[2];
                 */
                Matrix ElasticMatrix(4,4);
                CalculateElasticMatrix(ElasticMatrix, mE, mNU);
                noalias(StressVector) = prod(ElasticMatrix, StrainVector);
	}
	

//***********************************************************************************************
//***********************************************************************************************
	

 void Plasticity2D::CalculateStress( const Vector& StrainVector, Vector& StressVector)
            
	{
        
	KRATOS_TRY  
       
	double ElasticDomain           = 0.00;
        double delta_lamda             = 0.00;
        //double lamda                   = 0.00;

        array_1d<double,4>  StrainVector_aux;                //StrainVector_aux.resize(4, false); 
        array_1d<double,4>  StressVector_aux;         // StressVector_aux.resize(4, false);         
        array_1d<double,4>  ElasticStrain;             //ElasticStrain.resize(4, false);            
        array_1d<double,4>  Aux_1;                     //Aux_1.resize(4, false);               
        Vector DerivateFluencyCriteria;   DerivateFluencyCriteria.resize(6,false);      
        //Vector DerivatePotencialCriteria; DerivatePotencialCriteria.resize(6,false);           
	Matrix ElasticMatrix;             ElasticMatrix.resize(4,4, false);            
        //KRATOS_WATCH(mcurrent_plastic_strain)
        //         
         //calculating elastic strain
        StrainVector_aux[0] =  StrainVector[0];
        StrainVector_aux[1] =  StrainVector[1];
        StrainVector_aux[2] =  StrainVector[2];
        StrainVector_aux[3] =  -mNU/(1.00 - mNU)*(StrainVector_aux[0] + StrainVector_aux[1]) 
	+ (1.00/(1-mNU))*( mNU*(mcurrent_plastic_strain[0] + mcurrent_plastic_strain[1] - mcurrent_plastic_strain[3]) + mcurrent_plastic_strain[3] );
        
        
        noalias(ElasticStrain) = StrainVector_aux - mcurrent_plastic_strain;

        CalculateElasticStress(ElasticStrain, StressVector_aux); 
         
        // Tension elastica de prueba.
        StressVector[0] = StressVector_aux[0]; // Gxx
        StressVector[1] = StressVector_aux[1]; // Gyy
        StressVector[2] = StressVector_aux[2]; // Txy

       
        mpFluencyCriteria->CalculateEquivalentUniaxialStressViaInvariants(StressVector,ElasticDomain);
       
       if(ElasticDomain < 0.00)
          {
	     mlamda                           = 0.00;
             mcurrent_efective_plastic_strain = mefective_plastic_strain;
             noalias(mcurrent_plastic_strain) = mplastic_strain;
             ComputeTangentMatrix = false;
          }
       else
          {    
          ComputeTangentMatrix = true;
          
          double toler          =  1E-9;
          double aux_b          =  0.00;
          unsigned int max_iter =  1000;
          unsigned int iter     =  1;
          noalias(mDerivate_Fluency)   = ZeroVector(4);                
          CalculateElasticMatrix(ElasticMatrix, mE, mNU);
	  Vector Variables(1);
          Variables[0] = 0.00;

         // std::cout<<"Integrando Ecuacion Constitutiva"<<std::endl;
         // KRATOS_WATCH(Id())
          //std::cout<< "iteracion = "<< iter << " ** "<< "Elastic Domain = "<<ElasticDomain<<std::endl;
          while(fabs(ElasticDomain)>toler && iter<=max_iter)
              {

              mpFluencyCriteria->CalculateDerivateFluencyCriteria(StressVector, DerivateFluencyCriteria);
	      mDerivate_Fluency[0]      = DerivateFluencyCriteria[0];    
	      mDerivate_Fluency[1]      = DerivateFluencyCriteria[1];
	      mDerivate_Fluency[2]      = DerivateFluencyCriteria[3];
              mDerivate_Fluency[3]      = DerivateFluencyCriteria[2]; //deformacion plsatica ez
                        
              noalias(Aux_1) = prod(ElasticMatrix, mDerivate_Fluency);
              
              // Flujo no asociado
             switch (mpFluencyCriteria->mPotencialPlastic)
             { 
              case(Not_Associated):
              {
              std::cout<<"No potencial function has been added for plasticyty "<<std::endl;
              break;
              }

               case(Associated):
              {      
	      aux_b  = inner_prod(mDerivate_Fluency, Aux_1) + (*mpProperties)[PLASTIC_MODULUS];
              delta_lamda    = ElasticDomain/aux_b; 
              mlamda += delta_lamda;
	      noalias(mcurrent_plastic_strain)  +=  delta_lamda*mDerivate_Fluency;
              mcurrent_efective_plastic_strain  +=  sqrt((2.00/3.00)*inner_prod(delta_lamda*mDerivate_Fluency, delta_lamda*mDerivate_Fluency));
             
              //update stress
	      StrainVector_aux[0] =  StrainVector[0];
	      StrainVector_aux[1] =  StrainVector[1];
	      StrainVector_aux[2] =  StrainVector[2];
	      StrainVector_aux[3] =  -mNU/(1.00 - mNU)*(StrainVector_aux[0] + StrainVector_aux[1]) 
	      + (1.00/(1-mNU))*( mNU*(mcurrent_plastic_strain[0] + mcurrent_plastic_strain[1] - mcurrent_plastic_strain[3]) + mcurrent_plastic_strain[3]);
	
              noalias(ElasticStrain) = StrainVector_aux - mcurrent_plastic_strain;
              CalculateElasticStress(ElasticStrain, StressVector_aux);  // Tension elastica de prueba.
                    
              StressVector[0] = StressVector_aux[0];
              StressVector[1] = StressVector_aux[1];
              StressVector[2] = StressVector_aux[2];
           
              // updating internal variables
              //Updated_Internal_Variables(StressVector_aux, delta_lamda*mDerivate_Fluency);
              Variables[0] = mcurrent_efective_plastic_strain;
              mpFluencyCriteria->UpdateVariables(Variables);
              break;      
              }

              }
              mpFluencyCriteria->CalculateEquivalentUniaxialStressViaInvariants(StressVector, ElasticDomain);
              //std::cout<< "iteracion = "<< iter << " ** "<< "Elastic Domain = "<<ElasticDomain<<std::endl;
              if(ElasticDomain<0.00){ElasticDomain = 0.00;}               
              iter = iter + 1;
              
              if(iter>max_iter){std::cout<<" It has`nt reached the required convergence"<<std::endl;}
              }
        }        

    KRATOS_CATCH("")       
}


//***********************************************************************************************
//***********************************************************************************************	

  void Plasticity2D::CalculateCauchyStresses(
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
	 void Plasticity2D::CalculateStressAndTangentMatrix(Vector& StressVector,
                    const Vector& StrainVector,
                    Matrix& algorithmicTangent)
{
                Matrix ElasticMatrix(4,4);
                CalculateElasticMatrix(ElasticMatrix,mE,mNU);
                if(ComputeTangentMatrix==true)
                 {
                      // Inconsistent tangent Modular Matrix

                      double aux_b = 0.00;
                      Vector Aux_1(4);
                      Vector Aux_2(4);
                      Vector DerivateFluencyCriteria(6);
                      
                      mpFluencyCriteria->CalculateDerivateFluencyCriteria(StressVector, DerivateFluencyCriteria);
		      mDerivate_Fluency[0]      = DerivateFluencyCriteria[0];    
		      mDerivate_Fluency[1]      = DerivateFluencyCriteria[1];
		      mDerivate_Fluency[2]      = DerivateFluencyCriteria[3];
		      mDerivate_Fluency[3]      = DerivateFluencyCriteria[2]; //deformacion plsatica ez
			  

                      Aux_1 = prod(ElasticMatrix, mDerivate_Fluency);
                      
                      switch (mpFluencyCriteria->mPotencialPlastic) 
                       {   
                         case(Associated):
                          {                    
                          aux_b = inner_prod(mDerivate_Fluency, Aux_1) +  (*mpProperties)[PLASTIC_MODULUS];
                          noalias(ElasticMatrix) = ElasticMatrix - (outer_prod(Aux_1,Aux_1))/aux_b;
                          break;
                          }
                      case(Not_Associated): 
                        {         
                          break;         
                        }
                       }
                 }

                switch (mpFluencyCriteria->mState)
		  {
		    case Plane_Stress:
                      {
                      ComputeCondentationMatrix(ElasticMatrix,algorithmicTangent);              
                      break;
                      }
                    case Plane_Strain:
                      {
                      for(unsigned int i = 0; i<3; i ++)
		      {
                            for(unsigned int j= 0;j<3; j++)
                               {
                                  algorithmicTangent(i,j) = ElasticMatrix(i,j);  
                               }
                       }
                      }
                    }
}
                        
                         


//***********************************************************************************************
//***********************************************************************************************	

void Plasticity2D::UpdateMaterial( const Vector& StrainVector,
                                     const Properties& props,
                                     const GeometryType& geom,
                                     const Vector& ShapeFunctionsValues,
                                     const ProcessInfo& CurrentProcessInfo)
{
        mcurrent_efective_plastic_strain = mefective_plastic_strain;
        noalias(mcurrent_plastic_strain) = mplastic_strain;
        //if (CurrentProcessInfo[NL_ITERATION_NUMBER] != 1.00){Iniatilize = true;} 

       
}


//***********************************************************************************************
//***********************************************************************************************	


void Plasticity2D::Calculate(const Variable<double>& rVariable, 
                                    double& Output, 
                                    const ProcessInfo& rCurrentProcessInfo)
   {
   }

//***********************************************************************************************
//***********************************************************************************************	

void  Plasticity2D::Comprobate_State_Vector(Vector& Result)
		    {
                          for (unsigned int i = 0.00; i<Result.size(); i++){ 
                             if(fabs(Result(i))< 1E-9){
				    Result(i) = 0.00;}
		          } 
                    }


void Plasticity2D::ComputeCondentationMatrix(const Matrix& TangentMatrix, Matrix& Result)
{
     Matrix k11(3,3);
     Vector k12(3);
     Vector k21(3);
     Matrix k22(1,1);

     Result.resize(3,3, false);
     
     for(unsigned int i = 0; i<3; i++){
	      k12[i] = TangentMatrix(i,3);
              k21[i] = TangentMatrix(3,i);
          for(unsigned int j = 0; j<3; j++){
              k11(i,j)  = TangentMatrix(i,j);
             }
      }

    k22(0,0) =  TangentMatrix(3,3);
    noalias(Result) = k11  - outer_prod(k12,k21)/k22(0,0); 

}

}


