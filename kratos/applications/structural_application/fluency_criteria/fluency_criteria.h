/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

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
*   Last Modified by:    $Author: Nelson Lafontaine $
*   Date:                $Date: 26-06-2009 $
*   Revision:            $Revision: 1.0 $
*
* ***********************************************************/

#if !defined( KRATOS_FLUENCY_CRITERIA)
#define KRATOS_FLUENCY_CRITERIA

/* System includes */

/* External includes */
#include "boost/smart_ptr.hpp"

/* Project includes */
//#include "includes/define.h"
//#include "includes/variables.h"
//#include "utilities/math_utils.h"
//#include "custom_utilities/sd_math_utils.h"
#include "structural_application.h"
#include "custom_utilities/tensor_utils.h"
#include "includes/ublas_interface.h"
#include "includes/properties.h"
#include <cmath>
#include <string>
#include <iostream>



namespace Kratos
{

     enum myState{Plane_Stress, Plane_Strain, Tri_D };
     enum myPotencialPlastic{Not_Associated, Associated}; 
    
      class FluencyCriteria
	  {
        public:

	            double mSigma_e;          // Esfuerzo efectivo
	            double mSigma_y;          // Esfuerzo Resistencia de Comparacion
                    double mElasticDomain;
 
                    ///* Multisurface Platicity
                    Vector mMultisurface_Platicity_Yield;
                    Vector mMultisurface_Platicity_Sigma;  
                    //array_1d<double,6> mMorh_Coulomb_Yield;   
                    //array_1d<double,6> mMorh_Coulomb_Sigma_e; 
                    //array_1d<double,3> mRankine_Yield;   
                    //array_1d<double,3> mRankine_Sigma_e; 


		    typedef boost::numeric::ublas::vector<Vector> Second_Order_Tensor; // dos opciones: un tensor de segundo orden y/o un vector que almacena un vector		  
		    typedef boost::numeric::ublas::vector<Second_Order_Tensor> Third_Order_Tensor;

		    typedef boost::numeric::ublas::vector<boost::numeric::ublas::vector<Matrix> > Fourth_Order_Tensor;

		    typedef matrix<Second_Order_Tensor> Matrix_Second_Tensor; // Acumulo un tensor de 2 orden en una matriz.

		    typedef FluencyCriteria FluencyCriteriaType;

                    typedef Properties::Pointer PropertiesPointer;
	            
		    FluencyCriteria(){}

		    virtual ~FluencyCriteria(){}

                    KRATOS_CLASS_POINTER_DEFINITION( FluencyCriteria );

                     virtual boost::shared_ptr<FluencyCriteria> Clone() const
                      {
			    boost::shared_ptr<FluencyCriteria> p_clone(new FluencyCriteria ());
			    return p_clone;
		      }



		    virtual void InitializeMaterial(const Properties& props)
		      {
			  KRATOS_ERROR(std::logic_error,  "Called the virtual function for InitializeMaterial" , "");
		      }

		    virtual void CalculateEquivalentUniaxialStress(  
		    const Vector& StressVector, const Vector& StrainVector, const Matrix& Other, double& Result)
		    {
                        KRATOS_ERROR(std::logic_error,  "Called the virtual function for CalculateEquivalentUniaxialStress" , "");
		    }


		    virtual void CalculateEquivalentUniaxialStress(  
		    const Vector& StressVector,double& Result)
		    {
                        KRATOS_ERROR(std::logic_error,  "Called the virtual function for CalculateEquivalentUniaxialStress" , "");
		    }

		    virtual void CalculateEquivalentUniaxialStressViaPrincipalStress(  
		    const Vector& StressVector,double& Result)
		    {
	              KRATOS_ERROR(std::logic_error,  "Called the virtual function for CalculateEquivalentUniaxialStressViaPrincipalStress" , "");
		    }


		    virtual void CalculateEquivalentUniaxialStressViaInvariants( 
		    const Vector& StressVector,double& Result)

		    {
	              KRATOS_ERROR(std::logic_error,  "Called the virtual function for CalculateEquivalentUniaxialStressViaInvariants", "");
		    }

		    virtual void CalculateEquivalentUniaxialStressViaCilindricalCoordinate( 
		    const Vector& StressVector,double& Result)
		    {
	              KRATOS_ERROR(std::logic_error,  "Called the virtual function for CalculateEquivalentUniaxialStressViaCilindricalCoordinate", "");
		    }


		    virtual void CalculateEquivalentUniaxialStressMultiSurface(  
		    const Vector& StressVector,Vector& Result)
		    {
                        KRATOS_ERROR(std::logic_error,  "Called the virtual function for CalculateEquivalentUniaxialStressMultiSurface" , "");
		    }

                    

		    virtual void CalculateDerivateFluencyCriteria(const Vector& StressVector, Vector& DerivateFluencyCriteria)
		    {
	              KRATOS_ERROR(std::logic_error,  "Called the virtual function for CalculateDerivateFluencyCriteria", "");
		     }

		    virtual void CalculateDerivatePotencialFlowCriteria(const Vector& StressVector, Vector& DerivateFluencyCriteria)
		    {
	              KRATOS_ERROR(std::logic_error,  "Called the virtual function for CalculateDerivatePotencialFlowCriteria", "");
		     }

                    virtual void CalculateDerivateFluencyCriteriaMultiSurface(const Vector& StressVector, vector<Vector>& DerivateFluencyCriteria)
		    {
	              KRATOS_ERROR(std::logic_error,  "Called the virtual function for CalculateDerivateFluencyCriteria", "");
		    }

                    virtual void CalculateDerivatePotencialFlowCriteriaMultiSurface(const Vector& StressVector, vector<Vector>& DerivatePotencialFlow)
		    {
	              KRATOS_ERROR(std::logic_error,  "Called the virtual function for CalculateDerivateFluencyCriteria", "");
		    }
 

		    virtual void UpdateVariables( const Vector& Variables)
		    {
	              KRATOS_ERROR(std::logic_error,  "Called the virtual function for UpdateVariables", "");
		    }


		    /**
		    * Calculates the value of a given variable on the current integration point
		    * @param Trial_Stress = trial stress 
		    * @param delta_lamda  = The plastic Multiplicater 
		    * @param Result       = The corrected principal stress
		    */
                    
                    virtual void ReturnMapping(const Vector& StressVector, 
                    const Vector& StrainVector, 
                    Vector& delta_lamda,
                    array_1d<double,3>& Result)
		    {
	              KRATOS_ERROR(std::logic_error,  "Called the virtual function for ReturnMapping", "");
		    }

  
                    virtual void GetValue(double Result)
		    {
	              KRATOS_ERROR(std::logic_error,  "GetValue", "");
		    }

                    virtual void GetValue(Vector& Result)
		    {
	              KRATOS_ERROR(std::logic_error,  "GetValue", "");
		    }
                  

                    virtual void Finalize()
		     {
	              KRATOS_ERROR(std::logic_error,  "Finalize", "");
		     }  
 
                        
                        

		    //protected:


                   //
                    void CalculatePrincipalStressVector(const Vector& StressVector, array_1d<double,3>& Principal_Stress_Vector) 
                   {
		      int    iter      = 1000;
		      double zero      = 1.0E-12;
		      Matrix EigenVectors            = ZeroMatrix(3,3);
		      Matrix StressTensor            = ZeroMatrix(3,3);
		      Principal_Stress_Vector        = ZeroVector(3);
                      Vector Prin_Stress_Vector_Aux  = ZeroVector(3);
		      this->State_Tensor(StressVector,StressTensor);

		      SD_MathUtils<double>::EigenVectors(StressTensor, EigenVectors, Prin_Stress_Vector_Aux, zero, iter);

		      ///*sigma_1 >  sigma_2 > sigma_3 
		      sort (Prin_Stress_Vector_Aux.begin(), Prin_Stress_Vector_Aux.end()); 
		      reverse(Prin_Stress_Vector_Aux.begin(),Prin_Stress_Vector_Aux.end()); 
                      noalias(Principal_Stress_Vector) = Prin_Stress_Vector_Aux;   

                   }


		    void  Comprobate_State_Tensor(Matrix& StressTensor, const Vector& StressVector)
		    {
		      // Necesario para calcular eigen valores con subrutina de Jacobi, NO ACEPTA TERMINOS NULOS en diagonal principal.
		      if (fabs(StressTensor(0,0))<1E-10){StressTensor(0,0) = 1E-10; }
		      if (fabs(StressTensor(1,1))<1E-10){StressTensor(1,1) = 1E-10; }
		      if (fabs(StressTensor(2,2))<1E-10){StressTensor(2,2) = 1E-10; }       
		     }

                   void State_Tensor( const Vector& StressVector, Matrix& StressTensor)
                           {    
				
                                StressTensor.resize(3,3, false);
				switch (mState)
				  {
				    case Plane_Stress:
                                          {
					  StressTensor (0,0) = StressVector[0]; StressTensor (0,1) = StressVector[2]; StressTensor (0,2) = 0.00;
					  StressTensor (1,0) = StressVector[2]; StressTensor (1,1) = StressVector[1]; StressTensor (1,2) = 0.00;
					  StressTensor (2,0) = 0.00;            StressTensor (2,1) = 0.00;            StressTensor (2,2) = 1E-10;
					  break;
                                          }
				
				      case Plane_Strain:
                                          {
                                          double sigma_z = 0.00; 
					  sigma_z = (*mprops)[POISSON_RATIO]*(StressVector[0]+StressVector[1]);
					  StressTensor (0,0) = StressVector[0]; StressTensor (0,1) = StressVector[2]; StressTensor (0,2) = 0.00;
					  StressTensor (1,0) = StressVector[2]; StressTensor (1,1) = StressVector[1]; StressTensor (1,2) = 0.00;	
					  StressTensor (2,0) = 0.00;            StressTensor (2,1) = 0.00;            StressTensor (2,2) = sigma_z;
					  break;
                                          }
				      
				      case Tri_D:
					  StressTensor = MathUtils<double>::StressVectorToTensor(StressVector);
					  break;
                                      
                                      default:
				      std::cout<<"Warning: State not valid"<<std::endl;
				  }

                             Comprobate_State_Tensor(StressTensor, StressVector);
			  }


                     

                     void CalculateVectorFlowDerivate(const Vector& StressVector, Second_Order_Tensor& a)
                          {

			    matrix<double> StressTensor         = ZeroMatrix(3,3);
			    matrix<double> Aux_Tensor           = ZeroMatrix(3,3);
			    matrix<double> SphericComponent     = IdentityMatrix(3,3);
			    matrix<double> DesviatoricComponent = ZeroMatrix(3,3);

                            Vector  I(3);
                            Vector  J(3); 
                            Vector  J_des(3); 

			    a.resize(3, false);
			    a[0].resize(6, false); a[0] = ZeroVector(6);
			    a[1].resize(6, false); a[1] = ZeroVector(6);
			    a[2].resize(6, false); a[2] = ZeroVector(6);
                            
                            // Calculando a1
                            // Parcial(J1)/Parcial(sigma)
			    a[0](0)= 1.00;
			    a[0](1)= 1.00;
			    a[0](2)= 1.00;
			    a[0](3)= 0.00;
			    a[0](4)= 0.00;
			    a[0](5)= 0.00;
              
                            // calculando a2
                            State_Tensor(StressVector, StressTensor);
                            Tensor_Utils<double>::TensorialInvariants(StressTensor, I, J, J_des);
                            //KRATOS_WATCH(I)
                            //KRATOS_WATCH(StressVector)
                            //KRATOS_WATCH(StressTensor)


                           double factor_A =  1.00/(2.00*sqrt(J_des[1])); 
                           if (fabs(J_des[1])<=1E-8) 
			  {
			     mtetha_Lode = 0.00;
                             factor_A    = 0.00;                               
			  }
			  else
			  {  
			  mtetha_Lode = -(3.00*sqrt(3.00)*J_des(2))/(2.00*pow(J_des[1], 1.50));
			  if(fabs(mtetha_Lode) > 1.00){mtetha_Lode = 1.00; }
			  mtetha_Lode = asin(mtetha_Lode)/3.00; 
			  }



                            
	                    noalias(SphericComponent)     =  (I(0)/3.00)*SphericComponent;
	                    noalias(DesviatoricComponent) =  StressTensor - SphericComponent;   

			    a[1](0) = DesviatoricComponent(0,0);
			    a[1](1) = DesviatoricComponent(1,1);
			    a[1](2) = DesviatoricComponent(2,2);
			    a[1](3) = 2.00*StressTensor(0,1);       // xy
			    a[1](4) = 2.00*StressTensor(1,2);       // yz
			    a[1](5) = 2.00*StressTensor(0,2);       // xz
                            a[1]    =  factor_A*a[1];

                            // calculando a3
                            a[2](0) =  DesviatoricComponent(1,1)*DesviatoricComponent(2,2) - StressTensor(1,2)*StressTensor(1,2) + J_des[1]/3.00; 
                            a[2](1) =  DesviatoricComponent(0,0)*DesviatoricComponent(2,2) - StressTensor(0,2)*StressTensor(0,2) + J_des[1]/3.00;
                            a[2](2) =  DesviatoricComponent(0,0)*DesviatoricComponent(1,1) - StressTensor(0,1)*StressTensor(0,1) + J_des[1]/3.00;

                            a[2](3) =  2.00*(StressTensor(1,2)*StressTensor(0,2) - DesviatoricComponent(2,2)*StressTensor(0,1)); 
		            a[2](4) =  2.00*(StressTensor(0,1)*StressTensor(1,2) - DesviatoricComponent(1,1)*StressTensor(0,2));
                            a[2](5) =  2.00*(StressTensor(0,2)*StressTensor(0,1) - DesviatoricComponent(0,0)*StressTensor(1,2));    
                            //KRATOS_WATCH(a)

			  }  


                    void CalculateDerivatePrincipalStress(const Vector& StressVector, Second_Order_Tensor& DerivateStress)
                      {

                          // WARNING = For Revision
                          Second_Order_Tensor c;
                          Second_Order_Tensor a; 
			  c.resize(3, false);
                          c[0].resize(3, false); c[0] = ZeroVector(3);
                          c[1].resize(3, false); c[1] = ZeroVector(3);
                          c[2].resize(3, false); c[2] = ZeroVector(3);

                         CalculateVectorFlowDerivate(StressVector,  a);

                         c[0](0) = 1.00/3.00;  
                         c[0](1) = 1.00/3.00; 
                         c[0](2) = 1.00/3.00;
                    
                         if (fabs(mtetha_Lode) <= 0.50614548307836) // angulos comprendidos entre (-29 < tetha < +29  )
                         {    
 
                         double fact  = 0.00;
                         double denom = (sqrt(3.00) * cos(3.00 * mtetha_Lode)) ;
                         if(denom == 0.00) {fact = 1E14;}    
                         else {fact  = 2.00 / denom;}

                         c[1](0) = sin(mtetha_Lode - 3.00 * mtetha_Lode + 2.00 * PI / 3.00 );  
                         c[1](1) = sin(mtetha_Lode - 3.00 * mtetha_Lode ); 
                         c[1](2) = sin(mtetha_Lode - 3.00 * mtetha_Lode - 2.00 * PI / 3.00 );
                         noalias(c[1]) =  fact * (c[1]);

                        matrix<double> StressTensor = ZeroMatrix(3,3); 
                        Vector  I(3);
                        Vector  J(3); 
                        Vector  J_des(3); 
                        State_Tensor(StressVector, StressTensor);
                        Tensor_Utils<double>::TensorialInvariants(StressTensor, I, J, J_des);                        
                         
                        double fact2  = 0.00;
                        double denom2 = J_des[1] * cos(3.00 * cos(3.00 * mtetha_Lode)) ;
                        if(denom2 == 0.00) {fact2 = 1E14;}    
                        else {fact2  = 2.00 / denom2;}

                        c[2](0) = cos(mtetha_Lode  + 2.00 * PI / 3.00 );  
                        c[2](1) = cos(mtetha_Lode); 
                        c[2](2) = cos(mtetha_Lode - 2.00 * PI / 3.00 );
                        noalias(c[2]) =  fact2 * (c[2]);
                        }
                        else
                        {
                         double fact  = 2.00 / sqrt(3.00);
                         c[1](0) = sin(mtetha_Lode  + 2.00 * PI / 3.00 );  
                         c[1](1) = sin(mtetha_Lode); 
                         c[1](2) = sin(mtetha_Lode  - 2.00 * PI / 3.00 );
                         noalias(c[1]) =  fact * (c[1]); 

                        c[2](0) = 0.00;  
                        c[2](1) = 0.00; 
                        c[2](2) = 0.00;
                           
                        }

                       //KRATOS_WATCH(c)
                       //KRATOS_WATCH(a) 

                          DerivateStress.resize(3, false);
                          DerivateStress[0].resize(6, false); DerivateStress[0] = ZeroVector(6);
                          DerivateStress[1].resize(6, false); DerivateStress[1] = ZeroVector(6);
                          DerivateStress[2].resize(6, false); DerivateStress[2] = ZeroVector(6);

                        for(unsigned int i =0; i<3; i++)
                          {
                            for(unsigned int j =0; j<3; j++)
                            { 
                               noalias(DerivateStress[i]) = DerivateStress[i] +  c[j](i) * a[j]; 
                            }
                          }
                        //KRATOS_WATCH(DerivateStress)   

                      } 
                  
                        
		    public:
		    
		    const Properties *mprops;
                    myState mState;
                    myPotencialPlastic mPotencialPlastic;
                    double maccumulated_plastic_strain_current;   
                    double maccumulated_plastic_strain_old;

                    private:
                    double mtetha_Lode;
                    
                      

                 
    }; /* Class FluencyCriteria */
    

    /**
     * definition of CONSTITUTIVE_LAW variable
     */
}  /* namespace Kratos.*/
#endif /* FLUENCY_CRITERIA defined */
