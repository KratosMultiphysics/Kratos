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

	            double mSigma_e;  // Esfuerzo efectivo
	            double mSigma_y;  // Esfuerzo Resistencia de Comparacion
                    double mElasticDomain;


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


		    virtual void CalculateDerivateFluencyCriteria(const Vector& StressVector, Vector& DerivateFluencyCriteria)
		    {
	              KRATOS_ERROR(std::logic_error,  "Called the virtual function for CalculateDerivateFluencyCriteria", "");
		     }

		    virtual void CalculateDerivatePotencialFlowCriteria(const Vector& StressVector, Vector& DerivateFluencyCriteria)
		    {
	              KRATOS_ERROR(std::logic_error,  "Called the virtual function for CalculateDerivatePotencialFlowCriteria", "");
		     }

		    virtual void UpdateVariables( const Vector& Variables)
		    {
	              KRATOS_ERROR(std::logic_error,  "Called the virtual function for UpdateVariables", "");
		     }

		    //protected:

		    void  Comprobate_State_Tensor(Matrix& StressTensor, const Vector& StressVector)
		    {
		      // Necesario para calcular eigen valores con subrutina de Jacobi, NO ACEPTA TERMINOS NULOS en diagonal principal.
		      if (fabs(StressTensor(0,0))<1E-10){StressTensor(0,0) = 1E-10; }
		      if (fabs(StressTensor(1,1))<1E-10){StressTensor(1,1) = 1E-10; }
		      if (fabs(StressTensor(2,2))<1E-10){StressTensor(2,2) = 1E-10; }       
		     }

                   void State_Tensor( const Vector& StressVector, Matrix& StressTensor)
                           {

				    
				double sigma_z = 0.00;
                                StressTensor.resize(3,3, false);
				switch (mState)
				  {
				    case Plane_Stress:
                                          {
					  StressTensor (0,0) = StressVector(0); StressTensor (0,1) = StressVector(2); StressTensor (0,2) = 0.00;
					  StressTensor (1,0) = StressVector(2); StressTensor (1,1) = StressVector(1); StressTensor (1,2) = 0.00;
					  StressTensor (2,0) = 0.00;            StressTensor (2,1) = 0.00;            StressTensor (2,2) = 1E-10;
                                          
					  break;
                                          }
				
				      case Plane_Strain:
					  sigma_z = (*mprops)[POISSON_RATIO]*(StressVector(0)+StressVector(1));
					    
					  {
					  StressTensor (0,0) = StressVector(0); StressTensor (0,1) = StressVector(2); StressTensor (0,2) = 0.00;
					  StressTensor (1,0) = StressVector(2); StressTensor (1,1) = StressVector(1); StressTensor (1,2) = 0.00;	
					  StressTensor (2,0) = 0.00;            StressTensor (2,1) = 0.00;            StressTensor (2,2) =sigma_z;
                                          }
					  break;
				      
				      case Tri_D:
					  StressTensor = MathUtils<double>::StressVectorToTensor(StressVector);
					  break;
                                      
                                      default:
				      std::cout<<"Warning: State not valid"<<std::endl;
				  }
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
                            //KRATOS_WATCH(J_des)


                            double factor_A  =  (1.00/(2.00*sqrt(J_des(1)))); 
	                    noalias(SphericComponent)     =  (I(0)/3.00)*SphericComponent;
	                    noalias(DesviatoricComponent) =  StressTensor - SphericComponent;
			   

			    a[1](0) = DesviatoricComponent(0,0);
			    a[1](1) = DesviatoricComponent(1,1);
			    a[1](2) = DesviatoricComponent(2,2);
			    a[1](3) = 2.00*StressTensor(0,1);
			    a[1](4) = 2.00*StressTensor(1,2);
			    a[1](5) = 2.00*StressTensor(0,2);
                            a[1]    =  factor_A*a[1];

                            // calculando a3
                            a[2](0) =  DesviatoricComponent(1,1)*DesviatoricComponent(2,2) - StressTensor(1,2)*StressTensor(1,2) + J_des(1)/3.00; 
                            a[2](1) =  DesviatoricComponent(0,0)*DesviatoricComponent(2,2) - StressTensor(0,2)*StressTensor(0,2) + J_des(1)/3.00;
                            a[2](2) =  DesviatoricComponent(0,0)*DesviatoricComponent(1,1) - StressTensor(0,1)*StressTensor(0,1) + J_des(1)/3.00;
                            a[2](3) =  2.00*(StressTensor(1,2)*StressTensor(0,2) - DesviatoricComponent(2,2)*StressTensor(0,1)); 
                            a[2](4) =  2.00*(StressTensor(0,2)*StressTensor(0,1) - DesviatoricComponent(0,0)*StressTensor(1,2)); 
                            a[2](5) =  2.00*(StressTensor(0,1)*StressTensor(1,2) - DesviatoricComponent(1,1)*StressTensor(0,2));

			  }  
                        
		    public:
		    
		    const Properties *mprops;
                    myState mState;
                    myPotencialPlastic mPotencialPlastic;

                      

                 
    }; /* Class FluencyCriteria */
    

    /**
     * definition of CONSTITUTIVE_LAW variable
     */
}  /* namespace Kratos.*/
#endif /* FLUENCY_CRITERIA defined */
