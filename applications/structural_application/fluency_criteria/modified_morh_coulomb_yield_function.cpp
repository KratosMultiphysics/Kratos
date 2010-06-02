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


//#include "utilities/math_utils.h"
//#include "custom_utilities/sd_math_utils.h"
//#include "custom_utilities/tensor_utils.h"
#include "fluency_criteria/modified_morh_coulomb_yield_function.h"
//#include <cmath>



namespace Kratos
  {

  
	    typedef boost::numeric::ublas::vector<Vector> Second_Order_Tensor; // dos opciones: un tensor de segundo orden y/o un vector que almacena un vector		  
	    typedef boost::numeric::ublas::vector<Second_Order_Tensor> Third_Order_Tensor;
			  
            typedef boost::numeric::ublas::vector<boost::numeric::ublas::vector<Matrix> > Fourth_Order_Tensor;
			  
	    typedef matrix<Second_Order_Tensor> Matrix_Second_Tensor; // Acumulo un tensor de 2 orden en una matri
           
            Modified_Morh_Coulomb_Yield_Function::Modified_Morh_Coulomb_Yield_Function(myState State, myPotencialPlastic PotencialPlastic )
	    :FluencyCriteria()
	    {
              mState            = State;
              mPotencialPlastic = PotencialPlastic;
	    }

             Modified_Morh_Coulomb_Yield_Function::~Modified_Morh_Coulomb_Yield_Function() {}

    
//***********************************************************************
//***********************************************************************

		    void Modified_Morh_Coulomb_Yield_Function::InitializeMaterial(const Properties& props) 
                   {
			mprops = &props;
			double maxfriction_angle  = (*mprops)[MAX_FRICTION_INTERNAL_ANGLE]*PI/180.00;
			mEta = 2.00*tan(maxfriction_angle/2.00 + PI/4.00); 
                   }
 
		     

		    void Modified_Morh_Coulomb_Yield_Function:: CalculateEquivalentUniaxialStress(
		    const Vector& StressVector,double& Result) 
                     { 
                         CalculateEquivalentUniaxialStressViaPrincipalStress(StressVector, Result);
                         //CalculateEquivalentUniaxialStressViaInvariants(StressVector, Result);
                     }


		    void Modified_Morh_Coulomb_Yield_Function::CalculateEquivalentUniaxialStressViaPrincipalStress(
		    const Vector& StressVector,double& Result)
                       {     
                       int    iter      = 50;
                       double zero      = 1.0E-9;
                       Matrix EigenVectors    = ZeroMatrix(3,3);
                       Matrix StressTensor    = ZeroMatrix(3,3);
                       Vector PrincipalStress = ZeroVector(3);
                       this->State_Tensor(StressVector,StressTensor);
                     
                      SD_MathUtils<double>::EigenVectors(StressTensor, EigenVectors,PrincipalStress, zero, iter);
                      
                      ///*sigma_1 >  sigma_2 > sigma_3 
                      sort (PrincipalStress.begin(), PrincipalStress.end()); 
                      reverse(PrincipalStress.begin(),PrincipalStress.end()); 
                      
                      //KRATOS_WATCH(PrincipalStress) 
    
                      double  frictional_internal = (*mprops)[MAX_FRICTION_INTERNAL_ANGLE]*PI/180.00;       
                      double  sin_fric_inter      = sin(frictional_internal);
                      double  cos_fric_inter      = cos(frictional_internal);
                       

                      mMultisurface_Platicity_Sigma  = ZeroVector(6);
                      mMultisurface_Platicity_Yield    = ZeroVector(6);

                      ///* Multisurface Representation                        
 
                      mMultisurface_Platicity_Sigma[0] =   (PrincipalStress[0] - PrincipalStress[2]) +  (PrincipalStress[0] + PrincipalStress[2]) *  sin_fric_inter; 
                      mMultisurface_Platicity_Yield[0]   =  mMultisurface_Platicity_Sigma[0] - 2.00 * mcohesion*cos_fric_inter; 
                        

                      mMultisurface_Platicity_Sigma[1] =  (PrincipalStress[1] - PrincipalStress[2]) + (PrincipalStress[1] + PrincipalStress[2]) *  sin_fric_inter;
                      mMultisurface_Platicity_Yield[1]   =  mMultisurface_Platicity_Sigma[1] - 2.00*mcohesion*cos_fric_inter; 


                      mMultisurface_Platicity_Sigma[2] =  (PrincipalStress[1] - PrincipalStress[0]) + (PrincipalStress[1] + PrincipalStress[0]) *  sin_fric_inter;
                      mMultisurface_Platicity_Yield[2]   =  mMultisurface_Platicity_Sigma[2] - 2.00*mcohesion*cos_fric_inter; 


                      mMultisurface_Platicity_Sigma[3] =  (PrincipalStress[2] - PrincipalStress[0]) + (PrincipalStress[2] + PrincipalStress[0]) *  sin_fric_inter;
                      mMultisurface_Platicity_Yield[3]   =  mMultisurface_Platicity_Sigma[3] - 2.00*mcohesion*cos_fric_inter; 


                      mMultisurface_Platicity_Sigma[4] =  (PrincipalStress[2] - PrincipalStress[1]) +  (PrincipalStress[2] + PrincipalStress[1]) *  sin_fric_inter;
                      mMultisurface_Platicity_Yield[4]   =  mMultisurface_Platicity_Sigma[4] - 2.00*mcohesion*cos_fric_inter; 


                      mMultisurface_Platicity_Sigma[5] =  (PrincipalStress[0] - PrincipalStress[1]) + (PrincipalStress[0] + PrincipalStress[1]) *  sin_fric_inter;
                      mMultisurface_Platicity_Yield[5]   =  mMultisurface_Platicity_Sigma[5] - 2.00*mcohesion*cos_fric_inter; 

                            
		      /*KRATOS_WATCH(mMultisurface_Platicity_Sigma[0]) 
		      KRATOS_WATCH(mMultisurface_Platicity_Sigma[1])
		      KRATOS_WATCH(mMultisurface_Platicity_Sigma[2])
		      KRATOS_WATCH(mMultisurface_Platicity_Sigma[3])
                      KRATOS_WATCH(mMultisurface_Platicity_Sigma[4])
                      KRATOS_WATCH(mMultisurface_Platicity_Sigma[5]) */
                      Result   =  mMultisurface_Platicity_Yield[0];
                      //Result   = (*std::max_element(mMorh_Coulomb_Yield.begin(), mMorh_Coulomb_Yield.end()));                                 
                      }



void Modified_Morh_Coulomb_Yield_Function::CalculateEquivalentUniaxialStressViaInvariants(
const Vector& StressVector,double& Result)
{

      double tetha_Lode          = 0.00;
      double frictional_internal = 0.00;
      double R_morh              = 0.00;
      double Alfa                = 0.00;
      double K_one               = 0.00;
      double K_two               = 0.00;
      double K_three             = 0.00; 
      double sin_fric_inter      = 0.00;
      double cos_fric_inter      = 0.00;
      double sin_tetha_lode      = 0.00;
      double cos_tetha_lode      = 0.00;

      Vector I        = ZeroVector(3);
      Vector J        = ZeroVector(3);
      Vector J_des    = ZeroVector(3);		      

      Matrix StressTensor     = ZeroMatrix(3,3);
      this->State_Tensor(StressVector,StressTensor);
      Tensor_Utils<double>::TensorialInvariants(StressTensor, I, J, J_des);
      
      if (J_des[1]==0.00) 
        {
         tetha_Lode = 0.00;                              
        }
      else
        {  
	tetha_Lode = -(3.00*sqrt(3.00)*J_des[2])/(2.00*pow(J_des[1], 1.50));
	if(fabs(tetha_Lode) > 1.00){tetha_Lode = 1.00; }
	tetha_Lode = asin(tetha_Lode)/3.00; 
	}

      //PrincipalStress[0] = I[0]/3.00 + sqrt(J_des[1]) * (cos(tetha_Lode) - sin(tetha_Lode)/sqrt(3.00));   
      //PrincipalStress[1] = I[0]/3.00 + 2.00 * sqrt(J_des[1]/3.00) * (sin(tetha_Lode));
      //PrincipalStress[2] = I[0]/3.00 + sqrt(J_des[1]) * (-cos(tetha_Lode) - sin(tetha_Lode)/sqrt(3.00));
      //KRATOS_WATCH(PrincipalStress) 
      frictional_internal = (*mprops)[MAX_FRICTION_INTERNAL_ANGLE]*PI/180.00;


      R_morh = tan(frictional_internal/2.00 + PI/4.00);
      R_morh = R_morh*R_morh;

      Alfa = ((*mprops)[FC]/(*mprops)[FT])/R_morh;
      sin_fric_inter = sin(frictional_internal);
      cos_fric_inter = cos(frictional_internal);
      sin_tetha_lode = sin(tetha_Lode);
      cos_tetha_lode = cos(tetha_Lode);

      K_one    = 1.00;                //0.50*((1.00 + Alfa)-(1.00-Alfa)*sin_fric_inter);
      K_two    = 1.00;               // 1E15;
      if(sin_fric_inter!=0.00)  
      {K_two   = 1.00;}              //0.50*((1.00 + Alfa)-(1.00-Alfa)/sin_fric_inter);}
      K_three  = sin_fric_inter;    // 0.50*((1.00 + Alfa)*sin_fric_inter-(1.00-Alfa));

      mMultisurface_Platicity_Sigma[0] =  I[0]*K_three/3.00 + sqrt(J_des[1])*(K_one * cos_tetha_lode - K_two * sin_tetha_lode * sin_fric_inter/sqrt(3.00));
      mMultisurface_Platicity_Yield[0]   =  mMultisurface_Platicity_Sigma[0] - mcohesion*cos_fric_inter; 

      
      mMultisurface_Platicity_Sigma[1] =  I[0]*K_three/3.00 + 0.50 * sqrt(J_des[1]) * ( K_one * cos_tetha_lode * (1.00 + sin_fric_inter) +  K_two * sin_tetha_lode * ( sin_fric_inter - 3.00)/sqrt(3.00) ) ;  
      mMultisurface_Platicity_Yield[1]   =  mMultisurface_Platicity_Sigma[1] - mcohesion*cos_fric_inter; 


      mMultisurface_Platicity_Sigma[2] =  I[0]*K_three/3.00  +  0.50 * sqrt(J_des[1]) * (K_one * cos_tetha_lode * (1.00  - sin_fric_inter) +  K_two * sin_tetha_lode * ( (sin_fric_inter + 3.00)/sqrt(3.00) ) ) ;
      mMultisurface_Platicity_Yield[2]   = mMultisurface_Platicity_Sigma[2] - mcohesion*cos_fric_inter; 

      mMultisurface_Platicity_Sigma[3] =  I[0]*K_three/3.00 - sqrt(J_des[1])*(K_one*cos_tetha_lode + K_two * sin_tetha_lode * sin_fric_inter/sqrt(3.00));
      mMultisurface_Platicity_Yield[3]  = mMultisurface_Platicity_Sigma[3] - mcohesion*cos_fric_inter; 

      //KRATOS_WATCH(mMorh_Coulomb_Yield[0]) 
      //KRATOS_WATCH(mMorh_Coulomb_Yield[1])
      //KRATOS_WATCH(mMorh_Coulomb_Yield[2])
      //KRATOS_WATCH(mMorh_Coulomb_Yield[3])
      //KRATOS_WATCH("***********************") 
      Result   =  mMultisurface_Platicity_Yield[0];
      //Result   = (*std::max_element(mMorh_Coulomb_Yield.begin(), mMorh_Coulomb_Yield.end()));  
}


void Modified_Morh_Coulomb_Yield_Function::CalculateEquivalentUniaxialStressViaCilindricalCoordinate(
const Vector& StressVector,double& Result){}

// solo respecto de las tensiones
void Modified_Morh_Coulomb_Yield_Function::CalculateDerivateFluencyCriteria(const Vector& StressVector, Vector& DerivateFluencyCriteria)
{


}

void Modified_Morh_Coulomb_Yield_Function::CalculateDerivatePotencialFlowCriteria(const Vector& StressVector, Vector& DerivatePotencialFlow)
{

}

void Modified_Morh_Coulomb_Yield_Function::UpdateVariables(const Vector& Variables)
{
 mcohesion        = Variables(0);
}


void Modified_Morh_Coulomb_Yield_Function::CalculateDerivateFluencyCriteriaMultiSurface(const Vector& StressVector,  vector<Vector>& DerivateFluencyCriteria)
{
      KRATOS_TRY  

      
      double tetha_Lode      = 0.00;
      double R_morh          = 0.00;
      double Alfa            = 0.00;
      double K_one           = 0.00;
      double K_two           = 0.00;
      double K_three         = 0.00; 
      double sin_fric_inter  = 0.00;
      double cos_fric_inter  = 0.00;
      double sin_tetha_lode  = 0.00;
      double cos_tetha_lode  = 0.00;
      double frictional_internal = 0.00;

      Second_Order_Tensor a;  
      Second_Order_Tensor C;  
      Matrix StressTensor      = ZeroMatrix(3,3);
      array_1d<double,6 > A; A = ZeroVector(6);    

      DerivateFluencyCriteria.resize(4);
      DerivateFluencyCriteria[0].resize(6);
      DerivateFluencyCriteria[1].resize(6);
      DerivateFluencyCriteria[2].resize(6);
      DerivateFluencyCriteria[3].resize(6);

      ///* Vector de Flujo
      a.resize(3, false);  
      a[0].resize(6, false); a[0] = ZeroVector(6);
      a[1].resize(6, false); a[1] = ZeroVector(6);
      a[2].resize(6, false); a[2] = ZeroVector(6);

      ///* Superficies de morh-colulomb = 4;
      C.resize(4, false);
      C[0].resize(3, false); C[0] = ZeroVector(3);
      C[1].resize(3, false); C[1] = ZeroVector(3);
      C[2].resize(3, false); C[2] = ZeroVector(3);
      C[3].resize(3, false); C[3] = ZeroVector(3);

      frictional_internal = (*mprops)[MAX_FRICTION_INTERNAL_ANGLE]*PI/180.00;                                   

      Vector I          = ZeroVector(3);
      Vector J          = ZeroVector(3);
      Vector J_des      = ZeroVector(3);      				  

      this->State_Tensor(StressVector,StressTensor);
      Tensor_Utils<double>::TensorialInvariants(StressTensor, I, J, J_des);

      ///* Calculating Lode angle
      if (J_des[1]==0.00) 
      {
      tetha_Lode = 0.00;                              
      }
      else
      {  
      tetha_Lode = -(3.00*sqrt(3.00)*J_des(2))/(2.00*pow(J_des[1], 1.50));
      if(fabs(tetha_Lode) > 1.00){tetha_Lode = 1.00; }
      tetha_Lode = asin(tetha_Lode)/3.00; 
      }

      R_morh = tan(frictional_internal/2.00 + PI/4.00);
      R_morh = R_morh*R_morh;

      Alfa = ((*mprops)[FC]/(*mprops)[FT])/R_morh;
      sin_fric_inter = sin(frictional_internal);
      cos_fric_inter = cos(frictional_internal);
      sin_tetha_lode = sin(tetha_Lode);
      cos_tetha_lode = cos(tetha_Lode);

      ///*WARNING = This is the standard  morh coolomb 
      K_one    = 1.00; // 0.50*((1.00 + Alfa)-(1.00-Alfa)*sin_fric_inter);
      K_two    = 1.00; // 1E15;
      if(sin_fric_inter!=0.00)  
      {K_two   = 1.00;} // 0.50*((1.00 + Alfa)-(1.00-Alfa)/sin_fric_inter);}
      K_three  = sin_fric_inter; // 0.50*((1.00 + Alfa)*sin_fric_inter-(1.00-Alfa));

      C(0)[0] =  K_three/3.00;  ///*OK
      C(1)[0] =  K_three/3.00;  ///*OK
      C(2)[0] =  K_three/3.00; ///*OK
      C(3)[0] =  K_three/3.00;  ///*OK
      
      double& s  = sin_tetha_lode; ///*OK 
      double& c  = cos_tetha_lode; ///*OK       
      double& sf = sin_fric_inter; ///*OK
      

      double dF0_J =  (K_one  * c - (K_two * s * sf) /sqrt(3.00) );  ///*OK
      double dF1_J =   0.50 * ( K_one * c * (1.00 + sf) +  K_two * s * ( sf - 3.00)/sqrt(3.00));
      double dF2_J =   0.50 * ( K_one * c * (1.00 - sf) +  K_two * s * ( sf + 3.00)/sqrt(3.00));
      double dF3_J =  - ( K_one  * c  +  (K_two * s * sf)/sqrt(3.00) ); ///*OK 
 

      if (fabs(tetha_Lode) <= 0.50614548307836) // angulos comprendidos entre (-29 < tetha < +29  )
      {    
     
      double dF0_T  = sqrt(J_des[1]) * ( -K_one * s - K_two * c * sf/sqrt(3.00) ); ///*OK  
      double dF1_T  = 0.50 * sqrt(J_des[1]) * ( -K_one * s * (1.00 + sf ) + K_two * c * (sf - 3.00)/sqrt(3.00)); ///*OK
      double dF2_T  = 0.50 * sqrt(J_des[1]) * ( -K_one * s * (1.00 - sf ) + K_two * c * (sf + 3.00)/sqrt(3.00));  
      double dF3_T  = -sqrt(J_des[1]) * ( -K_one * s + K_two * c * sf/sqrt(3.00));


      double fact_1 = tan(3.00 * tetha_Lode) / sqrt(J_des[1]); ///*OK  
      double fact_2 = -sqrt(3.00) / (2.00 * cos(3.00 * tetha_Lode) * J_des[1] * sqrt(J_des[1])); ///*OK  

      C(0)[1] = dF0_J - dF0_T * fact_1;     C(0)[2] = dF0_T * fact_2;  
      C(1)[1] = dF1_J - dF1_T * fact_1;     C(1)[2] = dF1_T * fact_2;
      C(2)[1] = dF2_J - dF2_T * fact_1;     C(2)[2] = dF2_T * fact_2;
      C(3)[1] = dF3_J - dF3_T * fact_1;     C(3)[2] = dF3_T * fact_2; 
    
       }
      
      else  // para angulos (-30 < tehta < -29  and  30 > tehta > 29 )
      {    
      C(0)[1] =  dF0_J; 
      C(1)[1] =  dF1_J; 
      C(2)[1] =  dF2_J;  
      C(3)[1] =  dF3_J;  

      C(0)[2] =  0.00; 
      C(1)[2] =  0.00;
      C(2)[2] =  0.00;
      C(3)[2] =  0.00; 
      }  
      
      this->CalculateVectorFlowDerivate(StressVector, a);

      for(unsigned  int i=0; i<4; i++)
       { for(unsigned int j =0; j<3; j++)
           {
             noalias(A) += a[j]*C[i][j]; 
           }
           for(unsigned int k = 0; k<6; k++)
              { 
               if (fabs(A[k] ) < 1E-14) {A[k] = 0.00; }
              } 
           noalias(DerivateFluencyCriteria[i]) = A;
           A = ZeroVector(6);
       }  
      //KRATOS_WATCH(DerivateFluencyCriteria)
      //CalculateDerivatePrincipalStress_Fluency(StressVector,DerivateFluencyCriteria);
      KRATOS_CATCH("")

}

void Modified_Morh_Coulomb_Yield_Function::CalculateDerivatePotencialFlowCriteriaMultiSurface(const Vector& StressVector, vector<Vector>& DerivatePotencialFlow)
{
KRATOS_TRY  
      
      double tetha_Lode      = 0.00;
      double R_morh          = 0.00;
      double Alfa            = 0.00;
      double K_one           = 0.00;
      double K_two           = 0.00;
      double K_three         = 0.00; 
      double sin_fric_inter  = 0.00;
      double cos_fric_inter  = 0.00;
      double sin_tetha_lode  = 0.00;
      double cos_tetha_lode  = 0.00;
      double frictional_internal = 0.00;

      Second_Order_Tensor a;  
      Second_Order_Tensor C;  
      Matrix StressTensor      = ZeroMatrix(3,3);
      array_1d<double,6 > A; A = ZeroVector(6);    

      DerivatePotencialFlow.resize(4);
      DerivatePotencialFlow[0].resize(6);
      DerivatePotencialFlow[1].resize(6);
      DerivatePotencialFlow[2].resize(6);
      DerivatePotencialFlow[3].resize(6);
      ///* Vector de Flujo
      a.resize(3, false);  
      a[0].resize(6, false); a[0] = ZeroVector(6);
      a[1].resize(6, false); a[1] = ZeroVector(6);
      a[2].resize(6, false); a[2] = ZeroVector(6);

      ///* Superficies de morh-colulomb = 4;
      C.resize(4, false);
      C[0].resize(3, false); C[0] = ZeroVector(3);
      C[1].resize(3, false); C[1] = ZeroVector(3);
      C[2].resize(3, false); C[2] = ZeroVector(3);
      C[3].resize(3, false); C[3] = ZeroVector(3);

      frictional_internal = (*mprops)[MAX_DILATANCY_ANGLE]*PI/180.00;                                   

      Vector I          = ZeroVector(3);
      Vector J          = ZeroVector(3);
      Vector J_des      = ZeroVector(3);      				  

      this->State_Tensor(StressVector,StressTensor);
      Tensor_Utils<double>::TensorialInvariants(StressTensor, I, J, J_des);

      ///* Calculating Lode angle
      if (J_des[1]==0.00) 
      {
      tetha_Lode = 0.00;                              
      }
      else
      {  
      tetha_Lode = -(3.00*sqrt(3.00)*J_des(2))/(2.00*pow(J_des[1], 1.50));
      if(fabs(tetha_Lode) > 1.00){tetha_Lode = 1.00; }
      tetha_Lode = asin(tetha_Lode)/3.00; 
      }

      R_morh = tan(frictional_internal/2.00 + PI/4.00);
      R_morh = R_morh*R_morh;

      Alfa = ((*mprops)[FC]/(*mprops)[FT])/R_morh;
      sin_fric_inter = sin(frictional_internal);
      cos_fric_inter = cos(frictional_internal);
      sin_tetha_lode = sin(tetha_Lode);
      cos_tetha_lode = cos(tetha_Lode);

      ///*WARNING = This is the standard  morh coolomb 
      K_one    = 1.00; //0.50*((1.00 + Alfa)-(1.00-Alfa)*sin_fric_inter);
      K_two    = 1.00; // 1E15;
      if(sin_fric_inter!=0.00)  
      {K_two   =  1.00;} //0.50*((1.00 + Alfa)-(1.00-Alfa)/sin_fric_inter);}
      K_three  =  sin_fric_inter; //0.50*((1.00 + Alfa)*sin_fric_inter-(1.00-Alfa));

      C(0)[0] =  K_three/3.00;  ///*OK
      C(1)[0] =  K_three/3.00;  ///*OK
      C(2)[0] =  K_three/3.00; ///*OK
      C(3)[0] =  K_three/3.00;  ///*OK
      
      double& s  = sin_tetha_lode; ///*OK 
      double& c  = cos_tetha_lode; ///*OK       
      double& sf = sin_fric_inter; ///*OK
      

      double dF0_J =  (K_one  * c - (K_two * s * sf) /sqrt(3.00) );  ///*OK
      double dF1_J =   0.50 * ( K_one * c * (1.00 + sf) +  K_two * s * ( sf - 3.00)/sqrt(3.00));
      double dF2_J =   0.50 * ( K_one * c * (1.00 - sf) +  K_two * s * ( sf + 3.00)/sqrt(3.00));
      double dF3_J =  - ( K_one  * c  +  (K_two * s * sf)/sqrt(3.00) ); ///*OK 
 

      if (fabs(tetha_Lode) <= 0.50614548307836) // angulos comprendidos entre (-29 < tetha < +29  )
      {    
     
      double dF0_T  = sqrt(J_des[1]) * ( -K_one * s - K_two * c * sf/sqrt(3.00) ); ///*OK  
      double dF1_T  = 0.50 * sqrt(J_des[1]) * ( -K_one * s * (1.00 + sf ) + K_two * c * (sf - 3.00)/sqrt(3.00)); ///*OK
      double dF2_T  = 0.50 * sqrt(J_des[1]) * ( -K_one * s * (1.00 - sf ) + K_two * c * (sf + 3.00)/sqrt(3.00));  
      double dF3_T  = -sqrt(J_des[1]) * ( -K_one * s + K_two * c * sf/sqrt(3.00));


      double fact_1 = tan(3.00 * tetha_Lode) / sqrt(J_des[1]); ///*OK  
      double fact_2 = -sqrt(3.00) / (2.00 * cos(3.00 * tetha_Lode) * J_des[1] * sqrt(J_des[1])); ///*OK  

      C(0)[1] = dF0_J - dF0_T * fact_1;     C(0)[2] = dF0_T * fact_2;  
      C(1)[1] = dF1_J - dF1_T * fact_1;     C(1)[2] = dF1_T * fact_2;
      C(2)[1] = dF2_J - dF2_T * fact_1;     C(2)[2] = dF2_T * fact_2;
      C(3)[1] = dF3_J - dF3_T * fact_1;     C(3)[2] = dF3_T * fact_2; 
    
       }
      
      else  // para angulos (-30 < tehta < -29  and  30 > tehta > 29 )
      {    
      C(0)[1] =  dF0_J; 
      C(1)[1] =  dF1_J; 
      C(2)[1] =  dF2_J;  
      C(3)[1] =  dF3_J;  

      C(0)[2] =  0.00; 
      C(1)[2] =  0.00;
      C(2)[2] =  0.00;
      C(3)[2] =  0.00; 
      }  
      
      this->CalculateVectorFlowDerivate(StressVector, a);

      for(unsigned  int i=0; i<4; i++)
       { for(unsigned int j =0; j<3; j++)
           {
             noalias(A) += C[i][j]*a[j]; 
           }
           for(unsigned int k = 0; k<6; k++)
              { 
               if (fabs(A[k] ) < 1E-14) {A[k] = 0.00; }
              } 
           noalias(DerivatePotencialFlow[i]) = A;
           A = ZeroVector(6);
       }  
     
      //KRATOS_WATCH(DerivatePotencialFlow);
      //CalculateDerivatePrincipalStress_Flow(StressVector, DerivatePotencialFlow);
      //KRATOS_WATCH(DerivatePotencialFlow);
      //std::cout<<"-------------------------------------"<<std::endl;   

      KRATOS_CATCH("")

}

void Modified_Morh_Coulomb_Yield_Function::CalculateDerivatePrincipalStress_Fluency(const Vector& StressVector, vector<Vector>& DerivateFluencyCriteria)
{


vector<Vector> DerivateFluency; DerivateFluency.resize(4);
DerivateFluencyCriteria.resize(6);
DerivateFluency[0].resize(3);
DerivateFluency[1].resize(3);
DerivateFluency[2].resize(3);
DerivateFluency[3].resize(3);

DerivateFluencyCriteria.resize(4);
DerivateFluencyCriteria[0].resize(6); DerivateFluencyCriteria[0] = ZeroVector(6);
DerivateFluencyCriteria[1].resize(6); DerivateFluencyCriteria[1] = ZeroVector(6);
DerivateFluencyCriteria[2].resize(6); DerivateFluencyCriteria[2] = ZeroVector(6);
DerivateFluencyCriteria[3].resize(6); DerivateFluencyCriteria[3] = ZeroVector(6);

double  frictional_internal = (*mprops)[MAX_FRICTION_INTERNAL_ANGLE]*PI/180.00;      
double  s = sin(frictional_internal);

DerivateFluency[0](0) =  0.5 * ( s + 1.00);  
DerivateFluency[0](1) =  0.00;
DerivateFluency[0](2) =  0.5 * ( s - 1.00);

DerivateFluency[1](0) =  0.5 * ( s + 1.00);  
DerivateFluency[1](1) =  0.5 * ( s - 1.00);
DerivateFluency[1](2) =  0.00;

DerivateFluency[2](0) =  0.00;  
DerivateFluency[2](1) =  0.5 * ( s + 1.00);
DerivateFluency[2](2) =  0.5 * ( s - 1.00);

DerivateFluency[3](0) =  0.5 * ( s - 1.00);  
DerivateFluency[3](1) =  0.00;
DerivateFluency[3](2) =  0.5 * ( s + 1.00);

Second_Order_Tensor DerivateStress;
this-> CalculateDerivatePrincipalStress(StressVector, DerivateStress);

for(unsigned int i=0; i<4; i++)
{
      for (unsigned int j=0; j<3; j++)
      {noalias(DerivateFluencyCriteria[i]) = DerivateFluencyCriteria[i] +  DerivateFluency[i](j) * DerivateStress[j];}     
}

}

void Modified_Morh_Coulomb_Yield_Function::CalculateDerivatePrincipalStress_Flow(const Vector& StressVector, vector<Vector>& DerivatePotencialFlow)
{

vector<Vector> DerivatePotencial; DerivatePotencial.resize(4);
DerivatePotencial[0].resize(3);
DerivatePotencial[1].resize(3);
DerivatePotencial[2].resize(3);
DerivatePotencial[3].resize(3);

DerivatePotencialFlow.resize(4);
DerivatePotencialFlow[0].resize(6); DerivatePotencialFlow[0] = ZeroVector(6);
DerivatePotencialFlow[1].resize(6); DerivatePotencialFlow[0] = ZeroVector(6);
DerivatePotencialFlow[2].resize(6); DerivatePotencialFlow[0] = ZeroVector(6);
DerivatePotencialFlow[3].resize(6); DerivatePotencialFlow[0] = ZeroVector(6);

double  dilatancy = (*mprops)[MAX_DILATANCY_ANGLE]*PI/180.00;      
double  s         = sin(dilatancy);

DerivatePotencial[0](0) =  0.5 * ( s + 1.00);  
DerivatePotencial[0](1) =  0.00;
DerivatePotencial[0](2) =  0.5 * ( s - 1.00);

DerivatePotencial[1](0) =  0.5 * ( s + 1.00);  
DerivatePotencial[1](1) =  0.5 * ( s - 1.00);
DerivatePotencial[1](2) =  0.00;

DerivatePotencial[2](0) =  0.00;  
DerivatePotencial[2](1) =  0.5 * ( s + 1.00);
DerivatePotencial[2](2) =  0.5 * ( s - 1.00);

DerivatePotencial[3](0) =  0.5 * ( s - 1.00);  
DerivatePotencial[3](1) =  0.00;
DerivatePotencial[3](2) =  0.5 * ( s + 1.00);

Second_Order_Tensor DerivateStress;
this-> CalculateDerivatePrincipalStress(StressVector, DerivateStress);

for(unsigned int i=0; i<4; i++)
{
      for (unsigned int j=0; j<3; j++)
      {noalias(DerivatePotencialFlow[i]) = DerivatePotencialFlow[i] +  DerivatePotencial[i](j) * DerivateStress[j];}     
}

 
}  




}




