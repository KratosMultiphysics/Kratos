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
#include "custom_utilities/tensor_utils.h"
#include "fluency_criteria/isotropic_rankine_yield_function.h"
#include <cmath>



namespace Kratos
  {

            Isotropic_Rankine_Yield_Function::Isotropic_Rankine_Yield_Function() {}       
    
            Isotropic_Rankine_Yield_Function::Isotropic_Rankine_Yield_Function(
            const SoftHardPointerType& SofteningBehavior,
            const myState& State)
	    :FluencyCriteria()
	    {
	      mpSofteningBehaviorFt  = SofteningBehavior; 
	      mState                 = State;   
	    }

             Isotropic_Rankine_Yield_Function::~Isotropic_Rankine_Yield_Function() {}

    
//***********************************************************************
//***********************************************************************

	  void Isotropic_Rankine_Yield_Function::InitializeMaterial(const Properties& props) 
	  {   
           
           mprops                                      =  &props;
	   mFt                                         =  (*mprops)[FT];
	   mcurrent_Ft                                 =  mFt;
           mrankine_accumulated_plastic_strain_current = 0.00;  
           mrankine_accumulated_plastic_strain_old     = 0.00;
	   mpastic_damage_old                          = 0.00;
           mpastic_damage_current                      = 0.00;

	   

	   int  size = 4;
	   if(mState==Tri_D)
	      size = 6;

	   mplastic_strain.resize(size, false); 
	   mplastic_strain_old.resize(size, false);

	   mplastic_strain                             =  ZeroVector(size);
	   mplastic_strain_old                         =  ZeroVector(size);   
           mElastic_strain                             =  ZeroVector(size); 
           mElastic_strain_old                         =  ZeroVector(size); 
	   mPrincipalPlasticStrain_current             =  ZeroVector(3);
	   mPrincipalPlasticStrain_old                 =  ZeroVector(3);
	   mpSofteningBehaviorFt->InitializeMaterial(props);  

	  }

	  void Isotropic_Rankine_Yield_Function::UpdateMaterial()
	  {
	   
	    	                                                  
           mpastic_damage_current                      = mpastic_damage_old;  
           mcurrent_Ft                                 = mFt;  
	   mrankine_accumulated_plastic_strain_current = mrankine_accumulated_plastic_strain_old;
	   noalias(mPrincipalPlasticStrain_current)    = mPrincipalPlasticStrain_old;
	   noalias(mplastic_strain)                    = mplastic_strain_old; 
	   noalias(mElastic_strain)                    = mElastic_strain_old;
	  }

   
          bool Isotropic_Rankine_Yield_Function::CheckPlasticAdmisibility(const Vector& Stress)
          { 
	     array_1d<double, 3> check = ZeroVector(3);
	     double toler = 1E-6 * mcurrent_Ft;
	     if(toler>0.10) toler = 0.10;
	     if(toler<1E-3) toler = 1E-3;
	     
	     check[0] = Stress[0] - mcurrent_Ft;
	     check[1] = Stress[1] - mcurrent_Ft;
	     check[2] = Stress[2] - mcurrent_Ft;
	     if(check[0]>= toler || check[1]>= toler ||  check[2]>= toler)
	        return true;
	     else 
	        return false;
	  }

          int Isotropic_Rankine_Yield_Function::CheckSurface(const Vector& Stress)
          {
	     int count = 0;
	     if(Stress[0]>mcurrent_Ft) count++;
	     if(Stress[1]>mcurrent_Ft) count++;
	     if(Stress[2]>mcurrent_Ft) count++;    
	     return count;
	  }

          void Isotropic_Rankine_Yield_Function::ReturnMapping(const Vector& StrainVector, const Vector& TrialStress, Vector& StressVector)
            {             
	      
	        if(mcurrent_Ft>0.00){
		//const double& Young   = (*mprops)[YOUNG_MODULUS];
		//const double& Poisson = (*mprops)[POISSON_RATIO];
		//const double Gmodu    = Young/(2.00 * (1.00 + Poisson) );
		//const double Bulk     = Young/(3.00 * (1.00-2.00*Poisson)); 

		array_1d<double,3> PrincipalStress = ZeroVector(3);
		array_1d<double,3> Sigma = ZeroVector(3);
		array_1d<array_1d < double,3 > ,3> EigenVectors;
		array_1d<unsigned int,3> Order;
		//const double d3  = 0.3333333333333333333; 

		const int dim =  TrialStress.size();
		Vector Stress(dim);
		Stress = ZeroVector(dim);
		/// computing the trial kirchooff strain
		SpectralDecomposition(TrialStress, PrincipalStress, EigenVectors);
		IdentifyMaximunAndMinumumPrincipalStres_CalculateOrder(PrincipalStress, Order);
		noalias(Sigma)  = PrincipalStress;
		noalias(Stress) = TrialStress;
		
		
                const bool check = CheckPlasticAdmisibility(PrincipalStress);  
		if(check==true){
		 Vector delta_lamda;		 
		 // return to main plane
                 bool check = One_Vector_Return_Mapping_To_Main_Plane(PrincipalStress, delta_lamda,  Sigma);    

		 
                 if(check==false)
                  {
                     //return to corner
		     UpdateMaterial();
                     check =  Two_Vector_Return_Mapping_To_Corner(PrincipalStress, delta_lamda, Sigma);  
                     if (check==false)
                          {
                            //return to apex
			      UpdateMaterial();
                              Three_Vector_Return_Mapping_To_Apex (PrincipalStress, delta_lamda, Sigma);   
                          }         
                  }

		 AssembleUpdateStressAndStrainTensor(Sigma,  EigenVectors,  Order, StrainVector, Stress);

		 CalculatePlasticDamage(Sigma);
		 Matrix PlasticTensor;  PlasticTensor.resize(3,3,false); 
		 PlasticTensor = ZeroMatrix(3,3);
		 /// Updating the  elastic plastic strain tensor
		 for(unsigned int i=0; i<3; i++)     
		     noalias(PlasticTensor) +=  mPrincipalPlasticStrain_current[i] * Matrix(outer_prod(EigenVectors[Order[i]],  EigenVectors[Order[i]]));

		/// WARNING = Plane Strain       
		mplastic_strain[0] = PlasticTensor(0,0);
		mplastic_strain[1] = PlasticTensor(1,1);
		mplastic_strain[2] = 2.00 * PlasticTensor(1,0);
		mplastic_strain[3] = PlasticTensor(2,2);
		 
                 //CalculatePlasticStrain(StrainVector, StressVector, mplastic_strain, Sigma[2]);
		} 

	      CalculateElasticStrain(Stress, mElastic_strain);
	     
	      /// WARNING = Plane Strain 
	      StressVector[0] = Stress[0];
	      StressVector[1] = Stress[1];
	      StressVector[2] = Stress[2];
	      //KRATOS_WATCH("----------------------------")
            }
	    }
	    
 
void Isotropic_Rankine_Yield_Function::FinalizeSolutionStep()
    {
          mpastic_damage_old                          =  mpastic_damage_current;
	  mFt                                         =  mcurrent_Ft;    
	  mrankine_accumulated_plastic_strain_old     =  mrankine_accumulated_plastic_strain_current;
	  noalias(mPrincipalPlasticStrain_old)        =  mPrincipalPlasticStrain_current; 
	  noalias(mplastic_strain_old)                =  mplastic_strain; 
	  noalias(mElastic_strain_old)                =  mElastic_strain;
 	  
    }


bool Isotropic_Rankine_Yield_Function::One_Vector_Return_Mapping_To_Main_Plane(const array_1d<double,3>& PrincipalStress, Vector& delta_lamda,  array_1d<double,3>& Sigma)
              {
              unsigned int iter       = 0;    
              double norma            = 1.00;   
	      double delta_lamda_a    = 0.00; 
	      double E                 = (*mprops)[YOUNG_MODULUS];
	      double NU                = (*mprops)[POISSON_RATIO];             
              double G                 = 0.5*E / (1.00 + NU);
              double K                 =  E / (3.00 * (1.00-2.00*NU));
              double H                 = 0.00;
	      double d                 = 0.00;  
	      double residual          = 0.00;
              const double toler       = 1E-3;
	      const double raiz_2_3    = 0.81649658092773;
	      unsigned int max         = 1000;
	      double Partial_Ep_gama_a = 0.00; 
     
	      
	      double Inc             = 0.00;
	      const double d3        = 0.3333333333333333;
	      const double raiz2d3   = 0.8164958092773;
	      double Ppvs = 0.00;           /// principal plastic volumetric strain  
	      array_1d<double,3> Ppds = ZeroVector(3);      /// principal plastic desviatoric strain
	      array_1d<double,3> Pps  = ZeroVector(3);       /// principal plastic  strain
	      array_1d<double,3> I; 
	      I[0] = 1.00;
	      I[1] = 1.00;
	      I[2] = 1.00;

              delta_lamda         =  ZeroVector(1); 
	      residual            =  PrincipalStress[0] - mcurrent_Ft;;
              norma               =  residual/mcurrent_Ft;
	      Vector Imput_Parameters(4);
	      Imput_Parameters    = ZeroVector(4);
	      Imput_Parameters[0] =  mhe; 
	      Imput_Parameters[1] =  mrankine_accumulated_plastic_strain_current;
	      	      
              while(iter++<=max && norma>= toler) 
	          {  
		    Partial_Ep_gama_a   = raiz_2_3; /// Acumulated Vond misses
		    H                   = mpSofteningBehaviorFt->FirstDerivateFunctionBehavior(Imput_Parameters);
                    d                   = -4.00 * G / 3.00 - K - H * Partial_Ep_gama_a;
                    delta_lamda[0]      = delta_lamda[0] - (residual / d);

		    
		    ///normal acuulated
                    ///mrankine_accumulated_plastic_strain_current= mrankine_accumulated_plastic_strain_old+ delta_lamda[0];
		    
		    /// von mises
		    mrankine_accumulated_plastic_strain_current = mrankine_accumulated_plastic_strain_old + raiz_2_3 * delta_lamda[0];   
		      /// computing Kp_punto.
		      /// volumetric and desviatoric plastic strain
		      Pps[0]        = delta_lamda_a;
		      Pps[1]        = 0.00;
		      Pps[2]        = 0.00;
		      Ppvs          = Pps[0] + Pps[1] + Pps[2]; 
		      noalias(Ppds) = Pps - d3 * Ppvs * I;
		      Inc           = raiz2d3* std::sqrt(inner_prod(Pps,Pps));
		      mrankine_accumulated_plastic_strain_current = mrankine_accumulated_plastic_strain_old + Inc;
		      ComputeActualStrees(Ppvs, Ppds, PrincipalStress, Sigma);
		      ComputeActualStrain(Pps);
		      CalculatePlasticDamage(Sigma);
		    
		    
		    ///* Updatinf mFt
		    Imput_Parameters[1] =  mrankine_accumulated_plastic_strain_current;
		    Imput_Parameters[2] =  mpastic_damage_old;
                    Imput_Parameters[3] =  mpastic_damage_current;
                    mcurrent_Ft         =  mpSofteningBehaviorFt->FunctionBehavior(Imput_Parameters);
                    delta_lamda_a       =  delta_lamda[0];   
		    
                    ///* comprobando si mft se cumplio   
		    if(mcurrent_Ft <= toler) 
                       {
		           mcurrent_Ft = toler;
                           break;         
                       } 
                    else
                    { 
                    ///* update teh current value             
		    residual    =  PrincipalStress[0] - mcurrent_Ft - delta_lamda_a * (4.00  * G / 3.00 + K);  
                    norma       =  fabs(residual/mcurrent_Ft);  
                    }
                  }   
                  
		  if(iter>=max){ //||  std::delta_lamda_a<0.00){  
		     //return false;
		     KRATOS_WATCH(norma) 
		     KRATOS_WATCH(PrincipalStress[0])
		     KRATOS_WATCH(PrincipalStress[1])
		     KRATOS_WATCH(PrincipalStress[2])
		     KRATOS_WATCH(delta_lamda_a)
		     KRATOS_WATCH(residual)
		     KRATOS_WATCH(mcurrent_Ft)
		     
		     std::cout<<  "RETURN MAPPING TO MAIN PLANE RANKINE  NOT CONVERGED" << std::endl;
		     KRATOS_ERROR(std::logic_error,  "RETURN MAPPING TO MAIN PLANE RANKINE  NOT CONVERGED" , "");
		  }
                  
                     ///* Updating Stress  
		    if(mcurrent_Ft <=toler) 
                    {
		    Sigma[0] = 0.00;
		    Sigma[1] = 0.00; 
		    Sigma[2] = 0.00;  
		    } 
                    else
                    { 
		    /// volumetric and desviatoric plastic strain
		    Pps[0]        = delta_lamda_a;
		    Pps[1]        = 0.00;
		    Pps[2]        = 0.00;
		    Ppvs          = Pps[0] + Pps[1] + Pps[2]; 
		    noalias(Ppds) = Pps - d3 * Ppvs * I;
	            
		    //Sigma[0] = PrincipalStress[0] - delta_lamda_a*(4.00  * G / 3.00 + K );    
		    //Sigma[1] = PrincipalStress[1] - delta_lamda_a*(-2.00 * G / 3.00 + K ); 
		    //Sigma[2] = PrincipalStress[2] - delta_lamda_a*(-2.00 * G / 3.00 + K );    
		     noalias(Sigma) = PrincipalStress - 2.00 * G * Ppds - K * Ppvs * I;
                    }  

                    bool check = CheckValidity(Sigma);
		    if(check==true) 
		     {
		       Vector PPS_bar(3);
		       PPS_bar = ZeroVector(3);
                       ComputePlasticStrainBar(mplastic_strain_old ,m_inv_DeltaF, PPS_bar);
		       //updating the correct principal pastic strain
		       mPrincipalPlasticStrain_current[0] = /*mPrincipalPlasticStrain_old[0]*/  PPS_bar[0] +  delta_lamda_a;
		       mPrincipalPlasticStrain_current[1] = /*mPrincipalPlasticStrain_old[1]*/  PPS_bar[1]; 
		       mPrincipalPlasticStrain_current[2] = /*mPrincipalPlasticStrain_old[2]*/  PPS_bar[2];
		       return true;
		     }
		   else
		     return false;
		   
		  
		
		}


    
    bool Isotropic_Rankine_Yield_Function::Two_Vector_Return_Mapping_To_Corner (const array_1d<double,3>& PrincipalStress, Vector& delta_lamda ,array_1d<double,3>& Sigma)
       {
	 
	      
              unsigned int iter    = 0;    
	      unsigned int max     = 1000; 
              int singular         = 0;  
	      const double E       = (*mprops)[YOUNG_MODULUS];
	      const double NU      = (*mprops)[POISSON_RATIO];             
              const double G       = 0.5*E / (1.00 + NU);
              const double K       =  E / (3.00 * (1.00-2.00*NU));
	      double Inc           = 0.00;
	      double norma         = 1.00;   
	      double delta_lamda_a = 0.00; 
              double delta_lamda_b = 0.00; 
              const  double toler  = 1E-3;
	      double H             = 0.00;
             
	      double Partial_Ep_gama_a = 0.00; 
	      double Partial_Ep_gama_b = 0.00; 
	      
	      Matrix d;             d.resize(2,2,false);      noalias(d)    = ZeroMatrix(2,2);  
              Matrix d_inv;         d_inv.resize(2,2,false);  noalias(d_inv)= ZeroMatrix(2,2); 
	      Vector residual      = ZeroVector(2);
              delta_lamda          = ZeroVector(2); 
	      d.resize(2,2);
	      d_inv.resize(2,2);
	      Vector Imput_Parameters(4);
	      Imput_Parameters    =  ZeroVector(4); 
	      Imput_Parameters[0] =  mhe; 
	      Imput_Parameters[1] =  mrankine_accumulated_plastic_strain_current;
	      residual[0] = PrincipalStress[0] - mcurrent_Ft;
	      residual[1] = PrincipalStress[1] - mcurrent_Ft;
	      H           = mpSofteningBehaviorFt->FirstDerivateFunctionBehavior(Imput_Parameters);
	      
	      
	      const double d3        = 0.3333333333333333;
	      const double raiz2d3    = 0.8164958092773;
	      double Ppvs = 0.00;           /// principal plastic volumetric strain  
	      array_1d<double,3> Ppds = ZeroVector(3);      /// principal plastic desviatoric strain
	      array_1d<double,3> Pps  = ZeroVector(3);       /// principal plastic  strain
	      array_1d<double,3> I; 
	      I[0] = 1.00;
	      I[1] = 1.00;
	      I[2] = 1.00;
	      
              while(iter++<=max && norma>= toler) 
		  {
		      
		      d(0,0) = -4.00 * G / 3.00 - K  - H*Partial_Ep_gama_a;  d(0,1)   =  2.00 * G / 3.00 - K  - H*Partial_Ep_gama_b;  
		      d(1,0) =  2.00 * G / 3.00 - K  - H*Partial_Ep_gama_a;  d(1,1)   = -4.00 * G / 3.00 - K  - H*Partial_Ep_gama_b;    
		      
		      singular      =  SD_MathUtils<double>::InvertMatrix(d, d_inv);    
		      delta_lamda   =  delta_lamda - Vector(prod(d_inv, residual)); 
		      delta_lamda_a =  delta_lamda[0];
		      delta_lamda_b =  delta_lamda[1];
                      
		      // normal
		      //mrankine_accumulated_plastic_strain_current= mrankine_accumulated_plastic_strain_old+ delta_lamda_a + delta_lamda_b;
		      
		      // von mises
		      Inc = 0.81649658092773 * norm_2(delta_lamda);
		      mrankine_accumulated_plastic_strain_current = mrankine_accumulated_plastic_strain_old + Inc;     
		      
		      /// computing Kp_punto.
		      /// volumetric and desviatoric plastic strain
		      Pps[0]        = delta_lamda_a;
		      Pps[1]        = delta_lamda_b;
		      Pps[2]        = 0.00;
		      Ppvs          = Pps[0] + Pps[1] + Pps[2]; 
		      noalias(Ppds) = Pps - d3 * Ppvs * I;
		      Inc           = raiz2d3* std::sqrt(inner_prod(Pps,Pps));
		      mrankine_accumulated_plastic_strain_current = mrankine_accumulated_plastic_strain_old + Inc;
		      ComputeActualStrees(Ppvs, Ppds, PrincipalStress, Sigma);
		      ComputeActualStrain(Pps);
		      CalculatePlasticDamage(Sigma);
		      
		      
                       // Updatinf mFt
		       Partial_Ep_gama_a   = (2.00/3.00) * delta_lamda_a/Inc;
		       Partial_Ep_gama_b   = (2.00/3.00) * delta_lamda_b/Inc; 
		       Imput_Parameters[1] =  mrankine_accumulated_plastic_strain_current;
		       Imput_Parameters[2] =  mpastic_damage_old;
                       Imput_Parameters[3] =  mpastic_damage_current;
                       mcurrent_Ft         =  mpSofteningBehaviorFt->FunctionBehavior(Imput_Parameters);  
		       H                   =  mpSofteningBehaviorFt->FirstDerivateFunctionBehavior(Imput_Parameters);
		       //mcurrent_Ft =  mFt - H * (mrankine_accumulated_plastic_strain_current- maccumulated_plastic_strain_old); 
                      
                       /// ft se anulan
                      if(mcurrent_Ft<=toler)
                        {
                            mcurrent_Ft = 0.00;
                            delta_lamda_a  = delta_lamda[0];  
                            delta_lamda_b  = delta_lamda[1];    
                            break;
                        }
                       else 
                       {
                                                 
                          
                 	  residual[0] = PrincipalStress[0] - mcurrent_Ft;
                          residual[1] = PrincipalStress[1] - mcurrent_Ft;
                       
                          residual[0]=  residual[0] - delta_lamda_a*( 4.00 * G / 3.00 + K ) - delta_lamda_b*( -2.00  * G / 3.00 + K );   
                          residual[1]=  residual[1] - delta_lamda_a*(-2.00 * G / 3.00 + K ) - delta_lamda_b*(  4.00  * G / 3.00 + K );   
                          norma      =  norm_2(residual); 

                        }
                      }

		  if(iter>=max){
		     KRATOS_WATCH(norma)
		     KRATOS_WATCH(PrincipalStress[0])
		     KRATOS_WATCH(PrincipalStress[1])
		     KRATOS_WATCH(PrincipalStress[2])
		     KRATOS_WATCH(delta_lamda)
		     std::cout<< "RETURN MAPPING SIGMA 1 AND 2 RANKINE  NOT CONVERGED" << std::endl;
		     KRATOS_ERROR(std::logic_error,  "RETURN MAPPING SIGMA 1 AND 2 RANKINE  NOT CONVERGED" , "");
		  }


                     if(mcurrent_Ft <= 0.00) 
                        {
                           Sigma[0] = 0.00;
                           Sigma[1] = 0.00;
                           Sigma[2] = 0.00;
                        }
                     else 
                      {
					    /// volumetric and desviatoric plastic strain
		    Pps[0]        = delta_lamda_a;
		    Pps[1]        = delta_lamda_b;
		    Pps[2]        = 0.00;
		    Ppvs          = Pps[0] + Pps[1] + Pps[2]; 
		    noalias(Ppds) = Pps - d3 * Ppvs * I;  
		    noalias(Sigma) = PrincipalStress - 2.00 * G * Ppds - K * Ppvs * I;
			
		    //Sigma[0] = PrincipalStress[0] - delta_lamda_a*( 4.00  * G / 3.00 + K ) - delta_lamda_b*( -2.00  * G / 3.00 + K );
		    //Sigma[1] = PrincipalStress[1] - delta_lamda_a*(-2.00 * G  / 3.00 + K ) - delta_lamda_b*( 4.00  * G / 3.00 + K ); 
		    //Sigma[2] = PrincipalStress[2] - (delta_lamda_a + delta_lamda_b) * (-2.00 * G / 3.00 + K );  
                     }
                      
		      bool check = CheckValidity(Sigma);
		      //bool check_2 = bool( delta_lamda_a > 0.00 && delta_lamda_b > 0.00); 
		      //bool check   = bool(checK_1==true && check_2==true);
		      if(check==true) 
		      {
			/// returns is consistent
			Vector PPS_bar(3);
			PPS_bar = ZeroVector(3);
                        ComputePlasticStrainBar(mplastic_strain_old ,m_inv_DeltaF, PPS_bar);
			//updating the correct principal pastic strain
			mPrincipalPlasticStrain_current[0] = /*mPrincipalPlasticStrain_old[0]*/ PPS_bar[0] + delta_lamda_a;
			mPrincipalPlasticStrain_current[1] = /*mPrincipalPlasticStrain_old[1]*/ PPS_bar[1] + delta_lamda_b; 
			mPrincipalPlasticStrain_current[2] = /*mPrincipalPlasticStrain_old[2]*/ PPS_bar[2] ;
			return true;
		      }
		      else
		        return false;
		              

              }

     
     void Isotropic_Rankine_Yield_Function::Three_Vector_Return_Mapping_To_Apex(const array_1d<double,3>& PrincipalStress, Vector& delta_lamda ,array_1d<double,3>& Sigma)
      {
          
              unsigned int iter    = 0;   
	      unsigned int max     = 10;
              int singular         = 0;  
              double norma         = 1.00;   
	      double delta_lamda_a = 0.00; 
              double delta_lamda_b = 0.00; 
              double delta_lamda_c = 0.00; 
              const  double toler  = 1E-3;
	      double E             = (*mprops)[YOUNG_MODULUS];
	      double NU            = (*mprops)[POISSON_RATIO];             
              double G             = 0.5*E / (1.00 + NU);
              double K             =  E / (3.00 * (1.00-2.00*NU));
              double H             =  mH;
	      Matrix d;             d.resize(3,3,false);      noalias(d)    = ZeroMatrix(3,3);  
              Matrix d_inv;         d_inv.resize(3,3,false);  noalias(d_inv)= ZeroMatrix(3,3); 
	      Vector residual      = ZeroVector(3);
              delta_lamda          = ZeroVector(3); 
	      
	      Vector Imput_Parameters(4);
	      Imput_Parameters    =  ZeroVector(4);
	      Imput_Parameters[0] =  mhe; 
	      Imput_Parameters[1] =  mrankine_accumulated_plastic_strain_current;
	      
	      
	      residual[0] = PrincipalStress[0] - mcurrent_Ft;
	      residual[1] = PrincipalStress[1] - mcurrent_Ft;
	      residual[2] = PrincipalStress[2] - mcurrent_Ft;
	      H           = mpSofteningBehaviorFt->FirstDerivateFunctionBehavior(Imput_Parameters);
	      
            
	      double Partial_Ep_gama_a = 0.00; 
	      double Partial_Ep_gama_b = 0.00; 
	      double Partial_Ep_gama_c = 0.00; 
	      
	      double prod_H_a = 0.00;
	      double prod_H_b = 0.00;
	      double prod_H_c = 0.00;
	      double Inc      = 0.00;
	      
	      
	      
              d.resize(3,3);
	      d_inv.resize(3,3);
	      const double raiz2d3    = 0.8164958092773;
 	      const double d3         = 0.3333333333333333;
	      double Ppvs = 0.00;                            /// principal plastic volumetric strain  
	      array_1d<double,3> Ppds = ZeroVector(3);       /// principal plastic desviatoric strain
	      array_1d<double,3> Pps  = ZeroVector(3);       /// principal plastic  strain
	      array_1d<double,3> I; 
	      I[0] = 1.00;
	      I[1] = 1.00;
	      I[2] = 1.00;
	      
              while(iter++<=max && norma>= toler) 
		  {
		    
		      prod_H_a = H * Partial_Ep_gama_a;
		      prod_H_b = H * Partial_Ep_gama_b;
		      prod_H_c = H * Partial_Ep_gama_c;
		      
		    
		      d(0,0) = -4.00 * G / 3.00 - K - prod_H_a;  d(0,1)  =  2.00 * G / 3.00 - K - prod_H_b;  d(0,2)   =   2.00 * G / 3.00 - K  - prod_H_c; 
		      d(1,0) =  2.00 * G / 3.00 - K - prod_H_a;  d(1,1)  = -4.00 * G / 3.00 - K - prod_H_b;  d(1,2)   =   2.00 * G / 3.00 - K  - prod_H_c; 
		      d(2,0) =  2.00 * G / 3.00 - K - prod_H_a;  d(2,1)  =  2.00 * G / 3.00 - K - prod_H_b;  d(2,2)   =  -4.00 * G / 3.00 - K  - prod_H_c;   
	      
	              singular =  SD_MathUtils<double>::InvertMatrix(d, d_inv);   
		      noalias(delta_lamda) =  delta_lamda - Vector(prod(d_inv, residual)); 
		
                        
		      delta_lamda_a = delta_lamda[0];
		      delta_lamda_b = delta_lamda[1];
                      delta_lamda_c = delta_lamda[2];
                     
		      /// normal
		      ///mrankine_accumulated_plastic_strain_current= mrankine_accumulated_plastic_strain_old+ delta_lamda_a + delta_lamda_b +  delta_lamda_c;
		     
		      ///von mises
		      Inc = 0.81649658092773 * norm_2(delta_lamda);
		      mrankine_accumulated_plastic_strain_current = mrankine_accumulated_plastic_strain_old + Inc;  
		      
		      
		      /// computing Kp_punto.
		      /// volumetric and desviatoric plastic strain
		      Pps[0]        = delta_lamda_a;
		      Pps[1]        = delta_lamda_b;
		      Pps[2]        = delta_lamda_c;
		      Ppvs          = Pps[0] + Pps[1] + Pps[2]; 
		      noalias(Ppds) = Pps - d3 * Ppvs * I;
		      Inc           = raiz2d3* std::sqrt(inner_prod(Pps,Pps));
		      mrankine_accumulated_plastic_strain_current = mrankine_accumulated_plastic_strain_old + Inc;
		      ComputeActualStrees(Ppvs, Ppds, PrincipalStress, Sigma);
		      ComputeActualStrain(Pps);
		      CalculatePlasticDamage(Sigma);
		      
		      
                       /// Updatinf mFt
		       Partial_Ep_gama_a   = (2.00/3.00) * delta_lamda_a/Inc;
		       Partial_Ep_gama_b   = (2.00/3.00) * delta_lamda_b/Inc; 
		       Partial_Ep_gama_b   = (2.00/3.00) * delta_lamda_c/Inc; 
		       
		       Imput_Parameters[1] =  mrankine_accumulated_plastic_strain_current;
		       Imput_Parameters[2] =  mpastic_damage_old;
                       Imput_Parameters[3] =  mpastic_damage_current;
		       mcurrent_Ft         =  mpSofteningBehaviorFt->FunctionBehavior(Imput_Parameters);  
		       H                   =  mpSofteningBehaviorFt->FirstDerivateFunctionBehavior(Imput_Parameters);  
                      
		       
                       /// ft se anulan
                      if(mcurrent_Ft<=toler)
                        {
                            mcurrent_Ft = 0.00;
                            //delta_lamda_a  = delta_lamda[0];  
                            //delta_lamda_b  = delta_lamda[1];    
                            //delta_lamda_c  = delta_lamda[2];                             
                            break;
                        }
                       else 
                       {
                                                
                  
	            residual[0] = PrincipalStress[0] - mcurrent_Ft;
                    residual[1] = PrincipalStress[1] - mcurrent_Ft;
	            residual[2] = PrincipalStress[2] - mcurrent_Ft;
                      
                    residual[0] =  residual[0] - delta_lamda_a*( 4.00 * G / 3.00 + K ) - delta_lamda_b*( -2.00 * G / 3.00 + K ) - delta_lamda_c*( -2.00 * G / 3.00 + K );   
                    residual[1] =  residual[1] - delta_lamda_a*(-2.00 * G / 3.00 + K ) - delta_lamda_b*( 4.00  * G / 3.00 + K ) - delta_lamda_c*( -2.00 * G / 3.00 + K );   
                    residual[2] =  residual[2] - delta_lamda_a*(-2.00 * G / 3.00 + K ) - delta_lamda_b*( -2.00 * G / 3.00 + K ) - delta_lamda_c*( 4.00  * G / 3.00 + K );         
                    norma       =  norm_2(residual); 

                        }
                      }

                   if(iter>=max || delta_lamda_a<0.0 || delta_lamda_b<0.0 || delta_lamda_c<0){
		      KRATOS_WATCH(iter)
		      KRATOS_WATCH(norma)
		      KRATOS_WATCH(mcurrent_Ft)
		      KRATOS_WATCH(Sigma)
		      KRATOS_WATCH(PrincipalStress)
		      KRATOS_WATCH(delta_lamda)
		      std::cout<< "RETURN MAPPING APEX RANKINE  NOT CONVERGED" << std::endl;
		      KRATOS_ERROR(std::logic_error,  "RETURN MAPPING APEX RANKINE  NOT CONVERGED" , "");
		   }

                     if(mcurrent_Ft <= 0.00) 
                        {
                           Sigma[0] = 0.00;
                           Sigma[1] = 0.00;
                           Sigma[2] = 0.00;
                        }
                     else 
                      {
			/// volumetric and desviatoric plastic strain
			Pps[0]        = delta_lamda_a;
			Pps[1]        = delta_lamda_b;
			Pps[2]        = delta_lamda_c;
			Ppvs          = Pps[0] + Pps[1] + Pps[2]; 
			noalias(Ppds) = Pps - d3 * Ppvs * I;  
			noalias(Sigma) = PrincipalStress - 2.00 * G * Ppds - K * Ppvs * I;

			
		      //Sigma[0] = PrincipalStress[0] - delta_lamda_a*( 4.00  * G / 3.00 + K ) - delta_lamda_b*( -2.00  * G / 3.00 + K ) - delta_lamda_c*( -2.00  * G / 3.00 + K );
		      //Sigma[1] = PrincipalStress[1] - delta_lamda_a*(-2.00 * G  / 3.00 + K ) - delta_lamda_b*( 4.00  * G / 3.00 + K )  - delta_lamda_c*( -2.00  * G / 3.00 + K );
                      //Sigma[2] = PrincipalStress[2] - delta_lamda_a*(-2.00 * G  / 3.00 + K ) - delta_lamda_b*( -2.00  * G / 3.00 + K )  - delta_lamda_c*( 4.00  * G / 3.00 + K ); 

		      }
                       
                     Vector PPS_bar(3);
		     Imput_Parameters = ZeroVector(3);
                     ComputePlasticStrainBar(mplastic_strain_old ,m_inv_DeltaF, PPS_bar);
                     mPrincipalPlasticStrain_current[0] = /*mPrincipalPlasticStrain_old[0]*/ PPS_bar[0] +  delta_lamda_a;
	             mPrincipalPlasticStrain_current[1] = /*mPrincipalPlasticStrain_old[0]*/ PPS_bar[1] +  delta_lamda_b; 
		     mPrincipalPlasticStrain_current[2] = /*mPrincipalPlasticStrain_old[0]*/ PPS_bar[2] +  delta_lamda_c; 
      }
     
         
bool Isotropic_Rankine_Yield_Function::CheckValidity( array_1d<double,3>& Sigma)
      {
	    bool check   = false;
	    array_1d<double,3> Aux_Sigma = ZeroVector(3);
	    Aux_Sigma[0] = fabs(Sigma[0]);
	    Aux_Sigma[1] = fabs(Sigma[1]);
	    Aux_Sigma[2] = fabs(Sigma[2]);
	    double delta = (*max_element(Aux_Sigma.begin(), Aux_Sigma.end())) * 1.00E-6; 
	    if( (Sigma[0] + delta) >= Sigma[1] && (Sigma[1] + delta ) >= Sigma[2]){ check = true;}
	    return check;
      }  

 void Isotropic_Rankine_Yield_Function::GetValue(const Variable<Matrix>& rVariable, Matrix& Result)
 {
     
   if (rVariable==GREEN_LAGRANGE_PLASTIC_STRAIN_TENSOR) 
   {
     int size = SizePlasticStrain();
     Result.resize(1, size, false);
     for(int i = 0; i< size; i++ )
        Result(0,i) = mplastic_strain(i); 
   }
    return;
   }
   
   void Isotropic_Rankine_Yield_Function::GetValue(double& Result)
   {
     mhe = Result;   
     return;
   }
   
   void Isotropic_Rankine_Yield_Function::GetValue(Matrix& Result)
   {
	m_inv_DeltaF.resize(3,3, false);
	noalias(m_inv_DeltaF) = ZeroMatrix(3,3);
	switch(mState)
         {
	  case Plane_Stress:
          {
	    KRATOS_ERROR(std::logic_error,  "PLANE STRESS NOT IMPLEMENTED" , "");
	    break;
          } 
          case Plane_Strain:
            {
	      m_inv_DeltaF(0,0)    = Result(0,0);
	      m_inv_DeltaF(0,1)    = Result(0,1);
	      m_inv_DeltaF(1,0)    = Result(1,0);
	      m_inv_DeltaF(1,1)    = Result(1,1);
	      m_inv_DeltaF(2,2)    = 1.00;
	      break;
            }
            
	  case Tri_D:
            {
	      noalias(m_inv_DeltaF) = Result;
	      break;
            }
           
	 }
    }
 
void Isotropic_Rankine_Yield_Function::GetValue(const Variable<Vector>& rVariable, Vector& Result)
    {
      const int size = mplastic_strain.size();
      if(rVariable==ALMANSI_PLASTIC_STRAIN){
	    Result.resize(size);
	    noalias(Result) = mplastic_strain;
          }
      if(rVariable==ALMANSI_ELASTIC_STRAIN){
	    Result.resize(size);
	    noalias(Result) = mElastic_strain ;
          }  
    }
   
   void Isotropic_Rankine_Yield_Function::GetValue(const Variable<double>& rVariable, double& Result)
      {
	
          //const double& Ft   = (*mprops)[FT];
	  //const double& Ec   = (*mprops)[YOUNG_MODULUS];
	  //const double& GE   = (*mprops)[FRACTURE_ENERGY];
	  //const double Eu    =  (2.00 * GE)/(Ft * mhe);
	  
	if(rVariable==COHESION)
	  Result = 0.00;
        
        if(rVariable==FT){
	   Result = mcurrent_Ft;
	}
	
	if(rVariable == DILATANCY_ANGLE) /// NOT APPLY
          Result = 0.00; 
	
	if(rVariable == INTERNAL_FRICTION_ANGLE) /// NOT APPLY
	  Result = 0.00;
	
        if(rVariable == DAMAGE)
	   Result = mpastic_damage_current; 
      }
      

    void Isotropic_Rankine_Yield_Function::CalculatePlasticDamage(const array_1d<double,3>& Sigma)
        {  
	  /*
	  //double shi             =  0.81649658092773 * (mrankine_accumulated_plastic_strain_current- mrankine_accumulated_plastic_strain_old);
	  array_1d<double, 3> DeltaPlasticStrain = ZeroVector(3);
	  noalias(DeltaPlasticStrain) = mPrincipalPlasticStrain_current - mPrincipalPlasticStrain_old;
	  double disipation      =  inner_prod(Sigma, DeltaPlasticStrain);  
	  double gf              =  (*mprops)[FRACTURE_ENERGY]/mhe;
	  double kp_punto        =  disipation; //UniaxialTension(Sigma) * shi; /// son igales fijate en el factor 1.22 es raiz 1.5
	  kp_punto               =  kp_punto/gf; 
	  mpastic_damage_current =  mpastic_damage_old + 0.816496581 * kp_punto;  
	  */
	  const double toler = 1E-10;
	  double teta_a     =  Tensor_Utils<double>::Mc_aully(Sigma);
          double teta_b     =  std::fabs(Sigma[0]) + std::fabs(Sigma[1]) + std::fabs(Sigma[2]);
	  double teta       =  0.00;
          array_1d<double, 3> DeltaPlasticStrain = ZeroVector(3);
          noalias(DeltaPlasticStrain) = mPrincipalPlasticStrain_current - mPrincipalPlasticStrain_old;
	  double disipation  =  inner_prod(Sigma, DeltaPlasticStrain);
	  
	  if (fabs(teta_b) > toler)
          {
	   teta = teta_a/teta_b;
	   // computing Kp_punto
	   double gf_p            = (*mprops)[FRACTURE_ENERGY]/mhe;
           double h               = (teta/gf_p);
	   double kp_punto        = h * disipation;    
	   if(disipation > 0.00)
	   {
	     mpastic_damage_current =  mpastic_damage_old + kp_punto;
	     if(mpastic_damage_current > 1.00)
	       mpastic_damage_current = 1.00;
	   }
	  }
	  
	}

//***********************************************************************
//************************************************* ********************** 

void Isotropic_Rankine_Yield_Function::ComputeActualStrees(const double& Ppvs, 
						      const array_1d<double,3>& Ppds,
						      const array_1d<double,3>& PrincipalStress,
						      array_1d<double,3>& Sigma)
{
   const double& Young   = (*mprops)[YOUNG_MODULUS];
   const double& Poisson = (*mprops)[POISSON_RATIO];
   const double G        = Young/(2.00 * (1.00 + Poisson) );
   const double K        = Young/(3.00 * (1.00-2.00*Poisson)); 
   array_1d<double,3> I; 
   I[0] = 1.00;
   I[1] = 1.00;
   I[2] = 1.00;
   noalias(Sigma) = PrincipalStress - 2.00 * G * Ppds - K * Ppvs * I;
}    

//***********************************************************************
//*********************************************************************** 

void Isotropic_Rankine_Yield_Function::ComputeActualStrain(const array_1d<double,3>& Pps)
{
  Vector PPS_bar(3);
  PPS_bar = ZeroVector(3);
  ComputePlasticStrainBar(mplastic_strain_old ,m_inv_DeltaF, PPS_bar);  
  noalias(mPrincipalPlasticStrain_current) = PPS_bar + Pps;
}

//***********************************************************************
//*********************************************************************** 
double Isotropic_Rankine_Yield_Function::UniaxialTension(const Vector& Stress)
        {
	   return Stress[0];
	}

}


