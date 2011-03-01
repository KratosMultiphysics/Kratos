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



#include "fluency_criteria/modified_morh_coulomb_yield_function.h"

namespace Kratos
  {


	    
            Modified_Morh_Coulomb_Yield_Function::Modified_Morh_Coulomb_Yield_Function()
            {
	    }
           
            Modified_Morh_Coulomb_Yield_Function::Modified_Morh_Coulomb_Yield_Function(
            myState State, 
            const MorhCoulombPointerType MorhCoulomb, const RankinePointerType Rankine)
	    :FluencyCriteria()
	    {
	       mState       = State; 
               mMorhCoulomb = MorhCoulomb;
               mRankine     = Rankine;
	      
	    }
              

             Modified_Morh_Coulomb_Yield_Function::~Modified_Morh_Coulomb_Yield_Function() {}

             void Modified_Morh_Coulomb_Yield_Function::InitializeMaterial(const Properties& props)
             {
	       mprops =  &props;
	       mMorhCoulomb->InitializeMaterial(props);
	       mRankine->InitializeMaterial(props);  
	       
	       m_modified_morh_coulomb_maccumulated_plastic_strain_old         = 0.00;
               m_modified_morh_coulomb_maccumulated_plastic_strain_current     = 0.00;
	       mpastic_damage_old                                              = 0.00;
               mpastic_damage_current                                          = 0.00;
	       mplastic_strain                                                 = ZeroVector(4);
	       mplastic_strain_old                                             = ZeroVector(4); 
	       mPrincipalPlasticStrain_current                                 = ZeroVector(3);
	       mPrincipalPlasticStrain_old                                     = ZeroVector(3);
	     }
             
              void Modified_Morh_Coulomb_Yield_Function::GetValue(double Result)
              {
		mlength = Result;
	      }
             
             void Modified_Morh_Coulomb_Yield_Function::ReturnMapping(const Vector& StrainVector, Vector& StressVector)
             {
	       
	             
		Vector Aux_Trial_Stress;
		Aux_Trial_Stress.resize(StressVector.size());
		noalias(Aux_Trial_Stress) =  StressVector;

		array_1d<double,3> PrincipalStress;
		array_1d<array_1d < double,3 > ,3> EigenVectors;
		array_1d<unsigned int,3> Order;
		
		SpectralDecomposition(StressVector, PrincipalStress, EigenVectors);
		IdentifyMaximunAndMinumumPrincipalStres_CalculateOrder(PrincipalStress, Order);
	        
		const array_1d<double,3> Aux_PrincipalStress              = PrincipalStress; 
		const array_1d<array_1d < double,3 > ,3> Aux_EigenVectors = EigenVectors; 
		const array_1d<unsigned int,3>  Aux_Order                 = Order;
		
		bool plastic_1 = mMorhCoulomb->CheckPlasticAdmisibility(PrincipalStress); 
		bool plastic_2 = mRankine->CheckPlasticAdmisibility(PrincipalStress); 
		        
                // elastic step
		if(plastic_1==false && plastic_2==false)
		  return;
		

		// standar morh coulomb
		if(plastic_1==true)
		{ 
		  mMorhCoulomb->ReturnMapping(StrainVector, StressVector);
		  
		  //check tensile consistency
		  SpectralDecomposition(StressVector, PrincipalStress, EigenVectors);
		  IdentifyMaximunAndMinumumPrincipalStres_CalculateOrder(PrincipalStress, Order);
		  
		  bool tensile_consistency =  mRankine->CheckPlasticAdmisibility(PrincipalStress);
		  if(tensile_consistency==true) // sigue siendo plastico 
		  { 
		    //update variables
		    mMorhCoulomb->UpdateMaterial();
	            mRankine->UpdateMaterial();
		    
		    // Cero value is false 1 is equal to one 
		    array_1d<double, 3> Sigma;  // Principal Stress updated
		    noalias(StressVector)    = Aux_Trial_Stress;
		    noalias(PrincipalStress) = Aux_PrincipalStress;
		    if(Return_Mapping_Intersection_Main_Plane_And_Sigma1_Tesile_Plane(PrincipalStress , Order,  Sigma)==false)
		    {
		      //return mapping to corner
		      mMorhCoulomb->UpdateMaterial();
	              mRankine->UpdateMaterial();
		      array_1d<double, 3> T;
		      double sinpsi  = std::sin(PI * mMorhCoulomb->mcurrent_dilatancy_angle  / 180.00);
		      const double a = 3.00 + sinpsi;
		      const double b = 3.00 - sinpsi;
		      T[0] =  2.00/a - 1.00/b;
		      T[1] = -1.00/a - 1.00/b;
		      T[2] = -1.00/a - 2.00/b;
		     
		      bool rigth_corner = false;  
		      if(inner_prod(T,Aux_PrincipalStress)>0)
			 rigth_corner = true;
		      
		      if( rigth_corner==true) //Main Plane Corner and Sigma 1
			Return_Mapping_Intersection_Main_Plane_Corner_And_Sigma1_Tesile_Plane(PrincipalStress, Order, Sigma);
		      else //Main Plane Corner and Sigma 1 and Sigma 2
			Return_Mapping_Intersection_Main_Plane_Corner_And_Sigma1_And_Sigma_2_Tesile_Plane(PrincipalStress, Order, Sigma);  
		    }
		    
		    AssembleUpdateStressAndStrainTensor(Sigma,  EigenVectors,  Order, StrainVector, StressVector);
		  }
		  
		  else
		    {
		      mPrincipalPlasticStrain_current = mPrincipalPlasticStrain_old + 
		      (mMorhCoulomb-> mPrincipalPlasticStrain_current - mMorhCoulomb-> mPrincipalPlasticStrain_old);
		      
		      m_modified_morh_coulomb_maccumulated_plastic_strain_current = m_modified_morh_coulomb_maccumulated_plastic_strain_old +
		      (mMorhCoulomb-> mmorh_coulomb_maccumulated_plastic_strain_current - mMorhCoulomb-> mmorh_coulomb_maccumulated_plastic_strain_old);
		    }
		}
		else
		{
		  mRankine->ReturnMapping(StrainVector, StressVector);   
		  mPrincipalPlasticStrain_current =  mPrincipalPlasticStrain_old 
		  + ( mRankine-> mPrincipalPlasticStrain_current - mRankine-> mPrincipalPlasticStrain_old);
		  
		  m_modified_morh_coulomb_maccumulated_plastic_strain_current =  m_modified_morh_coulomb_maccumulated_plastic_strain_old + 
		  (mRankine->mrankine_accumulated_plastic_strain_current - mRankine->mrankine_accumulated_plastic_strain_old);

		}
		
	        // Updating the plastic strain tensor
                CalculatePlasticStrain(StrainVector,StressVector);
		
		// Calculating plastic damage
	        SpectralDecomposition(StressVector, PrincipalStress, EigenVectors);
		IdentifyMaximunAndMinumumPrincipalStres_CalculateOrder(PrincipalStress, Order);
		CalculatePlasticDamage(PrincipalStress);
		
	     }
	     
             void Modified_Morh_Coulomb_Yield_Function::FinalizeSolutionStep()
             { 
	       mMorhCoulomb->FinalizeSolutionStep();
	       mRankine->FinalizeSolutionStep();
	       
	       m_modified_morh_coulomb_maccumulated_plastic_strain_old =  m_modified_morh_coulomb_maccumulated_plastic_strain_current;
	       noalias(mPrincipalPlasticStrain_old)                    =  mPrincipalPlasticStrain_current; 
	       noalias(mplastic_strain_old)                            =  mplastic_strain; 
	       mpastic_damage_old                                      =  mpastic_damage_current;
	       
	     }
	     
             void Modified_Morh_Coulomb_Yield_Function::UpdateMaterial()
             { 
	       mMorhCoulomb->UpdateMaterial();
	       mRankine->UpdateMaterial();
	       
               m_modified_morh_coulomb_maccumulated_plastic_strain_current =  m_modified_morh_coulomb_maccumulated_plastic_strain_old;  
	       noalias(mPrincipalPlasticStrain_current)                    =  mPrincipalPlasticStrain_old;  
	       noalias(mplastic_strain)                                    =  mplastic_strain_old; 
               mpastic_damage_current                                      =  mpastic_damage_old;  
	     } 
	     
	     
	     bool Modified_Morh_Coulomb_Yield_Function::Return_Mapping_Intersection_Main_Plane_And_Sigma1_Tesile_Plane(const array_1d<double,3>& PrincipalStress, const array_1d<unsigned int,3>& order, array_1d<double,3>& Sigma)
	     {
	        
		const double& E         =   (*mprops)[YOUNG_MODULUS];
		const double& NU        =   (*mprops)[POISSON_RATIO];
		const double G          =   0.5 * E / (1.00 + NU);
		const double K          =   E / (3.00 * (1.00-2.00 * NU) );
		const double toler      =   1E-8;
		const unsigned max      =   1000;
		unsigned int  iter      = 0;

		
		array_1d<double, 2> dgama    = ZeroVector(2);     
		array_1d<double, 2> ddgama   = ZeroVector(2); 
		array_1d<double, 2> residual = ZeroVector(2);
		Matrix d                    = ZeroMatrix(2,2);
                Matrix d_inv                = ZeroMatrix(2,2);
		
		

		
	        double sinphi             =   std::sin(PI * mMorhCoulomb->minternal_friction_angle  / 180.00);
		double cosphi             =   std::cos(PI * mMorhCoulomb->minternal_friction_angle  / 180.00);
		double sinpsi             =   std::sin(PI * mMorhCoulomb->mcurrent_dilatancy_angle  / 180.00);
		double cospsi             =   std::cos(PI * mMorhCoulomb->mcurrent_dilatancy_angle  / 180.00); 
		
		residual[0]   = (PrincipalStress[0] - PrincipalStress[2]) + (PrincipalStress[0] + PrincipalStress[2]) * sinphi - 2.00 *  cosphi * mMorhCoulomb->mcurrent_cohesion;
		residual[1]   =  PrincipalStress[0]  - mRankine->mcurrent_Ft;

		
		double Partial_Cohesion  = 0.00;   
		double Partial_Friction  = 0.00;
		double Partial_Dilatancy = 0.00;
		
		
 
	        double A = 0.00; 
                double B = 0.00; 
                double C = 0.00; 
                double D = 0.00; 
                
		
		double Partial_A_gama_a  = 0.00; 
                double Partial_A_gama_b  = 0.00;
	        double Partial_B_gama_a  = 0.00; 
                double Partial_B_gama_b  = 0.00;
                double Partial_C_gama_a  = 0.00; 
                double Partial_C_gama_b  = 0.00;
                double Partial_D_gama_a  = 0.00; 
                double Partial_D_gama_b  = 0.00;

		double Partial_Ep_gama_a = 0.00; 
                double Partial_Ep_gama_b = 0.00;

                double norma = norm_2(residual);
		int singular = 0;  
		const double raiz2d3 = 0.8164958092773; 
		double aux_1 = 0.00; 
		double aux_2 = 0.00;   
		double aux   = 1.00;
		double& gama_a = dgama[0];
		double& gama_b = dgama[1];
		
		while(fabs(norma)>toler && iter++ < max )
		{  

		Partial_Cohesion   = 0.00;   //  WARNING->deberia cambiar valor
		Partial_Friction   = 0.00;   //  WARNING->deberia cambiar valor
		Partial_Dilatancy  = 0.00;   //  WARNING->deberia cambiar valor
                 
		 
		Partial_Ep_gama_a  = (2.00/3.00) * ( gama_a * ( (1.00 + sinpsi) * (1.00 + sinpsi) + (-1.00 + sinpsi) * (-1.00 + sinpsi) )  + gama_b * (1.00 + sinpsi) ); 
		Partial_Ep_gama_a  = Partial_Ep_gama_a / aux;
		
		Partial_Ep_gama_b  = (2.00/3.00) * ( gama_a * (1.00 + sinpsi)  + gama_b ); 
		Partial_Ep_gama_b  = Partial_Ep_gama_b /  aux;


                Partial_A_gama_a   = 1.00 + (1.00/3.00) * sinpsi + (1.00/3.00) * cospsi * gama_a * Partial_Dilatancy * Partial_Ep_gama_a; 
		Partial_A_gama_b   = (1.00/3.00) * gama_a * cospsi * Partial_Dilatancy * Partial_Ep_gama_b + 2.00/3.00;
		
                Partial_B_gama_a   = (-2.00/3.00) * sinpsi - (2.00/3.00) * gama_a * cospsi * Partial_Dilatancy * Partial_Ep_gama_a;
		Partial_B_gama_b   = (-2.00/3.00) * gama_a *cospsi * Partial_Dilatancy * Partial_Ep_gama_b - 1.00/3.00;
		
		Partial_C_gama_a   = (1.00/3.00) * sinpsi -1.00 + (1.00/3.00) * gama_a * cospsi * Partial_Dilatancy * Partial_Ep_gama_a;
		Partial_C_gama_b   = (1.00/3.00) * gama_a * cospsi * Partial_Dilatancy * Partial_Ep_gama_b - 1.00/3.00;
		
	        Partial_D_gama_a   = 2.00 * sinpsi + 2.00 * gama_a * cospsi * Partial_Dilatancy * Partial_Ep_gama_a;
		Partial_D_gama_b   = 2.00 * gama_a * cospsi * Partial_Dilatancy * Partial_Ep_gama_a + 1.00;
		
		d(0,0) = 
		  
		  ( PrincipalStress[0] + PrincipalStress[2] ) * cosphi * Partial_Friction * Partial_Ep_gama_a
		  - 2.00 * cosphi * Partial_Cohesion * Partial_Ep_gama_a
		  + 2.00 * sinphi * (mMorhCoulomb->mcurrent_cohesion) * Partial_Friction * Partial_Ep_gama_a
		  + 2.00 * G * Partial_C_gama_a 
		  - 2.00 * G * sinphi *  Partial_A_gama_a 
		  - 2.00 * G * cosphi * A * Partial_Friction * Partial_Ep_gama_a
		  - 2.00 * K * cosphi * D * Partial_Friction * Partial_Ep_gama_a
		  - 2.00 * G * sinphi * Partial_C_gama_a 
		  - 2.00 * G * cosphi * C * Partial_Friction * Partial_A_gama_a 
		  - 2.00  * K * sinphi * Partial_D_gama_a 
		  - 2.00 * G * Partial_A_gama_a; 
		  
		 d(0,1) = 
		  
		  ( PrincipalStress[0] + PrincipalStress[2] ) * cosphi * Partial_Friction * Partial_Ep_gama_b
		  - 2.00 * cosphi * Partial_Cohesion * Partial_Ep_gama_b
		  + 2.00 * sinphi * (mMorhCoulomb->mcurrent_cohesion) * Partial_Friction * Partial_Ep_gama_b
		  + 2.00 * G * Partial_C_gama_b 
		  - 2.00 * G * sinphi *  Partial_A_gama_b 
		  - 2.00 * G * cosphi * A * Partial_Friction * Partial_Ep_gama_b
		  - 2.00 * K * cosphi * D * Partial_Friction * Partial_Ep_gama_b
		  - 2.00 * G * sinphi * Partial_C_gama_b 
		  - 2.00 * G * cosphi * C * Partial_Friction * Partial_A_gama_b 
		  - 2.00  * K * sinphi * Partial_D_gama_b - 2.00 * G * Partial_A_gama_b;   
		  
		  
		  d(1,0) = -2.00 * G * Partial_A_gama_a - K * Partial_D_gama_a;  
		  
		  d(1,1) = -2.00 * G * Partial_A_gama_b - K * Partial_D_gama_b;  
		  
		  singular       =  SD_MathUtils<double>::InvertMatrix(d, d_inv);
		  ddgama         =  -Vector(prod(d_inv, residual));
  
		  //Compute Newton-Raphson increment and update variables DGAMA and DGAMB
		  noalias(dgama) += ddgama;  
		  
		 //Updating Materiales
		 aux_1   = (sinpsi+1.00) * gama_a + gama_b; aux_1 = aux_1 * aux_1;  
		 aux_2   = (sinpsi-1.00) * gama_a;          aux_2 = aux_2 * aux_2;  
		 aux     = raiz2d3 * std::sqrt(aux_1 + aux_2);
		 
		 m_modified_morh_coulomb_maccumulated_plastic_strain_current = m_modified_morh_coulomb_maccumulated_plastic_strain_old + aux;

		 ///WARNING UPDATE INTERNAL VARIABLES
                 /*
                   cohesion;
		   dilatancy;
		   friction;
		   ft;                           
                 */
		 
		sinphi             =   std::sin(PI * mMorhCoulomb->minternal_friction_angle  / 180.00);
		cosphi             =   std::cos(PI * mMorhCoulomb->minternal_friction_angle  / 180.00);
		sinpsi             =   std::sin(PI * mMorhCoulomb->mcurrent_dilatancy_angle  / 180.00);
		cospsi             =   std::cos(PI * mMorhCoulomb->mcurrent_dilatancy_angle  / 180.00); 
		
		A = gama_a * (1.00 + (1.00/3.00) * sinpsi) + (2.00/3.00) * gama_b;
		B = -(2.00/3.00) * gama_a * sinpsi - (1.00/3.00) * gama_b;                 
		C = gama_a * ( (1.00/3.00) *sinpsi  - 1.00 ) - (1.00/3.00) * gama_b;
		D = 2.00 * gama_a * sinpsi + gama_b;
		 
		residual[0]   = (PrincipalStress[0] -  PrincipalStress[2]) + (PrincipalStress[0] + PrincipalStress[2]) * sinphi - 2.00 *  cosphi * mMorhCoulomb->mcurrent_cohesion;
		residual[1]   =  PrincipalStress[0]  - mRankine->mcurrent_Ft;
		
	        residual[0]   = residual[0]-2.00 * G * A + 2.00 * G * C - 2.00 * G * sinphi * A - 2.00 * G * sinphi * C - 2.00 * K * sinphi * D;  
		residual[1]   = residual[1]-2.00 * G * A - K * D; 
		norma         = norm_2(residual); 
		 
		
		if(iter>=max){  
		  KRATOS_WATCH(norma)
		  KRATOS_WATCH(residual)
		  KRATOS_WATCH(dgama) 
		  KRATOS_ERROR(std::logic_error,  "RETURN MAPPING OF MAIN PLANE AND SIGMA 1 NOT CONVERGED" , " "); 
		 }
		}


		Sigma[0] = PrincipalStress[0] -  2.00 * G * A - K * D; 
                Sigma[1] = PrincipalStress[1] -  2.00 * G * B - K * D; 
                Sigma[2] = PrincipalStress[2] -  2.00 * G * C - K * D; 
		
		
		
		bool cond_a = mRankine->CheckValidity(Sigma);
		bool cond_b = mRankine->CheckPlasticAdmisibility(Sigma);     // false when no plastic
		bool cond_c = mMorhCoulomb->CheckPlasticAdmisibility(Sigma); // false when no plastic
		
		bool cond_d = false;
// 		if( dgama[0]< 0.00 || dgama[1]<0.00)
// 		     cond_d = true;
		
		
		if(cond_a==true && cond_b==false && cond_c ==false && cond_d==false)
		{
		    //updating the correct principal pastic strain
		    mPrincipalPlasticStrain_current[0] = mPrincipalPlasticStrain_old[0] + gama_a * (1.00 + sinpsi) + gama_b;
		    mPrincipalPlasticStrain_current[1] = mPrincipalPlasticStrain_current[1]; 
		    mPrincipalPlasticStrain_current[2] = mPrincipalPlasticStrain_old[2] + gama_a * (sinpsi-1.00 ); 
		    return true;
		}
		else
		  return false;
		
		}
		
		
    bool Modified_Morh_Coulomb_Yield_Function::Return_Mapping_Intersection_Main_Plane_Corner_And_Sigma1_Tesile_Plane(const array_1d<double,3>& PrincipalStress, const array_1d<unsigned int,3>& order, array_1d<double,3>& Sigma)
    {
                const double& E              =   (*mprops)[YOUNG_MODULUS];
		const double& NU             =   (*mprops)[POISSON_RATIO];
		const double G               =   0.5 * E / (1.00 + NU);
		const double K               =   E / (3.00 * (1.00-2.00 * NU) );
		const double toler           =  1E-9;
		const unsigned max           =  100;
		
		array_1d<double, 3> dgama    = ZeroVector(3);     
		array_1d<double, 3> ddgama   = ZeroVector(3); 
		array_1d<double, 3> residual = ZeroVector(3);
		Matrix d                     = ZeroMatrix(3,3);
                Matrix d_inv                 = ZeroMatrix(3,3);
		
		unsigned iter                =  0;
	        double sinphi                =   std::sin(PI * mMorhCoulomb->minternal_friction_angle  / 180.00);
		double cosphi                =   std::cos(PI * mMorhCoulomb->minternal_friction_angle  / 180.00);
		double sinpsi                =   std::sin(PI * mMorhCoulomb->mcurrent_dilatancy_angle  / 180.00);
		double cospsi                =   std::cos(PI * mMorhCoulomb->mcurrent_dilatancy_angle  / 180.00); 
		
		residual[0]   = (PrincipalStress[0] - PrincipalStress[2]) + (PrincipalStress[0] + PrincipalStress[2]) * sinphi - 2.00 *  cosphi * mMorhCoulomb->mcurrent_cohesion;
		residual[1]   = (PrincipalStress[0] - PrincipalStress[1]) + (PrincipalStress[0] + PrincipalStress[1]) * sinphi - 2.00 *  cosphi * mMorhCoulomb->mcurrent_cohesion;
		residual[2]   =  PrincipalStress[0] - mRankine->mcurrent_Ft;
		
		double Partial_Cohesion  = 0.00;   
		double Partial_Friction  = 0.00;
		double Partial_Dilatancy = 0.00;
 
	        double A = 0.00; 
                double B = 0.00; 
                double C = 0.00; 
                double D = 0.00; 
                
		
		double Partial_A_gama_a  = 0.00; 
                double Partial_A_gama_b  = 0.00;
		double Partial_A_gama_c  = 0.00;
		
	        double Partial_B_gama_a  = 0.00; 
                double Partial_B_gama_b  = 0.00;
		double Partial_B_gama_c  = 0.00;
		
                double Partial_C_gama_a  = 0.00; 
                double Partial_C_gama_b  = 0.00;
		double Partial_C_gama_c  = 0.00;
		
                double Partial_D_gama_a  = 0.00; 
                double Partial_D_gama_b  = 0.00;
		double Partial_D_gama_c  = 0.00;

		double Partial_Ep_gama_a = 0.00; 
                double Partial_Ep_gama_b = 0.00;
		double Partial_Ep_gama_c = 0.00;

                double norma = norm_2(residual);
		int singular = 0;  
		const double raiz2d3 = 0.8164958092773; 
		double aux_1 = 0.00; 
		double aux_2 = 0.00;   
		double aux   = 1.00;
		double& gama_a = dgama[0];
		double& gama_b = dgama[1];
		double& gama_c = dgama[2];
		
		while(fabs(norma)>toler && iter++ < max )
		{
		  
		  Partial_Cohesion   = 0.00;   //  WARNING->deberia cambiar valor
		  Partial_Friction   = 0.00;   //  WARNING->deberia cambiar valor
		  Partial_Dilatancy  = 0.00;   //  WARNING->deberia cambiar valor

		  sinphi             =   std::sin(PI * mMorhCoulomb->minternal_friction_angle  / 180.00);
		  cosphi             =   std::cos(PI * mMorhCoulomb->minternal_friction_angle  / 180.00);
		  sinpsi             =   std::sin(PI * mMorhCoulomb->mcurrent_dilatancy_angle  / 180.00);
		  cospsi             =   std::cos(PI * mMorhCoulomb->mcurrent_dilatancy_angle  / 180.00); 


		  Partial_Ep_gama_a  = (2.00/3.00) *( ( ( gama_a + gama_b )*(1.00 + sinpsi) + gama_c ) * (1.00 + sinpsi) +  gama_a * ( -1.00 + sinpsi) * (-1.00 + sinpsi) ) ; 
		  Partial_Ep_gama_a  = Partial_Ep_gama_a / aux;

		  Partial_Ep_gama_b  = (2.00/3.00) *( ( ( gama_a + gama_b )*(1.00 + sinpsi) + gama_c ) * (1.00 + sinpsi) +  gama_b * ( -1.00 + sinpsi) * (-1.00 + sinpsi) ) ; 
		  Partial_Ep_gama_b  = Partial_Ep_gama_b /  aux;
		  
		  Partial_Ep_gama_c  = (2.00/3.00) * ( ( gama_a + gama_b )  * (1.00 + sinpsi)  + gama_b ); 
		  Partial_Ep_gama_c  = Partial_Ep_gama_c /  aux;

		  
                  //********************************************************  
		  
		  
		  Partial_A_gama_a  = 1.00 + (1.00/3.00) * sinpsi + (1.00/3.00) * cospsi * gama_a * Partial_Dilatancy * Partial_Ep_gama_a 
		                     + (1.00/3.00) * cospsi * gama_b * Partial_Dilatancy * Partial_Ep_gama_a;
		
		  Partial_A_gama_b   = 1.00 + (1.00/3.00) * sinpsi + (1.00/3.00) * cospsi * gama_a * Partial_Dilatancy * Partial_Ep_gama_b
		                      + (1.00/3.00) * cospsi * gama_b * Partial_Dilatancy * Partial_Ep_gama_b;
		  
		  Partial_A_gama_c   = (1.00/3.00) * gama_a * cospsi * Partial_Dilatancy * Partial_Ep_gama_c
		                       + (1.00/3.00) * gama_b * cospsi * Partial_Dilatancy * Partial_Ep_gama_c + 2.00/3.00;

		  //******************************************************** 
		  
		  
		  Partial_B_gama_a   = (-2.00/3.00) * sinpsi - (2.00/3.00) * gama_a * cospsi * Partial_Dilatancy * Partial_Ep_gama_a
		                      + (1.00/3.00) * gama_b * cospsi * Partial_Dilatancy * Partial_Ep_gama_a;
		  

		  Partial_B_gama_b   = (1.00/3.00) * sinpsi - (2.00/3.00) * gama_a * cospsi * Partial_Dilatancy * Partial_Ep_gama_b
		                      + (1.00/3.00) * gama_b * cospsi * Partial_Dilatancy * Partial_Ep_gama_b - 1.00;
      
		 
		  Partial_B_gama_c   = (-2.00/3.00) * gama_a *cospsi * Partial_Dilatancy * Partial_Ep_gama_c 
		                      +(1.00/3.00) * gama_b *cospsi * Partial_Dilatancy * Partial_Ep_gama_c
		                      - 1.00/3.00;
		  
		  //********************************************************  


		  Partial_C_gama_a   = (1.00/3.00) * sinpsi -1.00 + (1.00/3.00) * gama_a * cospsi * Partial_Dilatancy * Partial_Ep_gama_a
		                       -(2.00/3.00) * gama_b * cospsi * Partial_Dilatancy * Partial_Ep_gama_a;
		  
		  Partial_C_gama_b   = (-2.00/3.00) * sinpsi + (1.00/3.00) * gama_a * cospsi * Partial_Dilatancy * Partial_Ep_gama_b
		                       -(2.00/3.00) * gama_b * cospsi * Partial_Dilatancy * Partial_Ep_gama_b;
				       
		  Partial_C_gama_c   = (1.00/3.00) * gama_a * cospsi * Partial_Dilatancy * Partial_Ep_gama_c 
		                       -(2.00/3.00) * gama_b * cospsi * Partial_Dilatancy * Partial_Ep_gama_c 
		                       - 1.00/3.00;
		  

		  //******************************************************** 
		  
		  
		  Partial_D_gama_a   = 2.00 * sinpsi + 2.00 * gama_a * cospsi * Partial_Dilatancy * Partial_Ep_gama_a
		                       + 2.00 * gama_b * cospsi * Partial_Dilatancy * Partial_Ep_gama_a;
		  
		  Partial_D_gama_b   = 2.00 * sinpsi + 2.00 * gama_a * cospsi * Partial_Dilatancy * Partial_Ep_gama_b
		                       + 2.00 * gama_b * cospsi * Partial_Dilatancy * Partial_Ep_gama_b;
		  
		  Partial_D_gama_c   = 2.00 * gama_a * cospsi * Partial_Dilatancy * Partial_Ep_gama_c 
		                       +2.00 * gama_b * cospsi * Partial_Dilatancy * Partial_Ep_gama_c + 1.00;

		   d(0,0) = 
		  ( PrincipalStress[0] + PrincipalStress[2] ) * cosphi * Partial_Friction * Partial_Ep_gama_a
		  - 2.00 * cosphi * Partial_Cohesion * Partial_Ep_gama_a
		  + 2.00 * sinphi * (mMorhCoulomb->mcurrent_cohesion) * Partial_Friction * Partial_Ep_gama_a
		  + 2.00 * G * Partial_C_gama_a 
		  - 2.00 * G * sinphi *  Partial_A_gama_a 
		  - 2.00 * G * cosphi * A * Partial_Friction * Partial_Ep_gama_a
		  - 2.00 * K * cosphi * D * Partial_Friction * Partial_Ep_gama_a
		  - 2.00 * G * sinphi * Partial_C_gama_a 
		  - 2.00 * G * cosphi * C * Partial_Friction * Partial_A_gama_a 
		  - 2.00  * K * sinphi * Partial_D_gama_a 
		  - 2.00 * G * Partial_A_gama_a;    
		  
		  d(0,1) = 
		  ( PrincipalStress[0] + PrincipalStress[2] ) * cosphi * Partial_Friction * Partial_Ep_gama_b
		  - 2.00 * cosphi * Partial_Cohesion * Partial_Ep_gama_b
		  + 2.00 * sinphi * (mMorhCoulomb->mcurrent_cohesion) * Partial_Friction * Partial_Ep_gama_b
		  + 2.00 * G * Partial_C_gama_b 
		  - 2.00 * G * sinphi *  Partial_A_gama_b 
		  - 2.00 * G * cosphi * A * Partial_Friction * Partial_Ep_gama_b
		  - 2.00 * K * cosphi * D * Partial_Friction * Partial_Ep_gama_b
		  - 2.00 * G * sinphi * Partial_C_gama_b 
		  - 2.00 * G * cosphi * C * Partial_Friction * Partial_A_gama_b 
		  - 2.00  * K * sinphi * Partial_D_gama_b 
		  - 2.00 * G * Partial_A_gama_b;  
		  
		  d(0,2) = 
		  ( PrincipalStress[0] + PrincipalStress[2] ) * cosphi * Partial_Friction * Partial_Ep_gama_c
		  - 2.00 * cosphi * Partial_Cohesion * Partial_Ep_gama_c
		  + 2.00 * sinphi * (mMorhCoulomb->mcurrent_cohesion) * Partial_Friction * Partial_Ep_gama_c
		  + 2.00 * G * Partial_C_gama_c 
		  - 2.00 * G * sinphi *  Partial_A_gama_c 
		  - 2.00 * G * cosphi * A * Partial_Friction * Partial_Ep_gama_c
		  - 2.00 * K * cosphi * D * Partial_Friction * Partial_Ep_gama_c
		  - 2.00 * G * sinphi * Partial_C_gama_c 
		  - 2.00 * G * cosphi * C * Partial_Friction * Partial_A_gama_c 
		  - 2.00  * K * sinphi * Partial_D_gama_c 
		  - 2.00 * G * Partial_A_gama_c;  
				       
	        
		  d(1,0) = 
		  ( PrincipalStress[0] + PrincipalStress[1] ) * cosphi * Partial_Friction * Partial_Ep_gama_a
		  - 2.00 * cosphi * Partial_Cohesion * Partial_Ep_gama_a
		  + 2.00 * sinphi * (mMorhCoulomb->mcurrent_cohesion) * Partial_Friction * Partial_Ep_gama_a
		  - 2.00 * G * Partial_A_gama_a 
		  + 2.00 * G * Partial_B_gama_a  
		  - 2.00 * G * cosphi * A * Partial_Friction * Partial_Ep_gama_a 
		  - 2.00 * G * sinphi *  Partial_A_gama_a 
		  - 2.00 * G * cosphi * B * Partial_Friction * Partial_Ep_gama_a 
		  - 2.00 * G * sinphi * Partial_B_gama_a  
		  - 2.00 * K * cosphi * D * Partial_Friction * Partial_Ep_gama_a
		  - 2.00 * K * sinphi * Partial_D_gama_a;
		 
		  
		  d(1,1) = 
		  ( PrincipalStress[0] + PrincipalStress[1] ) * cosphi * Partial_Friction * Partial_Ep_gama_b
		  - 2.00 * cosphi * Partial_Cohesion * Partial_Ep_gama_b
		  + 2.00 * sinphi * (mMorhCoulomb->mcurrent_cohesion) * Partial_Friction * Partial_Ep_gama_b
		  - 2.00 * G * Partial_A_gama_b 
		  + 2.00 * G * Partial_B_gama_b  
		  - 2.00 * G * cosphi * A * Partial_Friction * Partial_Ep_gama_b 
		  - 2.00 * G * sinphi *  Partial_A_gama_b 
		  - 2.00 * G * cosphi * B * Partial_Friction * Partial_Ep_gama_b 
		  - 2.00 * G * sinphi * Partial_B_gama_b  
		  - 2.00 * K * cosphi * D * Partial_Friction * Partial_Ep_gama_b
		  - 2.00 * K * sinphi * Partial_D_gama_b;
		  
		  d(1,2) = 
		  ( PrincipalStress[0] + PrincipalStress[1] ) * cosphi * Partial_Friction * Partial_Ep_gama_c
		  - 2.00 * cosphi * Partial_Cohesion * Partial_Ep_gama_c
		  + 2.00 * sinphi * (mMorhCoulomb->mcurrent_cohesion) * Partial_Friction * Partial_Ep_gama_c
		  - 2.00 * G * Partial_A_gama_c 
		  + 2.00 * G * Partial_B_gama_c  
		  - 2.00 * G * cosphi * A * Partial_Friction * Partial_Ep_gama_c 
		  - 2.00 * G * sinphi *  Partial_A_gama_c 
		  - 2.00 * G * cosphi * B * Partial_Friction * Partial_Ep_gama_c 
		  - 2.00 * G * sinphi * Partial_B_gama_c  
		  - 2.00 * K * cosphi * D * Partial_Friction * Partial_Ep_gama_c
		  - 2.00 * K * sinphi * Partial_D_gama_c;

		  
		  d(2,0) = -2.00 * G * Partial_A_gama_a - K * Partial_D_gama_a;
		  d(2,1) = -2.00 * G * Partial_A_gama_b - K * Partial_D_gama_b;
		  d(2,2) = -2.00 * G * Partial_A_gama_c - K * Partial_D_gama_c;
		  

		  singular       =  SD_MathUtils<double>::InvertMatrix(d, d_inv);
		  ddgama         =  -Vector(prod(d_inv, residual));
  
		  //Compute Newton-Raphson increment and update variables DGAMA and DGAMB
		  noalias(dgama) += ddgama;  
		  
		 //Updating Materiales
		 aux_1   = (sinpsi+1.00) * (gama_a + gama_b) + gama_c  ; aux_1 = aux_1 * aux_1;  
		 aux_2   = (sinpsi-1.00) *(sinpsi-1.00) * ( gama_a * gama_a + gama_b*gama_b) ;  
		 aux     = raiz2d3 * std::sqrt(aux_1 + aux_2);
		 
		 m_modified_morh_coulomb_maccumulated_plastic_strain_current = m_modified_morh_coulomb_maccumulated_plastic_strain_old + aux;

		 ///WARNING UPDATE INTERNAL VARIABLES
                 /*
                   cohesion;
		   dilatancy;
		   friction;
		   ft;                           
                 */
		 
		 
		  sinphi             =   std::sin(PI * mMorhCoulomb->minternal_friction_angle  / 180.00);
		  cosphi             =   std::cos(PI * mMorhCoulomb->minternal_friction_angle  / 180.00);
		  sinpsi             =   std::sin(PI * mMorhCoulomb->mcurrent_dilatancy_angle  / 180.00);
		  cospsi             =   std::cos(PI * mMorhCoulomb->mcurrent_dilatancy_angle  / 180.00); 
		    
		  A =  (gama_a + gama_b) * (1.00 + (1.00/3.00) * sinpsi) + (2.00/3.00) * gama_c; 
		  B = -(2.00/3.00) * gama_a * sinpsi + gama_b * ((1.00/3.00)*sinpsi - 1.00) - (1.00/3.00) * gama_c;                 
		  C = gama_a * ( (1.00/3.00) *sinpsi  - 1.00 ) -(2.00/3.00) * gama_b * sinpsi   - (1.00/3.00) * gama_c;
		  D = 2.00  * sinpsi * (gama_a + gama_b) + gama_c;
		  
		  
		  residual[0]   = (PrincipalStress[0] - PrincipalStress[2]) + (PrincipalStress[0] + PrincipalStress[2]) * sinphi - 2.00 *  cosphi * mMorhCoulomb->mcurrent_cohesion;
		  residual[1]   = (PrincipalStress[0] - PrincipalStress[1]) + (PrincipalStress[0] + PrincipalStress[1]) * sinphi - 2.00 *  cosphi * mMorhCoulomb->mcurrent_cohesion;
		  residual[2]   =  PrincipalStress[0] - mRankine->mcurrent_Ft;
		  
	          residual[0]   = residual[0]-2.00 * G * A + 2.00 * G * C - 2.00 * G * sinphi * A - 2.00 * G * sinphi * C - 2.00 * K * sinphi * D;
		  residual[1]   = residual[1]-2.00 * G * A + 2.00 * G * B - 2.00 * G * sinphi * A - 2.00 * G * sinphi * B - 2.00 * K * sinphi * D;
		  residual[2]   = residual[2]-2.00 * G * A - K * D; 
		
		  norma         = norm_2(residual); 
		
		if(iter>=max){  
		  KRATOS_WATCH(d)
		  KRATOS_WATCH(residual)
		  KRATOS_WATCH(dgama) 
		  KRATOS_WATCH("RETURN MAPPING OF MAIN PLANE, CORNER AND SIGMA 1 NOT CONVERGED")
		  KRATOS_ERROR(std::logic_error,  "RETURN MAPPING OF MAIN PLANE AND SIGMA 1 NOT CONVERGED" , ""); 
		 }
		}

		Sigma[0] = PrincipalStress[0] -  2.00 * G * A - K * D; 
                Sigma[1] = PrincipalStress[1] -  2.00 * G * B - K * D; 
                Sigma[2] = PrincipalStress[2] -  2.00 * G * C - K * D; 
		
		bool cond_a = mRankine->CheckValidity(Sigma);
		bool cond_b = mRankine->CheckPlasticAdmisibility(Sigma); // false when no plastic
		bool cond_3 = false;
// 		if( dgama[0]< 0.00 || dgama[1]<0.00)
// 		     cond_3 = true;
		
		
		if(cond_a==true && cond_b==false && cond_3 ==false)
		  {
		    //updating the correct principal pastic strain
		    mPrincipalPlasticStrain_current[0] = mPrincipalPlasticStrain_old[0] +  (gama_a + gama_b ) * (1.00 + sinpsi) + gama_c;
		    mPrincipalPlasticStrain_current[1] = mPrincipalPlasticStrain_old[1] +  (sinpsi -1.00) * gama_b ;
		    mPrincipalPlasticStrain_current[2] = mPrincipalPlasticStrain_old[2] +  (sinpsi -1.00) * gama_a;
		    return true;
		  }
		else
		  return false;
    }
    
    
    bool Modified_Morh_Coulomb_Yield_Function::Return_Mapping_Intersection_Main_Plane_Corner_And_Sigma1_And_Sigma_2_Tesile_Plane(const array_1d<double,3>& PrincipalStress, const array_1d<unsigned int,3>& order, array_1d<double,3>& Sigma)
    {
      
	const double& E              =   (*mprops)[YOUNG_MODULUS];
	const double& NU             =   (*mprops)[POISSON_RATIO];
	const double G               =   0.5 * E / (1.00 + NU);
	const double K               =   E / (3.00 * (1.00-2.00 * NU) );
	const double toler           =   1E-9;
	const unsigned max           =   100;

	array_1d<double, 4> dgama    =   ZeroVector(4);     
	array_1d<double, 4> ddgama   =   ZeroVector(4); 
	array_1d<double, 4> residual =   ZeroVector(4);
	Matrix d                     =   ZeroMatrix(4,4);
	Matrix d_inv                 =   ZeroMatrix(4,4);

	unsigned iter                =   0;
	double sinphi                =   std::sin(PI * mMorhCoulomb->minternal_friction_angle  / 180.00);
	double cosphi                =   std::cos(PI * mMorhCoulomb->minternal_friction_angle  / 180.00);
	double sinpsi                =   std::sin(PI * mMorhCoulomb->mcurrent_dilatancy_angle  / 180.00);
	double cospsi                =   std::cos(PI * mMorhCoulomb->mcurrent_dilatancy_angle  / 180.00); 

	
	residual[0]   = (PrincipalStress[0] - PrincipalStress[2]) + (PrincipalStress[0] + PrincipalStress[2]) * sinphi - 2.00 *  cosphi * mMorhCoulomb->mcurrent_cohesion;
	residual[1]   = (PrincipalStress[1] - PrincipalStress[2]) + (PrincipalStress[1] + PrincipalStress[2]) * sinphi - 2.00 *  cosphi * mMorhCoulomb->mcurrent_cohesion;
	residual[2]   =  PrincipalStress[0] - mRankine->mcurrent_Ft;
	residual[3]   =  PrincipalStress[1] - mRankine->mcurrent_Ft;
       

	
	double Partial_Cohesion  = 0.00;   
	double Partial_Friction  = 0.00;
	double Partial_Dilatancy = 0.00;

	double A = 0.00; 
	double B = 0.00; 
	double C = 0.00; 
	double D = 0.00; 


	double Partial_A_gama_a  = 0.00; 
	double Partial_A_gama_b  = 0.00;
	double Partial_A_gama_c  = 0.00;
        double Partial_A_gama_d  = 0.00; 
	
	double Partial_B_gama_a  = 0.00; 
	double Partial_B_gama_b  = 0.00;
	double Partial_B_gama_c  = 0.00;
        double Partial_B_gama_d  = 0.00;
	
	double Partial_C_gama_a  = 0.00; 
	double Partial_C_gama_b  = 0.00;
	double Partial_C_gama_c  = 0.00;
        double Partial_C_gama_d  = 0.00;
	
	double Partial_D_gama_a  = 0.00; 
	double Partial_D_gama_b  = 0.00;
	double Partial_D_gama_c  = 0.00;
	double Partial_D_gama_d  = 0.00;

	double Partial_Ep_gama_a = 0.00; 
	double Partial_Ep_gama_b = 0.00;
	double Partial_Ep_gama_c = 0.00;
	double Partial_Ep_gama_d = 0.00;
	

	double norma         = norm_2(residual);
	int singular         = 0;  
	const double raiz2d3 = 0.8164958092773; 
	double aux_1         = 0.00; 
	double aux_2         = 0.00;   
	double aux_3         = 0.00;  
	double aux           = 1.00;
	
	double& gama_a = dgama[0];
	double& gama_b = dgama[1];
	double& gama_c = dgama[2];
	double& gama_d = dgama[3];      
	
	while(fabs(norma)>toler && iter++ < max )
		{
		  
		  Partial_Cohesion   = 0.00;   //  WARNING->deberia cambiar valor
		  Partial_Friction   = 0.00;   //  WARNING->deberia cambiar valor
		  Partial_Dilatancy  = 0.00;   //  WARNING->deberia cambiar valor

		  Partial_Ep_gama_a  = (2.00/3.00) *( ( (gama_a)*(1.00 + sinpsi) + gama_c ) * (1.00 + sinpsi) +  (gama_a + gama_b)  * (sinpsi-1.00) * (sinpsi - 1.00 ) ) ; 
		  Partial_Ep_gama_a  = Partial_Ep_gama_a / aux;

		  Partial_Ep_gama_b  = (2.00/3.00) *( ( (gama_b)*(1.00 + sinpsi) + gama_d ) * (1.00 + sinpsi) +  ( gama_a + gama_b) * (sinpsi-1.00) * (sinpsi-1.00 ) ) ; 
		  Partial_Ep_gama_b  = Partial_Ep_gama_b /  aux;
		  
		  Partial_Ep_gama_c  = (2.00/3.00) * ( ( gama_a )  * (1.00 + sinpsi)  + gama_c ); 
		  Partial_Ep_gama_c  = Partial_Ep_gama_c /  aux;
		  
		  Partial_Ep_gama_d  = (2.00/3.00) * ( (gama_b )  * (1.00 + sinpsi)  + gama_d ); 
		  Partial_Ep_gama_d  = Partial_Ep_gama_d /  aux;

		  
                  //********************************************************  
		  
		  
		  Partial_A_gama_a  = 1.00 + (1.00/3.00) * sinpsi + (1.00/3.00) * cospsi * gama_a * Partial_Dilatancy * Partial_Ep_gama_a 
		                     - (2.00/3.00) * cospsi * gama_b * Partial_Dilatancy * Partial_Ep_gama_a;
		
		  Partial_A_gama_b   = - (2.00/3.00) * sinpsi + (1.00/3.00) * cospsi * gama_a * Partial_Dilatancy * Partial_Ep_gama_b
		                      - (2.00/3.00) * cospsi * gama_b * Partial_Dilatancy * Partial_Ep_gama_b;
		  
		  Partial_A_gama_c   = (1.00/3.00) * gama_a * cospsi * Partial_Dilatancy * Partial_Ep_gama_c
		                       - (2.00/3.00) * gama_b * cospsi * Partial_Dilatancy * Partial_Ep_gama_c + 2.00/3.00;
		
		  Partial_A_gama_d   = (1.00/3.00) * gama_a * cospsi * Partial_Dilatancy * Partial_Ep_gama_d
		                       - (2.00/3.00) * gama_b * cospsi * Partial_Dilatancy * Partial_Ep_gama_d - 1.00/3.00;

		  //******************************************************** 
		  
		  
		  Partial_B_gama_a   = (-2.00/3.00) * sinpsi - (2.00/3.00) * gama_a * cospsi * Partial_Dilatancy * Partial_Ep_gama_a
		                      + (1.00/3.00) * gama_b * cospsi * Partial_Dilatancy * Partial_Ep_gama_a;
		  

		  Partial_B_gama_b   = 1.00 + (1.00/3.00) * sinpsi - (2.00/3.00) * gama_a * cospsi * Partial_Dilatancy * Partial_Ep_gama_b
		                      + (1.00/3.00) * gama_b * cospsi * Partial_Dilatancy * Partial_Ep_gama_b;
      
		 
		  Partial_B_gama_c   = (-2.00/3.00) * gama_a *cospsi * Partial_Dilatancy * Partial_Ep_gama_c 
		                      +(1.00/3.00) * gama_b *cospsi * Partial_Dilatancy * Partial_Ep_gama_c
		                      - 1.00/3.00;
				      
		  Partial_B_gama_d   = (-2.00/3.00) * gama_a *cospsi * Partial_Dilatancy * Partial_Ep_gama_d 
		                      +(1.00/3.00) * gama_b *cospsi * Partial_Dilatancy * Partial_Ep_gama_d
		                      + 2.00/3.00;
		  
		  //********************************************************  


		  Partial_C_gama_a   = (1.00/3.00) * sinpsi + (1.00/3.00) * gama_a * cospsi * Partial_Dilatancy * Partial_Ep_gama_a
		                        + (1.00/3.00) * gama_b * cospsi * Partial_Dilatancy * Partial_Ep_gama_a
		                        -1.00;
		  
		  Partial_C_gama_b   = (1.00/3.00) * sinpsi + (1.00/3.00) * gama_a * cospsi * Partial_Dilatancy * Partial_Ep_gama_b
		                       +(1.00/3.00) * gama_b * cospsi * Partial_Dilatancy * Partial_Ep_gama_b 
		                       -1.00;
				       
		  Partial_C_gama_c   = (1.00/3.00) * gama_a * cospsi * Partial_Dilatancy * Partial_Ep_gama_c 
		                       +(1.00/3.00) * gama_b * cospsi * Partial_Dilatancy * Partial_Ep_gama_c 
		                       - 1.00/3.00;
		  
	          Partial_C_gama_d   = (1.00/3.00) * gama_a * cospsi * Partial_Dilatancy * Partial_Ep_gama_d 
		                       +(1.00/3.00) * gama_b * cospsi * Partial_Dilatancy * Partial_Ep_gama_d 
		                       - 1.00/3.00;

		  //******************************************************** 
		  
		  
		  Partial_D_gama_a   = 2.00 * sinpsi + 2.00 * gama_a * cospsi * Partial_Dilatancy * Partial_Ep_gama_a
		                       + 2.00 * gama_b * cospsi * Partial_Dilatancy * Partial_Ep_gama_a;
		  
		  Partial_D_gama_b   = 2.00 * sinpsi + 2.00 * gama_a * cospsi * Partial_Dilatancy * Partial_Ep_gama_b
		                       + 2.00 * gama_b * cospsi * Partial_Dilatancy * Partial_Ep_gama_b;
		  
		  Partial_D_gama_c   = 2.00 * gama_a * cospsi * Partial_Dilatancy * Partial_Ep_gama_c 
		                       +2.00 * gama_b * cospsi * Partial_Dilatancy * Partial_Ep_gama_c + 1.00;
				       
		  Partial_D_gama_d   = 2.00 * gama_a * cospsi * Partial_Dilatancy * Partial_Ep_gama_d 
		                       +2.00 * gama_b * cospsi * Partial_Dilatancy * Partial_Ep_gama_d + 1.00;
	
		
		  
				       
		d(0,0) =   
		  ( PrincipalStress[0] + PrincipalStress[2] ) * cosphi * Partial_Friction * Partial_Ep_gama_a
		  - 2.00 * cosphi * Partial_Cohesion * Partial_Ep_gama_a
		  + 2.00 * sinphi * (mMorhCoulomb->mcurrent_cohesion) * Partial_Friction * Partial_Ep_gama_a
		  + 2.00 * G * Partial_C_gama_a 
		  - 2.00 * G * sinphi *  Partial_A_gama_a 
		  - 2.00 * G * cosphi * A * Partial_Friction * Partial_Ep_gama_a
		  - 2.00 * K * cosphi * D * Partial_Friction * Partial_Ep_gama_a
		  - 2.00 * G * sinphi * Partial_C_gama_a 
		  - 2.00 * G * cosphi * C * Partial_Friction * Partial_A_gama_a 
		  - 2.00  * K * sinphi * Partial_D_gama_a 
		  - 2.00 * G * Partial_A_gama_a; 
		  
				       
		  d(0,1) = 
		  ( PrincipalStress[0] + PrincipalStress[2] ) * cosphi * Partial_Friction * Partial_Ep_gama_b
		  - 2.00 * cosphi * Partial_Cohesion * Partial_Ep_gama_b
		  + 2.00 * sinphi * (mMorhCoulomb->mcurrent_cohesion) * Partial_Friction * Partial_Ep_gama_b
		  + 2.00 * G * Partial_C_gama_b 
		  - 2.00 * G * sinphi *  Partial_A_gama_b 
		  - 2.00 * G * cosphi * A * Partial_Friction * Partial_Ep_gama_b
		  - 2.00 * K * cosphi * D * Partial_Friction * Partial_Ep_gama_b
		  - 2.00 * G * sinphi * Partial_C_gama_b 
		  - 2.00 * G * cosphi * C * Partial_Friction * Partial_A_gama_b 
		  - 2.00  * K * sinphi * Partial_D_gama_b
		  - 2.00 * G * Partial_A_gama_b;    
		  
		  
		  d(0,2) = 
		  ( PrincipalStress[0] + PrincipalStress[2] ) * cosphi * Partial_Friction * Partial_Ep_gama_c
		  - 2.00 * cosphi * Partial_Cohesion * Partial_Ep_gama_c
		  + 2.00 * sinphi * (mMorhCoulomb->mcurrent_cohesion) * Partial_Friction * Partial_Ep_gama_c
		  + 2.00 * G * Partial_C_gama_c 
		  - 2.00 * G * sinphi *  Partial_A_gama_c 
		  - 2.00 * G * cosphi * A * Partial_Friction * Partial_Ep_gama_c
		  - 2.00 * K * cosphi * D * Partial_Friction * Partial_Ep_gama_c
		  - 2.00 * G * sinphi * Partial_C_gama_c 
		  - 2.00 * G * cosphi * C * Partial_Friction * Partial_A_gama_c 
		  - 2.00  * K * sinphi * Partial_D_gama_c
		  - 2.00 * G * Partial_A_gama_c;    
		  
		  
		  d(0,3) = 
		  ( PrincipalStress[0] + PrincipalStress[2] ) * cosphi * Partial_Friction * Partial_Ep_gama_d
		  - 2.00 * cosphi * Partial_Cohesion * Partial_Ep_gama_d
		  + 2.00 * sinphi * (mMorhCoulomb->mcurrent_cohesion) * Partial_Friction * Partial_Ep_gama_d
		  + 2.00 * G * Partial_C_gama_d 
		  - 2.00 * G * sinphi *  Partial_A_gama_d 
		  - 2.00 * G * cosphi * A * Partial_Friction * Partial_Ep_gama_d
		  - 2.00 * K * cosphi * D * Partial_Friction * Partial_Ep_gama_d
		  - 2.00 * G * sinphi * Partial_C_gama_d 
		  - 2.00 * G * cosphi * C * Partial_Friction * Partial_A_gama_d 
		  - 2.00  * K * sinphi * Partial_D_gama_d
		  - 2.00 * G * Partial_A_gama_d;    
		  
		  
                   //************************************************************************************************
		  
				       
		       
		d(1,0) = 
		  ( PrincipalStress[1] + PrincipalStress[2] ) * cosphi * Partial_Friction * Partial_Ep_gama_a
		  - 2.00 * cosphi * Partial_Cohesion * Partial_Ep_gama_a
		  + 2.00 * sinphi * (mMorhCoulomb->mcurrent_cohesion) * Partial_Friction * Partial_Ep_gama_a
		  - 2.00 * G * Partial_B_gama_a 
		  + 2.00 * G * Partial_C_gama_a  
		  - 2.00 * G * cosphi * B * Partial_Friction * Partial_Ep_gama_a 
		  - 2.00 * G * sinphi * Partial_B_gama_a 
		  - 2.00 * G * cosphi * C * Partial_Friction * Partial_Ep_gama_a 
		  - 2.00 * G * sinphi * Partial_C_gama_a  
		  - 2.00 * K * cosphi * D * Partial_Friction * Partial_Ep_gama_a
		  - 2.00 * K * sinphi * Partial_D_gama_a;
		 

		  d(1,1) = 
		  ( PrincipalStress[1] + PrincipalStress[2] ) * cosphi * Partial_Friction * Partial_Ep_gama_b
		  - 2.00 * cosphi * Partial_Cohesion * Partial_Ep_gama_b
		  + 2.00 * sinphi * (mMorhCoulomb->mcurrent_cohesion) * Partial_Friction * Partial_Ep_gama_b
		  - 2.00 * G * Partial_B_gama_b 
		  + 2.00 * G * Partial_C_gama_b  
		  - 2.00 * G * cosphi * B * Partial_Friction * Partial_Ep_gama_b 
		  - 2.00 * G * sinphi * Partial_B_gama_b 
		  - 2.00 * G * cosphi * C * Partial_Friction * Partial_Ep_gama_b 
		  - 2.00 * G * sinphi * Partial_C_gama_b  
		  - 2.00 * K * cosphi * D * Partial_Friction * Partial_Ep_gama_b
		  - 2.00 * K * sinphi * Partial_D_gama_b;
		  
		  
		  d(1,2) = 
		  ( PrincipalStress[1] + PrincipalStress[2] ) * cosphi * Partial_Friction * Partial_Ep_gama_c
		  - 2.00 * cosphi * Partial_Cohesion * Partial_Ep_gama_c
		  + 2.00 * sinphi * (mMorhCoulomb->mcurrent_cohesion) * Partial_Friction * Partial_Ep_gama_c
		  - 2.00 * G * Partial_B_gama_c 
		  + 2.00 * G * Partial_C_gama_c  
		  - 2.00 * G * cosphi * B * Partial_Friction * Partial_Ep_gama_c 
		  - 2.00 * G * sinphi * Partial_B_gama_c 
		  - 2.00 * G * cosphi * C * Partial_Friction * Partial_Ep_gama_c 
		  - 2.00 * G * sinphi * Partial_C_gama_c  
		  - 2.00 * K * cosphi * D * Partial_Friction * Partial_Ep_gama_c
		  - 2.00 * K * sinphi * Partial_D_gama_c;
		  
		  
		  d(1,3) = 
		  ( PrincipalStress[1] + PrincipalStress[2] ) * cosphi * Partial_Friction * Partial_Ep_gama_d
		  - 2.00 * cosphi * Partial_Cohesion * Partial_Ep_gama_d
		  + 2.00 * sinphi * (mMorhCoulomb->mcurrent_cohesion) * Partial_Friction * Partial_Ep_gama_d
		  - 2.00 * G * Partial_B_gama_d 
		  + 2.00 * G * Partial_C_gama_d  
		  - 2.00 * G * cosphi * B * Partial_Friction * Partial_Ep_gama_d 
		  - 2.00 * G * sinphi * Partial_B_gama_d 
		  - 2.00 * G * cosphi * C * Partial_Friction * Partial_Ep_gama_d 
		  - 2.00 * G * sinphi * Partial_C_gama_d  
		  - 2.00 * K * cosphi * D * Partial_Friction * Partial_Ep_gama_d
		  - 2.00 * K * sinphi * Partial_D_gama_d;
		  
		  
		  d(2,0) = -2.00 * G * Partial_A_gama_a - K * Partial_D_gama_a;
		  d(2,1) = -2.00 * G * Partial_A_gama_b - K * Partial_D_gama_b;
		  d(2,2) = -2.00 * G * Partial_A_gama_c - K * Partial_D_gama_c;
		  d(2,3) = -2.00 * G * Partial_A_gama_d - K * Partial_D_gama_d;
		  
		  d(3,0) = -2.00 * G * Partial_B_gama_a - K * Partial_D_gama_a;
		  d(3,1) = -2.00 * G * Partial_B_gama_b - K * Partial_D_gama_b;
		  d(3,2) = -2.00 * G * Partial_B_gama_c - K * Partial_D_gama_c;
		  d(3,3) = -2.00 * G * Partial_B_gama_d - K * Partial_D_gama_d;
		  
		  
		  // Trunca los termnos de esta matriz. Lo hice con el proposito de que esta matriz tubieta 
		  // inversa. Al parecer no tiene inversa cuando utilizo la forma d como tal.
		  for( unsigned int i = 0; i< 4; i++){
		      for( unsigned int j = 0; j< 4; j++)
		      {
			aux    = d(i,j);
			d(i,j) = round(aux, 3);
		      }
		  }
		    
	          aux            =   0.00;   
		  singular       =   SD_MathUtils<double>::InvertMatrix(d, d_inv);
		  ddgama         =  -Vector(prod(d_inv, residual));
  
		  //Compute Newton-Raphson increment and update variables DGAMA and DGAMB
		  noalias(dgama) += ddgama;  
		  
		 //Updating Materiales
		 aux_1   = (sinpsi+1.00) *gama_a + gama_c  ;  aux_1 = aux_1 * aux_1;  
		 aux_2   = (sinpsi+1.00) *gama_b + gama_d  ;  aux_2 = aux_2 * aux_2;
		 aux_3   = (sinpsi-1.00) *( gama_a + gama_b); aux_3 = aux_3 * aux_3;
		 aux     = raiz2d3 * std::sqrt(aux_1 + aux_2 + aux_3);
		 
		 m_modified_morh_coulomb_maccumulated_plastic_strain_current = m_modified_morh_coulomb_maccumulated_plastic_strain_old + aux;

		 ///WARNING UPDATE INTERNAL VARIABLES
                 /*
                   cohesion;
		   dilatancy;
		   friction;
		   ft;                           
                 */
		  
		  sinphi             =   std::sin(PI * mMorhCoulomb->minternal_friction_angle  / 180.00);
		  cosphi             =   std::cos(PI * mMorhCoulomb->minternal_friction_angle  / 180.00);
		  sinpsi             =   std::sin(PI * mMorhCoulomb->mcurrent_dilatancy_angle  / 180.00);
		  cospsi             =   std::cos(PI * mMorhCoulomb->mcurrent_dilatancy_angle  / 180.00); 
		 
		 
		  D =   2.00  * sinpsi  * (gama_a + gama_b) + gama_c + gama_d;
		  A =  (1.00 + sinpsi) *  gama_a + gama_c  - D/3.00; 
		  B =  (1.00 + sinpsi) *  gama_b + gama_d  - D/3.00;               
		  C =  (sinpsi -1.00)  * (gama_a + gama_b) - D/3.00;  
		  
		  
		  
		  
		  residual[0]   = (PrincipalStress[0] - PrincipalStress[2]) + (PrincipalStress[0] + PrincipalStress[2]) * sinphi - 2.00 *  cosphi * mMorhCoulomb->mcurrent_cohesion;
		  residual[1]   = (PrincipalStress[1] - PrincipalStress[2]) + (PrincipalStress[1] + PrincipalStress[2]) * sinphi - 2.00 *  cosphi * mMorhCoulomb->mcurrent_cohesion;
		  residual[2]   =  PrincipalStress[0] - mRankine->mcurrent_Ft;
		  residual[3]   =  PrincipalStress[1] - mRankine->mcurrent_Ft;
		  
	          residual[0]   = residual[0]-2.00 * G * A + 2.00 * G * C - 2.00 * G * sinphi * A - 2.00 * G * sinphi * C - 2.00 * K * sinphi * D;
		  residual[1]   = residual[1]-2.00 * G * B + 2.00 * G * C - 2.00 * G * sinphi * B - 2.00 * G * sinphi * C - 2.00 * K * sinphi * D;
		  residual[2]   = residual[2]-2.00 * G * A - K * D; 
		  residual[3]   = residual[3]-2.00 * G * B - K * D; 
		
		  norma = norm_2(residual);
	          
		if(iter>=max){
		  KRATOS_WATCH(dgama)
		  KRATOS_WATCH(norma)
		  KRATOS_WATCH("RETURN MAPPING OF MAIN PLANE AND SIGMA 1 AND SIGMA 2 NOT CONVERGED") 
		  KRATOS_ERROR(std::logic_error,  "RETURN MAPPING OF MAIN PLANE AND SIGMA 1 AND SIGMA 2 NOT CONVERGED" , ""); 
		 }
		}

		Sigma[0] = PrincipalStress[0] -  2.00 * G * A - K * D; 
                Sigma[1] = PrincipalStress[1] -  2.00 * G * B - K * D; 
                Sigma[2] = PrincipalStress[2] -  2.00 * G * C - K * D; 
		
		bool cond_a = mRankine->CheckValidity(Sigma);
		bool cond_b = mRankine->CheckPlasticAdmisibility(Sigma); // false when no plastic
		bool cond_3 = false;
// 		if( dgama[0]< 0.00 || dgama[1]<0.00)
// 		     cond_3 = true;
		
		
		if(cond_a==true && cond_b==false && cond_3 ==false)
		   {
		    //updating the correct principal pastic strain
		    mPrincipalPlasticStrain_current[0] = mPrincipalPlasticStrain_old[0] + gama_a  * (1.00 + sinpsi) + gama_c;
		    mPrincipalPlasticStrain_current[1] = mPrincipalPlasticStrain_old[1] + gama_b  * (1.00 + sinpsi) + gama_d;
		    mPrincipalPlasticStrain_current[2] = mPrincipalPlasticStrain_old[2] + (sinpsi -1.00) * ( gama_a + gama_d);
		    return true;
		   }
		else{
		  return false;
		  KRATOS_ERROR(std::logic_error,  "RETURN MAPPING OF MAIN PLANE AND SIGMA 1 AND SIGMA 2 IS FALSE" , ""); 
		}
	
      
    }
    

    void Modified_Morh_Coulomb_Yield_Function::GetValue(const Variable<double>& rVariable, double& Result)
      {
	
	if(rVariable==COHESION)
	  Result = mMorhCoulomb->mcurrent_cohesion;
        
	if(rVariable == DILATANCY_ANGLE)
          Result = mMorhCoulomb->mcurrent_dilatancy_angle; 
	
	if(rVariable == INTERNAL_FRICTION_ANGLE)
	  Result = mMorhCoulomb->mcurrent_minternal_friction_angle;
	
        if(rVariable == DAMAGE)
	   Result = mpastic_damage_current;
	
      }
      
     void Modified_Morh_Coulomb_Yield_Function::GetValue(const Variable<Matrix>& rVariable, Matrix& Result)
	     {
	       unsigned int size = 0;  
	       if(this->mState==Plane_Strain)
		  size = 4;
	       else
		  size = 6;
	         
		if (rVariable==GREEN_LAGRANGE_PLASTIC_STRAIN_TENSOR)
		{
		Result.resize(1, size);
		for(unsigned int i = 0; i< size; i++ )
		     Result(0,i) = mplastic_strain(i); 
		
                }  
            return; 
            }
            
        void Modified_Morh_Coulomb_Yield_Function::CalculatePlasticDamage(const array_1d<double,3>& Sigma)
        {
	  array_1d<double, 3> DeltaPlasticStrain = mPrincipalPlasticStrain_current - mPrincipalPlasticStrain_old;
	  double disipation = inner_prod(Sigma, DeltaPlasticStrain);  
	  double teta_a     =  Tensor_Utils<double>::Mc_aully(Sigma);
          double teta_b     =  norm_1(Sigma);
	  double teta       = 0.00;
          if (fabs(teta_b) < 1E-10)
              {teta = 0.00;}
          else
              {teta = teta_a/teta_b;}

	  // updating plastic_damage
	  // computing Kp_punto
	  double gc_p     = (*mprops)[CRUSHING_ENERGY]/mlength; 
	  double gf_p     = (*mprops)[FRACTURE_ENERGY]/mlength;

	  double kp_punto           = ( teta/gf_p + (1.00-teta)/(gc_p) ) * disipation;
	  mpastic_damage_current    =  mpastic_damage_old + kp_punto;                 
	  if(mpastic_damage_current>=1.00) 
	      mpastic_damage_current = 1.00; 
        }
      

  }
       




