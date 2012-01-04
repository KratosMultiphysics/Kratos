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


///  WARNING = mirar los casos de switch
///  WARNING = mirar las funciones de los modelos
///  WARNING = mirar matrix


#include "soft_hard_behavior/softening_hardening_criteria.h"
#include "fluency_criteria/modified_morh_coulomb_yield_function.h"

namespace Kratos
  {

            typedef Morh_Coulomb_Yield_Function MorhCoulomb;
	    typedef Isotropic_Rankine_Yield_Function Rankine;
     
            Modified_Morh_Coulomb_Yield_Function::Modified_Morh_Coulomb_Yield_Function()
            {
	    }
           
            ///WARNING = THE STRUCTURE COULD BE DANGEROUS 
            Modified_Morh_Coulomb_Yield_Function::Modified_Morh_Coulomb_Yield_Function(
            myState State, 
            const  MorhCoulombPointerType MorhCoulomb, RankinePointerType Rankine)
	    {
	       mState       = State; 
               mMorhCoulomb = MorhCoulomb;
               mRankine     = Rankine;
	    }

             /*
	     Modified_Morh_Coulomb_Yield_Function::Modified_Morh_Coulomb_Yield_Function(
	     const SoftHardPointerType& SofteningBehavior,           /// For Ft
	     const SoftHardPointerType& SofteningBehaviorCohesion,   /// For Cohesion
	     const SoftHardPointerType& SofteningBehaviorFriction,   /// For Friction Angle
             const SoftHardPointerType& SofteningBehaviorDilatancy,  /// For Dilatancy
	     const myState& State, 
             const myPotencialPlastic& PotencialPlastic):
             FluencyCriteria(),
             Isotropic_Rankine_Yield_Function(SofteningBehavior,State), 
             Morh_Coulomb_Yield_Function(SofteningBehaviorCohesion, SofteningBehaviorFriction, SofteningBehaviorDilatancy, State, PotencialPlastic)
             {
	       mState = State; 
	     }
	     */

             Modified_Morh_Coulomb_Yield_Function::~Modified_Morh_Coulomb_Yield_Function() {}

             void Modified_Morh_Coulomb_Yield_Function::InitializeMaterial(const Properties& props)
             {
	       mprops =  &props;
	       mMorhCoulomb->InitializeMaterial(props);
	       mRankine->InitializeMaterial(props);  
	      
	       m_modified_morh_coulomb_maccumulated_plastic_strain_old         = 0.00;
               m_modified_morh_coulomb_maccumulated_plastic_strain_current     = 0.00;
               //m_modified_morh_coulomb_maccumulated_plastic_strain_current_pos = 0.00;
               //m_modified_morh_coulomb_maccumulated_plastic_strain_current_neg = 0.00;
	       mpastic_damage_old                                              = 0.00;
               mpastic_damage_current                                          = 0.00;
	       
	       int  size = 4;
	       if(mState==Tri_D)
	          size = 6;

	       mplastic_strain.resize(size, false); 
	       mplastic_strain_old.resize(size, false);

	       mplastic_strain                                   = ZeroVector(size);
	       mplastic_strain_old                               = ZeroVector(size); 
               mElastic_strain                                   = ZeroVector(size);
               mElastic_strain_old                               = ZeroVector(size);
	       mPrincipalPlasticStrain_current                   = ZeroVector(3);
	       mPrincipalPlasticStrain_old                       = ZeroVector(3);
	     }
             
              void Modified_Morh_Coulomb_Yield_Function::GetValue(double& Result)
              {
		mlength = Result;
		mRankine->GetValue(Result);
		mMorhCoulomb->GetValue(Result);
	      }
             
             void Modified_Morh_Coulomb_Yield_Function::ReturnMapping(const Vector& StrainVector, const Vector& TrialStress,  Vector& StressVector)
             {
	       
		noalias(StressVector) = ZeroVector(3);
		const double& Young   = (*mprops)[YOUNG_MODULUS];
		const double& Poisson = (*mprops)[POISSON_RATIO];
		//const double Gmodu    = Young/(2.00 * (1.00 + Poisson) );
		//const double Bulk     = Young/(3.00 * (1.00-2.00*Poisson)); 
	             
		Vector Aux_Trial_Stress;
		Aux_Trial_Stress.resize(StressVector.size());
		noalias(Aux_Trial_Stress) =  StressVector;

		array_1d<double,3> Sigma = ZeroVector(3);
		array_1d<double,3> PrincipalStress = ZeroVector(3);
		array_1d<array_1d < double,3 > ,3> EigenVectors;
		array_1d<unsigned int,3> Order;
		
		
		Vector Stress(TrialStress.size());
		Stress = ZeroVector(TrialStress.size());
		
		/// computing the trial kirchooff strain
		SpectralDecomposition(TrialStress, PrincipalStress, EigenVectors);
		IdentifyMaximunAndMinumumPrincipalStres_CalculateOrder(PrincipalStress, Order);
		noalias(Sigma)  = PrincipalStress;
		noalias(Stress) = TrialStress;
	        
		const array_1d<double,3> Aux_PrincipalStress              = PrincipalStress; 
		const array_1d<array_1d < double,3 > ,3> Aux_EigenVectors = EigenVectors; 
		const array_1d<unsigned int,3>  Aux_Order                 = Order;
		
		bool plastic_1 = mMorhCoulomb->CheckPlasticAdmisibility(PrincipalStress); 
		bool plastic_2 = mRankine->CheckPlasticAdmisibility(PrincipalStress); 
		      
                // elastic step
		if(plastic_1==false && plastic_2==false)
		 {
		    /// computing the pressure
		    CalculateElasticStrain(Stress, mElastic_strain);
		    //mpressure       =  d3 * (Sigma[0] + Sigma[1] + Sigma[2]);
		    StressVector[0] = Stress[0];
		    StressVector[1] = Stress[1];
		    StressVector[2] = Stress[2];
		    return;
		 }
		 
                else{
		// standar morh coulomb
		if(plastic_1==true)
		{ 
		  mMorhCoulomb->ReturnMapping(StrainVector, TrialStress, StressVector);
		  
		  /// For plane Strain
		  //check tensile consistency
		  Stress[0] = StressVector[0];
		  Stress[1] = StressVector[1];
		  Stress[2] = StressVector[2];
		  Stress[3] = mMorhCoulomb->msigma_z;

		  SpectralDecomposition(Stress, PrincipalStress, EigenVectors);
		  IdentifyMaximunAndMinumumPrincipalStres_CalculateOrder(PrincipalStress, Order);  
		  /// sino se cumple pasamos a la forma combinada de rankine y morh
		  bool tensile_consistency =  mRankine->CheckPlasticAdmisibility(PrincipalStress);
		  if(tensile_consistency==true) // sigue siendo plastico 
		  { 		   
		    //update variables
		    mMorhCoulomb->UpdateMaterial();
	            mRankine->UpdateMaterial();
		    // Cero value is false 1 is equal to one 
		    Sigma = ZeroVector(3);    // Principal Stress updated
		    noalias(StressVector)    = Aux_Trial_Stress;
		    noalias(PrincipalStress) = Aux_PrincipalStress;
		    //KRATOS_WATCH("AAAAAAAAAAAAAA")
		    if(Return_Mapping_Intersection_Main_Plane_And_Sigma1_Tesile_Plane(PrincipalStress , Order,  Sigma)==false)
		    {
		      //return mapping to corner
		      mMorhCoulomb->UpdateMaterial();
	              mRankine->UpdateMaterial();
		      /*
		      array_1d<double, 3> T;
		      double sinpsi  = std::sin(PI * mMorhCoulomb->mcurrent_dilatancy_angle  / 180.00);
		      const double a = 3.00 + sinpsi;
		      const double b = 3.00 - sinpsi;
		      T[0] =  2.00/a - 1.00/b;
		      T[1] = -1.00/a - 1.00/b;
		      T[2] = -1.00/a - 2.00/b;
		      */
		      //bool rigth_corner = mMorhCoulomb->ReturnToEdges(PrincipalStress);
		      /*
		      bool rigth_corner = false;  
		      if(inner_prod(T,Aux_PrincipalStress)>0)
			 rigth_corner = true;
		      */
		      //if( rigth_corner==true){//Main Plane Corner and Sigma 1
		       //KRATOS_WATCH("CCCCCCCCCCCCCCCC")
		       if(Return_Mapping_Intersection_Main_Plane_Corner_And_Sigma1_Tesile_Plane(PrincipalStress, Order, Sigma)==false){
		       //}
		       //else{//Main Plane Corner and Sigma 1 and Sigma 2
			//KRATOS_WATCH("DDDDDDDDDDDDDDDDDDD")
			Return_Mapping_Intersection_Main_Plane_Corner_And_Sigma1_And_Sigma_2_Tesile_Plane(PrincipalStress, Order, Sigma); 
		      }
		    }
		    
		    CalculatePlasticDamage(Sigma);
		    AssembleUpdateStressAndStrainTensor(Sigma,  EigenVectors,  Order, StrainVector, Stress);
		    AssembleStress(Stress, StressVector);
		  }
		  
		   else
		    {
		      /// standar morh coulomb
		      mPrincipalPlasticStrain_current = mMorhCoulomb->mPrincipalPlasticStrain_current; 
		      m_modified_morh_coulomb_maccumulated_plastic_strain_current = mMorhCoulomb->mmorh_coulomb_maccumulated_plastic_strain_current;
		      noalias(mElastic_strain) = mMorhCoulomb->mElastic_strain; 
		      mpastic_damage_current   = mMorhCoulomb->mpastic_damage_current;
		      
		      /// Updating Rankine
		      Vector Imput_Parameters_R; 
		      Imput_Parameters_R.resize(4); Imput_Parameters_R = ZeroVector(4);
		      Imput_Parameters_R[0] =  mlength;
		      Imput_Parameters_R[1] =  m_modified_morh_coulomb_maccumulated_plastic_strain_current; 
		      Imput_Parameters_R[2] =  mpastic_damage_old;
                      Imput_Parameters_R[3] =  mpastic_damage_current;
		      
		      mRankine->mcurrent_Ft =  mRankine->mpSofteningBehaviorFt->FunctionBehavior(Imput_Parameters_R);
		    }    
		}
		
		else 
		{
		 
		  mRankine->ReturnMapping(StrainVector, TrialStress, StressVector); 
		  mPrincipalPlasticStrain_current = mRankine->mPrincipalPlasticStrain_current; 
		  m_modified_morh_coulomb_maccumulated_plastic_strain_current = mRankine->mrankine_accumulated_plastic_strain_current;  
		  noalias(mElastic_strain) = mRankine->mElastic_strain; 
		  mpastic_damage_current   = mRankine->mpastic_damage_current;
		  
		  ///variable Compute Internal Variables of Morh Coulomb Models
		  Vector Imput(3);
		  Imput = ZeroVector(3);
		  Vector Imput_Parameters(4);
		  Imput_Parameters = ZeroVector(4); 
		  Imput[0] = mpastic_damage_current;
		  Imput[1] = mpastic_damage_old;
		  Imput[2] = mMorhCoulomb->mcohesion; ///old cohesion
		  mMorhCoulomb->mcurrent_cohesion =  mMorhCoulomb->mpSofteningBehavior_Cohesion->EvolucionLaws(Imput, Sigma);   
		  Imput_Parameters[0] =  mlength;
		  Imput_Parameters[1] =  m_modified_morh_coulomb_maccumulated_plastic_strain_current;
		  Imput_Parameters[3] =  mpastic_damage_current;
		  mMorhCoulomb->mcurrent_minternal_friction_angle =  mMorhCoulomb->mpSofteningBehavior_Friction->FunctionBehavior(Imput_Parameters);   
		  Imput_Parameters[2]                    =  mMorhCoulomb->mcurrent_minternal_friction_angle;
		  mMorhCoulomb->mcurrent_dilatancy_angle =  mMorhCoulomb->mpSofteningBehavior_Dilatancy->FunctionBehavior(Imput_Parameters);
		  }
		  
		 CalculatePlasticStrain(mPrincipalPlasticStrain_current, EigenVectors, Order, mplastic_strain);
		}

                double i = inner_prod(Stress, Stress);
		if(i!=i){
                  KRATOS_WATCH(Stress)
                  KRATOS_WATCH("AKIIIIIIIIIIII")
                  KRATOS_ERROR(std::logic_error,  "MODIFIEDDDDDD  " , " "); 
		}
                
		CalculateElasticStrain(Stress, mElastic_strain);
	     }
	     
	     
             void Modified_Morh_Coulomb_Yield_Function::FinalizeSolutionStep()
             { 

	       mRankine->mrankine_accumulated_plastic_strain_current           = m_modified_morh_coulomb_maccumulated_plastic_strain_current;
	       mMorhCoulomb->mmorh_coulomb_maccumulated_plastic_strain_current = m_modified_morh_coulomb_maccumulated_plastic_strain_current;	       
	       m_modified_morh_coulomb_maccumulated_plastic_strain_old         = m_modified_morh_coulomb_maccumulated_plastic_strain_current;
	       noalias(mPrincipalPlasticStrain_old)                            = mPrincipalPlasticStrain_current; 
	       noalias(mplastic_strain_old)                                    = mplastic_strain; 
	       noalias(mElastic_strain_old)                                    = mElastic_strain;
	       mpastic_damage_old                                              = mpastic_damage_current;
	       
	       /// Asumiendo las mismas deformaciones plasticas para ambos modelos
	       mRankine->mPrincipalPlasticStrain_current     = mPrincipalPlasticStrain_current;
	       mRankine->mplastic_strain                     = mplastic_strain; 
	       mRankine->mElastic_strain                     = mElastic_strain;
	       
	       mMorhCoulomb->mPrincipalPlasticStrain_current = mPrincipalPlasticStrain_current;
	       mMorhCoulomb->mplastic_strain                 = mplastic_strain; 
	       mMorhCoulomb->mElastic_strain                 = mElastic_strain;

	       
	       mMorhCoulomb->mpastic_damage_current          = mpastic_damage_current;
	       mRankine->mpastic_damage_current              = mpastic_damage_current;
	       
	       
	       mMorhCoulomb->FinalizeSolutionStep();
	       mRankine->FinalizeSolutionStep();
	       
	     }
	     
             void Modified_Morh_Coulomb_Yield_Function::UpdateMaterial()
             { 
	       mMorhCoulomb->UpdateMaterial();
	       mRankine->UpdateMaterial();
	       
               m_modified_morh_coulomb_maccumulated_plastic_strain_current =  m_modified_morh_coulomb_maccumulated_plastic_strain_old;  
	       noalias(mPrincipalPlasticStrain_current)                    =  mPrincipalPlasticStrain_old;  
	       noalias(mplastic_strain)                                    =  mplastic_strain_old; 
               mpastic_damage_current                                      =  mpastic_damage_old; 
	       noalias(mElastic_strain)                                    = mElastic_strain_old;
	     } 

	     bool Modified_Morh_Coulomb_Yield_Function::Return_Mapping_Intersection_Main_Plane_And_Sigma1_Tesile_Plane(const array_1d<double,3>& PrincipalStress, const array_1d<unsigned int,3>& order, array_1d<double,3>& Sigma)
	     {
	        
	       
	        UpdateMaterial();
		noalias(Sigma)          =   ZeroVector(3);
		const double d3         =   0.3333333333333333333; 
		const double& E         =   (*mprops)[YOUNG_MODULUS];
		const double& NU        =   (*mprops)[POISSON_RATIO];
		const double G          =   0.5 * E / (1.00 + NU);
		const double K          =   E / (3.00 * (1.00-2.00 * NU) );
		const double toler      =   1E-2;
		const unsigned max      =   1000;
		unsigned int  iter      =   0;

		
		array_1d<double, 2> dgama    = ZeroVector(2);     
		array_1d<double, 2> ddgama   = ZeroVector(2); 
		array_1d<double, 2> residual = ZeroVector(2);
		Matrix d                     = ZeroMatrix(2,2);
                Matrix d_inv                 = ZeroMatrix(2,2);
		
		

		
	        double sinphi             =   std::sin(PI * mMorhCoulomb->mcurrent_minternal_friction_angle  / 180.00);
		double cosphi             =   std::cos(PI * mMorhCoulomb->mcurrent_minternal_friction_angle  / 180.00);
		double sinpsi             =   std::sin(PI * mMorhCoulomb->mcurrent_dilatancy_angle  / 180.00);
		double cospsi             =   std::cos(PI * mMorhCoulomb->mcurrent_dilatancy_angle  / 180.00); 
		
		residual[0]   = (PrincipalStress[0] - PrincipalStress[2]) + (PrincipalStress[0] + PrincipalStress[2]) * sinphi - 2.00 *  cosphi * mMorhCoulomb->mcurrent_cohesion;
		residual[1]   =  PrincipalStress[0] - mRankine->mcurrent_Ft;

		
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
		double H     = 0.00;
		double& gama_a = dgama[0];
		double& gama_b = dgama[1];
		
		Vector Imput_Parameters_M;
		Vector Imput_Parameters_R;
		
		Imput_Parameters_M.resize(4); Imput_Parameters_M = ZeroVector(4);
		Imput_Parameters_R.resize(4); Imput_Parameters_R = ZeroVector(4);
		
		
		Imput_Parameters_M[0] =  mlength;
		Imput_Parameters_M[1] =  m_modified_morh_coulomb_maccumulated_plastic_strain_current;
		Imput_Parameters_M[2] =  mMorhCoulomb->mcurrent_minternal_friction_angle;
                Imput_Parameters_M[3] =  mpastic_damage_current;
		
		Imput_Parameters_R[0] =  mlength;
		Imput_Parameters_R[1] =  m_modified_morh_coulomb_maccumulated_plastic_strain_current; 
		Imput_Parameters_R[2] =  mpastic_damage_old;
		Imput_Parameters_R[3] =  mpastic_damage_current;
		
		Partial_Cohesion      = (mMorhCoulomb->mpSofteningBehavior_Cohesion)->FirstDerivateFunctionBehavior(Imput_Parameters_M);   
		Partial_Friction      = (mMorhCoulomb->mpSofteningBehavior_Friction)->FirstDerivateFunctionBehavior(Imput_Parameters_M);   
		Partial_Dilatancy     = (mMorhCoulomb->mpSofteningBehavior_Dilatancy)->FirstDerivateFunctionBehavior(Imput_Parameters_M);  
		H                     = (mRankine->mpSofteningBehaviorFt)->FirstDerivateFunctionBehavior(Imput_Parameters_R);
		
	        sinphi                =   std::sin(PI * mMorhCoulomb->mcurrent_minternal_friction_angle  / 180.00);
		cosphi                =   std::cos(PI * mMorhCoulomb->mcurrent_minternal_friction_angle  / 180.00);
		sinpsi                =   std::sin(PI * mMorhCoulomb->mcurrent_dilatancy_angle  / 180.00);
		cospsi                =   std::cos(PI * mMorhCoulomb->mcurrent_dilatancy_angle  / 180.00); 
		
		Vector Imput(3);
		Imput = ZeroVector(3);
		double Ppvs = 0.00;           /// principal plastic volumetric strain  
		array_1d<double,3> Ppds = ZeroVector(3);      /// principal plastic desviatoric strain
		array_1d<double,3> Pps  = ZeroVector(3);       /// principal plastic  strain
		array_1d<double,3> I; 
		I[0] = 1.00;
		I[1] = 1.00;
		I[2] = 1.00;
		
		while(fabs(norma)>toler && iter++ < max )
		{  
		  
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
		  
		  
		  d(1,0) = -2.00 * G * Partial_A_gama_a - K * Partial_D_gama_a - H * Partial_Ep_gama_a;  
		  
		  d(1,1) = -2.00 * G * Partial_A_gama_b - K * Partial_D_gama_b - H * Partial_Ep_gama_b; 
		  
		  singular       =  SD_MathUtils<double>::InvertMatrix(d, d_inv);
		  ddgama         =  -Vector(prod(d_inv, residual));
  
		  //Compute Newton-Raphson increment and update variables DGAMA and DGAMB
		  noalias(dgama) += ddgama;  
		  
		  /// volumetric and desviatoric plastic strain
		  Pps[0]        = (1.00 + sinpsi) * gama_a +  gama_b;
		  Pps[1]        = 0.00;
		  Pps[2]        = (sinpsi-1.00) * gama_a; 
		  Ppvs          = Pps[0] + Pps[1] + Pps[2]; 
		  noalias(Ppds) = Pps - d3 * Ppvs * I;
		  
		 //Updating Materiales
		 //aux_1   = (sinpsi+1.00) * gama_a + gama_b; aux_1 = aux_1 * aux_1;  
		 //aux_2   = (sinpsi-1.00) * gama_a;          aux_2 = aux_2 * aux_2;  
		 //aux     = raiz2d3 * std::sqrt(aux_1 + aux_2);
		 
		 aux  = raiz2d3* std::sqrt(inner_prod(Pps,Pps));
		 m_modified_morh_coulomb_maccumulated_plastic_strain_current = m_modified_morh_coulomb_maccumulated_plastic_strain_old + aux;
		 ComputeActualStrees(Ppvs, Ppds, PrincipalStress, Sigma);
		 ComputeActualStrain(Pps);
		 CalculatePlasticDamage(Sigma);
		
		 ///variable morh coulomb 
		Imput[0] = mpastic_damage_current;
		Imput[1] = mpastic_damage_old;
		Imput[2] = mMorhCoulomb->mcohesion; ///old cohesion
		
		mMorhCoulomb->mcurrent_cohesion = mMorhCoulomb->mpSofteningBehavior_Cohesion->EvolucionLaws(Imput, Sigma);

		// Compute Internal Variables   
		Imput_Parameters_M[1]             =  m_modified_morh_coulomb_maccumulated_plastic_strain_current;
		Imput_Parameters_M[3]             =  mpastic_damage_current;
		//mcurrent_cohesion               =  mpSofteningBehavior_Cohesion->FunctionBehavior(Imput_Parameters); 
		mMorhCoulomb->mcurrent_minternal_friction_angle =  mMorhCoulomb->mpSofteningBehavior_Friction->FunctionBehavior(Imput_Parameters_M);   
		Imput_Parameters_M[2]                           =  mMorhCoulomb->mcurrent_minternal_friction_angle;
		mMorhCoulomb->mcurrent_dilatancy_angle          =  mMorhCoulomb->mpSofteningBehavior_Dilatancy->FunctionBehavior(Imput_Parameters_M);
		 
		/// Variable Rankine
		///* Updatinf mFt
		Imput_Parameters_R[1] =  m_modified_morh_coulomb_maccumulated_plastic_strain_current;
		Imput_Parameters_R[2] =  mpastic_damage_old;
                Imput_Parameters_R[3] =  mpastic_damage_current;
		mRankine->mcurrent_Ft =  mRankine->mpSofteningBehaviorFt->FunctionBehavior(Imput_Parameters_R);
		
		
		/* 
	        Imput_Parameters_M[1]                           =  m_modified_morh_coulomb_maccumulated_plastic_strain_current;
		Imput_Parameters_R[1]                           =  m_modified_morh_coulomb_maccumulated_plastic_strain_current;
		mMorhCoulomb->mcurrent_cohesion                 =  mMorhCoulomb->mpSofteningBehavior_Cohesion->FunctionBehavior(Imput_Parameters_M);
		mMorhCoulomb->mcurrent_minternal_friction_angle =  mMorhCoulomb->mpSofteningBehavior_Friction->FunctionBehavior(Imput_Parameters_M); 
		Imput_Parameters_M[2]                           =  mMorhCoulomb->mcurrent_minternal_friction_angle;
		mMorhCoulomb->mcurrent_dilatancy_angle          =  mMorhCoulomb->mpSofteningBehavior_Dilatancy->FunctionBehavior(Imput_Parameters_M);
	        mRankine->mcurrent_Ft                           =  (mRankine->mpSofteningBehaviorFt)->FunctionBehavior(Imput_Parameters_R);
	       */
		
		Partial_Cohesion                                = (mMorhCoulomb->mpSofteningBehavior_Cohesion)->FirstDerivateFunctionBehavior(Imput_Parameters_M);   
		Partial_Friction                                = (mMorhCoulomb->mpSofteningBehavior_Friction)->FirstDerivateFunctionBehavior(Imput_Parameters_M);   
		Partial_Dilatancy                               = (mMorhCoulomb->mpSofteningBehavior_Dilatancy)->FirstDerivateFunctionBehavior(Imput_Parameters_M);  
		H                                               = (mRankine->mpSofteningBehaviorFt)->FirstDerivateFunctionBehavior(Imput_Parameters_R);
		 
		sinphi                                          =   std::sin(PI * mMorhCoulomb->mcurrent_minternal_friction_angle  / 180.00);
		cosphi                                          =   std::cos(PI * mMorhCoulomb->mcurrent_minternal_friction_angle  / 180.00);
		sinpsi                                          =   std::sin(PI * mMorhCoulomb->mcurrent_dilatancy_angle  / 180.00);
		cospsi                                          =   std::cos(PI * mMorhCoulomb->mcurrent_dilatancy_angle  / 180.00); 
		
		A = gama_a * (1.00 + (1.00/3.00) * sinpsi) + (2.00/3.00) * gama_b;
		B = -(2.00/3.00) * gama_a * sinpsi - (1.00/3.00) * gama_b;                 
		C = gama_a * ( (1.00/3.00) *sinpsi  - 1.00 ) - (1.00/3.00) * gama_b;
		D = 2.00 * gama_a * sinpsi + gama_b;
		 
		residual[0]   = (PrincipalStress[0] -  PrincipalStress[2]) + (PrincipalStress[0] + PrincipalStress[2]) * sinphi - 2.00 *  cosphi * mMorhCoulomb->mcurrent_cohesion;
		residual[1]   =  PrincipalStress[0]  - mRankine->mcurrent_Ft;
		
	        residual[0]   = residual[0]-2.00 * G * A + 2.00 * G * C - 2.00 * G * sinphi * A - 2.00 * G * sinphi * C - 2.00 * K * sinphi * D;  
		residual[1]   = residual[1]-2.00 * G * A - K * D; 
		
		
		norma         = norm_2(residual); 
		//std::cout<< "iter = "<< iter << std::endl;
		//KRATOS_WATCH(residual)
		//KRATOS_WATCH(Sigma)
                //mMorhCoulomb->CheckPlasticAdmisibility(Sigma);
		
		if(iter>=max){  

		  KRATOS_WATCH(norma) 
		  KRATOS_WATCH(residual)
		  KRATOS_WATCH(PrincipalStress)
		  KRATOS_WATCH(Sigma)
		  std::cout<< "RETURN MAPPING OF MAIN PLANE AND SIGMA 1 NOT CONVERGED"<< std::endl;
		  break;
		  //KRATOS_ERROR(std::logic_error,  "RETURN MAPPING OF MAIN PLANE AND SIGMA 1 NOT CONVERGED" , " "); 
		 }
		}

		/// volumetric and desviatoric plastic strain update again
		Pps[0]        = (1.00 + sinpsi) * gama_a +  gama_b;
		Pps[1]        = 0.00;
		Pps[2]        = (sinpsi-1.00) * gama_a; 
		Ppvs          = Pps[0] + Pps[1] + Pps[2]; 
		noalias(Ppds) = Pps - d3 * Ppvs * I;

		//Sigma[0] = PrincipalStress[0] -  2.00 * G * A - K * D; 
                //Sigma[1] = PrincipalStress[1] -  2.00 * G * B - K * D; 
                //Sigma[2] = PrincipalStress[2] -  2.00 * G * C - K * D; 
		noalias(Sigma) = PrincipalStress - 2.00 * G * Ppds - K * Ppvs * I;
		bool cond_a = mRankine->CheckValidity(Sigma);
		bool cond_b = mRankine->CheckPlasticAdmisibility(Sigma);     // false when no plastic
		bool cond_c = mMorhCoulomb->CheckPlasticAdmisibility(Sigma); // false when no plastic
		//std::cout<< "saliendo = "<< std::endl;
		//KRATOS_WATCH(cond_a)
		//KRATOS_WATCH(cond_b)
		//KRATOS_WATCH(cond_c)
		//KRATOS_WATCH(PrincipalStress)
		//KRATOS_WATCH(Sigma)
		//KRATOS_WATCH(mRankine->mcurrent_Ft)
		//CheckPlasticAdmisibility(Sigma);
		//KRATOS_ERROR(std::logic_error,  "RETURN MAPPING OF MAIN PLANE AND SIGMA 1 NOT CONVERGED" , " ");
		
		if(cond_a==true && cond_b==false && cond_c ==false)
		{
		    //KRATOS_WATCH("AAAAAAAAA")
		    Vector PPS_bar(3);
		    PPS_bar = ZeroVector(3);
                    ComputePlasticStrainBar(mplastic_strain_old ,m_inv_DeltaF, PPS_bar);
		    //updating the correct principal pastic strain
		    mPrincipalPlasticStrain_current[0] = /*mPrincipalPlasticStrain_old[0]*/ PPS_bar[0] + gama_a * (1.00 + sinpsi) + gama_b;
		    mPrincipalPlasticStrain_current[1] = /*mPrincipalPlasticStrain_current[1]*/ + PPS_bar[1]; 
		    mPrincipalPlasticStrain_current[2] = /*mPrincipalPlasticStrain_old[2]*/  PPS_bar[2] + gama_a * (sinpsi-1.00 ); 
		    return true;
		}
		else
		  return false;
		
		}
		
		
    bool Modified_Morh_Coulomb_Yield_Function::Return_Mapping_Intersection_Main_Plane_Corner_And_Sigma1_Tesile_Plane(const array_1d<double,3>& PrincipalStress, const array_1d<unsigned int,3>& order, array_1d<double,3>& Sigma)
    {
      
                UpdateMaterial();
		noalias(Sigma)               =   ZeroVector(3);
                const double& E              =   (*mprops)[YOUNG_MODULUS];
		const double d3              =   0.3333333333333333333; 
		const double& NU             =   (*mprops)[POISSON_RATIO];
		const double G               =   0.5 * E / (1.00 + NU);
		const double K               =   E / (3.00 * (1.00-2.00 * NU) );
		const double toler           =  1E-3;
		const unsigned max           =  1000;
		
		array_1d<double, 3> dgama    = ZeroVector(3);     
		array_1d<double, 3> ddgama   = ZeroVector(3); 
		array_1d<double, 3> residual = ZeroVector(3);
		Matrix d                     = ZeroMatrix(3,3);
                Matrix d_inv                 = ZeroMatrix(3,3);
		
		unsigned iter                =  0;
	        double sinphi                =   std::sin(PI * mMorhCoulomb->mcurrent_minternal_friction_angle  / 180.00);
		double cosphi                =   std::cos(PI * mMorhCoulomb->mcurrent_minternal_friction_angle  / 180.00);
		double sinpsi                =   std::sin(PI * mMorhCoulomb->mcurrent_dilatancy_angle  / 180.00);
		double cospsi                =   std::cos(PI * mMorhCoulomb->mcurrent_dilatancy_angle  / 180.00); 
		
		residual[0]   = (PrincipalStress[0] - PrincipalStress[2]) + (PrincipalStress[0] + PrincipalStress[2]) * sinphi - 2.00 *  cosphi * mMorhCoulomb->mcurrent_cohesion;
		residual[1]   = (PrincipalStress[0] - PrincipalStress[1]) + (PrincipalStress[0] + PrincipalStress[1]) * sinphi - 2.00 *  cosphi * mMorhCoulomb->mcurrent_cohesion;
		residual[2]   =  PrincipalStress[0] - mRankine->mcurrent_Ft;
		
		//KRATOS_WATCH(residual)
		
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
		double aux_1   = 0.00; 
		double aux_2   = 0.00;   
		double aux     = 1.00;
		double& gama_a = dgama[0];
		double& gama_b = dgama[1];
		double& gama_c = dgama[2];
		
		double H = 0.00;
		
		Vector Imput(3);
		Imput = ZeroVector(3);
		double Ppvs = 0.00;           /// principal plastic volumetric strain  
		array_1d<double,3> Ppds = ZeroVector(3);      /// principal plastic desviatoric strain
		array_1d<double,3> Pps  = ZeroVector(3);       /// principal plastic  strain
		array_1d<double,3> I; 
		I[0] = 1.00;
		I[1] = 1.00;
		I[2] = 1.00;
		
		
		Vector Imput_Parameters_M;
		Vector Imput_Parameters_R;

		Imput_Parameters_M.resize(4); Imput_Parameters_M = ZeroVector(4);
		Imput_Parameters_R.resize(4); Imput_Parameters_R = ZeroVector(4);


		Imput_Parameters_M[0] =  mlength;
		Imput_Parameters_M[1] =  m_modified_morh_coulomb_maccumulated_plastic_strain_current;
		Imput_Parameters_M[2] =  mMorhCoulomb->mcurrent_minternal_friction_angle;
                Imput_Parameters_M[3] =  mpastic_damage_current;
		
		Imput_Parameters_R[0] =  mlength;
		Imput_Parameters_R[1] =  m_modified_morh_coulomb_maccumulated_plastic_strain_current; 
		Imput_Parameters_R[2] =  mpastic_damage_old;
		Imput_Parameters_R[3] =  mpastic_damage_current;

		Partial_Cohesion      = (mMorhCoulomb->mpSofteningBehavior_Cohesion)->FirstDerivateFunctionBehavior(Imput_Parameters_M);   
		Partial_Friction      = (mMorhCoulomb->mpSofteningBehavior_Friction)->FirstDerivateFunctionBehavior(Imput_Parameters_M);   
		Partial_Dilatancy     = (mMorhCoulomb->mpSofteningBehavior_Dilatancy)->FirstDerivateFunctionBehavior(Imput_Parameters_M);  
		H                     = (mRankine->mpSofteningBehaviorFt)->FirstDerivateFunctionBehavior(Imput_Parameters_R);
		
		sinphi             =   std::sin(PI * mMorhCoulomb->mcurrent_minternal_friction_angle  / 180.00);
		cosphi             =   std::cos(PI * mMorhCoulomb->mcurrent_minternal_friction_angle  / 180.00);
		sinpsi             =   std::sin(PI * mMorhCoulomb->mcurrent_dilatancy_angle  / 180.00);
		cospsi             =   std::cos(PI * mMorhCoulomb->mcurrent_dilatancy_angle  / 180.00); 
		
		while(fabs(norma)>toler && iter++ < max )
		{

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
		 // aux_1   = (sinpsi+1.00) * (gama_a + gama_b) + gama_c  ; aux_1 = aux_1 * aux_1;  
		 //aux_2   = (sinpsi-1.00) *(sinpsi-1.00) * ( gama_a * gama_a + gama_b*gama_b) ;  
		 //aux     = raiz2d3 * std::sqrt(aux_1 + aux_2);
		 //m_modified_morh_coulomb_maccumulated_plastic_strain_current = m_modified_morh_coulomb_maccumulated_plastic_strain_old + aux;
		  
		 /// volumetric and desviatoric plastic strain
		  
		  Pps[0]        = (1.00 + sinpsi) * (gama_a + gama_b) + gama_c;
		  Pps[1]        = (sinpsi-1.00) * gama_b;
		  Pps[2]        = (sinpsi-1.00) * gama_a;
		  Ppvs          = Pps[0] + Pps[1] + Pps[2]; 
		  noalias(Ppds) = Pps - d3 * Ppvs * I;
		  
		 
		 aux  = raiz2d3* std::sqrt(inner_prod(Pps,Pps));
		 m_modified_morh_coulomb_maccumulated_plastic_strain_current = m_modified_morh_coulomb_maccumulated_plastic_strain_old + aux;
		 ComputeActualStrees(Ppvs, Ppds, PrincipalStress, Sigma);
		 ComputeActualStrain(Pps);
		 CalculatePlasticDamage(Sigma);
		
		 ///variable morh coulomb 
		Imput[0] = mpastic_damage_current;
		Imput[1] = mpastic_damage_old;
		Imput[2] = mMorhCoulomb->mcohesion; ///old cohesion
		
		mMorhCoulomb->mcurrent_cohesion = mMorhCoulomb->mpSofteningBehavior_Cohesion->EvolucionLaws(Imput, Sigma);

		// Compute Internal Variables   
		Imput_Parameters_M[1]             =  m_modified_morh_coulomb_maccumulated_plastic_strain_current;
		Imput_Parameters_M[3]             =  mpastic_damage_current;
		//mcurrent_cohesion               =  mpSofteningBehavior_Cohesion->FunctionBehavior(Imput_Parameters); 
		mMorhCoulomb->mcurrent_minternal_friction_angle =  mMorhCoulomb->mpSofteningBehavior_Friction->FunctionBehavior(Imput_Parameters_M);   
		Imput_Parameters_M[2]             =  mMorhCoulomb->mcurrent_minternal_friction_angle;
		mMorhCoulomb->mcurrent_dilatancy_angle          =  mMorhCoulomb->mpSofteningBehavior_Dilatancy->FunctionBehavior(Imput_Parameters_M);
		 
		/// Variable Rankine
		///* Updatinf mFt
		Imput_Parameters_R[1] =  m_modified_morh_coulomb_maccumulated_plastic_strain_current;
	        Imput_Parameters_R[2] =  mpastic_damage_old;
                Imput_Parameters_R[3] =  mpastic_damage_current;
		mRankine->mcurrent_Ft =  mRankine->mpSofteningBehaviorFt->FunctionBehavior(Imput_Parameters_R);
		
		  /*
		  Imput_Parameters_M[1]                           =  m_modified_morh_coulomb_maccumulated_plastic_strain_current;
		  Imput_Parameters_R[1]                           =  m_modified_morh_coulomb_maccumulated_plastic_strain_current;
		  mMorhCoulomb->mcurrent_cohesion                 =  mMorhCoulomb->mpSofteningBehavior_Cohesion->FunctionBehavior(Imput_Parameters_M);
		  mMorhCoulomb->mcurrent_minternal_friction_angle =  mMorhCoulomb->mpSofteningBehavior_Friction->FunctionBehavior(Imput_Parameters_M); 
		  Imput_Parameters_M[2]                           =  mMorhCoulomb->mcurrent_minternal_friction_angle;
		  mMorhCoulomb->mcurrent_dilatancy_angle          =  mMorhCoulomb->mpSofteningBehavior_Dilatancy->FunctionBehavior(Imput_Parameters_M);
		  mRankine->mcurrent_Ft                           =  (mRankine->mpSofteningBehaviorFt)->FunctionBehavior(Imput_Parameters_R);
		  */
		  
		  Partial_Cohesion                  = (mMorhCoulomb->mpSofteningBehavior_Cohesion)->FirstDerivateFunctionBehavior(Imput_Parameters_M);   
		  Partial_Friction                  = (mMorhCoulomb->mpSofteningBehavior_Friction)->FirstDerivateFunctionBehavior(Imput_Parameters_M);   
		  Partial_Dilatancy                 = (mMorhCoulomb->mpSofteningBehavior_Dilatancy)->FirstDerivateFunctionBehavior(Imput_Parameters_M);  
		  H                                 = (mRankine->mpSofteningBehaviorFt)->FirstDerivateFunctionBehavior(Imput_Parameters_R);
		 
		  sinphi             =   std::sin(PI * mMorhCoulomb->mcurrent_minternal_friction_angle  / 180.00);
		  cosphi             =   std::cos(PI * mMorhCoulomb->mcurrent_minternal_friction_angle  / 180.00);
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

		  KRATOS_WATCH(norma)
		  KRATOS_WATCH(iter)
		  std::cout<< "RETURN MAPPING OF MAIN PLANE CORNER AND SIGMA 1 NOT CONVERGED" << std::endl; 
		  break;
		  //KRATOS_ERROR(std::logic_error,  "RETURN MAPPING OF MAIN PLANE CORNER AND SIGMA 1 NOT CONVERGED" , ""); 
		 }
		}
		
		//KRATOS_WATCH(norma)
		//KRATOS_WATCH(Sigma)
		/// volumetric and desviatoric plastic strain update again
		  Pps[0]        = (1.00 + sinpsi) * (gama_a + gama_b) + gama_c;
		  Pps[1]        = (sinpsi-1.00) * gama_b;
		  Pps[2]        = (sinpsi-1.00) * gama_a;
		  Ppvs          = Pps[0] + Pps[1] + Pps[2]; 
		  noalias(Ppds) = Pps - d3 * Ppvs * I;

		//Sigma[0] = PrincipalStress[0] -  2.00 * G * A - K * D; 
                //Sigma[1] = PrincipalStress[1] -  2.00 * G * B - K * D; 
                //Sigma[2] = PrincipalStress[2] -  2.00 * G * C - K * D; 
		noalias(Sigma) = PrincipalStress - 2.00 * G * Ppds - K * Ppvs * I;
		bool cond_a = mRankine->CheckValidity(Sigma);
		bool cond_b = mRankine->CheckPlasticAdmisibility(Sigma); // false when no plastic
		bool cond_c = mMorhCoulomb->CheckPlasticAdmisibility(Sigma);

		//KRATOS_WATCH(PrincipalStress)
		//KRATOS_WATCH(Sigma)
		//KRATOS_WATCH(norma)
		//KRATOS_WATCH(Sigma[0]-mRankine->mcurrent_Ft)
		//KRATOS_WATCH(Sigma[1]-mRankine->mcurrent_Ft)
		//KRATOS_WATCH(Sigma[2]-mRankine->mcurrent_Ft)
		//CheckPlasticAdmisibility(Sigma);
		
		
		if(cond_a==true && cond_b==false && cond_c ==false)
		  {
		    //updating the correct principal pastic strain
		    Vector PPS_bar(3);
		    PPS_bar = ZeroVector(3);
                    ComputePlasticStrainBar(mplastic_strain_old ,m_inv_DeltaF, PPS_bar);
		    mPrincipalPlasticStrain_current[0] = /*mPrincipalPlasticStrain_old[0]*/ PPS_bar[0] +  (gama_a + gama_b ) * (1.00 + sinpsi) + gama_c;
		    mPrincipalPlasticStrain_current[1] = /*mPrincipalPlasticStrain_old[1]*/ PPS_bar[1] +  (sinpsi -1.00) * gama_b ;
		    mPrincipalPlasticStrain_current[2] = /*mPrincipalPlasticStrain_old[2]*/ PPS_bar[2] +  (sinpsi -1.00) * gama_a;
		    return true;
		  }
		else
		  return false;
    }
    
    
    bool Modified_Morh_Coulomb_Yield_Function::Return_Mapping_Intersection_Main_Plane_Corner_And_Sigma1_And_Sigma_2_Tesile_Plane(const array_1d<double,3>& PrincipalStress, const array_1d<unsigned int,3>& order, array_1d<double,3>& Sigma)
    {
      
	
      	UpdateMaterial();
        noalias(Sigma)               =   ZeroVector(3);
        const double& E              =   (*mprops)[YOUNG_MODULUS];
	const double& NU             =   (*mprops)[POISSON_RATIO];
	const double G               =   0.5 * E / (1.00 + NU);
	const double K               =   E / (3.00 * (1.00-2.00 * NU) );
	const double toler           =   1.00;
	const unsigned max           =   1000;
	const double d3              =   0.33333333333333333;

	array_1d<double, 4> dgama    =   ZeroVector(4);     
	array_1d<double, 4> ddgama   =   ZeroVector(4); 
	array_1d<double, 4> residual =   ZeroVector(4);
        
	Matrix d; d.resize(4,4, false);
	Matrix d_inv; d_inv.resize(4,4,false); 
	noalias(d)                   =   ZeroMatrix(4,4);
	noalias(d_inv)               =   ZeroMatrix(4,4);

	unsigned iter                =   0;
	double sinphi                =   std::sin(PI * mMorhCoulomb->mcurrent_minternal_friction_angle  / 180.00);
	double cosphi                =   std::cos(PI * mMorhCoulomb->mcurrent_minternal_friction_angle  / 180.00);
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
        double H       = 0.00;

	Vector Imput_Parameters_M;
	Vector Imput_Parameters_R;
		
	Imput_Parameters_M.resize(4);Imput_Parameters_M = ZeroVector(4);
	Imput_Parameters_R.resize(4);Imput_Parameters_R = ZeroVector(4);


	Imput_Parameters_M[0] =  mlength;
	Imput_Parameters_M[1] =  m_modified_morh_coulomb_maccumulated_plastic_strain_current;
	Imput_Parameters_M[2] =  mMorhCoulomb->mcurrent_minternal_friction_angle;
	Imput_Parameters_M[3] =  mpastic_damage_current;

	Imput_Parameters_R[0] =  mlength;
	Imput_Parameters_R[1] =  m_modified_morh_coulomb_maccumulated_plastic_strain_current; 
	Imput_Parameters_R[2] =  mpastic_damage_old;
	Imput_Parameters_R[3] =  mpastic_damage_current;

	Partial_Cohesion      = (mMorhCoulomb->mpSofteningBehavior_Cohesion)->FirstDerivateFunctionBehavior(Imput_Parameters_M);   
	Partial_Friction      = (mMorhCoulomb->mpSofteningBehavior_Friction)->FirstDerivateFunctionBehavior(Imput_Parameters_M);   
	Partial_Dilatancy     = (mMorhCoulomb->mpSofteningBehavior_Dilatancy)->FirstDerivateFunctionBehavior(Imput_Parameters_M);  
	H                     = (mRankine->mpSofteningBehaviorFt)->FirstDerivateFunctionBehavior(Imput_Parameters_R);
	
	sinphi             =   std::sin(PI * mMorhCoulomb->mcurrent_minternal_friction_angle  / 180.00);
	cosphi             =   std::cos(PI * mMorhCoulomb->mcurrent_minternal_friction_angle  / 180.00);
	sinpsi             =   std::sin(PI * mMorhCoulomb->mcurrent_dilatancy_angle  / 180.00);
	cospsi             =   std::cos(PI * mMorhCoulomb->mcurrent_dilatancy_angle  / 180.00); 
	
	Vector Imput(3);
	Imput = ZeroVector(3);
	double Ppvs = 0.00;           /// principal plastic volumetric strain  
	array_1d<double,3> Ppds = ZeroVector(3);      /// principal plastic desviatoric strain
	array_1d<double,3> Pps  = ZeroVector(3);       /// principal plastic  strain
	array_1d<double,3> I; 
	I[0] = 1.00;
	I[1] = 1.00;
	I[2] = 1.00;
	
	int count = 0;
	while(fabs(norma)>toler && iter++ < max )
		{
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
	
		
		  
		//if(count==0)
		{    
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
		  
		  
		  d(2,0) = -2.00 * G * Partial_A_gama_a - K * Partial_D_gama_a - H * Partial_D_gama_a;
		  d(2,1) = -2.00 * G * Partial_A_gama_b - K * Partial_D_gama_b - H * Partial_D_gama_b;
		  d(2,2) = -2.00 * G * Partial_A_gama_c - K * Partial_D_gama_c - H * Partial_D_gama_c;
		  d(2,3) = -2.00 * G * Partial_A_gama_d - K * Partial_D_gama_d - H * Partial_D_gama_d;
		  
		  d(3,0) = -2.00 * G * Partial_B_gama_a - K * Partial_D_gama_a- H * Partial_D_gama_a;
		  d(3,1) = -2.00 * G * Partial_B_gama_b - K * Partial_D_gama_b- H * Partial_D_gama_b;
		  d(3,2) = -2.00 * G * Partial_B_gama_c - K * Partial_D_gama_c- H * Partial_D_gama_c;
		  d(3,3) = -2.00 * G * Partial_B_gama_d - K * Partial_D_gama_d- H * Partial_D_gama_d;
		  
		  /*
		  // Trunca los termnos de esta matriz. Lo hice con el proposito de que esta matriz tubieta 
		  // inversa. Al parecer no tiene inversa cuando utilizo la forma d como tal.
		  for( unsigned int i = 0; i< 4; i++){
		      for( unsigned int j = 0; j< 4; j++)
		      {
			aux    = d(i,j);
			d(i,j) = round(aux, 2);
		      }
		  }
		   count++;
		   */
		   singular       =   SD_MathUtils<double>::InvertMatrix(d, d_inv);
		  }
		  

		  ddgama         =  -Vector(prod(d_inv, residual));
                  
		  //Compute Newton-Raphson increment and update variables DGAMA and DGAMB
		  noalias(dgama) += ddgama;  
		  
		 //Updating Materiales
		 //aux_1   = (sinpsi+1.00) *gama_a + gama_c  ;  aux_1 = aux_1 * aux_1;  
		 //aux_2   = (sinpsi+1.00) *gama_b + gama_d  ;  aux_2 = aux_2 * aux_2;
		 //aux_3   = (sinpsi-1.00) *( gama_a + gama_b); aux_3 = aux_3 * aux_3;
		 //aux     = raiz2d3 * std::sqrt(aux_1 + aux_2 + aux_3);		 
		 //m_modified_morh_coulomb_maccumulated_plastic_strain_current = m_modified_morh_coulomb_maccumulated_plastic_strain_old + aux;

		 /// volumetric and desviatoric plastic strain
		  Pps[0]        = (1.00 + sinpsi) * (gama_a) + gama_c;
		  Pps[1]        = (1.00 + sinpsi) * (gama_b) + gama_d;
		  Pps[2]        = (sinpsi-1.00) * ( gama_a + gama_b);
		  Ppvs          = Pps[0] + Pps[1] + Pps[2]; 
		  noalias(Ppds) = Pps - d3 * Ppvs * I;

		  
		  aux  = raiz2d3* std::sqrt(inner_prod(Pps,Pps));
		  m_modified_morh_coulomb_maccumulated_plastic_strain_current = m_modified_morh_coulomb_maccumulated_plastic_strain_old + aux;
		  ComputeActualStrees(Ppvs, Ppds, PrincipalStress, Sigma);
		  ComputeActualStrain(Pps);
		  CalculatePlasticDamage(Sigma);
		  
		  //KRATOS_WATCH(d)
		  //KRATOS_WATCH(d_inv)
		  //KRATOS_WATCH(PrincipalStress)
		  //KRATOS_WATCH(Sigma)
		  //KRATOS_WATCH(mpastic_damage_current)
		  //KRATOS_WATCH(mpastic_damage_old)
		  //KRATOS_WATCH("------------------------")
		  
		  ///variable morh coulomb 
		  Imput[0] = mpastic_damage_current;
		  Imput[1] = mpastic_damage_old;
		  Imput[2] = mMorhCoulomb->mcohesion; ///old cohesion

		  mMorhCoulomb->mcurrent_cohesion = mMorhCoulomb->mpSofteningBehavior_Cohesion->EvolucionLaws(Imput, Sigma);

		  // Compute Internal Variables   
		  Imput_Parameters_M[1]             =  m_modified_morh_coulomb_maccumulated_plastic_strain_current;
		  Imput_Parameters_M[3]             =  mpastic_damage_current;
		  //mcurrent_cohesion               =  mpSofteningBehavior_Cohesion->FunctionBehavior(Imput_Parameters); 
		  mMorhCoulomb->mcurrent_minternal_friction_angle =  mMorhCoulomb->mpSofteningBehavior_Friction->FunctionBehavior(Imput_Parameters_M);   
		  Imput_Parameters_M[2]             =  mMorhCoulomb->mcurrent_minternal_friction_angle;
		  mMorhCoulomb->mcurrent_dilatancy_angle          =  mMorhCoulomb->mpSofteningBehavior_Dilatancy->FunctionBehavior(Imput_Parameters_M);

		  /// Variable Rankine
		  ///* Updatinf mFt
		  Imput_Parameters_R[1]           =  m_modified_morh_coulomb_maccumulated_plastic_strain_current;
		  Imput_Parameters_R[2]           =  mpastic_damage_old;
                  Imput_Parameters_R[3]           =  mpastic_damage_current;
		  mRankine->mcurrent_Ft           =  mRankine->mpSofteningBehaviorFt->FunctionBehavior(Imput_Parameters_R);
		 
		  /*
		  Imput_Parameters_M[1]                           =  m_modified_morh_coulomb_maccumulated_plastic_strain_current;
		  Imput_Parameters_R[1]                           =  m_modified_morh_coulomb_maccumulated_plastic_strain_current;
		  mMorhCoulomb->mcurrent_cohesion                 =  mMorhCoulomb->mpSofteningBehavior_Cohesion->FunctionBehavior(Imput_Parameters_M);
		  mMorhCoulomb->mcurrent_minternal_friction_angle =  mMorhCoulomb->mpSofteningBehavior_Friction->FunctionBehavior(Imput_Parameters_M); 
		  Imput_Parameters_M[2]                           = mMorhCoulomb-> mcurrent_minternal_friction_angle;
		  mMorhCoulomb->mcurrent_dilatancy_angle          =  mMorhCoulomb->mpSofteningBehavior_Dilatancy->FunctionBehavior(Imput_Parameters_M);
		  mRankine->mcurrent_Ft                           =  (mRankine->mpSofteningBehaviorFt)->FunctionBehavior(Imput_Parameters_R);
		  */
		  
		  Partial_Cohesion                  = (mMorhCoulomb->mpSofteningBehavior_Cohesion)->FirstDerivateFunctionBehavior(Imput_Parameters_M);   
		  Partial_Friction                  = (mMorhCoulomb->mpSofteningBehavior_Friction)->FirstDerivateFunctionBehavior(Imput_Parameters_M);   
		  Partial_Dilatancy                 = (mMorhCoulomb->mpSofteningBehavior_Dilatancy)->FirstDerivateFunctionBehavior(Imput_Parameters_M);  
		  H                                 = (mRankine->mpSofteningBehaviorFt)->FirstDerivateFunctionBehavior(Imput_Parameters_R);
		 
		  sinphi             =   std::sin(PI * mMorhCoulomb->mcurrent_minternal_friction_angle  / 180.00);
		  cosphi             =   std::cos(PI * mMorhCoulomb->mcurrent_minternal_friction_angle  / 180.00);
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
	          
		if(iter>=max){ //|| gama_a < 0.00 || gama_b < 0.00 || gama_c < 0.00 || gama_d < 0.00 ){
		  KRATOS_WATCH(Sigma)
		  KRATOS_WATCH(norma)
		  KRATOS_WATCH(dgama)
		  std::cout<< "RETURN MAPPING OF MAIN PLANE AND SIGMA 1 AND SIGMA 2 NOT CONVERGED" << std::endl;
		  break;
		  //KRATOS_ERROR(std::logic_error,  "RETURN MAPPING OF MAIN PLANE AND SIGMA 1 AND SIGMA 2 NOT CONVERGED" , ""); 
		 }
		}

		//Sigma[0] = PrincipalStress[0] -  2.00 * G * A - K * D; 
                //Sigma[1] = PrincipalStress[1] -  2.00 * G * B - K * D; 
                //Sigma[2] = PrincipalStress[2] -  2.00 * G * C - K * D; 
		/// volumetric and desviatoric plastic strain update again
		Pps[0]        = (1.00 + sinpsi) * (gama_a) + gama_c;
		Pps[1]        = (1.00 + sinpsi) * (gama_b) + gama_d;
		Pps[2]        = (sinpsi-1.00) * ( gama_a + gama_b);
		Ppvs          = Pps[0] + Pps[1] + Pps[2]; 
		noalias(Ppds) = Pps - d3 * Ppvs * I;
		
		noalias(Sigma) = PrincipalStress - 2.00 * G * Ppds - K * Ppvs * I;
		
		bool cond_a = mRankine->CheckValidity(Sigma);
		bool cond_b = mRankine->CheckPlasticAdmisibility(Sigma); // false when no plastic
		bool cond_c = mMorhCoulomb->CheckPlasticAdmisibility(Sigma);
		
		//KRATOS_WATCH(PrincipalStress)
		//KRATOS_WATCH(Sigma)
		//KRATOS_WATCH(Sigma[0]-mRankine->mcurrent_Ft)
		//KRATOS_WATCH(Sigma[1]-mRankine->mcurrent_Ft)
		//KRATOS_WATCH(Sigma[2]-mRankine->mcurrent_Ft)
		//CheckPlasticAdmisibility(Sigma);
		
		if(cond_a==true && cond_b==false && cond_c ==false)
		   {
		    //updating the correct principal pastic strain
		    Vector PPS_bar(3);
		    PPS_bar = ZeroVector(3);
                    ComputePlasticStrainBar(mplastic_strain_old ,m_inv_DeltaF, PPS_bar);
		    mPrincipalPlasticStrain_current[0] = /*mPrincipalPlasticStrain_old[0]*/ PPS_bar[0] + gama_a  * (1.00 + sinpsi) + gama_c;
		    mPrincipalPlasticStrain_current[1] = /*mPrincipalPlasticStrain_old[1]*/ PPS_bar[1] + gama_b  * (1.00 + sinpsi) + gama_d;
		    mPrincipalPlasticStrain_current[2] = /*mPrincipalPlasticStrain_old[2]*/ PPS_bar[2] + (sinpsi -1.00) * ( gama_a + gama_d);
		    return true;
		   }
		else{
		  KRATOS_WATCH(PrincipalStress)
		  KRATOS_WATCH(Sigma)
		  KRATOS_ERROR(std::logic_error,  "RETURN MAPPING OF MAIN PLANE AND SIGMA 1 AND SIGMA 2 IS FALSE" , ""); 
		  return false;
		}
	
      
    }
    

    void Modified_Morh_Coulomb_Yield_Function::GetValue(const Variable<double>& rVariable, double& Result)
      {
	
	if(rVariable==COHESION)
	{Result = mMorhCoulomb->mcurrent_cohesion;}
        
        if(rVariable==FT){ 
	  Result = mRankine->mcurrent_Ft;}
	
	if(rVariable == DILATANCY_ANGLE)
	{Result = mMorhCoulomb->mcurrent_dilatancy_angle; }
	
	if(rVariable == INTERNAL_FRICTION_ANGLE)
	{Result = mMorhCoulomb->mcurrent_minternal_friction_angle;}
	
        if(rVariable == DAMAGE)
	{Result = mpastic_damage_current;}
	
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
            
        
        bool Modified_Morh_Coulomb_Yield_Function::CheckPlasticAdmisibility(const Vector& Stress)
        { 
	  	const double& friction  = mMorhCoulomb->mcurrent_minternal_friction_angle; 
	        const double& cohe      = mMorhCoulomb->mcurrent_cohesion;
	        const double sinphi     = std::sin(PI * friction  / 180.00);
	        const double cosphi     = std::cos(PI * friction  / 180.00);

	        // Check plastic admissibility
	        double sigma_ef = (Stress[0] - Stress[2]) + (Stress[0] + Stress[2]) * sinphi;
	        double phia     = sigma_ef - 2.00 *  cosphi * cohe;
                return false;
	  /*
	  bool plastic_1 = mMorhCoulomb->CheckPlasticAdmisibility(Stress); 
	  bool plastic_2 = mRankine->CheckPlasticAdmisibility(Stress); 
	  if(plastic_1== false && plastic_2==false)
	     return false; 
	  else
	    return true;
	  */
	}
	
	
      void Modified_Morh_Coulomb_Yield_Function::GetValue(const Variable<Vector>& rVariable, Vector& Result)
      {
      const int size = mplastic_strain.size();
      if(rVariable==ALMANSI_PLASTIC_STRAIN){
	    Result.resize(size);
	    noalias(Result) = mplastic_strain;
          }
      if(rVariable==ALMANSI_ELASTIC_STRAIN){
	    Result.resize(size);
	    noalias(Result) = mElastic_strain;
          }
      }
     
      void Modified_Morh_Coulomb_Yield_Function::GetValue(Matrix& Result)
      {
	m_inv_DeltaF;
	m_inv_DeltaF.resize(3,3, false);
	noalias(m_inv_DeltaF) = ZeroMatrix(3,3);
	switch(mState)
         {
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
	 
	  mRankine->GetValue(m_inv_DeltaF);
	  mMorhCoulomb->GetValue(m_inv_DeltaF);
      }
	
            
        void Modified_Morh_Coulomb_Yield_Function::CalculatePlasticDamage(const array_1d<double,3>& Sigma)
        {
	 
	  /*
	  const double& Ft   = (*mprops)[FT];
	  const double& Ec   = (*mprops)[YOUNG_MODULUS];
	  const double& GE   = (*mprops)[FRACTURE_ENERGY];
	  const double Eu    =  (2.00 * GE)/(Ft * mlength);
	  mpastic_damage_current = m_modified_morh_coulomb_maccumulated_plastic_strain_current/Eu;
	  if(mpastic_damage_current>=1.00) 
	      mpastic_damage_current = 1.00; 
	  */
	  
	            const double toler = 1E-6;
	  double teta_a     =  Tensor_Utils<double>::Mc_aully(Sigma);
          double teta_b     =  std::fabs(Sigma[0]) + std::fabs(Sigma[1]) + std::fabs(Sigma[2]);
	  double teta       =  0.00;
          array_1d<double, 3> DeltaPlasticStrain = mPrincipalPlasticStrain_current - mPrincipalPlasticStrain_old;
	  double disipation =  inner_prod(Sigma, DeltaPlasticStrain);
	  
	  
	  if (fabs(teta_b) > 1E-8)
          {
	   teta = teta_a/teta_b;
	   // computing Kp_punto
	   double gc_p            = (*mprops)[CRUSHING_ENERGY]/mlength; 
	   double gf_p            = (*mprops)[FRACTURE_ENERGY]/mlength;
           double h               = (teta/gf_p + (1.00-teta)/gc_p);
	   double kp_punto        = h * disipation;
	   if(disipation > 0.00)
	   {
	    mpastic_damage_current =  mpastic_damage_old + kp_punto;
	    if(mpastic_damage_current > 1.00)
	      mpastic_damage_current = 1.00;
	    }  
	  }
	}
      

void Modified_Morh_Coulomb_Yield_Function::ComputeActualStrees(const double& Ppvs, 
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

void Modified_Morh_Coulomb_Yield_Function::ComputeActualStrain(const array_1d<double,3>& Pps)
{
  Vector PPS_bar(3);
  PPS_bar =ZeroVector(3);
  ComputePlasticStrainBar(mplastic_strain_old ,m_inv_DeltaF, PPS_bar);  
  noalias(mPrincipalPlasticStrain_current) = PPS_bar + Pps;
}


//***********************************************************************
//***********************************************************************
      

  }
       




