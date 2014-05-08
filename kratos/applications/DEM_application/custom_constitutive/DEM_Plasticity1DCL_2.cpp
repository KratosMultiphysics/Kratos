
#include "DEM_constitutive_law.h"
#include "DEM_Plasticity1DCL.h"

namespace Kratos
{
                                                                                                      
      void CalculateAll(array_1d<double, 3>& rContactForce, array_1d<double, 3>& rContactMoment,array_1d<double, 3>& rElasticForce, array_1d<double, 3>& rInitialRotaMoment,ProcessInfo& rCurrentProcessInfo){


   KRATOS_TRY

        const double dt = rCurrentProcessInfo[DELTA_TIME];        
        const double dt_i = 1 / dt; 
        const int time_steps = rCurrentProcessInfo[TIME_STEPS];
                
        /* Initializations */
                          
        const array_1d<double, 3>& vel          = this->GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY);
        const array_1d<double, 3>& delta_displ  = this->GetGeometry()(0)->FastGetSolutionStepValue(DELTA_DISPLACEMENT);
        const array_1d<double, 3>& ang_vel      = this->GetGeometry()(0)->FastGetSolutionStepValue(ANGULAR_VELOCITY);
        double& rRepresentative_Volume          = this->GetGeometry()(0)->FastGetSolutionStepValue(REPRESENTATIVE_VOLUME);
        const double moment_of_inertia          = this->GetGeometry()(0)->FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA);
        double RotaAcc[3]                       = {0.0};

        // if (mRotationOption){
        if (this->Is(DEMFlags::HAS_ROTATION) ){    
            RotaAcc[0]                         = ang_vel[0] * dt_i;
            RotaAcc[1]                         = ang_vel[1] * dt_i;
            RotaAcc[2]                         = ang_vel[2] * dt_i;

            rInitialRotaMoment[0] = RotaAcc[0] * moment_of_inertia;       
            rInitialRotaMoment[1] = RotaAcc[1] * moment_of_inertia;
            rInitialRotaMoment[2] = RotaAcc[2] * moment_of_inertia;

        }        

        //size_t i_neighbour_count = 0;       
        
        //r(ParticleWeakIteratorType neighbour_iterator = mrNeighbours.begin(); neighbour_iterator != mrNeighbours.end(); neighbour_iterator++) {
        for( unsigned int i_neighbour_count = 0; i_neighbour_count < mNeighbourElements.size(); i_neighbour_count++) {
            SphericParticle* neighbour_iterator = mNeighbourElements[i_neighbour_count];
            //if(mNeighbourElements[i]->Id() != neighbour_iterator->Id() ) KRATOS_ERROR(std::logic_error,"ALAAAAAAAAAA",this->Id());
        
            array_1d<double,3> other_to_me_vect   = this->GetGeometry()(0)->Coordinates() - neighbour_iterator->GetGeometry()(0)->Coordinates();
            //const double &other_radius                  = neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(RADIUS);
            const double &other_radius                  = neighbour_iterator->GetRadius();
            //const double &other_sqrt_of_mass            = neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(SQRT_OF_MASS);           
            const double &other_sqrt_of_mass            = neighbour_iterator->GetSqrtOfRealMass();    
 
            double distance                       = sqrt(other_to_me_vect[0] * other_to_me_vect[0] +
                                                          other_to_me_vect[1] * other_to_me_vect[1] +
                                                          other_to_me_vect[2] * other_to_me_vect[2]);
            double radius_sum                     = mRadius + other_radius;
            double radius_sum_i                   = 1.0 / radius_sum;
            double equiv_radius                   = 2.0 * mRadius * other_radius * radius_sum_i; 
            
            double initial_delta                  = mNeighbourDelta[i_neighbour_count]; //*
            double initial_dist                   = (radius_sum - initial_delta);
            double initial_dist_i                 = 1.0 / initial_dist;
            double indentation                    = initial_dist - distance;   //#1
            double equiv_area                     = 0.25*M_PI * equiv_radius * equiv_radius; //#2 
            double calculation_area               = equiv_area;
            double equiv_mass                     = mSqrtOfRealMass * other_sqrt_of_mass;
            double myYoung                        = GetYoung();
            double myPoisson                      = GetPoisson();
            double myLnOfRestitCoeff              = GetLnOfRestitCoeff();
            double myTgOfFrictionAngle            = GetTgOfFrictionAngle();

            double equiv_young;
            double equiv_poisson;
            double equiv_visco_damp_coeff_normal;
            double equiv_visco_damp_coeff_tangential;
            double equiv_ln_of_restit_coeff;
            double kn_el;
            double kt_el;
            double equiv_tg_of_fri_ang;

            double DeltDisp[3]                    = {0.0};
            double RelVel[3]                      = {0.0};

            double LocalCoordSystem[3][3]         = {{0.0}, {0.0}, {0.0}};
            double OldLocalCoordSystem[3][3]      = {{0.0}, {0.0}, {0.0}};

            bool sliding = false;

            int mapping_new_ini = mMapping_New_Ini[i_neighbour_count]; //*
            int mapping_new_cont =mMapping_New_Cont[i_neighbour_count];

            double contact_tau = 0.0;
            double contact_sigma = 0.0;
            double failure_criterion_state = 0.0; 
            double acumulated_damage = 0.0; 
            
            unsigned int neighbour_iterator_id = neighbour_iterator->Id();

            // Getting neighbour properties
            double other_young               = neighbour_iterator->GetYoung();
            double other_poisson             = neighbour_iterator->GetPoisson();
            double other_ln_of_restit_coeff  = neighbour_iterator->GetLnOfRestitCoeff();
            double other_tg_of_fri_angle     = neighbour_iterator->GetTgOfFrictionAngle();

            equiv_young                       = 2.0 * myYoung * other_young / (myYoung + other_young);
            equiv_poisson                     = 2.0 * myPoisson * other_poisson / (myPoisson + other_poisson);
            equiv_ln_of_restit_coeff          = 0.5 * (myLnOfRestitCoeff + other_ln_of_restit_coeff);
            equiv_tg_of_fri_ang               = 0.5 * (myTgOfFrictionAngle + other_tg_of_fri_angle);
        
            double aux_norm_to_tang = 0.0;


              double rmin = mRadius;
              if(other_radius<mRadius) rmin = other_radius;
              
              calculation_area = M_PI*rmin*rmin;
              double equiv_shear = equiv_young/(2.0*(1+equiv_poisson));
              
              kn_el = equiv_young*calculation_area*initial_dist_i;
              kt_el = equiv_shear*calculation_area*initial_dist_i;

            

        
            //if (mCriticalTimeOption){
            if ( this->Is(DEMFlags::HAS_CRITICAL_TIME) ){
                double historic = rCurrentProcessInfo[HISTORICAL_MIN_K];

                if ((kn_el < historic) || (kt_el < historic)){
                    historic = std::min(kn_el, kt_el);
                }

            }

            EvaluateDeltaDisplacement(DeltDisp, RelVel, LocalCoordSystem, OldLocalCoordSystem, other_to_me_vect, vel, delta_displ, neighbour_iterator, distance);

            DisplacementDueToRotation(DeltDisp, OldLocalCoordSystem, other_radius, dt, ang_vel, neighbour_iterator);

            double LocalDeltDisp[3] = {0.0};
            double LocalElasticContactForce[3]  = {0.0}; // 0: first tangential, // 1: second tangential, // 2: normal force
            double GlobalElasticContactForce[3] = {0.0};

            GlobalElasticContactForce[0] = mOldNeighbourContactForces[i_neighbour_count][0];  
            GlobalElasticContactForce[1] = mOldNeighbourContactForces[i_neighbour_count][1];
            GlobalElasticContactForce[2] = mOldNeighbourContactForces[i_neighbour_count][2];
            
            GeometryFunctions::VectorGlobal2Local(OldLocalCoordSystem, GlobalElasticContactForce, LocalElasticContactForce); 

            GeometryFunctions::VectorGlobal2Local(OldLocalCoordSystem, DeltDisp, LocalDeltDisp);
      

            /* Translational Forces */

            if  (indentation > 0.0 || (mNeighbourFailureId[i_neighbour_count] == 0) )
            {                                                                                                
              
                //Normal Forces
                
                    if (mapping_new_cont!=-1)
                    {
               
                      PlasticityAndDamage1D(LocalElasticContactForce, kn_el, equiv_young, indentation, calculation_area,radius_sum_i, failure_criterion_state, acumulated_damage, i_neighbour_count,mapping_new_cont, mapping_new_ini, rCurrentProcessInfo );
                                        
                    }
                    

            } //if compression cohesive contact
               
     
            //Tangential. With degradation:
            
            double degradation = 1.0;

            if(mapping_new_cont!= -1)
            {
             
              if(indentation >= 0.0 ) //COMPRESSION
              {
              
                degradation = mHistory[mapping_new_cont][3];
               
              }
              else
              {
               
                degradation = (1.0 -  mHistory[mapping_new_cont][2]);
               
              }
              
            }
                 
            LocalElasticContactForce[0] += - degradation*kt_el * LocalDeltDisp[0];  // 0: first tangential
            LocalElasticContactForce[1] += - degradation*kt_el * LocalDeltDisp[1];  // 1: second tangential  
              
                                         
            double ShearForceNow = sqrt(LocalElasticContactForce[0] * LocalElasticContactForce[0]
                                    +   LocalElasticContactForce[1] * LocalElasticContactForce[1]); 
                  
                /* Evaluating Failure for the continuum contacts */
          
                if(mNeighbourFailureId[i_neighbour_count] == 0)
                {                
                  /*
                  mNeighbourFailureId[i_neighbour_count] = 2; //shear in compression
                  mNeighbourFailureId[i_neighbour_count] = 3;  //shear failure tension
                  mNeighbourFailureId[i_neighbour_count] = 4; //tension failure
                  mNeighbourFailureId[i_neighbour_count] = 12; //both shear and tension
                  */
                  EvaluateFailureCriteria(LocalElasticContactForce,ShearForceNow,calculation_area,i_neighbour_count,contact_sigma,contact_tau, failure_criterion_state, sliding, mapping_new_ini);
                }
                   
                // if(*mpActivateSearch == 0)
                if(rCurrentProcessInfo[ACTIVATE_SEARCH] == 0)
                {
                    if(mNeighbourFailureId[i_neighbour_count]!=0)
                    {
                        rCurrentProcessInfo[ACTIVATE_SEARCH_VECTOR][OpenMPUtils::ThisThread()]=1;
                    }
                }
               
                // Tangential Friction for broken bonds //dempack and kdem do the same.
                
                if ( mNeighbourFailureId[i_neighbour_count] != 0 )  
                {
                   double Frictional_ShearForceMax = equiv_tg_of_fri_ang * LocalElasticContactForce[2];
                
                if (Frictional_ShearForceMax < 0.0)
                {
                  Frictional_ShearForceMax = 0.0;
                  
                }
                                    
                  failure_criterion_state = 1.0;
                                                        
                  if( (ShearForceNow >  Frictional_ShearForceMax) && (ShearForceNow != 0.0) ) 
                  {
                      LocalElasticContactForce[0] = (Frictional_ShearForceMax / ShearForceNow) * LocalElasticContactForce[0];
                      LocalElasticContactForce[1] = (Frictional_ShearForceMax / ShearForceNow )* LocalElasticContactForce[1];
                      sliding = true;
                                  
                  }
                  
                }

                   
               /* Viscodamping (applyied locally)*/ 
               
                 //DAMPING:

              double ViscoDampingLocalContactForce[3] = {0.0};
              
              if  (indentation > 0.0 || (mNeighbourFailureId[i_neighbour_count] == 0) )
              {  
                
                double LocalRelVel[3] = {0.0};
                
                GeometryFunctions::VectorGlobal2Local(OldLocalCoordSystem, RelVel, LocalRelVel);
                

                  equiv_visco_damp_coeff_normal = mDempack_damping*2.0*sqrt(kn_el/(mSqrtOfRealMass * mSqrtOfRealMass + other_sqrt_of_mass * other_sqrt_of_mass))*equiv_mass;   // := 2d0* sqrt ( kn_el*(m1*m2)/(m1+m2) )

                  equiv_visco_damp_coeff_tangential = equiv_visco_damp_coeff_normal * aux_norm_to_tang;
                
               
                CalculateViscoDamping(LocalRelVel,ViscoDampingLocalContactForce,indentation,equiv_visco_damp_coeff_normal,equiv_visco_damp_coeff_tangential,sliding);
                  
              }
                
            // Transforming to global forces and adding up
            double LocalContactForce[3] =                 {0.0};
            double ViscoDampingGlobalContactForce[3] =    {0.0}; 
            double GlobalContactForce[3] =                {0.0};
            
              
            AddUpForcesAndProject(LocalCoordSystem, LocalContactForce,LocalElasticContactForce,GlobalContactForce,
                                  GlobalElasticContactForce,ViscoDampingLocalContactForce,ViscoDampingGlobalContactForce,rContactForce,rElasticForce,
                                  i_neighbour_count);
            
     
            //AddPoissonContribution(LocalCoordSystem, GlobalContactForce, GlobalElasticContactForce, ViscoDampingGlobalContactForce, rContactForce, damp_forces); //MSIMSI 10

            //if(mRotationOption){
            if (this->Is(DEMFlags::HAS_ROTATION) ){
              ComputeMoments(LocalElasticContactForce,GlobalElasticContactForce,rInitialRotaMoment,LocalCoordSystem,other_radius,rContactMoment,neighbour_iterator);
            }

            if(mContactMeshOption==1 && (mapping_new_cont !=-1) && this->Id() < neighbour_iterator_id) 
            {

              CalculateOnContactElements( neighbour_iterator_id ,i_neighbour_count, mapping_new_cont, LocalElasticContactForce, contact_sigma, contact_tau, failure_criterion_state, acumulated_damage, time_steps);
              
            }

            if(mStressStrainOption)
            {
              
        
            StressTensorOperations(mStressTensor,GlobalElasticContactForce,
                                   other_to_me_vect,distance,radius_sum,calculation_area,
                                   neighbour_iterator,rCurrentProcessInfo,rRepresentative_Volume);
            


            }

            //i_neighbour_count++;

        }//for each neighbour
        
       
        KRATOS_CATCH("")         

}





    void PlasticityAndDamage(double LocalElasticContactForce[3], double kn_el, double equiv_young, double indentation, double      calculation_area, double radius_sum_i, double& failure_criterion_state, double& acumulated_damage, int i_neighbour_count, int mapping_new_cont, int mapping_new_ini, int time_steps){
       

        
      //VARIABLES
    
      // a guardar:
  
      double kn_b = kn_el / mN1;
      double kn_c = kn_el / mN2;
      double kn_d = kn_el / mN3;
      double kp_el = mYoungPlastic/equiv_young * kn_el;
      double Yields_el = mPlasticityLimit * calculation_area;
        
      double Ncstr1_el = mC1 * calculation_area;
      double Ncstr2_el = mC2 * calculation_area;
      double Ncstr3_el = mC3 * calculation_area;
      double Ntstr_el  = mTensionLimit * calculation_area;
      double u_max = mHistory[mapping_new_cont][0];
      
      double& fn = LocalElasticContactForce[2]; //[2] means 'normal' contact force
                
      
      if( indentation >= 0.0 ) //COMPRESSION
      {
        
          fn = kn_el * indentation;
          
          double u_ela1 = Ncstr1_el/kn_el;;
          double u_ela2 = u_ela1 + (Ncstr2_el-Ncstr1_el)/(kn_b);
          double u_ela3 = u_ela2 + (Ncstr3_el-Ncstr2_el)/(kn_c);

          if ( ( indentation > u_max ) || ( time_steps <= 1) )//maximum historical intentation OR first step  MSIMSI 0
            
          {

            mHistory[mapping_new_cont][0]  = indentation;             // Guarda el treshold del màxim desplaçament
            
            
            if (indentation > u_ela3) //4rt tram
            {
              
              fn = Ncstr3_el + ( indentation - u_ela3 )*kn_d;
              mHistory[mapping_new_cont][3] = kn_d/kn_el;
              
            }
            else if (indentation > u_ela2) //3r tram
            {
              
              fn = Ncstr2_el + ( indentation - u_ela2 )*kn_c;
              mHistory[mapping_new_cont][3] = kn_c/kn_el;
              
            }
            else
            {    
              if( indentation > u_ela1) //2n tram
              {
                fn = Ncstr1_el + (indentation - u_ela1)*kn_b;
                mHistory[mapping_new_cont][3] = kn_b/kn_el;
              
              }
              
            }
          
          mHistory[mapping_new_cont][1] = fn; //actualitzar la força màxima a compressió.
          
          }
          
          else //Per sota del màxim.
          {

              if(mHistory[mapping_new_cont][1] > 0.0)  //Màxim en compressió. 
              {
 
                  double u_plas;        //MSIMSI 2 akesta operació de saber quant val la u_plastica es fa cada pas de temps i en realitat es fixe sempre.

                  if(Yields_el <= Ncstr1_el) //si el punt de plastificació està en la primera rama elastica.
                  {
                      u_plas = Yields_el/kn_el;

                  }
                  else
                  {  
                    if(Yields_el <= Ncstr2_el) //si està en la segona...
                    {
                        u_plas = u_ela1 + (Yields_el-Ncstr1_el)/(kn_b);

                    }
                     else if(Yields_el <= Ncstr3_el) //si està en la tercera...
                    {
                        u_plas = u_ela2 + (Yields_el-Ncstr2_el)/(kn_c);

                    }
                    
                    else //en la quarta
                    {

                      u_plas = u_ela3 + (Yields_el-Ncstr3_el)/(kn_d);
                    }
                    
                  }

                  
                  if ( u_plas < u_max ) //si nosaltres estem per sota del maxim pero ja estem plastificant 
                  {
                    fn = mHistory[mapping_new_cont][1] - kp_el*(u_max - indentation); // Esta en zona de descarga plastica (pot estar en carga/descarga)
                    mHistory[mapping_new_cont][3] = kp_el/kn_el;
                    
                    
                  }
                  else                                   // Esta en zona descarga elastica, ens despreocupem de la plasticitat
                  {

                    if ( indentation > u_ela3)  //en la 4a ramma
                    {
                      fn = Ncstr3_el + (indentation - u_ela3)*kn_d;
                      
                    }
                    
                    else if ( indentation > u_ela2)  //en la 3a ramma
                    {
                      fn = Ncstr2_el + (indentation - u_ela2)*kn_c;
                      
                    }
                    
                    else
                    {
                      if(indentation > u_ela1)  //en la 2a rama
                      {
                        fn = Ncstr1_el + (indentation-u_ela1)*kn_b;
                      }
            
                    }
            
                  }

              } //si tenim precàrrega en compressió.
              
          }//Per sota del màxim.


      } //Compression
              

      else //tension
      {
        fn = kn_el * indentation; 
        
        double u1 = Ntstr_el / kn_el;

        double u2 = u1*(1+ mDamageMaxDisplacementFactor);
 
          if(fabs(indentation) > u2)                  // FULL DAMAGE 
          {
            mNeighbourFailureId[i_neighbour_count] = 4; //tension failure
            mIniNeighbourFailureId[ mapping_new_ini ] = 4;
            acumulated_damage = 1.0;
            fn = 0.0;
          }
          else
          {
            if (fabs(indentation) > u1)  
            {
              double u_frac = (fabs(indentation) - u1)/(u2 - u1);
              //failure_criterion_state = fabs(indentation)/u2;
              acumulated_damage = u_frac;
                            
              if (u_frac > mHistory[mapping_new_cont][2])  
              {
                mHistory[mapping_new_cont][2] = u_frac;
              }
              
            }
            
            fn = indentation * kn_el*(1.0 -  mHistory[mapping_new_cont][2]);  // normal adhesive force (gap +)
            
          }
      
        }//Tension

    }

} /* namespace Kratos.*/
