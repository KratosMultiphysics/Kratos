
//   Project Name:        Kratos
//   Last Modified by:    $Author: M.Santasusana $
//   Date:                $Date: 2006-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


// System includes
#include <string>
#include <iostream>

// External includes


// Project includes
#include "includes/define.h"
#include "spheric_particle.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "DEM_application.h"
#include "utilities/openmp_utils.h"

//TIMER....................
#include "utilities/timer.h"

#define CUSTOMTIMER 0  // ACTIVATES AND DISABLES ::TIMER:::::

#ifdef CUSTOMTIMER
#define KRATOS_TIMER_START(t) Timer::Start(t);
#define KRATOS_TIMER_STOP(t) Timer::Stop(t);
#else
#define KRATOS_TIMER_START(t)
#define KRATOS_TIMER_STOP(t)
#endif

//......................
//std::cout<<print...<<std::endl;

namespace Kratos
{

      SphericContinuumParticle::SphericContinuumParticle() : SphericParticle(){}

      SphericContinuumParticle::SphericContinuumParticle( IndexType NewId, GeometryType::Pointer pGeometry) : SphericParticle(NewId, pGeometry){}

      SphericContinuumParticle::SphericContinuumParticle( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
      : SphericParticle(NewId, pGeometry, pProperties){}

      SphericContinuumParticle::SphericContinuumParticle(IndexType NewId, NodesArrayType const& ThisNodes)
      : SphericParticle(NewId, ThisNodes){}

      Element::Pointer SphericContinuumParticle::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
      {
         return SphericParticle::Pointer(new SphericContinuumParticle(NewId, GetGeometry().Create( ThisNodes ), pProperties) );
      }

      /// Destructor.
      SphericContinuumParticle::~SphericContinuumParticle(){}

  
      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void SphericContinuumParticle::SetInitialContacts( const ProcessInfo& r_process_info  ) //vull ficar que sigui zero si no son veins cohesius.
      {   
                
          /*
            * 
            * ELEMENT_NEIGHBOURS / NEIGHBOURS_IDS: the ones important for calculating forces.
            * INI_NEIGHBOURS_IDS: the ones to be treated specially due to initial delta or continuum case.
            * INI_CONTINUUM_NEIGHBOURS_IDS: only the ones that are continuum at 0 step and we should treat the possible detachment.
            * 
            * These 3 classes do NOT coincide at t=0!

          */

            // DEFINING THE REFERENCES TO THE MAIN PARAMETERS

            ParticleWeakVectorType& mrNeighbours                  = this->GetValue(NEIGHBOUR_ELEMENTS);            
            ParticleWeakVectorType& r_continuum_ini_neighbours    = this->GetValue(CONTINUUM_INI_NEIGHBOUR_ELEMENTS);
            vector<int>& r_continuum_ini_neighbours_ids          = this->GetValue(CONTINUUM_INI_NEIGHBOURS_IDS);
            
            r_continuum_ini_neighbours.clear();
            r_continuum_ini_neighbours_ids.clear();
                        
            size_t ini_size = 0;
            size_t continuum_ini_size =0;

                       
            //default
            mFailureId=1;

            size_t cont_ini_mapping_index= 0;
            
            //SAVING THE INICIAL NEIGHBOURS, THE DELTAS AND THE FAILURE ID
                              
            for(ParticleWeakIteratorType_ptr ineighbour = mrNeighbours.ptr_begin();  //loop over the temp neighbours and store into a initial_neighbours vector.
            ineighbour != mrNeighbours.ptr_end(); ineighbour++)
            {
               
              
                array_1d<double,3> other_to_me_vect = this->GetGeometry()(0)->Coordinates() - ((*ineighbour).lock())->GetGeometry()(0)->Coordinates();
                double distance                     = sqrt(other_to_me_vect[0] * other_to_me_vect[0] +
                                                      other_to_me_vect[1] * other_to_me_vect[1] +
                                                      other_to_me_vect[2] * other_to_me_vect[2]);

                double radius_sum                   = mRadius + ((*ineighbour).lock())->GetGeometry()(0)->FastGetSolutionStepValue(RADIUS);
                double initial_delta                = radius_sum - distance;

                int r_other_continuum_group         = ((*ineighbour).lock())->GetGeometry()(0)->FastGetSolutionStepValue(PARTICLE_CONTINUUM);
   
                if( ( (r_other_continuum_group == mContinuumGroup) && (mContinuumGroup != 0) ) || ( fabs(initial_delta)>1.0e-6 ) ) // which would be the appropiate tolerance?
                
                //for the ones that need to be stored... #5
                {
                    ini_size++;
                    
                    mIniNeighbourIds.resize(ini_size); 
                    mIniNeighbourToIniContinuum.resize(ini_size);
                    mIniNeighbourDelta.resize(ini_size);
                    mIniNeighbourFailureId.resize(ini_size);           
                    mMapping_New_Ini.resize(ini_size); 
        
                    mIniNeighbourIds[ini_size - 1]          = ((*ineighbour).lock())->Id();
                    mIniNeighbourDelta[ini_size - 1]        = 0.0;
                    mIniNeighbourFailureId[ini_size - 1]    = 1;
                    mMapping_New_Ini[ini_size - 1]          = ini_size - 1;
                    
                    mIniNeighbourToIniContinuum[ini_size-1] = -1; //-1 is initial but not continuum.             
                    
                    
                    if (mDeltaOption == true)
                    {            
                        mIniNeighbourDelta[ini_size - 1]    = initial_delta;

                                                
                    }

                    if (mContinuumSimulationOption == true)
                    {
                        
                        if ( (r_other_continuum_group == mContinuumGroup) && (mContinuumGroup != 0) )
                        {
                            
                            mIniNeighbourToIniContinuum[ini_size-1] = cont_ini_mapping_index;
                        
                            mIniNeighbourFailureId[ini_size - 1]=0;
                           
                            mFailureId=0; // if a cohesive contact exist, the FailureId becomes 0. 
                            
                            continuum_ini_size++;
                            cont_ini_mapping_index++;
                                
                            r_continuum_ini_neighbours.push_back(*ineighbour);
                              
                            r_continuum_ini_neighbours_ids.resize(continuum_ini_size);
                            
                            mMapping_New_Cont.resize(continuum_ini_size);
                            mMapping_New_Cont[continuum_ini_size - 1] = -1;
                            
                            mHistory.resize(continuum_ini_size);
                            
                            mHistory[continuum_ini_size - 1][0] = 0.0; //maximum indentation reached
                            mHistory[continuum_ini_size - 1][1] = 0.0; //maximum force reached
                            mHistory[continuum_ini_size - 1][2] = 0.0;
                            
                            r_continuum_ini_neighbours_ids[continuum_ini_size - 1] = ((*ineighbour).lock())->Id();
                                                                                    
                           if(mContactMeshOption)
                            {

                                (this->GetGeometry()(0))->GetValue(NODE_TO_NEIGH_ELEMENT_POINTER).resize(continuum_ini_size);
                    
                            } //if(mContactMeshOption) 
                            
                        }//if ( (r_other_continuum_group == mContinuumGroup) && (mContinuumGroup != 0) )

                    }//for mContinuumSimulationOption      

                } // FOR THE CASES THAT NEED STORING INITIAL NEIGHBOURS
              
            } //end for: ParticleWeakIteratorType ineighbour
            

            if (mContinuumSimulationOption == true)
            {
              
                if(mDimension == 3)
                {
                  ContactAreaWeighting3D(r_process_info);
                }
                else if (mDimension ==2)
                {
                  
                  ContactAreaWeighting2D(r_process_info);
                  
                }  
    
            } 
            
         
        }//SetInitialContacts
        
        void SphericContinuumParticle::NeighNeighMapping( ProcessInfo& rCurrentProcessInfo  ) 
        {
          
//         ParticleWeakVectorType& r_continuum_ini_neighbours    = this->GetValue(CONTINUUM_INI_NEIGHBOUR_ELEMENTS);
// 
//         int my_id = this->Id();
//         
//         vector<int> my_index_from_my_neigh;
//         
//         my_index_from_my_neigh.resize(r_continuum_ini_neighbours.size());
//         
//         KRATOS_WATCH("   ")
//         KRATOS_WATCH(this->Id())
//                  
//         int cont_ini_neighbour_counter = 0;
//                 
//         for(ParticleWeakIteratorType cont_ini_neighbour_iterator = r_continuum_ini_neighbours.begin();
//             cont_ini_neighbour_iterator != r_continuum_ini_neighbours.end(); cont_ini_neighbour_iterator++)
//         {
//              
//               KRATOS_WATCH(cont_ini_neighbour_iterator->Id())
//               
//           
//               ParticleWeakVectorType& r_continuum_ini_neighbours_of_my_continuum_ini_neighbour    = cont_ini_neighbour_iterator->GetValue(CONTINUUM_INI_NEIGHBOUR_ELEMENTS);
// 
//               //KRATOS_WATCH(r_continuum_ini_neighbours_of_my_continuum_ini_neighbour.size())
//               /*
//               int index_me_from_my_neigh = 0;
//               
//                 for(ParticleWeakIteratorType cont_ini_neighbour_of_my_cont_ini_neigh_iterator = r_continuum_ini_neighbours_of_my_continuum_ini_neighbour.begin();
//                   cont_ini_neighbour_of_my_cont_ini_neigh_iterator != r_continuum_ini_neighbours_of_my_continuum_ini_neighbour.end(); cont_ini_neighbour_of_my_cont_ini_neigh_iterator++)
//               {
//           
//                 if( cont_ini_neighbour_of_my_cont_ini_neigh_iterator->Id() == my_id)
//                 {
//                                   
//                 my_index_from_my_neigh[cont_ini_neighbour_counter] = index_me_from_my_neigh;
//                                 
//                 break;  
//                 }
//                 
//                 index_me_from_my_neigh++;
//                
//               } //neighbours of my neighbour
//              */ 
//         } //my neighbours      
//         
//         
//         //checks:
//         //loop veins ini cont
//           //vei->trencat?
//           
//           //vei->vei(jo)->trencat?
//         
//         
//         
//         KRATOS_WATCH(my_index_from_my_neigh)
//         
        
      }// 
      

      
      void SphericContinuumParticle::ContactAreaWeighting3D(const ProcessInfo& r_process_info) //MISMI 10: POOYAN this could be done by calculating on the bars. not looking at the neighbous of my neighbours.
      { 

          double alpha = 1.0;
          double external_sphere_area = 4*M_PI*mRadius*mRadius;  
          
          double total_equiv_area = 0.0;

          ParticleWeakVectorType r_continuum_ini_neighbours    = this->GetValue(CONTINUUM_INI_NEIGHBOUR_ELEMENTS);
          int cont_ini_neighbours_size                         = r_continuum_ini_neighbours.size();
          
          mcont_ini_neigh_area.resize(cont_ini_neighbours_size);
          
          //computing the total equivalent area
          
          size_t index = 0;
          
          for(ParticleWeakIteratorType ini_cont_neighbour_iterator = r_continuum_ini_neighbours.begin();     // MSIMSI 99:Could this loop be done during the bar creation in the strategy and so avoid another repetition?
              ini_cont_neighbour_iterator != r_continuum_ini_neighbours.end(); ini_cont_neighbour_iterator++)
          {   
              double other_radius     = ini_cont_neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(RADIUS);
              double equiv_radius     = 2*mRadius * other_radius / (mRadius + other_radius);        
              double equiv_area       = (0.25)*M_PI * equiv_radius * equiv_radius; //we now take 1/2 of the efective mRadius.
              total_equiv_area       += equiv_area;
          
              mcont_ini_neigh_area[index] = equiv_area; //*
              index++; //*
              
              
          } //for every neighbour
       
            if (cont_ini_neighbours_size >= 4) //more than 3 neigbours. 
            {
                if(!*mSkinSphere)
                {
                
                  AuxiliaryFunctions::CalculateAlphaFactor3D(cont_ini_neighbours_size, external_sphere_area, total_equiv_area, alpha); 
                  
                  size_t not_skin_index = 0;
              
                  for(ParticleWeakIteratorType ini_cont_neighbour_iterator = r_continuum_ini_neighbours.begin();
                      ini_cont_neighbour_iterator != r_continuum_ini_neighbours.end(); ini_cont_neighbour_iterator++)
                      
                      {      
                          mcont_ini_neigh_area[not_skin_index] = alpha*mcont_ini_neigh_area[not_skin_index];
                          not_skin_index++;  
                          
                      } //for every neighbour

                }
                
                else //skin sphere 
                {
                  
                  size_t skin_index = 0;
              
                  for(ParticleWeakIteratorType ini_cont_neighbour_iterator = r_continuum_ini_neighbours.begin();
                      ini_cont_neighbour_iterator != r_continuum_ini_neighbours.end(); ini_cont_neighbour_iterator++)
                      
                      {  
                  
                      alpha  = 1.00*(1.40727)*(external_sphere_area/total_equiv_area)*((double(cont_ini_neighbours_size))/11);
                      mcont_ini_neigh_area[skin_index] = alpha*mcont_ini_neigh_area[skin_index];
                  
                      skin_index++;
                      }
                    
                }//skin particles.
              
            }//if more than 3 neighbours
   
      } //Contact Area Weighting           
                              
   
      /**
       * Calculates all particle's ball-to-ball forces based on its neighbours
       * @param rContactForce
       * @param rContactMoment
       * @param rCurrentProcessInfo
       **/
                                                                                                                      
      void SphericContinuumParticle::ComputeBallToBallContactForce(array_1d<double, 3>& rContactForce, array_1d<double, 3>& rContactMoment, 
                                                                   array_1d<double, 3>& rElasticForce, array_1d<double, 3>& rInitialRotaMoment, 
                                                                   ProcessInfo& rCurrentProcessInfo)
      {
                 
        KRATOS_TRY

        ParticleWeakVectorType& mrNeighbours         = this->GetValue(NEIGHBOUR_ELEMENTS); 

        double dt = rCurrentProcessInfo[DELTA_TIME];
        double dt_i = 1 / dt; 
                
        /* Initializations */
                          
        const array_1d<double, 3>& vel         = this->GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY);
        const array_1d<double, 3>& delta_displ = this->GetGeometry()(0)->FastGetSolutionStepValue(DELTA_DISPLACEMENT);
        const array_1d<double, 3>& ang_vel     = this->GetGeometry()(0)->FastGetSolutionStepValue(ANGULAR_VELOCITY);
        double RotaAcc[3]                      = {0.0};
        double InitialRotaMoment[3]            = {0.0};

        if (mRotationOption){
            RotaAcc[0]                         = ang_vel[0] * dt_i;
            RotaAcc[1]                         = ang_vel[1] * dt_i;
            RotaAcc[2]                         = ang_vel[2] * dt_i;

            InitialRotaMoment[0] = RotaAcc[0] * mMomentOfInertia;       
            InitialRotaMoment[1] = RotaAcc[1] * mMomentOfInertia;
            InitialRotaMoment[2] = RotaAcc[2] * mMomentOfInertia;

        }        

        size_t i_neighbour_count = 0;

        for(ParticleWeakIteratorType neighbour_iterator = mrNeighbours.begin();
            neighbour_iterator != mrNeighbours.end(); neighbour_iterator++)
        {

            array_1d<double,3> other_to_me_vect   = this->GetGeometry()(0)->Coordinates() - neighbour_iterator->GetGeometry()(0)->Coordinates();
            const double &other_radius                  = neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(RADIUS);
            const double &other_sqrt_of_mass            = neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(SQRT_OF_MASS);
            const double &other_ln_of_restit_coeff  = neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(LN_OF_RESTITUTION_COEFF);

            double distance                       = sqrt(other_to_me_vect[0] * other_to_me_vect[0] +
                                                          other_to_me_vect[1] * other_to_me_vect[1] +
                                                          other_to_me_vect[2] * other_to_me_vect[2]);
            double radius_sum                     = mRadius + other_radius;
            double radius_sum_i                   = 1 / radius_sum;
            double equiv_radius                   = 2 * mRadius * other_radius * radius_sum_i; 
            
            double initial_delta                  = mNeighbourDelta[i_neighbour_count]; //*
            double initial_dist                   = (radius_sum - initial_delta);
            double indentation                    = initial_dist - distance;   //#1
            double equiv_area                     = 0.25*M_PI * equiv_radius * equiv_radius; //#2 
            double calculation_area               = equiv_area;
            double equiv_mass                     = mSqrtOfRealMass * other_sqrt_of_mass;

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

            double NormalDir[3]                   = {0.0};
            double OldNormalDir[3]                = {0.0};

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
          

            if (mUniformMaterialOption){
                equiv_radius                      = mRadius;
                equiv_young                       = mYoung;
                equiv_poisson                     = mPoisson;
                equiv_ln_of_restit_coeff          = mLnOfRestitCoeff;
                equiv_tg_of_fri_ang               = mTgOfFrictionAngle;
            }

            else {
                // Getting neighbour properties
                double &other_young               = neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(YOUNG_MODULUS);
                double &other_poisson             = neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(POISSON_RATIO);
                double &other_ln_of_restit_coeff  = neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(LN_OF_RESTITUTION_COEFF);
                double &other_tg_of_fri_angle     = neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(PARTICLE_FRICTION);

                equiv_young                       = 2 * mYoung * other_young / (mYoung + other_young);
                equiv_poisson                     = 2 * mPoisson * other_poisson / (mPoisson + other_poisson);
                equiv_ln_of_restit_coeff          = 0.5 * (mLnOfRestitCoeff + other_ln_of_restit_coeff);
                equiv_tg_of_fri_ang               = 0.5 * (mTgOfFrictionAngle + other_tg_of_fri_angle);
            }

            // Globally defined parameters
                double aux_norm_to_tang = 0.0;

            if (mGlobalVariablesOption){
              
                kn_el                                = mGlobalKn;
                kt_el                                = mGlobalKt;
                aux_norm_to_tang                     = mGlobalAuxNormToTang;

            }

            else if(mDempack){
                
                double rad_squared = mRadius * other_radius;
                
                calculation_area = 4.0*M_PI*(rad_squared*rad_squared)*radius_sum_i*radius_sum_i;
                
                
                double equiv_shear = equiv_young/(2*(1+equiv_poisson));
                kn_el = equiv_young*calculation_area*radius_sum_i;
                kt_el = equiv_shear*calculation_area*radius_sum_i;
              

            }
            
            else
            {

                if(mContinuumSimulationOption==1 && (mapping_new_ini !=-1))
                {
                  calculation_area = mcont_ini_neigh_area[mapping_new_ini];                            
                }
              
                kn_el              = mMagicFactor * equiv_young * calculation_area * radius_sum_i; //MSIMSI 1: initial gap? aki dividim nomes per suma de radis.
                kt_el              = mMagicFactorPoisson * kn_el/(2.0 + equiv_poisson + equiv_poisson);
                aux_norm_to_tang   = sqrt(kt_el / kn_el);


            }
          
            if (mCriticalTimeOption){
                double historic = rCurrentProcessInfo[HISTORICAL_MIN_K];

                if ((kn_el < historic) || (kt_el < historic)){
                    historic = fmin(kn_el, kt_el);
                }

            }

            //DAMPING:
            
            if(mDempack){
            
            equiv_visco_damp_coeff_normal     = mDempack_damping*2.0*sqrt(kn_el/(mRealMass+other_sqrt_of_mass*other_sqrt_of_mass))*equiv_mass;   // := 2d0* sqrt ( kn_el*(m1*m2)/(m1+m2) )
            equiv_visco_damp_coeff_tangential = equiv_visco_damp_coeff_normal; // dempack no l'utilitza...
            
            }
              
            else{ //KDEM
            
              if (mLnOfRestitCoeff > 0.0 || other_ln_of_restit_coeff > 0.0){
                  
                  equiv_visco_damp_coeff_normal     = 2 * sqrt(equiv_mass * kn_el);
                  equiv_visco_damp_coeff_tangential = equiv_visco_damp_coeff_normal * aux_norm_to_tang; 
              }

              else {
                  
                  equiv_visco_damp_coeff_normal     = - 2 * equiv_ln_of_restit_coeff * sqrt(equiv_mass * kn_el / (equiv_ln_of_restit_coeff * equiv_ln_of_restit_coeff + M_PI * M_PI));
                  equiv_visco_damp_coeff_tangential = equiv_visco_damp_coeff_normal * aux_norm_to_tang; 
              }
              
            }
            
            
            EvaluateDeltaDisplacement(DeltDisp, RelVel, NormalDir, OldNormalDir, LocalCoordSystem, OldLocalCoordSystem, other_to_me_vect, vel, delta_displ, neighbour_iterator);

            DisplacementDueToRotation(DeltDisp, OldNormalDir, OldLocalCoordSystem, other_radius, dt, ang_vel, neighbour_iterator);
            
            double LocalDeltDisp[3] = {0.0};
            double LocalElasticContactForce[3]  = {0.0}; // 0: first tangential, // 1: second tangential, // 2: normal force
            double GlobalElasticContactForce[3] = {0.0};
            double LocalRelVel[3] = {0.0};
            double ViscoDampingLocalContactForce[3]    = {0.0};

            GlobalElasticContactForce[0] = mOldNeighbourContactForces[i_neighbour_count][0];  
            GlobalElasticContactForce[1] = mOldNeighbourContactForces[i_neighbour_count][1];
            GlobalElasticContactForce[2] = mOldNeighbourContactForces[i_neighbour_count][2];
            
            GeometryFunctions::VectorGlobal2Local(OldLocalCoordSystem, GlobalElasticContactForce, LocalElasticContactForce); 
            //we recover this way the old local forces projected in the new coordinates in the way they were in the old ones; Now they will be increased if its the necessary
            GeometryFunctions::VectorGlobal2Local(OldLocalCoordSystem, DeltDisp, LocalDeltDisp);
            GeometryFunctions::VectorGlobal2Local(OldLocalCoordSystem, RelVel, LocalRelVel);
            
            /* Translational Forces */
            
            //Tangential Forces   #4


            LocalElasticContactForce[0] += - kt_el * LocalDeltDisp[0];  // 0: first tangential
            LocalElasticContactForce[1] += - kt_el * LocalDeltDisp[1];  // 1: second tangential
                              
            double ShearForceNow = sqrt(LocalElasticContactForce[0] * LocalElasticContactForce[0]
                                    +   LocalElasticContactForce[1] * LocalElasticContactForce[1]); 
        
            double Frictional_ShearForceMax = equiv_tg_of_fri_ang * LocalElasticContactForce[2];
            
            if (Frictional_ShearForceMax < 0.0)
            {
              Frictional_ShearForceMax = 0.0;
              
            }
            
               
            if  (indentation > 0.0 || (mNeighbourFailureId[i_neighbour_count] == 0) )//*  //#3
            {                                                                                                
              
                //Normal Forces
                
                if(mElasticityType < 2)
                {
                      NormalForceCalculation(LocalElasticContactForce,kn_el,indentation);

                }
                
                else if(mElasticityType==2){
                
               
                    if (mapping_new_cont!=-1)
                    {
                    
                  
                      
                      PlasticityAndDamage1D(LocalElasticContactForce, kn_el, indentation, calculation_area,radius_sum_i, failure_criterion_state, acumulated_damage, i_neighbour_count,mapping_new_cont, mapping_new_ini );
                    
                      
                    }
                    
                    else
                    
                    {
                      NormalForceCalculation(LocalElasticContactForce,kn_el,indentation);
                      
                    }
                    
                }//plasticity and damage for the initial continuum contacts only.

            } //if compression or cohesive contact
               

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
                   
                if(*mpActivateSearch == 0)
                {
                    if(mNeighbourFailureId[i_neighbour_count]!=0)
                    {
                        rCurrentProcessInfo[ACTIVATE_SEARCH_VECTOR][OpenMPUtils::ThisThread()]=1;
                    }
                  
                }
                
                       /*   MSIMSI 10 Activar la busqueda quan un peti. 
          if(mNeighbourFailureId[i_neighbour_count] != 0 && rCurrentProcessInfo[ACTIVATE_SEARCH]==0)
          {
              rCurrentProcessInfo.SetValue(ACTIVATE_SEARCH, 1);

              KRATOS_WATCH(" ")
              KRATOS_WATCH("-------->From now on, searching neighbours, some contacs have failed<-------")
              KRATOS_WATCH("Time step:")
              KRATOS_WATCH(*mpTimeStep)
              KRATOS_WATCH("Particle_1")
              KRATOS_WATCH(this->Id())
              KRATOS_WATCH("Particle_2")
              KRATOS_WATCH(neighbour_iterator->Id())      
              KRATOS_WATCH(" ")
          

           }
      */
                
                
                /* Tangential Friction for broken bonds */  //dempack and kdem do the same.
                
                if ( mNeighbourFailureId[i_neighbour_count] != 0 ) //*   //degut als canvis de DEMPACK hi ha hagut una modificació, ara despres de trencar es fa akest maping de maxima tangencial que és correcte!
                {
                  failure_criterion_state = 1.0;
                                                        
                  if( (ShearForceNow >  Frictional_ShearForceMax) && (ShearForceNow != 0.0) ) 
                  {
                      LocalElasticContactForce[0] = (Frictional_ShearForceMax / ShearForceNow) * LocalElasticContactForce[0];
                      LocalElasticContactForce[1] = (Frictional_ShearForceMax / ShearForceNow )* LocalElasticContactForce[1];
                      sliding = true;
                                  
                  }
                  
                }

                   
               /* Viscodamping (applyied locally)*/ 

                  if  (indentation > 0.0 || (mNeighbourFailureId[i_neighbour_count] == 0) )//*  //#3
                  {  
                    
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

            
            ComputeMoments(LocalElasticContactForce,GlobalElasticContactForce,InitialRotaMoment,LocalCoordSystem,other_radius,rContactMoment,neighbour_iterator);
            
            //StressTensorOperations(mStressTensor,GlobalElasticContactForce,other_to_me_vect,distance,radius_sum,calculation_area,neighbour_iterator,rCurrentProcessInfo); //MSISI 10
  
            if(mContactMeshOption==1 && (mapping_new_cont !=-1)) 
            {

              CalculateOnContactElements( neighbour_iterator_id ,i_neighbour_count, mapping_new_cont, LocalElasticContactForce, contact_sigma, contact_tau, failure_criterion_state, acumulated_damage);

              
            }

            i_neighbour_count++;

        }//for each neighbour
        
        rInitialRotaMoment [0] = InitialRotaMoment [0];
        rInitialRotaMoment [1] = InitialRotaMoment [1];
        rInitialRotaMoment [2] = InitialRotaMoment [2];


        //ComputeStressStrain(mStressTensor, rCurrentProcessInfo);  //MSIMSI 10
        
        KRATOS_CATCH("")         
            
        
      }//ComputeBallToBallContactForce


      void SphericContinuumParticle::NonlinearNormalForceCalculation(double LocalElasticContactForce[3], double kn1, double kn2, double distance, double max_dist, double initial_dist)
      {

          LocalElasticContactForce[2] = kn1 * (initial_dist - max_dist) + kn2 * (max_dist - distance);
      }

    
      void SphericContinuumParticle::ApplyLocalMomentsDamping(const ProcessInfo& rCurrentProcessInfo )
      {

          KRATOS_TRY
      
          array_1d<double, 3 > & RotaMoment       = this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT);
          double RotaDampRatio                    = this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_ROTATION_DAMP_RATIO);

          // LOCAL DAMPING OPTION FOR THE UNBALANCED FORCES (IN GLOBAL CORDINATES).
        
          for (int iDof = 0; iDof < 3; iDof++)
          {
              if (this->GetGeometry()(0)->FastGetSolutionStepValue(ANGULAR_VELOCITY)[iDof] > 0.0)
              {
                  RotaMoment[iDof] = RotaMoment[iDof] - RotaDampRatio * fabs(RotaMoment[iDof]);
              }
              else
              {
                  RotaMoment[iDof] = RotaMoment[iDof] + RotaDampRatio * fabs(RotaMoment[iDof]);
              }
          }

          KRATOS_CATCH("")

      } //ApplyLocalMomentsDamping

      void SphericContinuumParticle::CharacteristicParticleFailureId(const ProcessInfo& rCurrentProcessInfo ) //MSIMSI 9: aixo pot ser %veins trencats...
      {  

          KRATOS_TRY

          // SETTING THE CARACHTERISTIC FAILURE TYPE OF THE PARTICLE (THE DOMINANT ONE) FOR VISUALITZATION.

            /* ContactFailureId = 2 (partially detached) it's never considerer in a contact, but it can be a desired output on the visualitzation.
             *
             * If the particle has no neighbour, but it is not from ContinuumGroup = 0 its becouse it has lost all the neighbours,
             * we are not interested in its global failure type but in its neighbours failure type. See: //// DETECTING A COHESIVE PARTICLE THAT HAS BEEN COMPLETELY DETACHED. in AfterForceCalculation.
             *
             * If one of the failure types completely dominates over the others this is the FailureId choosen for the particle
             *
             * If its partially detached by some type(s) of failure but the dominant contact state is "still atached", so failureId =0, we chose FailureId = 2.
             *
             * If the particle is not simulating the continuum, the failure type is 1 (default set from the particle whose ContinuumGroup is 0); It won't be printed as output.
             *
             */

            /*
             *   mFailureId values:
             *      0 := all neighbours attached
             *      1 := General detachment (all neighbours lost or no initial continuum case)
             *      2 := Partially detached and no dominance
             *      3 := tensile dominating on the contact detachment
             *      4 := shear dominating on the contact detachment
             */


          int tempType[5] = {0,0,0,0,0};

          if (mFailureId != 1)  // for mFailureId == 1 there's no failure to represent, the particle is not a continuum-simulating one or has been completelly detached already.
          {
              vector<int>& r_initial_neighbours_id          = this->GetValue(INI_NEIGHBOURS_IDS);

              int ini_neighbours_size = r_initial_neighbours_id.size();
     
              for(int index = 0; index < ini_neighbours_size; index++)    //loop over neighbours, just to run the index.
              {    
                  if( this->GetValue(PARTICLE_INITIAL_FAILURE_ID)[index] == 0)
                  {
                      tempType[0]++;
                  }
                  else if( this->GetValue(PARTICLE_INITIAL_FAILURE_ID)[index] == 1)
                  {
                      tempType[1]++;
                  }

                  // mContactFailureId == 2 intentionally skipped!

                  else if( this->GetValue(PARTICLE_INITIAL_FAILURE_ID)[index] == 3)
                  {
                      tempType[3]++;
                  }
                  else if( this->GetValue(PARTICLE_INITIAL_FAILURE_ID)[index] == 4)
                  {
                      tempType[4]++;
                  }
              }

              if ( tempType[0] == 0)  //no neighbour is attached
              {
                  mFailureId = 1;
              }   // no one neighbour is attached (but maybe still contacting).
              else if( (tempType[3] > tempType[4]) ) //some neighbour attached but failure 3 dominates over 4.
              {
                  mFailureId = 3;
              }
              else if( (tempType[4] > tempType[3]) ) // the same but 4 dominates over 3.
              {
                  mFailureId = 4;
              }
              else if ( (tempType[4] > 0) || (tempType[3] > 0) ) // no 3 neither 4 dominates but one of them may exist.
              {
                  mFailureId = 2;  // Partially detached / mix case.
              }
              else
              {
                  mFailureId = 0;  // last option: no one detached.
              }

          }// if (mFailureId != 1)

          KRATOS_CATCH("")

      } //CharacteristicParticleFailureId

     
      void SphericContinuumParticle::SymmetrizeTensor(const ProcessInfo& rCurrentProcessInfo)   //MSIMSI10
      {

        /*
        
        for (int i=0; i<3; i++)
        {
            for (int j=0; j<3; j++)
            {
                mSymmStressTensor[i][j] = 0.5*(mStressTensor[i][j]+mStressTensor[j][i]);
                
            }

        }

            //composició 
        
          ParticleWeakVectorType& mrNeighbours                = this->GetValue(NEIGHBOUR_ELEMENTS); //MSIMSI 10 aixo no sha de treure de un getvalue sino del membre.
        
            
            for(ParticleWeakIteratorType neighbour_iterator = mrNeighbours.begin();
                neighbour_iterator != mrNeighbours.end(); neighbour_iterator++)
            {            
                // GETTING NEIGHBOUR PROPERTIES

                //searching for the area
                double other_radius     = neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(RADIUS);
                double other_poisson    = neighbour_iterator->GetGeometry()[0].FastGetSolutionStepValue(POISSON_RATIO);
                double equiv_radius     = 2*mRadius * other_radius / (mRadius + other_radius);
                int size_ini_cont_neigh = this->GetValue(CONTINUUM_INI_NEIGHBOURS_IDS).size();
                double equiv_area       = (0.25)*M_PI * equiv_radius * equiv_radius; // 0.25 is becouse we take only the half of the equivalent mRadius, corresponding to the case of one sphere with mRadius Requivalent and other = mRadius 0.
                double equiv_poisson    = 2* mPoisson * other_poisson / (mPoisson + other_poisson);
                //double equiv_young      = 2 * young * other_young / (young + other_young);
                //bool is_continuum       = false;
                  double calculation_area = equiv_area;
                  
                  bool found = false;
                  
                  for (int index_area=0; index_area<size_ini_cont_neigh; index_area++)
                  {

                        if ( this->GetValue(CONTINUUM_INI_NEIGHBOURS_IDS)[index_area] == int(neighbour_iterator->Id()) ) 
                        {
                            
                          //MPI_CARLOS_ en MPI no pot funcionar lo de symmetrize de moment

                                Element::Pointer lock_p_weak = (this->GetGeometry()[0].GetValue(NODE_TO_NEIGH_ELEMENT_POINTER)(index_area)).lock();

                                calculation_area = lock_p_weak->GetValue(MEAN_CONTACT_AREA);

                                  array_1d<double,3> other_to_me_vect = this->GetGeometry()(0)->Coordinates() - neighbour_iterator->GetGeometry()(0)->Coordinates();
                
                                double NormalDir[3]           = {0.0};
                                double LocalCoordSystem[3][3] = {{0.0}, {0.0}, {0.0}};
                                NormalDir[0] = other_to_me_vect[0];  
                                NormalDir[1] = other_to_me_vect[1];
                                NormalDir[2] = other_to_me_vect[2];
                                
                                
                                double Auxiliar[3][2] = {{0.0}, {0.0}, {0.0}};
                                double stress_projection[2] = {0.0};
                                GeometryFunctions::ComputeContactLocalCoordSystem(NormalDir, LocalCoordSystem);  
                                
                                    
                                //assumció de que la suma de sigmaX i sigmaZ en un pla sempre dona el mateix valor.
                                              
                                //vector ortogonal 1 = LocalCoordSystem[0]
                                //vector ortogonal 2 = LocalCoordSystem[1]
                                
                                //NOTE:puc fer un clean si un valor es massa petit
                                
                                for (int i=0;i<2;i++)//only for 0 and 1, dos auxiliars
                                { 
                                            
                                    for (int j=0;j<3;j++)//for 0,1,2. Component dels auxiliars
                                    {
                                      for (int u=0;u<3;u++)
                                      {
                      
                                      
                                      Auxiliar[j][i] += mSymmStressTensor[j][u]*LocalCoordSystem[i][u];
                                    
                                      }

                                        for (int k=0;k<3;k++)//for 0,1,2.
                                        {
                                                        
                                            stress_projection[i] += Auxiliar[k][i]*LocalCoordSystem[i][k];
                                      
                                        } 
                                    }  
                                    
                                }
                                
                                double Normal_Contact_Contribution = -1.0*calculation_area*equiv_poisson*(stress_projection[0]+stress_projection[1]); 
                                
                                //storing the value in the bar and doing the mean
                              
                                if(this->Id() < neighbour_iterator->Id())
                                {
                                  
                                  lock_p_weak->GetValue(LOW_POISSON_FORCE) = Normal_Contact_Contribution; //crec que ja té la direcció cap a on ha danar el veí.
                                    
                                    if(*mpTimeStep==500)
                                    {
                                    KRATOS_WATCH(lock_p_weak->GetValue(LOW_POISSON_FORCE))             
                                    }
                                  
                                    
                                }
                                else
                                {
                                  
                                  lock_p_weak->GetValue(HIGH_POISSON_FORCE) = Normal_Contact_Contribution;
                                    if(*mpTimeStep==500)
                                    {
                                    KRATOS_WATCH(lock_p_weak->GetValue(HIGH_POISSON_FORCE))   
                                    
                                    }            
                                }
                                                
                                found = true;
                                
                                break;
                                
                
                          
                        }// if ( this->GetValue(CONTINUUM_INI_NEIGHBOURS_IDS)[index_area] == int(neighbour_iterator->Id()) ) 
                                            
                    }//for every ini_neighbour      
                
                
                if(found == false)
                {
                  
                  KRATOS_WATCH("ERROR!!!!!!!!!!!!!   THIS IS ONLY VALID FOR KNOWN INITIAL NEIGHBOURS, NO NEW ONES")KRATOS_WATCH(*mpTimeStep)
                  
                }
                
                //NOTE: this area is provisionally obtained diferent from each side of the contact. bars are not yet ready for mpi.
                
                
                
                  
                  //this->GetValue(POISSON_NORMAL_FORCE) += 0.5*Normal_Contact_Contribution*LocalCoordSystem[2]; //seria el mateix que other to me normalitzat oi?
                  //this->GetValue(POISSON_NORMAL_FORCE) += 0.5*Normal_Contact_Contribution*LocalCoordSystem[2];
                  
                  //neighbour_iterator->GetValue(POISSON_NORMAL_FORCE) += 0.5*Normal_Contact_Contribution;
                  
                  
            
                
            } // for every neighbour    
            
            */
        } //SymmetrizeTensor
              
  
    
void SphericContinuumParticle::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) //Note: this function is only called once by now in the continuum strategy
      { 
   
        const ProcessInfo& r_process_info = rCurrentProcessInfo;


         MemberDeclarationFirstStep(r_process_info); //MSI 1: el cases va fer un custom memberdeclarationfirststep per ficar tot aixo...
           
         //DEM_CONTINUUM
         
         if (!mInitializedVariablesFlag){

         mDempack = r_process_info[DEMPACK_OPTION]; 
         
         mpCurrentTime =&(r_process_info[TIME]); 
                 
         if(mDempack)
         {
           
           mDempack_damping = r_process_info[DEMPACK_DAMPING]; 
           mDempack_global_damping = r_process_info[DEMPACK_GLOBAL_DAMPING]; 

         }
         
         mpActivateSearch = &(r_process_info[ACTIVATE_SEARCH]);

         mGamma1 = r_process_info[DONZE_G1];
         mGamma2 = r_process_info[DONZE_G2];
         mGamma3 = r_process_info[DONZE_G3];
         mMaxDef = r_process_info[DONZE_MAX_DEF];

         mMagicFactorPoisson  = r_process_info[DEM_MAGIC_FACTOR_POISSON];
         
         mSkinSphere = &(this->GetValue(SKIN_SPHERE));
           

         mFinalSimulationTime = r_process_info[FINAL_SIMULATION_TIME];
         
         mContinuumGroup        = this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_CONTINUUM);             

         mpCaseOption                   = r_process_info[CASE_OPTION]; 
         
         mContactMeshOption             = r_process_info[CONTACT_MESH_OPTION];
         
         mTriaxialOption                = r_process_info[TRIAXIAL_TEST_OPTION];

         if (mTriaxialOption )
         {
             mFinalPressureTime      = 0.01*r_process_info[TIME_INCREASING_RATIO] * mFinalSimulationTime; 
         }
         
         AuxiliaryFunctions::SwitchCase(mpCaseOption, mDeltaOption, mContinuumSimulationOption);
         
         mContactInternalFriccion       = r_process_info[CONTACT_INTERNAL_FRICC]*M_PI/180;
         
         mSinContactInternalFriccion    = sin(mContactInternalFriccion);
         mCosContactInternalFriccion    = cos(mContactInternalFriccion);
         mTanContactInternalFriccion    = tan(mContactInternalFriccion);
         
         mFailureCriterionOption        = r_process_info[FAILURE_CRITERION_OPTION];

         mTensionLimit                  = r_process_info[CONTACT_SIGMA_MIN]*1e6; //N/m2
         mCompressionLimit              = r_process_info[CONTACT_SIGMA_MAX]*1e6;
         mTauZero                       = r_process_info[CONTACT_TAU_ZERO]*1e6;           
         
         if( mpCaseOption !=0 ) 
          {
            
            SetInitialContacts( r_process_info);
            //NeighNeighMapping( r_process_info);//MSIMSI DEBUG
            
          }
   
         mInitializedVariablesFlag = 1;
         
         
         //nonlinear parameters:
         
         if(mElasticityType == 2)
         {
           
            mN1 = r_process_info[SLOPE_FRACTION_N1];
            mN2 = r_process_info[SLOPE_FRACTION_N2];
            mN3 = r_process_info[SLOPE_FRACTION_N3];
            mC1 = r_process_info[SLOPE_LIMIT_COEFF_C1];
            mC2 = r_process_info[SLOPE_LIMIT_COEFF_C2];
            mC3 = r_process_info[SLOPE_LIMIT_COEFF_C2]; 
            mYoungPlastic = r_process_info[YOUNG_MODULUS_PLASTIC];
            mPlasticityLimit = r_process_info[PLASTIC_YIELD_STRESS];
            mDamageMaxDisplacementFactor = r_process_info[DAMAGE_FACTOR];
           
         }
         
           
        }// if (!mInitializedVariablesFlag)
         
   
          //#C6: Initialize Volume... Stresstraintensors...MSIMSI 7

      }//void SphericContinuumParticle::InitializeSolutionStep(ProcessInfo& r_process_info)
   
   
    void SphericContinuumParticle::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo)   
      {

          //this->GetGeometry()[0].FastGetSolutionStepValue(EXPORT_PARTICLE_FAILURE_ID) = double(this->GetValue(PARTICLE_FAILURE_ID)); //temporarily unused
          
  
          if(rCurrentProcessInfo[PRINT_SKIN_SPHERE]==1)
          {
            
            this->GetGeometry()[0].FastGetSolutionStepValue(EXPORT_SKIN_SPHERE) = double(*mSkinSphere);  
            
          }
          
          if (rCurrentProcessInfo[PRINT_GROUP_ID] == 1)
          {
              this->GetGeometry()[0].FastGetSolutionStepValue(EXPORT_GROUP_ID) = double(this->GetGeometry()[0].FastGetSolutionStepValue(GROUP_ID));
          }
          
         
        
          //this->GetGeometry()[0].FastGetSolutionStepValue(NUM_OF_NEIGH) = this->GetValue(NEIGHBOUR_ELEMENTS).size();
         /*
          if( mContactMeshOption ==1 && rCurrentProcessInfo[STRESS_STRAIN_OPTION] )
          {
      
          this->GetGeometry()[0].FastGetSolutionStepValue(DEM_STRESS_XX) =  mStressTensor[0][0];
          this->GetGeometry()[0].FastGetSolutionStepValue(DEM_STRESS_XY) =  mStressTensor[0][1];
          this->GetGeometry()[0].FastGetSolutionStepValue(DEM_STRESS_XZ) =  mStressTensor[0][2];
          this->GetGeometry()[0].FastGetSolutionStepValue(DEM_STRESS_YX) =  mStressTensor[1][0];
          this->GetGeometry()[0].FastGetSolutionStepValue(DEM_STRESS_YY) =  mStressTensor[1][1];
          this->GetGeometry()[0].FastGetSolutionStepValue(DEM_STRESS_YZ) =  mStressTensor[1][2];
          this->GetGeometry()[0].FastGetSolutionStepValue(DEM_STRESS_ZX) =  mStressTensor[2][0];
          this->GetGeometry()[0].FastGetSolutionStepValue(DEM_STRESS_ZY) =  mStressTensor[2][1];
          this->GetGeometry()[0].FastGetSolutionStepValue(DEM_STRESS_ZZ) =  mStressTensor[2][2];
          
          }
          */
           if(rCurrentProcessInfo[PRINT_RADIAL_DISPLACEMENT]==1)
          {
            
            double X = this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_X);
            double Z = this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Z);
            
            this->GetGeometry()[0].FastGetSolutionStepValue(RADIAL_DISPLACEMENT) = sqrt(X*X+Z*Z);
            
          }
           // the elemental variable is copied to a nodal variable in order to export the results onto GiD Post. Also a casting to double is necessary for GiD interpretation.
      }
  
  //VELL:
   
  void SphericContinuumParticle::ComputeNewNeighboursHistoricalData() 
  {

  ParticleWeakVectorType& TempNeighbours = mTempNeighbours;
  ParticleWeakVectorType& neighbour_elements = this->GetValue(NEIGHBOUR_ELEMENTS);
  
  TempNeighbours.swap(neighbour_elements); 
  
  unsigned int temp_size = TempNeighbours.size();
  
  neighbour_elements.clear(); 
        
  unsigned int neighbour_counter       = 0;
      
  std::vector<int>&                  temp_neighbours_ids = mTempNeighboursIds;
  std::vector<double>&               temp_neighbours_delta = mTempNeighboursDelta;
  std::vector<int>&                  temp_neighbours_failure_id = mTempNeighboursFailureId;
  std::vector<array_1d<double, 3> >& temp_neighbours_contact_forces = mTempNeighboursContactForces;
  std::vector<int>&                  temp_neighbours_mapping = mTempNeighboursMapping;
  std::vector<int>&                  temp_cont_neighbours_mapping = mTempContNeighboursMapping;
  
  temp_neighbours_ids.resize(temp_size);
  temp_neighbours_delta.resize(temp_size);
  temp_neighbours_failure_id.resize(temp_size);
  temp_neighbours_contact_forces.resize(temp_size);
  temp_neighbours_mapping.resize(temp_size);
  temp_cont_neighbours_mapping.resize(temp_size);
  

  array_1d<double, 3> vector_of_zeros;
  vector_of_zeros[0]                   = 0.0;
  vector_of_zeros[1]                   = 0.0;
  vector_of_zeros[2]                   = 0.0;
  
  for (ParticleWeakIteratorType i = TempNeighbours.begin(); i != TempNeighbours.end(); i++)
  
  {

    double                ini_delta           = 0.0;
    int                   failure_id          = 1;
    array_1d<double, 3>   neigh_forces        = vector_of_zeros;
    double                mapping_new_ini     = -1;  
    double                mapping_new_cont    = -1;

    //Loop Over Initial Neighbours
      //unsigned int start_searching_here = 0; //only to be used if neighbours are already sorted
      
      
      //for (unsigned int k = start_searching_here; k != mIniNeighbourIds.size(); k++) //only to be used if neighbours are already sorted
    for (unsigned int k = 0; k != mIniNeighbourIds.size(); k++) 
      {
        //if (static_cast<int>((i)->Id()) < mIniNeighbourIds[k])  break;         //theoretically useful but it loses a lot of time   
        if (  (i)->Id() == mIniNeighbourIds[k]) //****
        {                               
          ini_delta  = mIniNeighbourDelta[k];
          failure_id = mIniNeighbourFailureId[k];
          mapping_new_ini = k; 
          mapping_new_cont = mIniNeighbourToIniContinuum[k];
          //start_searching_here = k + 1;      //only to be used if neighbours are already sorted           
          break;
        }
      }
                
    //Loop Over Last time-step Neighbours
      //start_searching_here = 0;     //only to be used if neighbours are already sorted       
      //for (unsigned int j = start_searching_here; j != mOldNeighbourIds.size(); j++) //only to be used if neighbours are already sorted
      for (unsigned int j = 0; j != mOldNeighbourIds.size(); j++)
      {
        //if (static_cast<int>(i->Id()) < mOldNeighbourIds[j]) break;  //theoretically useful but it loses a lot of time    
        if ( i->Id() == mOldNeighbourIds[j])
        {
          neigh_forces = mOldNeighbourContactForces[j];
          //start_searching_here = j + 1; //only to be used if neighbours are already sorted
          break;
        }
      }
      
      //Judge if its neighbour            
      double other_radius                 = i->GetGeometry()[0].FastGetSolutionStepValue(RADIUS);
      double radius_sum                   = mRadius + other_radius;
      array_1d<double,3> other_to_me_vect = this->GetGeometry()(0)->Coordinates() - i->GetGeometry()(0)->Coordinates();
      double distance                     = sqrt(other_to_me_vect[0] * other_to_me_vect[0] + other_to_me_vect[1] * other_to_me_vect[1] + other_to_me_vect[2] * other_to_me_vect[2]);
      double indentation                  = radius_sum - distance - ini_delta;
      
      if ( indentation > 0.0 || failure_id == 0 )  //WE NEED TO SET A NUMERICAL TOLERANCE FUNCTION OF THE RADIUS.  MSIMSI 10
      {
      
          neighbour_elements.push_back(*(i.base()));                
          
          //temp_neighbours_ids[neighbour_counter]              = static_cast<int>((i)->Id());
          temp_neighbours_ids[neighbour_counter]              = ((i)->Id());
          temp_neighbours_mapping[neighbour_counter]          = mapping_new_ini;
          temp_cont_neighbours_mapping[neighbour_counter]     = mapping_new_cont;                
          temp_neighbours_delta[neighbour_counter]            = ini_delta;
          temp_neighbours_failure_id[neighbour_counter]       = failure_id;
          temp_neighbours_contact_forces[neighbour_counter]   = neigh_forces;
          
          neighbour_counter++;
          
      }

    }//for ParticleWeakIteratorType i
    
    int final_size = neighbour_elements.size();
    temp_neighbours_ids.resize(final_size);
    temp_neighbours_delta.resize(final_size);
    temp_neighbours_failure_id.resize(final_size);
    temp_neighbours_contact_forces.resize(final_size);
    temp_neighbours_mapping.resize(final_size);
    temp_cont_neighbours_mapping.resize(final_size);
    
    mMapping_New_Ini.swap(temp_neighbours_mapping);
    mMapping_New_Cont.swap(temp_cont_neighbours_mapping);
    mOldNeighbourIds.swap(temp_neighbours_ids);
    mNeighbourDelta.swap(temp_neighbours_delta);
    mNeighbourFailureId.swap(temp_neighbours_failure_id);
    mOldNeighbourContactForces.swap(temp_neighbours_contact_forces);
    
  } //ComputeNewNeighboursHistoricalData
  
 /*
  //RIC!!!!!
  void SphericContinuumParticle::ComputeNewNeighboursHistoricalData() //NOTA: LOOP SOBRE TOTS ELS VEINS PROVISIONALS, TEN KEDERAS UNS QUANTS FENT PUSHBACK. ALS VECTORS DELTA ETC.. HI HAS DE POSAR
     //LA POSICIÓ DELS QUE SON DEFINITIUS.
     {

        ParticleWeakVectorType& TempNeighbours = this->GetValue(NEIGHBOUR_ELEMENTS);

        unsigned int neighbour_counter       = 0;
        unsigned int temp_neighbour_counter  = 0;

        unsigned int temp_size = TempNeighbours.size();

        //vector<int>                  temp_neighbours_mapping(temp_size);
        //vector<int>                  temp_cont_neighbours_mapping(temp_size);
        //vector<int>                  temp_neighbours_ids(temp_size);
        //vector<double>               temp_neighbours_delta(temp_size,0.0);
        //vector<int>                  temp_neighbours_failure_id(temp_size,1);
        //vector<array_1d<double, 3> > temp_neighbours_contact_forces;
        //temp_neighbours_contact_forces.resize(temp_size);              //TODO::Ric. No em deixa fer el //boost::numeric::ublas::matrix<double, check_counter, 3 >
  
        std::vector<int>&                  temp_neighbours_ids = mTempNeighboursIds;
        std::vector<double>&               temp_neighbours_delta = mTempNeighboursDelta;
        std::vector<int>&                  temp_neighbours_failure_id = mTempNeighboursFailureId;
        std::vector<array_1d<double, 3> >& temp_neighbours_contact_forces = mTempNeighboursContactForces;
        std::vector<int>&                  temp_neighbours_mapping = mTempNeighboursMapping;
        std::vector<int>&                  temp_cont_neighbours_mapping = mTempContNeighboursMapping;   
        
        temp_neighbours_ids.resize(temp_size);
        temp_neighbours_delta.resize(temp_size);
        temp_neighbours_failure_id.resize(temp_size);
        temp_neighbours_contact_forces.resize(temp_size);
        temp_neighbours_mapping.resize(temp_size);
        temp_cont_neighbours_mapping.resize(temp_size);
        
        double                ini_delta           = 0.0;
        int                   failure_id          = 1;
        array_1d<double, 3>   neigh_forces        (3,0.0);
        double                mapping_new_ini     = -1;  
        double                mapping_new_cont    = -1;
               

        for (ParticleWeakIteratorType i = TempNeighbours.begin(); i != TempNeighbours.end(); i++)
        
        {

          ini_delta        = 0.0;
          failure_id       = 1;
          neigh_forces[0]  = 0.0;
          neigh_forces[1]  = 0.0;
          neigh_forces[2]  = 0.0;
          mapping_new_ini  = -1;  
          mapping_new_cont = -1;

          //Loop Over Initial Neighbours

          for (unsigned int k = 0; k != mIniNeighbourIds.size(); k++)
          {
                        
            if (static_cast<int>((i)->Id()) == mIniNeighbourIds[k])
            {               
              
              ini_delta  = mIniNeighbourDelta[k];
              failure_id = mIniNeighbourFailureId[k];
              mapping_new_ini = k; 
              mapping_new_cont = mIniNeighbourToIniContinuum[k];
              
              break;
            }

          }
    
          //Judge if its neighbour
          
          double other_radius                 = i->GetGeometry()[0].FastGetSolutionStepValue(RADIUS);
          double radius_sum                   = mRadius + other_radius;
          array_1d<double,3> other_to_me_vect = this->GetGeometry()(0)->Coordinates() - i->GetGeometry()(0)->Coordinates();
          double distance                     = sqrt(other_to_me_vect[0] * other_to_me_vect[0] + other_to_me_vect[1] * other_to_me_vect[1] + other_to_me_vect[2] * other_to_me_vect[2]);
          double indentation                  = radius_sum - distance - ini_delta;
          
          if ( indentation > 0.0 || failure_id == 0 )  //WE NEED TO SET A NUMERICAL TOLERANCE FUNCTION OF THE RADIUS.  MSIMSI 10
          {
        
            //Loop Over Last time-step Neighbours
          
            for (unsigned int j = 0; j != mOldNeighbourIds.size(); j++)
            {
              if (static_cast<int>(i->Id()) == mOldNeighbourIds[j])
              {
                neigh_forces = mOldNeighbourContactForces[j];
                break;
              }

            }
            
            if(neighbour_counter != temp_neighbour_counter)
            {
              
              (*(TempNeighbours.ptr_begin() + neighbour_counter)).swap( (*(TempNeighbours.ptr_begin() + temp_neighbour_counter))) ;     
                              
            }
            
            temp_neighbours_mapping[neighbour_counter]          = mapping_new_ini;
            temp_cont_neighbours_mapping[neighbour_counter]     = mapping_new_cont;
            temp_neighbours_ids[neighbour_counter]              = static_cast<int>((i)->Id());
            temp_neighbours_delta[neighbour_counter]            = ini_delta;
            temp_neighbours_failure_id[neighbour_counter]       = failure_id;
            temp_neighbours_contact_forces[neighbour_counter]   = neigh_forces;
            
            neighbour_counter++;
            
          }
            
            temp_neighbour_counter++;

        }
        
        TempNeighbours.erase(TempNeighbours.begin()+neighbour_counter, TempNeighbours.end());
     
        if(mMapping_New_Ini.size() != neighbour_counter)   //si en comptes de fer tot aixo fes un resize del mMapping_New_... etc... no quedaria tallada la part ke no minteressa i hagues pogut ferlo servir ampliat i ja sta. fer resize i quedarme amb lo bo.
        {
            mMapping_New_Ini.resize(temp_size);
            mMapping_New_Cont.resize(temp_size);
            mOldNeighbourIds.resize(temp_size);
            mNeighbourDelta.resize(temp_size);
            mNeighbourFailureId.resize(temp_size);
            mOldNeighbourContactForces.resize(temp_size);
            
            for(unsigned int w=0; w<neighbour_counter; w++)
            {
              mMapping_New_Ini[w]           = temp_neighbours_mapping[w];
              mMapping_New_Cont[w]          = temp_cont_neighbours_mapping[w];
              mOldNeighbourIds[w]           = temp_neighbours_ids[w];
              mNeighbourDelta[w]            = temp_neighbours_delta[w];
              mNeighbourFailureId[w]        = temp_neighbours_failure_id[w];
              mOldNeighbourContactForces[w] = temp_neighbours_contact_forces[w];
            }
        }
          
        else
        {
            mMapping_New_Ini.swap(temp_neighbours_mapping);
            mMapping_New_Cont.swap(temp_cont_neighbours_mapping);
            mOldNeighbourIds.swap(temp_neighbours_ids);
            mNeighbourDelta.swap(temp_neighbours_delta);
            mNeighbourFailureId.swap(temp_neighbours_failure_id);
            mOldNeighbourContactForces.swap(temp_neighbours_contact_forces);
        }

      } //ComputeNewNeighboursHistoricalData*/

  
   
      void SphericContinuumParticle::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo)
      {
        //KRATOS_WATCH(rVariable)
        
        KRATOS_TRY

        if (rVariable == DELTA_TIME)
        {

            double coeff = rCurrentProcessInfo[NODAL_MASS_COEFF];

            if (coeff>1.0)
            {
                KRATOS_ERROR(std::runtime_error,"The coefficient assigned for vitual mass is larger than one, virtual_mass_coeff= ",coeff)
            }

            else if ((coeff==1.0) && (rCurrentProcessInfo[VIRTUAL_MASS_OPTION]))
            {
                Output = 9.0E09;
            }

            else
            {

                if (rCurrentProcessInfo[VIRTUAL_MASS_OPTION])
                {
                    mRealMass = mRealMass/(1-coeff);
                }

                double K = mYoung * M_PI * this->GetGeometry()(0)->FastGetSolutionStepValue(RADIUS); //M. Error, should be the same that the local definition.

                if (rCurrentProcessInfo[GLOBAL_VARIABLES_OPTION]==1)
                    K = rCurrentProcessInfo[GLOBAL_KN];

                Output = 0.34 * sqrt( mRealMass / K);

                if(rCurrentProcessInfo[ROTATION_OPTION] == 1)
                {
                    Output = Output * 0.5; //factor for critical time step when rotation is allowed.
                }
            }
            
        }//CRITICAL DELTA CALCULATION

        if (rVariable == PARTICLE_ROTATION_DAMP_RATIO)
        {
            //ApplyLocalMomentsDamping( rCurrentProcessInfo ); MSIMSI
        } //DAMPING

        if (rVariable == DEM_STRESS_XX)  //operations with the stress_strain tensors
        {
          
            SymmetrizeTensor( rCurrentProcessInfo );
            
          
        } //EULER_ANGLES
        
        if (rVariable == DUMMY_DEBUG_DOUBLE) //Dummy variable for debugging  MSIMSI DEBUG
        {
          
          CheckPairWiseBreaking();
        
        }
        
        if (rVariable == MEAN_CONTACT_AREA)
        {

            int my_id = this->Id();
            bool im_skin = bool(this->GetValue(SKIN_SPHERE));
            ParticleWeakVectorType& r_continuum_ini_neighbours    = this->GetValue(CONTINUUM_INI_NEIGHBOUR_ELEMENTS);
            
            size_t cont_neigh_index = 0; 
            
            for(ParticleWeakIteratorType ini_cont_neighbour_iterator = r_continuum_ini_neighbours.begin();
              ini_cont_neighbour_iterator != r_continuum_ini_neighbours.end(); ini_cont_neighbour_iterator++)
              
            {
                
                  Element::Pointer lock_p_weak = (this->GetGeometry()[0].GetValue(NODE_TO_NEIGH_ELEMENT_POINTER)(cont_neigh_index)).lock();
                 
                  if(!rCurrentProcessInfo[AREA_CALCULATED_FLAG])
                  {
                  
                    bool neigh_is_skin = bool(ini_cont_neighbour_iterator->GetValue(SKIN_SPHERE));
                    
                    int neigh_id = ini_cont_neighbour_iterator->Id();
                                        
                    if ( (im_skin && neigh_is_skin) || ( !im_skin && !neigh_is_skin ) )
                    {
                                             
                      if( my_id < neigh_id )
                        {

                          lock_p_weak -> GetValue(LOCAL_CONTACT_AREA_LOW) = mcont_ini_neigh_area[cont_neigh_index];
                                                    
                        } // if my id < neigh id
                        
                        else
                        {

                          lock_p_weak -> GetValue(LOCAL_CONTACT_AREA_HIGH) = mcont_ini_neigh_area[cont_neigh_index];
                          
                        }
                        
                      
                    } //both skin or both inner.
                    
                    else if( !im_skin && neigh_is_skin ) //we will store both the same only comming from the inner to the skin.
                    {
                      
                
                      lock_p_weak -> GetValue(LOCAL_CONTACT_AREA_HIGH) = mcont_ini_neigh_area[cont_neigh_index];
                      lock_p_weak -> GetValue(LOCAL_CONTACT_AREA_LOW) = mcont_ini_neigh_area[cont_neigh_index];
                      
                                          
                    } //neigh skin

                      
                }//if(first_time)
                
                else //last operation
                {
                  
                  mcont_ini_neigh_area[cont_neigh_index] = lock_p_weak->GetValue(MEAN_CONTACT_AREA);
                
                }
                
                cont_neigh_index++;
            
            }//loop neigh.
          
          
        } //MEAN_CONTACT_AREA
        
        if (rVariable == LOCAL_CONTACT_AREA_HIGH)
        {
            
            Output = AreaDebugging( rCurrentProcessInfo);
        }
          
           if (rVariable == DEMPACK_DAMPING)
        {
                  
             array_1d<double, 3>& total_force = this->GetGeometry()(0)->FastGetSolutionStepValue(TOTAL_FORCES); //Includes all elastic, damping, but not external (gravity)
             array_1d<double, 3>& velocity = this->GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY);
                           

            for (int i = 0; i<3; i++)
            {
                if ( this->GetGeometry()(0)->pGetDof(VELOCITY_Y)->IsFixed() == false ){
                    total_force[i] = total_force[i] - mDempack_global_damping*fabs(total_force[i])*GeometryFunctions::sign(velocity[i]);                                
                }
            }
        } 
            
    
          
          
          
     
      KRATOS_CATCH("")
      
      }//calculate
      
           
     void SphericContinuumParticle::Calculate(const Variable<Vector >& rVariable, Vector& Output,
                           const ProcessInfo& rCurrentProcessInfo)
    {
       
      if (rVariable == DEM_AREA_VECTOR)  //weighting area.
          {
            
              Output = mcont_ini_neigh_area;
                   
          } //EULER_ANGLES
      
      
    }//calculate Output vector.
      

         
     void SphericContinuumParticle::ComputeAdditionalForces(array_1d<double, 3>& contact_force, array_1d<double, 3>& contact_moment,
                                                            array_1d<double, 3>& additionally_applied_force, array_1d<double, 3>& additionally_applied_moment, ProcessInfo& rCurrentProcessInfo)
    {
          const array_1d<double,3>& gravity         = rCurrentProcessInfo[GRAVITY];


          if(mTriaxialOption && *mSkinSphere) //could be applified to selected particles.
          {
            
            ComputePressureForces(additionally_applied_force, rCurrentProcessInfo);
            
          }
                
          if( mRotationOption != 0 && mRotationSpringOption != 0 )
          {
              //ComputeParticleRotationSpring(); MSI: #C2
          }
          
          //CharacteristicParticleFailureId(rCurrentProcessInfo);
          
          //noalias(additionally_applied_force) += mRealMass * gravity;  //MSIMSI 1: Test that this works
          additionally_applied_force[0] += mRealMass * gravity[0];
          additionally_applied_force[1] += mRealMass * gravity[1];
          additionally_applied_force[2] += mRealMass * gravity[2];
      }
 
      void SphericContinuumParticle::CustomInitialize()
      {         
          
          double& mSectionalInertia         = this->GetGeometry()(0)->FastGetSolutionStepValue(PARTICLE_INERTIA);   
          mSectionalInertia                 = 0.25 * M_PI * mRadius * mRadius * mRadius  * mRadius ;    
          
          double& mRepresentative_Volume    = this->GetGeometry()(0)->FastGetSolutionStepValue(REPRESENTATIVE_VOLUME);             
          mRepresentative_Volume            = 0.0;        
          
      }

      
      void SphericContinuumParticle::EvaluateFailureCriteria(double LocalElasticContactForce[3],double ShearForceNow,double calculation_area, int i_neighbour_count,double& contact_sigma, double& contact_tau,double& failure_criterion_state, bool& sliding, int mapping_new_ini)
      {             

             //(1) MOHR-COULOMB FAILURE: (we don't consider rotational spring!!!!! here) need to be thought.
      
              if (mFailureCriterionOption==1)  //MOHR-COULOMB
              {   
                  contact_tau = ShearForceNow/(calculation_area);
                  contact_sigma = LocalElasticContactForce[2]/(calculation_area);

                  double sigma_max, sigma_min;

                  if (LocalElasticContactForce[2]>=0)
                  {
                      sigma_max = contact_sigma;
                      sigma_min = 0;
                  }
                  else 
                  {
                      sigma_max = 0;
                      sigma_min = contact_sigma;
                  }

                  //change into principal stresses

                  double centre = 0.5*(sigma_max + sigma_min);
                  double radius = sqrt( (sigma_max - centre)*(sigma_max - centre) + contact_tau*contact_tau   ) ;

                  double sigma_I = centre + radius;
                  double sigma_II = centre - radius;
                 
                  // Check:
                 
                  double distance_to_failure = ( mTauZero/(mTanContactInternalFriccion) + centre )*mSinContactInternalFriccion;
              
                  failure_criterion_state = radius/distance_to_failure;
    
                  
                if ( sigma_I - sigma_II >= 2*mTauZero*mCosContactInternalFriccion + (sigma_I + sigma_II)*mSinContactInternalFriccion )
                  {
                      //breaks

                      mNeighbourFailureId[i_neighbour_count] = 5; //mohr coulomb   
                      mIniNeighbourFailureId[ mapping_new_ini ] = 5;
                      failure_criterion_state = 1.0;
                      sliding = true ;

                  }
  
                  
              } //MOHR-COULOMB
        
              ///(2) UNCOUPLED FRACTURE
                      
              if (mFailureCriterionOption==2)//UNCOUPLED FRACTURE and DEMPACK
              {    
 
                  contact_tau = ShearForceNow/(calculation_area);
                  
                  contact_sigma = LocalElasticContactForce[2]/(calculation_area);

                  //double mTauZero = 0.5*sqrt(mCompressionLimit*mTensionLimit); 
                               
                  if (LocalElasticContactForce[2]>=0)
                  {
                      double tau_strength = mTauZero+mTanContactInternalFriccion*contact_sigma;
                      
                      failure_criterion_state = contact_tau/tau_strength;
                                            
                      if(contact_tau>tau_strength)
                      {
                          mNeighbourFailureId[i_neighbour_count] = 2; //shear in compression
                          mIniNeighbourFailureId[ mapping_new_ini ] = 2;
                          failure_criterion_state = 1.0;
                          sliding = true;
                      }
                  } //positive sigmas
                  
                  else //negative sigmas
                  {
  
                        double tau_strength = mTauZero;

                        failure_criterion_state = GeometryFunctions::max(contact_tau/tau_strength, -contact_sigma/mTensionLimit) ;
                        
                        if(contact_tau > tau_strength)
                        {
                            mNeighbourFailureId[i_neighbour_count] = 3;  //shear failure tension
                            mIniNeighbourFailureId[ mapping_new_ini ] = 3;
                            sliding = true;
                            failure_criterion_state = 1.0;
                      
                          //Amb Dempack la fractura tracció és el limit del dany i es mira al calcul de forces...
                          
                            
                         /*  
                            
                            if(contact_sigma<-mTensionLimit && mElasticityType<2)
                            {
                                mNeighbourFailureId[i_neighbour_count] = 12; //both shear and tension
                                mIniNeighbourFailureId[ mapping_new_ini ] = 12;
                                failure_criterion_state = 1.0;
                            } //both shear and tension
                          */
                          
                          
                        }
                        
                      /*
                        else if (contact_sigma<-mTensionLimit && mElasticityType<2)
                        {
                            mNeighbourFailureId[i_neighbour_count] = 4; //tension failure
                            mIniNeighbourFailureId[ mapping_new_ini ] = 4;
                            sliding = true;
                            failure_criterion_state = 1.0;
                            
                        }*/
                     
                      
                  } //negative values of sigma              
          
              } //UNCOUPLED FRACTURE
              
    
      }
      
      void SphericContinuumParticle::CalculateOnContactElements(unsigned int neighbour_iterator_id, size_t i_neighbour_count, int mapping_new_cont, double LocalElasticContactForce[3], 
                                                          double  contact_sigma, double  contact_tau, double failure_criterion_state, double acumulated_damage)
      {
      KRATOS_TRY
       //obtaining pointer to contact element.

       Element::Pointer lock_p_weak = (this->GetGeometry()[0].GetValue(NODE_TO_NEIGH_ELEMENT_POINTER)(mapping_new_cont)).lock();
      
       if ( (this->GetGeometry()[0].GetValue(NODE_TO_NEIGH_ELEMENT_POINTER)(mapping_new_cont)).expired() ) 
           KRATOS_WATCH("You are using a pointer that points to nowhere. Something has to be Fix!!!")
                  
       if( this->Id() < neighbour_iterator_id )  // Since areas are the same, the values are the same and we only store from lower ids.
        {
            //COPY VARIABLES LOW
                                        
            //storing values:
                
            //HIGH-LOW variables
                  
            lock_p_weak->GetValue(LOCAL_CONTACT_FORCE)[0] = LocalElasticContactForce[0];
            lock_p_weak->GetValue(LOCAL_CONTACT_FORCE)[1] = LocalElasticContactForce[1];
            lock_p_weak->GetValue(LOCAL_CONTACT_FORCE)[2] = LocalElasticContactForce[2];
                   
  
            lock_p_weak->GetValue(CONTACT_SIGMA) = contact_sigma;
            lock_p_weak->GetValue(CONTACT_TAU)   = contact_tau;

            lock_p_weak->GetValue(CONTACT_FAILURE) = (mNeighbourFailureId[i_neighbour_count]);                                        
            lock_p_weak->GetValue(FAILURE_CRITERION_STATE) = failure_criterion_state;
            if( ( acumulated_damage > lock_p_weak->GetValue(UNIDIMENSIONAL_DAMAGE) ) || (*mpTimeStep == 0) )
             { lock_p_weak->GetValue(UNIDIMENSIONAL_DAMAGE) = acumulated_damage; }
              
        
            
                  
        } // if Target Id < Neigh Id
//         else   
//         {
//             //COPY VARIABLES HIGH 
//                   
//             //HIGH-LOW variables
//                   
//             lock_p_weak->GetValue(LOCAL_CONTACT_FORCE_HIGH)[0] = LocalElasticContactForce[0];
//             lock_p_weak->GetValue(LOCAL_CONTACT_FORCE_HIGH)[1] = LocalElasticContactForce[1];
//             lock_p_weak->GetValue(LOCAL_CONTACT_FORCE_HIGH)[2] = LocalElasticContactForce[2];
//                                                   
//             //COMBINED MEAN       
//   
//             lock_p_weak->GetValue(CONTACT_TAU)                  = contact_sigma;
//             
//                                   
//         }
        

        ///////////////////////////////////////////////////////////////////////////////////////////////      
        //QUESTION!::: M.S.I.     
        //what happens with the initial continuum contacs that now are not found becouse they are broken....
        //should be assured that they become 0 when they break and this value keeps.
        ///////////////////////////////////////////////////////////////////////////////////////////////

      KRATOS_CATCH("")
      
      }//CalculateOnContactElements
      
      
      
      
      
      /**
       * ComputeStressStrain
       * @param mStressTensor StressTensor matrix
       * @param Representative_Volume NO_SE_QUE_ES
       **/
      void SphericContinuumParticle::ComputeStressStrain(double mStressTensor[3][3],ProcessInfo& rCurrentProcessInfo)
      {
          if(rCurrentProcessInfo[STRESS_STRAIN_OPTION] == 1) // if stress_strain_options ON 
          {
              double& Representative_Volume = this->GetGeometry()[0].FastGetSolutionStepValue(REPRESENTATIVE_VOLUME);
          
            
              if ( ( Representative_Volume <= 0.0 ))// && ( this->GetValue(SKIN_SPHERE) == 0 ) )
              {
                  this->GetGeometry()(0)->FastGetSolutionStepValue(GROUP_ID) = 15;
                  KRATOS_WATCH(this->Id())
                  KRATOS_WATCH("Negatiu volume")
                  KRATOS_WATCH(*mpTimeStep)
              }
                
              else
              {
                  for (int i=0; i<3; i++)
                  {
                      for (int j=0; j<3; j++)
                      {
                          //KRATOS_WATCH(Representative_Volume)
                          //KRATOS_WATCH(mStressTensor[i][j])
                          mStressTensor[i][j] = (1/Representative_Volume)*0.5*(mStressTensor [i][j] + mStressTensor[j][i]);  //THIS WAY THE TENSOR BECOMES SYMMETRIC
                      }
                  }       
              }  
          }//if stress_strain_options
      }

      /**
       * Calculates Stress Tensor
       * @param mStressTensor StressTensor matrix
       * @param GlobalElasticContactForce GlobalElasticContactForce vector
       * @param other_to_me_vect NO_SE_QUE_ES
       * @param distance distance with the neighbour
       * @param radius_sum neighbour's radius plus our radius
       * @param calculation_area cirrected contact area between particles
       * @param Representative_Volume NO_SE_QUE_ES
       * @param neighbour_iterator neighbour pointer
       **/
      void SphericContinuumParticle::StressTensorOperations(double mStressTensor[3][3],
                                                                    double GlobalElasticContactForce[3],
                                                                    array_1d<double,3> &other_to_me_vect,
                                                                    const double &distance,
                                                                    const double &radius_sum,
                                                                    const double &calculation_area,
                                                                    ParticleWeakIteratorType neighbour_iterator, ProcessInfo& rCurrentProcessInfo)
      {
          if(rCurrentProcessInfo[STRESS_STRAIN_OPTION]==1) //TODO: Change this with class members or flags
          {
              double gap                  = distance - radius_sum;
            
              array_1d<double,3> normal_vector_on_contact =  -1 * other_to_me_vect; //outwards
            
              double Dummy_Dummy = 0.0;
              GeometryFunctions::normalize(normal_vector_on_contact,Dummy_Dummy); // Normalize to unitary module

              array_1d<double,3> coord_target    = this->GetGeometry()(0)->Coordinates();
              array_1d<double,3> coord_neigh     = neighbour_iterator->GetGeometry()(0)->Coordinates();
              //double neigh_radius                = neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(RADIUS);
              double target_radius               = this->GetGeometry()(0)->FastGetSolutionStepValue(RADIUS);
              array_1d<double,3> x_centroid      = (target_radius + 0.5*gap) * normal_vector_on_contact;
            
              //KRATOS_WATCH(kn_el)
              double result_product = GeometryFunctions::DotProduct(x_centroid,normal_vector_on_contact);
            
              double& Representative_Volume = this->GetGeometry()[0].FastGetSolutionStepValue(REPRESENTATIVE_VOLUME);
          
              Representative_Volume = Representative_Volume + 0.33333333333333 * (result_product * calculation_area);
            
              for (int i=0; i<3; i++)
              {
                  for (int j=0; j<3; j++)
                  {   
                      mStressTensor[i][j] += (x_centroid[j]) * GlobalElasticContactForce[i]; //ref: Katalin Bagi 1995 Mean stress tensor           
                  }
              }
          }
      }
      
      /**
       * Calculate Adds Posisson Contribution
       * @param LocalCoordSystem DEFINITION HERE
       * @param GlobalContactForce DEFINITION HERE
       * @param GlobalElasticContactForce DEFINITION HERE
       * @param ViscoDampingGlobalContactForce DEFINITION HERE
       * @param rContactForce DEFINITION HERE
       * @param damp_forces DEFINITION HERE
       **/
      
      
      void SphericContinuumParticle::AddPoissonContribution(double LocalCoordSystem[3][3], //TODO: array 1d i boost::numeric::ublas::bounded_matrix o Matrix(,) si no saps el tamany a alocar.
                                                                    double GlobalContactForce[3],
                                                                    double GlobalElasticContactForce[3],
                                                                    double ViscoDampingGlobalContactForce[3],
                                                                    array_1d<double, 3>& rContactForce,
                                                                    array_1d<double,3>& damp_forces)
      {     
          double PoissonContactForce[3]        = {0.0};
          double GlobalContactPoissonForce[3]  = {0.0};
          
          PoissonContactForce [2] = 0.0; //neigh_poisson_contribution; //MSI: cambiar a que tingui poisson de la matriu membre
          GeometryFunctions::VectorLocal2Global(LocalCoordSystem, PoissonContactForce, GlobalContactPoissonForce);

          
          // TODO:noalias(rContactForce) += GlobalContactForce PER PODER FER AIXO AHN DE SER TOTS ARRAY1D
          rContactForce[0] += GlobalContactForce[0] + GlobalContactPoissonForce[0];
          rContactForce[1] += GlobalContactForce[1] + GlobalContactPoissonForce[1];
          rContactForce[2] += GlobalContactForce[2] + GlobalContactPoissonForce[2];
       
      }     
      
      
      
   void SphericContinuumParticle::ComputePressureForces(array_1d<double, 3>& externally_applied_force, ProcessInfo& rCurrentProcessInfo) 
        
     {
 
        array_1d<double, 3> total_externally_applied_force  = this->GetGeometry()[0].FastGetSolutionStepValue(EXTERNAL_APPLIED_FORCE); 
        
//         double magnitude = 0.0;
//         
//         GeometryFunctions::module(total_externally_applied_force, magnitude);
//         
//      
        
        double time_now = rCurrentProcessInfo[TIME]; //MSIMSI 1 I tried to do a *mpTIME

    
        if( mFinalPressureTime <= 1e-10 )
        {
          
            KRATOS_WATCH("WARNING: SIMULATION TIME TO CLOSE TO ZERO")
          
        }
        
        else if ( (mFinalPressureTime > 1e-10) && (time_now < mFinalPressureTime) )
        {
          
        externally_applied_force = AuxiliaryFunctions::LinearTimeIncreasingFunction(total_externally_applied_force, time_now, mFinalPressureTime);
   
//               if(this->Id() ==10118)
//               {
//                       double magnitude;
//                       GeometryFunctions::module(externally_applied_force, magnitude);
//                       KRATOS_WATCH(time_now)
//                   KRATOS_WATCH(magnitude)
//               }
        }       
        else
        {
          externally_applied_force = total_externally_applied_force;
          
        }
        
     } //SphericContinuumParticle::ComputePressureForces
    

    void SphericContinuumParticle::PlasticityAndDamage1D(double LocalElasticContactForce[3], double kn_el, double indentation, double calculation_area, double radius_sum_i, double& failure_criterion_state, double& acumulated_damage, int i_neighbour_count, double mapping_new_cont, double mapping_new_ini )
    {

      //VARIABLES
      
      // a guardar:
  
      double kn_b = kn_el / mN1;
      double kn_c = kn_el / mN2;
      double kn_d = kn_el / mN3;
      double kp_el = mYoungPlastic * kn_el;

      double CompressionLimit_1 = mC1 * mCompressionLimit;
      double CompressionLimit_2 = mC2 * mCompressionLimit;
      double CompressionLimit_3 = mC3 * mCompressionLimit;
      
      double Yields_el = mPlasticityLimit * mCompressionLimit * calculation_area;
        
      double Ncstr1_el = CompressionLimit_1 * calculation_area;
      double Ncstr2_el = CompressionLimit_2 * calculation_area;
      double Ncstr3_el = CompressionLimit_3 * calculation_area;
      double Ntstr_el  = mTensionLimit * calculation_area;
      
      //double sigma_a = (kn_el * indentation)/(calculation_area);
      //double sigma_b = mCompressionLimit_1 + kn_b*(indentation/calculation_area - mCompressionLimit_1/kn_el);

      double u_max = mHistory[mapping_new_cont][0];
      
      double& fn = LocalElasticContactForce[2]; //[2] means 'normal' contact force
                
      
      if( indentation >= 0.0 ) //COMPRESSION
      {
        
          fn = kn_el * indentation;
          
          double u_ela1 = Ncstr1_el/kn_el;;
          double u_ela2 = u_ela1 + (Ncstr2_el-Ncstr1_el)/(kn_b);
          double u_ela3 = u_ela2 + (Ncstr3_el-Ncstr2_el)/(kn_c);

          if ( ( indentation > u_max ) || (*mpTimeStep <= 1) )//maximum historical intentation OR first step
            
          {

            mHistory[mapping_new_cont][0]  = indentation;             // Guarda el treshold del màxim desplaçament
            
            
            if (indentation > u_ela3) //4rt tram
            {
              
              fn = Ncstr3_el + ( indentation - u_ela3 )*kn_d;
              
            }
            else if (indentation > u_ela2) //3r tram
            {
              
              fn = Ncstr2_el + ( indentation - u_ela2 )*kn_c;
              
            }
            else
            {    
              if( indentation > u_ela1) //2n tram
              {
                fn = Ncstr1_el + (indentation - u_ela1)*kn_b;
              
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
            
            fn = indentation * kn_el*(1 -  mHistory[mapping_new_cont][2]);  // normal adhesive force (gap +)
            
          }
      
        }//Tension
     
    }//PlasticityAndDamage1D

      
      
     void SphericContinuumParticle::CheckPairWiseBreaking()  //MSIMSI DEBUG
     
      {              
      
        //MSIMSI DEBUG
      }
      
       double SphericContinuumParticle::AreaDebugging(const ProcessInfo& rCurrentProcessInfo)  //MSIMSI DEBUG
     
      { 
         //DEBUG MEDICIÓ
           
          double area_vertical_centre = 0.0;
          
          if (this->GetGeometry()[0].FastGetSolutionStepValue(GROUP_ID)==5)
                
          {

            
           size_t i_neighbour_count = 0;

           ParticleWeakVectorType& r_continuum_ini_neighbours    = this->GetValue(CONTINUUM_INI_NEIGHBOUR_ELEMENTS);
            
           
            for(ParticleWeakIteratorType neighbour_iterator = r_continuum_ini_neighbours.begin();
              neighbour_iterator != r_continuum_ini_neighbours.end(); neighbour_iterator++)
            {
            
                if ( neighbour_iterator->GetGeometry()[0].FastGetSolutionStepValue(GROUP_ID) != 5 )
                  
                {
                    
                                  
                  array_1d<double,3> other_to_me_vect = this->GetGeometry()(0)->Coordinates() - neighbour_iterator->GetGeometry()(0)->Coordinates();
                  array_1d<double,3> normal_vector_on_contact =  -1 * other_to_me_vect; //outwards     
                            
                  double Dummy_Dummy = 0.0;
                  GeometryFunctions::normalize(normal_vector_on_contact,Dummy_Dummy); // Normalize to unitary module
                        
                  
                  area_vertical_centre += mcont_ini_neigh_area[i_neighbour_count]*fabs(normal_vector_on_contact[1]); //X(0), Y(1), Z(2)

                }                   
                
                i_neighbour_count++;
                
            }  // loop neigh
      
          }  //group_id == 5
          
          
          return area_vertical_centre;

      } //AreaDebugging
      
      
    
      
      
      
      
      
      
    void SphericContinuumParticle::Calculate(const Variable<array_1d<double, 3 > >& rVariable, array_1d<double, 3 > & Output, const ProcessInfo& rCurrentProcessInfo){}
    void SphericContinuumParticle::Calculate(const Variable<Matrix >& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo){}

      
}  // namespace Kratos.



//NOTE::
/*
 * #1: Here, Initial_delta is expected to be positive if it is embeding and negative if there's a separation.
 * #2: 0.25 is becouse we take only the half of the equivalent mRadius, corresponding to the case of one sphere with mRadius Requivalent and other = mRadius 0.
 * #3: For detached particles we enter only if the indentation is > 0. For attached particles we enter only if the particle is still attached.
 * #4: we use incremental calculation. YADE develops a complicated "absolute method"
 * #5: We only store in the initial neighbours array the ones that are cohesive or the ones that have possitive or negative initial indentation. In other words,
 *     the non-cohesive ones with 0 indentation (<some tolerance) don't have to be stored since we can treat it indistinctly from other new neighbours that the particle in stury would meet.
*/






           
//#C2: ComputeParticleRotationSpring;
           
//           void SphericContinuumParticle::ComputeParticleRotationSpring() // shan de corregir areees etc etc etc...
//       {
//         //double dt                           = rCurrentProcessInfo[DELTA_TIME]; //C.F.: neew
//         /*
//                     c=objecte_contacte(particula,vei)
// 
//             força=(c.RADI)*3;  //M: idea: to create a class contact, create objects of contacts with all the paramaters. easier...
//                                 /no puc amb MPI oi? pk hauria de passar punters...
//           */
//           ParticleWeakVectorType& mrNeighbours                = this->GetValue(NEIGHBOUR_ELEMENTS);
// 
//           array_1d<double, 3 > & mRota_Moment = GetGeometry()(0)->FastGetSolutionStepValue(PARTICLE_MOMENT);
// 
// 
//           Vector & mRotaSpringFailureType  = this->GetValue(PARTICLE_ROTATE_SPRING_FAILURE_TYPE);
// 
//           size_t i_neighbour_count = 0;
// 
//           for(ParticleWeakIteratorType ineighbour = mrNeighbours.begin(); ineighbour != mrNeighbours.end(); ineighbour++)
//           {
// 
//               //if(mIfInitalContact[i_neighbour_count] == 1 && mRotaSpringFailureType[i_neighbour_count] == 0) ///M.S:NEWWWW,  if THE SPRING BRAKES... NO MORE CONTRIBUION.
//               //if( mRotaSpringFailureType[i_neighbour_count] == 0) //M.S: CAL FICAR A INITIALIZE QUE SIGUI 1 I DESPRES INITIAL CONTACTS POSAR 0 SI NECESITEN, IGUAL QUE FAILURE NORMAL.
//               //mmm.. what about the other failure types? if a contact is broken due to shear or tensile, it cant be a bending
//               {
//                 
//                 
//                   array_1d<double, 3 > & mRotaSpringMoment  = this->GetValue(PARTICLE_ROTATE_SPRING_MOMENT)[ i_neighbour_count ];
// 
//                   double other_radius    = ineighbour->GetGeometry()(0)->FastGetSolutionStepValue(RADIUS);
//                   double other_young     = ineighbour->GetGeometry()[0].FastGetSolutionStepValue(YOUNG_MODULUS);
//                   double other_poisson   = ineighbour->GetGeometry()[0].FastGetSolutionStepValue(POISSON_RATIO);
//                   double other_tension   = ineighbour->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_TENSION); //MSI: these variables are eliminated.
//                   double other_cohesion  = ineighbour->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_COHESION);
//                   double other_inertia   = ineighbour->GetGeometry()(0)->FastGetSolutionStepValue(PARTICLE_INERTIA);
// 
//                   double equiv_tension  = (mTension  + other_tension ) * 0.5;
//                   double equiv_cohesion = (mCohesion + other_cohesion) * 0.5;
// 
//                   double equiv_radius     = (mRadius + other_radius) * 0.5 ;
//                   double equiv_area       = M_PI * equiv_radius * equiv_radius;
//                   double equiv_poisson    = (mPoisson + other_poisson) * 0.5 ;
//                   double equiv_young      = (mYoung  + other_young)  * 0.5;
// 
//                   double kn_el               = mMagicFactor*equiv_young * equiv_area / (2.0 * equiv_radius);
//                   double ks               = kn_el / (2.0 * (1.0 + equiv_poisson));
// 
//                   array_1d<double,3> other_to_me_vect = GetGeometry()(0)->Coordinates() - ineighbour->GetGeometry()(0)->Coordinates();
// 
//                   /////Cfeng: Forming the Local Contact Coordinate system
//                   double NormalDir[3]           = {0.0};
//                   double LocalCoordSystem[3][3] = {{0.0}, {0.0}, {0.0}};
//                   NormalDir[0] = other_to_me_vect[0];
//                   NormalDir[1] = other_to_me_vect[1];
//                   NormalDir[2] = other_to_me_vect[2];
//                   GeometryFunctions::ComputeContactLocalCoordSystem(NormalDir, LocalCoordSystem);
// 
//                   double LocalRotaSpringMoment[3]     = {0.0};
//                   double GlobalRotaSpringMoment[3]    = {0.0};
//                   double GlobalRotaSpringMomentOld[3] = {0.0};
// 
//                   double DeltRotaDisp[3] = {0.0};
//                   //double DeltRotaDisp2[3] = {0.0};
// 
//                   double LocalDeltRotaDisp[3] = {0.0};
// 
//                   double TargetDeltRotaDist[3] = {0.0};
//                   double NeighbourDeltRotaDist[3] = {0.0};
//                 
//                   for (int i=0;i<3;i++)
//                   {
//                       TargetDeltRotaDist[i] = this->GetGeometry()(0)->FastGetSolutionStepValue(DELTA_ROTA_DISPLACEMENT)[i];
//                       NeighbourDeltRotaDist[i] = ineighbour->GetGeometry()(0)->FastGetSolutionStepValue(DELTA_ROTA_DISPLACEMENT)[i];
//                       DeltRotaDisp[i] =  - ( TargetDeltRotaDist[i] - NeighbourDeltRotaDist[i] );
//                   }
//               
//                   GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, DeltRotaDisp, LocalDeltRotaDisp);
// 
//                   GlobalRotaSpringMomentOld[0] = mRotaSpringMoment[ 0 ];
// 
//                   GlobalRotaSpringMomentOld[1] = mRotaSpringMoment[ 1 ];
// 
//                   GlobalRotaSpringMomentOld[2] = mRotaSpringMoment[ 2 ];
// 
//                   GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, GlobalRotaSpringMomentOld, LocalRotaSpringMoment);
// 
// 
//                   double Inertia_I = (mSectionalInertia + other_inertia) * 0.5;
// 
//                   double Inertia_J = Inertia_I * 2.0;
// 
// 
//                   LocalRotaSpringMoment[0] +=  - Inertia_I * LocalDeltRotaDisp[0] * kn_el / equiv_area;
// 
//                   LocalRotaSpringMoment[1] +=  - Inertia_I * LocalDeltRotaDisp[1] * kn_el / equiv_area;
// 
//                   LocalRotaSpringMoment[2] +=  - Inertia_J * LocalDeltRotaDisp[2] * ks / equiv_area;
// 
// 
//                   ////Judge if the rotate spring is broken or not
//                   double GlobalElasticContactForce[3]  = {0.0};
//                   double LocalElasticContactForce [3]  = {0.0};
// 
//                   GlobalElasticContactForce[0] = this->GetValue(PARTICLE_CONTACT_FORCES)[ i_neighbour_count ][ 0 ];
//                   GlobalElasticContactForce[1] = this->GetValue(PARTICLE_CONTACT_FORCES)[ i_neighbour_count ][ 1 ];
//                   GlobalElasticContactForce[2] = this->GetValue(PARTICLE_CONTACT_FORCES)[ i_neighbour_count ][ 2 ]; //GlobalElasticContactForce[2] = mContactForces[3 * i_neighbour_count  + 2 ];
//                   GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, GlobalElasticContactForce, LocalElasticContactForce);
// 
//                   double ForceN  = LocalElasticContactForce[2];
//                   double ForceS  = sqrt( LocalElasticContactForce[0] * LocalElasticContactForce[0] + LocalElasticContactForce[1] * LocalElasticContactForce[1]);
//                   double MomentS = sqrt(LocalRotaSpringMoment[0] * LocalRotaSpringMoment[0] + LocalRotaSpringMoment[1] * LocalRotaSpringMoment[1]);
//                   double MomentN = LocalRotaSpringMoment[2];
// 
//                   //////bending stress and axial stress add together, use edge of the bar will failure first
//                   double TensiMax = -ForceN / equiv_area + MomentS        / Inertia_I * equiv_radius;
//                   double ShearMax = ForceS  / equiv_area + fabs(MomentN)  / Inertia_J * equiv_radius;
// 
//                   if(TensiMax > equiv_tension || ShearMax > equiv_cohesion)
//                   {
//                       mRotaSpringFailureType[i_neighbour_count] = 1;
// 
//                       LocalRotaSpringMoment[0] = 0.0;
//                       LocalRotaSpringMoment[1] = 0.0;
//                       LocalRotaSpringMoment[2] = 0.0;
//                   }
// 
//                   GeometryFunctions::VectorLocal2Global(LocalCoordSystem, LocalRotaSpringMoment, GlobalRotaSpringMoment);
// 
//                   mRotaSpringMoment[ 0 ] = GlobalRotaSpringMoment[0];
//                   mRotaSpringMoment[ 1 ] = GlobalRotaSpringMoment[1];
//                   mRotaSpringMoment[ 2 ] = GlobalRotaSpringMoment[2];
// 
//                   ////feedback, contact moment----induce by rotation spring
//                   mRota_Moment[0] -= GlobalRotaSpringMoment[0];
//                   mRota_Moment[1] -= GlobalRotaSpringMoment[1];
//                   mRota_Moment[2] -= GlobalRotaSpringMoment[2];
//               }
// 
//               i_neighbour_count++;
//           }
//       }//ComputeParticleRotationSpring
//      
/*         
 * 
 * 
 */

// #C3: InitializeContactElements : akesta funcio comensava abans de evaluatedeltadisplacement per poisson etc pero no crec ke faci falta ara.
    
    
                  

         
    
            
   /*
   //#C6 : initalizesolutionstep del volume strain stress tensor...  
         
          double& Representative_Volume = this->GetGeometry()[0].FastGetSolutionStepValue(REPRESENTATIVE_VOLUME);
         
          Representative_Volume = 0.0;
          
          for (int i=0; i<3; i++)
          {
          
              for (int j=0; j<3; j++)
              {
                   mStressTensor[i][j] = 0.0;
                   mSymmStressTensor[i][j] = 0.0;
              }
          
           }
         //MSIMSI
           
*/
   
   /* // MSI #C4
              else if(mElasticityType==3)
              {
                  double current_def            = (initial_dist - distance)/initial_dist ;
                  double kn1                    = kn_el;
                  double kn2                    = kn1 * mGamma3 + kn1 * mGamma1 * exp(mGamma2 * (current_def - mMaxDef));
                         if (kn2 > kn1) {
                             kn2                = kn1;
                         };
                  double max_dist               = initial_dist * (1 - mMaxDef);
                  
                  initial_dist
                  double kn_plas                = kn1;  // modificable en el futuro con un input

                  LocalElasticContactForce[2] = kn1 * (initial_dist - distance); // = kn1 * indentation; fuerza en parte lineal

                  if (distance < mHistDist) // vemos si esta en carga, comparando la distancia actual entre particulas con la maxima historica
                  {
                      mHistDist = distance;
                      
                      if( distance < max_dist )  //  se supera el limite para le cambio de pendiente.
                      {
                          NonlinearNormalForceCalculation (LocalElasticContactForce, kn1, kn2, distance, max_dist, initial_dist);                       
                      }
                      else
                      {
                          LocalElasticContactForce[2] = kn1 * (initial_dist - distance); // fuerza en parte lineal
                      }
                   }
                   double hist_fn = LocalElasticContactForce[2]; //guardamos el maximo historico de fuerza fn
                   else  // esta en descarga, la distancia entre particulas aumenta respecto la historica, current_dist > mHistDist
                   {
                      if (hist_fn > 0 );   //  fuerza normal esta en el rango de compresion de la curva
                      {
                          plast_dist = max_dist; // initial_dist*(1-plast_def) distancia associada al valor de plast_def impuesto, por ahora coincide con el cambio de pendiente
                          if (plast_dist > mHistDist ) // mientras se este por encima de la maxima historica, estamos en plasticidad.
                          {
                              fn = hist_fn + kn_el_plas*(mHistDist - distance); // en descarga: 500 - kn_plas(10 - 12). en carga 500 - kn_el(10 - 11) pero con distance > mHistDistance
                          }
                          else  // esta en descarga pero no en la zona plastica, descarga por la linea elastica
                          {
                              if( distance < max_dist )  // se supera el limite para le cambio de pendiente, mientras plast_dist=max_dist nunca pasara, nunca descargara por kn2
                              {
                                  NonlinearNormalForceCalculation (LocalElasticContactForce, kn1, kn2, distance, max_dist, initial_dist);
                              }
                              else // descarga por la primera rama elastica
                              {
                                  LocalElasticContactForce[2] = kn1 * (initial_dist - distance); // fuerza en parte lineal
                              }
                          }
                      }
                   }

               }
*/
              
      
                   
