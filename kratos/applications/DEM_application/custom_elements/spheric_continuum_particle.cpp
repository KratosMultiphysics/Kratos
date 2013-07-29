
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

      void SphericContinuumParticle::SetInitialContacts( ProcessInfo& rCurrentProcessInfo  ) //vull ficar que sigui zero si no son veins cohesius.
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

            //SAVING THE INICIAL NEIGHBOURS, THE DELTAS AND THE FAILURE ID
                  
            for(ParticleWeakIteratorType_ptr ineighbour = mrNeighbours.ptr_begin();  //loop over the temp neighbours and store into a initial_neighbours vector.
            ineighbour != mrNeighbours.ptr_end(); ineighbour++)
            {
               
              
                array_1d<double,3> other_to_me_vect = this->GetGeometry()(0)->Coordinates() - ((*ineighbour).lock())->GetGeometry()(0)->Coordinates();
                double distance                     = sqrt(other_to_me_vect[0] * other_to_me_vect[0] +
                                                      other_to_me_vect[1] * other_to_me_vect[1] +
                                                      other_to_me_vect[2] * other_to_me_vect[2]);

                double radius_sum                   = mRadius + ((*ineighbour).lock())->GetGeometry()(0)->GetSolutionStepValue(RADIUS);
                double initial_delta                = radius_sum - distance;

                int r_other_continuum_group         = ((*ineighbour).lock())->GetGeometry()(0)->GetSolutionStepValue(PARTICLE_CONTINUUM);
   
                if( ( (r_other_continuum_group == mContinuumGroup) && (mContinuumGroup != 0) ) || ( fabs(initial_delta)>1.0e-6 ) ) // which would be the appropiate tolerance?
                
                //for the ones that need to be stored... #5
                {
                    ini_size++;
                    
                    mIniNeighbourIds.resize(ini_size);         
                    mIniNeighbourDelta.resize(ini_size);
                    mIniNeighbourFailureId.resize(ini_size);           
                    mMapping_New_Ini.resize(ini_size); 
        
                    mIniNeighbourIds[ini_size - 1]          = ((*ineighbour).lock())->Id();
                    mIniNeighbourDelta[ini_size - 1]        = 0.0;
                    mIniNeighbourFailureId[ini_size - 1]    = 1;
                    mMapping_New_Ini[ini_size - 1]          = ini_size - 1;
                                 
                    if (mDeltaOption == true)
                    {            
                        mIniNeighbourDelta[ini_size - 1]    = initial_delta;

                                                
                    }

                    if (mContinuumSimulationOption == true)
                    {
                        
                        if ( (r_other_continuum_group == mContinuumGroup) && (mContinuumGroup != 0) )
                        {
                            
                            mIniNeighbourFailureId[ini_size - 1]=0;
                           
                            mFailureId=0; // if a cohesive contact exist, the FailureId becomes 0.          
                            continuum_ini_size++;
                                
                            r_continuum_ini_neighbours.push_back(*ineighbour);
                              
                            r_continuum_ini_neighbours_ids.resize(continuum_ini_size);                               
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
              
                ContactAreaWeighting(rCurrentProcessInfo);
                  
            } 
            
         
        }//SetInitialContacts
 

      void SphericContinuumParticle::ContactAreaWeighting(const ProcessInfo& rCurrentProcessInfo) //MISMI 10: POOYAN this could be done by calculating on the bars. not looking at the neighbous of my neighbours.
      { 

          double alpha = 1.0;
          double external_sphere_area = 4*M_PI*mRadius*mRadius;  
          
          mtotal_equiv_area = 0.0;

          ParticleWeakVectorType r_continuum_ini_neighbours    = this->GetValue(CONTINUUM_INI_NEIGHBOUR_ELEMENTS);
          int cont_ini_neighbours_size                         = r_continuum_ini_neighbours.size();
          
          mcont_ini_neigh_area.resize(cont_ini_neighbours_size);
          
          //computing the total equivalent area
          
          size_t index = 0;
          
          for(ParticleWeakIteratorType ini_cont_neighbour_iterator = r_continuum_ini_neighbours.begin();     // MSIMSI 99:Could this loop be done during the bar creation in the strategy and so avoid another repetition?
              ini_cont_neighbour_iterator != r_continuum_ini_neighbours.end(); ini_cont_neighbour_iterator++)
          {   
              double other_radius     = ini_cont_neighbour_iterator->GetGeometry()(0)->GetSolutionStepValue(RADIUS);
              double equiv_radius     = 2*mRadius * other_radius / (mRadius + other_radius);        
              double equiv_area       = (0.25)*M_PI * equiv_radius * equiv_radius; //we now take 1/2 of the efective mRadius.
              mtotal_equiv_area       += equiv_area;
          
              mcont_ini_neigh_area[index] = equiv_area; //*
              index++; //*
              
              
          } //for every neighbour
       
       
          if(!mSkinSphere)
          {
          
            AuxiliaryFunctions::CalculateAlphaFactor(cont_ini_neighbours_size, external_sphere_area, mtotal_equiv_area, alpha); 
            
            size_t not_skin_index = 0;
        
            for(ParticleWeakIteratorType ini_cont_neighbour_iterator = r_continuum_ini_neighbours.begin();
                ini_cont_neighbour_iterator != r_continuum_ini_neighbours.end(); ini_cont_neighbour_iterator++)
                
                {      
                    mcont_ini_neigh_area[not_skin_index] = alpha*mcont_ini_neigh_area[not_skin_index];
                    not_skin_index++;  
                    
                } //for every neighbour

          }
          
          else //skin sphere //AIXO HAURIA DE CANVIAR PER L'AREA COPIADA BONA DEL VEI KE NO ES 
          {
          
               
              size_t skin_index = 0; 
              
              for(ParticleWeakIteratorType ini_cont_neighbour_iterator = r_continuum_ini_neighbours.begin();
                ini_cont_neighbour_iterator != r_continuum_ini_neighbours.end(); ini_cont_neighbour_iterator++)
                
                {
                     if(ini_cont_neighbour_iterator->GetValue(SKIN_SPHERE))
                     {
                        alpha            = 1.0*(1.40727)*(external_sphere_area/mtotal_equiv_area)*((double(cont_ini_neighbours_size))/11);
                        mcont_ini_neigh_area[skin_index] = alpha*mcont_ini_neigh_area[skin_index];
                     }
                     
                     else
                     {
                        
                       double neigh_area = 0.0;
                       Vector vector_neigh_area;
                       
                       ParticleWeakVectorType r_cont_neigh_of_my_cont_neigh = ini_cont_neighbour_iterator->GetValue(CONTINUUM_INI_NEIGHBOUR_ELEMENTS);  
                         
                       size_t skin_skin_index = 0;
                         
                       for(ParticleWeakIteratorType r_cont_neigh_of_my_cont_neigh_iterator = r_cont_neigh_of_my_cont_neigh.begin();            //MSIMSI 99: if there would be more loop on the neighbours of my neighbours would be nice to have a mapping to avoid this on the solve loop
                           r_cont_neigh_of_my_cont_neigh_iterator != r_cont_neigh_of_my_cont_neigh.end(); r_cont_neigh_of_my_cont_neigh_iterator++)
                
                          {
                              
                            if(r_cont_neigh_of_my_cont_neigh_iterator->Id() == this->Id() )  //MSIMSI 99:Pooyan, would be better mId than this->Id()?
                            {
 
                              r_cont_neigh_of_my_cont_neigh_iterator->Calculate(DEM_AREA_VECTOR,vector_neigh_area,rCurrentProcessInfo);
                              neigh_area = vector_neigh_area[skin_skin_index];
                               
                              break;
                            }
                                                            
                            skin_skin_index++;
                          }//loop of the cont neighs of my cont neighbour.
  
  
                        mcont_ini_neigh_area[skin_index] = neigh_area;
  
                     }//not skin neighbours of skin particles
            
                skin_index++;
                
                }//loop on cont neighs       
                
          }//skin particles.
   
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
            double indentation                    = (radius_sum - initial_delta) - distance;   //#1
            double equiv_area                     = 0.25*M_PI * equiv_radius * equiv_radius; //#2 
            double corrected_area                 = equiv_area;
            double equiv_mass                     = mSqrtOfRealMass * other_sqrt_of_mass;

            double equiv_young;
            double equiv_poisson;
            double equiv_visco_damp_coeff_normal;
            double equiv_visco_damp_coeff_tangential;
            double equiv_ln_of_restit_coeff;
            double kn;
            double kt;
            double equiv_tg_of_fri_ang;

            double DeltDisp[3]                    = {0.0};
            double RelVel[3]                      = {0.0};

            double NormalDir[3]                   = {0.0};
            double OldNormalDir[3]                = {0.0};

            double LocalCoordSystem[3][3]         = {{0.0}, {0.0}, {0.0}};
            double OldLocalCoordSystem[3][3]      = {{0.0}, {0.0}, {0.0}};

            bool sliding = false;
            
            int mapping = mMapping_New_Ini[i_neighbour_count]; //*
                  
            if(mContactMeshOption==1 && (mapping !=-1))
            {
              corrected_area = mcont_ini_neigh_area[mapping];                            
            }
            
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
                double aux_norm_to_tang;
                
            if (mGlobalVariablesOption){
                kn                                = mGlobalKn;
                kt                                = mGlobalKt;
                aux_norm_to_tang                  = mGlobalAuxNormToTang;
            }

            else {
              
                kn                                = mMagicFactor * equiv_young * corrected_area * radius_sum_i;
                kt                                = kn / (2.0 + equiv_poisson + equiv_poisson);
                aux_norm_to_tang                  = sqrt(kt / kn);
            }
          
            if (mCriticalTimeOption){
                double historic = rCurrentProcessInfo[HISTORICAL_MIN_K];

                if ((kn < historic) || (kt < historic)){
                    historic = fmin(kn, kt);
                }

            }

            if (mLnOfRestitCoeff > 0.0 || other_ln_of_restit_coeff > 0.0){
                equiv_visco_damp_coeff_normal     = 2 * sqrt(equiv_mass * kn);
                equiv_visco_damp_coeff_tangential = equiv_visco_damp_coeff_normal * aux_norm_to_tang; 
            }

            else {
                equiv_visco_damp_coeff_normal     = - 2 * equiv_ln_of_restit_coeff * sqrt(equiv_mass * kn / (equiv_ln_of_restit_coeff * equiv_ln_of_restit_coeff + M_PI * M_PI));

                equiv_visco_damp_coeff_tangential = equiv_visco_damp_coeff_normal * aux_norm_to_tang; 
            }

            EvaluateDeltaDisplacement(DeltDisp, RelVel, NormalDir, OldNormalDir, LocalCoordSystem, OldLocalCoordSystem, other_to_me_vect, vel, delta_displ, neighbour_iterator);

            DisplacementDueToRotation(DeltDisp, OldNormalDir, OldLocalCoordSystem, other_radius, dt, ang_vel, neighbour_iterator);
            
            double LocalDeltDisp[3] = {0.0};
            double LocalElasticContactForce[3]  = {0.0}; // 0: first tangential, // 1: second tangential, // 2: normal force
            double GlobalElasticContactForce[3] = {0.0};
            double LocalRelVel[3] = {0.0};

            GlobalElasticContactForce[0] = mOldNeighbourContactForces[i_neighbour_count][0];  
            GlobalElasticContactForce[1] = mOldNeighbourContactForces[i_neighbour_count][1];
            GlobalElasticContactForce[2] = mOldNeighbourContactForces[i_neighbour_count][2];
            
            GeometryFunctions::VectorGlobal2Local(OldLocalCoordSystem, GlobalElasticContactForce, LocalElasticContactForce); 
            //we recover this way the old local forces projected in the new coordinates in the way they were in the old ones; Now they will be increased if its the necessary
            GeometryFunctions::VectorGlobal2Local(OldLocalCoordSystem, DeltDisp, LocalDeltDisp);
            GeometryFunctions::VectorGlobal2Local(OldLocalCoordSystem, RelVel, LocalRelVel);
            
            /* Translational Forces */
          
            if  (indentation > 0.0 || (mNeighbourFailureId[i_neighbour_count] == 0) )//*  //#3
            {                                                                                                
              
              //Normal Forces
              NormalForceCalculation(LocalElasticContactForce,kn,indentation,mElasticityType);
                              
              //Nonlinear...() MSI #C4
              
              //Tangential Forces   #4
              LocalElasticContactForce[0] += - kt * LocalDeltDisp[0];  // 0: first tangential
              LocalElasticContactForce[1] += - kt * LocalDeltDisp[1];  // 1: second tangential
              
            }

            /* Evaluating Failure for the continuum contacts */
            
              double failure_criterion_state = 0.0;   
      
              double ShearForceNow = sqrt(LocalElasticContactForce[0] * LocalElasticContactForce[0]
                                    +      LocalElasticContactForce[1] * LocalElasticContactForce[1]); 
              
              double Frictional_ShearForceMax = equiv_tg_of_fri_ang * LocalElasticContactForce[2];
              
              if (Frictional_ShearForceMax < 0.0){
                    Frictional_ShearForceMax = 0.0;
                  }
              
              if ( mNeighbourFailureId[i_neighbour_count] != 0 ) //*
                {
                    failure_criterion_state = 1.0;
                                                          
                    if( (ShearForceNow >  Frictional_ShearForceMax) && (ShearForceNow != 0.0) ) 
                    {
                        LocalElasticContactForce[0] = (Frictional_ShearForceMax / ShearForceNow) * LocalElasticContactForce[0];
                        LocalElasticContactForce[1] = (Frictional_ShearForceMax / ShearForceNow )* LocalElasticContactForce[1];
                        sliding = true;
                                    
                    }
                    
                }
            
            if(corrected_area <1e-09) {KRATOS_WATCH(corrected_area) KRATOS_WATCH(this->Id())} // MSIMSI 10

            double contact_tau = 0.0;
            double contact_sigma = 0.0;

            if(mNeighbourFailureId[i_neighbour_count] == 0)
            {                
              EvaluateFailureCriteria(LocalElasticContactForce,ShearForceNow,corrected_area,i_neighbour_count,contact_sigma,contact_tau, failure_criterion_state, sliding, mapping);
            }
      
            // VISCODAMPING (applyied locally)
            
            double ViscoDampingLocalContactForce[3]    = {0.0};

            if  (indentation > 0.0 || (mNeighbourFailureId[i_neighbour_count] == 0) )//*  //#3
            {  
            
            CalculateViscoDamping(LocalRelVel,ViscoDampingLocalContactForce,indentation,equiv_visco_damp_coeff_normal,equiv_visco_damp_coeff_tangential,sliding);
            
            }
  
            // Transforming to global forces and adding up
            double LocalContactForce[3] =                 {0.0};
            double ViscoDampingGlobalContactForce[3] =    {0.0}; 
            double GlobalContactForce[3] =                {0.0};
            
              
            AddUpForcesAndProject(LocalCoordSystem,mOldNeighbourContactForces,LocalContactForce,LocalElasticContactForce,GlobalContactForce,
                                  GlobalElasticContactForce,ViscoDampingLocalContactForce,ViscoDampingGlobalContactForce,rContactForce,rElasticForce,
                                  i_neighbour_count);

            //AddPoissonContribution(LocalCoordSystem, GlobalContactForce, GlobalElasticContactForce, ViscoDampingGlobalContactForce, rContactForce, damp_forces); //MSIMSI 10

            
            ComputeMoments(LocalElasticContactForce,GlobalElasticContactForce,InitialRotaMoment,LocalCoordSystem,other_radius,rContactMoment,neighbour_iterator);
            
            //StressTensorOperations(mStressTensor,GlobalElasticContactForce,other_to_me_vect,distance,radius_sum,corrected_area,neighbour_iterator,rCurrentProcessInfo); //MSISI 10
  
            if(mContactMeshOption==1 && (mapping !=-1)) 
            {

              CalculateOnContactElements( neighbour_iterator ,i_neighbour_count, mapping, LocalElasticContactForce, corrected_area, contact_sigma, contact_tau, failure_criterion_state);
  
            }

            i_neighbour_count++;

        }//for each neighbour
        
        rInitialRotaMoment [0] = InitialRotaMoment [0];
        rInitialRotaMoment [1] = InitialRotaMoment [1];
        rInitialRotaMoment [2] = InitialRotaMoment [2];


        //ComputeStressStrain(mStressTensor, rCurrentProcessInfo);  //MSIMSI 10
        
        KRATOS_CATCH("")         
            
        
      }//ComputeBallToBallContactForce

    
      void SphericContinuumParticle::ApplyLocalMomentsDamping(const ProcessInfo& rCurrentProcessInfo )
      {

          KRATOS_TRY
      
          array_1d<double, 3 > & RotaMoment       = this->GetGeometry()[0].GetSolutionStepValue(PARTICLE_MOMENT);
          double RotaDampRatio                    = this->GetGeometry()[0].GetSolutionStepValue(PARTICLE_ROTATION_DAMP_RATIO);

          // LOCAL DAMPING OPTION FOR THE UNBALANCED FORCES (IN GLOBAL CORDINATES).
        
          for (int iDof = 0; iDof < 3; iDof++)
          {
              if (this->GetGeometry()(0)->GetSolutionStepValue(ANGULAR_VELOCITY)[iDof] > 0.0)
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

     
      void SphericContinuumParticle::SymmetrizeTensor(const ProcessInfo& rCurrentProcessInfo)
      {

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
                double other_radius     = neighbour_iterator->GetGeometry()(0)->GetSolutionStepValue(RADIUS);
                double other_poisson    = neighbour_iterator->GetGeometry()[0].GetSolutionStepValue(POISSON_RATIO);
                double equiv_radius     = 2*mRadius * other_radius / (mRadius + other_radius);
                int size_ini_cont_neigh = this->GetValue(CONTINUUM_INI_NEIGHBOURS_IDS).size();
                double equiv_area       = (0.25)*M_PI * equiv_radius * equiv_radius; // 0.25 is becouse we take only the half of the equivalent mRadius, corresponding to the case of one sphere with mRadius Requivalent and other = mRadius 0.
                double equiv_poisson    = 2* mPoisson * other_poisson / (mPoisson + other_poisson);
                //double equiv_young      = 2 * young * other_young / (young + other_young);
                //bool is_continuum       = false;
                  double corrected_area = equiv_area;
                  
                  bool found = false;
                  
                  for (int index_area=0; index_area<size_ini_cont_neigh; index_area++)
                  {

                        if ( this->GetValue(CONTINUUM_INI_NEIGHBOURS_IDS)[index_area] == int(neighbour_iterator->Id()) ) 
                        {
                            
                          //MPI_CARLOS_ en MPI no pot funcionar lo de symmetrize de moment

                                Element::Pointer lock_p_weak = (this->GetGeometry()[0].GetValue(NODE_TO_NEIGH_ELEMENT_POINTER)(index_area)).lock();

                                corrected_area = lock_p_weak->GetValue(MEAN_CONTACT_AREA);

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
                                
                                double Normal_Contact_Contribution = -1.0*corrected_area*equiv_poisson*(stress_projection[0]+stress_projection[1]); 
                                
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
        }
              
  
    

      void SphericContinuumParticle::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
      { 
        
           MemberDeclarationFirstStep(rCurrentProcessInfo);
           
         //DEM_CONTINUUM
         
         if (!mInitializedVariablesFlag){
         
         
         mSkinSphere = this->GetValue(SKIN_SPHERE);
           
         mFinalSimulationTime = rCurrentProcessInfo[FINAL_SIMULATION_TIME];
         
         mContinuumGroup        = this->GetGeometry()[0].GetSolutionStepValue(PARTICLE_CONTINUUM);             

         
         mpCaseOption                   = &(rCurrentProcessInfo[CASE_OPTION]);   //NOTE: pointer
         
         mContactMeshOption             = rCurrentProcessInfo[CONTACT_MESH_OPTION];
         
         mTriaxialOption                = rCurrentProcessInfo[TRIAXIAL_TEST_OPTION];

         if (mTriaxialOption )
         {
             mInitialPressureTime    = rCurrentProcessInfo[INITIAL_PRESSURE_TIME];
             mFinalPressureTime      = 0.01*rCurrentProcessInfo[TIME_INCREASING_RATIO] * mFinalSimulationTime; 
         }
         
         AuxiliaryFunctions::SwitchCase(*mpCaseOption, mDeltaOption, mContinuumSimulationOption);
         
         mContactInternalFriccion       = rCurrentProcessInfo[CONTACT_INTERNAL_FRICC]*M_PI/180;
         
         mSinContactInternalFriccion    = sin(mContactInternalFriccion);
         mCosContactInternalFriccion    = cos(mContactInternalFriccion);
         mTanContactInternalFriccion    = tan(mContactInternalFriccion);
         
         mFailureCriterionOption        = rCurrentProcessInfo[FAILURE_CRITERION_OPTION];

         mTensionLimit                  = rCurrentProcessInfo[CONTACT_SIGMA_MIN]*1e6; //N/m2
         mCompressionLimit              = rCurrentProcessInfo[CONTACT_SIGMA_MAX]*1e6;
         mTauZero                       = rCurrentProcessInfo[CONTACT_TAU_ZERO]*1e6;           
         
         if( (*mpCaseOption!=0) )
          {
            
            SetInitialContacts( rCurrentProcessInfo);
            
          }
          
         mSwitchPressure = &(rCurrentProcessInfo[SWITCH_PRESSURE]); 
         mInitializedVariablesFlag = 1;
           
        }// if (!mInitializedVariablesFlag)
         
   
          //#C6: Initialize Volume... Stresstraintensors...MSIMSI 7

      }//void SphericContinuumParticle::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
   
   
    void SphericContinuumParticle::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo)   
      {

          //this->GetGeometry()[0].FastGetSolutionStepValue(EXPORT_PARTICLE_FAILURE_ID) = double(this->GetValue(PARTICLE_FAILURE_ID)); //temporarily unused
          
  
          if(rCurrentProcessInfo[PRINT_SKIN_SPHERE]==1)
          {
            
            this->GetGeometry()[0].FastGetSolutionStepValue(EXPORT_SKIN_SPHERE) = double(mSkinSphere);  
            
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
            
            double X = this->GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT_X);
            double Z = this->GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT_Z);
            
            this->GetGeometry()[0].GetSolutionStepValue(RADIAL_DISPLACEMENT) = sqrt(X*X+Z*Z);
            
          }
           // the elemental variable is copied to a nodal variable in order to export the results onto GiD Post. Also a casting to double is necessary for GiD interpretation.
      }

   
   
     void SphericContinuumParticle::ComputeNewNeighboursHistoricalData() //NOTA: LOOP SOBRE TOTS ELS VEINS PROVISIONALS, TEN KEDERAS UNS QUANTS FENT PUSHBACK. ALS VECTORS DELTA ETC.. HI HAS DE POSAR
     //LA POSICIÓ DELS QUE SON DEFINITIUS.
     {
   
       ParticleWeakVectorType TempNeighbours;
       TempNeighbours.swap(this->GetValue(NEIGHBOUR_ELEMENTS)); //GetValue is needed becouse this information comes from the strategy (the search function)
       
       this->GetValue(NEIGHBOUR_ELEMENTS).clear(); 
              
       unsigned int neighbour_counter       = 0;
       
       vector<int>                  temp_neighbours_ids;
       vector<double>               temp_neighbours_delta;
       vector<int>                  temp_neighbours_failure_id;
       vector<array_1d<double, 3> > temp_neighbours_contact_forces;
       vector<int>                  temp_neighbours_mapping;

       array_1d<double, 3> vector_of_zeros;
       vector_of_zeros[0]                   = 0.0;
       vector_of_zeros[1]                   = 0.0;
       vector_of_zeros[2]                   = 0.0;
   
       for (ParticleWeakIteratorType_ptr i = TempNeighbours.ptr_begin(); i != TempNeighbours.ptr_end(); i++)
       
       {

         //neigh_added = false;
          
          double                ini_delta       = 0.0;
          int                   failure_id      = 1;
          array_1d<double, 3>   neigh_forces    = vector_of_zeros;
          double                mapping            = -1;    

          //Loop Over Initial Neighbours

            for (unsigned int k = 0; k != mIniNeighbourIds.size(); k++)
            {
                           
              if (static_cast<int>((*i).lock()->Id()) == mIniNeighbourIds[k])
              {               
                
                ini_delta  = mIniNeighbourDelta[k];
                failure_id = mIniNeighbourFailureId[k];
                mapping = k; 
                
                break;
              }

            }
                           
          //Loop Over Last time-step Neighbours
          
            for (unsigned int j = 0; j != mOldNeighbourIds.size(); j++)
            {
              if (static_cast<int>((*i).lock()->Id()) == mOldNeighbourIds[j])
              {
                neigh_forces = mOldNeighbourContactForces[j];
                break;
              }

            }
            
            //Judge if its neighbour
            
            double other_radius                 = (*i).lock()->GetGeometry()[0].GetSolutionStepValue(RADIUS);
            double radius_sum                   = mRadius + other_radius;
            array_1d<double,3> other_to_me_vect = this->GetGeometry()(0)->Coordinates() - (*i).lock()->GetGeometry()(0)->Coordinates();
            double distance                     = sqrt(other_to_me_vect[0] * other_to_me_vect[0] + other_to_me_vect[1] * other_to_me_vect[1] + other_to_me_vect[2] * other_to_me_vect[2]);
            double indentation                  = radius_sum - distance - ini_delta;
            
            if ( indentation > 0.0 || (indentation < 1.0e-6 && failure_id == 0 ) )  //WE NEED TO SET A NUMERICAL TOLERANCE FUNCTION OF THE RADIUS.  MSIMSI 10
            {
           
                this->GetValue(NEIGHBOUR_ELEMENTS).push_back(*i);
                size_t size = this->GetValue(NEIGHBOUR_ELEMENTS).size();
                
                temp_neighbours_ids.resize(size);
                temp_neighbours_delta.resize(size);
                temp_neighbours_failure_id.resize(size);
                temp_neighbours_contact_forces.resize(size);
                temp_neighbours_mapping.resize(size);
                
                temp_neighbours_ids[neighbour_counter]              = static_cast<int>((*i).lock()->Id());
                temp_neighbours_mapping[neighbour_counter]          = mapping;
                temp_neighbours_delta[neighbour_counter]            = ini_delta;
                temp_neighbours_failure_id[neighbour_counter]       = failure_id;
                temp_neighbours_contact_forces[neighbour_counter]   = neigh_forces;
                
                neighbour_counter++;
                
            }

        }
        
        mMapping_New_Ini.swap(temp_neighbours_mapping);
        mOldNeighbourIds.swap(temp_neighbours_ids);
        mNeighbourDelta.swap(temp_neighbours_delta);
        mNeighbourFailureId.swap(temp_neighbours_failure_id);
        mOldNeighbourContactForces.swap(temp_neighbours_contact_forces);

      } //ComputeNewNeighboursHistoricalData

   
      void SphericContinuumParticle::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo)
      {
          //KRATOS_WATCH(rVariable)

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
          

      }//calculate
      
           
     void SphericContinuumParticle::Calculate(const Variable<Vector >& rVariable, Vector& Output,
                           const ProcessInfo& rCurrentProcessInfo)
    {
       
      if (rVariable == DEM_AREA_VECTOR)  //weighting area.
          {
            
              Output = mcont_ini_neigh_area;
                   
          } //EULER_ANGLES
      
      
    }//calculate Output vector.
      
      
   
     void SphericContinuumParticle::CustomCalculateRightHandSide(array_1d<double, 3>& contact_force, array_1d<double, 3>& contact_moment, 
                                                        array_1d<double, 3>& exterally_applied_force, ProcessInfo& rCurrentProcessInfo)
     {
          if(mTriaxialOption && mSkinSphere) //could be applified to selected particles.
          {
            
            ComputePressureForces(exterally_applied_force, rCurrentProcessInfo);
            
            KRATOS_WATCH(exterally_applied_force)
            
          }
        
        
          if( mRotationOption != 0 && mRotationSpringOption != 0 )
          {
              //ComputeParticleRotationSpring(); MSI: #C2
          }
          
          //CharacteristicParticleFailureId(rCurrentProcessInfo);
      }
 
      void SphericContinuumParticle::CustomInitialize()
      {
         
          double& mSectionalInertia         = GetGeometry()(0)->FastGetSolutionStepValue(PARTICLE_INERTIA);   
          mSectionalInertia                 = 0.25 * M_PI * mRadius * mRadius * mRadius  * mRadius ;    
          
          double& mRepresentative_Volume = this->GetGeometry()[0].GetSolutionStepValue(REPRESENTATIVE_VOLUME);   
          
          mRepresentative_Volume = 0.0;
      
      }

      
      void SphericContinuumParticle::EvaluateFailureCriteria(double LocalElasticContactForce[3],double ShearForceNow,double corrected_area, int i_neighbour_count,double& contact_sigma, double& contact_tau,double& failure_criterion_state, bool& sliding, int mapping)
      {             

             //(1) MOHR-COULOMB FAILURE: (we don't consider rotational spring!!!!! here) need to be thought.
      
              if (mFailureCriterionOption==1)  //MOHR-COULOMB
              {   
                  contact_tau = ShearForceNow/(corrected_area);
                  contact_sigma = LocalElasticContactForce[2]/(corrected_area);

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
                      mIniNeighbourFailureId[ mapping ] = 5;
                      sliding = true ;

                  }
  
                  
              } //MOHR-COULOMB
        
              ///(2) UNCOUPLED FRACTURE
        
              ///vam decidir amb miguel angel de no fer el mapping de les shear fins al pas seguent.. esta correcte? afecta quan trenca?

       
              if (mFailureCriterionOption==2)//UNCOUPLED FRACTURE
              {    
       
                  contact_tau = ShearForceNow/(corrected_area);
                  contact_sigma = LocalElasticContactForce[2]/(corrected_area);

                  //double mTauZero = 0.5*sqrt(mCompressionLimit*mTensionLimit); 
                               
                  if (LocalElasticContactForce[2]>=0)
                  {
                      double tau_strength = mTauZero+mTanContactInternalFriccion*contact_sigma;
                      
                      failure_criterion_state = contact_tau/tau_strength;
                                            
                      if(contact_tau>tau_strength)
                      {
                          mNeighbourFailureId[i_neighbour_count] = 2;
                          mIniNeighbourFailureId[ mapping ] = 2;
                          
                          sliding = true;
                      }
                  } //positive sigmas
                  
                  else //negative sigmas
                  {
                      double tau_strength = mTauZero;
            
                      failure_criterion_state = GeometryFunctions::max(contact_tau/tau_strength, -contact_sigma/mTensionLimit) ;

                      if(contact_tau > tau_strength)
                      {
                          mNeighbourFailureId[i_neighbour_count] = 3;  //shear failure
                          mIniNeighbourFailureId[ mapping ] = 3;
                          sliding = true;
              
                          if(contact_sigma<-mTensionLimit)
                          {
                              mNeighbourFailureId[i_neighbour_count] = 12;
                              mIniNeighbourFailureId[ mapping ] = 12;
                          } //both shear and tension
              
                      }
                   
                      else if(contact_sigma<-mTensionLimit)
                      {
                          mNeighbourFailureId[i_neighbour_count] = 4; //tension failure
                          mIniNeighbourFailureId[ mapping ] = 4;
                          sliding = true;
                           
                      }
                      
                  } //negative values of sigma              
          
              } //UNCOUPLED FRACTURE
              
    
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
      }
      
      void SphericContinuumParticle::CalculateOnContactElements(ParticleWeakIteratorType neighbour_iterator, size_t i_neighbour_count, int mapping, double LocalElasticContactForce[3], 
                                                          double corrected_area, double  contact_sigma, double  contact_tau, double failure_criterion_state)
      {
   
       //obtaining pointer to contact element.

       Element::Pointer lock_p_weak = (this->GetGeometry()[0].GetValue(NODE_TO_NEIGH_ELEMENT_POINTER)(mapping)).lock();
                      
       if( this->Id() < neighbour_iterator->Id() )  // if id pequeña
        {
            //COPY VARIABLES LOW
                                        
            //storing values:
                
            //HIGH-LOW variables
                  
            lock_p_weak->GetValue(LOCAL_CONTACT_FORCE_LOW)[0] = LocalElasticContactForce[0];
            lock_p_weak->GetValue(LOCAL_CONTACT_FORCE_LOW)[1] = LocalElasticContactForce[1];
            lock_p_weak->GetValue(LOCAL_CONTACT_FORCE_LOW)[2] = LocalElasticContactForce[2];
            
            if(*mpTimeStep==0)
            {
            lock_p_weak->GetValue(LOCAL_CONTACT_AREA_LOW) = corrected_area;
            
            }
            
            //COMBINED MEAN          
  
            lock_p_weak->GetValue(CONTACT_SIGMA) += 0.5*contact_sigma;
            lock_p_weak->GetValue(CONTACT_TAU)   += 0.5*contact_tau;
                                                                      
            //UNIQUE VALUES
            
            //1) failure
                lock_p_weak->GetValue(CONTACT_FAILURE) = (mNeighbourFailureId[i_neighbour_count]);                                        
                      
                if(failure_criterion_state<=1.0)
                {
                    lock_p_weak->GetValue(FAILURE_CRITERION_STATE) = failure_criterion_state;
                   
                }
                
                else
                {
                    //KRATOS_WATCH (failure_criterion_state )
                }   
                                
        } // if Target Id < Neigh Id
        else   
        {
            //COPY VARIABLES HIGH 
                  
            //HIGH-LOW variables
                  
            lock_p_weak->GetValue(LOCAL_CONTACT_FORCE_HIGH)[0] = LocalElasticContactForce[0];
            lock_p_weak->GetValue(LOCAL_CONTACT_FORCE_HIGH)[1] = LocalElasticContactForce[1];
            lock_p_weak->GetValue(LOCAL_CONTACT_FORCE_HIGH)[2] = LocalElasticContactForce[2];
            
            if(*mpTimeStep==0)
            {

            lock_p_weak->GetValue(LOCAL_CONTACT_AREA_HIGH) = corrected_area;

            }
            
                                                  
            //COMBINED MEAN       
  
            lock_p_weak->GetValue(CONTACT_SIGMA)                += 0.5*contact_sigma;
            lock_p_weak->GetValue(CONTACT_TAU)                  += 0.5*contact_tau;
            
                                  
        }
        
        //CONTACT AREA
        
          if ( ( *mpTimeStep==0 ) && (!mSkinSphere ) && ( neighbour_iterator->GetValue(SKIN_SPHERE)==0 ) )
              {
                                            
                lock_p_weak->GetValue(MEAN_CONTACT_AREA)   += 0.5*corrected_area;

              }
              
              else if ( ( *mpTimeStep==0 ) && ( mSkinSphere ) && ( neighbour_iterator->GetValue(SKIN_SPHERE)==1 ) )
              {
          
                lock_p_weak->GetValue(MEAN_CONTACT_AREA)   += 0.5*corrected_area;
                                            
              }
              
                      
              
              else if ( ( *mpTimeStep==0 ) && ( !mSkinSphere ) && ( neighbour_iterator->GetValue(SKIN_SPHERE)==1 ) )
              {
              
              lock_p_weak->GetValue(MEAN_CONTACT_AREA)   = corrected_area;
                
              }
    
        ///////////////////////////////////////////////////////////////////////////////////////////////      
        //QUESTION!::: M.S.I.     
        //what happens with the initial continuum contacs that now are not found becouse they are broken....
        //should be assured that they become 0 when they break and this value keeps.
        ///////////////////////////////////////////////////////////////////////////////////////////////


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
              double& Representative_Volume = this->GetGeometry()[0].GetSolutionStepValue(REPRESENTATIVE_VOLUME);
          
            
              if ( ( Representative_Volume <= 0.0 ))// && ( this->GetValue(SKIN_SPHERE) == 0 ) )
              {
                  this->GetGeometry()(0)->GetSolutionStepValue(GROUP_ID) = 15;
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
       * @param corrected_area cirrected contact area between particles
       * @param Representative_Volume NO_SE_QUE_ES
       * @param neighbour_iterator neighbour pointer
       **/
      void SphericContinuumParticle::StressTensorOperations(double mStressTensor[3][3],
                                                                    double GlobalElasticContactForce[3],
                                                                    array_1d<double,3> &other_to_me_vect,
                                                                    const double &distance,
                                                                    const double &radius_sum,
                                                                    const double &corrected_area,
                                                                    ParticleWeakIteratorType neighbour_iterator, ProcessInfo& rCurrentProcessInfo)
      {
          if(rCurrentProcessInfo[STRESS_STRAIN_OPTION]==1) //TODO: Change this with class members or flags
          {
              double gap                  = distance - radius_sum;
            
              array_1d<double,3> normal_vector_on_contact =  -1 * other_to_me_vect; //outwards
            
              double Dummy_Dummy = 0.0;
              GeometryFunctions::norm(normal_vector_on_contact,Dummy_Dummy); // Normalize to unitary module

              array_1d<double,3> coord_target    = this->GetGeometry()(0)->Coordinates();
              array_1d<double,3> coord_neigh     = neighbour_iterator->GetGeometry()(0)->Coordinates();
              //double neigh_radius                = neighbour_iterator->GetGeometry()(0)->GetSolutionStepValue(RADIUS);
              double target_radius               = this->GetGeometry()(0)->GetSolutionStepValue(RADIUS);
              array_1d<double,3> x_centroid      = (target_radius + 0.5*gap) * normal_vector_on_contact;
            
              //KRATOS_WATCH(kn)
              double result_product = GeometryFunctions::DotProduct(x_centroid,normal_vector_on_contact);
            
              double& Representative_Volume = this->GetGeometry()[0].GetSolutionStepValue(REPRESENTATIVE_VOLUME);
          
              Representative_Volume = Representative_Volume + 0.33333333333333 * (result_product * corrected_area);
            
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
      
      
      void SphericContinuumParticle::AddPoissonContribution(double LocalCoordSystem[3][3],
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

          rContactForce[0] += GlobalContactForce[0] + GlobalContactPoissonForce[0];
          rContactForce[1] += GlobalContactForce[1] + GlobalContactPoissonForce[1];
          rContactForce[2] += GlobalContactForce[2] + GlobalContactPoissonForce[2];
          
          //array_1d<double,3>& elastiqueueueueueue; 
          //elastiqueueueueueue[0] += GlobalElasticContactForce[0]; //
          //elastiqueueueueueue[1] += GlobalElasticContactForce[1];
          //elastiqueueueueueue[2] += GlobalElasticContactForce[2];           //MSI: si vols exportar les forces elastiques les has de repondre

          //damp_forces[0] += ViscoDampingGlobalContactForce[0];
          //damp_forces[1] += ViscoDampingGlobalContactForce[1];            //MSI: si les vols guardar... pos ho has dactivar
          //damp_forces[2] += ViscoDampingGlobalContactForce[2];          
      }     
      
      
      
   void SphericContinuumParticle::ComputePressureForces(array_1d<double, 3>& externally_applied_force, ProcessInfo& rCurrentProcessInfo) 
        
     {
 
        array_1d<double, 3> total_externally_applied_force  = this->GetGeometry()[0].GetSolutionStepValue(EXTERNAL_APPLIED_FORCE); 

        double current_time     = rCurrentProcessInfo[TIME];     
  
        if( mFinalPressureTime <= 1e-10 )
        {
          
            KRATOS_WATCH("WARNING: SIMULATION TIME TO CLOSE TO ZERO")
          
        }
        
        else if ( (mFinalPressureTime > 1e-10) && (current_time < mFinalPressureTime) )
        {
          
        externally_applied_force = AuxiliaryFunctions::LinearTimeIncreasingFunction(total_externally_applied_force,mInitialPressureTime,current_time,mFinalPressureTime);
        
        if(mSkinSphere)
        {
          
          
        }
        
        }
        
        else
        {
          externally_applied_force = total_externally_applied_force;
          
          if(*mSwitchPressure==0)
          {
            
            KRATOS_WATCH("Confinement application finished at time :")
            KRATOS_WATCH(current_time)
            *mSwitchPressure = 1;
          }  
        }
        
     } //SphericContinuumParticle::ComputePressureForces
    

      
      
      
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
//               //if(mIfInitalContact[i_neighbour_count] == 1 && mRotaSpringFailureType[i_neighbour_count] == 0) ///M.S:NEWWWW, IF THE SPRING BRAKES... NO MORE CONTRIBUION.
//               //if( mRotaSpringFailureType[i_neighbour_count] == 0) //M.S: CAL FICAR A INITIALIZE QUE SIGUI 1 I DESPRES INITIAL CONTACTS POSAR 0 SI NECESITEN, IGUAL QUE FAILURE NORMAL.
//               //mmm.. what about the other failure types? if a contact is broken due to shear or tensile, it cant be a bending
//               {
//                 
//                 
//                   array_1d<double, 3 > & mRotaSpringMoment  = this->GetValue(PARTICLE_ROTATE_SPRING_MOMENT)[ i_neighbour_count ];
// 
//                   double other_radius    = ineighbour->GetGeometry()(0)->FastGetSolutionStepValue(RADIUS);
//                   double other_young     = ineighbour->GetGeometry()[0].GetSolutionStepValue(YOUNG_MODULUS);
//                   double other_poisson   = ineighbour->GetGeometry()[0].GetSolutionStepValue(POISSON_RATIO);
//                   double other_tension   = ineighbour->GetGeometry()[0].GetSolutionStepValue(PARTICLE_TENSION); //MSI: these variables are eliminated.
//                   double other_cohesion  = ineighbour->GetGeometry()[0].GetSolutionStepValue(PARTICLE_COHESION);
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
//                   double kn               = mMagicFactor*equiv_young * equiv_area / (2.0 * equiv_radius);
//                   double ks               = kn / (2.0 * (1.0 + equiv_poisson));
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
//                   LocalRotaSpringMoment[0] +=  - Inertia_I * LocalDeltRotaDisp[0] * kn / equiv_area;
// 
//                   LocalRotaSpringMoment[1] +=  - Inertia_I * LocalDeltRotaDisp[1] * kn / equiv_area;
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
    
    
                  
//#C4 //nonlinear().....        
                  
                  
            /* MSIMSI aixo avans estava dintre de el case 0.
            if(rCurrentProcessInfo[NON_LINEAR_OPTION])
                            {
                             
                              
                                if (this->Id() == 1 && *mpTimeStep == 4 ) KRATOS_WATCH( "MUST BE IMPROVED, THE CALCULATION OF STRESS IN THE NON_LINEAR_OPTION ONLY TAKES INTO ACCOUNT THE INDENTATION, IS IT OKAY??" )
                                  
                                double kn_b = kn/rCurrentProcessInfo[SLOPE_FRACTION_N1];
                                double kn_c = kn/rCurrentProcessInfo[SLOPE_FRACTION_N2];
                                
                                
                                double mCompressionLimit_1 = rCurrentProcessInfo[SLOPE_LIMIT_COEFF_C1]*rCurrentProcessInfo[CONTACT_SIGMA_MAX]*1e6;
                                double mCompressionLimit_2 = rCurrentProcessInfo[SLOPE_LIMIT_COEFF_C2]*rCurrentProcessInfo[CONTACT_SIGMA_MAX]*1e6;
                                
                                double sigma_a = (kn * indentation)/(corrected_area);
                             
                                double sigma_b = mCompressionLimit_1 + kn_b*(indentation/corrected_area - mCompressionLimit_1/kn);
                                
                                
                              
                                if( (indentation >= 0.0) && (sigma_a < mCompressionLimit_1) ) 
                                {
                                    LocalElasticContactForce[2]= kn * indentation;
                                    
                                    if(mContactMeshOption)
                                    
                                    {
                                      //lock_p_weak->GetValue(NON_ELASTIC_STAGE) = 1.0;
                                    }
                                    
                                }
                                
                                else if( (indentation >= 0.0) && (sigma_a >= mCompressionLimit_1) && ( sigma_b < mCompressionLimit_2 ) ) 
                                {
                                    LocalElasticContactForce[2]= mCompressionLimit_1*corrected_area + kn_b*(indentation - mCompressionLimit_1*corrected_area/kn);
                                     
                                    if(mContactMeshOption)
                                    
                                    {
                                    //lock_p_weak->GetValue(NON_ELASTIC_STAGE) = 2.0;
                                    
                                    }
                                }
                                
                                else if ( indentation >= 0.0 ) 
                                {
                                    
                                    LocalElasticContactForce[2]= mCompressionLimit_2*corrected_area + kn_c*(indentation - corrected_area*(mCompressionLimit_1/kn + mCompressionLimit_2/kn_b - mCompressionLimit_1/kn_b));
  
                                }
                                
                                else {LocalElasticContactForce[2]= kn * indentation; }
                            }
            
            
   
   //#C6 : initalizesolutionstep del volume strain stress tensor...  
         /*
          double& Representative_Volume = this->GetGeometry()[0].GetSolutionStepValue(REPRESENTATIVE_VOLUME);
         
          Representative_Volume = 0.0;
          
          for (int i=0; i<3; i++)
          {
          
              for (int j=0; j<3; j++)
              {
                   mStressTensor[i][j] = 0.0;
                   mSymmStressTensor[i][j] = 0.0;
              }
          
           }
        */ //MSIMSI
           

      
                   
