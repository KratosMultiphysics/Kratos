//
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
#include "includes/kratos_flags.h"



namespace Kratos
{
     // using namespace GeometryFunctions;

      SphericParticle::SphericParticle()
        : DiscreteElement(), mParticleId(-1), mSqrtOfRealMass(0)/*, mInitializedVariablesFlag(0)*/ {
          mRadius = 0;
          mSqrtOfRealMass = 0;
        }

      SphericParticle::SphericParticle(IndexType NewId, GeometryType::Pointer pGeometry)
        : DiscreteElement(NewId, pGeometry), mParticleId(NewId), mSqrtOfRealMass(0)/*, mInitializedVariablesFlag(0)*/ {
          mRadius = 0;
          mSqrtOfRealMass = 0;          
        }

      SphericParticle::SphericParticle(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
      : DiscreteElement(NewId, pGeometry, pProperties), mParticleId(NewId), mSqrtOfRealMass(0)/*, mInitializedVariablesFlag(0)*/ {
          mRadius = 0;
          mSqrtOfRealMass = 0;             
      }

      SphericParticle::SphericParticle(IndexType NewId, NodesArrayType const& ThisNodes)
      : DiscreteElement(NewId, ThisNodes), mParticleId(NewId), mSqrtOfRealMass(0)/*, mInitializedVariablesFlag(0)*/ {
          mRadius = 0;
          mSqrtOfRealMass = 0;             
      }

      Element::Pointer SphericParticle::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
      {
           return Element::Pointer(new SphericParticle(NewId, GetGeometry().Create(ThisNodes), pProperties));

      }      

      /// Destructor.
      SphericParticle::~SphericParticle(){}

      void SphericParticle::Initialize()
      {
          KRATOS_TRY 

          mDimension                = 3;
          mRadius                   = GetGeometry()[0].FastGetSolutionStepValue(RADIUS);
          double density            = GetDensity();          
          double& sqrt_of_mass      = GetGeometry()[0].FastGetSolutionStepValue(SQRT_OF_MASS);          
          
          //double mass               = 4.0 / 3.0 * KRATOS_M_PI * density * mRadius * mRadius * mRadius;
          double mass               = 1.33333333333333333 * KRATOS_M_PI * density * mRadius * mRadius * mRadius;
          sqrt_of_mass              = sqrt(mass);          
          mSqrtOfRealMass           = sqrt_of_mass;

          if (this->Is(DEMFlags::HAS_ROTATION) ){
            double moment_of_inertia = 0.4 * mass * mRadius * mRadius;   
            GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA) = moment_of_inertia;
          }
          else{
            array_1d<double, 3>& angular_velocity = GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
            angular_velocity[0] = 0.0;
            angular_velocity[1] = 0.0;
            angular_velocity[2] = 0.0;
          }

          if ( GetGeometry()[0].GetDof(VELOCITY_X).IsFixed() )         GetGeometry()[0].Set(DEMFlags::FIXED_VEL_X,true);
          else                                                          GetGeometry()[0].Set(DEMFlags::FIXED_VEL_X,false);
          if ( GetGeometry()[0].GetDof(VELOCITY_Y).IsFixed() )         GetGeometry()[0].Set(DEMFlags::FIXED_VEL_Y,true);
          else                                                          GetGeometry()[0].Set(DEMFlags::FIXED_VEL_Y,false);
          if ( GetGeometry()[0].GetDof(VELOCITY_Z).IsFixed() )         GetGeometry()[0].Set(DEMFlags::FIXED_VEL_Z,true);
          else                                                          GetGeometry()[0].Set(DEMFlags::FIXED_VEL_Z,false);
          if ( GetGeometry()[0].GetDof(ANGULAR_VELOCITY_X).IsFixed() ) GetGeometry()[0].Set(DEMFlags::FIXED_ANG_VEL_X,true);
          else                                                          GetGeometry()[0].Set(DEMFlags::FIXED_ANG_VEL_X,false);
          if ( GetGeometry()[0].GetDof(ANGULAR_VELOCITY_Y).IsFixed() ) GetGeometry()[0].Set(DEMFlags::FIXED_ANG_VEL_Y,true);
          else                                                          GetGeometry()[0].Set(DEMFlags::FIXED_ANG_VEL_Y,false);
          if ( GetGeometry()[0].GetDof(ANGULAR_VELOCITY_Z).IsFixed() ) GetGeometry()[0].Set(DEMFlags::FIXED_ANG_VEL_Z,true);
          else                                                          GetGeometry()[0].Set(DEMFlags::FIXED_ANG_VEL_Z,false);

          //mDiscontinuumConstitutiveLaw = GetProperties()[DEM_DISCONTINUUM_CONSTITUTIVE_LAW_POINTER]->Clone();
          //mDiscontinuumConstitutiveLaw->Initialize();
          //mDiscontinuumConstitutiveLaw = GetProperties()[DEM_DISCONTINUUM_CONSTITUTIVE_LAW_POINTER]->Clone();
          //mDiscontinuumConstitutiveLaw->Initialize();
          
          CustomInitialize();

          KRATOS_CATCH( "" )
      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void SphericParticle::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
      {
          KRATOS_TRY

          //array_1d<double, 3> contact_force;
          //array_1d<double, 3> contact_moment;
          array_1d<double, 3> additionally_applied_force;
          array_1d<double, 3> additionally_applied_moment;
          array_1d<double, 3> initial_rotation_moment;     
          array_1d<double, 3>& elastic_force = this->GetGeometry()[0].FastGetSolutionStepValue(ELASTIC_FORCES);

          //contact_force.clear();
          mContactForce.clear();
          //contact_moment.clear();
          mContactMoment.clear();
          additionally_applied_force.clear();
          additionally_applied_moment.clear();
          initial_rotation_moment.clear();
          elastic_force.clear();

          bool multi_stage_RHS = false;

          ComputeBallToBallContactForce(/*contact_force, contact_moment, */elastic_force, initial_rotation_moment, rCurrentProcessInfo, multi_stage_RHS );

          //Cfeng,RigidFace
          if( mFemOldNeighbourIds.size() > 0)
          {
            ComputeBallToRigidFaceContactForce(/*contact_force, contact_moment,*/ elastic_force, initial_rotation_moment, rCurrentProcessInfo);
          }
              
          ComputeAdditionalForces(additionally_applied_force, additionally_applied_moment, rCurrentProcessInfo);

          /*rRightHandSideVector[0] = contact_force[0]  + additionally_applied_force[0];
          rRightHandSideVector[1] = contact_force[1]  + additionally_applied_force[1];
          rRightHandSideVector[2] = contact_force[2]  + additionally_applied_force[2];
          rRightHandSideVector[3] = contact_moment[0] + additionally_applied_moment[0];
          rRightHandSideVector[4] = contact_moment[1] + additionally_applied_moment[1];
          rRightHandSideVector[5] = contact_moment[2] + additionally_applied_moment[2];*/
          
          rRightHandSideVector[0] = mContactForce[0]  + additionally_applied_force[0];
          rRightHandSideVector[1] = mContactForce[1]  + additionally_applied_force[1];
          rRightHandSideVector[2] = mContactForce[2]  + additionally_applied_force[2];
          rRightHandSideVector[3] = mContactMoment[0] + additionally_applied_moment[0];
          rRightHandSideVector[4] = mContactMoment[1] + additionally_applied_moment[1];
          rRightHandSideVector[5] = mContactMoment[2] + additionally_applied_moment[2];
          
          array_1d<double,3>& total_forces  = this->GetGeometry()[0].FastGetSolutionStepValue(TOTAL_FORCES);
          array_1d<double,3>& total_moment = this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT);

          for (int i = 0; i < 3; i++){
              total_forces[i] = rRightHandSideVector[i];
              total_moment[i] = rRightHandSideVector[3 + i];
          }	  

          KRATOS_CATCH( "" )
      }

      
      void SphericParticle::FirstCalculateRightHandSide(ProcessInfo& rCurrentProcessInfo)
      {
          KRATOS_TRY

          array_1d<double, 3> initial_rotation_moment;     
          array_1d<double, 3>& elastic_force = this->GetGeometry()[0].FastGetSolutionStepValue(ELASTIC_FORCES);

          mContactForce.clear();
          mContactMoment.clear();
          initial_rotation_moment.clear();
          elastic_force.clear();
          
          bool multi_stage_RHS = true;

          ComputeBallToBallContactForce(elastic_force, initial_rotation_moment, rCurrentProcessInfo, multi_stage_RHS );

          
          
          //ConditionWeakVectorType& rNeighbours    = this->GetValue(NEIGHBOUR_RIGID_FACES);
          std::vector<DEMWall*>& rNeighbours = this->mNeighbourRigidFaces;
        
          //for(ConditionWeakIteratorType ineighbour = rNeighbours.begin(); ineighbour != rNeighbours.end(); ineighbour++)
          for (unsigned int i=0; i<rNeighbours.size(); i++) 
          {             
              //Vector& neighbour_rigid_faces_elastic_contact_force = this->GetValue(NEIGHBOUR_RIGID_FACES_ELASTIC_CONTACT_FORCE);
              //Vector& neighbour_rigid_faces_total_contact_force = this->GetValue(NEIGHBOUR_RIGID_FACES_TOTAL_CONTACT_FORCE);
              std::vector<double>& neighbour_rigid_faces_elastic_contact_force = this->mNeighbourRigidFacesElasticContactForce;
              std::vector<double>& neighbour_rigid_faces_total_contact_force = this->mNeighbourRigidFacesTotalContactForce;
              
              neighbour_rigid_faces_elastic_contact_force[3 * i + 0] = 0.0;
              neighbour_rigid_faces_elastic_contact_force[3 * i + 1] = 0.0;
              neighbour_rigid_faces_elastic_contact_force[3 * i + 2] = 0.0;

              neighbour_rigid_faces_total_contact_force[3 * i + 0] = 0.0;
              neighbour_rigid_faces_total_contact_force[3 * i + 1] = 0.0;
              neighbour_rigid_faces_total_contact_force[3 * i + 2] = 0.0;              
          }
          
          
          //Cfeng,RigidFace
          if( mFemOldNeighbourIds.size() > 0)
          {
            ComputeBallToRigidFaceContactForce(elastic_force, initial_rotation_moment, rCurrentProcessInfo);
          }

          KRATOS_CATCH( "" )
      }
      
      void SphericParticle::CollectCalculateRightHandSide(ProcessInfo& rCurrentProcessInfo)
      {
          KRATOS_TRY  
          
          //LOOP OVER NEIGHBOURS BEGINS:          
          for( unsigned int i = 0; i < mNeighbourElements.size(); i++)           
          {
              SphericParticle* neighbour_iterator = mNeighbourElements[i];

              if( this->Is(NEW_ENTITY) && neighbour_iterator->Is(NEW_ENTITY)) continue;
              if( this->Id() < neighbour_iterator->Id() ) continue;
              
              for( unsigned int j = 0; j < neighbour_iterator->mNeighbourElements.size(); j++)  //loop to find the neighbour of the neighbours which is me
              {
                  SphericParticle* is_that_me = neighbour_iterator->mNeighbourElements[j];
                  if ( is_that_me->Id() == this->Id() ){ 
                    mOldNeighbourElasticContactForces[i][0] = -1.0 * neighbour_iterator->mOldNeighbourElasticContactForces[j][0];
                    mOldNeighbourElasticContactForces[i][1] = -1.0 * neighbour_iterator->mOldNeighbourElasticContactForces[j][1];
                    mOldNeighbourElasticContactForces[i][2] = -1.0 * neighbour_iterator->mOldNeighbourElasticContactForces[j][2];

                    mOldNeighbourTotalContactForces[i][0] = -1.0 * neighbour_iterator->mOldNeighbourTotalContactForces[j][0];
                    mOldNeighbourTotalContactForces[i][1] = -1.0 * neighbour_iterator->mOldNeighbourTotalContactForces[j][1];
                    mOldNeighbourTotalContactForces[i][2] = -1.0 * neighbour_iterator->mOldNeighbourTotalContactForces[j][2];      

                    mContactForce[0] += mOldNeighbourTotalContactForces[i][0];
                    mContactForce[1] += mOldNeighbourTotalContactForces[i][1];
                    mContactForce[2] += mOldNeighbourTotalContactForces[i][2];
                    
                    array_1d<double, 3>& rElasticForce = this->GetGeometry()[0].FastGetSolutionStepValue(ELASTIC_FORCES);
                    rElasticForce[0] += mOldNeighbourElasticContactForces[i][0];
                    rElasticForce[1] += mOldNeighbourElasticContactForces[i][1];
                    rElasticForce[2] += mOldNeighbourElasticContactForces[i][2];                            
                    
                    break;
                  }
                  //if(j == neighbour_iterator->mNeighbourElements.size()-1) KRATOS_WATCH("ERRORRRRRRRRRRRRRRRRR")
              }              
              
          }
          

          KRATOS_CATCH( "" )
      }
      
      void SphericParticle::FinalCalculateRightHandSide(ProcessInfo& rCurrentProcessInfo)
      {
          KRATOS_TRY
          
          
          if( this->Is(DEMFlags::HAS_ROTATION) ) {      
                          
                //LOOP OVER NEIGHBOURS BEGINS:          
                for( unsigned int i = 0; i < mNeighbourElements.size(); i++)           
                {
                    SphericParticle* neighbour_iterator = mNeighbourElements[i];

                    if( this->Is(NEW_ENTITY) && neighbour_iterator->Is(NEW_ENTITY)) continue;                                            
                       
                          
                    array_1d<double, 3> other_to_me_vect    = this->GetGeometry()[0].Coordinates() - neighbour_iterator->GetGeometry()[0].Coordinates();
                    double distance                         = sqrt(other_to_me_vect[0] * other_to_me_vect[0] + other_to_me_vect[1] * other_to_me_vect[1] + other_to_me_vect[2] * other_to_me_vect[2]);
                    double inv_distance                     = 1.0/distance;
                    double other_to_me_vect_unitary[3]      = {0.0};
                    other_to_me_vect_unitary[0]             = other_to_me_vect[0] * inv_distance;
                    other_to_me_vect_unitary[1]             = other_to_me_vect[1] * inv_distance;
                    other_to_me_vect_unitary[2]             = other_to_me_vect[2] * inv_distance;


                    double projection_to_local_axis2 = mOldNeighbourElasticContactForces[i][0] * other_to_me_vect_unitary[0] 
                                                     + mOldNeighbourElasticContactForces[i][1] * other_to_me_vect_unitary[1]
                                                     + mOldNeighbourElasticContactForces[i][2] * other_to_me_vect_unitary[2];
                  
                    
                    double RotaAcc[3]                      = {0.0};
                    const array_1d<double, 3> ang_vel = this->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
                    double dt = rCurrentProcessInfo[DELTA_TIME];
                    double dt_i = 1.0 / dt;
                    RotaAcc[0]                         = ang_vel[0] * dt_i;
                    RotaAcc[1]                         = ang_vel[1] * dt_i;
                    RotaAcc[2]                         = ang_vel[2] * dt_i;

                    const double moment_of_inertia         = this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA);
                    array_1d<double, 3> initial_rotation_moment; 
                    initial_rotation_moment[0]               = RotaAcc[0] * moment_of_inertia;
                    initial_rotation_moment[1]               = RotaAcc[1] * moment_of_inertia;
                    initial_rotation_moment[2]               = RotaAcc[2] * moment_of_inertia;
                    if (this->Is(DEMFlags::HAS_ROTATION) ) {
                        ComputeMoments(projection_to_local_axis2, mOldNeighbourElasticContactForces[i], initial_rotation_moment, other_to_me_vect_unitary, neighbour_iterator);
                    }
              
                }
          } //if( this->Is(DEMFlags::HAS_ROTATION) ) 

          array_1d<double, 3> additionally_applied_force;
          array_1d<double, 3> additionally_applied_moment;
          additionally_applied_force.clear();
          additionally_applied_moment.clear();
          
          ComputeAdditionalForces(additionally_applied_force, additionally_applied_moment, rCurrentProcessInfo);
          
          array_1d<double,3>& total_forces  = this->GetGeometry()[0].FastGetSolutionStepValue(TOTAL_FORCES);
          array_1d<double,3>& total_moment = this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT);  
                  
          total_forces[0] = mContactForce[0]  + additionally_applied_force[0];
          total_forces[1] = mContactForce[1]  + additionally_applied_force[1];
          total_forces[2] = mContactForce[2]  + additionally_applied_force[2];
          total_moment[0] = mContactMoment[0] + additionally_applied_moment[0];
          total_moment[1] = mContactMoment[1] + additionally_applied_moment[1];
          total_moment[2] = mContactMoment[2] + additionally_applied_moment[2];
          
                  	  

          KRATOS_CATCH( "" )
      }



      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void SphericParticle::CalculateMaxIndentation(double& rCurrentMaxIndentation, const double& rTolerance)
      {

          if (rCurrentMaxIndentation > rTolerance){
              //ParticleWeakVectorType& rNeighbours = this->GetValue(NEIGHBOUR_ELEMENTS);
              //double& radius                      = GetGeometry()[0].FastGetSolutionStepValue(RADIUS); // cannot use the member variable mRadius because it may not be initialized
              rCurrentMaxIndentation              = 0.0;

              //for (ParticleWeakIteratorType i = rNeighbours.begin(); i != rNeighbours.end(); i++){
              for (unsigned int j = 0; j < mNeighbourElements.size(); j++) {
                  SphericParticle* i = mNeighbourElements[j];
                  
                  array_1d<double, 3> other_to_me_vect  = this->GetGeometry()[0].Coordinates() - i->GetGeometry()[0].Coordinates();
                  //double &other_radius                  = i->GetGeometry()[0].FastGetSolutionStepValue(RADIUS);
                  double other_radius                  = i->GetRadius();
                  double distance                       = sqrt(other_to_me_vect[0] * other_to_me_vect[0] +
                                                               other_to_me_vect[1] * other_to_me_vect[1] +
                                                               other_to_me_vect[2] * other_to_me_vect[2]);
                  double radius_sum                     = mRadius + other_radius;
                  double indentation                    = radius_sum - distance;

                  rCurrentMaxIndentation = (indentation > rCurrentMaxIndentation) ? indentation : rCurrentMaxIndentation;

              }

          }

      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void SphericParticle::CalculateKineticEnergy(double& rKineticEnergy)
      {
          const array_1d<double, 3>& vel    = this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
          const array_1d<double, 3> ang_vel = this->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
          const double moment_of_inertia         = this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA);
          double square_of_celerity         = vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2];
          double square_of_angular_celerity = ang_vel[0] * ang_vel[0] + ang_vel[1] * ang_vel[1] + ang_vel[2] * ang_vel[2];          

          rKineticEnergy = 0.5 * (mSqrtOfRealMass * mSqrtOfRealMass * square_of_celerity + moment_of_inertia * square_of_angular_celerity);
      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void SphericParticle::CalculateElasticEnergyOfContacts(double& rElasticEnergy) // Calculates the elastic energy stored in the sum of all the contacts shared by the particle and all its neighbours
      {
          //ParticleWeakVectorType& rNeighbours         = this->GetValue(NEIGHBOUR_ELEMENTS);
          double added_potential_energy_of_contacts   = 0.0;
          size_t i_neighbour_count                    = 0;
          double myYoung = GetYoung();
          double myPoisson = GetPoisson();

          
          std::vector<int> mTempNeighboursIds;
          std::vector<array_1d<double, 3> > mTempNeighbourElasticContactForces;
          std::vector<array_1d<double, 3> > mTempNeighbourTotalContactForces;
          ComputeNewNeighboursHistoricalData(mTempNeighboursIds, mTempNeighbourElasticContactForces, mTempNeighbourTotalContactForces);

          //for (ParticleWeakIteratorType neighbour_iterator = rNeighbours.begin(); neighbour_iterator != rNeighbours.end(); neighbour_iterator++){
          for( unsigned int i = 0; i < mNeighbourElements.size(); i++) {
              SphericParticle* neighbour_iterator = mNeighbourElements[i];
              
              //const double &other_radius              = neighbour_iterator->GetGeometry()[0].FastGetSolutionStepValue(RADIUS);
              const double &other_radius              = neighbour_iterator->GetRadius();
              double radius_sum                       = mRadius + other_radius;
              double radius_sum_i                     = 1.0 / radius_sum;
              double equiv_radius                     = 2.0 * mRadius * other_radius * radius_sum_i;
              double equiv_area                       = 0.25 * KRATOS_M_PI * equiv_radius * equiv_radius; // 0.25 is becouse we take only the half of the equivalent radius, corresponding to the case of one ball with radius Requivalent and other = radius 0.
              double equiv_young;
              double equiv_poisson;
              double kn;
              double kt;


              // Getting neighbour properties
              //const double &other_young           = neighbour_iterator->GetGeometry()[0].FastGetSolutionStepValue(YOUNG_MODULUS);
              //const double &other_poisson         = neighbour_iterator->GetGeometry()[0].FastGetSolutionStepValue(POISSON_RATIO);
              const double other_young           = neighbour_iterator->GetYoung();
              const double other_poisson         = neighbour_iterator->GetPoisson();

              equiv_young                         = 2.0 * myYoung * other_young / (myYoung + other_young);
              equiv_poisson                       = 2.0 * myPoisson * other_poisson / (myPoisson + other_poisson);
          

       
              kn                                  = equiv_young * equiv_area * radius_sum_i; //KRATOS_M_PI * 0.5 * equiv_young * equiv_radius; //M: CANET FORMULA
              kt                                  = kn / (2.0 + equiv_poisson + equiv_poisson);
           

              // Normal contribution

              double aux_power_of_contact_i_normal_force;

              switch (mElasticityType){ //  0 ---linear compression & tension ; 1 --- Hertzian (non-linear compression, linear tension)
                   case 0:

                       aux_power_of_contact_i_normal_force = mOldNeighbourElasticContactForces[i_neighbour_count][2] * mOldNeighbourElasticContactForces[i_neighbour_count][2];
                       added_potential_energy_of_contacts  += 0.5 * aux_power_of_contact_i_normal_force / kn;

                   break;

                   case 1:
                        aux_power_of_contact_i_normal_force = pow(fabs(mOldNeighbourElasticContactForces[i_neighbour_count][2]), 5 / 3); //error: substitute divisions by known result!!!
                        added_potential_energy_of_contacts  += 0.4 * aux_power_of_contact_i_normal_force / pow(kn, 2 / 3); //error: substitute divisions by known result!!!

                   break;

                   default:

                       aux_power_of_contact_i_normal_force = mOldNeighbourElasticContactForces[i_neighbour_count][2] * mOldNeighbourElasticContactForces[i_neighbour_count][2];
                       added_potential_energy_of_contacts  += 0.5 * aux_power_of_contact_i_normal_force / kn;

                  break;

               }//switch

              // Tangential Contribution

              double aux_power_of_contact_i_tang_force = mOldNeighbourElasticContactForces[i_neighbour_count][0] * mOldNeighbourElasticContactForces[i_neighbour_count][0] + mOldNeighbourElasticContactForces[i_neighbour_count][1] * mOldNeighbourElasticContactForces[i_neighbour_count][1];
              added_potential_energy_of_contacts       += 0.5 * aux_power_of_contact_i_tang_force / kt;

              i_neighbour_count ++;

          }

          rElasticEnergy = added_potential_energy_of_contacts;
      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void SphericParticle::CalculateMomentum(array_1d<double, 3>& rMomentum)
      {
          const array_1d<double, 3>& vel = this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
          rMomentum = mSqrtOfRealMass * mSqrtOfRealMass * vel;
      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void SphericParticle::CalculateLocalAngularMomentum(array_1d<double, 3>& rAngularMomentum)
      {
          const array_1d<double, 3> ang_vel  = this->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
          const double moment_of_inertia         = this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA);
          rAngularMomentum = moment_of_inertia * ang_vel;
      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

     void SphericParticle::ComputeNewNeighboursHistoricalData(std::vector<int>& mTempNeighboursIds, std::vector<array_1d<double, 3> >& mTempNeighbourElasticContactForces,
                                                       std::vector<array_1d<double, 3> >& mTempNeighbourTotalContactForces)
     {
       unsigned int new_size                = mNeighbourElements.size();
       
       mTempNeighboursIds.resize(new_size);       
       
       mTempNeighbourElasticContactForces.resize(new_size);
       mTempNeighbourTotalContactForces.resize(new_size);
       array_1d<double, 3> vector_of_zeros;
       vector_of_zeros[0]                   = 0.0;
       vector_of_zeros[1]                   = 0.0;
       vector_of_zeros[2]                   = 0.0;

       for (unsigned int i = 0; i < mNeighbourElements.size(); i++) {
           SphericParticle* i_neighbour = mNeighbourElements[i];

           mTempNeighboursIds[i] = static_cast<int>(i_neighbour->Id());
           mTempNeighbourElasticContactForces[i] = vector_of_zeros;
           mTempNeighbourTotalContactForces[i] = vector_of_zeros;

           for (unsigned int j = 0; j != mOldNeighbourIds.size(); j++){

               if (static_cast<int>(i_neighbour->Id()) == mOldNeighbourIds[j]){
                   mTempNeighbourElasticContactForces[i] = mOldNeighbourElasticContactForces[j];
                   mTempNeighbourTotalContactForces[i]   = mOldNeighbourTotalContactForces[j];
                   break;
               }

            }

        }

        mOldNeighbourIds.swap(mTempNeighboursIds);
        mOldNeighbourElasticContactForces.swap(mTempNeighbourElasticContactForces);
        mOldNeighbourTotalContactForces.swap(mTempNeighbourTotalContactForces);

      }

     void SphericParticle::ComputeNewRigidFaceNeighboursHistoricalData()
     {

       //ConditionWeakVectorType& rNeighbours  = this->GetValue(NEIGHBOUR_RIGID_FACES);
       std::vector<DEMWall*>& rNeighbours = this->mNeighbourRigidFaces;
       unsigned int new_size                = rNeighbours.size();
       std::vector<int> temp_neighbours_ids(new_size); //these two temporal vectors are very small, saving them as a member of the particle loses time.
       std::vector<array_1d<double, 3> > temp_neighbours_contact_forces(new_size);       
       
       array_1d<double, 3> vector_of_zeros;
       vector_of_zeros[0]                   = 0.0;
       vector_of_zeros[1]                   = 0.0;
       vector_of_zeros[2]                   = 0.0;

       //for (ConditionWeakIteratorType i = rNeighbours.begin(); i != rNeighbours.end(); i++){
       for (unsigned int i=0; i<rNeighbours.size(); i++){    

           temp_neighbours_ids[i] = static_cast<int>(rNeighbours[i]->Id());
           temp_neighbours_contact_forces[i] = vector_of_zeros;

           for (unsigned int j = 0; j != mFemOldNeighbourIds.size(); j++) {

               if (static_cast<int>(rNeighbours[i]->Id()) == mFemOldNeighbourIds[j]) {
                   temp_neighbours_contact_forces[i] = mFemOldNeighbourContactForces[j];
                   break;
               }
            }
        }

        mFemOldNeighbourIds.swap(temp_neighbours_ids);
        mFemOldNeighbourContactForces.swap(temp_neighbours_contact_forces);

      }


      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void SphericParticle::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo){}

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void SphericParticle::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
      {
          rMassMatrix(0,0) = mSqrtOfRealMass * mSqrtOfRealMass;
      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void SphericParticle::EvaluateDeltaDisplacement(double DeltDisp[3],
                                                      double RelVel[3],                                                      
                                                      double LocalCoordSystem[3][3],
                                                      double OldLocalCoordSystem[3][3],
                                                      array_1d<double, 3>& other_to_me_vect,
                                                      const array_1d<double, 3>& vel,
                                                      const array_1d<double, 3>& delta_displ,
                                                      //ParticleWeakIteratorType neighbour_iterator,
                                                      SphericParticle* neighbour_iterator,
                                                      double& distance)
      {
          // FORMING LOCAL CORDINATES

          //Notes: Since we will normally inherit the mesh from GiD, we respect the global system X,Y,Z [0],[1],[2]
          //In the local coordinates we will define the normal direction of the contact as the [2] component!!!!!
          //the way the normal direction is defined (other_to_me_vect) compression is positive!!!
          GeometryFunctions::ComputeContactLocalCoordSystem(other_to_me_vect, distance, LocalCoordSystem); //new Local Coord System (normalizes other_to_me_vect)

          // FORMING OLD LOCAL CORDINATES
          const array_1d<double,3> old_coord_target       = this->GetGeometry()[0].Coordinates() - delta_displ;
          const array_1d<double, 3 > & other_delta_displ  = neighbour_iterator->GetGeometry()[0].FastGetSolutionStepValue(DELTA_DISPLACEMENT);
          const array_1d<double,3> old_coord_neigh        = neighbour_iterator->GetGeometry()[0].Coordinates() - other_delta_displ;
          
          array_1d<double,3> Old_other_to_me_vect = old_coord_target - old_coord_neigh;                    

          const double old_distance = sqrt(Old_other_to_me_vect[0]*Old_other_to_me_vect[0] + Old_other_to_me_vect[1]*Old_other_to_me_vect[1] + Old_other_to_me_vect[2]*Old_other_to_me_vect[2]);
          
          GeometryFunctions::ComputeContactLocalCoordSystem(Old_other_to_me_vect, old_distance, OldLocalCoordSystem); //Old Local Coord System

          // VELOCITIES AND DISPLACEMENTS
          const array_1d<double, 3 >& other_vel          = neighbour_iterator->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
          
          RelVel[0] = (vel[0] - other_vel[0]);
          RelVel[1] = (vel[1] - other_vel[1]);
          RelVel[2] = (vel[2] - other_vel[2]);

          //DeltDisp in global cordinates
          DeltDisp[0] = (delta_displ[0] - other_delta_displ[0]);
          DeltDisp[1] = (delta_displ[1] - other_delta_displ[1]);
          DeltDisp[2] = (delta_displ[2] - other_delta_displ[2]);
      }

      void SphericParticle::DisplacementDueToRotation(double DeltDisp[3],
                                                      double OldLocalCoordSystem[3][3],
                                                      const double& other_radius,
                                                      const double& dt,
                                                      const array_1d<double, 3>& ang_vel,
                                                      SphericParticle* neighbour_iterator)
      {
          array_1d<double, 3> other_ang_vel     = neighbour_iterator->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
          double velA[3]                        = {0.0};
          double velB[3]                        = {0.0};
          double dRotaDisp[3]                   = {0.0};

          double VelTemp[3]                     = {      ang_vel[0],       ang_vel[1],       ang_vel[2]};
          double OtherVelTemp[3]                = {other_ang_vel[0], other_ang_vel[1], other_ang_vel[2]};
          GeometryFunctions::CrossProduct(     VelTemp, OldLocalCoordSystem[2], velA); //it was Local Coordinate system, now we do OLD.
          GeometryFunctions::CrossProduct(OtherVelTemp, OldLocalCoordSystem[2], velB);

          dRotaDisp[0] = - velA[0] * mRadius - velB[0] * other_radius;
          dRotaDisp[1] = - velA[1] * mRadius - velB[1] * other_radius;
          dRotaDisp[2] = - velA[2] * mRadius - velB[2] * other_radius;
          // Contribution of the rotation velocity
          DeltDisp[0] += dRotaDisp[0] * dt;
          DeltDisp[1] += dRotaDisp[1] * dt;
          DeltDisp[2] += dRotaDisp[2] * dt;
      }

      void SphericParticle::ComputeMoments(double normalLocalElasticContactForce,
                                           array_1d<double, 3>& GlobalElasticContactForce,
                                           array_1d<double, 3>& rInitialRotaMoment,
                                           double LocalCoordSystem_2[3],
                                           SphericParticle* neighbour_iterator)
      {
          double MA[3]         = {0.0};      

          GeometryFunctions::CrossProduct(LocalCoordSystem_2, GlobalElasticContactForce, MA);
          
          mContactMoment[0] -= MA[0] * mRadius; 
          mContactMoment[1] -= MA[1] * mRadius;
          mContactMoment[2] -= MA[2] * mRadius;

          // ROLLING FRICTION
          if (this->Is(DEMFlags::HAS_ROLLING_FRICTION) ){    
              double rolling_friction_coeff            = GetRollingFriction() * mRadius;
              const double other_rolling_friction = neighbour_iterator->GetRollingFriction();
              double other_rolling_friction_coeff  = other_rolling_friction * neighbour_iterator->GetRadius();
              double equiv_rolling_friction_coeff  = std::min(rolling_friction_coeff,other_rolling_friction_coeff);

              if (equiv_rolling_friction_coeff != 0.0){
                  double MaxRotaMoment[3]      = {0.0};
                  double CoordSystemMoment1[3] = {0.0};
                  double CoordSystemMoment2[3] = {0.0};
                  double MR[3]                 = {0.0};

                  MaxRotaMoment[0] = rInitialRotaMoment[0] + mContactMoment[0];
                  MaxRotaMoment[1] = rInitialRotaMoment[1] + mContactMoment[1];
                  MaxRotaMoment[2] = rInitialRotaMoment[2] + mContactMoment[2];
                             
                  GeometryFunctions::CrossProduct(LocalCoordSystem_2, MaxRotaMoment, CoordSystemMoment1);
                  double det_coor_sys_moment_i_1 = 1 / sqrt(CoordSystemMoment1[0] * CoordSystemMoment1[0] + CoordSystemMoment1[1] * CoordSystemMoment1[1] + CoordSystemMoment1[2] * CoordSystemMoment1[2]);                             
                  CoordSystemMoment1[0] *= det_coor_sys_moment_i_1;
                  CoordSystemMoment1[1] *= det_coor_sys_moment_i_1;
                  CoordSystemMoment1[2] *= det_coor_sys_moment_i_1;
                                                          
                  GeometryFunctions::CrossProduct(MaxRotaMoment, CoordSystemMoment1, CoordSystemMoment2);
                  double det_coor_sys_moment_i_2 = 1 / sqrt(CoordSystemMoment2[0] * CoordSystemMoment2[0] + CoordSystemMoment2[1] * CoordSystemMoment2[1] + CoordSystemMoment2[2] * CoordSystemMoment2[2]);

                  CoordSystemMoment2[0] *= det_coor_sys_moment_i_2;
                  CoordSystemMoment2[1] *= det_coor_sys_moment_i_2;
                  CoordSystemMoment2[2] *= det_coor_sys_moment_i_2;
                             
                  GeometryFunctions::CrossProduct(CoordSystemMoment2, CoordSystemMoment1, MR);
                  MR[0] *= fabs(normalLocalElasticContactForce);
                  MR[1] *= fabs(normalLocalElasticContactForce);
                  MR[2] *= fabs(normalLocalElasticContactForce);
                             
                  double det_MR = sqrt(MR[0] * MR[0] + MR[1] * MR[1] + MR[2] * MR[2]);
                  double MR_now = det_MR * equiv_rolling_friction_coeff;
                  double MR_max = sqrt(MaxRotaMoment[0] * MaxRotaMoment[0] + MaxRotaMoment[1] * MaxRotaMoment[1] + MaxRotaMoment[2] * MaxRotaMoment[2]);

                  if (MR_max > MR_now){
                       mContactMoment[0] += MR[0] * equiv_rolling_friction_coeff;
                       mContactMoment[1] += MR[1] * equiv_rolling_friction_coeff;
                       mContactMoment[2] += MR[2] * equiv_rolling_friction_coeff;
                   }

                   else {
                       mContactMoment[0] = - rInitialRotaMoment[0];
                       mContactMoment[1] = - rInitialRotaMoment[1];
                       mContactMoment[2] = - rInitialRotaMoment[2];
                   }

               } // if (equiv_rolling_friction_coeff != 0.0)

            } // if (mRotationDampType == 2)            

      }

      void SphericParticle::ComputeBallToBallContactForce(array_1d<double, 3>& rElasticForce,
                                                          array_1d<double, 3>& rInitialRotaMoment,
                                                          ProcessInfo& rCurrentProcessInfo,
                                                          const bool multi_stage_RHS)
      {
          KRATOS_TRY

          // KINEMATICS
          double dt = rCurrentProcessInfo[DELTA_TIME];
          double dt_i = 1 / dt;

          const array_1d<double, 3>& vel         = this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
          const array_1d<double, 3>& delta_displ = this->GetGeometry()[0].FastGetSolutionStepValue(DELTA_DISPLACEMENT);
          const array_1d<double, 3>& ang_vel     = this->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
                    
          if (this->Is(DEMFlags::HAS_ROTATION) && !multi_stage_RHS){
              double RotaAcc[3]                      = {0.0};
              RotaAcc[0]                         = ang_vel[0] * dt_i;
              RotaAcc[1]                         = ang_vel[1] * dt_i;
              RotaAcc[2]                         = ang_vel[2] * dt_i;

              const double moment_of_inertia         = this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA);
              rInitialRotaMoment[0]               = RotaAcc[0] * moment_of_inertia;
              rInitialRotaMoment[1]               = RotaAcc[1] * moment_of_inertia;
              rInitialRotaMoment[2]               = RotaAcc[2] * moment_of_inertia;
          }                    
          
          double kn;
	  double kt;
          double equiv_visco_damp_coeff_normal;
          double equiv_visco_damp_coeff_tangential;
          double equiv_tg_of_fri_ang;
          
          double LocalCoordSystem[3][3]            = {{0.0}, {0.0}, {0.0}};
          double OldLocalCoordSystem[3][3]         = {{0.0}, {0.0}, {0.0}};
          double DeltDisp[3]                       = {0.0};
          double LocalDeltDisp[3]                  = {0.0};
          double RelVel[3]                         = {0.0};
          double LocalRelVel[3]                    = {0.0};
          
          //LOOP OVER NEIGHBOURS BEGINS:
          for( unsigned int i = 0; i < mNeighbourElements.size(); i++)           
          {
              SphericParticle* neighbour_iterator = mNeighbourElements[i];

              if( this->Is(NEW_ENTITY) && neighbour_iterator->Is(NEW_ENTITY)) continue;
              if( multi_stage_RHS  &&  this->Id() > neighbour_iterator->Id()) continue;
              
              
              // BASIC CALCULATIONS              
              array_1d<double, 3> other_to_me_vect    = this->GetGeometry()[0].Coordinates() - neighbour_iterator->GetGeometry()[0].Coordinates();
              const double &other_radius              = neighbour_iterator->GetRadius();
              double distance                         = sqrt(other_to_me_vect[0] * other_to_me_vect[0] + other_to_me_vect[1] * other_to_me_vect[1] + other_to_me_vect[2] * other_to_me_vect[2]);
              
              double radius_sum                       = mRadius + other_radius;
              double indentation                      = radius_sum - distance;
	          

              CalculateEquivalentConstitutiveParameters(other_to_me_vect, other_radius, radius_sum, kn, kt, equiv_visco_damp_coeff_normal, equiv_visco_damp_coeff_tangential, equiv_tg_of_fri_ang, neighbour_iterator);
                            
              DeltDisp[0]      = 0.0; DeltDisp[1]      = 0.0; DeltDisp[2]      = 0.0;
              LocalDeltDisp[0] = 0.0; LocalDeltDisp[1] = 0.0; LocalDeltDisp[2] = 0.0;
              RelVel[0]        = 0.0; RelVel[1]        = 0.0; RelVel[2]        = 0.0;
              LocalRelVel[0]   = 0.0; LocalRelVel[1]   = 0.0; LocalRelVel[2]   = 0.0;                                         
              LocalCoordSystem[0][0]   = 0.0; LocalCoordSystem[0][1]   = 0.0; LocalCoordSystem[0][2]   = 0.0;
              LocalCoordSystem[1][0]   = 0.0; LocalCoordSystem[1][1]   = 0.0; LocalCoordSystem[1][2]   = 0.0;
              LocalCoordSystem[2][0]   = 0.0; LocalCoordSystem[2][1]   = 0.0; LocalCoordSystem[2][2]   = 0.0;
              OldLocalCoordSystem[0][0]= 0.0; OldLocalCoordSystem[0][1]= 0.0; OldLocalCoordSystem[0][2]= 0.0;
              OldLocalCoordSystem[1][0]= 0.0; OldLocalCoordSystem[1][1]= 0.0; OldLocalCoordSystem[1][2]= 0.0;
              OldLocalCoordSystem[2][0]= 0.0; OldLocalCoordSystem[2][1]= 0.0; OldLocalCoordSystem[2][2]= 0.0;

              EvaluateDeltaDisplacement(DeltDisp, RelVel, LocalCoordSystem, OldLocalCoordSystem, other_to_me_vect, vel, delta_displ, neighbour_iterator, distance);

              if (this->Is(DEMFlags::HAS_ROTATION) ){    
                  DisplacementDueToRotation(DeltDisp, OldLocalCoordSystem, other_radius, dt, ang_vel, neighbour_iterator);
              }

              double normal_force                      = 0.0;
              double LocalContactForce[3]              = {0.0};
              double GlobalContactForce[3]             = {0.0};
              double LocalElasticContactForce[3]       = {0.0};
              double GlobalElasticContactForce[3]      = {0.0};
              double ViscoDampingLocalContactForce[3]  = {0.0};
              //double ViscoDampingGlobalContactForce[3] = {0.0};

              GlobalElasticContactForce[0]             = mOldNeighbourElasticContactForces[i][0];
              GlobalElasticContactForce[1]             = mOldNeighbourElasticContactForces[i][1];
              GlobalElasticContactForce[2]             = mOldNeighbourElasticContactForces[i][2];

              GeometryFunctions::VectorGlobal2Local(OldLocalCoordSystem, GlobalElasticContactForce, LocalElasticContactForce); // Here we recover the old local forces projected in the new coordinates in the way they were in the old ones; Now they will be increased if its the necessary
              GeometryFunctions::VectorGlobal2Local(OldLocalCoordSystem, DeltDisp, LocalDeltDisp);
              GeometryFunctions::VectorGlobal2Local(OldLocalCoordSystem, RelVel, LocalRelVel);

              // TRANSLATION FORCES

              bool sliding = false;

              if (indentation > 0.0){
                  NormalForceCalculation(normal_force, kn, indentation);  //ERROR here! The normal force should be added to the local NEW axis!!

                  // TANGENTIAL FORCE. Incremental calculation. An "absolute method" could also be used (YADE?)
                  LocalElasticContactForce[2] = 0.0;
                  TangentialForceCalculation(normal_force, LocalElasticContactForce, LocalDeltDisp, kt, equiv_tg_of_fri_ang, sliding);

                  // VISCODAMPING (applyied locally)
                  CalculateViscoDamping(LocalRelVel, ViscoDampingLocalContactForce, indentation, equiv_visco_damp_coeff_normal, equiv_visco_damp_coeff_tangential, sliding);
              }
              
              // Transforming to global forces and adding up
              AddUpForcesAndProject(OldLocalCoordSystem, LocalCoordSystem, normal_force, LocalContactForce, LocalElasticContactForce, GlobalContactForce, GlobalElasticContactForce, ViscoDampingLocalContactForce, /*ViscoDampingGlobalContactForce, rContactForce, */rElasticForce, i);

              // ROTATION FORCES
              if (this->Is(DEMFlags::HAS_ROTATION) && !multi_stage_RHS){    
                  ComputeMoments(normal_force, mOldNeighbourElasticContactForces[i], rInitialRotaMoment, LocalCoordSystem[2], /*other_radius, rContactMoment,*/ neighbour_iterator);
              }             

          }// for each neighbour

          KRATOS_CATCH("")
      }// ComputeBallToBallContactForce

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************


///////////////////////////////////////////////////Cfeng,RigidFace Contact Force calculation//////////////////////////////
//////////*******************************************************07,Oct,2013*******************************////////////


void SphericParticle::ComputeRigidFaceToMeVelocity(DEMWall* rObj_2, std::size_t ino, 
                             double LocalCoordSystem[3][3], double & DistPToB, array_1d<double, 3 > &other_to_me_vel, int & ContactType)
{
	 KRATOS_TRY
	 
	
	 double Weight[4] = {0.0};
	 
	 //Vector & RF_Pram= this->GetValue(NEIGHBOUR_RIGID_FACES_PRAM);
         std::vector<double>& RF_Pram = this->mNeighbourRigidFacesPram;
	 
	
     int ino1 = ino * 16;
	
	 LocalCoordSystem[0][0] = RF_Pram[ino1 + 0];
	 LocalCoordSystem[0][1] = RF_Pram[ino1 + 1];
	 LocalCoordSystem[0][2] = RF_Pram[ino1 + 2];
	 LocalCoordSystem[1][0] = RF_Pram[ino1 + 3];
	 LocalCoordSystem[1][1] = RF_Pram[ino1 + 4];
	 LocalCoordSystem[1][2] = RF_Pram[ino1 + 5];
	 LocalCoordSystem[2][0] = RF_Pram[ino1 + 6];
	 LocalCoordSystem[2][1] = RF_Pram[ino1 + 7];
	 LocalCoordSystem[2][2] = RF_Pram[ino1 + 8];
	 DistPToB               = RF_Pram[ino1 + 9];
	 Weight[0]              = RF_Pram[ino1 + 10];
	 Weight[1]              = RF_Pram[ino1 + 11];
	 Weight[2]              = RF_Pram[ino1 + 12];
	 Weight[3]              = RF_Pram[ino1 + 13];
	 int iNeighborID        = static_cast<int> (RF_Pram[ino1 + 14]);
     ContactType            = RF_Pram[ino1 + 15];
	 
	 if(iNeighborID == static_cast<int>(rObj_2->Id()) )
	 {
		for(std::size_t inode = 0; inode < rObj_2->GetGeometry().size(); inode++)
		{
		   other_to_me_vel += rObj_2->GetGeometry()(inode)->FastGetSolutionStepValue(VELOCITY) * Weight[inode];
		}		 
	 }
	
	 KRATOS_CATCH("")
}



    void SphericParticle::ComputeBallToRigidFaceContactForce(array_1d<double, 3>& rElasticForce,
                                                            array_1d<double, 3>& rInitialRotaMoment,
                                                            ProcessInfo& rCurrentProcessInfo)
    {
        
      KRATOS_TRY
          
      //ConditionWeakVectorType& rNeighbours    = this->GetValue(NEIGHBOUR_RIGID_FACES);
      std::vector<DEMWall*>& rNeighbours = this->mNeighbourRigidFaces;

      double mTimeStep        = rCurrentProcessInfo[DELTA_TIME];
      double myYoung          = GetYoung();
      double myPoisson        = GetPoisson();
      double area             = KRATOS_M_PI * mRadius * mRadius;

      array_1d<double, 3 > vel = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);

      //for(ConditionWeakIteratorType ineighbour = rNeighbours.begin(); ineighbour != rNeighbours.end(); ineighbour++)
      for (unsigned int i=0; i<rNeighbours.size(); i++)    
      {
        //Condition* p_neighbour_condition = &(*ineighbour);        
        //RigidFace3D* cast_neighbour = dynamic_cast<RigidFace3D*>( p_neighbour_condition );
        DEMWall* cast_neighbour = rNeighbours[i];

        double WallBallFriction = cast_neighbour->mTgOfFrictionAngle;
        
        double LocalElasticContactForce[3]  = {0.0}; 
        double GlobalElasticContactForce[3] = {0.0};        

        array_1d<double, 3 > other_to_me_vel;
        noalias(other_to_me_vel) = ZeroVector(3);

        double LocalCoordSystem[3][3] = {{0.0}, {0.0}, {0.0}};			
        
        array_1d<double, 3> node_coor_array = this->GetGeometry()[0].Coordinates();
        
        double node_coor[3];
        
        node_coor[0]=node_coor_array[0]; //MSIMSI 1 can be optimized.
        node_coor[1]=node_coor_array[1];
        node_coor[2]=node_coor_array[2];
        
        double ini_delta = GetInitialDelta(i);
        double kn_el;


        switch (mElasticityType){ //  0 ---linear compression & tension ; 1 --- Hertzian (non-linear compression, linear tension)
        
            case 0:
                kn_el    = 1.0 * myYoung * area / (2.0 * mRadius - ini_delta);
            break;
            
            case 1:
                kn_el    = (4/3) * (myYoung / 2 * (1 - myPoisson * myPoisson)) * sqrt(mRadius - ini_delta);
            break;
            
            default:
                kn_el    = 1.0 * myYoung * area / (2.0 * mRadius - ini_delta);
            break;
            
        }//switch
        
        double ks_el = kn_el / (2.0 * (1.0 + myPoisson));
        double aux_norm_to_tang = sqrt(ks_el / kn_el);
        double DistPToB = 0.0;

        int ContactType = -1;

        ComputeRigidFaceToMeVelocity(rNeighbours[i], i, LocalCoordSystem, DistPToB, other_to_me_vel, ContactType);

        //MSI: Renew Distance for steps in between searches.
        // The flag ContactType takes value 0 for a plane contact, 1 for line and 2 for point. The optimized function is created for Plane contact, not for other cases yet.
        

        if(ContactType == 0)
        {
         double Coord[4][3] = { {0.0},{0.0},{0.0},{0.0} };

          // Triangle
          
          Coord[0][0] = rNeighbours[i]->GetGeometry()[0].Coordinates()[0];    //MSIMSI 1 can be optimized with vector access.
          Coord[0][1] = rNeighbours[i]->GetGeometry()[0].Coordinates()[1];
          Coord[0][2] = rNeighbours[i]->GetGeometry()[0].Coordinates()[2];

          Coord[1][0] = rNeighbours[i]->GetGeometry()[1].Coordinates()[0];
          Coord[1][1] = rNeighbours[i]->GetGeometry()[1].Coordinates()[1];
          Coord[1][2] = rNeighbours[i]->GetGeometry()[1].Coordinates()[2];

          Coord[2][0] = rNeighbours[i]->GetGeometry()[2].Coordinates()[0];
          Coord[2][1] = rNeighbours[i]->GetGeometry()[2].Coordinates()[1];
          Coord[2][2] = rNeighbours[i]->GetGeometry()[2].Coordinates()[2];

          if(rNeighbours[i]->GetGeometry().size() == 4)
          {
            Coord[3][0] = rNeighbours[i]->GetGeometry()[3].Coordinates()[0];
            Coord[3][1] = rNeighbours[i]->GetGeometry()[3].Coordinates()[1];
            Coord[3][2] = rNeighbours[i]->GetGeometry()[3].Coordinates()[2];
          }
            
            GeometryFunctions::QuickDistanceForAKnownNeighbour(Coord , node_coor, mRadius, DistPToB);

        }
        
        double indentation = -(DistPToB - mRadius) - ini_delta;
        
        double DeltDisp[3] = {0.0};
        double DeltVel [3] = {0.0};

        DeltVel[0] = (vel[0] - other_to_me_vel[0]);
        DeltVel[1] = (vel[1] - other_to_me_vel[1]);
        DeltVel[2] = (vel[2] - other_to_me_vel[2]);

        // For translation movement delt displacement
        DeltDisp[0] = DeltVel[0] * mTimeStep;
        DeltDisp[1] = DeltVel[1] * mTimeStep;
        DeltDisp[2] = DeltVel[2] * mTimeStep;

        if (this->Is(DEMFlags::HAS_ROTATION) )
        {
          double velA[3]   = {0.0};
          double dRotaDisp[3] = {0.0};

          array_1d<double, 3 > AngularVel= GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
          double Vel_Temp[3] = { AngularVel[0], AngularVel[1], AngularVel[2]};
          GeometryFunctions::CrossProduct(Vel_Temp, LocalCoordSystem[2], velA);

          dRotaDisp[0] = -velA[0] * mRadius;
          dRotaDisp[1] = -velA[1] * mRadius;
          dRotaDisp[2] = -velA[2] * mRadius;

          //////contribution of the rotation vel
          DeltDisp[0] += dRotaDisp[0] * mTimeStep;
          DeltDisp[1] += dRotaDisp[1] * mTimeStep;
          DeltDisp[2] += dRotaDisp[2] * mTimeStep;
        }

        double LocalDeltDisp[3] = {0.0};
        GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, DeltDisp, LocalDeltDisp);

        //////120323,for global storage

        GlobalElasticContactForce[0] = mFemOldNeighbourContactForces[i][0];
        GlobalElasticContactForce[1] = mFemOldNeighbourContactForces[i][1];
        GlobalElasticContactForce[2] = mFemOldNeighbourContactForces[i][2];

        GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, GlobalElasticContactForce, LocalElasticContactForce);
        LocalElasticContactForce[0] +=  - ks_el * LocalDeltDisp[0];
        LocalElasticContactForce[1] +=  - ks_el * LocalDeltDisp[1];


        if(indentation > 0.0)
        {
          LocalElasticContactForce[2] =  kn_el * indentation;
        }
        else
        {
          LocalElasticContactForce[2] = 0.0;
        }


        bool If_sliding = false;

        if (-LocalElasticContactForce[2] > 0.0)
        {
          LocalElasticContactForce[0] = 0.0;
          LocalElasticContactForce[1]  = 0.0;
          LocalElasticContactForce[2]  = 0.0;
        }
        else
        {
          double ShearForceMax = LocalElasticContactForce[2] * WallBallFriction;
          double ShearForceNow = sqrt(LocalElasticContactForce[0] * LocalElasticContactForce[0]
          +      LocalElasticContactForce[1] * LocalElasticContactForce[1]);


          //Cfeng: for shear failure
          if(ShearForceMax == 0.0)
          {
            LocalElasticContactForce[0] = 0.0;
            LocalElasticContactForce[1] = 0.0;
          }
          else if(ShearForceNow > ShearForceMax)
          {
            double inv_ShearForceNow = 1.0 / ShearForceNow;
            LocalElasticContactForce[0] = ShearForceMax * inv_ShearForceNow * LocalElasticContactForce[0];
            LocalElasticContactForce[1] = ShearForceMax * inv_ShearForceNow * LocalElasticContactForce[1];

            If_sliding = true;
          }
        }
        
        //---------------------------------------------DAMPING-FORCE-CALCULATION------------------------------------------------------------//

        double ViscoDampingLocalContactForce[3] = {0.0};
        
        if(indentation > 0.0) // Compression , use visc damp
        {			
          double LocalRelVel[3]                   = {0.0};
          GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, DeltVel, LocalRelVel);
          double equiv_visco_damp_coeff_normal = 0.0;
          double equiv_visco_damp_coeff_tangential = 0.0;
          double mass = mSqrtOfRealMass * mSqrtOfRealMass;
          
          if(rCurrentProcessInfo[DEMPACK_OPTION])  //MSIMSI: DEMPACK
          {
            
            equiv_visco_damp_coeff_normal     = rCurrentProcessInfo[DEMPACK_DAMPING]*2.0*sqrt(kn_el/(mass+mass))*mass;//other_sqrt_of_mass*other_sqrt_of_mass))*equiv_mass;   // := 2d0* sqrt ( kn_el*(m1*m2)/(m1+m2) )  //MSIMSI: DEMPACK
            equiv_visco_damp_coeff_tangential = equiv_visco_damp_coeff_normal * aux_norm_to_tang; // dempack no l'utilitza...
          
          }
          
          else
          {
            
            if ( GetLnOfRestitCoeff() > 0.0)
            {
              equiv_visco_damp_coeff_normal = 2.0 * sqrt(mass * kn_el);
            }
            else 
            {
              equiv_visco_damp_coeff_normal = -2.0 * GetLnOfRestitCoeff() * sqrt(mass * kn_el / (GetLnOfRestitCoeff() * GetLnOfRestitCoeff() + KRATOS_M_PI * KRATOS_M_PI));
            }

            equiv_visco_damp_coeff_tangential = equiv_visco_damp_coeff_normal * sqrt(ks_el / kn_el);
            
          }
          
          CalculateViscoDamping(LocalRelVel, ViscoDampingLocalContactForce, indentation, equiv_visco_damp_coeff_normal,
          equiv_visco_damp_coeff_tangential, If_sliding);

        }
        

         double LocalContactForce[3] =                 {0.0};
         double GlobalContactForce[3] =                {0.0};
        
         AddUpFEMForcesAndProject(LocalCoordSystem, LocalContactForce,LocalElasticContactForce,GlobalContactForce,
                                  GlobalElasticContactForce,ViscoDampingLocalContactForce,rElasticForce, i);
        
        if ( this->Is(DEMFlags::HAS_ROTATION) ){    
                  ComputeMoments(LocalElasticContactForce[2], mFemOldNeighbourContactForces[i], rInitialRotaMoment, LocalCoordSystem[2], this); //WARNING: sending itself as the neighbour!!
        }  
        /*if (this->Is(DEMFlags::HAS_ROTATION) )
        {
          double MA[3]                     = {0.0};
          double RotaMoment[3]             = {0.0};

          RotaMoment[0] = mContactMoment[0];
          RotaMoment[1] = mContactMoment[1];
          RotaMoment[2] = mContactMoment[2];

          GeometryFunctions::CrossProduct(LocalCoordSystem[2], GlobalElasticContactForce, MA);

          RotaMoment[0] -= MA[0] * mRadius;
          RotaMoment[1] -= MA[1] * mRadius;
          RotaMoment[2] -= MA[2] * mRadius;

          if (this->Is(DEMFlags::HAS_ROLLING_FRICTION) )
          {  // Rolling friction type              
            
            double rolling_friction_coeff       = GetRollingFriction() * mRadius;

            if (rolling_friction_coeff != 0.0)
            {
              double MaxRotaMoment[3]      = {0.0};
              double CoordSystemMoment1[3] = {0.0};
              double CoordSystemMoment2[3] = {0.0};
              double MR[3]                 = {0.0};

              MaxRotaMoment[0] = rInitialRotaMoment[0] + RotaMoment[0];
              MaxRotaMoment[1] = rInitialRotaMoment[1] + RotaMoment[1];
              MaxRotaMoment[2] = rInitialRotaMoment[2] + RotaMoment[2];

              GeometryFunctions::CrossProduct(LocalCoordSystem[2], MaxRotaMoment, CoordSystemMoment1);
              double det_coor_sys_moment_i_1 = 1 / sqrt(CoordSystemMoment1[0] * CoordSystemMoment1[0] + CoordSystemMoment1[1] * CoordSystemMoment1[1] + CoordSystemMoment1[2] * CoordSystemMoment1[2]);                             
              CoordSystemMoment1[0] *= det_coor_sys_moment_i_1;
              CoordSystemMoment1[1] *= det_coor_sys_moment_i_1;
              CoordSystemMoment1[2] *= det_coor_sys_moment_i_1;

              GeometryFunctions::CrossProduct(MaxRotaMoment, CoordSystemMoment1, CoordSystemMoment2);
              double det_coor_sys_moment_i_2 = 1 / sqrt(CoordSystemMoment2[0] * CoordSystemMoment2[0] + CoordSystemMoment2[1] * CoordSystemMoment2[1] + CoordSystemMoment2[2] * CoordSystemMoment2[2]);

              CoordSystemMoment2[0] *= det_coor_sys_moment_i_2;
              CoordSystemMoment2[1] *= det_coor_sys_moment_i_2;
              CoordSystemMoment2[2] *= det_coor_sys_moment_i_2;

              GeometryFunctions::CrossProduct(CoordSystemMoment2, CoordSystemMoment1, MR);
              MR[0] *= fabs(LocalElasticContactForce[2]);
              MR[1] *= fabs(LocalElasticContactForce[2]);
              MR[2] *= fabs(LocalElasticContactForce[2]);

              double det_MR = sqrt(MR[0] * MR[0] + MR[1] * MR[1] + MR[2] * MR[2]);
              double MR_now = det_MR * rolling_friction_coeff;
              double MR_max = sqrt(MaxRotaMoment[0] * MaxRotaMoment[0] + MaxRotaMoment[1] * MaxRotaMoment[1] + MaxRotaMoment[2] * MaxRotaMoment[2]);

              if (MR_max > MR_now){
              RotaMoment[0]    += MR[0] * rolling_friction_coeff;
              RotaMoment[1]    += MR[1] * rolling_friction_coeff;
              RotaMoment[2]    += MR[2] * rolling_friction_coeff;

              MaxRotaMoment[0] += MR[0] * rolling_friction_coeff;
              MaxRotaMoment[1] += MR[1] * rolling_friction_coeff;
              MaxRotaMoment[2] += MR[2] * rolling_friction_coeff;
              }

              else 
              {
                RotaMoment[0]     = -rInitialRotaMoment[0];
                RotaMoment[1]     = -rInitialRotaMoment[1];
                RotaMoment[2]     = -rInitialRotaMoment[2];
              }

            } // if (rolling_friction_coeff != 0.0)

          } // if (mRollingFrictionOption)

          mContactMoment[0] = RotaMoment[0];
          mContactMoment[1] = RotaMoment[1];
          mContactMoment[2] = RotaMoment[2];

        } //if (mRotationOption)*/


      }

      KRATOS_CATCH("")
        
    }// ComputeBallToRigidFaceContactForce

//////////////////////////////**********Finish***********//////////////////////////////////////////////////



      void SphericParticle::ComputeBallToSurfaceContactForce(array_1d<double, 3>& rContactForce,
                                                             array_1d<double, 3>& rContactMoment,
                                                             array_1d<double, 3>& rInitialRotaMoment,
                                                             int surface_num,
                                                             ProcessInfo& rCurrentProcessInfo)
      {
          KRATOS_TRY

          // CONTACT WITH A PLANE

          double dt                            = rCurrentProcessInfo[DELTA_TIME];

          // INITIALIZATIONS

          const array_1d<double, 3> vel         = this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
          const array_1d<double, 3> delta_displ = this->GetGeometry()[0].FastGetSolutionStepValue(DELTA_DISPLACEMENT);
          const array_1d<double, 3> ang_vel     = this->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
          double InitialRotaMoment[3]           = {0.0};
          double visco_damp_coeff_normal;
          double visco_damp_coeff_tangential;
          //const double &mLnOfRestitCoeff        = this->GetGeometry()[0].FastGetSolutionStepValue(LN_OF_RESTITUTION_COEFF);
          double myYoung = GetYoung();
          double myPoisson = GetPoisson();

          InitialRotaMoment [0] = rInitialRotaMoment [0];
          InitialRotaMoment [1] = rInitialRotaMoment [1];
          InitialRotaMoment [2] = rInitialRotaMoment [2];
             
          // SLIDING

          bool sliding = false;

          // BASIC CALCULATIONS

          array_1d<double,3> surface_normal_dir = rCurrentProcessInfo[SURFACE_NORMAL_DIR_1];
          array_1d<double,3> surface_point_coor = rCurrentProcessInfo[SURFACE_POINT_COOR_1];
          double surface_friction               = rCurrentProcessInfo[SURFACE_FRICTION_1];          

          if (surface_num == 1){
              surface_normal_dir = rCurrentProcessInfo[SURFACE_NORMAL_DIR_2];
              surface_point_coor = rCurrentProcessInfo[SURFACE_POINT_COOR_2];
              surface_friction   = rCurrentProcessInfo[SURFACE_FRICTION_2];
		  }
		  
          if (surface_num == 2){
              surface_normal_dir = rCurrentProcessInfo[SURFACE_NORMAL_DIR_3];
              surface_point_coor = rCurrentProcessInfo[SURFACE_POINT_COOR_3];
              surface_friction   = rCurrentProcessInfo[SURFACE_FRICTION_3];
		  }
		  
          if (surface_num == 3){
              surface_normal_dir = rCurrentProcessInfo[SURFACE_NORMAL_DIR_4];
              surface_point_coor = rCurrentProcessInfo[SURFACE_POINT_COOR_4];
              surface_friction   = rCurrentProcessInfo[SURFACE_FRICTION_4];
		  }

          if (surface_num == 4){
              surface_normal_dir = rCurrentProcessInfo[SURFACE_NORMAL_DIR_5];
              surface_point_coor = rCurrentProcessInfo[SURFACE_POINT_COOR_5];
              surface_friction   = rCurrentProcessInfo[SURFACE_FRICTION_5];
		  }
          
          array_1d<double, 3>& GlobalSurfContactForce = this->GetValue(PARTICLE_SURFACE_CONTACT_FORCES_1);
          if (surface_num == 1){GlobalSurfContactForce = this->GetValue(PARTICLE_SURFACE_CONTACT_FORCES_2);}
          if (surface_num == 2){GlobalSurfContactForce = this->GetValue(PARTICLE_SURFACE_CONTACT_FORCES_3);}
          if (surface_num == 3){GlobalSurfContactForce = this->GetValue(PARTICLE_SURFACE_CONTACT_FORCES_4);}
          if (surface_num == 4){GlobalSurfContactForce = this->GetValue(PARTICLE_SURFACE_CONTACT_FORCES_5);}
          
          array_1d<double, 3> point_coor = this->GetGeometry()[0].Coordinates();

          //Calculate surface equation

          double surface_ecuation[4] = {0.0};
          surface_ecuation[0] = surface_normal_dir[0];
          surface_ecuation[1] = surface_normal_dir[1];
          surface_ecuation[2] = surface_normal_dir[2];
          surface_ecuation[3] = -(surface_normal_dir[0] * surface_point_coor[0] + surface_normal_dir[1] * surface_point_coor[1] + surface_normal_dir[2] * surface_point_coor[2]);

          // Calculate indentation

          double distance =  fabs (surface_ecuation[0] * point_coor[0]       + surface_ecuation[1] * point_coor[1]       + surface_ecuation[2] * point_coor[2] + surface_ecuation[3])
                             / sqrt(surface_ecuation[0] * surface_ecuation[0] + surface_ecuation[1] * surface_ecuation[1] + surface_ecuation[2] * surface_ecuation[2]);

          double indentation = mRadius - distance; //M: Here, Initial_delta is expected to be positive if it is embeding and negative if it's separation.

          if (indentation <= 0.0){
                 GlobalSurfContactForce[0] = 0.0;  // 0: first tangential
                 GlobalSurfContactForce[1] = 0.0;  // 1: second tangential
                 GlobalSurfContactForce[2] = 0.0;  // 2: normal force
          }

          if (indentation > 0.0){
                     // MACRO PARAMETERS

                     double kn;
                     double kt;
                     double equiv_young;

                     switch (mElasticityType){ //  0 ---linear compression & tension ; 1 --- Hertzian (non-linear compression, linear tension)
                
                         case 0:
                             kn = KRATOS_M_PI * 0.5 * myYoung * mRadius; //KRATOS_M_PI * 0.5 * equiv_young * equiv_radius; //M: CANET FORMULA
                             kt = kn / (2.0 * (1.0 + myPoisson));
                        
                         break;

                         case 1:
                             equiv_young = myYoung / (1- myPoisson * myPoisson);
                             kn                 = 1.3333333333 * equiv_young * sqrt(mRadius);
                             kt                 = 2.0 * kn * (1 - myPoisson * myPoisson) / ((2.0 - myPoisson) * (1 + myPoisson));               
                 
                         break;

                         default:
                             kn = KRATOS_M_PI * 0.5 * myYoung * mRadius; //KRATOS_M_PI * 0.5 * equiv_young * equiv_radius; //M: CANET FORMULA
                             kt = kn / (2.0 * (1.0 + myPoisson));

                         break;

                     }//switch


                 // Historical minimun K for the critical time:
                 //if (mCriticalTimeOption){
                 if ( this->Is(DEMFlags::HAS_CRITICAL_TIME) ){
                     double historic = rCurrentProcessInfo[HISTORICAL_MIN_K];

                     if ((kn < historic) || (kt < historic)){
                         historic = std::min(kn, kt);
                     }

                 }   
                 double aux_norm_to_tang = sqrt(kt / kn);
                 double mass = mSqrtOfRealMass * mSqrtOfRealMass;

                 if (GetLnOfRestitCoeff() > 0.0){
                     visco_damp_coeff_normal     = 2 * sqrt(mass * kn);
                     visco_damp_coeff_tangential = visco_damp_coeff_normal * aux_norm_to_tang; // 2 * sqrt(mass * kt);
                 }

          else {
                     visco_damp_coeff_normal     = - 2 * GetLnOfRestitCoeff() * sqrt(mass * kn / (GetLnOfRestitCoeff() * GetLnOfRestitCoeff() + KRATOS_M_PI * KRATOS_M_PI));
                     visco_damp_coeff_tangential = visco_damp_coeff_normal * aux_norm_to_tang; //= -(2 * log(restitution_coeff) * sqrt(mass * kt)) / (sqrt((log(restitution_coeff) * log(restitution_coeff)) + (KRATOS_M_PI * KRATOS_M_PI)));
          }

                 // FORMING LOCAL CORDINATES

                 // Notes: Since we will normally inherit the mesh from GiD, we respect the global system X,Y,Z [0],[1],[2]
                 // In the local coordinates we will define the normal direction of the contact as the [2] component!!!!!
                 // the way the normal direction is defined compression is positive

                 //double NormalDir[3]            = {0.0};
                 array_1d<double, 3> NormalDir;
                 double LocalCoordSystem[3][3]  = {{0.0}, {0.0}, {0.0}};
                 double norm_surface_normal_dir = sqrt(surface_normal_dir[0] * surface_normal_dir[0] + surface_normal_dir[1] * surface_normal_dir[1] + surface_normal_dir[2] * surface_normal_dir[2]);
                 NormalDir[0] = surface_normal_dir[0] / norm_surface_normal_dir;
                 NormalDir[1] = surface_normal_dir[1] / norm_surface_normal_dir;
                 NormalDir[2] = surface_normal_dir[2] / norm_surface_normal_dir;

                 GeometryFunctions::ComputeContactLocalCoordSystem(NormalDir,1.0, LocalCoordSystem); //new Local Coord System

                 // VELOCITIES AND DISPLACEMENTS

                 double DeltDisp[3] = {0.0};
                 double RelVel  [3] = {0.0};

                 RelVel[0] = vel[0];
                 RelVel[1] = vel[1];
                 RelVel[2] = vel[2];

                 // DeltDisp in global cordinates

                 DeltDisp[0] = delta_displ[0];
                 DeltDisp[1] = delta_displ[1];
                 DeltDisp[2] = delta_displ[2];

                 //if (mRotationOption){
                 if (this->Is(DEMFlags::HAS_ROTATION) ){    
                     double velA[3]      = {0.0};
                     double dRotaDisp[3] = {0.0};
                     double Vel_Temp[3] = {ang_vel[0], ang_vel[1], ang_vel[2]};
                     GeometryFunctions::CrossProduct(Vel_Temp, LocalCoordSystem[2], velA);

                     dRotaDisp[0] = -velA[0] * mRadius;
                     dRotaDisp[1] = -velA[1] * mRadius;
                     dRotaDisp[2] = -velA[2] * mRadius;
                     // Contribution of the rotation vel
                     DeltDisp[0] += dRotaDisp[0] * dt;
                     DeltDisp[1] += dRotaDisp[1] * dt;
                     DeltDisp[2] += dRotaDisp[2] * dt;
                 }// if (mRotationOption)

                 double LocalDeltDisp[3]             = {0.0};
                 double LocalElasticContactForce[3]  = {0.0};
                 double GlobalElasticContactForce[3] = {0.0};
                 double LocalRelVel[3]               = {0.0};

                 GlobalElasticContactForce[0] = GlobalSurfContactForce[0];   // GlobalSurfContactForce saved in a container PARTICLE_SURFACE_CONTACT_FORCES
                 GlobalElasticContactForce[1] = GlobalSurfContactForce[1];
                 GlobalElasticContactForce[2] = GlobalSurfContactForce[2];

                 GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, GlobalElasticContactForce, LocalElasticContactForce); //we recover this way the old local forces projected in the new coordinates in the way they were in the old ones; Now they will be increased if its the necessary
                 GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, DeltDisp, LocalDeltDisp);
                 GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, RelVel, LocalRelVel);

                 // FORCES

                 // NORMAL FORCE

                 switch (mElasticityType) //  0---linear comp ; 1 --- Hertzian
                 {
                     case 0:
                         LocalElasticContactForce[2] = kn * indentation;
                         break;

                     case 1:
                         LocalElasticContactForce[2] = kn * pow(indentation, 1.5);
                         break;
                 }

                 // TANGENTIAL FORCE

                 LocalElasticContactForce[0] += - kt * LocalDeltDisp[0];  // 0: first tangential
                 LocalElasticContactForce[1] += - kt * LocalDeltDisp[1];  // 1: second tangential
                 
                 double dyn_friction_angle =  surface_friction * KRATOS_M_PI / 180;
                 double ShearForceNow = sqrt(LocalElasticContactForce[0] * LocalElasticContactForce[0] + LocalElasticContactForce[1] * LocalElasticContactForce[1]);
                 double Frictional_ShearForceMax = tan(dyn_friction_angle) * LocalElasticContactForce[2];                

                 if (Frictional_ShearForceMax < 0.0){
                     Frictional_ShearForceMax = 0.0;
                 }

                 if ((ShearForceNow > Frictional_ShearForceMax) && (ShearForceNow != 0.0)){
                     LocalElasticContactForce[0] = (Frictional_ShearForceMax / ShearForceNow) * LocalElasticContactForce[0];
                     LocalElasticContactForce[1] = (Frictional_ShearForceMax / ShearForceNow) * LocalElasticContactForce[1];
                     sliding = true;
                 }

                 // VISCODAMPING (applyied locally)

                 double ViscoDampingLocalContactForce[3]    = {0.0};

                 if ((mDampType > 0) && ((indentation > 0.0))){

                     if (mDampType == 11 || mDampType == 10){
                         ViscoDampingLocalContactForce[2] = - visco_damp_coeff_normal * LocalRelVel[2];
                     }

                     if (sliding == false && (mDampType == 1 || mDampType == 11)){ //only applied when no sliding to help to the regularized friction law or the spring convergence
                         ViscoDampingLocalContactForce[0] = - visco_damp_coeff_tangential * LocalRelVel[0];
                         ViscoDampingLocalContactForce[1] = - visco_damp_coeff_tangential * LocalRelVel[1];
                     }
                 }

                 // Transforming to global forces and adding up

                 double LocalContactForce[3]              = {0.0};
                 //double ViscoDampingGlobalContactForce[3] = {0.0};
                 double GlobalContactForce[3]             = {0.0};

                 for (unsigned int index = 0; index < 3; index++){
                     LocalContactForce[index] = LocalElasticContactForce[index]  + ViscoDampingLocalContactForce[index];
                 }

                 GeometryFunctions::VectorLocal2Global(LocalCoordSystem, LocalElasticContactForce, GlobalElasticContactForce);
                 //GeometryFunctions::VectorLocal2Global(LocalCoordSystem, ViscoDampingLocalContactForce, ViscoDampingGlobalContactForce);
                 GeometryFunctions::VectorLocal2Global(LocalCoordSystem, LocalContactForce, GlobalContactForce);

                 rContactForce[0] += GlobalContactForce[0];
                 rContactForce[1] += GlobalContactForce[1];
                 rContactForce[2] += GlobalContactForce[2];

                 // Saving contact forces

                    if (surface_num == 0){
                         this->GetValue(PARTICLE_SURFACE_CONTACT_FORCES_1)[0] = GlobalElasticContactForce[0];
                         this->GetValue(PARTICLE_SURFACE_CONTACT_FORCES_1)[1] = GlobalElasticContactForce[1];
                         this->GetValue(PARTICLE_SURFACE_CONTACT_FORCES_1)[2] = GlobalElasticContactForce[2];
			         }
                     if (surface_num == 1){
                         this->GetValue(PARTICLE_SURFACE_CONTACT_FORCES_2)[0] = GlobalElasticContactForce[0];
                         this->GetValue(PARTICLE_SURFACE_CONTACT_FORCES_2)[1] = GlobalElasticContactForce[1];
                         this->GetValue(PARTICLE_SURFACE_CONTACT_FORCES_2)[2] = GlobalElasticContactForce[2];
			         }
                     if (surface_num == 2){
                         this->GetValue(PARTICLE_SURFACE_CONTACT_FORCES_3)[0] = GlobalElasticContactForce[0];
                         this->GetValue(PARTICLE_SURFACE_CONTACT_FORCES_3)[1] = GlobalElasticContactForce[1];
                         this->GetValue(PARTICLE_SURFACE_CONTACT_FORCES_3)[2] = GlobalElasticContactForce[2];
			         }
                     if (surface_num == 3){
                         this->GetValue(PARTICLE_SURFACE_CONTACT_FORCES_4)[0] = GlobalElasticContactForce[0];
                         this->GetValue(PARTICLE_SURFACE_CONTACT_FORCES_4)[1] = GlobalElasticContactForce[1];
                         this->GetValue(PARTICLE_SURFACE_CONTACT_FORCES_4)[2] = GlobalElasticContactForce[2];
			         }
                     if (surface_num == 4){
                         this->GetValue(PARTICLE_SURFACE_CONTACT_FORCES_5)[0] = GlobalElasticContactForce[0];
                         this->GetValue(PARTICLE_SURFACE_CONTACT_FORCES_5)[1] = GlobalElasticContactForce[1];
                         this->GetValue(PARTICLE_SURFACE_CONTACT_FORCES_5)[2] = GlobalElasticContactForce[2];
			         }

                 //if (mRotationOption){
                 if (this->Is(DEMFlags::HAS_ROTATION) ){    
                     double MA[3]                     = {0.0};
                     double RotaMoment[3]             = {0.0};

                     RotaMoment[0] = rContactMoment[0];
                     RotaMoment[1] = rContactMoment[1];
                     RotaMoment[2] = rContactMoment[2];

                     GeometryFunctions::CrossProduct(LocalCoordSystem[2], GlobalElasticContactForce, MA);

                     RotaMoment[0] -= MA[0] * mRadius;
                     RotaMoment[1] -= MA[1] * mRadius;
                     RotaMoment[2] -= MA[2] * mRadius;

                     //if (mRollingFrictionOption){  // Rolling friction type
                     if (this->Is(DEMFlags::HAS_ROLLING_FRICTION) ){
                         //double rolling_friction             = this->GetGeometry()[0].FastGetSolutionStepValue(ROLLING_FRICTION);
                         
                         double rolling_friction_coeff       = GetRollingFriction() * mRadius;

                         if (rolling_friction_coeff != 0.0){
							 double MaxRotaMoment[3]      = {0.0};
							 double CoordSystemMoment1[3] = {0.0};
                             double CoordSystemMoment2[3] = {0.0};
                             double MR[3]                 = {0.0};

                             MaxRotaMoment[0] = InitialRotaMoment[0] + RotaMoment[0];
                             MaxRotaMoment[1] = InitialRotaMoment[1] + RotaMoment[1];
                             MaxRotaMoment[2] = InitialRotaMoment[2] + RotaMoment[2];
                             
                             GeometryFunctions::CrossProduct(LocalCoordSystem[2], MaxRotaMoment, CoordSystemMoment1);
                             double det_coor_sys_moment_i_1 = 1 / sqrt(CoordSystemMoment1[0] * CoordSystemMoment1[0] + CoordSystemMoment1[1] * CoordSystemMoment1[1] + CoordSystemMoment1[2] * CoordSystemMoment1[2]);                             
                             CoordSystemMoment1[0] *= det_coor_sys_moment_i_1;
                             CoordSystemMoment1[1] *= det_coor_sys_moment_i_1;
                             CoordSystemMoment1[2] *= det_coor_sys_moment_i_1;
                                                          
                             GeometryFunctions::CrossProduct(MaxRotaMoment, CoordSystemMoment1, CoordSystemMoment2);
                             double det_coor_sys_moment_i_2 = 1 / sqrt(CoordSystemMoment2[0] * CoordSystemMoment2[0] + CoordSystemMoment2[1] * CoordSystemMoment2[1] + CoordSystemMoment2[2] * CoordSystemMoment2[2]);

                             CoordSystemMoment2[0] *= det_coor_sys_moment_i_2;
                             CoordSystemMoment2[1] *= det_coor_sys_moment_i_2;
                             CoordSystemMoment2[2] *= det_coor_sys_moment_i_2;
                             
                             GeometryFunctions::CrossProduct(CoordSystemMoment2, CoordSystemMoment1, MR);
                             MR[0] *= fabs(LocalElasticContactForce[2]);
                             MR[1] *= fabs(LocalElasticContactForce[2]);
                             MR[2] *= fabs(LocalElasticContactForce[2]);

                             double det_MR = sqrt(MR[0] * MR[0] + MR[1] * MR[1] + MR[2] * MR[2]);
                             double MR_now = det_MR * rolling_friction_coeff;
                             double MR_max = sqrt(MaxRotaMoment[0] * MaxRotaMoment[0] + MaxRotaMoment[1] * MaxRotaMoment[1] + MaxRotaMoment[2] * MaxRotaMoment[2]);

                             if (MR_max > MR_now){
                                 RotaMoment[0]    += MR[0] * rolling_friction_coeff;
                                 RotaMoment[1]    += MR[1] * rolling_friction_coeff;
                                 RotaMoment[2]    += MR[2] * rolling_friction_coeff;

                                 MaxRotaMoment[0] += MR[0] * rolling_friction_coeff;
                                 MaxRotaMoment[1] += MR[1] * rolling_friction_coeff;
                                 MaxRotaMoment[2] += MR[2] * rolling_friction_coeff;
                             }

                             else 
							{
                                 RotaMoment[0]     = -InitialRotaMoment[0];
                                 RotaMoment[1]     = -InitialRotaMoment[1];
                                 RotaMoment[2]     = -InitialRotaMoment[2];
                             }

                         } // if (RollingFrictionCoeff != 0.0)

                     } // if (mRotationDampType == 2)
            
                     rContactMoment[0] = RotaMoment[0];
                     rContactMoment[1] = RotaMoment[1];
                     rContactMoment[2] = RotaMoment[2];

                 } //if (mRotationOption)

             }  //if (indentation >= 0.0)

         KRATOS_CATCH("")
      }//ComputeBallToSurfaceContactForce

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void SphericParticle::ComputeBallToCylinderContactForce(array_1d<double, 3>& rContactForce,
                                                              array_1d<double, 3>& rContactMoment,
                                                              array_1d<double, 3>& rInitialRotaMoment,
                                                              int cylinder_num,
                                                              ProcessInfo& rCurrentProcessInfo)
      {
          KRATOS_TRY 

          double dt             = rCurrentProcessInfo[DELTA_TIME];
          int time_step         = rCurrentProcessInfo[TIME_STEPS];
          
          double myYoung = GetYoung();
          double myPoisson = GetPoisson();
          double myLnOfRestitCoeff = GetLnOfRestitCoeff();

          array_1d<double,3> CylinderAxisDir           = rCurrentProcessInfo[CYLINDER_AXIS_DIR_1];
          array_1d<double,3> InitialBaseCylinderCentre = rCurrentProcessInfo[INITIAL_BASE_CYLINDER_CENTRE_1];
          double CylinderRadius                        = rCurrentProcessInfo[CYLINDER_RADIUS_1];
		  double CylinderVelocity                      = rCurrentProcessInfo[CYLINDER_VELOCITY_1];
          double CylinderAngularVelocity               = rCurrentProcessInfo[CYLINDER_ANGULAR_VELOCITY_1];          
          double CylinderFriction                      = rCurrentProcessInfo[CYLINDER_FRICTION_1];  
          //const double &mLnOfRestitCoeff        = this->GetGeometry()[0].FastGetSolutionStepValue(LN_OF_RESTITUTION_COEFF);

          if (cylinder_num == 1){
			  CylinderAxisDir           = rCurrentProcessInfo[CYLINDER_AXIS_DIR_2];
			  InitialBaseCylinderCentre = rCurrentProcessInfo[INITIAL_BASE_CYLINDER_CENTRE_2];
              CylinderRadius            = rCurrentProcessInfo[CYLINDER_RADIUS_2];
              CylinderVelocity          = rCurrentProcessInfo[CYLINDER_VELOCITY_2];
              CylinderAngularVelocity   = rCurrentProcessInfo[CYLINDER_ANGULAR_VELOCITY_2];
              CylinderFriction          = rCurrentProcessInfo[CYLINDER_FRICTION_2];
		  }
		  
          if (cylinder_num == 2){
			  CylinderAxisDir           = rCurrentProcessInfo[CYLINDER_AXIS_DIR_3];
			  InitialBaseCylinderCentre = rCurrentProcessInfo[INITIAL_BASE_CYLINDER_CENTRE_3];
              CylinderRadius            = rCurrentProcessInfo[CYLINDER_RADIUS_3];
              CylinderVelocity          = rCurrentProcessInfo[CYLINDER_VELOCITY_3];
              CylinderAngularVelocity   = rCurrentProcessInfo[CYLINDER_ANGULAR_VELOCITY_3];
              CylinderFriction          = rCurrentProcessInfo[CYLINDER_FRICTION_3];
		  }
		  
          if (cylinder_num == 3){
			  CylinderAxisDir           = rCurrentProcessInfo[CYLINDER_AXIS_DIR_4];
			  InitialBaseCylinderCentre = rCurrentProcessInfo[INITIAL_BASE_CYLINDER_CENTRE_4];
              CylinderRadius            = rCurrentProcessInfo[CYLINDER_RADIUS_4];
              CylinderVelocity          = rCurrentProcessInfo[CYLINDER_VELOCITY_4];
              CylinderAngularVelocity   = rCurrentProcessInfo[CYLINDER_ANGULAR_VELOCITY_4];
              CylinderFriction          = rCurrentProcessInfo[CYLINDER_FRICTION_4];
		  }

          if (cylinder_num == 4){
			  CylinderAxisDir           = rCurrentProcessInfo[CYLINDER_AXIS_DIR_5];
			  InitialBaseCylinderCentre = rCurrentProcessInfo[INITIAL_BASE_CYLINDER_CENTRE_5];
              CylinderRadius            = rCurrentProcessInfo[CYLINDER_RADIUS_5];
              CylinderVelocity          = rCurrentProcessInfo[CYLINDER_VELOCITY_5];
              CylinderAngularVelocity   = rCurrentProcessInfo[CYLINDER_ANGULAR_VELOCITY_5];
              CylinderFriction          = rCurrentProcessInfo[CYLINDER_FRICTION_5];
		  }
          
          // CONTACT WITH A CYLINDER

             // INITIALIZATIONS

             const array_1d<double, 3> vel         = this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
             const array_1d<double, 3> delta_displ = this->GetGeometry()[0].FastGetSolutionStepValue(DELTA_DISPLACEMENT);
             const array_1d<double, 3> ang_vel     = this->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
             double InitialRotaMoment[3]           = {0.0};
             double visco_damp_coeff_normal;
             double visco_damp_coeff_tangential;

             InitialRotaMoment [0] = rInitialRotaMoment [0];
             InitialRotaMoment [1] = rInitialRotaMoment [1];
             InitialRotaMoment [2] = rInitialRotaMoment [2];

             // SLIDING

             bool sliding = false;

             // BASIC CALCULATIONS

             array_1d<double, 3>& GlobalCylinderContactForce   = this->GetValue(PARTICLE_CYLINDER_CONTACT_FORCES_1);
             if (cylinder_num == 1){GlobalCylinderContactForce = this->GetValue(PARTICLE_CYLINDER_CONTACT_FORCES_2);}
             if (cylinder_num == 2){GlobalCylinderContactForce = this->GetValue(PARTICLE_CYLINDER_CONTACT_FORCES_3);}
             if (cylinder_num == 3){GlobalCylinderContactForce = this->GetValue(PARTICLE_CYLINDER_CONTACT_FORCES_4);}
             if (cylinder_num == 4){GlobalCylinderContactForce = this->GetValue(PARTICLE_CYLINDER_CONTACT_FORCES_5);}
             
             array_1d<double, 3> point_coord                   = this->GetGeometry()[0].Coordinates();

             //Calculate surface equation-> Surface perpendicular to cylinder axis and containing the point centre of the ball

             double cylinder_axis_dir[3] = {0.0};
             cylinder_axis_dir[0] = CylinderAxisDir[0];
             cylinder_axis_dir[1] = CylinderAxisDir[1];
             cylinder_axis_dir[2] = CylinderAxisDir[2];   

             double det_cylinder_axis_dir = sqrt( cylinder_axis_dir[0] + cylinder_axis_dir[0] * cylinder_axis_dir[1] + cylinder_axis_dir[1] * cylinder_axis_dir[2] + cylinder_axis_dir[2] );

             if(det_cylinder_axis_dir != 0.0){
				 cylinder_axis_dir[0] /= det_cylinder_axis_dir;
				 cylinder_axis_dir[1] /= det_cylinder_axis_dir;
				 cylinder_axis_dir[2] /= det_cylinder_axis_dir;
			 }

             double surface_ecuation[4] = {0.0};
             surface_ecuation[0] = cylinder_axis_dir[0];
             surface_ecuation[1] = cylinder_axis_dir[1];
             surface_ecuation[2] = cylinder_axis_dir[2];
             surface_ecuation[3] = -(cylinder_axis_dir[0] * point_coord[0] + cylinder_axis_dir[1] * point_coord[1] + cylinder_axis_dir[2] * point_coord[2]);

             //Calculate intersection point between cylinder axis and surface
             
             double lambda = -(surface_ecuation[0] * InitialBaseCylinderCentre[0] + surface_ecuation[1] * InitialBaseCylinderCentre[1] + surface_ecuation[2] * InitialBaseCylinderCentre[2] + surface_ecuation[3]) / (surface_ecuation[0] * cylinder_axis_dir[0] + surface_ecuation[1] * cylinder_axis_dir[1] + surface_ecuation[2] * cylinder_axis_dir[2]);
             
             double axis_point_coord[3] = {0.0};
             axis_point_coord[0] = InitialBaseCylinderCentre[0] + cylinder_axis_dir[0]*lambda;
             axis_point_coord[1] = InitialBaseCylinderCentre[1] + cylinder_axis_dir[1]*lambda;
             axis_point_coord[2] = InitialBaseCylinderCentre[2] + cylinder_axis_dir[2]*lambda;

             //Calculate normal direction and distance

             double contact_normal_dir[3] = {0.0};
             contact_normal_dir[0] = axis_point_coord[0] - point_coord[0];
             contact_normal_dir[1] = axis_point_coord[1] - point_coord[1];
             contact_normal_dir[2] = axis_point_coord[2] - point_coord[2];

             double det_normal_vect = sqrt( contact_normal_dir[0] * contact_normal_dir[0] + contact_normal_dir[1] * contact_normal_dir[1] + contact_normal_dir[2] * contact_normal_dir[2] );

             if (det_normal_vect == 0.0){
                 GlobalCylinderContactForce[0] = 0.0;  // 0: first tangential
                 GlobalCylinderContactForce[1] = 0.0;  // 1: second tangential
                 GlobalCylinderContactForce[2] = 0.0;  // 2: normal force
             }
             
             if (det_normal_vect > CylinderRadius){
                 contact_normal_dir[0] = -contact_normal_dir[0];
                 contact_normal_dir[1] = -contact_normal_dir[1];
                 contact_normal_dir[2] = -contact_normal_dir[2];
			 }				 

             if(det_normal_vect != 0.0){
		 //double contact_normal_dir_unit[3] = {0.0};
                 array_1d<double, 3> contact_normal_dir_unit;
                 contact_normal_dir_unit[0] = contact_normal_dir[0]/det_normal_vect;
                 contact_normal_dir_unit[1] = contact_normal_dir[1]/det_normal_vect;
                 contact_normal_dir_unit[2] = contact_normal_dir[2]/det_normal_vect;
                 
                 double sphere_contact_point[3] = {0.0};
                 
                 if (det_normal_vect > CylinderRadius){                 
                     sphere_contact_point[0] = - contact_normal_dir[0] + contact_normal_dir_unit[0] * mRadius;
                     sphere_contact_point[1] = - contact_normal_dir[1] + contact_normal_dir_unit[1] * mRadius;
                     sphere_contact_point[2] = - contact_normal_dir[2] + contact_normal_dir_unit[2] * mRadius;
				 }
				 
				 else {
                     sphere_contact_point[0] = - contact_normal_dir[0] - contact_normal_dir_unit[0] * mRadius;
                     sphere_contact_point[1] = - contact_normal_dir[1] - contact_normal_dir_unit[1] * mRadius;
                     sphere_contact_point[2] = - contact_normal_dir[2] - contact_normal_dir_unit[2] * mRadius;
				 }

                 double sphere_distance = sqrt( sphere_contact_point[0] * sphere_contact_point[0] + sphere_contact_point[1] * sphere_contact_point[1] + sphere_contact_point[2] * sphere_contact_point[2] );                
                 
                 // Calculate indentation
                 
                 double indentation = 0.0;

                 if (det_normal_vect <= CylinderRadius && det_normal_vect != 0.0){
                     indentation = sphere_distance - CylinderRadius; //M: Here, Initial_delta is expected to be positive if it is embeding and negative if it's separation.
			     }
             
                 if (det_normal_vect > CylinderRadius && det_normal_vect != 0.0){
                     indentation = CylinderRadius - sphere_distance; //M: Here, Initial_delta is expected to be positive if it is embeding and negative if it's separation.
			     }
			     
                 // Calculate cylinder outlet
                 
                 double outlet_point_coord[3] = {0.0};
                 outlet_point_coord[0] = InitialBaseCylinderCentre[0] + (cylinder_axis_dir[0] * CylinderVelocity * dt * time_step);
                 outlet_point_coord[1] = InitialBaseCylinderCentre[1] + (cylinder_axis_dir[1] * CylinderVelocity * dt * time_step);
                 outlet_point_coord[2] = InitialBaseCylinderCentre[2] + (cylinder_axis_dir[2] * CylinderVelocity * dt * time_step);
                 
                 double cylinder_outlet_surface_ecuation[4] = {0.0};
                 cylinder_outlet_surface_ecuation[0] = cylinder_axis_dir[0];
                 cylinder_outlet_surface_ecuation[1] = cylinder_axis_dir[1];
                 cylinder_outlet_surface_ecuation[2] = cylinder_axis_dir[2];
                 cylinder_outlet_surface_ecuation[3] = -(cylinder_axis_dir[0] * outlet_point_coord[0] + cylinder_axis_dir[1] * outlet_point_coord[1] + cylinder_axis_dir[2] * outlet_point_coord[2]);
                 
                 double position_ball_outlet = cylinder_outlet_surface_ecuation[0] * point_coord[0] + cylinder_outlet_surface_ecuation[1] * point_coord[1] + cylinder_outlet_surface_ecuation[2] * point_coord[2] + cylinder_outlet_surface_ecuation[3];                 			     

                 if (indentation <= 0.0 || position_ball_outlet <= 0.0){
                     GlobalCylinderContactForce[0] = 0.0;  // 0: first tangential
                     GlobalCylinderContactForce[1] = 0.0;  // 1: second tangential
                     GlobalCylinderContactForce[2] = 0.0;  // 2: normal force
                 }

                 if (indentation > 0.0 && position_ball_outlet > 0.0){
                     // MACRO PARAMETERS

                     double kn;
                     double kt;
                     double equiv_young;
                     double effective_radius;

                     switch (mElasticityType){ //  0 ---linear compression & tension ; 1 --- Hertzian (non-linear compression, linear tension)
                
                         case 0:
                             kn = KRATOS_M_PI * 0.5 * myYoung * mRadius; //KRATOS_M_PI * 0.5 * equiv_young * equiv_radius; //M: CANET FORMULA
                             kt = kn / (2.0 * (1.0 + myPoisson));
                        
                         break;

                         case 1:
                             equiv_young      = myYoung / (1- myPoisson * myPoisson);
                             effective_radius = mRadius * CylinderRadius / (CylinderRadius - mRadius);
                             kn                      = 1.3333333333 * equiv_young * sqrt(effective_radius);
                             kt                      = 2.0 * kn * (1 - myPoisson * myPoisson) / ((2.0 - myPoisson) * (1 + myPoisson));               
                 
                         break;

                         default:
                             kn = KRATOS_M_PI * 0.5 * myYoung * mRadius; //KRATOS_M_PI * 0.5 * equiv_young * equiv_radius; //M: CANET FORMULA
                             kt = kn / (2.0 * (1.0 + myPoisson));

                         break;

                     }//switch
                    
        
                     // Historical minimun K for the critical time:
                     //if (mCriticalTimeOption){
                     if ( this->Is(DEMFlags::HAS_CRITICAL_TIME) ){
                         double historic = rCurrentProcessInfo[HISTORICAL_MIN_K];
    
                         if ((kn < historic) || (kt < historic)){
                             historic = std::min(kn, kt);
                         }

                     }
                     double aux_norm_to_tang = sqrt(kt / kn);
                     double mass = mSqrtOfRealMass * mSqrtOfRealMass;

                     if (myLnOfRestitCoeff > 0.0){
                         visco_damp_coeff_normal     = 2 * sqrt(mass * kn);
                         visco_damp_coeff_tangential = visco_damp_coeff_normal * aux_norm_to_tang; // 2 * sqrt(mass * kt);
                     }
    
                     else {
                         visco_damp_coeff_normal     = - 2 * myLnOfRestitCoeff * sqrt(mass * kn / (myLnOfRestitCoeff * myLnOfRestitCoeff + KRATOS_M_PI * KRATOS_M_PI));
                         visco_damp_coeff_tangential = visco_damp_coeff_normal * aux_norm_to_tang; //= -(2 * log(restitution_coeff) * sqrt(mass * kt)) / (sqrt((log(restitution_coeff) * log(restitution_coeff)) + (KRATOS_M_PI * KRATOS_M_PI)));
                     }

                     // FORMING LOCAL CORDINATES

                     // Notes: Since we will normally inherit the mesh from GiD, we respect the global system X,Y,Z [0],[1],[2]
                     // In the local coordinates we will define the normal direction of the contact as the [2] component!!!!!
                     // the way the normal direction is defined compression is positive

                     double LocalCoordSystem[3][3]  = {{0.0}, {0.0}, {0.0}};
                 
                     GeometryFunctions::ComputeContactLocalCoordSystem(contact_normal_dir_unit, 1.0, LocalCoordSystem);

                     // VELOCITIES AND DISPLACEMENTS

                     double DeltDisp[3] = {0.0};
                     double RelVel  [3] = {0.0};

                     RelVel[0] = (vel[0] - CylinderVelocity * cylinder_axis_dir[0]);
                     RelVel[1] = (vel[1] - CylinderVelocity * cylinder_axis_dir[1]);
                     RelVel[2] = (vel[2] - CylinderVelocity * cylinder_axis_dir[2]);

                     // DeltDisp in global cordinates

                     DeltDisp[0] = (delta_displ[0] - CylinderVelocity * cylinder_axis_dir[0] * dt);
                     DeltDisp[1] = (delta_displ[1] - CylinderVelocity * cylinder_axis_dir[1] * dt);
                     DeltDisp[2] = (delta_displ[2] - CylinderVelocity * cylinder_axis_dir[2] * dt);

                     //if (mRotationOption){
                     if (this->Is(DEMFlags::HAS_ROTATION) ){                     
                         double velA[3]      = {0.0};
                         double dRotaDisp[3] = {0.0};
                         double Vel_Temp[3] = {ang_vel[0], ang_vel[1], ang_vel[2]};
                         GeometryFunctions::CrossProduct(Vel_Temp, LocalCoordSystem[2], velA);

                         dRotaDisp[0] = -velA[0] * mRadius;
                         dRotaDisp[1] = -velA[1] * mRadius;
                         dRotaDisp[2] = -velA[2] * mRadius;
                         // Contribution of the rotation vel
                         DeltDisp[0] += dRotaDisp[0] * dt;
                         DeltDisp[1] += dRotaDisp[1] * dt;
                         DeltDisp[2] += dRotaDisp[2] * dt;
                     }// if (mRotationOption)
                     
                     if (CylinderAngularVelocity != 0.0){
						 //double cylinder_angular_velocity[3] = {0.0};
                         array_1d<double, 3> cylinder_angular_velocity;
                         cylinder_angular_velocity[0] = CylinderAngularVelocity * cylinder_axis_dir[0];
                         cylinder_angular_velocity[1] = CylinderAngularVelocity * cylinder_axis_dir[1];
                         cylinder_angular_velocity[2] = CylinderAngularVelocity * cylinder_axis_dir[2];
                         //double vel_cyl[3]      = {0.0};
                         array_1d<double, 3> vel_cyl;
                         GeometryFunctions::CrossProduct(cylinder_angular_velocity, contact_normal_dir_unit, vel_cyl);
                         
                         if (det_normal_vect > CylinderRadius){                 
                             DeltDisp[0] -= vel_cyl[0] * CylinderRadius * dt;
                             DeltDisp[1] -= vel_cyl[1] * CylinderRadius * dt;
                             DeltDisp[2] -= vel_cyl[2] * CylinderRadius * dt;
				         }
				 
				         else {
                             DeltDisp[0] += vel_cyl[0] * CylinderRadius * dt;
                             DeltDisp[1] += vel_cyl[1] * CylinderRadius * dt;
                             DeltDisp[2] += vel_cyl[2] * CylinderRadius * dt;
				         }                         
					 }

                     double LocalDeltDisp[3]             = {0.0};
                     double LocalElasticContactForce[3]  = {0.0};
                     double GlobalElasticContactForce[3] = {0.0};
                     double LocalRelVel[3]               = {0.0};

                     GlobalElasticContactForce[0] = GlobalCylinderContactForce[0];   // GlobalCylinderContactForce saved in a container PARTICLE_CYLINDER_CONTACT_FORCES
                     GlobalElasticContactForce[1] = GlobalCylinderContactForce[1];
                     GlobalElasticContactForce[2] = GlobalCylinderContactForce[2];

                     GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, GlobalElasticContactForce, LocalElasticContactForce); //we recover this way the old local forces projected in the new coordinates in the way they were in the old ones; Now they will be increased if its the necessary
                     GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, DeltDisp, LocalDeltDisp);
                     GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, RelVel, LocalRelVel);

                     // FORCES
    
                     // NORMAL FORCE

                     switch (mElasticityType) //  0---linear comp ; 1 --- Hertzian
                     {
                         case 0:
                             LocalElasticContactForce[2] = kn * indentation;
                             break;

                         case 1:
                             LocalElasticContactForce[2] = kn * pow(indentation, 1.5);
                             break;
                     }

                     // TANGENTIAL FORCE

                     LocalElasticContactForce[0] += - kt * LocalDeltDisp[0];  // 0: first tangential
                     LocalElasticContactForce[1] += - kt * LocalDeltDisp[1];  // 1: second tangential

                     double dyn_friction_angle =  CylinderFriction * KRATOS_M_PI / 180;
                     double ShearForceNow = sqrt(LocalElasticContactForce[0] * LocalElasticContactForce[0] + LocalElasticContactForce[1] * LocalElasticContactForce[1]);
                     double Frictional_ShearForceMax = tan(dyn_friction_angle) * LocalElasticContactForce[2];

                     if (Frictional_ShearForceMax < 0.0){
                         Frictional_ShearForceMax = 0.0;
                     }

                     if ((ShearForceNow > Frictional_ShearForceMax) && (ShearForceNow != 0.0)){
                         LocalElasticContactForce[0] = (Frictional_ShearForceMax / ShearForceNow) * LocalElasticContactForce[0];
                         LocalElasticContactForce[1] = (Frictional_ShearForceMax / ShearForceNow) * LocalElasticContactForce[1];
                         sliding = true;
                     }

                     // VISCODAMPING (applyied locally)

                     double ViscoDampingLocalContactForce[3]    = {0.0};

                     if ((mDampType > 0) && ((indentation > 0.0))){

                         if (mDampType == 11 || mDampType == 10){
                             ViscoDampingLocalContactForce[2] = - visco_damp_coeff_normal * LocalRelVel[2];
                         }

                         if (sliding == false && (mDampType == 1 || mDampType == 11)){ //only applied when no sliding to help to the regularized friction law or the spring convergence
                             ViscoDampingLocalContactForce[0] = - visco_damp_coeff_tangential * LocalRelVel[0];
                             ViscoDampingLocalContactForce[1] = - visco_damp_coeff_tangential * LocalRelVel[1];
                         }
                     }

                     // Transforming to global forces and adding up

                     double LocalContactForce[3]              = {0.0};
                     //double ViscoDampingGlobalContactForce[3] = {0.0};
                     double GlobalContactForce[3]             = {0.0};

                     for (unsigned int index = 0; index < 3; index++){
                         LocalContactForce[index] = LocalElasticContactForce[index]  + ViscoDampingLocalContactForce[index];
                     }

                     GeometryFunctions::VectorLocal2Global(LocalCoordSystem, LocalElasticContactForce, GlobalElasticContactForce);
                     //GeometryFunctions::VectorLocal2Global(LocalCoordSystem, ViscoDampingLocalContactForce, ViscoDampingGlobalContactForce);
                     GeometryFunctions::VectorLocal2Global(LocalCoordSystem, LocalContactForce, GlobalContactForce);

                     rContactForce[0] += GlobalContactForce[0];
                     rContactForce[1] += GlobalContactForce[1];
                     rContactForce[2] += GlobalContactForce[2];

                     // Saving contact forces
                     
                     if (cylinder_num == 0){
                         this->GetValue(PARTICLE_CYLINDER_CONTACT_FORCES_1)[0] = GlobalElasticContactForce[0];
                         this->GetValue(PARTICLE_CYLINDER_CONTACT_FORCES_1)[1] = GlobalElasticContactForce[1];
                         this->GetValue(PARTICLE_CYLINDER_CONTACT_FORCES_1)[2] = GlobalElasticContactForce[2];
			         }
                     if (cylinder_num == 1){
                         this->GetValue(PARTICLE_CYLINDER_CONTACT_FORCES_2)[0] = GlobalElasticContactForce[0];
                         this->GetValue(PARTICLE_CYLINDER_CONTACT_FORCES_2)[1] = GlobalElasticContactForce[1];
                         this->GetValue(PARTICLE_CYLINDER_CONTACT_FORCES_2)[2] = GlobalElasticContactForce[2];
			         }
                     if (cylinder_num == 2){
                         this->GetValue(PARTICLE_CYLINDER_CONTACT_FORCES_3)[0] = GlobalElasticContactForce[0];
                         this->GetValue(PARTICLE_CYLINDER_CONTACT_FORCES_3)[1] = GlobalElasticContactForce[1];
                         this->GetValue(PARTICLE_CYLINDER_CONTACT_FORCES_3)[2] = GlobalElasticContactForce[2];
			         }
                     if (cylinder_num == 3){
                         this->GetValue(PARTICLE_CYLINDER_CONTACT_FORCES_4)[0] = GlobalElasticContactForce[0];
                         this->GetValue(PARTICLE_CYLINDER_CONTACT_FORCES_4)[1] = GlobalElasticContactForce[1];
                         this->GetValue(PARTICLE_CYLINDER_CONTACT_FORCES_4)[2] = GlobalElasticContactForce[2];
			         }
                     if (cylinder_num == 4){
                         this->GetValue(PARTICLE_CYLINDER_CONTACT_FORCES_5)[0] = GlobalElasticContactForce[0];
                         this->GetValue(PARTICLE_CYLINDER_CONTACT_FORCES_5)[1] = GlobalElasticContactForce[1];
                         this->GetValue(PARTICLE_CYLINDER_CONTACT_FORCES_5)[2] = GlobalElasticContactForce[2];
			         }

                     //if (mRotationOption){
                     if (this->Is(DEMFlags::HAS_ROTATION) ){    
                         double MA[3]                     = {0.0};
                         double RotaMoment[3]             = {0.0};
                         
                         RotaMoment[0] = rContactMoment[0];
                         RotaMoment[1] = rContactMoment[1];
                         RotaMoment[2] = rContactMoment[2];                         

                         GeometryFunctions::CrossProduct(LocalCoordSystem[2], GlobalElasticContactForce, MA);

                         RotaMoment[0] -= MA[0] * mRadius;
                         RotaMoment[1] -= MA[1] * mRadius;
                         RotaMoment[2] -= MA[2] * mRadius;

                         //if (mRollingFrictionOption){  // Rolling friction 
                         if (this->Is(DEMFlags::HAS_ROLLING_FRICTION) ){
                             //double rolling_friction             = this->GetGeometry()[0].FastGetSolutionStepValue(ROLLING_FRICTION);
                             double rolling_friction_coeff       = GetRollingFriction() * mRadius;

                             if (rolling_friction_coeff != 0.0){
								 double MaxRotaMoment[3]      = {0.0};
							     double CoordSystemMoment1[3] = {0.0};
                                 double CoordSystemMoment2[3] = {0.0};
                                 double MR[3]                 = {0.0};

                                 MaxRotaMoment[0] = InitialRotaMoment[0] + RotaMoment[0];
                                 MaxRotaMoment[1] = InitialRotaMoment[1] + RotaMoment[1];
                                 MaxRotaMoment[2] = InitialRotaMoment[2] + RotaMoment[2];
                             
                                 GeometryFunctions::CrossProduct(LocalCoordSystem[2], MaxRotaMoment, CoordSystemMoment1);
                                 double det_coor_sys_moment_i_1 = 1 / sqrt(CoordSystemMoment1[0] * CoordSystemMoment1[0] + CoordSystemMoment1[1] * CoordSystemMoment1[1] + CoordSystemMoment1[2] * CoordSystemMoment1[2]);                             
                                 CoordSystemMoment1[0] *= det_coor_sys_moment_i_1;
                                 CoordSystemMoment1[1] *= det_coor_sys_moment_i_1;
                                 CoordSystemMoment1[2] *= det_coor_sys_moment_i_1;
                                                          
                                 GeometryFunctions::CrossProduct(MaxRotaMoment, CoordSystemMoment1, CoordSystemMoment2);
                                 double det_coor_sys_moment_i_2 = 1 / sqrt(CoordSystemMoment2[0] * CoordSystemMoment2[0] + CoordSystemMoment2[1] * CoordSystemMoment2[1] + CoordSystemMoment2[2] * CoordSystemMoment2[2]);

                                 CoordSystemMoment2[0] *= det_coor_sys_moment_i_2;
                                 CoordSystemMoment2[1] *= det_coor_sys_moment_i_2;
                                 CoordSystemMoment2[2] *= det_coor_sys_moment_i_2;
                             
                                 GeometryFunctions::CrossProduct(CoordSystemMoment2, CoordSystemMoment1, MR);
                                 MR[0] *= fabs(LocalElasticContactForce[2]);
                                 MR[1] *= fabs(LocalElasticContactForce[2]);
                                 MR[2] *= fabs(LocalElasticContactForce[2]);

                                 double det_MR = sqrt(MR[0] * MR[0] + MR[1] * MR[1] + MR[2] * MR[2]);
                                 double MR_now = det_MR * rolling_friction_coeff;
                                 double MR_max = sqrt(MaxRotaMoment[0] * MaxRotaMoment[0] + MaxRotaMoment[1] * MaxRotaMoment[1] + MaxRotaMoment[2] * MaxRotaMoment[2]);

                                 if (MR_max > MR_now){
                                     RotaMoment[0]    += MR[0] * rolling_friction_coeff;
                                     RotaMoment[1]    += MR[1] * rolling_friction_coeff;
                                     RotaMoment[2]    += MR[2] * rolling_friction_coeff;
                                 }

                                 else {
                                     RotaMoment[0]     = -InitialRotaMoment[0];
                                     RotaMoment[1]     = -InitialRotaMoment[1];
                                     RotaMoment[2]     = -InitialRotaMoment[2];                            
                                 }

                             } // if (RollingFrictionCoeff != 0.0)

                         } // if (mRotationDampType == 2)
            
                         rContactMoment[0] = RotaMoment[0];
                         rContactMoment[1] = RotaMoment[1];
                         rContactMoment[2] = RotaMoment[2];

                         rInitialRotaMoment [0] = InitialRotaMoment [0];
                         rInitialRotaMoment [1] = InitialRotaMoment [1];
                         rInitialRotaMoment [2] = InitialRotaMoment [2];

                     } //if (mRotationOption)

                 }  //if (indentation >= 0.0)
             
             }  //if(det_normal_vect != 0.0)
             
             KRATOS_CATCH("")
      }//ComputeBallToCylinderContactForce

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void SphericParticle::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo){}

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void SphericParticle::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
      {
          KRATOS_TRY

          ElementalDofList.resize(0);

          for (unsigned int i = 0; i < GetGeometry().size(); i++){
              ElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_X));
              ElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_Y));
              if (GetGeometry().WorkingSpaceDimension() == 3){
                  ElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_Z));
              }

              ElementalDofList.push_back(GetGeometry()[i].pGetDof(ANGULAR_VELOCITY_X));
              ElementalDofList.push_back(GetGeometry()[i].pGetDof(ANGULAR_VELOCITY_Y));
              if (GetGeometry().WorkingSpaceDimension() == 3){
                  ElementalDofList.push_back(GetGeometry()[i].pGetDof(ANGULAR_VELOCITY_Z));
              }

          }

          KRATOS_CATCH("")
      }


      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void SphericParticle::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
      {
          KRATOS_TRY
          
          const ProcessInfo& r_process_info = rCurrentProcessInfo;

          MemberDeclarationFirstStep(r_process_info);
          //mInitializedVariablesFlag = 1;
          this->Set(DEMFlags::HAS_INITIALIZED_VARIABLES);

          KRATOS_CATCH("")

      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void SphericParticle::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo)
      {
        

      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void SphericParticle::CustomInitialize(){}

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void SphericParticle::CalculateEquivalentConstitutiveParameters(array_1d<double, 3>& other_to_me_vect,
                                                                      const double& other_radius,
                                                                      const double& radius_sum,
                                                                      double& kn,
                                                                      double& kt,
                                                                      double& equiv_visco_damp_coeff_normal,
                                                                      double& equiv_visco_damp_coeff_tangential,
                                                                      double& equiv_tg_of_fri_ang,
                                                                      SphericParticle* neighbour_iterator)
      {
          
        double other_sqrt_of_mass = neighbour_iterator->GetSqrtOfRealMass();        
        double radius_sum_i = 1 / radius_sum;
        double equiv_radius = 2 * mRadius * other_radius * radius_sum_i;
        double equiv_area = 0.25 * KRATOS_M_PI * equiv_radius * equiv_radius; // 0.25 is because we take only the half of the equivalent radius, corresponding to the case of one ball with radius Requivalent and other = radius 0.
        double equiv_mass = mSqrtOfRealMass * other_sqrt_of_mass;
         
        double myYoung = GetYoung();
        double myPoisson = GetPoisson();
        double myLnOfRestitCoeff = GetLnOfRestitCoeff();  
        double myTgOfFrictionAngle = GetTgOfFrictionAngle();
        
        double other_young;
        double other_poisson;
        double equiv_ln_of_restit_coeff;  
        double equiv_young;
        double equiv_poisson;
        double other_ln_of_restit_coeff;
          
        if( mFastProperties == neighbour_iterator->GetFastProperties() ) {
            other_young = myYoung;
            
            other_poisson = myPoisson;
            equiv_poisson = myPoisson;
            
            other_ln_of_restit_coeff = myLnOfRestitCoeff;
            equiv_ln_of_restit_coeff = myLnOfRestitCoeff;
            
            equiv_tg_of_fri_ang = myTgOfFrictionAngle;
            
        }
        else{
            other_young = neighbour_iterator->GetYoung();            
            
            other_poisson = neighbour_iterator->GetPoisson();
            equiv_poisson = 2.0 * myPoisson * other_poisson / (myPoisson + other_poisson); 
            
            other_ln_of_restit_coeff = neighbour_iterator->GetLnOfRestitCoeff();
            equiv_ln_of_restit_coeff        = 0.5 * (myLnOfRestitCoeff + other_ln_of_restit_coeff);
            
            const double other_tg_of_fri_angle = neighbour_iterator->GetTgOfFrictionAngle();                        
            equiv_tg_of_fri_ang             = 0.5 * (myTgOfFrictionAngle + other_tg_of_fri_angle);                        
        }   
                
                
        switch (mElasticityType){ //  0 ---linear compression & tension ; 1 --- Hertzian (non-linear compression, linear tension)
        
            case 0:
                equiv_young                     = 2.0 * myYoung * other_young / (myYoung + other_young);                
                kn                              = equiv_young * equiv_area * radius_sum_i; //KRATOS_M_PI * 0.5 * equiv_young * equiv_radius; //M: CANET FORMULA
                kt                              = kn / (2.0 + equiv_poisson + equiv_poisson);                
                
            break;

            case 1:
                equiv_young                     = myYoung * other_young / (other_young * (1.0 - myPoisson * myPoisson) + myYoung * (1.0 - other_poisson * other_poisson));
                kn                              = 1.3333333333333 * equiv_young * sqrt(0.5 * equiv_radius);
                kt                              = 2.0 * kn * (1.0 - equiv_poisson * equiv_poisson) / ((2.0 - equiv_poisson) * (1.0 + equiv_poisson));
          
            break;

            default:
                equiv_young                     = 2.0 * myYoung * other_young / (myYoung + other_young);
                kn                              = equiv_young * equiv_area * radius_sum_i; //KRATOS_M_PI * 0.5 * equiv_young * equiv_radius; //M: CANET FORMULA
                kt                              = kn / (2.0 + equiv_poisson + equiv_poisson);
            break;

        }//switch


        if (GetLnOfRestitCoeff() > 0.0 || other_ln_of_restit_coeff > 0.0){ // Limit expressions when the restitution coefficient tends to 0. Variable lnRestitCoeff is set to 1.0 (instead of minus infinite) by the problem type.
            equiv_visco_damp_coeff_normal = 2.0 * sqrt(equiv_mass * kn);

        }

        else {

            equiv_visco_damp_coeff_normal = - 2.0 * equiv_ln_of_restit_coeff * sqrt(equiv_mass * kn / (equiv_ln_of_restit_coeff * equiv_ln_of_restit_coeff + KRATOS_M_PI * KRATOS_M_PI));
        }
        
        double aux_norm_to_tang = sqrt(kt / kn);
        equiv_visco_damp_coeff_tangential = equiv_visco_damp_coeff_normal * aux_norm_to_tang;
		

      }
	  

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void SphericParticle::ComputeAdditionalForces(array_1d<double, 3>& externally_applied_force, array_1d<double, 3>& externally_applied_moment, ProcessInfo& rCurrentProcessInfo)
      {

        const array_1d<double,3>& gravity = rCurrentProcessInfo[GRAVITY];
        double mass = mSqrtOfRealMass * mSqrtOfRealMass;

        externally_applied_force[0] = mass * gravity[0];
        externally_applied_force[1] = mass * gravity[1];
        externally_applied_force[2] = mass * gravity[2];

      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void SphericParticle::CalculateViscoDamping(double LocalRelVel[3],
                                                  double ViscoDampingLocalContactForce[3],
                                                  double indentation,
                                                  double equiv_visco_damp_coeff_normal,
                                                  double equiv_visco_damp_coeff_tangential,
                                                  bool sliding)
      {

          //*** The compbrobation is component-wise since localContactForce and RelVel have in principle no relationship.
          // The visco force can be higher than the contact force only if they go to the same direction. (in my opinion)
          // But in oposite direction the visco damping can't overpass the force...

          if (mDampType > 0){

              if (mDampType == 11 || mDampType == 10){
                  ViscoDampingLocalContactForce[2] = - equiv_visco_damp_coeff_normal * LocalRelVel[2];
              }

              if (sliding == false && (mDampType == 1 || mDampType == 11)){ //only applied when no sliding to help to the regularized friction law or the spring convergence
                  ViscoDampingLocalContactForce[0] = - equiv_visco_damp_coeff_tangential * LocalRelVel[0];
                  ViscoDampingLocalContactForce[1] = - equiv_visco_damp_coeff_tangential * LocalRelVel[1];
              }
          }
      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void SphericParticle::AddUpForcesAndProject(double OldCoordSystem[3][3],
                                                  double LocalCoordSystem[3][3], 
                                                  double normal_force,
                                                  double LocalContactForce[3],
                                                  double LocalElasticContactForce[3],
                                                  double GlobalContactForce[3],
                                                  double GlobalElasticContactForce[3],
                                                  double ViscoDampingLocalContactForce[3],
                                                  //double ViscoDampingGlobalContactForce[3],
                                                  //array_1d<double, 3> &rContactForce,
                                                  array_1d<double, 3> &rElasticForce,
                                                  const double &i_neighbour_count)
      {
          for (unsigned int index = 0; index < 3; index++)
          {
              LocalContactForce[index] = LocalElasticContactForce[index] + ViscoDampingLocalContactForce[index];
          }

          GeometryFunctions::VectorLocal2Global(OldCoordSystem, LocalElasticContactForce, GlobalElasticContactForce);
          //GeometryFunctions::VectorLocal2Global(LocalCoordSystem, ViscoDampingLocalContactForce, ViscoDampingGlobalContactForce); //is this line necessary?????
          GeometryFunctions::VectorLocal2Global(OldCoordSystem, LocalContactForce, GlobalContactForce);
          
          double vect_normal_force[3] = {0.0};
          double global_normal_force[3] = {0.0};
          vect_normal_force[2] = normal_force;
          GeometryFunctions::VectorLocal2Global(LocalCoordSystem, vect_normal_force, global_normal_force);
          
          GlobalElasticContactForce[0] += global_normal_force[0];
          GlobalElasticContactForce[1] += global_normal_force[1];
          GlobalElasticContactForce[2] += global_normal_force[2];
          GlobalContactForce[0]        += global_normal_force[0];
          GlobalContactForce[1]        += global_normal_force[1];
          GlobalContactForce[2]        += global_normal_force[2];

          // Saving contact forces (We need to, since tangential elastic force is history-dependent)
          mOldNeighbourElasticContactForces[i_neighbour_count][0] = GlobalElasticContactForce[0];
          mOldNeighbourElasticContactForces[i_neighbour_count][1] = GlobalElasticContactForce[1];
          mOldNeighbourElasticContactForces[i_neighbour_count][2] = GlobalElasticContactForce[2];
          
          mOldNeighbourTotalContactForces[i_neighbour_count][0] = GlobalContactForce[0];
          mOldNeighbourTotalContactForces[i_neighbour_count][1] = GlobalContactForce[1];
          mOldNeighbourTotalContactForces[i_neighbour_count][2] = GlobalContactForce[2];
          

          //rContactForce[0] += GlobalContactForce[0];
          //rContactForce[1] += GlobalContactForce[1];
          //rContactForce[2] += GlobalContactForce[2];
          
          mContactForce[0] += GlobalContactForce[0];
          mContactForce[1] += GlobalContactForce[1];
          mContactForce[2] += GlobalContactForce[2];
  
          rElasticForce[0] += GlobalElasticContactForce[0];
          rElasticForce[1] += GlobalElasticContactForce[1];
          rElasticForce[2] += GlobalElasticContactForce[2];
      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************
 
      void SphericParticle::AddUpFEMForcesAndProject(double LocalCoordSystem[3][3],
                                                  double LocalContactForce[3],
                                                  double LocalElasticContactForce[3],
                                                  double GlobalContactForce[3],
                                                  double GlobalElasticContactForce[3],
                                                  double ViscoDampingLocalContactForce[3],
                                                  //double ViscoDampingGlobalContactForce[3],
                                                  //array_1d<double, 3> &rContactForce,
                                                  array_1d<double, 3> &rElasticForce,
                                                  const double &iRigidFaceNeighbour)
      {
          for (unsigned int index = 0; index < 3; index++)
          {
              LocalContactForce[index] = LocalElasticContactForce[index] + ViscoDampingLocalContactForce[index];
          }

          GeometryFunctions::VectorLocal2Global(LocalCoordSystem, LocalElasticContactForce, GlobalElasticContactForce);
          //GeometryFunctions::VectorLocal2Global(LocalCoordSystem, ViscoDampingLocalContactForce, ViscoDampingGlobalContactForce);
          GeometryFunctions::VectorLocal2Global(LocalCoordSystem, LocalContactForce, GlobalContactForce);

          // Saving contact forces (We need to, since tangential elastic force is history-dependent)
          mFemOldNeighbourContactForces[iRigidFaceNeighbour][0] = GlobalElasticContactForce[0];
          mFemOldNeighbourContactForces[iRigidFaceNeighbour][1] = GlobalElasticContactForce[1];
          mFemOldNeighbourContactForces[iRigidFaceNeighbour][2] = GlobalElasticContactForce[2];

          mContactForce[0] += GlobalContactForce[0];
          mContactForce[1] += GlobalContactForce[1];
          mContactForce[2] += GlobalContactForce[2];
  
          rElasticForce[0] += GlobalElasticContactForce[0];
          rElasticForce[1] += GlobalElasticContactForce[1];
          rElasticForce[2] += GlobalElasticContactForce[2];
          
          //----------------------------------------------------------------------------------------------------------------------------------//

          ///Global stored contact force between rigid face and particle, used by fem elements
          std::vector<double>& neighbour_rigid_faces_elastic_contact_force = this->mNeighbourRigidFacesElasticContactForce;
          std::vector<double>& neighbour_rigid_faces_total_contact_force = this->mNeighbourRigidFacesTotalContactForce;
          
          neighbour_rigid_faces_elastic_contact_force[3 * iRigidFaceNeighbour + 0] = GlobalElasticContactForce[0];
          neighbour_rigid_faces_elastic_contact_force[3 * iRigidFaceNeighbour + 1] = GlobalElasticContactForce[1];
          neighbour_rigid_faces_elastic_contact_force[3 * iRigidFaceNeighbour + 2] = GlobalElasticContactForce[2];

          neighbour_rigid_faces_total_contact_force[3 * iRigidFaceNeighbour + 0] = GlobalContactForce[0];
          neighbour_rigid_faces_total_contact_force[3 * iRigidFaceNeighbour + 1] = GlobalContactForce[1];
          neighbour_rigid_faces_total_contact_force[3 * iRigidFaceNeighbour + 2] = GlobalContactForce[2];

          //----------------------------------------------------------------------------------------------------------------------------------//
  
      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************
      
      void SphericParticle::MemberDeclarationFirstStep(const ProcessInfo& r_process_info)

      {

          // Paso al nodo la id del elemento cuando inicializo al mismo
          
          //if (!mInitializedVariablesFlag){
          if (this->IsNot(DEMFlags::HAS_INITIALIZED_VARIABLES) ){

               if (r_process_info[PRINT_EXPORT_ID] == 1){
                   this->GetGeometry()[0].FastGetSolutionStepValue(EXPORT_ID) = double(this->Id());
               }
                          

              mDampType                                    = r_process_info[DAMP_TYPE];
              mElasticityType                              = r_process_info[FORCE_CALCULATION_TYPE];
              if (r_process_info[ROTATION_OPTION])         this->Set(DEMFlags::HAS_ROTATION,true);
              else                                         this->Set(DEMFlags::HAS_ROTATION,false);
              if (r_process_info[ROLLING_FRICTION_OPTION]) this->Set(DEMFlags::HAS_ROLLING_FRICTION,true);
              else                                         this->Set(DEMFlags::HAS_ROLLING_FRICTION,false);
              if (r_process_info[CRITICAL_TIME_OPTION])    this->Set(DEMFlags::HAS_CRITICAL_TIME,true);
              else                                         this->Set(DEMFlags::HAS_CRITICAL_TIME,false);
                                                           this->Set(DEMFlags::HAS_ROTATION_SPRING,false);              

              AdditionalMemberDeclarationFirstStep(r_process_info);
          }

      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************
      
      void SphericParticle::NormalForceCalculation(double& normal_force, double kn, double indentation)
              
      {

         switch (mElasticityType){ //  0 ---linear compression & tension ; 1 --- Hertzian (non-linear compression, linear tension)
              case 0:
                  normal_force = kn * indentation;

              break;

              case 1:

                  normal_force = kn * pow(indentation, 1.5);
                  
                 
              break;

              default:

                  normal_force = kn * indentation;

              break;

          }//switch

          
      } //NormalForceCalculation
      
      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void SphericParticle::TangentialForceCalculation(const double normal_force, double LocalElasticContactForce[3], double LocalDeltDisp[3], const double& kt, const double& equiv_tg_of_fri_ang, bool& sliding)
      {

        LocalElasticContactForce[0] += - kt * LocalDeltDisp[0];  // 0: first tangential
        LocalElasticContactForce[1] += - kt * LocalDeltDisp[1];  // 1: second tangential

        double ShearForceNow = sqrt(LocalElasticContactForce[0] * LocalElasticContactForce[0] + LocalElasticContactForce[1] * LocalElasticContactForce[1]);
        double Frictional_ShearForceMax = equiv_tg_of_fri_ang * normal_force;

        if (Frictional_ShearForceMax < 0.0){
            Frictional_ShearForceMax = 0.0;
        }

        if ((ShearForceNow > Frictional_ShearForceMax) && (ShearForceNow != 0.0)){
            double ShearForceNowMaxRatio = Frictional_ShearForceMax / ShearForceNow;
            LocalElasticContactForce[0] = ShearForceNowMaxRatio * LocalElasticContactForce[0];
            LocalElasticContactForce[1] = ShearForceNowMaxRatio * LocalElasticContactForce[1];
            sliding = true;
        }

      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      double SphericParticle::GetInitialDelta(int index)
      
      {
	double delta = 0.0; //only available in continuum_particle

	return delta;
	
      }
      
      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      
      void SphericParticle::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo)
      {
          
        KRATOS_TRY

          if (rVariable == PARTICLE_ID){
              Output = mParticleId; // (NOT YET ACTIVE!!)
              return;
          }
          //CRITICAL DELTA CALCULATION

          if (rVariable == DELTA_TIME){
 
              double mass = mSqrtOfRealMass * mSqrtOfRealMass;
              double coeff = rCurrentProcessInfo[NODAL_MASS_COEFF];

              if (coeff > 1.0){
                  KRATOS_ERROR(std::runtime_error, "The coefficient assigned for vitual mass is larger than one, virtual_mass_coeff= ", coeff);
              }

              else if ((coeff == 1.0) && (rCurrentProcessInfo[VIRTUAL_MASS_OPTION])){
                  Output = 9.0E09;
              }

              else {

                  if (rCurrentProcessInfo[VIRTUAL_MASS_OPTION]){
                      mass = mass / (1 - coeff);
                  }

                  double K = KRATOS_M_PI * GetYoung() * mRadius; //M. Error, should be the same that the local definition.

                  Output = 0.34 * sqrt(mass / K);

                  //if (mRotationOption){
                  if (this->Is(DEMFlags::HAS_ROTATION) ){
                      Output *= 0.5; //factor for critical time step when rotation is allowed.
                  }

             }
              return;
          }

          if (rVariable == MAX_INDENTATION){
              CalculateMaxIndentation(Output, rCurrentProcessInfo[DISTANCE_TOLERANCE]);
              return;
          }

          if (rVariable == KINETIC_ENERGY){
              CalculateKineticEnergy(Output);
              return;
          }

          if (rVariable == ELASTIC_ENERGY_OF_CONTACTS){
              CalculateElasticEnergyOfContacts(Output);
              return;
          }                   
          
          if (rVariable == CALCULATE_COMPUTE_NEW_RIGID_FACE_NEIGHBOURS_HISTORICAL_DATA){
              ComputeNewRigidFaceNeighboursHistoricalData();
              return;
          }
          
          //AdditionalCalculate(rVariable, Output, rCurrentProcessInfo);

          KRATOS_CATCH("")

      }// Calculate

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void SphericParticle::Calculate(const Variable<array_1d<double, 3> >& rVariable, array_1d<double, 3>& Output, const ProcessInfo& rCurrentProcessInfo)
      {

          if (rVariable == MOMENTUM){
              CalculateMomentum(Output);
          }

          else if (rVariable == ANGULAR_MOMENTUM){
              CalculateLocalAngularMomentum(Output);
          }
      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void SphericParticle::Calculate(const Variable<Vector >& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo){}
      void SphericParticle::Calculate(const Variable<Matrix >& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo){}

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void SphericParticle::AdditionalMemberDeclarationFirstStep(const ProcessInfo& r_process_info){}
      void SphericParticle::AdditionalCalculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo){}

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************
      double SphericParticle::GetRadius()                                                      { return mRadius;                                                                        }
      void   SphericParticle::SetRadius(double radius)                                         { mRadius = radius;                                                                      }
      double SphericParticle::GetSqrtOfRealMass()                                              { return mSqrtOfRealMass;                                                                }
      void   SphericParticle::SetSqrtOfRealMass(double sqrt_of_real_mass)                      { mSqrtOfRealMass = sqrt_of_real_mass;                                                   }
      
      /*double SphericParticle::GetYoung()                                                       { return mYoung;                                                                         }
      double SphericParticle::GetRollingFriction()                                             { return mRollingFriction;                                                               }
      double SphericParticle::GetPoisson()                                                     { return mPoisson;                                                                       }
      double SphericParticle::GetTgOfFrictionAngle()                                           { return mTgOfFrictionAngle;                                                             }
      double SphericParticle::GetLnOfRestitCoeff()                                             { return mLnOfRestitCoeff;                                                               }*/

      /*void   SphericParticle::SetYoungFromProperties(double* young)                            { mYoung = *young;                                                                       }
      void   SphericParticle::SetPoissonFromProperties(double* poisson)                        { mPoisson = *poisson;                                                                   }
      void   SphericParticle::SetRollingFrictionFromProperties(double* rolling_friction)       { mRollingFriction = *rolling_friction;                                                  }
      void   SphericParticle::SetTgOfFrictionAngleFromProperties(double* tg_of_friction_angle) { mTgOfFrictionAngle = *tg_of_friction_angle;                                            }
      void   SphericParticle::SetLnOfRestitCoeffFromProperties(double* ln_of_restit_coeff)     { mLnOfRestitCoeff = *ln_of_restit_coeff;                                                }*/
      
      double SphericParticle::GetYoung()                                                       { return GetFastProperties()->GetYoung();                                                }
      double SphericParticle::GetRollingFriction()                                             { return GetFastProperties()->GetRollingFriction();                                      }
      double SphericParticle::GetPoisson()                                                     { return GetFastProperties()->GetPoisson();                                              }
      double SphericParticle::GetTgOfFrictionAngle()                                           { return GetFastProperties()->GetTgOfFrictionAngle() ;                                   }
      double SphericParticle::GetLnOfRestitCoeff()                                             { return GetFastProperties()->GetLnOfRestitCoeff();                                      }
      double SphericParticle::GetDensity()                                                     { return GetFastProperties()->GetDensity();                                              }
      
      void   SphericParticle::SetYoungFromProperties(double* young)                            { GetFastProperties()->SetYoungFromProperties( young);                                   }
      void   SphericParticle::SetRollingFrictionFromProperties(double* rolling_friction)       { GetFastProperties()->SetRollingFrictionFromProperties( rolling_friction);              }
      void   SphericParticle::SetPoissonFromProperties(double* poisson)                        { GetFastProperties()->SetPoissonFromProperties( poisson);                               }
      void   SphericParticle::SetTgOfFrictionAngleFromProperties(double* tg_of_friction_angle) { GetFastProperties()->SetTgOfFrictionAngleFromProperties( tg_of_friction_angle);        }
      void   SphericParticle::SetLnOfRestitCoeffFromProperties(double* ln_of_restit_coeff)     { GetFastProperties()->SetLnOfRestitCoeffFromProperties( ln_of_restit_coeff);            }  
      void   SphericParticle::SetDensityFromProperties(double* density)                        { GetFastProperties()->SetDensityFromProperties( density);                               }

      PropertiesProxy* SphericParticle::GetFastProperties()                                    { return mFastProperties;                                                                }
      void   SphericParticle::SetFastProperties(PropertiesProxy* pProps)                       { mFastProperties = pProps;                                                              }

}  // namespace Kratos.

