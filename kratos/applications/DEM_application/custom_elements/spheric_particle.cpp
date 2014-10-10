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

#define DEM_COPY_FIRST_TO_SECOND_3(a, b) a[0] = b[0]; a[1] = b[1]; a[2] = b[2];
#define DEM_SET_COMPONENTS_TO_ZERO_3(a)  a[0] = 0.0;  a[1] = 0.0;  a[2] = 0.0;

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

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void SphericParticle::Initialize()
      {
          KRATOS_TRY 

          mDimension                = 3;
          mRadius                   = GetGeometry()[0].FastGetSolutionStepValue(RADIUS);
          double density            = GetDensity();          
          double& sqrt_of_mass      = GetGeometry()[0].FastGetSolutionStepValue(SQRT_OF_MASS);          
          double mass               = 1.33333333333333333 * KRATOS_M_PI * density * mRadius * mRadius * mRadius;
          sqrt_of_mass              = sqrt(mass);          
          mSqrtOfRealMass           = sqrt_of_mass;

          if (this->Is(DEMFlags::HAS_ROTATION) ){
              double moment_of_inertia = 0.4 * mass * mRadius * mRadius;
              GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA) = moment_of_inertia;
          }

          else {
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

      void SphericParticle::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, double dt, const array_1d<double,3>& gravity)
      {
          KRATOS_TRY

          array_1d<double, 3> additional_forces;
          array_1d<double, 3> additionally_applied_moment;
          array_1d<double, 3> initial_rotation_moment;     
          array_1d<double, 3>& elastic_force = this->GetGeometry()[0].FastGetSolutionStepValue(ELASTIC_FORCES);

          mContactForce.clear();
          mContactMoment.clear();
          additional_forces.clear();
          additionally_applied_moment.clear();
          initial_rotation_moment.clear();
          elastic_force.clear();

          bool multi_stage_RHS = false;

          ComputeBallToBallContactForce(elastic_force, initial_rotation_moment, rCurrentProcessInfo, dt, multi_stage_RHS);

          if (mFemOldNeighbourIds.size() > 0){
              ComputeBallToRigidFaceContactForce(elastic_force, initial_rotation_moment, rCurrentProcessInfo);
          }
              
          ComputeAdditionalForces(additional_forces, additionally_applied_moment, rCurrentProcessInfo, gravity);

          rRightHandSideVector[0] = mContactForce[0]  + additional_forces[0];
          rRightHandSideVector[1] = mContactForce[1]  + additional_forces[1];
          rRightHandSideVector[2] = mContactForce[2]  + additional_forces[2];
          rRightHandSideVector[3] = mContactMoment[0] + additionally_applied_moment[0];
          rRightHandSideVector[4] = mContactMoment[1] + additionally_applied_moment[1];
          rRightHandSideVector[5] = mContactMoment[2] + additionally_applied_moment[2];
          
          array_1d<double,3>& total_forces  = this->GetGeometry()[0].FastGetSolutionStepValue(TOTAL_FORCES);
          array_1d<double,3>& total_moment  = this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT);

          for (int i = 0; i < 3; i++){
              total_forces[i] = rRightHandSideVector[i];
              total_moment[i] = rRightHandSideVector[3 + i];
          }	  

          KRATOS_CATCH( "" )
      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void SphericParticle::FirstCalculateRightHandSide(ProcessInfo& rCurrentProcessInfo, double dt)
      {
          KRATOS_TRY

          array_1d<double, 3> initial_rotation_moment;     
          array_1d<double, 3>& elastic_force = this->GetGeometry()[0].FastGetSolutionStepValue(ELASTIC_FORCES);

          mContactForce.clear();
          mContactMoment.clear();
          initial_rotation_moment.clear();
          elastic_force.clear();
          
          bool multi_stage_RHS = true;

          ComputeBallToBallContactForce(elastic_force, initial_rotation_moment, rCurrentProcessInfo, dt, multi_stage_RHS );

          std::vector<double>& neighbour_rigid_faces_elastic_contact_force = this->mNeighbourRigidFacesElasticContactForce;
          std::vector<double>& neighbour_rigid_faces_total_contact_force = this->mNeighbourRigidFacesTotalContactForce;
          std::fill(neighbour_rigid_faces_elastic_contact_force.begin(), neighbour_rigid_faces_elastic_contact_force.end(), 0.0);
          std::fill(neighbour_rigid_faces_total_contact_force.begin(), neighbour_rigid_faces_total_contact_force.end(), 0.0);
          
          if (mFemOldNeighbourIds.size() > 0){
              ComputeBallToRigidFaceContactForce(elastic_force, initial_rotation_moment, rCurrentProcessInfo);
          }

          KRATOS_CATCH( "" )
      }
      
      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void SphericParticle::CollectCalculateRightHandSide(ProcessInfo& rCurrentProcessInfo)
      {
          KRATOS_TRY  
          
          //LOOP OVER NEIGHBOURS BEGINS:          
          for (unsigned int i = 0; i < mNeighbourElements.size(); i++){
              SphericParticle* neighbour_iterator = mNeighbourElements[i];

              if (this->Is(NEW_ENTITY) && neighbour_iterator->Is(NEW_ENTITY)) continue;
              if (this->Id() < neighbour_iterator->Id()) continue;
              
              for (unsigned int j = 0; j < neighbour_iterator->mNeighbourElements.size(); j++){  //loop to find the neighbour of the neighbours which is me
                  SphericParticle* is_that_me = neighbour_iterator->mNeighbourElements[j];

                  if (is_that_me->Id() == this->Id()){
                      mOldNeighbourElasticContactForces[i] = - neighbour_iterator->mOldNeighbourElasticContactForces[j];
                      mOldNeighbourTotalContactForces[i]   = - neighbour_iterator->mOldNeighbourTotalContactForces[j];
                      mContactForce += mOldNeighbourTotalContactForces[i];
                      array_1d<double, 3>& rElasticForce = this->GetGeometry()[0].FastGetSolutionStepValue(ELASTIC_FORCES);
                      rElasticForce += mOldNeighbourElasticContactForces[i];
                      break;
                  }

              }              
              
          }          

          KRATOS_CATCH( "" )
      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void SphericParticle::FinalCalculateRightHandSide(ProcessInfo& rCurrentProcessInfo, double dt, const array_1d<double,3>& gravity)
      {
          KRATOS_TRY
                    
          if (this->Is(DEMFlags::HAS_ROTATION)){
                          
              //LOOP OVER NEIGHBOURS BEGINS:
              for( unsigned int i = 0; i < mNeighbourElements.size(); i++){
                  SphericParticle* neighbour_iterator = mNeighbourElements[i];

                  if (this->Is(NEW_ENTITY) && neighbour_iterator->Is(NEW_ENTITY)) continue;

                  array_1d<double, 3> other_to_me_vect = this->GetGeometry()[0].Coordinates() - neighbour_iterator->GetGeometry()[0].Coordinates();
                  double distance                      = sqrt(other_to_me_vect[0] * other_to_me_vect[0] + other_to_me_vect[1] * other_to_me_vect[1] + other_to_me_vect[2] * other_to_me_vect[2]);
                  double inv_distance                  = 1.0/distance;
                  double other_to_me_vect_unitary[3]   = {0.0};
                  other_to_me_vect_unitary[0]          = other_to_me_vect[0] * inv_distance;
                  other_to_me_vect_unitary[1]          = other_to_me_vect[1] * inv_distance;
                  other_to_me_vect_unitary[2]          = other_to_me_vect[2] * inv_distance;

                  double projection_to_local_axis2 = mOldNeighbourElasticContactForces[i][0] * other_to_me_vect_unitary[0]
                                                   + mOldNeighbourElasticContactForces[i][1] * other_to_me_vect_unitary[1]
                                                   + mOldNeighbourElasticContactForces[i][2] * other_to_me_vect_unitary[2];

                  const array_1d<double, 3> ang_vel = this->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
                  const double coeff_acc            = this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA) / dt;
                  array_1d<double, 3> initial_rotation_moment;
                  initial_rotation_moment           = coeff_acc * ang_vel; // the moment needed to stop the spin in one time step

                  if (this->Is(DEMFlags::HAS_ROTATION)){
                      ComputeMoments(projection_to_local_axis2, mOldNeighbourElasticContactForces[i], initial_rotation_moment, other_to_me_vect_unitary, neighbour_iterator);
                  }

              }

          } //if( this->Is(DEMFlags::HAS_ROTATION) ) 

          array_1d<double, 3> additional_forces;
          array_1d<double, 3> additionally_applied_moment;
          additional_forces.clear();
          additionally_applied_moment.clear();
          
          ComputeAdditionalForces(additional_forces, additionally_applied_moment, rCurrentProcessInfo, gravity);
          
          array_1d<double,3>& total_forces = this->GetGeometry()[0].FastGetSolutionStepValue(TOTAL_FORCES);
          array_1d<double,3>& total_moment = this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT);
                  
          total_forces = mContactForce  + additional_forces;
          total_moment = mContactMoment + additionally_applied_moment;

          KRATOS_CATCH( "" )
      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void SphericParticle::CalculateMaxBallToBallIndentation(double& rCurrentMaxIndentation)
      {
          rCurrentMaxIndentation = - std::numeric_limits<double>::max();

          for (unsigned int j = 0; j < mNeighbourElements.size(); j++){
              SphericParticle* i = mNeighbourElements[j];

              array_1d<double, 3> other_to_me_vect  = this->GetGeometry()[0].Coordinates() - i->GetGeometry()[0].Coordinates();
              double other_radius                   = i->GetRadius();
              double distance                       = sqrt(other_to_me_vect[0] * other_to_me_vect[0] +
                                                           other_to_me_vect[1] * other_to_me_vect[1] +
                                                           other_to_me_vect[2] * other_to_me_vect[2]);
              double radius_sum                     = mRadius + other_radius;
              double indentation                    = radius_sum - distance;

              rCurrentMaxIndentation = (indentation > rCurrentMaxIndentation) ? indentation : rCurrentMaxIndentation;
           }

      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void SphericParticle::CalculateMaxBallToFaceIndentation(double& rCurrentMaxIndentation)
      {

          rCurrentMaxIndentation = - std::numeric_limits<double>::max();
          array_1d<double, 3> node_coor_array = this->GetGeometry()[0].Coordinates();

          double node_coor[3];

          DEM_COPY_FIRST_TO_SECOND_3(node_coor, node_coor_array) //MSIMSI 1 can be optimized.

          for (unsigned int j = 0; j < mNeighbourRigidFaces.size(); j++){
              double DistPToB;
              DEMWall* pNeighbour = mNeighbourRigidFaces[j];
              double Coord[4][3] = {{0.0}, {0.0}, {0.0}, {0.0}};

              // Triangle
              DEM_COPY_FIRST_TO_SECOND_3(Coord[0], pNeighbour->GetGeometry()[0].Coordinates()) //MSIMSI 1 can be optimized with vector access.
              DEM_COPY_FIRST_TO_SECOND_3(Coord[1], pNeighbour->GetGeometry()[1].Coordinates())
              DEM_COPY_FIRST_TO_SECOND_3(Coord[2], pNeighbour->GetGeometry()[2].Coordinates())

              if (pNeighbour->GetGeometry().size() == 4){
                  DEM_COPY_FIRST_TO_SECOND_3(Coord[3], pNeighbour->GetGeometry()[3].Coordinates())
              }

              GeometryFunctions::QuickDistanceForAKnownNeighbour(Coord , node_coor, mRadius, DistPToB);
              double indentation = mRadius - DistPToB;
              rCurrentMaxIndentation = (indentation > rCurrentMaxIndentation) ? indentation : rCurrentMaxIndentation;

            }

      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void SphericParticle::CalculateKineticEnergy(double& rKineticEnergy)
      {
          const array_1d<double, 3>& vel    = this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
          const array_1d<double, 3> ang_vel = this->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
          const double moment_of_inertia    = this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA);
          double square_of_celerity         = vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2];
          double square_of_angular_celerity = ang_vel[0] * ang_vel[0] + ang_vel[1] * ang_vel[1] + ang_vel[2] * ang_vel[2];          

          rKineticEnergy = 0.5 * (mSqrtOfRealMass * mSqrtOfRealMass * square_of_celerity + moment_of_inertia * square_of_angular_celerity);
      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void SphericParticle::CalculateElasticEnergyOfContacts(double& rElasticEnergy) // Calculates the elastic energy stored in the sum of all the contacts shared by the particle and all its neighbours
      {
          double added_potential_energy_of_contacts   = 0.0;
          size_t i_neighbour_count                    = 0;
          double myYoung = GetYoung();
          double myPoisson = GetPoisson();
          
          std::vector<int> mTempNeighboursIds;
          std::vector<array_1d<double, 3> > mTempNeighbourElasticContactForces;
          std::vector<array_1d<double, 3> > mTempNeighbourTotalContactForces;
          ComputeNewNeighboursHistoricalData(mTempNeighboursIds, mTempNeighbourElasticContactForces, mTempNeighbourTotalContactForces);

          //for (ParticleWeakIteratorType neighbour_iterator = rNeighbours.begin(); neighbour_iterator != rNeighbours.end(); neighbour_iterator++){
          for (unsigned int i = 0; i < mNeighbourElements.size(); i++) {
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
          const double moment_of_inertia     = this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA);
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
         array_1d<double, 3> vector_of_zeros = ZeroVector(3);

         for (unsigned int i = 0; i < mNeighbourElements.size(); i++){
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

     //**************************************************************************************************************************************************
     //**************************************************************************************************************************************************

     void SphericParticle::ComputeNewRigidFaceNeighboursHistoricalData()
     {
       //ConditionWeakVectorType& rNeighbours  = this->GetValue(NEIGHBOUR_RIGID_FACES);
       std::vector<DEMWall*>& rNeighbours = this->mNeighbourRigidFaces;
       unsigned int new_size              = rNeighbours.size();
       std::vector<int> temp_neighbours_ids(new_size); //these two temporal vectors are very small, saving them as a member of the particle loses time.
       std::vector<array_1d<double, 3> > temp_neighbours_contact_forces(new_size);

       array_1d<double, 3> vector_of_zeros = ZeroVector(3);

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
          const array_1d<double,3> old_coord_target      = this->GetGeometry()[0].Coordinates() - delta_displ;
          const array_1d<double, 3 > & other_delta_displ = neighbour_iterator->GetGeometry()[0].FastGetSolutionStepValue(DELTA_DISPLACEMENT);
          const array_1d<double,3> old_coord_neigh       = neighbour_iterator->GetGeometry()[0].Coordinates() - other_delta_displ;
          
          array_1d<double,3> old_other_to_me_vect = old_coord_target - old_coord_neigh;

          const double old_distance = sqrt(old_other_to_me_vect[0]*old_other_to_me_vect[0] + old_other_to_me_vect[1]*old_other_to_me_vect[1] + old_other_to_me_vect[2]*old_other_to_me_vect[2]);
          
          GeometryFunctions::ComputeContactLocalCoordSystem(old_other_to_me_vect, old_distance, OldLocalCoordSystem); //Old Local Coord System

          // VELOCITIES AND DISPLACEMENTS
          const array_1d<double, 3 >& other_vel = neighbour_iterator->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
          
          RelVel[0] = (vel[0] - other_vel[0]);
          RelVel[1] = (vel[1] - other_vel[1]);
          RelVel[2] = (vel[2] - other_vel[2]);

          //DeltDisp in global cordinates
          DeltDisp[0] = (delta_displ[0] - other_delta_displ[0]);
          DeltDisp[1] = (delta_displ[1] - other_delta_displ[1]);
          DeltDisp[2] = (delta_displ[2] - other_delta_displ[2]);          
      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void SphericParticle::DisplacementDueToRotation(double DeltDisp[3],
                                                      double OldLocalCoordSystem[3][3],
                                                      const double& other_radius,
                                                      const double& dt,
                                                      const array_1d<double, 3>& ang_vel,
                                                      SphericParticle* neighbour_iterator)
      {
          array_1d<double, 3> other_ang_vel = neighbour_iterator->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
          double velA[3]                    = {0.0};
          double velB[3]                    = {0.0};
          double dRotaDisp[3]               = {0.0};

          double VelTemp[3]                 = {      ang_vel[0],       ang_vel[1],       ang_vel[2]};
          double OtherVelTemp[3]            = {other_ang_vel[0], other_ang_vel[1], other_ang_vel[2]};

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

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

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
              double rolling_friction_coeff       = GetRollingFriction() * mRadius;
              const double other_rolling_friction = neighbour_iterator->GetRollingFriction();
              double other_rolling_friction_coeff = other_rolling_friction * neighbour_iterator->GetRadius();
              double equiv_rolling_friction_coeff = std::min(rolling_friction_coeff,other_rolling_friction_coeff);

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
                       mContactMoment = - rInitialRotaMoment;
                   }

               } // if (equiv_rolling_friction_coeff != 0.0)

            } // if (mRotationDampType == 2)            

      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void SphericParticle::ComputeBallToBallContactForce(array_1d<double, 3>& rElasticForce,
                                                          array_1d<double, 3>& rInitialRotaMoment,
                                                          ProcessInfo& rCurrentProcessInfo,
                                                          double dt,
                                                          const bool multi_stage_RHS)
      {
          KRATOS_TRY

          // KINEMATICS

          const array_1d<double, 3>& vel         = this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
          const array_1d<double, 3>& delta_displ = this->GetGeometry()[0].FastGetSolutionStepValue(DELTA_DISPLACEMENT);
          const array_1d<double, 3>& ang_vel     = this->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
                    
          if (this->Is(DEMFlags::HAS_ROTATION) && !multi_stage_RHS){
              const double coeff_acc             = this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA) / dt;
              rInitialRotaMoment                 = coeff_acc * ang_vel; // the moment needed to stop the spin in one time step
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
          for (unsigned int i = 0; i < mNeighbourElements.size(); i++){
              SphericParticle* neighbour_iterator = mNeighbourElements[i];

              if (this->Is(NEW_ENTITY) && neighbour_iterator->Is(NEW_ENTITY)) continue;
              if (multi_stage_RHS  &&  this->Id() > neighbour_iterator->Id()) continue;
              
              // BASIC CALCULATIONS              
              array_1d<double, 3> other_to_me_vect    = this->GetGeometry()[0].Coordinates() - neighbour_iterator->GetGeometry()[0].Coordinates();
              const double &other_radius              = neighbour_iterator->GetRadius();
              double distance                         = sqrt(other_to_me_vect[0] * other_to_me_vect[0] + other_to_me_vect[1] * other_to_me_vect[1] + other_to_me_vect[2] * other_to_me_vect[2]);
              
              double radius_sum                       = mRadius + other_radius;
              double indentation                      = radius_sum - distance;
	          
              CalculateEquivalentConstitutiveParameters(other_to_me_vect, other_radius, radius_sum, kn, kt, equiv_visco_damp_coeff_normal, equiv_visco_damp_coeff_tangential, equiv_tg_of_fri_ang, neighbour_iterator);

              DEM_SET_COMPONENTS_TO_ZERO_3(DeltDisp)
              DEM_SET_COMPONENTS_TO_ZERO_3(LocalDeltDisp)
              DEM_SET_COMPONENTS_TO_ZERO_3(RelVel)
              DEM_SET_COMPONENTS_TO_ZERO_3(LocalRelVel)

              DEM_SET_COMPONENTS_TO_ZERO_3(LocalCoordSystem[0])
              DEM_SET_COMPONENTS_TO_ZERO_3(LocalCoordSystem[1])
              DEM_SET_COMPONENTS_TO_ZERO_3(LocalCoordSystem[2])

              DEM_SET_COMPONENTS_TO_ZERO_3(OldLocalCoordSystem[0])
              DEM_SET_COMPONENTS_TO_ZERO_3(OldLocalCoordSystem[1])
              DEM_SET_COMPONENTS_TO_ZERO_3(OldLocalCoordSystem[2])

              EvaluateDeltaDisplacement(DeltDisp, RelVel, LocalCoordSystem, OldLocalCoordSystem, other_to_me_vect, vel, delta_displ, neighbour_iterator, distance);

              if (this->Is(DEMFlags::HAS_ROTATION)){
                  DisplacementDueToRotation(DeltDisp, OldLocalCoordSystem, other_radius, dt, ang_vel, neighbour_iterator);
                }

              double normal_force                      = 0.0;
              double LocalContactForce[3]              = {0.0};
              double GlobalContactForce[3]             = {0.0};
              double LocalElasticContactForce[3]       = {0.0};
              double GlobalElasticContactForce[3]      = {0.0};
              double ViscoDampingLocalContactForce[3]  = {0.0};

              DEM_COPY_FIRST_TO_SECOND_3(GlobalElasticContactForce, mOldNeighbourElasticContactForces[i])

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

          if (iNeighborID == static_cast<int>(rObj_2->Id())){

              for(std::size_t inode = 0; inode < rObj_2->GetGeometry().size(); inode++){
                 other_to_me_vel += rObj_2->GetGeometry()(inode)->FastGetSolutionStepValue(VELOCITY) * Weight[inode];
                }

            }

          KRATOS_CATCH("")
      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void SphericParticle::ComputeBallToRigidFaceContactForce(array_1d<double, 3>& rElasticForce,
                                                               array_1d<double, 3>& rInitialRotaMoment,
                                                               ProcessInfo& rCurrentProcessInfo)
      {

        KRATOS_TRY

        //ConditionWeakVectorType& rNeighbours    = this->GetValue(NEIGHBOUR_RIGID_FACES);
        std::vector<DEMWall*>& rNeighbours = this->mNeighbourRigidFaces;

        double mTimeStep        = rCurrentProcessInfo[DELTA_TIME]; //TODO: to be removed and sent as an argument
        double myYoung          = GetYoung();
        double myPoisson        = GetPoisson();
        double area             = KRATOS_M_PI * mRadius * mRadius;

        array_1d<double, 3 > vel = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);

        for (unsigned int i=0; i<rNeighbours.size(); i++){
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
          DEM_COPY_FIRST_TO_SECOND_3(node_coor, node_coor_array) //MSIMSI 1 can be optimized.

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

          if (ContactType == 0){
              double Coord[4][3] = { {0.0},{0.0},{0.0},{0.0} };

              // Triangle

              DEM_COPY_FIRST_TO_SECOND_3(Coord[0], rNeighbours[i]->GetGeometry()[0].Coordinates()) //MSIMSI 1 can be optimized with vector access.
              DEM_COPY_FIRST_TO_SECOND_3(Coord[1], rNeighbours[i]->GetGeometry()[1].Coordinates())
              DEM_COPY_FIRST_TO_SECOND_3(Coord[2], rNeighbours[i]->GetGeometry()[2].Coordinates())

              if (rNeighbours[i]->GetGeometry().size() == 4){
                  DEM_COPY_FIRST_TO_SECOND_3(Coord[3], rNeighbours[i]->GetGeometry()[3].Coordinates())
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

          if (this->Is(DEMFlags::HAS_ROTATION)){
              double velA[3]      = {0.0};
              double dRotaDisp[3] = {0.0};

              array_1d<double, 3 > AngularVel= GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
              double Vel_Temp[3] = {AngularVel[0], AngularVel[1], AngularVel[2]};
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

          DEM_COPY_FIRST_TO_SECOND_3(GlobalElasticContactForce, mFemOldNeighbourContactForces[i])

          GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, GlobalElasticContactForce, LocalElasticContactForce);
          LocalElasticContactForce[0] +=  - ks_el * LocalDeltDisp[0];
          LocalElasticContactForce[1] +=  - ks_el * LocalDeltDisp[1];

          if (indentation > 0.0){
              LocalElasticContactForce[2] =  kn_el * indentation;
          }

          else {
              LocalElasticContactForce[2] = 0.0;
          }

          bool If_sliding = false;

          if (-LocalElasticContactForce[2] > 0.0){
              DEM_SET_COMPONENTS_TO_ZERO_3(LocalElasticContactForce)
          }

          else {
              double ShearForceMax = LocalElasticContactForce[2] * WallBallFriction;
              double ShearForceNow = sqrt(LocalElasticContactForce[0] * LocalElasticContactForce[0]
              +      LocalElasticContactForce[1] * LocalElasticContactForce[1]);

              if (ShearForceMax == 0.0){
                  LocalElasticContactForce[0] = 0.0;
                  LocalElasticContactForce[1] = 0.0;
              }

              else if (ShearForceNow > ShearForceMax){
                  double inv_ShearForceNow = 1.0 / ShearForceNow;
                  LocalElasticContactForce[0] = ShearForceMax * inv_ShearForceNow * LocalElasticContactForce[0];
                  LocalElasticContactForce[1] = ShearForceMax * inv_ShearForceNow * LocalElasticContactForce[1];

                  If_sliding = true;
              }
          }

          //---------------------------------------------DAMPING-FORCE-CALCULATION------------------------------------------------------------//

          double ViscoDampingLocalContactForce[3] = {0.0};

          if (indentation > 0.0){ // Compression , use visc damp
              double LocalRelVel[3] = {0.0};
              GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, DeltVel, LocalRelVel);
              double equiv_visco_damp_coeff_normal     = 0.0;
              double equiv_visco_damp_coeff_tangential = 0.0;
              double mass = mSqrtOfRealMass * mSqrtOfRealMass;

              if (rCurrentProcessInfo[DEMPACK_OPTION]){  //MSIMSI: DEMPACK
                  equiv_visco_damp_coeff_normal     = rCurrentProcessInfo[DEMPACK_DAMPING]*2.0*sqrt(kn_el/(mass+mass))*mass;//other_sqrt_of_mass*other_sqrt_of_mass))*equiv_mass;   // := 2d0* sqrt ( kn_el*(m1*m2)/(m1+m2) )  //MSIMSI: DEMPACK
                  equiv_visco_damp_coeff_tangential = equiv_visco_damp_coeff_normal * aux_norm_to_tang; // dempack no l'utilitza...
              }

              else {

                  if ( GetLnOfRestitCoeff() > 0.0){
                      equiv_visco_damp_coeff_normal = 2.0 * sqrt(mass * kn_el);
                  }

                  else {
                      equiv_visco_damp_coeff_normal = -2.0 * GetLnOfRestitCoeff() * sqrt(mass * kn_el / (GetLnOfRestitCoeff() * GetLnOfRestitCoeff() + KRATOS_M_PI * KRATOS_M_PI));
                  }

                  equiv_visco_damp_coeff_tangential = equiv_visco_damp_coeff_normal * sqrt(ks_el / kn_el);
              }

              CalculateViscoDamping(LocalRelVel, ViscoDampingLocalContactForce, indentation, equiv_visco_damp_coeff_normal,
              equiv_visco_damp_coeff_tangential, If_sliding);
            }

          double LocalContactForce[3]  = {0.0};
          double GlobalContactForce[3] = {0.0};

          AddUpFEMForcesAndProject(LocalCoordSystem, LocalContactForce,LocalElasticContactForce,GlobalContactForce,
                                  GlobalElasticContactForce,ViscoDampingLocalContactForce,rElasticForce, i);

          if (this->Is(DEMFlags::HAS_ROTATION)){
              ComputeMoments(LocalElasticContactForce[2], mFemOldNeighbourContactForces[i], rInitialRotaMoment, LocalCoordSystem[2], this); //WARNING: sending itself as the neighbour!!
            }

        }

        KRATOS_CATCH("")

      }// ComputeBallToRigidFaceContactForce

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

          this->Set(DEMFlags::HAS_INITIALIZED_VARIABLES);

          KRATOS_CATCH("")

      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void SphericParticle::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo){}

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
          
        if (mFastProperties == neighbour_iterator->GetFastProperties()){
            other_young   = myYoung;
            other_poisson = myPoisson;
            equiv_poisson = myPoisson;            
            other_ln_of_restit_coeff = myLnOfRestitCoeff;
            equiv_ln_of_restit_coeff = myLnOfRestitCoeff;            
            equiv_tg_of_fri_ang      = myTgOfFrictionAngle;
          }

        else {
            other_young = neighbour_iterator->GetYoung();                        
            other_poisson = neighbour_iterator->GetPoisson();
            equiv_poisson = 2.0 * myPoisson * other_poisson / (myPoisson + other_poisson);
            other_ln_of_restit_coeff = neighbour_iterator->GetLnOfRestitCoeff();
            equiv_ln_of_restit_coeff = 0.5 * (myLnOfRestitCoeff + other_ln_of_restit_coeff);
            const double other_tg_of_fri_angle = neighbour_iterator->GetTgOfFrictionAngle();                        
            equiv_tg_of_fri_ang      = 0.5 * (myTgOfFrictionAngle + other_tg_of_fri_angle);
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

      void SphericParticle::ComputeAdditionalForces(array_1d<double, 3>& externally_applied_force, array_1d<double, 3>& externally_applied_moment, ProcessInfo& rCurrentProcessInfo, const array_1d<double,3>& gravity)
      {
        double mass = mSqrtOfRealMass * mSqrtOfRealMass;
        externally_applied_force = mass * gravity;
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
                                                  array_1d<double, 3> &rElasticForce,
                                                  const double &i_neighbour_count)
      {

          for (unsigned int index = 0; index < 3; index++){
              LocalContactForce[index] = LocalElasticContactForce[index] + ViscoDampingLocalContactForce[index];
          }

          GeometryFunctions::VectorLocal2Global(OldCoordSystem, LocalElasticContactForce, GlobalElasticContactForce);
          GeometryFunctions::VectorLocal2Global(OldCoordSystem, LocalContactForce, GlobalContactForce);

          double global_normal_force[3] = {0.0};
          double vect_normal_force[3]   = {0.0};
          vect_normal_force[2] = normal_force;

          GeometryFunctions::VectorLocal2Global(LocalCoordSystem, vect_normal_force, global_normal_force);
          
          GlobalElasticContactForce[0] += global_normal_force[0];
          GlobalElasticContactForce[1] += global_normal_force[1];
          GlobalElasticContactForce[2] += global_normal_force[2];
          GlobalContactForce[0]        += global_normal_force[0];
          GlobalContactForce[1]        += global_normal_force[1];
          GlobalContactForce[2]        += global_normal_force[2];

          // Saving contact forces (We need to, since tangential elastic force is history-dependent)
          DEM_COPY_FIRST_TO_SECOND_3(mOldNeighbourElasticContactForces[i_neighbour_count], GlobalElasticContactForce)
          DEM_COPY_FIRST_TO_SECOND_3(mOldNeighbourTotalContactForces[i_neighbour_count], GlobalContactForce)

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
                                                  array_1d<double, 3> &rElasticForce,
                                                  const double &iRigidFaceNeighbour)
      {

          for (unsigned int index = 0; index < 3; index++){
              LocalContactForce[index] = LocalElasticContactForce[index] + ViscoDampingLocalContactForce[index];
          }

          GeometryFunctions::VectorLocal2Global(LocalCoordSystem, LocalElasticContactForce, GlobalElasticContactForce);
          GeometryFunctions::VectorLocal2Global(LocalCoordSystem, LocalContactForce, GlobalContactForce);

          // Saving contact forces (We need to, since tangential elastic force is history-dependent)
          std::copy(GlobalElasticContactForce, GlobalElasticContactForce + 3, mFemOldNeighbourContactForces[iRigidFaceNeighbour].begin());

          mContactForce[0] += GlobalContactForce[0];
          mContactForce[1] += GlobalContactForce[1];
          mContactForce[2] += GlobalContactForce[2];
  
          rElasticForce[0] += GlobalElasticContactForce[0];
          rElasticForce[1] += GlobalElasticContactForce[1];
          rElasticForce[2] += GlobalElasticContactForce[2];
          
          ///Global stored contact force between rigid face and particle, used by fem elements
          std::vector<double>& neighbour_rigid_faces_elastic_contact_force = this->mNeighbourRigidFacesElasticContactForce;
          std::vector<double>& neighbour_rigid_faces_total_contact_force = this->mNeighbourRigidFacesTotalContactForce;
          
          for (unsigned int i = 0; i < 3; ++i){
              neighbour_rigid_faces_elastic_contact_force[3 * iRigidFaceNeighbour + i] = GlobalElasticContactForce[i];
              neighbour_rigid_faces_total_contact_force[3 * iRigidFaceNeighbour + i]   = GlobalContactForce[i];
          }
  
      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************
      
      void SphericParticle::MemberDeclarationFirstStep(const ProcessInfo& r_process_info)

      {
          // Passing the element id to the node upon initialization
          
          if (this->IsNot(DEMFlags::HAS_INITIALIZED_VARIABLES)){

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

                  if (this->Is(DEMFlags::HAS_ROTATION) ){
                      Output *= 0.5; //factor for critical time step when rotation is allowed.
                  }

             }
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
      void   SphericParticle::SetFastProperties(PropertiesProxy* pProps)                       { mFastProperties = pProps;
                                                                                               }
      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

}  // namespace Kratos.

