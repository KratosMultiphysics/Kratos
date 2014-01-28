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
        : DiscreteElement(), mParticleId(-1), mInitializedVariablesFlag(0){}

      SphericParticle::SphericParticle(IndexType NewId, GeometryType::Pointer pGeometry)
        : DiscreteElement(NewId, pGeometry), mParticleId(NewId), mInitializedVariablesFlag(0){}

      SphericParticle::SphericParticle(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
      : DiscreteElement(NewId, pGeometry, pProperties), mParticleId(NewId), mInitializedVariablesFlag(0){}

      SphericParticle::SphericParticle(IndexType NewId, NodesArrayType const& ThisNodes)
      : DiscreteElement(NewId, ThisNodes), mParticleId(NewId), mInitializedVariablesFlag(0){}

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
          mRadius                   = GetGeometry()(0)->FastGetSolutionStepValue(RADIUS);
          mYoung                    = GetGeometry()(0)->FastGetSolutionStepValue(YOUNG_MODULUS);         
          mPoisson                  = GetGeometry()(0)->FastGetSolutionStepValue(POISSON_RATIO);
          mTgOfFrictionAngle        = GetGeometry()(0)->FastGetSolutionStepValue(PARTICLE_FRICTION);
          mLnOfRestitCoeff          = GetGeometry()(0)->FastGetSolutionStepValue(LN_OF_RESTITUTION_COEFF);
          double& density           = GetGeometry()(0)->FastGetSolutionStepValue(PARTICLE_DENSITY);
          //double& mass              = GetGeometry()(0)->FastGetSolutionStepValue(NODAL_MASS);
          double& sqrt_of_mass      = GetGeometry()(0)->FastGetSolutionStepValue(SQRT_OF_MASS);          

          double mass                      = 4.0 / 3.0 * M_PI * density * mRadius * mRadius * mRadius;
          sqrt_of_mass              = sqrt(mass);
          double moment_of_inertia         = 0.4 * mass * mRadius * mRadius; 
          
          mRealMass                 = mass;                    
          mSqrtOfRealMass           = sqrt_of_mass;
          mMomentOfInertia          = moment_of_inertia;

          if(mRotationOption){
            GetGeometry()(0)->FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA) = moment_of_inertia;
          }
          
          //OPTIMIZATION
          GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY_X_DOF_POS) = GetGeometry()[0].GetDofPosition(VELOCITY_X); 
          GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY_Y_DOF_POS) = GetGeometry()[0].GetDofPosition(VELOCITY_Y); 
          GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY_Z_DOF_POS) = GetGeometry()[0].GetDofPosition(VELOCITY_Z); 
          GetGeometry()(0)->FastGetSolutionStepValue(ANGULAR_VELOCITY_X_DOF_POS) = GetGeometry()[0].GetDofPosition(ANGULAR_VELOCITY_X); 
          GetGeometry()(0)->FastGetSolutionStepValue(ANGULAR_VELOCITY_Y_DOF_POS) = GetGeometry()[0].GetDofPosition(ANGULAR_VELOCITY_Y); 
          GetGeometry()(0)->FastGetSolutionStepValue(ANGULAR_VELOCITY_Z_DOF_POS) = GetGeometry()[0].GetDofPosition(ANGULAR_VELOCITY_Z);
          GetGeometry()(0)->FastGetSolutionStepValue(OLD_COORDINATES) = GetGeometry()(0)->Coordinates();

          CustomInitialize();

          KRATOS_CATCH( "" )
      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void SphericParticle::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
      {
          KRATOS_TRY

          array_1d<double, 3> contact_force;
          array_1d<double, 3> contact_moment;
          array_1d<double, 3> additionally_applied_force;
          array_1d<double, 3> additionally_applied_moment;
          array_1d<double, 3> initial_rotation_moment;     
           array_1d<double, 3>& elastic_force = this->GetGeometry()[0].FastGetSolutionStepValue(ELASTIC_FORCES);

          contact_force.clear();
          contact_moment.clear();
          additionally_applied_force.clear();
          additionally_applied_moment.clear();
          initial_rotation_moment.clear();
          elastic_force.clear();
          

          ComputeBallToBallContactForce(contact_force, contact_moment, elastic_force, initial_rotation_moment, rCurrentProcessInfo);

          if (mLimitSurfaceOption > 0){

              for (int surface_num = 0; surface_num < mLimitSurfaceOption; surface_num++){
                  ComputeBallToSurfaceContactForce(contact_force, contact_moment, initial_rotation_moment, surface_num, rCurrentProcessInfo);
              }
          }
          
          if (mLimitCylinderOption > 0){

              for (int cylinder_num = 0; cylinder_num < mLimitCylinderOption; cylinder_num++){
                  ComputeBallToCylinderContactForce(contact_force, contact_moment, initial_rotation_moment, cylinder_num, rCurrentProcessInfo);
              }
          }
		  
		  //Cfeng,RigidFace
		  if( mOldRigidFaceNeighbourIds.size() > 0)
		  {
			  ComputeBallToRigidFaceContactForce(contact_force, contact_moment, elastic_force, initial_rotation_moment, rCurrentProcessInfo);
		  }
          
          ComputeAdditionalForces(contact_force, contact_moment, additionally_applied_force, additionally_applied_moment, rCurrentProcessInfo);

          rRightHandSideVector[0] = contact_force[0]  + additionally_applied_force[0];
          rRightHandSideVector[1] = contact_force[1]  + additionally_applied_force[1];
          rRightHandSideVector[2] = contact_force[2]  + additionally_applied_force[2];
          rRightHandSideVector[3] = contact_moment[0] + additionally_applied_moment[0];
          rRightHandSideVector[4] = contact_moment[1] + additionally_applied_moment[1];
          rRightHandSideVector[5] = contact_moment[2] + additionally_applied_moment[2];
		  

          KRATOS_CATCH( "" )
      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void SphericParticle::CalculateMaxIndentation(double& rCurrentMaxIndentation, const double& rTolerance)
      {

          if (rCurrentMaxIndentation > rTolerance){
              ParticleWeakVectorType& rNeighbours = this->GetValue(NEIGHBOUR_ELEMENTS);
              double& radius                      = GetGeometry()(0)->FastGetSolutionStepValue(RADIUS); // cannot use the member variable mRadius because it may not be initialized
              rCurrentMaxIndentation              = 0.0;

              for (ParticleWeakIteratorType i = rNeighbours.begin(); i != rNeighbours.end(); i++){
                  array_1d<double, 3> other_to_me_vect  = this->GetGeometry()(0)->Coordinates() - i->GetGeometry()(0)->Coordinates();
                  double &other_radius                  = i->GetGeometry()(0)->FastGetSolutionStepValue(RADIUS);
                  double distance                       = sqrt(other_to_me_vect[0] * other_to_me_vect[0] +
                                                               other_to_me_vect[1] * other_to_me_vect[1] +
                                                               other_to_me_vect[2] * other_to_me_vect[2]);
                  double radius_sum                     = radius + other_radius;
                  double indentation                    = radius_sum - distance;

                  rCurrentMaxIndentation = (indentation > rCurrentMaxIndentation) ? indentation : rCurrentMaxIndentation;

              }

          }

      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void SphericParticle::CalculateKineticEnergy(double& rKineticEnergy)
      {
          const array_1d<double, 3>& vel    = this->GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY);
          const array_1d<double, 3> ang_vel = this->GetGeometry()(0)->FastGetSolutionStepValue(ANGULAR_VELOCITY);
          double square_of_celerity         = vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2];
          double square_of_angular_celerity = ang_vel[0] * ang_vel[0] + ang_vel[1] * ang_vel[1] + ang_vel[2] * ang_vel[2];

          rKineticEnergy = 0.5 * (mRealMass * square_of_celerity + mMomentOfInertia * square_of_angular_celerity);
      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void SphericParticle::CalculateElasticEnergyOfContacts(double& rElasticEnergy) // Calculates the elastic energy stored in the sum of all the contacts shared by the particle and all its neighbours
      {
          ParticleWeakVectorType& rNeighbours         = this->GetValue(NEIGHBOUR_ELEMENTS);
          double added_potential_energy_of_contacts   = 0.0;
          size_t i_neighbour_count                    = 0;

          ComputeNewNeighboursHistoricalData();

          for (ParticleWeakIteratorType neighbour_iterator = rNeighbours.begin();
              neighbour_iterator != rNeighbours.end(); neighbour_iterator++){
              const double &other_radius              = neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(RADIUS);
              double radius_sum                       = mRadius + other_radius;
              double radius_sum_i                     = 1 / radius_sum;
              double equiv_radius                     = 2 * mRadius * other_radius * radius_sum_i;
              double equiv_area                       = 0.25 * M_PI * equiv_radius * equiv_radius; // 0.25 is becouse we take only the half of the equivalent radius, corresponding to the case of one ball with radius Requivalent and other = radius 0.
              double corrected_area                   = equiv_area;
              double equiv_young;
              double equiv_poisson;
              double kn;
              double kt;


              // Getting neighbour properties
              const double &other_young           = neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(YOUNG_MODULUS);
              const double &other_poisson         = neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(POISSON_RATIO);

              equiv_young                         = 2 * mYoung * other_young / (mYoung + other_young);
              equiv_poisson                       = 2 * mPoisson * other_poisson / (mPoisson + other_poisson);
          

       
              kn                                  = equiv_young * corrected_area * radius_sum_i; //M_PI * 0.5 * equiv_young * equiv_radius; //M: CANET FORMULA
              kt                                  = kn / (2.0 + equiv_poisson + equiv_poisson);
           

              // Normal contribution

              double aux_power_of_contact_i_normal_force;

              switch (mElasticityType){ //  0 ---linear compression & tension ; 1 --- Hertzian (non-linear compression, linear tension)
                   case 0:

                       aux_power_of_contact_i_normal_force = mOldNeighbourContactForces[i_neighbour_count][2] * mOldNeighbourContactForces[i_neighbour_count][2];
                       added_potential_energy_of_contacts  += 0.5 * aux_power_of_contact_i_normal_force / kn;

                   break;

                   case 1:
                        aux_power_of_contact_i_normal_force = pow(fabs(mOldNeighbourContactForces[i_neighbour_count][2]), 5 / 3);
                        added_potential_energy_of_contacts  += 0.4 * aux_power_of_contact_i_normal_force / pow(kn, 2 / 3);

                   break;

                   default:

                       aux_power_of_contact_i_normal_force = mOldNeighbourContactForces[i_neighbour_count][2] * mOldNeighbourContactForces[i_neighbour_count][2];
                       added_potential_energy_of_contacts  += 0.5 * aux_power_of_contact_i_normal_force / kn;

                  break;

               }//switch

              // Tangential Contribution

              double aux_power_of_contact_i_tang_force = mOldNeighbourContactForces[i_neighbour_count][0] * mOldNeighbourContactForces[i_neighbour_count][0] + mOldNeighbourContactForces[i_neighbour_count][1] * mOldNeighbourContactForces[i_neighbour_count][1];
              added_potential_energy_of_contacts       += 0.5 * aux_power_of_contact_i_tang_force / kt;

              i_neighbour_count ++;

          }

          rElasticEnergy = added_potential_energy_of_contacts;
      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void SphericParticle::CalculateMomentum(array_1d<double, 3>& rMomentum)
      {
          const array_1d<double, 3>& vel = this->GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY);
          rMomentum = mRealMass * vel;
      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void SphericParticle::CalculateLocalAngularMomentum(array_1d<double, 3>& rAngularMomentum)
      {
          const array_1d<double, 3> ang_vel  = this->GetGeometry()(0)->FastGetSolutionStepValue(ANGULAR_VELOCITY);
          rAngularMomentum = mMomentOfInertia * ang_vel;
      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

     void SphericParticle::ComputeNewNeighboursHistoricalData()
     {
      
       ParticleWeakVectorType& rNeighbours  = this->GetValue(NEIGHBOUR_ELEMENTS);
       unsigned int new_size                = rNeighbours.size();
       unsigned int neighbour_counter       = 0;
       //vector<int> temp_neighbours_ids(new_size);
       std::vector<int>& temp_neighbours_ids = mTempNeighboursIds;
       temp_neighbours_ids.resize(new_size);
       //vector<array_1d<double, 3> > temp_neighbours_contact_forces(new_size);
       std::vector<array_1d<double, 3> >& temp_neighbours_contact_forces = mTempNeighboursContactForces;
       temp_neighbours_contact_forces.resize(new_size);
       array_1d<double, 3> vector_of_zeros;
       vector_of_zeros[0]                   = 0.0;
       vector_of_zeros[1]                   = 0.0;
       vector_of_zeros[2]                   = 0.0;

       for (ParticleWeakIteratorType i = rNeighbours.begin(); i != rNeighbours.end(); i++){

           temp_neighbours_ids[neighbour_counter] = static_cast<int>(i->Id());
           temp_neighbours_contact_forces[neighbour_counter] = vector_of_zeros;

           for (unsigned int j = 0; j != mOldNeighbourIds.size(); j++){

               if (static_cast<int>(i->Id()) == mOldNeighbourIds[j]){
                   temp_neighbours_contact_forces[neighbour_counter] = mOldNeighbourContactForces[j];
                   //temp_neighbours_contact_forces[neighbour_counter] = mOldNeighbourContactForces[j]; continuum
                   break;
               }

            }

            neighbour_counter++;
        }

        mOldNeighbourIds.swap(temp_neighbours_ids);
        mOldNeighbourContactForces.swap(temp_neighbours_contact_forces);

      }

	//////Cfeng,RigidFace
	void SphericParticle::ComputeNewRigidFaceNeighboursHistoricalData()
     {

       ConditionWeakVectorType& rNeighbours  = this->GetValue(NEIGHBOUR_RIGID_FACES);
       unsigned int new_size                = rNeighbours.size();
       unsigned int neighbour_counter       = 0;
       vector<int> temp_neighbours_ids(new_size); //these two temporal vectors are very small, saving them as a member of the particle loses time.
       vector<array_1d<double, 3> > temp_neighbours_contact_forces(new_size);       
       
       array_1d<double, 3> vector_of_zeros;
       vector_of_zeros[0]                   = 0.0;
       vector_of_zeros[1]                   = 0.0;
       vector_of_zeros[2]                   = 0.0;

       for (ConditionWeakIteratorType i = rNeighbours.begin(); i != rNeighbours.end(); i++){

           temp_neighbours_ids[neighbour_counter] = static_cast<int>(i->Id());
           temp_neighbours_contact_forces[neighbour_counter] = vector_of_zeros;

           for (unsigned int j = 0; j != mOldRigidFaceNeighbourIds.size(); j++)
			{

               if (static_cast<int>(i->Id()) == mOldRigidFaceNeighbourIds[j])
			   {
                   temp_neighbours_contact_forces[neighbour_counter] = mOldRigidFaceNeighbourContactForces[j];
                   break;
               }

            }

            neighbour_counter++;
        }

        mOldRigidFaceNeighbourIds.swap(temp_neighbours_ids);
        mOldRigidFaceNeighbourContactForces.swap(temp_neighbours_contact_forces);

      }


      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void SphericParticle::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo){}

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void SphericParticle::MassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
      {
          rMassMatrix(0,0) = mRealMass;
      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void SphericParticle::EvaluateDeltaDisplacement(double DeltDisp[3],
                                                      double RelVel[3],
                                                      //double NormalDir[3],
                                                      //double OldNormalDir[3],
                                                      double LocalCoordSystem[3][3],
                                                      double OldLocalCoordSystem[3][3],
                                                      array_1d<double, 3>& other_to_me_vect,
                                                      const array_1d<double, 3>& vel,
                                                      const array_1d<double, 3>& delta_displ,
                                                      ParticleWeakIteratorType neighbour_iterator,
                                                      double& distance)
      {


          // FORMING LOCAL CORDINATES

          //Notes: Since we will normally inherit the mesh from GiD, we respect the global system X,Y,Z [0],[1],[2]
          //In the local coordinates we will define the normal direction of the contact as the [2] component!!!!!
          //the way the normal direction is defined (other_to_me_vect) compression is positive

          //NormalDir[0] = other_to_me_vect[0];
          //NormalDir[1] = other_to_me_vect[1];
          //NormalDir[2] = other_to_me_vect[2];
          
          GeometryFunctions::ComputeContactLocalCoordSystem(other_to_me_vect, distance, LocalCoordSystem); //new Local Coord System (normalizes other_to_me_vect)

          // FORMING OLD LOCAL CORDINATES

          //array_1d<double,3> old_coord_target     = this->GetGeometry()(0)->GetInitialPosition() + this->GetGeometry()(0)->FastGetSolutionStepValue(DISPLACEMENT/*,1*/);          
          array_1d<double,3> old_coord_target     = this->GetGeometry()(0)->FastGetSolutionStepValue(OLD_COORDINATES);
          //array_1d<double,3> old_coord_neigh      = neighbour_iterator->GetGeometry()(0)->GetInitialPosition()+neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(DISPLACEMENT/*,1*/);
          array_1d<double,3> old_coord_neigh      = neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(OLD_COORDINATES);
          array_1d<double,3> Old_other_to_me_vect = old_coord_target - old_coord_neigh;                    

          //OldNormalDir[0] = Old_other_to_me_vect[0];   // M. this way the compresion is positive.
          //OldNormalDir[1] = Old_other_to_me_vect[1];
          //OldNormalDir[2] = Old_other_to_me_vect[2];
          
          //GeometryFunctions::ComputeContactLocalCoordSystem(OldNormalDir, OldLocalCoordSystem); //Old Local Coord System
          double old_distance = sqrt(Old_other_to_me_vect[0]*Old_other_to_me_vect[0] + Old_other_to_me_vect[1]*Old_other_to_me_vect[1] + Old_other_to_me_vect[2]*Old_other_to_me_vect[2]);
          
          GeometryFunctions::ComputeContactLocalCoordSystem(Old_other_to_me_vect, old_distance, OldLocalCoordSystem); //Old Local Coord System

          // VELOCITIES AND DISPLACEMENTS
          array_1d<double, 3 > other_vel          = neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY);
          array_1d<double, 3 > other_delta_displ  = neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(DELTA_DISPLACEMENT);

          RelVel[0] = (vel[0] - other_vel[0]);
          RelVel[1] = (vel[1] - other_vel[1]);
          RelVel[2] = (vel[2] - other_vel[2]);

          //DeltDisp in global cordinates
          DeltDisp[0] = (delta_displ[0] - other_delta_displ[0]);
          DeltDisp[1] = (delta_displ[1] - other_delta_displ[1]);
          DeltDisp[2] = (delta_displ[2] - other_delta_displ[2]);
      }

      void SphericParticle::DisplacementDueToRotation(double DeltDisp[3],
                                                      //double OldNormalDir[3],
                                                      double OldLocalCoordSystem[3][3],
                                                      const double& other_radius,
                                                      const double& dt,
                                                      const array_1d<double, 3>& ang_vel,
                                                      ParticleWeakIteratorType neighbour_iterator)
      {
          array_1d<double, 3> other_ang_vel     = neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(ANGULAR_VELOCITY);
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

      void SphericParticle::ComputeMoments(double LocalElasticContactForce[3],
                                           double GlobalElasticContactForce[3],
                                           //double InitialRotaMoment[3],
                                           array_1d<double, 3>& rInitialRotaMoment,
                                           double LocalCoordSystem[3][3],
                                           const double& other_radius,
                                           array_1d<double, 3>& rContactMoment,
                                           ParticleWeakIteratorType neighbour_iterator)
      {
          double MA[3]         = {0.0};
          //double RotaMoment[3] = {0.0};
          
          /*RotaMoment[0] = rContactMoment[0];
          RotaMoment[1] = rContactMoment[1];
          RotaMoment[2] = rContactMoment[2]; */         

          GeometryFunctions::CrossProduct(LocalCoordSystem[2], GlobalElasticContactForce, MA);

          /*RotaMoment[0] -= MA[0] * mRadius;
          RotaMoment[1] -= MA[1] * mRadius;
          RotaMoment[2] -= MA[2] * mRadius;*/
          rContactMoment[0] -= MA[0] * mRadius;
          rContactMoment[1] -= MA[1] * mRadius;
          rContactMoment[2] -= MA[2] * mRadius;

          // ROLLING FRICTION

          if (mRollingFrictionOption){  // rolling friction 
              double rolling_friction_coeff            = mRollingFriction * mRadius;
              double equiv_rolling_friction_coeff;


              const double& other_rolling_friction = neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(ROLLING_FRICTION);
              double other_rolling_friction_coeff  = other_rolling_friction * other_radius;
              equiv_rolling_friction_coeff         = 2 * rolling_friction_coeff * other_rolling_friction_coeff / (rolling_friction_coeff + other_rolling_friction_coeff);
              

              if (equiv_rolling_friction_coeff != 0.0){
                  double MaxRotaMoment[3]      = {0.0};
                  double CoordSystemMoment1[3] = {0.0};
                  double CoordSystemMoment2[3] = {0.0};
                  double MR[3]                 = {0.0};

                  /*MaxRotaMoment[0] = InitialRotaMoment[0] + RotaMoment[0];
                  MaxRotaMoment[1] = InitialRotaMoment[1] + RotaMoment[1];
                  MaxRotaMoment[2] = InitialRotaMoment[2] + RotaMoment[2];*/
                  MaxRotaMoment[0] = rInitialRotaMoment[0] + rContactMoment[0];
                  MaxRotaMoment[1] = rInitialRotaMoment[1] + rContactMoment[1];
                  MaxRotaMoment[2] = rInitialRotaMoment[2] + rContactMoment[2];
                             
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
                  double MR_now = det_MR * equiv_rolling_friction_coeff;
                  double MR_max = sqrt(MaxRotaMoment[0] * MaxRotaMoment[0] + MaxRotaMoment[1] * MaxRotaMoment[1] + MaxRotaMoment[2] * MaxRotaMoment[2]);

                  if (MR_max > MR_now){
                       /*RotaMoment[0] += MR[0] * equiv_rolling_friction_coeff;
                       RotaMoment[1] += MR[1] * equiv_rolling_friction_coeff;
                       RotaMoment[2] += MR[2] * equiv_rolling_friction_coeff;*/
                       rContactMoment[0] += MR[0] * equiv_rolling_friction_coeff;
                       rContactMoment[1] += MR[1] * equiv_rolling_friction_coeff;
                       rContactMoment[2] += MR[2] * equiv_rolling_friction_coeff;
                   }

                   else {
                       /*RotaMoment[0] = - InitialRotaMoment[0];
                       RotaMoment[1] = - InitialRotaMoment[1];
                       RotaMoment[2] = - InitialRotaMoment[2];*/
                       rContactMoment[0] = - rInitialRotaMoment[0];
                       rContactMoment[1] = - rInitialRotaMoment[1];
                       rContactMoment[2] = - rInitialRotaMoment[2];
                   }

               } // if (equiv_rolling_friction_coeff != 0.0)

            } // if (mRotationDampType == 2)
            
            /*rContactMoment[0] = RotaMoment[0];
            rContactMoment[1] = RotaMoment[1];
            rContactMoment[2] = RotaMoment[2];      */      

      }

      void SphericParticle::ComputeBallToBallContactForce(array_1d<double, 3>& rContactForce,
                                                          array_1d<double, 3>& rContactMoment,
                                                          array_1d<double, 3>& rElasticForce,
                                                          array_1d<double, 3>& rInitialRotaMoment,
                                                          ProcessInfo& rCurrentProcessInfo)
      {
          KRATOS_TRY

          ParticleWeakVectorType& rNeighbours    = this->GetValue(NEIGHBOUR_ELEMENTS);

          // KINEMATICS

          double dt = rCurrentProcessInfo[DELTA_TIME];
          double dt_i = 1 / dt;

          const array_1d<double, 3>& vel         = this->GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY);
          const array_1d<double, 3>& delta_displ = this->GetGeometry()(0)->FastGetSolutionStepValue(DELTA_DISPLACEMENT);
          const array_1d<double, 3>& ang_vel     = this->GetGeometry()(0)->FastGetSolutionStepValue(ANGULAR_VELOCITY);
          double RotaAcc[3]                      = {0.0};
          //double InitialRotaMoment[3]            = {0.0};

          if (mRotationOption){
              RotaAcc[0]                         = ang_vel[0] * dt_i;
              RotaAcc[1]                         = ang_vel[1] * dt_i;
              RotaAcc[2]                         = ang_vel[2] * dt_i;

              /*InitialRotaMoment[0]               = RotaAcc[0] * mMomentOfInertia;
              InitialRotaMoment[1]               = RotaAcc[1] * mMomentOfInertia;
              InitialRotaMoment[2]               = RotaAcc[2] * mMomentOfInertia;*/
              rInitialRotaMoment[0]               = RotaAcc[0] * mMomentOfInertia;
              rInitialRotaMoment[1]               = RotaAcc[1] * mMomentOfInertia;
              rInitialRotaMoment[2]               = RotaAcc[2] * mMomentOfInertia;
          }

          //LOOP OVER NEIGHBOURS BEGINS:
          size_t i_neighbour_count = 0;

          for (ParticleWeakIteratorType neighbour_iterator = rNeighbours.begin(); neighbour_iterator != rNeighbours.end(); neighbour_iterator++){

              if( this->Is(NEW_ENTITY) && neighbour_iterator->Is(NEW_ENTITY)) continue;
              
              // BASIC CALCULATIONS              
              array_1d<double, 3> other_to_me_vect    = this->GetGeometry()(0)->Coordinates() - neighbour_iterator->GetGeometry()(0)->Coordinates();
              const double &other_radius              = neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(RADIUS);
              double distance                         = sqrt(other_to_me_vect[0] * other_to_me_vect[0] + other_to_me_vect[1] * other_to_me_vect[1] + other_to_me_vect[2] * other_to_me_vect[2]);
              double radius_sum                       = mRadius + other_radius;
              double indentation                      = radius_sum - distance;
	      double kn;
	      double kt;
              double equiv_visco_damp_coeff_normal;
              double equiv_visco_damp_coeff_tangential;
              double equiv_tg_of_fri_ang;

              CalculateEquivalentConstitutiveParameters(other_to_me_vect, other_radius, radius_sum, kn, kt, equiv_visco_damp_coeff_normal, equiv_visco_damp_coeff_tangential, equiv_tg_of_fri_ang, neighbour_iterator);
              
              double DeltDisp[3]                       = {0.0};
              double LocalDeltDisp[3]                  = {0.0};
              double RelVel[3]                         = {0.0};
              double LocalRelVel[3]                    = {0.0};
              //double NormalDir[3]                      = {0.0};
              //double OldNormalDir[3]                   = {0.0};
              double LocalCoordSystem[3][3]            = {{0.0}, {0.0}, {0.0}};
              double OldLocalCoordSystem[3][3]         = {{0.0}, {0.0}, {0.0}};

              EvaluateDeltaDisplacement(DeltDisp, RelVel, /*NormalDir, OldNormalDir,*/ LocalCoordSystem, OldLocalCoordSystem, other_to_me_vect, vel, delta_displ, neighbour_iterator, distance);

              if (mRotationOption){
                  DisplacementDueToRotation(DeltDisp, /*OldNormalDir,*/ OldLocalCoordSystem, other_radius, dt, ang_vel, neighbour_iterator);
              }

              double LocalContactForce[3]              = {0.0};
              double GlobalContactForce[3]             = {0.0};
              double LocalElasticContactForce[3]       = {0.0};
              double GlobalElasticContactForce[3]      = {0.0};
              double ViscoDampingLocalContactForce[3]  = {0.0};
              double ViscoDampingGlobalContactForce[3] = {0.0};

              GlobalElasticContactForce[0]             = mOldNeighbourContactForces[i_neighbour_count][0];
              GlobalElasticContactForce[1]             = mOldNeighbourContactForces[i_neighbour_count][1];
              GlobalElasticContactForce[2]             = mOldNeighbourContactForces[i_neighbour_count][2];

              GeometryFunctions::VectorGlobal2Local(OldLocalCoordSystem, GlobalElasticContactForce, LocalElasticContactForce); // Here we recover the old local forces projected in the new coordinates in the way they were in the old ones; Now they will be increased if its the necessary
              GeometryFunctions::VectorGlobal2Local(OldLocalCoordSystem, DeltDisp, LocalDeltDisp);
              GeometryFunctions::VectorGlobal2Local(OldLocalCoordSystem, RelVel, LocalRelVel);

              // TRANSLATION FORCES

              bool sliding = false;

              if (indentation > 0.0){
                  NormalForceCalculation(LocalElasticContactForce, kn, indentation);

                  // TANGENTIAL FORCE. Incremental calculation. YADE develops a complicated "absolute method"

                  TangentialForceCalculation(LocalElasticContactForce, LocalDeltDisp, kt, equiv_tg_of_fri_ang, sliding);

                  // VISCODAMPING (applyied locally)

                  CalculateViscoDamping(LocalRelVel, ViscoDampingLocalContactForce, indentation, equiv_visco_damp_coeff_normal, equiv_visco_damp_coeff_tangential, sliding);
              }
              
              // Transforming to global forces and adding up

              AddUpForcesAndProject(LocalCoordSystem, LocalContactForce, LocalElasticContactForce, GlobalContactForce, GlobalElasticContactForce, ViscoDampingLocalContactForce, ViscoDampingGlobalContactForce, rContactForce, rElasticForce, i_neighbour_count);

              // ROTATION FORCES

              if (mRotationOption){
                  ComputeMoments(LocalElasticContactForce, GlobalElasticContactForce, /*InitialRotaMoment, */rInitialRotaMoment, LocalCoordSystem, other_radius, rContactMoment, neighbour_iterator);
              }

              i_neighbour_count++;

          }// for each neighbour

          /*rInitialRotaMoment[0] = InitialRotaMoment[0];
          rInitialRotaMoment[1] = InitialRotaMoment[1];
          rInitialRotaMoment[2] = InitialRotaMoment[2];*/

          KRATOS_CATCH("")
      }// ComputeBallToBallContactForce

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************


///////////////////////////////////////////////////Cfeng,RigidFace Contact Force calculation//////////////////////////////
//////////*******************************************************07,Oct,2013*******************************////////////


void SphericParticle::ComputeRigidFaceToMeVelocity(ConditionWeakIteratorType rObj_2, std::size_t ino, 
                             double LocalCoordSystem[3][3], double & DistPToB, array_1d<double, 3 > &other_to_me_vel)
{
	 KRATOS_TRY
	 
	
	 double Weight[4] = {0.0};
	 
	 Vector & RF_Pram= this->GetValue(NEIGHBOUR_RIGID_FACES_PRAM);
	 
	
	 int ino1 = ino * 15;
	
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
	 
	 if(iNeighborID == static_cast<int>(rObj_2->Id()) )
	 {
		for(std::size_t inode = 0; inode < rObj_2->GetGeometry().size(); inode++)
		{
		   other_to_me_vel += rObj_2->GetGeometry()(inode)->FastGetSolutionStepValue(VELOCITY) * Weight[inode];
		}		 
	 }
	
	 KRATOS_CATCH("")
}



  void SphericParticle::ComputeBallToRigidFaceContactForce(array_1d<double, 3>& rContactForce,
                                                          array_1d<double, 3>& rContactMoment,
                                                          array_1d<double, 3>& rElasticForce,
                                                          array_1d<double, 3>& rInitialRotaMoment,
                                                          ProcessInfo& rCurrentProcessInfo)
	{
		  
		  KRATOS_TRY
		  

          ConditionWeakVectorType& rNeighbours    = this->GetValue(NEIGHBOUR_RIGID_FACES);
          
		
		double mTimeStep    = rCurrentProcessInfo[DELTA_TIME];
		/////int CalRotateOption = rCurrentProcessInfo[RIGID_FACE_FLAG];

		double Friction       = mTgOfFrictionAngle;
		double young          = mYoung;
		double poisson        = mPoisson;
		double radius         = GetGeometry()(0)->FastGetSolutionStepValue(RADIUS);
		double area           = M_PI * radius * radius;
		double kn             = 1.0*young * area / (2.0 * radius);      
		double ks             = kn / (2.0 * (1.0 + poisson));
                //const double &mLnOfRestitCoeff        = this->GetGeometry()(0)->FastGetSolutionStepValue(LN_OF_RESTITUTION_COEFF);
		array_1d<double, 3 > vel = GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY);
                
        std::size_t iRigidFaceNeighbour = 0;

        for(ConditionWeakIteratorType ineighbour = rNeighbours.begin(); ineighbour != rNeighbours.end(); ineighbour++)
        {
			
            double LocalContactForce[3]  = {0.0};
            double GlobalContactForce[3] = {0.0};
            double GlobalContactForceOld[3] = {0.0};			            

            array_1d<double, 3 > other_to_me_vel;
            noalias(other_to_me_vel) = ZeroVector(3);

            double LocalCoordSystem[3][3] = {{0.0}, {0.0}, {0.0}};			
            double DistPToB = 0.0;

	    ComputeRigidFaceToMeVelocity(ineighbour, iRigidFaceNeighbour, LocalCoordSystem, DistPToB, other_to_me_vel);						

            double DeltDisp[3] = {0.0};
            double DeltVel [3] = {0.0};

            DeltVel[0] = (vel[0] - other_to_me_vel[0]);
            DeltVel[1] = (vel[1] - other_to_me_vel[1]);
            DeltVel[2] = (vel[2] - other_to_me_vel[2]);

            // For translation movement delt displacement
            DeltDisp[0] = DeltVel[0] * mTimeStep;
            DeltDisp[1] = DeltVel[1] * mTimeStep;
            DeltDisp[2] = DeltVel[2] * mTimeStep;


            if (mRotationOption)
            {
                double velA[3]   = {0.0};
                double dRotaDisp[3] = {0.0};

                array_1d<double, 3 > AngularVel= GetGeometry()(0)->FastGetSolutionStepValue(ANGULAR_VELOCITY);
                double Vel_Temp[3] = { AngularVel[0], AngularVel[1], AngularVel[2]};
                GeometryFunctions::CrossProduct(Vel_Temp, LocalCoordSystem[2], velA);

                dRotaDisp[0] = -velA[0] * radius;
                dRotaDisp[1] = -velA[1] * radius;
                dRotaDisp[2] = -velA[2] * radius;

                //////contribution of the rotation vel
                DeltDisp[0] += dRotaDisp[0] * mTimeStep;
                DeltDisp[1] += dRotaDisp[1] * mTimeStep;
                DeltDisp[2] += dRotaDisp[2] * mTimeStep;
            }


            double LocalDeltDisp[3] = {0.0};
            GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, DeltDisp, LocalDeltDisp);
			
			

            //////120323,for global storage
			
		    GlobalContactForceOld[0] = mOldRigidFaceNeighbourContactForces[iRigidFaceNeighbour][0];
		    GlobalContactForceOld[1] = mOldRigidFaceNeighbourContactForces[iRigidFaceNeighbour][1];
		    GlobalContactForceOld[2] = mOldRigidFaceNeighbourContactForces[iRigidFaceNeighbour][2];
			
			
            GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, GlobalContactForceOld, LocalContactForce);
            LocalContactForce[0] +=  - ks * LocalDeltDisp[0];
            LocalContactForce[1] +=  - ks * LocalDeltDisp[1];
            //LocalContactForce[2] +=  - kn * LocalDeltDisp[2];
			
		
			 if(DistPToB < radius)
			{
				LocalContactForce[2] =  -kn * (DistPToB - radius);
			}
			else
			{
				LocalContactForce[2] = 0.0;
			}
// 			if(this->Id()==2)
//             {
//             KRATOS_WATCH(*mpTimeStep)
//             KRATOS_WATCH(rCurrentProcessInfo[TIME_STEPS])
//             }
// 			if(*mpTimeStep == 0)
//             {
// 			KRATOS_WATCH((DistPToB - radius))
//             }
				

            bool If_sliding = false;
			
            if (-LocalContactForce[2] > 0.0)
            {
                LocalContactForce[0] = 0.0;
                LocalContactForce[1]  = 0.0;
                LocalContactForce[2]  = 0.0;
            }
            else
            {

                double ShearForceMax = LocalContactForce[2] * Friction;
                double ShearForceNow = sqrt(LocalContactForce[0] * LocalContactForce[0]
                                     +      LocalContactForce[1] * LocalContactForce[1]);


                //Cfeng: for shear failure
                if(ShearForceMax == 0.0)
                {
                    LocalContactForce[0] = 0.0;
                    LocalContactForce[1] = 0.0;
                }
                else if(ShearForceNow > ShearForceMax)
                {
                    LocalContactForce[0] = ShearForceMax / ShearForceNow * LocalContactForce[0];
                    LocalContactForce[1] = ShearForceMax / ShearForceNow * LocalContactForce[1];
					
					If_sliding = true;
                }
            }
			

            GeometryFunctions::VectorLocal2Global(LocalCoordSystem, LocalContactForce, GlobalContactForce);

            mOldRigidFaceNeighbourContactForces[iRigidFaceNeighbour][0] = GlobalContactForce[0];
            mOldRigidFaceNeighbourContactForces[iRigidFaceNeighbour][1] = GlobalContactForce[1];
            mOldRigidFaceNeighbourContactForces[iRigidFaceNeighbour][2] = GlobalContactForce[2];
			
	    ///Global stored contact force between rigid face and particle, used by fem elements
            Vector& neighbour_rigid_faces_contact_force = this->GetValue(NEIGHBOUR_RIGID_FACES_CONTACT_FORCE);
	    neighbour_rigid_faces_contact_force[3 * iRigidFaceNeighbour + 0] = GlobalContactForce[0];
	    neighbour_rigid_faces_contact_force[3 * iRigidFaceNeighbour + 1] = GlobalContactForce[1];
	    neighbour_rigid_faces_contact_force[3 * iRigidFaceNeighbour + 2] = GlobalContactForce[2];
			

            rContactForce[0] += GlobalContactForce[0];
            rContactForce[1] += GlobalContactForce[1];
            rContactForce[2] += GlobalContactForce[2];
			
			/////////////////////////Damping Force Calculation//////////////////////////7
			if(LocalContactForce[2] > 0.0) // Compression , use visc damp
			{			
				double LocalRelVel[3]                   = {0.0};
				double ViscoDampingLocalContactForce[3] = {0.0};
				GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, DeltVel, LocalRelVel);
				double equiv_visco_damp_coeff_normal = 0.0;
				double equiv_visco_damp_coeff_tangential = 0.0;
				
				if (mLnOfRestitCoeff > 0.0)
				{
					equiv_visco_damp_coeff_normal = 2.0 * sqrt(mRealMass * kn);
				}
				else 
				{
					equiv_visco_damp_coeff_normal = -2.0 * mLnOfRestitCoeff * sqrt(mRealMass * kn / (mLnOfRestitCoeff * mLnOfRestitCoeff + M_PI * M_PI));
				}

				equiv_visco_damp_coeff_tangential = equiv_visco_damp_coeff_normal * sqrt(ks / kn);
			
				CalculateViscoDamping(LocalRelVel, ViscoDampingLocalContactForce,DistPToB, equiv_visco_damp_coeff_normal,
									  equiv_visco_damp_coeff_tangential, If_sliding);
				
				double GlobalDampForce[3]   = {0.0};
				GeometryFunctions::VectorLocal2Global(LocalCoordSystem, ViscoDampingLocalContactForce, GlobalDampForce);	
				
				rContactForce[0] += GlobalDampForce[0];
				rContactForce[1] += GlobalDampForce[1];
				rContactForce[2] += GlobalDampForce[2];
			}
			
			/////////////////////////////////////////////////////////////////////////////

            if ( mRotationOption)
            {
                double MA[3] = {0.0};
                GeometryFunctions::CrossProduct(LocalCoordSystem[2], GlobalContactForce, MA);
                rContactMoment[0] -= MA[0] * radius;
                rContactMoment[1] -= MA[1] * radius;
                rContactMoment[2] -= MA[2] * radius;
            }


            iRigidFaceNeighbour++;

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

          const array_1d<double, 3> vel         = this->GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY);
          const array_1d<double, 3> delta_displ = this->GetGeometry()(0)->FastGetSolutionStepValue(DELTA_DISPLACEMENT);
          const array_1d<double, 3> ang_vel     = this->GetGeometry()(0)->FastGetSolutionStepValue(ANGULAR_VELOCITY);
          double InitialRotaMoment[3]           = {0.0};
          double visco_damp_coeff_normal;
          double visco_damp_coeff_tangential;
          //const double &mLnOfRestitCoeff        = this->GetGeometry()(0)->FastGetSolutionStepValue(LN_OF_RESTITUTION_COEFF);

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
          
          array_1d<double, 3> point_coor = this->GetGeometry()(0)->Coordinates();

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

          if (indentation > 0.0 && point_coor[1] <= 0.25){
                     // MACRO PARAMETERS

                     double kn;
                     double kt;
                     double equiv_young;

                     switch (mElasticityType){ //  0 ---linear compression & tension ; 1 --- Hertzian (non-linear compression, linear tension)
                
                         case 0:
                             kn = M_PI * 0.5 * mYoung * mRadius; //M_PI * 0.5 * equiv_young * equiv_radius; //M: CANET FORMULA
                             kt = kn / (2.0 * (1.0 + mPoisson));
                        
                         break;

                         case 1:
                             equiv_young = mYoung / (1- mPoisson * mPoisson);
                             kn                 = (4/3) * equiv_young * sqrt(mRadius);
                             kt                 = 2.0 * kn * (1 - mPoisson * mPoisson) / ((2.0 - mPoisson) * (1 + mPoisson));               
                 
                         break;

                         default:
                             kn = M_PI * 0.5 * mYoung * mRadius; //M_PI * 0.5 * equiv_young * equiv_radius; //M: CANET FORMULA
                             kt = kn / (2.0 * (1.0 + mPoisson));

                         break;

                     }//switch


                 // Historical minimun K for the critical time:
                 if (mCriticalTimeOption){
                     double historic = rCurrentProcessInfo[HISTORICAL_MIN_K];

                     if ((kn < historic) || (kt < historic)){
                         historic = std::min(kn, kt);
                     }

                 }   
                 double aux_norm_to_tang = sqrt(kt / kn);

                 if (mLnOfRestitCoeff > 0.0){
                     visco_damp_coeff_normal     = 2 * sqrt(mRealMass * kn);
                     visco_damp_coeff_tangential = visco_damp_coeff_normal * aux_norm_to_tang; // 2 * sqrt(mass * kt);
                 }

          else {
                     visco_damp_coeff_normal     = - 2 * mLnOfRestitCoeff * sqrt(mRealMass * kn / (mLnOfRestitCoeff * mLnOfRestitCoeff + M_PI * M_PI));
                     visco_damp_coeff_tangential = visco_damp_coeff_normal * aux_norm_to_tang; //= -(2 * log(restitution_coeff) * sqrt(mass * kt)) / (sqrt((log(restitution_coeff) * log(restitution_coeff)) + (M_PI * M_PI)));
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

                 if (mRotationOption){
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
                 
                 double dyn_friction_angle =  surface_friction * M_PI / 180;
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
                 double ViscoDampingGlobalContactForce[3] = {0.0};
                 double GlobalContactForce[3]             = {0.0};

                 for (unsigned int index = 0; index < 3; index++){
                     LocalContactForce[index] = LocalElasticContactForce[index]  + ViscoDampingLocalContactForce[index];
                 }

                 GeometryFunctions::VectorLocal2Global(LocalCoordSystem, LocalElasticContactForce, GlobalElasticContactForce);
                 GeometryFunctions::VectorLocal2Global(LocalCoordSystem, ViscoDampingLocalContactForce, ViscoDampingGlobalContactForce);
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

                 if (mRotationOption){
                     double MA[3]                     = {0.0};
                     double RotaMoment[3]             = {0.0};

                     RotaMoment[0] = rContactMoment[0];
                     RotaMoment[1] = rContactMoment[1];
                     RotaMoment[2] = rContactMoment[2];

                     GeometryFunctions::CrossProduct(LocalCoordSystem[2], GlobalElasticContactForce, MA);

                     RotaMoment[0] -= MA[0] * mRadius;
                     RotaMoment[1] -= MA[1] * mRadius;
                     RotaMoment[2] -= MA[2] * mRadius;

                     if (mRollingFrictionOption){  // Rolling friction type
                         double rolling_friction             = this->GetGeometry()(0)->FastGetSolutionStepValue(ROLLING_FRICTION);
                         double rolling_friction_coeff       = rolling_friction * mRadius;

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

          array_1d<double,3> CylinderAxisDir           = rCurrentProcessInfo[CYLINDER_AXIS_DIR_1];
          array_1d<double,3> InitialBaseCylinderCentre = rCurrentProcessInfo[INITIAL_BASE_CYLINDER_CENTRE_1];
          double CylinderRadius                        = rCurrentProcessInfo[CYLINDER_RADIUS_1];
		  double CylinderVelocity                      = rCurrentProcessInfo[CYLINDER_VELOCITY_1];
          double CylinderAngularVelocity               = rCurrentProcessInfo[CYLINDER_ANGULAR_VELOCITY_1];          
          double CylinderFriction                      = rCurrentProcessInfo[CYLINDER_FRICTION_1];  
          //const double &mLnOfRestitCoeff        = this->GetGeometry()(0)->FastGetSolutionStepValue(LN_OF_RESTITUTION_COEFF);

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

             const array_1d<double, 3> vel         = this->GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY);
             const array_1d<double, 3> delta_displ = this->GetGeometry()(0)->FastGetSolutionStepValue(DELTA_DISPLACEMENT);
             const array_1d<double, 3> ang_vel     = this->GetGeometry()(0)->FastGetSolutionStepValue(ANGULAR_VELOCITY);
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
             
             array_1d<double, 3> point_coord                   = this->GetGeometry()(0)->Coordinates();

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
                             kn = M_PI * 0.5 * mYoung * mRadius; //M_PI * 0.5 * equiv_young * equiv_radius; //M: CANET FORMULA
                             kt = kn / (2.0 * (1.0 + mPoisson));
                        
                         break;

                         case 1:
                             equiv_young      = mYoung / (1- mPoisson * mPoisson);
                             effective_radius = mRadius * CylinderRadius / (CylinderRadius - mRadius);
                             kn                      = (4/3) * equiv_young * sqrt(effective_radius);
                             kt                      = 2.0 * kn * (1 - mPoisson * mPoisson) / ((2.0 - mPoisson) * (1 + mPoisson));               
                 
                         break;

                         default:
                             kn = M_PI * 0.5 * mYoung * mRadius; //M_PI * 0.5 * equiv_young * equiv_radius; //M: CANET FORMULA
                             kt = kn / (2.0 * (1.0 + mPoisson));

                         break;

                     }//switch
                    
        
                     // Historical minimun K for the critical time:
                     if (mCriticalTimeOption){
                         double historic = rCurrentProcessInfo[HISTORICAL_MIN_K];
    
                         if ((kn < historic) || (kt < historic)){
                             historic = std::min(kn, kt);
                         }

                     }
                     double aux_norm_to_tang = sqrt(kt / kn);

                     if (mLnOfRestitCoeff > 0.0){
                         visco_damp_coeff_normal     = 2 * sqrt(mRealMass * kn);
                         visco_damp_coeff_tangential = visco_damp_coeff_normal * aux_norm_to_tang; // 2 * sqrt(mass * kt);
                     }
    
                     else {
                         visco_damp_coeff_normal     = - 2 * mLnOfRestitCoeff * sqrt(mRealMass * kn / (mLnOfRestitCoeff * mLnOfRestitCoeff + M_PI * M_PI));
                         visco_damp_coeff_tangential = visco_damp_coeff_normal * aux_norm_to_tang; //= -(2 * log(restitution_coeff) * sqrt(mass * kt)) / (sqrt((log(restitution_coeff) * log(restitution_coeff)) + (M_PI * M_PI)));
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

                     if (mRotationOption){
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

                     double dyn_friction_angle =  CylinderFriction * M_PI / 180;
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
                     double ViscoDampingGlobalContactForce[3] = {0.0};
                     double GlobalContactForce[3]             = {0.0};

                     for (unsigned int index = 0; index < 3; index++){
                         LocalContactForce[index] = LocalElasticContactForce[index]  + ViscoDampingLocalContactForce[index];
                     }

                     GeometryFunctions::VectorLocal2Global(LocalCoordSystem, LocalElasticContactForce, GlobalElasticContactForce);
                     GeometryFunctions::VectorLocal2Global(LocalCoordSystem, ViscoDampingLocalContactForce, ViscoDampingGlobalContactForce);
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

                     if (mRotationOption){
                         double MA[3]                     = {0.0};
                         double RotaMoment[3]             = {0.0};
                         
                         RotaMoment[0] = rContactMoment[0];
                         RotaMoment[1] = rContactMoment[1];
                         RotaMoment[2] = rContactMoment[2];                         

                         GeometryFunctions::CrossProduct(LocalCoordSystem[2], GlobalElasticContactForce, MA);

                         RotaMoment[0] -= MA[0] * mRadius;
                         RotaMoment[1] -= MA[1] * mRadius;
                         RotaMoment[2] -= MA[2] * mRadius;

                         if (mRollingFrictionOption){  // Rolling friction 
                             double rolling_friction             = this->GetGeometry()(0)->FastGetSolutionStepValue(ROLLING_FRICTION);
                             double rolling_friction_coeff       = rolling_friction * mRadius;

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

      void SphericParticle::DampMatrix(MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo){}

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void SphericParticle::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
      {
          KRATOS_TRY

          ElementalDofList.resize(0);

          for (unsigned int i = 0; i < GetGeometry().size(); i++){
              ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
              ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));

              if (GetGeometry().WorkingSpaceDimension() == 3){
                  ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
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
          mInitializedVariablesFlag = 1;

          KRATOS_CATCH("")

      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void SphericParticle::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo){
        
          if (rCurrentProcessInfo[PRINT_GROUP_ID] == 1){
              this->GetGeometry()[0].FastGetSolutionStepValue(EXPORT_GROUP_ID) = double(this->GetGeometry()[0].FastGetSolutionStepValue(GROUP_ID));
          }
    
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
                                                                      ParticleWeakIteratorType neighbour_iterator)
      {
        const double &other_sqrt_of_mass        = neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(SQRT_OF_MASS);
        const double &other_ln_of_restit_coeff  = neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(LN_OF_RESTITUTION_COEFF);
        const double &other_tg_of_fri_angle     = neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(PARTICLE_FRICTION);
        //const double &mLnOfRestitCoeff          = this->GetGeometry()(0)->FastGetSolutionStepValue(LN_OF_RESTITUTION_COEFF);

        double radius_sum_i                     = 1 / radius_sum;
        double equiv_radius                     = 2 * mRadius * other_radius * radius_sum_i;
        double equiv_area                       = 0.25 * M_PI * equiv_radius * equiv_radius; // 0.25 is becouse we take only the half of the equivalent radius, corresponding to the case of one ball with radius Requivalent and other = radius 0.
        double equiv_mass                       = mSqrtOfRealMass * other_sqrt_of_mass;
        double equiv_ln_of_restit_coeff;
        double aux_norm_to_tang;

        // Globally defined parameters

        
        double equiv_young;
        double equiv_poisson;
        double corrected_area               = equiv_area;

        const double &other_young       = neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(YOUNG_MODULUS);
        const double &other_poisson     = neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(POISSON_RATIO);
        double effective_radius         = 0.5 * equiv_radius;

        switch (mElasticityType){ //  0 ---linear compression & tension ; 1 --- Hertzian (non-linear compression, linear tension)
        
            case 0:
                equiv_young                     = 2 * mYoung * other_young / (mYoung + other_young);
                equiv_poisson                   = 2 * mPoisson * other_poisson / (mPoisson + other_poisson);
                equiv_ln_of_restit_coeff        = 0.5 * (mLnOfRestitCoeff + other_ln_of_restit_coeff);
                equiv_tg_of_fri_ang             = 0.5 * (mTgOfFrictionAngle + other_tg_of_fri_angle);

                kn                              = equiv_young * corrected_area * radius_sum_i; //M_PI * 0.5 * equiv_young * equiv_radius; //M: CANET FORMULA
                kt                              = kn / (2.0 + equiv_poisson + equiv_poisson);
                aux_norm_to_tang                = sqrt(kt / kn);
                
            break;

            case 1:
                equiv_young                     = mYoung * other_young / (other_young * (1- mPoisson * mPoisson) + mYoung * (1- other_poisson * other_poisson));
                equiv_poisson                   = 2 * mPoisson * other_poisson / (mPoisson + other_poisson);
                equiv_ln_of_restit_coeff        = 0.5 * (mLnOfRestitCoeff + other_ln_of_restit_coeff);
                equiv_tg_of_fri_ang             = 0.5 * (mTgOfFrictionAngle + other_tg_of_fri_angle);

                kn                              = (4/3) * equiv_young * sqrt(effective_radius);
                kt                              = 2.0 * kn * (1 - equiv_poisson * equiv_poisson) / ((2.0 - equiv_poisson) * (1 + equiv_poisson));
                aux_norm_to_tang                = sqrt(kt / kn); 
          
            break;

            default:
                equiv_young                     = 2 * mYoung * other_young / (mYoung + other_young);
                equiv_poisson                   = 2 * mPoisson * other_poisson / (mPoisson + other_poisson);
                equiv_ln_of_restit_coeff        = 0.5 * (mLnOfRestitCoeff + other_ln_of_restit_coeff);
                equiv_tg_of_fri_ang             = 0.5 * (mTgOfFrictionAngle + other_tg_of_fri_angle);

                kn                              = equiv_young * corrected_area * radius_sum_i; //M_PI * 0.5 * equiv_young * equiv_radius; //M: CANET FORMULA
                kt                              = kn / (2.0 + equiv_poisson + equiv_poisson);
                aux_norm_to_tang                = sqrt(kt / kn);
            break;

        }//switch


        if (mLnOfRestitCoeff > 0.0 || other_ln_of_restit_coeff > 0.0){ // Limit expressions when the restitution coefficient tends to 0. Variable lnRestitCoeff is set to 1.0 (instead of minus infinite) by the problem type.
            equiv_visco_damp_coeff_normal = 2 * sqrt(equiv_mass * kn);

        }

        else {

            equiv_visco_damp_coeff_normal = - 2 * equiv_ln_of_restit_coeff * sqrt(equiv_mass * kn / (equiv_ln_of_restit_coeff * equiv_ln_of_restit_coeff + M_PI * M_PI));
        }

        equiv_visco_damp_coeff_tangential = equiv_visco_damp_coeff_normal * aux_norm_to_tang;
		

      }
	  

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void SphericParticle::ComputeAdditionalForces(array_1d<double, 3>& contact_force, array_1d<double, 3>& contact_moment,
                                                      array_1d<double, 3>& externally_applied_force, array_1d<double, 3>& externally_applied_moment, ProcessInfo& rCurrentProcessInfo)
      {

        const array_1d<double,3>& gravity = rCurrentProcessInfo[GRAVITY];

        externally_applied_force[0] = mRealMass * gravity[0];
        externally_applied_force[1] = mRealMass * gravity[1];
        externally_applied_force[2] = mRealMass * gravity[2];

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

      void SphericParticle::AddUpForcesAndProject(double LocalCoordSystem[3][3],
                                                  double LocalContactForce[3],
                                                  double LocalElasticContactForce[3],
                                                  double GlobalContactForce[3],
                                                  double GlobalElasticContactForce[3],
                                                  double ViscoDampingLocalContactForce[3],
                                                  double ViscoDampingGlobalContactForce[3],
                                                  array_1d<double, 3> &rContactForce,
                                                  array_1d<double, 3> &rElasticForce,
                                                  const double &i_neighbour_count)
      {
          for (unsigned int index = 0; index < 3; index++)
          {
              LocalContactForce[index] = LocalElasticContactForce[index] + ViscoDampingLocalContactForce[index];
          }

          GeometryFunctions::VectorLocal2Global(LocalCoordSystem, LocalElasticContactForce, GlobalElasticContactForce);
          GeometryFunctions::VectorLocal2Global(LocalCoordSystem, ViscoDampingLocalContactForce, ViscoDampingGlobalContactForce);
          GeometryFunctions::VectorLocal2Global(LocalCoordSystem, LocalContactForce, GlobalContactForce);

          // Saving contact forces (We need to, since tangential elastic force is history-dependent)
          mOldNeighbourContactForces[i_neighbour_count][0] = GlobalElasticContactForce[0];
          mOldNeighbourContactForces[i_neighbour_count][1] = GlobalElasticContactForce[1];
          mOldNeighbourContactForces[i_neighbour_count][2] = GlobalElasticContactForce[2];

          rContactForce[0] += GlobalContactForce[0];
          rContactForce[1] += GlobalContactForce[1];
          rContactForce[2] += GlobalContactForce[2];
  
          rElasticForce[0] += GlobalElasticContactForce[0];
          rElasticForce[1] += GlobalElasticContactForce[1];
          rElasticForce[2] += GlobalElasticContactForce[2];
      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************
      
      void SphericParticle::MemberDeclarationFirstStep(const ProcessInfo& r_process_info)

      {

          // Paso al nodo la id del elemento cuando inicializo al mismo
          
          if (!mInitializedVariablesFlag){

               if (r_process_info[PRINT_EXPORT_ID] == 1){
                   this->GetGeometry()(0)->FastGetSolutionStepValue(EXPORT_ID) = double(this->Id());
               }
                          

              mDampType                      = r_process_info[DAMP_TYPE];
              mElasticityType                = r_process_info[FORCE_CALCULATION_TYPE];
              mRotationOption                = r_process_info[ROTATION_OPTION]; //M:  it's 1/0, should be a boolean
              mRollingFrictionOption         = r_process_info[ROLLING_FRICTION_OPTION];
              mCriticalTimeOption            = r_process_info[CRITICAL_TIME_OPTION];
             
              mLimitSurfaceOption            = r_process_info[LIMIT_SURFACE_OPTION];
              mLimitCylinderOption           = r_process_info[LIMIT_CYLINDER_OPTION];
              
              mpTimeStep                     =  &(r_process_info[TIME_STEPS]); // reference. ****
              mpActivateSearch               =  &(r_process_info[ACTIVATE_SEARCH]);

              if (mRollingFrictionOption){
                  mRollingFriction           = this->GetGeometry()(0)->FastGetSolutionStepValue(ROLLING_FRICTION);
              }

              AdditionalMemberDeclarationFirstStep(r_process_info);
          }

      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************
      
      void SphericParticle::NormalForceCalculation(double LocalElasticContactForce[3], double kn, double indentation)
              
      {

         switch (mElasticityType){ //  0 ---linear compression & tension ; 1 --- Hertzian (non-linear compression, linear tension)
              case 0:
                  LocalElasticContactForce[2] = kn * indentation;

              break;

              case 1:

                  LocalElasticContactForce[2] = kn * pow(indentation, 1.5);
                  
                 
              break;

              default:

                  LocalElasticContactForce[2] = kn * indentation;

              break;

          }//switch

          
      } //NormalForceCalculation
      
      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void SphericParticle::TangentialForceCalculation(double LocalElasticContactForce[3], double LocalDeltDisp[3], const double& kt, const double& equiv_tg_of_fri_ang, bool& sliding)
      {

        LocalElasticContactForce[0] += - kt * LocalDeltDisp[0];  // 0: first tangential
        LocalElasticContactForce[1] += - kt * LocalDeltDisp[1];  // 1: second tangential

        double ShearForceNow = sqrt(LocalElasticContactForce[0] * LocalElasticContactForce[0] + LocalElasticContactForce[1] * LocalElasticContactForce[1]);
        double Frictional_ShearForceMax = equiv_tg_of_fri_ang * LocalElasticContactForce[2];

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

      void SphericParticle::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo)
      {
          
        KRATOS_TRY

          if (rVariable == PARTICLE_ID){
              Output = mParticleId; // (NOT YET ACTIVE!!)
              return;
          }
          //CRITICAL DELTA CALCULATION

          if (rVariable == DELTA_TIME){
 
              double mass  = mRealMass;
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

                  double K = M_PI * mYoung * mRadius; //M. Error, should be the same that the local definition.

                  Output = 0.34 * sqrt(mass / K);

                  if (mRotationOption){
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
          
          if (rVariable == CALCULATE_COMPUTE_NEW_NEIGHBOURS_HISTORICAL_DATA){
             ComputeNewNeighboursHistoricalData();
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


}  // namespace Kratos.

