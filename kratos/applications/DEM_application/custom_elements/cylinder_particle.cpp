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
#include "cylinder_particle.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "DEM_application.h"


namespace Kratos
{
     // using namespace GeometryFunctions;

      CylinderParticle::CylinderParticle()
      : DiscreteElement(){mInitializedVariablesFlag = 0;}

      CylinderParticle::CylinderParticle(IndexType NewId, GeometryType::Pointer pGeometry)
      : DiscreteElement(NewId, pGeometry){mInitializedVariablesFlag = 0;}

      CylinderParticle::CylinderParticle(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
      : DiscreteElement(NewId, pGeometry, pProperties){mInitializedVariablesFlag = 0;}

      CylinderParticle::CylinderParticle(IndexType NewId, NodesArrayType const& ThisNodes)
      : DiscreteElement(NewId, ThisNodes){mInitializedVariablesFlag = 0;}

      Element::Pointer CylinderParticle::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
      {
           return DiscreteElement::Pointer(new CylinderParticle(NewId, GetGeometry().Create(ThisNodes), pProperties));

      }

      /// Destructor.
      CylinderParticle::~CylinderParticle(){}

      void CylinderParticle::Initialize()
      {
          KRATOS_TRY

          mRadius                   = GetGeometry()(0)->FastGetSolutionStepValue(RADIUS);
          mYoung                    = GetGeometry()(0)->FastGetSolutionStepValue(YOUNG_MODULUS);         
          mPoisson                  = GetGeometry()(0)->FastGetSolutionStepValue(POISSON_RATIO);
          mTgOfFrictionAngle        = GetGeometry()(0)->FastGetSolutionStepValue(PARTICLE_FRICTION);
          mLnOfRestitCoeff          = GetGeometry()(0)->FastGetSolutionStepValue(LN_OF_RESTITUTION_COEFF);
          double& density           = GetGeometry()(0)->FastGetSolutionStepValue(PARTICLE_DENSITY);
          double& mass              = GetGeometry()(0)->FastGetSolutionStepValue(NODAL_MASS);
          double& sqrt_of_mass      = GetGeometry()(0)->FastGetSolutionStepValue(SQRT_OF_MASS);
          double& moment_of_inertia = GetGeometry()(0)->FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA);

          mass                      = M_PI * density * mRadius * mRadius * 1.0;
          sqrt_of_mass              = sqrt(mass);
          moment_of_inertia         = 0.5 * mass * mRadius * mRadius;
          mRealMass                 = mass;          
          mSqrtOfRealMass           = sqrt_of_mass;
          mMomentOfInertia          = moment_of_inertia;

          CustomInitialize();

          KRATOS_CATCH( "" )
      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void CylinderParticle::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
      {
          KRATOS_TRY

          const array_1d<double, 3>& gravity  = rCurrentProcessInfo[GRAVITY];

          //array_1d<double, 3>& applied_force  = this->GetGeometry()[0].GetSolutionStepValue(APPLIED_FORCE); //MSI: canviar nom a USER_DEFINED_FORCE

          array_1d<double, 3> contact_force;
          array_1d<double, 3> contact_moment;
          array_1d<double, 3> initial_rotation_moment;
          array_1d<double, 3> max_rotation_moment;

          contact_force[0]  = 0.0;
          contact_force[1]  = 0.0;
          contact_force[2]  = 0.0;
          contact_moment[0] = 0.0;
          contact_moment[1] = 0.0;
          contact_moment[2] = 0.0;
          initial_rotation_moment[0]  = 0.0;
          initial_rotation_moment[1]  = 0.0;
          initial_rotation_moment[2]  = 0.0;
          max_rotation_moment[0] = 0.0;
          max_rotation_moment[1] = 0.0;
          max_rotation_moment[2] = 0.0;  

          ComputeNewNeighboursHistoricalData();

          ComputeBallToBallContactForce(   contact_force, contact_moment, initial_rotation_moment, max_rotation_moment, rCurrentProcessInfo); //MSI: processInfo will be eliminated since all variables will be member

          if (mLimitSurfaceOption){
              ComputeBallToSurfaceContactForce(contact_force, contact_moment, initial_rotation_moment, max_rotation_moment, rCurrentProcessInfo); //MSI: eliminate processInfo
          }

          CustomCalculateRightHandSide(contact_force, contact_moment);

          rRightHandSideVector[0] = contact_force[0] + mRealMass * gravity[0];// + applied_force[0];
          rRightHandSideVector[1] = contact_force[1] + mRealMass * gravity[1];// + applied_force[1];
          rRightHandSideVector[2] = 0 ;
          rRightHandSideVector[3] = contact_moment[0];
          rRightHandSideVector[4] = contact_moment[1];
          rRightHandSideVector[5] = 0;

          KRATOS_CATCH( "" )
      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void CylinderParticle::CalculateMaxIndentation(double& rCurrentMaxIndentation, const double& rTolerance)
      {

          if (rCurrentMaxIndentation > rTolerance){
              ParticleWeakVectorType& rNeighbours = this->GetValue(NEIGHBOUR_ELEMENTS);
              double& radius                      = GetGeometry()(0)->FastGetSolutionStepValue(RADIUS); // cannot use the member variable mRadius because it may not be initialized
              rCurrentMaxIndentation              = 0.0;

              for (ParticleWeakIteratorType i = rNeighbours.begin(); i != rNeighbours.end(); i++){
                  array_1d<double, 3> other_to_me_vect  = this->GetGeometry()(0)->Coordinates() - i->GetGeometry()(0)->Coordinates();
                  double &other_radius                  = i->GetGeometry()(0)->FastGetSolutionStepValue(RADIUS);
                  double distance                       = sqrt(other_to_me_vect[0] * other_to_me_vect[0] +
                                                               other_to_me_vect[1] * other_to_me_vect[1]);
                  double radius_sum                     = radius + other_radius;
                  double indentation                    = radius_sum - distance;

                  rCurrentMaxIndentation = (indentation > rCurrentMaxIndentation) ? indentation : rCurrentMaxIndentation;

              }

          }

      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

     void CylinderParticle::ComputeNewNeighboursHistoricalData()
     {

       ParticleWeakVectorType& rNeighbours  = this->GetValue(NEIGHBOUR_ELEMENTS);
       unsigned int new_size                = rNeighbours.size();
       unsigned int neighbour_counter       = 0;
       vector<int> temp_neighbours_ids(new_size);
       vector<array_1d<double, 3> > temp_neighbours_contact_forces(new_size);
       array_1d<double, 3> vector_of_zeros;
       vector_of_zeros[0]                   = 0.0;
       vector_of_zeros[1]                   = 0.0;
       vector_of_zeros[2]                   = 0.0;

       for (ParticleWeakIteratorType i = rNeighbours.begin(); i != rNeighbours.end(); i++){

          temp_neighbours_ids[neighbour_counter] = static_cast<int>(i->Id());
          temp_neighbours_contact_forces[neighbour_counter] = vector_of_zeros;

            for (unsigned int j = 0; j != mOldNeighbourIds.size(); j++)
            {
              if (static_cast<int>(i->Id()) == mOldNeighbourIds[j])
              {
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


      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void CylinderParticle::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo){}

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void CylinderParticle::MassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
      {
          rMassMatrix(0,0) = mRealMass;
      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void CylinderParticle::EvaluateDeltaDisplacement(double DeltDisp[3],
                                                       double RelVel[3],
                                                       double NormalDir[3],
                                                       double OldNormalDir[3],
                                                       double LocalCoordSystem[3][3],
                                                       double OldLocalCoordSystem[3][3],
                                                       const array_1d<double, 3> &other_to_me_vect,
                                                       const array_1d<double, 3> &vel,
                                                       const array_1d<double, 3> &delta_displ,
                                                       ParticleWeakIteratorType neighbour_iterator)
      {


          // FORMING LOCAL CORDINATES

          //Notes: Since we will normally inherit the mesh from GiD, we respect the global system X,Y,Z [0],[1],[2]
          //In the local coordinates we will define the normal direction of the contact as the [2] component!!!!!
          //the way the normal direction is defined (other_to_me_vect) compression is positive

          NormalDir[0] = other_to_me_vect[0];
          NormalDir[1] = other_to_me_vect[1];
          NormalDir[2] = other_to_me_vect[2];

          GeometryFunctions::ComputeContactLocalCoordSystem(NormalDir, LocalCoordSystem); //new Local Coord System

          // FORMING OLD LOCAL CORDINATES

          array_1d<double,3> old_coord_target     = this->GetGeometry()(0)->GetInitialPosition() + this->GetGeometry()(0)->GetSolutionStepValue(DISPLACEMENT,1);
          array_1d<double,3> old_coord_neigh      = neighbour_iterator->GetGeometry()(0)->GetInitialPosition()+neighbour_iterator->GetGeometry()(0)->GetSolutionStepValue(DISPLACEMENT,1);
          array_1d<double,3> Old_other_to_me_vect = old_coord_target - old_coord_neigh;

          OldNormalDir[0] = Old_other_to_me_vect[0];   // M. this way the compresion is positive.
          OldNormalDir[1] = Old_other_to_me_vect[1];
          OldNormalDir[2] = Old_other_to_me_vect[2];
          GeometryFunctions::ComputeContactLocalCoordSystem(OldNormalDir, OldLocalCoordSystem); //Old Local Coord System

          // VELOCITIES AND DISPLACEMENTS
          array_1d<double, 3 > other_vel            = neighbour_iterator->GetGeometry()(0)->GetSolutionStepValue(VELOCITY);
          array_1d<double, 3 > other_delta_displ    = neighbour_iterator->GetGeometry()(0)->GetSolutionStepValue(DELTA_DISPLACEMENT);

          RelVel[0] = (vel[0] - other_vel[0]);
          RelVel[1] = (vel[1] - other_vel[1]);
          RelVel[2] = 0;

          //DeltDisp in global cordinates
          DeltDisp[0] = (delta_displ[0] - other_delta_displ[0]);
          DeltDisp[1] = (delta_displ[1] - other_delta_displ[1]);
          DeltDisp[2] = 0;
      }

      void CylinderParticle::DisplacementDueToRotation(double DeltDisp[3],
                                                      double OldNormalDir[3],
                                                      double OldLocalCoordSystem[3][3],
                                                      const double &other_radius,
                                                      const double &dt,
                                                      const array_1d<double, 3> &ang_vel,
                                                      ParticleWeakIteratorType neighbour_iterator)
      {
          if (mRotationOption == 1){
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
              dRotaDisp[2] = 0;
              // Contribution of the rotation velocity
              DeltDisp[0] += dRotaDisp[0] * dt;
              DeltDisp[1] += dRotaDisp[1] * dt;
              DeltDisp[2] += 0;
          }
      }

      void CylinderParticle::ComputeMoments(double LocalElasticContactForce[3],
                                           double GlobalElasticContactForce[3],
                                           double InitialRotaMoment[3],
                                           double MaxRotaMoment[3],
                                           double LocalCoordSystem[3][3],
                                           const double &other_radius,
                                           array_1d<double, 3>& rContactMoment,
                                           ParticleWeakIteratorType neighbour_iterator)
      {
          if (mRotationOption){

              double MA[3]         = {0.0};
              double RotaMoment[3] = {0.0};

              GeometryFunctions::CrossProduct(LocalCoordSystem[2], GlobalElasticContactForce, MA);

              RotaMoment[0] -= MA[0] * mRadius;
              RotaMoment[1] -= MA[1] * mRadius;
              RotaMoment[2] -= MA[2] * mRadius;

              // ROLLING FRICTION

              if (mRotationDampType == 2){  // rolling friccion type
                  double rolling_friction_coeff            = mRollingFriction * mRadius;
                  double equiv_rolling_friction_coeff;

                  if (mUniformMaterialOption){
                      equiv_rolling_friction_coeff         = rolling_friction_coeff;
                  }

                  else {
                      const double& other_rolling_friction = neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(ROLLING_FRICTION);
                      double other_rolling_friction_coeff  = other_rolling_friction * other_radius;
                      equiv_rolling_friction_coeff         = 2 * rolling_friction_coeff * other_rolling_friction_coeff / (rolling_friction_coeff + other_rolling_friction_coeff);
                  }

                  if (equiv_rolling_friction_coeff != 0.0){
                      double CoordSystemMoment[3] = {0.0};
                      double MR[3]                = {0.0};
                      double NormalForce[3]       = {0.0};

                      MaxRotaMoment[0] += RotaMoment[0];
                      MaxRotaMoment[1] += RotaMoment[1];
                      MaxRotaMoment[2] += RotaMoment[2];

                      NormalForce[0] = LocalCoordSystem[2][0] * fabs(LocalElasticContactForce[2]);
                      NormalForce[1] = LocalCoordSystem[2][1] * fabs(LocalElasticContactForce[2]);
                      NormalForce[2] = 0;

                      GeometryFunctions::CrossProduct(LocalCoordSystem[2], MaxRotaMoment, CoordSystemMoment);
                      double det_coor_sys_moment_i = 1 / sqrt(CoordSystemMoment[0] * CoordSystemMoment[0] + CoordSystemMoment[1] * CoordSystemMoment[1] + CoordSystemMoment[2] * CoordSystemMoment[2]);

                      CoordSystemMoment[0] *= det_coor_sys_moment_i;
                      CoordSystemMoment[1] *= det_coor_sys_moment_i;
                      CoordSystemMoment[2] *= det_coor_sys_moment_i;

                      GeometryFunctions::CrossProduct(NormalForce, CoordSystemMoment, MR);
                      double det_MR = sqrt(MR[0] * MR[0] + MR[1] * MR[1] + MR[2] * MR[2]);
                      double MR_now = det_MR * equiv_rolling_friction_coeff;
                      double MR_max = sqrt(MaxRotaMoment[0] * MaxRotaMoment[0] + MaxRotaMoment[1] * MaxRotaMoment[1] + MaxRotaMoment[2] * MaxRotaMoment[2]);

                      if (MR_max > MR_now){
                          RotaMoment[0]    += MR[0] * equiv_rolling_friction_coeff;
                          RotaMoment[1]    += MR[1] * equiv_rolling_friction_coeff;
                          RotaMoment[2]    += MR[2] * equiv_rolling_friction_coeff;

                          MaxRotaMoment[0] += MR[0] * equiv_rolling_friction_coeff;
                          MaxRotaMoment[1] += MR[1] * equiv_rolling_friction_coeff;
                          MaxRotaMoment[2] += MR[2] * equiv_rolling_friction_coeff;
                      }

                      else {
                          RotaMoment[0]     = - InitialRotaMoment[0];
                          RotaMoment[1]     = - InitialRotaMoment[1];
                          RotaMoment[2]     = - InitialRotaMoment[2];
                      }

                    } // if (equiv_rolling_friction_coeff != 0.0)

                } // if (mRotationDampType == 2)

              rContactMoment[0] += RotaMoment[0];
              rContactMoment[1] += RotaMoment[1];
              rContactMoment[2] += RotaMoment[2];

          } // if (mRotationOption)
      }

      void CylinderParticle::ComputeBallToBallContactForce(array_1d<double, 3>& rContactForce,
                                                          array_1d<double, 3>& rContactMoment,
                                                          array_1d<double, 3>& rInitialRotaMoment,
                                                          array_1d<double, 3>& rMaxRotaMoment,
                                                          ProcessInfo& rCurrentProcessInfo)
      {
          KRATOS_TRY

          
          // PROCESS INFO

          ParticleWeakVectorType& mrNeighbours         = this->GetValue(NEIGHBOUR_ELEMENTS);
          //vector<double>& r_VectorContactInitialDelta  = this->GetValue(PARTICLE_CONTACT_DELTA);  //MSI: must be changed in the same fashion as contactforces

          double dt = rCurrentProcessInfo[DELTA_TIME];
          double dt_i = 1 / dt;

          //INITIALIZATIONS

          const array_1d<double, 3>& vel         = this->GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY);
          const array_1d<double, 3>& delta_displ = this->GetGeometry()(0)->FastGetSolutionStepValue(DELTA_DISPLACEMENT);
          const array_1d<double, 3>& ang_vel     = this->GetGeometry()(0)->FastGetSolutionStepValue(ANGULAR_VELOCITY);
          double RotaAcc[3]                      = {0.0};
          double InitialRotaMoment[3]            = {0.0};
          double MaxRotaMoment[3]                = {0.0};

          if (mRotationOption){
              RotaAcc[0]                         = ang_vel[0] * dt_i;
              RotaAcc[1]                         = ang_vel[1] * dt_i;
              RotaAcc[2]                         = ang_vel[2] * dt_i;

              InitialRotaMoment[0]               = RotaAcc[0] * mMomentOfInertia;
              InitialRotaMoment[1]               = RotaAcc[1] * mMomentOfInertia;
              InitialRotaMoment[2]               = RotaAcc[2] * mMomentOfInertia;

              MaxRotaMoment[0]                   = InitialRotaMoment[0];
              MaxRotaMoment[1]                   = InitialRotaMoment[1];
              MaxRotaMoment[2]                   = InitialRotaMoment[2];
          }


          //LOOP OVER NEIGHBOURS BEGINS

          size_t i_neighbour_count = 0;

          for (ParticleWeakIteratorType neighbour_iterator = mrNeighbours.begin();
              neighbour_iterator != mrNeighbours.end(); neighbour_iterator++){

              // BASIC CALCULATIONS

              array_1d<double, 3> other_to_me_vect    = this->GetGeometry()(0)->Coordinates() - neighbour_iterator->GetGeometry()(0)->Coordinates();
              const double &other_radius              = neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(RADIUS);
              const double &other_sqrt_of_mass        = neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(SQRT_OF_MASS);
              const double &other_ln_of_restit_coeff  = neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(LN_OF_RESTITUTION_COEFF);
              double distance                         = sqrt(other_to_me_vect[0] * other_to_me_vect[0] +
                                                             other_to_me_vect[1] * other_to_me_vect[1] +
                                                             other_to_me_vect[2] * other_to_me_vect[2]);
              double radius_sum                       = mRadius + other_radius;
              double radius_sum_i                     = 1 / radius_sum;
              double equiv_radius                     = 2 * mRadius * other_radius * radius_sum_i;
              double indentation                      = radius_sum - distance;
              double equiv_area                       = 0.25 * M_PI * equiv_radius * equiv_radius; // 0.25 is becouse we take only the half of the equivalent radius, corresponding to the case of one ball with radius Requivalent and other = radius 0.
              double corrected_area                   = equiv_area;
              double equiv_mass                       = mSqrtOfRealMass * other_sqrt_of_mass;

              double equiv_young;
              double equiv_poisson;
              double equiv_visco_damp_coeff_normal;
              double equiv_visco_damp_coeff_tangential;
              double equiv_ln_of_restit_coeff;
              double kn;
              double kt;
              double equiv_tg_of_fri_ang;

              double DeltDisp[3]                      = {0.0};
              double RelVel[3]                        = {0.0};

              double NormalDir[3]                     = {0.0};
              double OldNormalDir[3]                  = {0.0};

              double LocalCoordSystem[3][3]           = {{0.0}, {0.0}, {0.0}};
              double OldLocalCoordSystem[3][3]        = {{0.0}, {0.0}, {0.0}};

              bool sliding = false;
              
              if (mUniformMaterialOption){
                  equiv_radius                       = mRadius;
                  equiv_young                        = mYoung;
                  equiv_poisson                      = mPoisson;
                  equiv_ln_of_restit_coeff           = mLnOfRestitCoeff;
                  equiv_tg_of_fri_ang                = mTgOfFrictionAngle;
              }

              else {
                  // Getting neighbour properties
                  const double &other_young           = neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(YOUNG_MODULUS);
                  const double &other_poisson         = neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(POISSON_RATIO);
                  const double &other_tg_of_fri_angle = neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(PARTICLE_FRICTION);

                  equiv_young                         = 2 * mYoung * other_young / (mYoung + other_young);
                  equiv_poisson                       = 2 * mPoisson * other_poisson / (mPoisson + other_poisson);
                  equiv_ln_of_restit_coeff            = 0.5 * (mLnOfRestitCoeff + other_ln_of_restit_coeff);
                  equiv_tg_of_fri_ang                 = 0.5 * (mTgOfFrictionAngle + other_tg_of_fri_angle);
              }

              // Globally defined parameters
                  double aux_norm_to_tang;

              if (mGlobalVariablesOption){
                  kn                                  = mGlobalKn;
                  kt                                  = mGlobalKt;
                  aux_norm_to_tang                    = mGlobalAuxNormToTang;
              }

              else {
                  kn                                  = mMagicFactor * equiv_young * corrected_area * radius_sum_i; //M_PI * 0.5 * equiv_young * equiv_radius; //M: CANET FORMULA
                  kt                                  = kn / (2.0 + equiv_poisson + equiv_poisson);
                  aux_norm_to_tang                    = sqrt(kt / kn);
              }

              // Historical minimun K for the critical time

              if (mCriticalTimeOption){
                  double historic = rCurrentProcessInfo[HISTORICAL_MIN_K];

                  if ((kn < historic) || (kt < historic)){
                      historic = fmin(kn, kt);
                  }

              }

            
              if (mLnOfRestitCoeff > 0.0 || other_ln_of_restit_coeff > 0.0){
                  equiv_visco_damp_coeff_normal     = 2 * sqrt(equiv_mass * kn);
                  equiv_visco_damp_coeff_tangential = equiv_visco_damp_coeff_normal * aux_norm_to_tang; // 2 * sqrt(equiv_mass * kt);
              }

              else {
                  equiv_visco_damp_coeff_normal     = - 2 * equiv_ln_of_restit_coeff * sqrt(equiv_mass * kn / (equiv_ln_of_restit_coeff * equiv_ln_of_restit_coeff + M_PI * M_PI));
                  //equiv_visco_damp_coeff_normal     = - 2 * log(equiv_restitution_coeff) * sqrt(equiv_mass * kn) / sqrt((log(equiv_restitution_coeff) * log(equiv_restitution_coeff)) + (M_PI * M_PI));
                  equiv_visco_damp_coeff_tangential = equiv_visco_damp_coeff_normal * aux_norm_to_tang; //= -(2 * log(equiv_restitution_coeff) * sqrt(equiv_mass * kt)) / (sqrt((log(equiv_restitution_coeff) * log(equiv_restitution_coeff)) + (M_PI * M_PI)));
              }

              EvaluateDeltaDisplacement(DeltDisp, RelVel, NormalDir, OldNormalDir, LocalCoordSystem, OldLocalCoordSystem, other_to_me_vect, vel, delta_displ, neighbour_iterator);
              DisplacementDueToRotation(DeltDisp, OldNormalDir, OldLocalCoordSystem, other_radius, dt, ang_vel, neighbour_iterator);

              double LocalDeltDisp[3]               = {0.0};
              double LocalElasticContactForce[3]    = {0.0};
              double GlobalElasticContactForce[3]   = {0.0};
              double LocalRelVel[3]                 = {0.0};

              GlobalElasticContactForce[0]          = mOldNeighbourContactForces[i_neighbour_count][0];   //M:aqui tenim guardades les del neighbour calculator.
              GlobalElasticContactForce[1]          = mOldNeighbourContactForces[i_neighbour_count][1];
              GlobalElasticContactForce[2]          = mOldNeighbourContactForces[i_neighbour_count][2];

              GeometryFunctions::VectorGlobal2Local(OldLocalCoordSystem, GlobalElasticContactForce, LocalElasticContactForce); //we recover this way the old local forces projected in the new coordinates in the way they were in the old ones; Now they will be increased if its the necessary
              GeometryFunctions::VectorGlobal2Local(OldLocalCoordSystem, DeltDisp, LocalDeltDisp);
              GeometryFunctions::VectorGlobal2Local(OldLocalCoordSystem, RelVel, LocalRelVel);

              // TRANSLATION FORCES

              if (indentation > 0.0){
        
              NormalForceCalculation(LocalElasticContactForce, kn, indentation, mElasticityType);
                            
                  // TANGENTIAL FORCE
                  // Incremental calculation. YADE develops a complicated "absolute method"

                  LocalElasticContactForce[0] += - kt * LocalDeltDisp[0];  // 0: first tangential
                  LocalElasticContactForce[1] += - kt * LocalDeltDisp[1];  // 1: second tangential

                  double ShearForceNow = sqrt(LocalElasticContactForce[0] * LocalElasticContactForce[0]
                                            + LocalElasticContactForce[1] * LocalElasticContactForce[1]);
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

              // VISCODAMPING (applyied locally)
              double ViscoDampingLocalContactForce[3] = {0.0};

              CalculateViscoDamping(LocalRelVel, ViscoDampingLocalContactForce, indentation, equiv_visco_damp_coeff_normal, equiv_visco_damp_coeff_tangential, sliding);

              // Transforming to global forces and adding up
              double LocalContactForce[3]              = {0.0};
              double ViscoDampingGlobalContactForce[3] = {0.0};
              double GlobalContactForce[3]             = {0.0};

              AddUpForcesAndProject(LocalCoordSystem, mOldNeighbourContactForces, LocalContactForce, LocalElasticContactForce, GlobalContactForce, GlobalElasticContactForce, ViscoDampingLocalContactForce, ViscoDampingGlobalContactForce, rContactForce, i_neighbour_count);

              // ROTATION FORCES
              ComputeMoments(LocalElasticContactForce, GlobalElasticContactForce, InitialRotaMoment, MaxRotaMoment, LocalCoordSystem, other_radius, rContactMoment, neighbour_iterator);

              i_neighbour_count++;

          }// for each neighbour
          
          rInitialRotaMoment [0] = InitialRotaMoment [0];
          rInitialRotaMoment [1] = InitialRotaMoment [1];
          rInitialRotaMoment [2] = InitialRotaMoment [2];
          rMaxRotaMoment [0]     = MaxRotaMoment [0];
          rMaxRotaMoment [1]     = MaxRotaMoment [1];
          rMaxRotaMoment [2]     = MaxRotaMoment [2];

          KRATOS_CATCH("")

      }// ComputeBallToBallContactForce

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void CylinderParticle::ComputeBallToSurfaceContactForce(array_1d<double, 3>& rContactForce,
                                                             array_1d<double, 3>& rContactMoment,
                                                             array_1d<double, 3>& rInitialRotaMoment,
                                                             array_1d<double, 3>& rMaxRotaMoment,
                                                             ProcessInfo& rCurrentProcessInfo)
      {
          KRATOS_TRY

          // CONTACT WITH A PLANE

           if (rCurrentProcessInfo[LIMIT_SURFACE_OPTION] != 0){
              double dt                            = rCurrentProcessInfo[DELTA_TIME];

             // INITIALIZATIONS

             const array_1d<double, 3> vel         = this->GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY);
             const array_1d<double, 3> delta_displ = this->GetGeometry()(0)->FastGetSolutionStepValue(DELTA_DISPLACEMENT);
             const array_1d<double, 3> ang_vel     = this->GetGeometry()(0)->FastGetSolutionStepValue(ANGULAR_VELOCITY);
             double InitialRotaMoment[3]           = {0.0};
             double MaxRotaMoment[3]               = {0.0};
             double visco_damp_coeff_normal;
             double visco_damp_coeff_tangential;

             InitialRotaMoment [0] = rInitialRotaMoment [0];
             InitialRotaMoment [1] = rInitialRotaMoment [1];
             InitialRotaMoment [2] = rInitialRotaMoment [2];
             MaxRotaMoment [0]     = rMaxRotaMoment [0];
             MaxRotaMoment [1]     = rMaxRotaMoment [1];
             MaxRotaMoment [2]     = rMaxRotaMoment [2];

             // SLIDING

             bool sliding = false;

             // BASIC CALCULATIONS

             const array_1d<double,3>& surface_normal_dir = rCurrentProcessInfo[SURFACE_NORMAL_DIR];
             const array_1d<double,3>& surface_point_coor = rCurrentProcessInfo[SURFACE_POINT_COOR];
             array_1d<double, 3> & GlobalSurfContactForce = this->GetValue(PARTICLE_SURFACE_CONTACT_FORCES);
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

             if (indentation > 0.0){
                 // MACRO PARAMETERS
                 double kn = M_PI * 0.5 * mYoung * mRadius; //M_PI * 0.5 * equiv_young * equiv_radius; //M: CANET FORMULA
                 double kt = kn / (2.0 * (1.0 + mPoisson));

                 if (mGlobalVariablesOption){ //globally defined parameters
                     kn = mGlobalKn;
                     kt = mGlobalKt;
                 }

                 // Historical minimun K for the critical time:
                 if (mCriticalTimeOption){
                     double historic = rCurrentProcessInfo[HISTORICAL_MIN_K];

                     if ((kn < historic) || (kt < historic)){
                         historic = fmin(kn, kt);
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

                 double NormalDir[3]            = {0.0};
                 double LocalCoordSystem[3][3]  = {{0.0}, {0.0}, {0.0}};
                 double norm_surface_normal_dir = sqrt(surface_normal_dir[0] * surface_normal_dir[0] + surface_normal_dir[1] * surface_normal_dir[1] + surface_normal_dir[2] * surface_normal_dir[2]);
                 NormalDir[0] = surface_normal_dir[0] / norm_surface_normal_dir;
                 NormalDir[1] = surface_normal_dir[1] / norm_surface_normal_dir;
                 NormalDir[2] = surface_normal_dir[2] / norm_surface_normal_dir;

                 GeometryFunctions::ComputeContactLocalCoordSystem(NormalDir, LocalCoordSystem); //new Local Coord System

                 // VELOCITIES AND DISPLACEMENTS

                 double DeltDisp[3] = {0.0};
                 double RelVel  [3] = {0.0};

                 RelVel[0] = vel[0];
                 RelVel[1] = vel[1];
                 RelVel[2] = 0;

                 // DeltDisp in global cordinates

                 DeltDisp[0] = delta_displ[0];
                 DeltDisp[1] = delta_displ[1];
                 DeltDisp[2] = 0;

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

                 double dyn_friction_angle =  rCurrentProcessInfo[SURFACE_FRICC] * M_PI / 180;
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

                     if (sliding == false && (mDampType == 1 || mDampType == 11)){ //only applied when no sliding to help to the regularized friccion law or the spring convergence
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

                 this->GetValue(PARTICLE_SURFACE_CONTACT_FORCES)[0] = GlobalElasticContactForce[0];
                 this->GetValue(PARTICLE_SURFACE_CONTACT_FORCES)[1] = GlobalElasticContactForce[1];
                 this->GetValue(PARTICLE_SURFACE_CONTACT_FORCES)[2] = GlobalElasticContactForce[2];

                 if (mRotationOption){
                     double MA[3]                     = {0.0};
                     double RotaMoment[3]             = {0.0};

                     GeometryFunctions::CrossProduct(LocalCoordSystem[2], GlobalElasticContactForce, MA);

                     RotaMoment[0] -= MA[0] * mRadius;
                     RotaMoment[1] -= MA[1] * mRadius;
                     RotaMoment[2] -= MA[2] * mRadius;

                     if (mRotationDampType == 2){  // Rolling friccion type
                         double rolling_friction             = this->GetGeometry()(0)->FastGetSolutionStepValue(ROLLING_FRICTION);
                         double rolling_friction_coeff       = rolling_friction * mRadius;

                         if (rolling_friction_coeff != 0.0){
                             double CoordSystemMoment[3] = {0.0};
                             double MR[3]                = {0.0};
                             double NormalForce[3]       = {0.0};

                             MaxRotaMoment[0] += RotaMoment[0];
                             MaxRotaMoment[1] += RotaMoment[1];
                             MaxRotaMoment[2] += RotaMoment[2];

                             NormalForce[0] = LocalCoordSystem[2][0] * fabs(LocalElasticContactForce[2]);
                             NormalForce[1] = LocalCoordSystem[2][1] * fabs(LocalElasticContactForce[2]);
                             NormalForce[2] = LocalCoordSystem[2][2] * fabs(LocalElasticContactForce[2]);

                             GeometryFunctions::CrossProduct(LocalCoordSystem[2], MaxRotaMoment, CoordSystemMoment);
                             double det_coor_sys_moment_i = 1 / sqrt(CoordSystemMoment[0] * CoordSystemMoment[0] + CoordSystemMoment[1] * CoordSystemMoment[1] + CoordSystemMoment[2] * CoordSystemMoment[2]);

                             CoordSystemMoment[0] *= det_coor_sys_moment_i;
                             CoordSystemMoment[1] *= det_coor_sys_moment_i;
                             CoordSystemMoment[2] *= det_coor_sys_moment_i;

                             GeometryFunctions::CrossProduct(NormalForce, CoordSystemMoment, MR);
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

                             else {
                                 RotaMoment[0]     = -InitialRotaMoment[0];
                                 RotaMoment[1]     = -InitialRotaMoment[1];
                                 RotaMoment[2]     = -InitialRotaMoment[2];
                             }

                         } // if (RollingFrictionCoeff != 0.0)

                     } // if (mRotationDampType == 2)

                     rContactMoment[0] += RotaMoment[0];
                     rContactMoment[1] += RotaMoment[1];
                     rContactMoment[2] += RotaMoment[2];

                     rInitialRotaMoment [0] = InitialRotaMoment [0];
                     rInitialRotaMoment [1] = InitialRotaMoment [1];
                     rInitialRotaMoment [2] = InitialRotaMoment [2];
                     rMaxRotaMoment [0]     = MaxRotaMoment [0];
                     rMaxRotaMoment [1]     = MaxRotaMoment [1];
                     rMaxRotaMoment [2]     = MaxRotaMoment [2];

                 } //if (mRotationOption)

             }  //if (indentation >= 0.0)

         } //if (LIMIT_SURFACE_OPTION)

         KRATOS_CATCH("")
      }//ComputeBallToSurfaceContactForce

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void CylinderParticle::DampMatrix(MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo){}

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void CylinderParticle::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
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

      void CylinderParticle::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
      {
          KRATOS_TRY

          MemberDeclarationFirstStep(rCurrentProcessInfo);
          mInitializedVariablesFlag = 1;

          KRATOS_CATCH("")

      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void CylinderParticle::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo){}

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void CylinderParticle::CustomInitialize()
      {

      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void CylinderParticle::CustomCalculateRightHandSide(array_1d<double, 3>& contact_force, array_1d<double, 3>& contact_moment)
      {
      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void CylinderParticle::CalculateViscoDamping(double LocalRelVel[3],
                                                  double ViscoDampingLocalContactForce[3],
                                                  double indentation,
                                                  double equiv_visco_damp_coeff_normal,
                                                  double equiv_visco_damp_coeff_tangential,
                                                  bool sliding)
      {


          //*** The compbrobation is component-wise since localContactForce and RelVel have in principle no relationship.
          // The visco force can be higher than the contact force only if they go to the same direction. (in my opinion)
          // But in oposite direction the visco damping can't overpass the force...

          if ((mDampType > 0) && (indentation > 0.0)){

              if (mDampType == 11 || mDampType == 10){
                  ViscoDampingLocalContactForce[2] = - equiv_visco_damp_coeff_normal * LocalRelVel[2];
              }

              if (sliding == false && (mDampType == 1 || mDampType == 11)){ //only applied when no sliding to help to the regularized friccion law or the spring convergence
                  ViscoDampingLocalContactForce[0] = - equiv_visco_damp_coeff_tangential * LocalRelVel[0];
                  ViscoDampingLocalContactForce[1] = - equiv_visco_damp_coeff_tangential * LocalRelVel[1];
              }
          }
      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void CylinderParticle::AddUpForcesAndProject(double LocalCoordSystem[3][3],
                                 VectorArray3Double &mOldNeighbourContactForces,
                                 double LocalContactForce[3],
                                 double LocalElasticContactForce[3],
                                 double GlobalContactForce[3],
                                 double GlobalElasticContactForce[3],
                                 double ViscoDampingLocalContactForce[3],
                                 double ViscoDampingGlobalContactForce[3],
                                 array_1d<double, 3> &rContactForce,
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
      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void CylinderParticle::MemberDeclarationFirstStep(ProcessInfo& rCurrentProcessInfo)
      {

          // Paso al nodo la id del elemento cuando inicializo al mismo

          if (!mInitializedVariablesFlag){

              if (rCurrentProcessInfo[INT_DUMMY_3] == 1){
                  this->GetGeometry()(0)->FastGetSolutionStepValue(EXPORT_ID) = double(this->Id());
              }

              mDampType                 = rCurrentProcessInfo[DAMP_TYPE];
              mElasticityType           = rCurrentProcessInfo[FORCE_CALCULATION_TYPE];
              mRotationOption           = rCurrentProcessInfo[ROTATION_OPTION]; //M:  it's 1/0, should be a boolean
              mRotationDampType         = rCurrentProcessInfo[ROTA_DAMP_TYPE];
              mGlobalVariablesOption    = rCurrentProcessInfo[GLOBAL_VARIABLES_OPTION]; //M:  it's 1/0, should be a boolean
              mCriticalTimeOption       = rCurrentProcessInfo[CRITICAL_TIME_OPTION];
              mUniformMaterialOption    = rCurrentProcessInfo[UNIFORM_MATERIAL_OPTION];
              mMagicFactor              = rCurrentProcessInfo[DEM_MAGIC_FACTOR];
              mGlobalKn                 = rCurrentProcessInfo[GLOBAL_KN];
              mGlobalKt                 = rCurrentProcessInfo[GLOBAL_KT];
              mLimitSurfaceOption       = rCurrentProcessInfo[LIMIT_SURFACE_OPTION];
              mRotationSpringOption     = rCurrentProcessInfo[ROTATION_SPRING_OPTION];

              if (mRotationOption){
                  mRollingFriction      = this->GetGeometry()(0)->FastGetSolutionStepValue(ROLLING_FRICTION);
              }

              if (mGlobalKn == 0.0){
                  mGlobalAuxNormToTang  = 0.0;
              }

              else {
                  mGlobalAuxNormToTang  = sqrt(mGlobalKt / mGlobalKn);;
              }

          }

      }
      
       
      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************
      
       void CylinderParticle::NormalForceCalculation(double LocalElasticContactForce[3], double kn, double indentation, int mElasticityType)
              
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

      void CylinderParticle::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo)
      {
          KRATOS_TRY

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

                  if (mGlobalVariablesOption == 1){
                      K = mGlobalKn;
                  }

                  Output = 0.34 * sqrt(mass / K);

                  if (mRotationOption == 1){
                      Output *= 0.5; //factor for critical time step when rotation is allowed.
                  }

             }

          }

          else if (rVariable == MAX_INDENTATION){
              CalculateMaxIndentation(Output, rCurrentProcessInfo[DISTANCE_TOLERANCE]);
          }

          KRATOS_CATCH("")

      }// Calculate

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void CylinderParticle::Calculate(const Variable<array_1d<double, 3> >& rVariable, array_1d<double, 3>& Output, const ProcessInfo& rCurrentProcessInfo){}
      void CylinderParticle::Calculate(const Variable<Vector >& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo){}
      void CylinderParticle::Calculate(const Variable<Matrix >& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo){}

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

}  // namespace Kratos.

