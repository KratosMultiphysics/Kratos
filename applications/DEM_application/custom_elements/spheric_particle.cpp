//
// Authors: 
// Miguel Angel Celigueta maceli@cimne.upc.edu
// Salvador Latorre latorre@cimne.upc.edu
// Miquel Santasusana msantasusana@cimne.upc.edu
// Guillermo Casas gcasas@cimne.upc.edu
//

// System includes
#include <string>
#include <iostream>
#include <cmath>

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
    : DiscreteElement(), mRealMass(0)
{
    mRadius = 0;
    mRealMass = 0;
    mStressTensor = NULL;
    mSymmStressTensor = NULL;
}

SphericParticle::SphericParticle(IndexType NewId, GeometryType::Pointer pGeometry)
    : DiscreteElement(NewId, pGeometry), mRealMass(0){
    mRadius = 0;
    mRealMass = 0;
    mStressTensor = NULL;
    mSymmStressTensor = NULL;
}

SphericParticle::SphericParticle(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : DiscreteElement(NewId, pGeometry, pProperties), mRealMass(0)
{
    mRadius = 0;
    mRealMass = 0;
    mStressTensor = NULL;
    mSymmStressTensor = NULL;
}

SphericParticle::SphericParticle(IndexType NewId, NodesArrayType const& ThisNodes)
    : DiscreteElement(NewId, ThisNodes), mRealMass(0)
{
    mRadius = 0;
    mRealMass = 0;
    mStressTensor = NULL;
    mSymmStressTensor = NULL;
}

Element::Pointer SphericParticle::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new SphericParticle(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

/// Destructor.
SphericParticle::~SphericParticle(){
    
    if(mStressTensor!=NULL) {
        delete mStressTensor;
        mStressTensor = NULL;
        delete mSymmStressTensor;
        mSymmStressTensor = NULL;
    }
}

void SphericParticle::Initialize()
{
    KRATOS_TRY

    mRadius                   = GetGeometry()[0].FastGetSolutionStepValue(RADIUS);
    double density            = GetDensity();
    double& mass              = GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS);
    mass                      = 4 *KRATOS_M_PI_3 * density * GetRadius() * GetRadius() * GetRadius();
    mRealMass                 = mass;
    
    if (this->IsNot(BLOCKED)) GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MATERIAL) = GetParticleMaterial();
    
    mClusterId = -1;

    if (this->Is(DEMFlags::HAS_ROTATION) ){
        double moment_of_inertia = 0.4 * mass * GetRadius() * GetRadius();
        GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA) = moment_of_inertia;
    }

    else {
        array_1d<double, 3>& angular_velocity = GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
        angular_velocity = ZeroVector(3);
    }

    if (GetGeometry()[0].GetDof(VELOCITY_X).IsFixed())         GetGeometry()[0].Set(DEMFlags::FIXED_VEL_X,true);
    else                                                          GetGeometry()[0].Set(DEMFlags::FIXED_VEL_X,false);
    if (GetGeometry()[0].GetDof(VELOCITY_Y).IsFixed())         GetGeometry()[0].Set(DEMFlags::FIXED_VEL_Y,true);
    else                                                          GetGeometry()[0].Set(DEMFlags::FIXED_VEL_Y,false);
    if (GetGeometry()[0].GetDof(VELOCITY_Z).IsFixed())         GetGeometry()[0].Set(DEMFlags::FIXED_VEL_Z,true);
    else                                                          GetGeometry()[0].Set(DEMFlags::FIXED_VEL_Z,false);
    if (GetGeometry()[0].GetDof(ANGULAR_VELOCITY_X).IsFixed()) GetGeometry()[0].Set(DEMFlags::FIXED_ANG_VEL_X,true);
    else                                                          GetGeometry()[0].Set(DEMFlags::FIXED_ANG_VEL_X,false);
    if (GetGeometry()[0].GetDof(ANGULAR_VELOCITY_Y).IsFixed()) GetGeometry()[0].Set(DEMFlags::FIXED_ANG_VEL_Y,true);
    else                                                          GetGeometry()[0].Set(DEMFlags::FIXED_ANG_VEL_Y,false);
    if (GetGeometry()[0].GetDof(ANGULAR_VELOCITY_Z).IsFixed()) GetGeometry()[0].Set(DEMFlags::FIXED_ANG_VEL_Z,true);
    else                                                          GetGeometry()[0].Set(DEMFlags::FIXED_ANG_VEL_Z,false);

    CustomInitialize();

    KRATOS_CATCH( "" )
}

void SphericParticle::FullInitialize(const ProcessInfo& r_process_info)
{
    KRATOS_TRY
    MemberDeclarationFirstStep(r_process_info);
    Initialize();   
    CreateDiscontinuumConstitutiveLaws(r_process_info);    
    KRATOS_CATCH( "" )
}

void SphericParticle::CalculateRightHandSide(VectorType& r_right_hand_side_vector, ProcessInfo& r_process_info,
                                             double dt, const array_1d<double,3>& gravity, int search_control)
{
    KRATOS_TRY           

    array_1d<double, 3> additional_forces;
    array_1d<double, 3> additionally_applied_moment;
    array_1d<double, 3> initial_rotation_moment;
    array_1d<double, 3>& elastic_force       = this->GetGeometry()[0].FastGetSolutionStepValue(ELASTIC_FORCES);
    array_1d<double, 3>& contact_force       = this->GetGeometry()[0].FastGetSolutionStepValue(CONTACT_FORCES);
    array_1d<double, 3>& rigid_element_force = this->GetGeometry()[0].FastGetSolutionStepValue(RIGID_ELEMENT_FORCE);

    mContactForce.clear();
    mContactMoment.clear();
    additional_forces.clear();
    additionally_applied_moment.clear();
    initial_rotation_moment.clear();
    elastic_force.clear();
    contact_force.clear();
    rigid_element_force.clear();

    bool multi_stage_RHS = false;

    ComputeBallToBallContactForce(elastic_force, contact_force, initial_rotation_moment, r_process_info, dt, multi_stage_RHS);
    
    ComputeBallToRigidFaceContactForce(elastic_force, contact_force, initial_rotation_moment, rigid_element_force, r_process_info, dt, search_control);    
    
    if (this->IsNot(DEMFlags::BELONGS_TO_A_CLUSTER)){
        ComputeAdditionalForces(additional_forces, additionally_applied_moment, r_process_info, gravity);
    }

    r_right_hand_side_vector[0] =  mContactForce[0] + additional_forces[0];
    r_right_hand_side_vector[1] =  mContactForce[1] + additional_forces[1];
    r_right_hand_side_vector[2] =  mContactForce[2] + additional_forces[2];
    r_right_hand_side_vector[3] = mContactMoment[0] + additionally_applied_moment[0];
    r_right_hand_side_vector[4] = mContactMoment[1] + additionally_applied_moment[1];
    r_right_hand_side_vector[5] = mContactMoment[2] + additionally_applied_moment[2];

    array_1d<double,3>& total_forces = this->GetGeometry()[0].FastGetSolutionStepValue(TOTAL_FORCES);
    array_1d<double,3>& total_moment = this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT);

    for (int i = 0; i < 3; i++){
        total_forces[i] = r_right_hand_side_vector[i];
        total_moment[i] = r_right_hand_side_vector[3 + i];
    }

    KRATOS_CATCH( "" )
}

void SphericParticle::FirstCalculateRightHandSide(ProcessInfo& r_process_info, double dt, int search_control)
{
    KRATOS_TRY

    /*array_1d<double, 3> initial_rotation_moment;
    array_1d<double, 3>& elastic_force       = this->GetGeometry()[0].FastGetSolutionStepValue(ELASTIC_FORCES);
    array_1d<double, 3>& contact_force       = this->GetGeometry()[0].FastGetSolutionStepValue(CONTACT_FORCES);
    array_1d<double, 3>& rigid_element_force = this->GetGeometry()[0].FastGetSolutionStepValue(RIGID_ELEMENT_FORCE);

    mContactForce.clear();
    mContactMoment.clear();
    initial_rotation_moment.clear();
    elastic_force.clear();
    contact_force.clear();
    rigid_element_force.clear();

    ComputeBallToBallContactForce(elastic_force, contact_force, initial_rotation_moment, r_process_info, dt, true); // The preset argument 'true' should be removed

    std::vector<double>& neighbour_rigid_faces_elastic_contact_force = this->mNeighbourRigidFacesElasticContactForce;
    std::vector<double>& neighbour_rigid_faces_total_contact_force = this->mNeighbourRigidFacesTotalContactForce;
    std::fill(neighbour_rigid_faces_elastic_contact_force.begin(), neighbour_rigid_faces_elastic_contact_force.end(), 0.0);
    std::fill(neighbour_rigid_faces_total_contact_force.begin(), neighbour_rigid_faces_total_contact_force.end(), 0.0);

    if (mFemOldNeighbourIds.size() > 0){
        ComputeBallToRigidFaceContactForce(elastic_force, contact_force, initial_rotation_moment, rigid_element_force, r_process_info, dt, search_control);
    }
*/
    KRATOS_CATCH( "" )
}

void SphericParticle::CollectCalculateRightHandSide(ProcessInfo& r_process_info)
{
    KRATOS_TRY

    for (unsigned int i = 0; i < mNeighbourElements.size(); i++){
        SphericParticle* ineighbour = mNeighbourElements[i];

        if (this->Is(NEW_ENTITY) && ineighbour->Is(NEW_ENTITY)) continue;
        if (this->Id() < ineighbour->Id())                      continue;

        for (unsigned int j = 0; j < ineighbour->mNeighbourElements.size(); j++){  //loop to find the neighbour of the neighbours which is me
            SphericParticle* is_that_me = ineighbour->mNeighbourElements[j];

            if (is_that_me->Id() == this->Id()){
                noalias(mNeighbourElasticContactForces[i]) = - ineighbour->mNeighbourElasticContactForces[j];
                noalias(mNeighbourTotalContactForces[i])  = - ineighbour->mNeighbourTotalContactForces[j];
                noalias(mContactForce) += mNeighbourTotalContactForces[i];
                array_1d<double, 3>& r_elastic_force = this->GetGeometry()[0].FastGetSolutionStepValue(ELASTIC_FORCES);
                noalias(r_elastic_force) += mNeighbourElasticContactForces[i];
                array_1d<double, 3>& r_contact_force = this->GetGeometry()[0].FastGetSolutionStepValue(CONTACT_FORCES);
                noalias(r_contact_force) += mNeighbourTotalContactForces[i];
                break;
            }
        }
    }

    KRATOS_CATCH( "" )
}

void SphericParticle::FinalCalculateRightHandSide(ProcessInfo& r_process_info, double dt, const array_1d<double,3>& gravity)
{
    KRATOS_TRY

    /*if (this->Is(DEMFlags::HAS_ROTATION)){
        const array_1d<double, 3> ang_vel = this->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
        const double coeff_acc            = this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA) / dt;
        array_1d<double, 3> initial_rotation_moment;
        noalias(initial_rotation_moment) = coeff_acc * ang_vel; // the moment needed to stop the spin in one time step

        for (unsigned int i = 0; i < mNeighbourElements.size(); i++){
            SphericParticle* ineighbour = mNeighbourElements[i];

            if (this->Is(NEW_ENTITY) && ineighbour->Is(NEW_ENTITY)) continue;

            array_1d<double, 3> other_to_me_vect;
            noalias(other_to_me_vect) = this->GetGeometry()[0].Coordinates() - ineighbour->GetGeometry()[0].Coordinates();
            double inv_distance                  = 1.0 / DEM_MODULUS_3(other_to_me_vect);
            double other_to_me_vect_unitary[3]   = {other_to_me_vect[0] * inv_distance, other_to_me_vect[1] * inv_distance, other_to_me_vect[2] * inv_distance};
            double projection_to_local_axis2     = DEM_INNER_PRODUCT_3(mNeighbourElasticContactForces[i], other_to_me_vect_unitary);

            if (this->Is(DEMFlags::HAS_ROTATION)){
                ComputeMoments(projection_to_local_axis2, mNeighbourElasticContactForces[i], initial_rotation_moment, other_to_me_vect_unitary, ineighbour, indentation);
            }
        }

    } //if( this->Is(DEMFlags::HAS_ROTATION) )

    array_1d<double, 3> additional_forces;
    array_1d<double, 3> additionally_applied_moment;
    additional_forces.clear();
    additionally_applied_moment.clear();

    ComputeAdditionalForces(additional_forces, additionally_applied_moment, r_process_info, gravity);

    array_1d<double,3>& total_forces = this->GetGeometry()[0].FastGetSolutionStepValue(TOTAL_FORCES);
    array_1d<double,3>& total_moment = this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT);

    noalias(total_forces) = mContactForce  + additional_forces;
    noalias(total_moment) = mContactMoment + additionally_applied_moment;
*/
    KRATOS_CATCH( "" )
}

void SphericParticle::CalculateMaxBallToBallIndentation(double& r_current_max_indentation)
{
    r_current_max_indentation = - std::numeric_limits<double>::max();

    for (unsigned int i = 0; i < mNeighbourElements.size(); i++){
        SphericParticle* ineighbour = mNeighbourElements[i];

        array_1d<double, 3> other_to_me_vect;
        noalias(other_to_me_vect) = this->GetGeometry()[0].Coordinates() - ineighbour->GetGeometry()[0].Coordinates();
        double other_radius                  = ineighbour->GetRadius();
        double distance                      = DEM_MODULUS_3(other_to_me_vect);
        double radius_sum                    = GetRadius() + other_radius;
        double indentation                   = radius_sum - distance;

        r_current_max_indentation = (indentation > r_current_max_indentation) ? indentation : r_current_max_indentation;
    }
}

void SphericParticle::CalculateMaxBallToFaceIndentation(double& r_current_max_indentation)
{
    r_current_max_indentation = - std::numeric_limits<double>::max();

    const std::vector<double>& RF_Pram = this->mNeighbourRigidFacesPram;

    for (unsigned int i = 0; i < mNeighbourRigidFaces.size(); i++){
        
        int ino1                = i * 16;
        double DistPToB         = RF_Pram[ino1 +  9];
        int ContactType         = RF_Pram[ino1 + 15];
        
        if(ContactType > 0){
            double indentation = GetRadius() - DistPToB;
            r_current_max_indentation = (indentation > r_current_max_indentation) ? indentation : r_current_max_indentation;    
            
        }
                
    } //for every rigidface neighbor
}
/*
void SphericParticle::CalculateBalltoBallElasticEnergy(double& total_normal_elastic_energy){
  
  KRATOS_TRY
   
  for (unsigned int i = 0; i < mNeighbourElements.size(); i++) {
        SphericParticle* ineighbour = mNeighbourElements[i];

        if (this->Is(NEW_ENTITY) && ineighbour->Is(NEW_ENTITY)) continue;        
        //if (multi_stage_RHS  &&  this->Id() > ineighbour->Id()) continue;
        
        array_1d<double, 3> other_to_me_vect;
        noalias(other_to_me_vect)    = this->GetGeometry()[0].Coordinates() - ineighbour->GetGeometry()[0].Coordinates();
        const double& other_radius   = ineighbour->GetRadius();
        double distance              = DEM_MODULUS_3(other_to_me_vect);
        double radius_sum            = GetRadius() + other_radius;
        double indentation           = radius_sum - distance;
        double cohesive_force        =  0.0;
        double normal_elastic_energy =  0.0;
        
        if (indentation > 0.0) {
            mDiscontinuumConstitutiveLaw->CalculateElasticEnergy(normal_elastic_energy, indentation, cohesive_force, this, ineighbour);
        }
        total_normal_elastic_energy += normal_elastic_energy;
        
  }//for every neighbor
  
  KRATOS_CATCH("")
}
*/

void SphericParticle::CalculateKinematicEnergy(double& r_kinematic_energy)
{
    const array_1d<double, 3>& vel    = this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
    const array_1d<double, 3> ang_vel = this->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
    const double moment_of_inertia    = this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA);
    double square_of_celerity         = vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2];
    double square_of_angular_celerity = ang_vel[0] * ang_vel[0] + ang_vel[1] * ang_vel[1] + ang_vel[2] * ang_vel[2];

    r_kinematic_energy = 0.5 * (mRealMass * square_of_celerity + moment_of_inertia * square_of_angular_celerity);
}

void SphericParticle::CalculateGravitationalEnergy(const array_1d<double,3>& gravity, double& r_gravitational_energy)
{
    
    double mass = this->GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS);
    double gh = 0.0;
    array_1d<double, 3> coord = this->GetGeometry()[0].Coordinates();

    for (unsigned int i=0;i<3;i++){
      
      gh += coord(i)*-gravity(i); //negative makes gravity introduce positive potential energy
    
    }

    r_gravitational_energy = mass*gh;

}//CalculateGravitationalEnergy


void SphericParticle::CalculateMomentum(array_1d<double, 3>& r_momentum)
{
    const array_1d<double, 3>& vel = this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
    noalias(r_momentum) = mRealMass * vel;
}

void SphericParticle::CalculateLocalAngularMomentum(array_1d<double, 3>& r_angular_momentum)
{
    const array_1d<double, 3> ang_vel  = this->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
    const double moment_of_inertia     = this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA);
    noalias(r_angular_momentum) = moment_of_inertia * ang_vel;
}

void SphericParticle::ComputeNewNeighboursHistoricalData(std::vector<unsigned int>& mTempNeighboursIds, std::vector<array_1d<double, 3> >& mTempNeighbourElasticContactForces,
                                                         std::vector<array_1d<double, 3> >& mTempNeighbourTotalContactForces)
{
    unsigned int new_size = mNeighbourElements.size();
    array_1d<double, 3> vector_of_zeros = ZeroVector(3);
    mTempNeighboursIds.resize(new_size);
    mTempNeighbourElasticContactForces.resize(new_size);
    mTempNeighbourTotalContactForces.resize(new_size);

    for (unsigned int i = 0; i < new_size; i++){
        SphericParticle* i_neighbour = mNeighbourElements[i];
        mTempNeighboursIds[i] = i_neighbour->Id();
        noalias(mTempNeighbourElasticContactForces[i]) = vector_of_zeros;
        noalias(mTempNeighbourTotalContactForces[i]) = vector_of_zeros;

        for (unsigned int j = 0; j < mOldNeighbourIds.size(); j++){
            if (i_neighbour->Id() == mOldNeighbourIds[j]){
                noalias(mTempNeighbourElasticContactForces[i]) = mNeighbourElasticContactForces[j];
                noalias(mTempNeighbourTotalContactForces[i])   = mNeighbourTotalContactForces[j];
                break;
            }
        }
    }

    mOldNeighbourIds.swap(mTempNeighboursIds);
    mNeighbourElasticContactForces.swap(mTempNeighbourElasticContactForces);
    mNeighbourTotalContactForces.swap(mTempNeighbourTotalContactForces);
}

void SphericParticle::ComputeNewRigidFaceNeighboursHistoricalData()
{
    array_1d<double, 3> vector_of_zeros = ZeroVector(3);
    std::vector<DEMWall*>& rNeighbours = this->mNeighbourRigidFaces;
    unsigned int new_size              = rNeighbours.size();
    std::vector<int> temp_neighbours_ids(new_size); //these two temporal vectors are very small, saving them as a member of the particle loses time (usually they consist on 1 member).
    std::vector<array_1d<double, 3> > temp_neighbours_elastic_contact_forces(new_size);
    std::vector<array_1d<double, 3> > temp_neighbours_contact_forces(new_size);

    for (unsigned int i = 0; i<rNeighbours.size(); i++){
        temp_neighbours_ids[i] = static_cast<int>(rNeighbours[i]->Id());
        noalias(temp_neighbours_contact_forces[i]) = vector_of_zeros;
        noalias(temp_neighbours_elastic_contact_forces[i]) = vector_of_zeros;

        for (unsigned int j = 0; j != mFemOldNeighbourIds.size(); j++) {
            if (static_cast<int>(rNeighbours[i]->Id()) == mFemOldNeighbourIds[j]) {
                noalias(temp_neighbours_elastic_contact_forces[i]) = mNeighbourRigidFacesElasticContactForce[j];
                noalias(temp_neighbours_contact_forces[i]) = mNeighbourRigidFacesTotalContactForce[j];
                break;
            }
        }
    }

    mFemOldNeighbourIds.swap(temp_neighbours_ids);
    mNeighbourRigidFacesElasticContactForce.swap(temp_neighbours_elastic_contact_forces);
    mNeighbourRigidFacesTotalContactForce.swap(temp_neighbours_contact_forces);
}

void SphericParticle::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& r_process_info){}

void SphericParticle::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& r_process_info)
{
    rMassMatrix(0,0) = mRealMass;
}

void SphericParticle::EvaluateDeltaDisplacement(double RelDeltDisp[3],
                                                double RelVel[3],
                                                double LocalCoordSystem[3][3],
                                                double OldLocalCoordSystem[3][3],
                                                array_1d<double, 3>& other_to_me_vect,
                                                const array_1d<double, 3>& vel,
                                                const array_1d<double, 3>& delta_displ,
                                                SphericParticle* p_neighbour,
                                                double& distance)
{
    // FORMING LOCAL COORDINATES

    //Notes: Since we will normally inherit the mesh from GiD, we respect the global system X,Y,Z [0],[1],[2]
    //In the local coordinates we will define the normal direction of the contact as the [2] component!!!!!
    //the way the normal direction is defined (other_to_me_vect) compression is positive!!!
    GeometryFunctions::ComputeContactLocalCoordSystem(other_to_me_vect, distance, LocalCoordSystem); //new Local Coordinate System (normalizes other_to_me_vect)

    // FORMING OLD LOCAL COORDINATES
    array_1d<double, 3> old_coord_target;
    noalias(old_coord_target) = this->GetGeometry()[0].Coordinates() - delta_displ;

    const array_1d<double, 3>& other_delta_displ = p_neighbour->GetGeometry()[0].FastGetSolutionStepValue(DELTA_DISPLACEMENT);
    array_1d<double, 3> old_coord_neigh;
    noalias(old_coord_neigh) = p_neighbour->GetGeometry()[0].Coordinates() - other_delta_displ;

    array_1d<double, 3> old_other_to_me_vect;
    noalias(old_other_to_me_vect) = old_coord_target - old_coord_neigh;


    const double old_distance = DEM_MODULUS_3(old_other_to_me_vect);

    GeometryFunctions::ComputeContactLocalCoordSystem(old_other_to_me_vect, old_distance, OldLocalCoordSystem); //Old Local Coordinate System

    // VELOCITIES AND DISPLACEMENTS
    const array_1d<double, 3 >& other_vel = p_neighbour->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);

    RelVel[0] = (vel[0] - other_vel[0]);
    RelVel[1] = (vel[1] - other_vel[1]);
    RelVel[2] = (vel[2] - other_vel[2]);

    // DeltDisp in global coordinates
    RelDeltDisp[0] = (delta_displ[0] - other_delta_displ[0]);
    RelDeltDisp[1] = (delta_displ[1] - other_delta_displ[1]);
    RelDeltDisp[2] = (delta_displ[2] - other_delta_displ[2]);
 
}

void SphericParticle::DisplacementDueToRotation(const double indentation,
                                                double RelDeltDisp[3],
                                                double RelVel[3],
                                                double OldLocalCoordSystem[3][3],
                                                const double& other_radius,
                                                const double& dt,
                                                const array_1d<double, 3>& my_ang_vel,
                                                SphericParticle* p_neighbour)
{
    const array_1d<double, 3>& other_ang_vel = p_neighbour->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
    const array_1d<double, 3>& my_delta_rotation = GetGeometry()[0].FastGetSolutionStepValue(DELTA_ROTATION);
    const array_1d<double, 3>& other_delta_rotation = p_neighbour->GetGeometry()[0].FastGetSolutionStepValue(DELTA_ROTATION);
    array_1d<double, 3> my_arm_vector;
    array_1d<double, 3> other_arm_vector;
    array_1d<double, 3> my_vel_at_contact_point_due_to_rotation;
    array_1d<double, 3> other_vel_at_contact_point_due_to_rotation;
    array_1d<double, 3> my_delta_disp_at_contact_point_due_to_rotation;
    array_1d<double, 3> other_delta_disp_at_contact_point_due_to_rotation;
    
    const double radius_sum        = GetRadius() + other_radius;
    const double inv_radius_sum    = 1.0 / radius_sum;
    const double my_arm_length     = GetRadius()  - indentation * GetRadius()  * inv_radius_sum;
    const double other_arm_length  = other_radius - indentation * other_radius * inv_radius_sum;
        
    my_arm_vector[0] = -OldLocalCoordSystem[2][0] * my_arm_length;
    my_arm_vector[1] = -OldLocalCoordSystem[2][1] * my_arm_length;
    my_arm_vector[2] = -OldLocalCoordSystem[2][2] * my_arm_length;          
    
    GeometryFunctions::CrossProduct(my_ang_vel, my_arm_vector, my_vel_at_contact_point_due_to_rotation);    
    
    other_arm_vector[0] = OldLocalCoordSystem[2][0] * other_arm_length;
    other_arm_vector[1] = OldLocalCoordSystem[2][1] * other_arm_length;
    other_arm_vector[2] = OldLocalCoordSystem[2][2] * other_arm_length;    
    
    GeometryFunctions::CrossProduct(other_ang_vel, other_arm_vector, other_vel_at_contact_point_due_to_rotation);
    
    RelVel[0] += my_vel_at_contact_point_due_to_rotation[0] - other_vel_at_contact_point_due_to_rotation[0];
    RelVel[1] += my_vel_at_contact_point_due_to_rotation[1] - other_vel_at_contact_point_due_to_rotation[1];
    RelVel[2] += my_vel_at_contact_point_due_to_rotation[2] - other_vel_at_contact_point_due_to_rotation[2];                                                                              
    
    GeometryFunctions::CrossProduct(my_delta_rotation,    my_arm_vector,    my_delta_disp_at_contact_point_due_to_rotation   );
    GeometryFunctions::CrossProduct(other_delta_rotation, other_arm_vector, other_delta_disp_at_contact_point_due_to_rotation);    

    // Contribution of the rotation 
    RelDeltDisp[0] += my_delta_disp_at_contact_point_due_to_rotation[0] - other_delta_disp_at_contact_point_due_to_rotation[0];
    RelDeltDisp[1] += my_delta_disp_at_contact_point_due_to_rotation[1] - other_delta_disp_at_contact_point_due_to_rotation[1];
    RelDeltDisp[2] += my_delta_disp_at_contact_point_due_to_rotation[2] - other_delta_disp_at_contact_point_due_to_rotation[2];
}


void SphericParticle::DisplacementDueToRotationMatrix(double DeltDisp[3],
                                                double RelVel[3],
                                                double OldLocalCoordSystem[3][3],
                                                const double& other_radius,
                                                const double& dt,
                                                const array_1d<double, 3>& angular_vel,
                                                SphericParticle* p_neighbour)
{    
        array_1d<double, 3>& neigh_angular_vel = p_neighbour->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY); 
        array_1d<double, 3> temp_angular_vel = angular_vel;
        array_1d<double, 3> angular_velocity = angular_vel;
        array_1d<double, 3> other_angular_velocity = neigh_angular_vel;
        array_1d<double, 3> temp_neigh_angular_vel = neigh_angular_vel;
        const double other_young = p_neighbour->GetYoung();
        const double my_young = GetYoung();
        
        DEM_MULTIPLY_BY_SCALAR_3(temp_angular_vel, dt); //THIS MACRO OPERATION CONVERTS THE temp_angular_vel VARIABLE INTO AN ANGLE, DESPITE THE NAME
        DEM_MULTIPLY_BY_SCALAR_3(temp_neigh_angular_vel, dt); // SAME COMMENT AS IN THE PREVIOUS LINE
        
        double my_rotated_angle = DEM_MODULUS_3(temp_angular_vel);
        double other_rotated_angle = DEM_MODULUS_3(temp_neigh_angular_vel);
        
        array_1d<double, 3> other_to_me_vect;
        noalias(other_to_me_vect)  = this->GetGeometry()[0].Coordinates() - p_neighbour->GetGeometry()[0].Coordinates();
        double distance            = DEM_MODULUS_3(other_to_me_vect);
        double radius_sum          = GetRadius() + other_radius;
        double indentation         = radius_sum - distance;
        
        double arm = GetRadius() - indentation * other_young / (other_young + my_young);
        array_1d<double, 3> e1;
        DEM_COPY_SECOND_TO_FIRST_3(e1, OldLocalCoordSystem[2]);
        DEM_MULTIPLY_BY_SCALAR_3(e1, -arm);
        array_1d<double, 3> new_axes1 = e1;
               
        double other_arm = other_radius - indentation * my_young / (other_young + my_young);
        array_1d<double, 3> e2;
        DEM_COPY_SECOND_TO_FIRST_3(e2, OldLocalCoordSystem[2]);
        DEM_MULTIPLY_BY_SCALAR_3(e2, other_arm);
        array_1d<double, 3> new_axes2 = e2;
        
        if (my_rotated_angle) {
            
            array_1d<double, 3> axis_1;
            axis_1[0] = temp_angular_vel[0] / my_rotated_angle;
            axis_1[1] = temp_angular_vel[1] / my_rotated_angle;
            axis_1[2] = temp_angular_vel[2] / my_rotated_angle;
            
            GeometryFunctions::RotateAVectorAGivenAngleAroundAUnitaryVector(e1, axis_1, my_rotated_angle, new_axes1);
        
        } //if my_rotated_angle

        if (other_rotated_angle) {
            
            array_1d<double, 3> axis_2;
            axis_2[0] = temp_neigh_angular_vel[0] / other_rotated_angle;
            axis_2[1] = temp_neigh_angular_vel[1] / other_rotated_angle;
            axis_2[2] = temp_neigh_angular_vel[2] / other_rotated_angle;
            
            GeometryFunctions::RotateAVectorAGivenAngleAroundAUnitaryVector(e2, axis_2, other_rotated_angle, new_axes2);     
   
        } //if other_rotated_angle
        
        array_1d<double, 3> radial_vector = - other_to_me_vect;
        array_1d<double, 3> other_radial_vector = other_to_me_vect;
        
        GeometryFunctions::normalize(radial_vector);
        GeometryFunctions::normalize(other_radial_vector);
        
        radial_vector *= arm;
        other_radial_vector *= other_arm;

        array_1d<double, 3> vel = ZeroVector(3);
        array_1d<double, 3> other_vel = ZeroVector(3);
        
        GeometryFunctions::CrossProduct(angular_velocity, radial_vector, vel);
        GeometryFunctions::CrossProduct(other_angular_velocity, other_radial_vector, other_vel);
                
        RelVel[0] += vel[0] - other_vel[0];
        RelVel[1] += vel[1] - other_vel[1];
        RelVel[2] += vel[2] - other_vel[2];
        
        // Contribution of the rotation velocity
        DeltDisp[0] += (new_axes1[0] - e1[0] - new_axes2[0] + e2[0]);
        DeltDisp[1] += (new_axes1[1] - e1[1] - new_axes2[1] + e2[1]);
        DeltDisp[2] += (new_axes1[2] - e1[2] - new_axes2[2] + e2[2]);
}

void SphericParticle::ComputeMoments(double NormalLocalElasticContactForce,
                                     array_1d<double, 3>& Force,
                                     array_1d<double, 3>& rInitialRotaMoment,
                                     double LocalCoordSystem2[3],
                                     SphericParticle* p_neighbour,
                                     double indentation,
                                     bool wall)
{
    double arm_length = GetRadius() - indentation;
    
    if (!wall) {
        const double other_young = p_neighbour->GetYoung();
        arm_length = GetRadius() - indentation * other_young / (other_young + GetYoung());
    }
        
    array_1d<double, 3> arm_vector;
    arm_vector[0] = -LocalCoordSystem2[0] * arm_length;
    arm_vector[1] = -LocalCoordSystem2[1] * arm_length;
    arm_vector[2] = -LocalCoordSystem2[2] * arm_length;
    
    array_1d<double, 3> moment_of_this_neighbour;
    GeometryFunctions::CrossProduct(arm_vector, Force, moment_of_this_neighbour);
    noalias(mContactMoment) += moment_of_this_neighbour;

    // ROLLING FRICTION
    if (this->Is(DEMFlags::HAS_ROLLING_FRICTION)) {
        
        double rolling_friction_coeff       = GetRollingFriction() * GetRadius();
        const double other_rolling_friction = p_neighbour->GetRollingFriction();
        double other_rolling_friction_coeff = other_rolling_friction * p_neighbour->GetRadius();
        double equiv_rolling_friction_coeff = std::min(rolling_friction_coeff, other_rolling_friction_coeff);

        if (equiv_rolling_friction_coeff != 0.0) {
            double MaxRotaMoment[3]       = {rInitialRotaMoment[0] + mContactMoment[0], rInitialRotaMoment[1] + mContactMoment[1], rInitialRotaMoment[2] + mContactMoment[2]};
            double CoordSystemMoment1[3]  = {0.0};
            double CoordSystemMoment2[3]  = {0.0};
            double MR[3]                  = {0.0};

            GeometryFunctions::CrossProduct(LocalCoordSystem2, MaxRotaMoment, CoordSystemMoment1);
            if (DEM_MODULUS_3(CoordSystemMoment1) > 0.0) {
                double det_coor_sys_moment_i_1 = 1.0 / DEM_MODULUS_3(CoordSystemMoment1);
                DEM_MULTIPLY_BY_SCALAR_3(CoordSystemMoment1, det_coor_sys_moment_i_1)
            }
                  
            GeometryFunctions::CrossProduct(MaxRotaMoment, CoordSystemMoment1, CoordSystemMoment2);
            if (DEM_MODULUS_3(CoordSystemMoment2) > 0.0) {
                double det_coor_sys_moment_i_2 = 1.0 / DEM_MODULUS_3(CoordSystemMoment2);
                DEM_MULTIPLY_BY_SCALAR_3(CoordSystemMoment2, det_coor_sys_moment_i_2)
            }
                             
            if (DEM_MODULUS_3(CoordSystemMoment1) > 0.0 && DEM_MODULUS_3(CoordSystemMoment2) > 0.0) GeometryFunctions::CrossProduct(CoordSystemMoment2, CoordSystemMoment1, MR);
            double AbsoluteNormalLocalElasticContactForce = fabs(NormalLocalElasticContactForce);
            DEM_MULTIPLY_BY_SCALAR_3(MR, AbsoluteNormalLocalElasticContactForce)
                  
            double MR_now = DEM_MODULUS_3(MR) * equiv_rolling_friction_coeff;
            double MR_max = DEM_MODULUS_3(MaxRotaMoment);
            
            if (MR_max > MR_now) {
                mContactMoment[0] += MR[0] * equiv_rolling_friction_coeff;
                mContactMoment[1] += MR[1] * equiv_rolling_friction_coeff;
                mContactMoment[2] += MR[2] * equiv_rolling_friction_coeff;
            }
            else {
                mContactMoment = - rInitialRotaMoment;
            }
        } // if (equiv_rolling_friction_coeff != 0.0)
    } // if (this->Is(DEMFlags::HAS_ROLLING_FRICTION) )
            
}

void SphericParticle::ComputeBallToBallContactForce(array_1d<double, 3>& r_elastic_force,
                                                    array_1d<double, 3>& r_contact_force,
                                                    array_1d<double, 3>& rInitialRotaMoment,
                                                    ProcessInfo& r_process_info,
                                                    double dt,
                                                    const bool multi_stage_RHS)
{
    KRATOS_TRY

    const array_1d<double, 3>& velocity     = this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
    const array_1d<double, 3>& delta_displ  = this->GetGeometry()[0].FastGetSolutionStepValue(DELTA_DISPLACEMENT);
    const array_1d<double, 3>& ang_velocity = this->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);      
    double LocalCoordSystem[3][3]    = {{0.0}, {0.0}, {0.0}};
    double OldLocalCoordSystem[3][3] = {{0.0}, {0.0}, {0.0}};
    double DeltDisp[3]               = {0.0};
    double LocalDeltDisp[3]          = {0.0};
    double RelVel[3]                 = {0.0};    
    bool sliding = false;
    
    //LOOP OVER NEIGHBORS:
    for (unsigned int i = 0; i < mNeighbourElements.size(); i++) {
        SphericParticle* ineighbour = mNeighbourElements[i];
        
        if (this->Is(NEW_ENTITY) && ineighbour->Is(NEW_ENTITY)) continue;        
        if (multi_stage_RHS  &&  this->Id() > ineighbour->Id()) continue;
        
        array_1d<double, 3> other_to_me_vect;
        noalias(other_to_me_vect)  = this->GetGeometry()[0].Coordinates() - ineighbour->GetGeometry()[0].Coordinates();
        const double& other_radius = ineighbour->GetRadius();
        double distance            = DEM_MODULUS_3(other_to_me_vect);
        double radius_sum          = GetRadius() + other_radius;
        double indentation         = radius_sum - distance;

        DEM_SET_COMPONENTS_TO_ZERO_3(DeltDisp)
        DEM_SET_COMPONENTS_TO_ZERO_3(LocalDeltDisp)
        DEM_SET_COMPONENTS_TO_ZERO_3(RelVel)        
        DEM_SET_COMPONENTS_TO_ZERO_3x3(LocalCoordSystem)
        DEM_SET_COMPONENTS_TO_ZERO_3x3(OldLocalCoordSystem)

        EvaluateDeltaDisplacement(DeltDisp, RelVel, LocalCoordSystem, OldLocalCoordSystem, other_to_me_vect, velocity, delta_displ, ineighbour, distance);

        if (this->Is(DEMFlags::HAS_ROTATION)) {
            DisplacementDueToRotation(indentation, DeltDisp, RelVel, OldLocalCoordSystem, other_radius, dt, ang_velocity, ineighbour);
        }

        double LocalContactForce[3]             = {0.0};
        double GlobalContactForce[3]            = {0.0};
        double LocalElasticContactForce[3]      = {0.0};
        double GlobalElasticContactForce[3]     = {0.0};
        double ViscoDampingLocalContactForce[3] = {0.0};
        double cohesive_force                   =  0.0;

        GeometryFunctions::VectorGlobal2Local(OldLocalCoordSystem, DeltDisp, LocalDeltDisp);   
        
        double OldLocalElasticContactForce[3] = {0.0};

        GeometryFunctions::VectorGlobal2Local(OldLocalCoordSystem, mNeighbourElasticContactForces[i], OldLocalElasticContactForce); // Here we recover the old local forces projected in the new coordinates in the way they were in the old ones; Now they will be increased if its the necessary

        const double previous_indentation = indentation + LocalDeltDisp[2];
        
        
        if (indentation > 0.0) {
            double LocalRelVel[3]            = {0.0};
            GeometryFunctions::VectorGlobal2Local(OldLocalCoordSystem, RelVel, LocalRelVel);            
            mDiscontinuumConstitutiveLaw->CalculateForces(r_process_info, OldLocalElasticContactForce, LocalElasticContactForce, LocalDeltDisp, LocalRelVel, indentation, previous_indentation, ViscoDampingLocalContactForce, cohesive_force, this, ineighbour, sliding);                      
        }
        
        // Transforming to global forces and adding up
        AddUpForcesAndProject(OldLocalCoordSystem, LocalCoordSystem, LocalContactForce, LocalElasticContactForce, GlobalContactForce, GlobalElasticContactForce, ViscoDampingLocalContactForce, cohesive_force, r_elastic_force, r_contact_force, i);
                
        // ROTATION FORCES
        if (this->Is(DEMFlags::HAS_ROTATION) && !multi_stage_RHS) {
            if (this->Is(DEMFlags::HAS_ROLLING_FRICTION) && !multi_stage_RHS) {
                const double coeff_acc      = this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA) / dt;
                noalias(rInitialRotaMoment) = coeff_acc * ang_velocity; // the moment needed to stop the spin in one time step
            }
            
            ComputeMoments(LocalElasticContactForce[2], mNeighbourTotalContactForces[i], rInitialRotaMoment, LocalCoordSystem[2], ineighbour, indentation);
        }
        
        if (r_process_info[STRESS_STRAIN_OPTION]) {
                AddNeighbourContributionToStressTensor(GlobalElasticContactForce, LocalCoordSystem[2], distance, radius_sum);
        }    
        
    }// for each neighbor

    KRATOS_CATCH("")
}// ComputeBallToBallContactForce


void SphericParticle::ComputeRigidFaceToMeVelocity(DEMWall* rObj_2, std::size_t ino,
                                                   double LocalCoordSystem[3][3], 
                                                   double& DistPToB,
                                                   double Weight[4],
                                                   array_1d<double,3>& wall_delta_disp_at_contact_point,
                                                   array_1d<double,3>& wall_velocity_at_contact_point, 
                                                   int& ContactType)
{
    KRATOS_TRY


    const std::vector<double>& RF_Pram = this->mNeighbourRigidFacesPram;

    int ino1 = ino * 16;

    LocalCoordSystem[0][0] = RF_Pram[ino1 +  0];
    LocalCoordSystem[0][1] = RF_Pram[ino1 +  1];
    LocalCoordSystem[0][2] = RF_Pram[ino1 +  2];
    LocalCoordSystem[1][0] = RF_Pram[ino1 +  3];
    LocalCoordSystem[1][1] = RF_Pram[ino1 +  4];
    LocalCoordSystem[1][2] = RF_Pram[ino1 +  5];
    LocalCoordSystem[2][0] = RF_Pram[ino1 +  6];
    LocalCoordSystem[2][1] = RF_Pram[ino1 +  7];
    LocalCoordSystem[2][2] = RF_Pram[ino1 +  8];
    DistPToB               = RF_Pram[ino1 +  9];
    Weight[0]              = RF_Pram[ino1 + 10];
    Weight[1]              = RF_Pram[ino1 + 11];
    Weight[2]              = RF_Pram[ino1 + 12];
    Weight[3]              = RF_Pram[ino1 + 13];
    ContactType            = RF_Pram[ino1 + 15];
    
    for (std::size_t inode = 0; inode < rObj_2->GetGeometry().size(); inode++){
        noalias(wall_velocity_at_contact_point)   += rObj_2->GetGeometry()[inode].FastGetSolutionStepValue(VELOCITY) * Weight[inode];
        noalias(wall_delta_disp_at_contact_point) += rObj_2->GetGeometry()[inode].FastGetSolutionStepValue(DELTA_DISPLACEMENT) * Weight[inode];
    }

    KRATOS_CATCH("")
}

void SphericParticle::UpdateRF_Pram(DEMWall* rObj_2, 
                                    const std::size_t ino,
                                    const double LocalCoordSystem[3][3], 
                                    const double DistPToB,
                                    const double Weight[4], 
                                    const int ContactType)
{
    KRATOS_TRY


    std::vector<double>& RF_Pram = this->mNeighbourRigidFacesPram;

    int ino1 = ino * 16;

    RF_Pram[ino1 +  0] = LocalCoordSystem[0][0];
    RF_Pram[ino1 +  1] = LocalCoordSystem[0][1];
    RF_Pram[ino1 +  2] = LocalCoordSystem[0][2];
    RF_Pram[ino1 +  3] = LocalCoordSystem[1][0];
    RF_Pram[ino1 +  4] = LocalCoordSystem[1][1];
    RF_Pram[ino1 +  5] = LocalCoordSystem[1][2];
    RF_Pram[ino1 +  6] = LocalCoordSystem[2][0];
    RF_Pram[ino1 +  7] = LocalCoordSystem[2][1];
    RF_Pram[ino1 +  8] = LocalCoordSystem[2][2];
    RF_Pram[ino1 +  9] = DistPToB;
    RF_Pram[ino1 + 10] = Weight[0];
    RF_Pram[ino1 + 11] = Weight[1];
    RF_Pram[ino1 + 12] = Weight[2];
    RF_Pram[ino1 + 13] = Weight[3];
    RF_Pram[ino1 + 15] = ContactType;


    KRATOS_CATCH("")
}

void SphericParticle::ComputeBallToRigidFaceContactForce(array_1d<double, 3>& r_elastic_force,
                                                         array_1d<double, 3>& r_contact_force,
                                                         array_1d<double, 3>& rInitialRotaMoment,
                                                         array_1d<double, 3>& rigid_element_force,
                                                         ProcessInfo& r_process_info,
                                                         double mTimeStep,
                                                         int search_control)
{
    KRATOS_TRY

    std::vector<DEMWall*>& rNeighbours   = this->mNeighbourRigidFaces;        
    array_1d<double, 3> vel              = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);            
    const array_1d<double,3>& AngularVel = GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
    
    for (unsigned int i = 0; i < rNeighbours.size(); i++) {
        
        DEMWall* wall = rNeighbours[i];        
        double LocalElasticContactForce[3]       = {0.0};
        double GlobalElasticContactForce[3]      = {0.0};        
        double ViscoDampingLocalContactForce[3]  = {0.0};
        double cohesive_force                    =  0.0;
        double LocalCoordSystem[3][3]            = {{0.0}, {0.0}, {0.0}};
        array_1d<double, 3> wall_velocity_at_contact_point = ZeroVector(3);
        array_1d<double, 3> wall_delta_disp_at_contact_point = ZeroVector(3);
        bool sliding = false;
        

        double ini_delta = GetInitialDeltaWithFEM(i);              
        double DistPToB = 0.0;

        int ContactType = -1;
        double Weight[4] = {0.0};                
        
        ComputeRigidFaceToMeVelocity(rNeighbours[i], i, LocalCoordSystem, DistPToB, Weight, wall_delta_disp_at_contact_point, wall_velocity_at_contact_point, ContactType);

        //Renew neighbour parameters for steps between searches
        if (ContactType == 1 || ContactType == 2 || ContactType == 3) {
            if (search_control == 1) { //Search active but not performed in this timestep
                UpdateDistanceToWall(rNeighbours[i], i, LocalCoordSystem, DistPToB, Weight, ContactType);                                
            }
        }
        if (ContactType == 1 || ContactType == 2 || ContactType == 3) {
            double indentation = -(DistPToB - GetRadius()) - ini_delta;
            double DeltDisp[3] = {0.0};
            double DeltVel [3] = {0.0};

            DeltVel[0] = vel[0] - wall_velocity_at_contact_point[0];
            DeltVel[1] = vel[1] - wall_velocity_at_contact_point[1];
            DeltVel[2] = vel[2] - wall_velocity_at_contact_point[2];

            // For translation movement delta displacement
            const array_1d<double, 3>& delta_displ  = this->GetGeometry()[0].FastGetSolutionStepValue(DELTA_DISPLACEMENT);
            DeltDisp[0] = delta_displ[0] - wall_delta_disp_at_contact_point[0];
            DeltDisp[1] = delta_displ[1] - wall_delta_disp_at_contact_point[1];
            DeltDisp[2] = delta_displ[2] - wall_delta_disp_at_contact_point[2];

            const double actual_distance_to_contact_point = GetRadius() - indentation;
            
            if (this->Is(DEMFlags::HAS_ROTATION)) {            
                const array_1d<double,3>& delta_rotation = GetGeometry()[0].FastGetSolutionStepValue(DELTA_ROTATION);
                                
                array_1d<double, 3> actual_arm_vector;
                actual_arm_vector[0] = -LocalCoordSystem[2][0] * actual_distance_to_contact_point;
                actual_arm_vector[1] = -LocalCoordSystem[2][1] * actual_distance_to_contact_point;
                actual_arm_vector[2] = -LocalCoordSystem[2][2] * actual_distance_to_contact_point;
                
                double tangential_vel[3]           = {0.0};
                double tangential_displacement_due_to_rotation[3]  = {0.0};
                GeometryFunctions::CrossProduct(AngularVel, actual_arm_vector, tangential_vel); 
                GeometryFunctions::CrossProduct(delta_rotation, actual_arm_vector, tangential_displacement_due_to_rotation); 

                DEM_ADD_SECOND_TO_FIRST(DeltVel, tangential_vel)
                DEM_ADD_SECOND_TO_FIRST(DeltDisp, tangential_displacement_due_to_rotation)
            }

            double LocalDeltDisp[3] = {0.0};
            GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, DeltDisp, LocalDeltDisp);
            
            double OldLocalElasticContactForce[3] = {0.0}; 

            GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, mNeighbourRigidFacesElasticContactForce[i], OldLocalElasticContactForce); 
            const double previous_indentation = indentation + LocalDeltDisp[2];
            
            if (indentation > 0.0) {
                double LocalRelVel[3]            = {0.0};
                GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, DeltVel, LocalRelVel);
                
                mDiscontinuumConstitutiveLaw->CalculateForcesWithFEM(r_process_info,OldLocalElasticContactForce, LocalElasticContactForce, LocalDeltDisp, LocalRelVel, indentation,
                                                                     previous_indentation, ViscoDampingLocalContactForce, cohesive_force, this, wall, sliding);         
			
            }

            double LocalContactForce[3]  = {0.0};
            double GlobalContactForce[3] = {0.0};

            AddUpFEMForcesAndProject(LocalCoordSystem, LocalContactForce, LocalElasticContactForce, GlobalContactForce,
                                     GlobalElasticContactForce, ViscoDampingLocalContactForce, cohesive_force, r_elastic_force, r_contact_force, i);
            
            rigid_element_force[0] -= GlobalContactForce[0];
            rigid_element_force[1] -= GlobalContactForce[1];
            rigid_element_force[2] -= GlobalContactForce[2];
            array_1d<double, 3> GlobalContactForce_array;
            DEM_COPY_SECOND_TO_FIRST_3(GlobalContactForce_array,GlobalContactForce)
                    
            if (this->Is(DEMFlags::HAS_ROTATION)) {
                if (this->Is(DEMFlags::HAS_ROLLING_FRICTION)) {
                    const double coeff_acc      = this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA) / mTimeStep;
                    noalias(rInitialRotaMoment) = coeff_acc * AngularVel; // the moment needed to stop the spin in one time step
                }
              
              ComputeMoments(LocalElasticContactForce[2], GlobalContactForce_array, rInitialRotaMoment, LocalCoordSystem[2], this, indentation, true); //WARNING: sending itself as the neighbor!!
            }
        
            //WEAR        
            if (wall->GetProperties()[PRINT_WEAR]) {                
                const double area                        = KRATOS_M_PI * GetRadius() * GetRadius();
                const double density                     = GetDensity();
                const double inverse_of_volume           = 1.0 / (4.0 * 0.333333333333333 * area * GetRadius());
                ComputeWear(LocalCoordSystem, vel, DeltVel, mTimeStep, density, sliding, inverse_of_volume, LocalElasticContactForce[2], wall);
            } //wall->GetProperties()[PRINT_WEAR] if 
            
            if (r_process_info[STRESS_STRAIN_OPTION]) {
                AddWallContributionToStressTensor(GlobalElasticContactForce, LocalCoordSystem[2], actual_distance_to_contact_point, 0.0);
            } 
        } //ContactType if          
    } //rNeighbours.size loop
    
    KRATOS_CATCH("")
}// ComputeBallToRigidFaceContactForce

void SphericParticle::UpdateDistanceToWall(DEMWall* const wall, 
                                            const int neighbour_index, 
                                            double LocalCoordSystem[3][3], 
                                            double& DistPToB, 
                                            double Weight[4], 
                                            int& ContactType){

    size_t FE_size = wall->GetGeometry().size(); 
    std::vector< array_1d <double,3> >Coord;
    Coord.resize(FE_size, array_1d<double,3>(3,0.0) );
    std::vector<double> TempWeight;
    TempWeight.resize(FE_size);
  
    double total_weight = 0.0;
    int points = 0;
    int inode1 = 0, inode2 = 0;

    for (unsigned int inode = 0; inode < FE_size; inode++) {

        if (Weight[inode] > 1.0e-6){
          
            for (unsigned int j = 0; j<3; j++)
            {
                Coord[inode][j] = wall->GetGeometry()[inode].Coordinates()[j];
            }
            total_weight = total_weight + Weight[inode];
            points++;
            if (points == 1) {inode1 = inode;}
            if (points == 2) {inode2 = inode;}
            
        }

        if (fabs(total_weight - 1.0) < 1.0e-6){
            break;
        }
    }

    bool contact_exists = true;
    array_1d<double, 3>& node_coordinates = this->GetGeometry()[0].Coordinates();
    double node_coor[3] = {0.0};
    DEM_COPY_SECOND_TO_FIRST_3(node_coor, node_coordinates)
    
    const double radius = this->GetRadius();

    if (points == 3 || points == 4) 
    {
      unsigned int dummy_current_edge_index;
      contact_exists = GeometryFunctions::FacetCheck(Coord, node_coor, radius, LocalCoordSystem, DistPToB, TempWeight, dummy_current_edge_index);
      Weight[0]=TempWeight[0];
      Weight[1]=TempWeight[1];
      Weight[2]=TempWeight[2];
      if(points == 4)
      {
        Weight[3] = TempWeight[3];
      }
      else
      {
       Weight[3] = 0.0; 
      }
    }
 
    if (points == 2) {
      
        double eta = 0.0;
        contact_exists = GeometryFunctions::EdgeCheck(Coord[inode1], Coord[inode2], node_coor, radius, LocalCoordSystem, DistPToB, eta);
     
        Weight[inode1] = 1-eta;
        Weight[inode2] = eta;
        
    }

    if (points == 1) {
        contact_exists = GeometryFunctions::VertexCheck(Coord[inode1], node_coor, radius, LocalCoordSystem, DistPToB);
        Weight[inode1] = 1.0;
    }

    if (contact_exists == false) {ContactType = -1;}

    UpdateRF_Pram(wall, neighbour_index, LocalCoordSystem, DistPToB, Weight, ContactType);

}

void SphericParticle::ComputeWear(double LocalCoordSystem[3][3], array_1d<double, 3>& vel, double tangential_vel[3],
                                  double mTimeStep, double density, bool sliding, double inverse_of_volume,
                                  double LocalElasticContactForce, DEMWall* wall) {
        
    array_1d<double, 3>& node_coor_array = this->GetGeometry()[0].Coordinates();
    array_1d<double, 3> local_vel;
    GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, vel, local_vel);
    double non_dim_volume_wear;
    double WallSeverityOfWear           = wall->GetProperties()[SEVERITY_OF_WEAR];
    double WallImpactSeverityOfWear     = wall->GetProperties()[IMPACT_WEAR_SEVERITY];
    double InverseOfWallBrinellHardness = 1.0 / (wall->GetProperties()[BRINELL_HARDNESS]);
    double vel_module = DEM_MODULUS_3(vel);
    double quotient_of_vels = fabs(local_vel[2]) / vel_module;
    double Sliding_0 = tangential_vel[0] * mTimeStep;
    double Sliding_1 = tangential_vel[1] * mTimeStep;
    double non_dim_impact_wear = 0.5 * WallImpactSeverityOfWear * InverseOfWallBrinellHardness * density * quotient_of_vels * quotient_of_vels
                                 * pow(1.0 - 4.0 * (1.0 - quotient_of_vels), 2) * vel_module * vel_module;

    if (sliding) {
    non_dim_volume_wear = WallSeverityOfWear * InverseOfWallBrinellHardness * inverse_of_volume * LocalElasticContactForce
                        * sqrt(Sliding_0 * Sliding_0 + Sliding_1 * Sliding_1);
    } else {
        non_dim_volume_wear = 0.0;
    }

    //COMPUTING THE PROJECTED POINT

    array_1d<double, 3> normal_to_wall;

    wall->CalculateNormal(normal_to_wall);

    array_1d<double, 3> relative_vector = wall->GetGeometry()[0].Coordinates() - node_coor_array; //We could have chosen [1] or [2], also.

    double dot_prod = DEM_INNER_PRODUCT_3(relative_vector, normal_to_wall);

    DEM_MULTIPLY_BY_SCALAR_3(normal_to_wall, dot_prod);

    array_1d<double, 3> inner_point = node_coor_array + normal_to_wall;

    array_1d<double, 3> relative_vector_0 = inner_point - wall->GetGeometry()[0].Coordinates();
    array_1d<double, 3> relative_vector_1 = inner_point - wall->GetGeometry()[1].Coordinates();
    array_1d<double, 3> relative_vector_2 = inner_point - wall->GetGeometry()[2].Coordinates();

    double distance_0 = DEM_MODULUS_3(relative_vector_0);
    double distance_1 = DEM_MODULUS_3(relative_vector_1);
    double distance_2 = DEM_MODULUS_3(relative_vector_2);
    double inverse_of_total_distance = 1.0 / (distance_0 + distance_1 + distance_2);

    double weight_0 = 1.0 - (distance_1 + distance_2) * inverse_of_total_distance;
    double weight_1 = 1.0 - (distance_2 + distance_0) * inverse_of_total_distance;
    double weight_2 = 1.0 - (distance_0 + distance_1) * inverse_of_total_distance; // It could be also: weight_2 = 1.0 - weight_0 - weight_1;
    
    
    
    wall->GetGeometry()[0].SetLock();
    wall->GetGeometry()[0].FastGetSolutionStepValue(NON_DIMENSIONAL_VOLUME_WEAR) += weight_0 * non_dim_volume_wear;
    wall->GetGeometry()[0].FastGetSolutionStepValue(IMPACT_WEAR) += weight_0 * non_dim_impact_wear;
    wall->GetGeometry()[0].UnSetLock();

    wall->GetGeometry()[1].SetLock();
    wall->GetGeometry()[1].FastGetSolutionStepValue(NON_DIMENSIONAL_VOLUME_WEAR) += weight_1 * non_dim_volume_wear;
    wall->GetGeometry()[1].FastGetSolutionStepValue(IMPACT_WEAR) += weight_1 * non_dim_impact_wear;
    wall->GetGeometry()[1].UnSetLock();

    wall->GetGeometry()[2].SetLock();
    wall->GetGeometry()[2].FastGetSolutionStepValue(NON_DIMENSIONAL_VOLUME_WEAR) += weight_2 * non_dim_volume_wear;
    wall->GetGeometry()[2].FastGetSolutionStepValue(IMPACT_WEAR) += weight_2 * non_dim_impact_wear;
    wall->GetGeometry()[2].UnSetLock();
            
} //ComputeWear

void SphericParticle::CreateDiscontinuumConstitutiveLaws(const ProcessInfo& r_process_info)
{

    mDiscontinuumConstitutiveLaw = GetProperties()[DEM_DISCONTINUUM_CONSTITUTIVE_LAW_POINTER]->Clone();
    mDiscontinuumConstitutiveLaw->Initialize(r_process_info);

}

void SphericParticle::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& r_process_info){}

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

void SphericParticle::InitializeSolutionStep(ProcessInfo& r_process_info)
{
    KRATOS_TRY
            
    mRadius = this->GetGeometry()[0].FastGetSolutionStepValue(RADIUS); //Just in case that someone is overwriting the Radius from Python
    
    if (r_process_info[STRESS_STRAIN_OPTION]) {
        this->GetGeometry()[0].FastGetSolutionStepValue(REPRESENTATIVE_VOLUME) = 0.0;
        for (int i=0; i<3; i++) {  for (int j=0; j<3; j++) { (*mStressTensor)(i,j) = 0.0; } }          
    }

    KRATOS_CATCH("")
}

void SphericParticle::AddNeighbourContributionToStressTensor(const double Force[3],
                                                            const double other_to_me_vect[3],
                                                            const double distance,
                                                            const double radius_sum) {

        KRATOS_TRY

        double gap = distance - radius_sum;
        double real_distance = GetRadius() + 0.5 * gap;
//        double& rRepresentative_Volume = this->GetGeometry()[0].FastGetSolutionStepValue(REPRESENTATIVE_VOLUME);
//        rRepresentative_Volume += 0.33333333333333 * (real_distance * contact_area);

        array_1d<double, 3> normal_vector_on_contact;
        normal_vector_on_contact[0] = -1 * other_to_me_vect[0]; //outwards
        normal_vector_on_contact[1] = -1 * other_to_me_vect[1]; //outwards
        normal_vector_on_contact[2] = -1 * other_to_me_vect[2]; //outwards
        
        array_1d<double, 3> x_centroid = real_distance * normal_vector_on_contact;

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                (*mStressTensor)(i,j) += (x_centroid[j]) * Force[i]; //ref: Katalin Bagi 1995 Mean stress tensor                                     
            }
        }        
        KRATOS_CATCH("")
}

void SphericParticle::AddContributionToRepresentativeVolume(const double distance,
                                            const double radius_sum,
                                            const double contact_area) {

        KRATOS_TRY

        double gap = distance - radius_sum;
        double real_distance = GetRadius() + 0.5 * gap;
        double& rRepresentative_Volume = this->GetGeometry()[0].FastGetSolutionStepValue(REPRESENTATIVE_VOLUME);
        rRepresentative_Volume += 0.33333333333333 * (real_distance * contact_area);

        KRATOS_CATCH("")
}

void SphericParticle::AddWallContributionToStressTensor(const double Force[3],
                                                        const double other_to_me_vect[3],
                                                        const double distance,
                                                        const double contact_area) {

        KRATOS_TRY

        double& rRepresentative_Volume = this->GetGeometry()[0].FastGetSolutionStepValue(REPRESENTATIVE_VOLUME);
        rRepresentative_Volume += 0.33333333333333 * (distance * contact_area);

        array_1d<double, 3> normal_vector_on_contact;
        normal_vector_on_contact[0] = -1 * other_to_me_vect[0]; //outwards
        normal_vector_on_contact[1] = -1 * other_to_me_vect[1]; //outwards
        normal_vector_on_contact[2] = -1 * other_to_me_vect[2]; //outwards
        
        array_1d<double, 3> x_centroid = distance * normal_vector_on_contact;

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                (*mStressTensor)(i,j) += (x_centroid[j]) * Force[i]; //ref: Katalin Bagi 1995 Mean stress tensor                                     
            }
        }
        
        KRATOS_CATCH("")
}

void SphericParticle::FinalizeSolutionStep(ProcessInfo& r_process_info){
    
    KRATOS_TRY

    double& rRepresentative_Volume = this->GetGeometry()[0].FastGetSolutionStepValue(REPRESENTATIVE_VOLUME);
    const double sphere_volume = 4 *KRATOS_M_PI_3 * GetRadius() * GetRadius() * GetRadius();
    if ((rRepresentative_Volume <= sphere_volume)) { //In case it gets 0.0 (discontinuum). Also sometimes the error can be too big. This puts some bound to the error for continuum.
        rRepresentative_Volume = sphere_volume;
    }

    //Update sphere mass and inertia taking into acount the real volume of the represented volume:
    mRealMass = rRepresentative_Volume * GetDensity();
    GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) = mRealMass;
    GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA) = 0.4 * mRealMass * GetRadius() * GetRadius();

    if (r_process_info[STRESS_STRAIN_OPTION]) {
        //Divide Stress Tensor by the total volume:
        for (int i = 0; i < 3; i++) {  for (int j = 0; j < 3; j++) {  (*mStressTensor)(i,j) /= rRepresentative_Volume;   }  }

        //The following operation symmetrizes the tensor. We will work with the symmetric stress tensor always, because the non-symmetric one is being filled while forces are being calculated
        for (int i = 0; i < 3; i++) {
            for (int j = i; j < 3; j++) {
                (*mSymmStressTensor)(i,j) = (*mSymmStressTensor)(j,i) = 0.5 * ((*mStressTensor)(i,j) + (*mStressTensor)(j,i));
            }
        }                           
    }      
    KRATOS_CATCH("")
}

void SphericParticle::PrepareForPrinting(ProcessInfo& r_process_info){
    
    if (r_process_info[STRESS_STRAIN_OPTION]) {
        this->GetGeometry()[0].FastGetSolutionStepValue(DEM_STRESS_XX) = (*mSymmStressTensor)(0,0);
        this->GetGeometry()[0].FastGetSolutionStepValue(DEM_STRESS_XY) = (*mSymmStressTensor)(0,1);
        this->GetGeometry()[0].FastGetSolutionStepValue(DEM_STRESS_XZ) = (*mSymmStressTensor)(0,2);
        this->GetGeometry()[0].FastGetSolutionStepValue(DEM_STRESS_YX) = (*mSymmStressTensor)(1,0);
        this->GetGeometry()[0].FastGetSolutionStepValue(DEM_STRESS_YY) = (*mSymmStressTensor)(1,1);
        this->GetGeometry()[0].FastGetSolutionStepValue(DEM_STRESS_YZ) = (*mSymmStressTensor)(1,2);
        this->GetGeometry()[0].FastGetSolutionStepValue(DEM_STRESS_ZX) = (*mSymmStressTensor)(2,0);
        this->GetGeometry()[0].FastGetSolutionStepValue(DEM_STRESS_ZY) = (*mSymmStressTensor)(2,1);
        this->GetGeometry()[0].FastGetSolutionStepValue(DEM_STRESS_ZZ) = (*mSymmStressTensor)(2,2);
    }
}

void SphericParticle::CustomInitialize(){}

void SphericParticle::ComputeAdditionalForces(array_1d<double, 3>& externally_applied_force,
                                              array_1d<double, 3>& externally_applied_moment,
                                              ProcessInfo& r_process_info, const array_1d<double,3>& gravity)
{
    noalias(externally_applied_force) += mRealMass * gravity;
}

void SphericParticle::AddUpForcesAndProject(double OldCoordSystem[3][3],
                                            double LocalCoordSystem[3][3],
                                            double LocalContactForce[3],
                                            double LocalElasticContactForce[3],
                                            double GlobalContactForce[3],
                                            double GlobalElasticContactForce[3],
                                            double ViscoDampingLocalContactForce[3],
                                            const double cohesive_force,
                                            array_1d<double, 3>& r_elastic_force,
                                            array_1d<double, 3>& r_contact_force,
                                            const unsigned int i_neighbour_count)
{
    for (unsigned int index = 0; index < 3; index++) {
        LocalContactForce[index] = LocalElasticContactForce[index] + ViscoDampingLocalContactForce[index];
    }
    LocalContactForce[2] -= cohesive_force;

    GeometryFunctions::VectorLocal2Global(OldCoordSystem, LocalElasticContactForce, GlobalElasticContactForce);
    GeometryFunctions::VectorLocal2Global(OldCoordSystem, LocalContactForce, GlobalContactForce);

    // Saving contact forces (We need to, since tangential elastic force is history-dependent)
    DEM_COPY_SECOND_TO_FIRST_3(mNeighbourElasticContactForces[i_neighbour_count], GlobalElasticContactForce)
    DEM_COPY_SECOND_TO_FIRST_3(mNeighbourTotalContactForces[i_neighbour_count], GlobalContactForce)
    DEM_ADD_SECOND_TO_FIRST(mContactForce, GlobalContactForce)
    DEM_ADD_SECOND_TO_FIRST(r_elastic_force, GlobalElasticContactForce)
    DEM_ADD_SECOND_TO_FIRST(r_contact_force, GlobalContactForce)
}

void SphericParticle::AddUpFEMForcesAndProject(double LocalCoordSystem[3][3],
                                               double LocalContactForce[3],
                                               double LocalElasticContactForce[3],
                                               double GlobalContactForce[3],
                                               double GlobalElasticContactForce[3],
                                               double ViscoDampingLocalContactForce[3],
                                               const double cohesive_force,
                                               array_1d<double, 3>& r_elastic_force,
                                               array_1d<double, 3>& r_contact_force,
                                               const unsigned int iRigidFaceNeighbour)
{
    for (unsigned int index = 0; index < 3; index++) {
        LocalContactForce[index] = LocalElasticContactForce[index] + ViscoDampingLocalContactForce[index];
    }
    LocalContactForce[2] -= cohesive_force;

    GeometryFunctions::VectorLocal2Global(LocalCoordSystem, LocalElasticContactForce, GlobalElasticContactForce);
    GeometryFunctions::VectorLocal2Global(LocalCoordSystem, LocalContactForce, GlobalContactForce);
    // Saving contact forces (We need to, since tangential elastic force is history-dependent)
    DEM_COPY_SECOND_TO_FIRST_3(mNeighbourRigidFacesElasticContactForce[iRigidFaceNeighbour],GlobalElasticContactForce)
    DEM_COPY_SECOND_TO_FIRST_3(mNeighbourRigidFacesTotalContactForce[iRigidFaceNeighbour],GlobalContactForce)
    DEM_ADD_SECOND_TO_FIRST(mContactForce, GlobalContactForce)
    DEM_ADD_SECOND_TO_FIRST(r_elastic_force, GlobalElasticContactForce)
    DEM_ADD_SECOND_TO_FIRST(r_contact_force, GlobalContactForce)
    
}

void SphericParticle::MemberDeclarationFirstStep(const ProcessInfo& r_process_info)

{
    // Passing the element id to the node upon initialization    
    if (r_process_info[PRINT_EXPORT_ID] == 1) {
        this->GetGeometry()[0].FastGetSolutionStepValue(EXPORT_ID) = double(this->Id());
    }

    if (r_process_info[ROTATION_OPTION])         this->Set(DEMFlags::HAS_ROTATION, true);
    else                                         this->Set(DEMFlags::HAS_ROTATION, false);
    if (r_process_info[ROLLING_FRICTION_OPTION]) this->Set(DEMFlags::HAS_ROLLING_FRICTION, true);
    else                                         this->Set(DEMFlags::HAS_ROLLING_FRICTION, false);
    if (r_process_info[CRITICAL_TIME_OPTION])    this->Set(DEMFlags::HAS_CRITICAL_TIME, true);
    else                                         this->Set(DEMFlags::HAS_CRITICAL_TIME, false);
    if (r_process_info[ROTATION_SPRING_OPTION])  this->Set(DEMFlags::HAS_ROTATION_SPRING, true);
    else                                         this->Set(DEMFlags::HAS_ROTATION_SPRING, false);
                                                 
    if (r_process_info[STRESS_STRAIN_OPTION]){        
        mStressTensor = new Matrix(3,3);
        *mStressTensor = ZeroMatrix(3,3);
        
        mSymmStressTensor = new Matrix(3,3);
        *mSymmStressTensor = ZeroMatrix(3,3);
    }
    else {
        mStressTensor = NULL;
        mSymmStressTensor = NULL;
    }

    AdditionalMemberDeclarationFirstStep(r_process_info);
    
}

double SphericParticle::GetInitialDeltaWithFEM(int index) //only available in continuum_particle
{
    return 0.0;
}

void SphericParticle::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& r_process_info)
{
    KRATOS_TRY
    
    //CRITICAL DELTA CALCULATION

    if (rVariable == DELTA_TIME){
        double mass = mRealMass;
        double coeff = r_process_info[NODAL_MASS_COEFF];

        if (coeff > 1.0){
            KRATOS_THROW_ERROR(std::runtime_error, "The coefficient assigned for virtual mass is larger than one. Virtual_mass_coeff is ", coeff);
        }

        else if ((coeff == 1.0) && (r_process_info[VIRTUAL_MASS_OPTION])){
            Output = 9.0E09;
        }

        else {

            if (r_process_info[VIRTUAL_MASS_OPTION]){
                mass = mass / (1 - coeff);
            }
            
            double eq_mass = 0.5*mass; //"mass" of the contact

            double kn = 0.0;
            double kt = 0.0;
            double ini_delta = 0.05*GetRadius(); // Hertz needs an initial Delta, linear ignores it
            
            mDiscontinuumConstitutiveLaw->GetContactStiffness(this, this, ini_delta, kn, kt);
             
            //double K = KRATOS_M_PI * GetYoung() * GetRadius(); //M. Error, should be the same that the local definition.

            Output = 0.34 * sqrt(eq_mass / kn);

            if (this->Is(DEMFlags::HAS_ROTATION) ){
                //Output *= 0.5; //factor for critical time step when rotation is allowed.
            }
        }

        return;
    }

    AdditionalCalculate(rVariable, Output, r_process_info);

    KRATOS_CATCH("")
}// Calculate

void SphericParticle::Calculate(const Variable<array_1d<double, 3> >& rVariable, array_1d<double, 3>& Output,
                                const ProcessInfo& r_process_info)
{
    if (rVariable == MOMENTUM){
        CalculateMomentum(Output);
    }

    else if (rVariable == ANGULAR_MOMENTUM){
        CalculateLocalAngularMomentum(Output);
    }
}

void SphericParticle::Calculate(const Variable<Vector >& rVariable, Vector& Output, const ProcessInfo& r_process_info){}
void SphericParticle::Calculate(const Variable<Matrix >& rVariable, Matrix& Output, const ProcessInfo& r_process_info){}

void SphericParticle::AdditionalMemberDeclarationFirstStep(const ProcessInfo& r_process_info){}
void SphericParticle::AdditionalCalculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& r_process_info){}

int    SphericParticle::GetClusterId()                                                   { return mClusterId;      }
void   SphericParticle::SetClusterId(int givenId)                                        { mClusterId = givenId;   }
double SphericParticle::GetRadius()                                                      { return mRadius;         }
void   SphericParticle::SetRadius(double radius)                                         { mRadius = radius;       }
double SphericParticle::GetSearchRadius()                                                { return GetRadius();     }
double SphericParticle::GetMass()                                                        { return mRealMass;       }
void   SphericParticle::SetMass(double real_mass)                                        { mRealMass = real_mass;  }

double SphericParticle::GetYoung()                                                       { return GetFastProperties()->GetYoung();                     }
double SphericParticle::GetRollingFriction()                                             { return GetFastProperties()->GetRollingFriction();           }
double SphericParticle::GetPoisson()                                                     { return GetFastProperties()->GetPoisson();                   }
double SphericParticle::GetTgOfFrictionAngle()                                           { return GetFastProperties()->GetTgOfFrictionAngle() ;        }
double SphericParticle::GetCoefficientOfRestitution()                                    { return GetFastProperties()->GetCoefficientOfRestitution();  }
double SphericParticle::GetLnOfRestitCoeff()                                             { return GetFastProperties()->GetLnOfRestitCoeff();           }
double SphericParticle::GetDensity()                                                     { return GetFastProperties()->GetDensity();                   }
int    SphericParticle::GetParticleMaterial()                                            { return GetFastProperties()->GetParticleMaterial();          }
double SphericParticle::GetParticleCohesion()                                            { return GetFastProperties()->GetParticleCohesion();          }

void   SphericParticle::SetYoungFromProperties(double* young)                            { GetFastProperties()->SetYoungFromProperties( young);                             }
void   SphericParticle::SetRollingFrictionFromProperties(double* rolling_friction)       { GetFastProperties()->SetRollingFrictionFromProperties( rolling_friction);        }
void   SphericParticle::SetPoissonFromProperties(double* poisson)                        { GetFastProperties()->SetPoissonFromProperties( poisson);                         }
void   SphericParticle::SetTgOfFrictionAngleFromProperties(double* tg_of_friction_angle) { GetFastProperties()->SetTgOfFrictionAngleFromProperties( tg_of_friction_angle);  }
void   SphericParticle::SetCoefficientOfRestitutionFromProperties(double* coefficient_of_restitution)     { GetFastProperties()->SetCoefficientOfRestitutionFromProperties( coefficient_of_restitution);      }
void   SphericParticle::SetLnOfRestitCoeffFromProperties(double* ln_of_restit_coeff)     { GetFastProperties()->SetLnOfRestitCoeffFromProperties( ln_of_restit_coeff);      }
void   SphericParticle::SetDensityFromProperties(double* density)                        { GetFastProperties()->SetDensityFromProperties( density);                         }
void   SphericParticle::SetParticleMaterialFromProperties(int* particle_material)        { GetFastProperties()->SetParticleMaterialFromProperties( particle_material);      }
void   SphericParticle::SetParticleCohesionFromProperties(double* particle_cohesion)     { GetFastProperties()->SetParticleCohesionFromProperties( particle_cohesion);      }

PropertiesProxy* SphericParticle::GetFastProperties()                                    { return mFastProperties;                                                          }
void   SphericParticle::SetFastProperties(PropertiesProxy* pProps)                       { mFastProperties = pProps;                                                        }
void   SphericParticle::SetFastProperties(std::vector<PropertiesProxy>& list_of_proxies) {
    for (unsigned int j=0;j<list_of_proxies.size(); j++){
        if ( list_of_proxies[j].GetId() == GetProperties().Id() ) {
            SetFastProperties(&list_of_proxies[j]);
            return;
        }
    }
}

double SphericParticle::SlowGetYoung()                                                   { return GetProperties()[YOUNG_MODULUS];                                           }
double SphericParticle::SlowGetRollingFriction()                                         { return GetProperties()[ROLLING_FRICTION_OPTION];                                 }
double SphericParticle::SlowGetPoisson()                                                 { return GetProperties()[POISSON_RATIO];                                           }
double SphericParticle::SlowGetTgOfFrictionAngle()                                       { return GetProperties()[PARTICLE_FRICTION];                                       }
double SphericParticle::SlowGetCoefficientOfRestitution()                                { return GetProperties()[COEFFICIENT_OF_RESTITUTION];                              }
double SphericParticle::SlowGetLnOfRestitCoeff()                                         { return GetProperties()[LN_OF_RESTITUTION_COEFF];                                 }
double SphericParticle::SlowGetDensity()                                                 { return GetProperties()[PARTICLE_DENSITY];                                        }
int    SphericParticle::SlowGetParticleMaterial()                                        { return GetProperties()[PARTICLE_MATERIAL];                                       }
double SphericParticle::SlowGetParticleCohesion()                                        { return GetProperties()[PARTICLE_COHESION];                                       }

}  // namespace Kratos.

