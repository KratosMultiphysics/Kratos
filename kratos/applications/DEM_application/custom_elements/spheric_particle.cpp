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

#include <fstream>

// External includes

// Project includes
#include "spheric_particle.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "custom_utilities/discrete_particle_configure.h"


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

    if (mStressTensor!=NULL) {
        delete mStressTensor;
        mStressTensor = NULL;
        delete mSymmStressTensor;
        mSymmStressTensor = NULL;
    }
}

void SphericParticle::Initialize(const ProcessInfo& r_process_info)
{
    KRATOS_TRY
    
    SetValue(NEIGHBOUR_IDS, boost::numeric::ublas::vector<int>());
    
    MemberDeclarationFirstStep(r_process_info);

    NodeType& node = GetGeometry()[0];

    SetRadius(node.GetSolutionStepValue(RADIUS));
    SetMass(GetDensity() * CalculateVolume());

    if (this->IsNot(BLOCKED)) node.GetSolutionStepValue(PARTICLE_MATERIAL) = GetParticleMaterial();

    mClusterId = -1;

    if (this->Is(DEMFlags::HAS_ROTATION)) {
        node.GetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA) = CalculateMomentOfInertia();
    }

    else {
        array_1d<double, 3>& angular_velocity = node.GetSolutionStepValue(ANGULAR_VELOCITY);
        angular_velocity = ZeroVector(3);
    }

    if (node.GetDof(VELOCITY_X).IsFixed())         {node.Set(DEMFlags::FIXED_VEL_X,true);}
    else                                           {node.Set(DEMFlags::FIXED_VEL_X,false);}
    if (node.GetDof(VELOCITY_Y).IsFixed())         {node.Set(DEMFlags::FIXED_VEL_Y,true);}
    else                                           {node.Set(DEMFlags::FIXED_VEL_Y,false);}
    if (node.GetDof(VELOCITY_Z).IsFixed())         {node.Set(DEMFlags::FIXED_VEL_Z,true);}
    else                                           {node.Set(DEMFlags::FIXED_VEL_Z,false);}
    if (node.GetDof(ANGULAR_VELOCITY_X).IsFixed()) {node.Set(DEMFlags::FIXED_ANG_VEL_X,true);}
    else                                           {node.Set(DEMFlags::FIXED_ANG_VEL_X,false);}
    if (node.GetDof(ANGULAR_VELOCITY_Y).IsFixed()) {node.Set(DEMFlags::FIXED_ANG_VEL_Y,true);}
    else                                           {node.Set(DEMFlags::FIXED_ANG_VEL_Y,false);}
    if (node.GetDof(ANGULAR_VELOCITY_Z).IsFixed()) {node.Set(DEMFlags::FIXED_ANG_VEL_Z,true);}
    else                                           {node.Set(DEMFlags::FIXED_ANG_VEL_Z,false);}
    
    double& elastic_energy = this->GetElasticEnergy();
    elastic_energy = 0.0;
    double& inelastic_frictional_energy = this->GetInelasticFrictionalEnergy();
    inelastic_frictional_energy = 0.0;
    double& inelastic_viscodamping_energy = this->GetInelasticViscodampingEnergy();
    inelastic_viscodamping_energy = 0.0;

    mBoundDeltaDispSq = 0.0;

    CreateDiscontinuumConstitutiveLaws(r_process_info);
    KRATOS_CATCH( "" )
}

void SphericParticle::CalculateRightHandSide(ProcessInfo& r_process_info, double dt, const array_1d<double,3>& gravity, int search_control)
{
    KRATOS_TRY
    array_1d<double, 3> additional_forces;
    array_1d<double, 3> additionally_applied_moment;
    array_1d<double, 3> initial_rotation_moment;
    array_1d<double, 3>& elastic_force       = this->GetGeometry()[0].FastGetSolutionStepValue(ELASTIC_FORCES);
    array_1d<double, 3>& contact_force       = this->GetGeometry()[0].FastGetSolutionStepValue(CONTACT_FORCES);
    array_1d<double, 3>& rigid_element_force = this->GetGeometry()[0].FastGetSolutionStepValue(RIGID_ELEMENT_FORCE);

    mContactMoment.clear();
    additional_forces.clear();
    additionally_applied_moment.clear();
    initial_rotation_moment.clear();
    elastic_force.clear();
    contact_force.clear();
    rigid_element_force.clear();
    
    if (this->Is(DEMFlags::HAS_ROTATION)) {
        if (this->Is(DEMFlags::HAS_ROLLING_FRICTION)) {
            array_1d<double, 3>& rolling_resistance_moment = this->GetGeometry()[0].FastGetSolutionStepValue(ROLLING_RESISTANCE_MOMENT);
            rolling_resistance_moment.clear();
        }
    }

    bool multi_stage_RHS = false;

    InitializeForceComputation(r_process_info);

    ComputeBallToBallContactForce(elastic_force, contact_force, initial_rotation_moment, r_process_info, dt, multi_stage_RHS);

    ComputeBallToRigidFaceContactForce(elastic_force, contact_force, initial_rotation_moment, rigid_element_force, r_process_info, dt, search_control);

    if (this->IsNot(DEMFlags::BELONGS_TO_A_CLUSTER)){
        ComputeAdditionalForces(additional_forces, additionally_applied_moment, r_process_info, gravity);
    }

    array_1d<double,3>& total_forces = this->GetGeometry()[0].FastGetSolutionStepValue(TOTAL_FORCES);
    array_1d<double,3>& total_moment = this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT);

    total_forces[0] =  contact_force[0] + additional_forces[0];
    total_forces[1] =  contact_force[1] + additional_forces[1];
    total_forces[2] =  contact_force[2] + additional_forces[2];

    total_moment[0] = mContactMoment[0] + additionally_applied_moment[0];
    total_moment[1] = mContactMoment[1] + additionally_applied_moment[1];
    total_moment[2] = mContactMoment[2] + additionally_applied_moment[2];

    KRATOS_CATCH("")
}

void SphericParticle::InitializeForceComputation(ProcessInfo& r_process_info){};

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

    /*for (unsigned int i = 0; i < mNeighbourElements.size(); i++){
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
*/
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

        if (!DiscreteParticleConfigure<3>::GetDomainPeriodicity()){ // default infinite-domain case
            noalias(other_to_me_vect) = this->GetGeometry()[0].Coordinates() - ineighbour->GetGeometry()[0].Coordinates();
        }

        else { // periodic domain
            double my_coors[3] = {this->GetGeometry()[0].Coordinates()[0], this->GetGeometry()[0].Coordinates()[1], this->GetGeometry()[0].Coordinates()[2]};
            double other_coors[3] = {ineighbour->GetGeometry()[0].Coordinates()[0], ineighbour->GetGeometry()[0].Coordinates()[1], ineighbour->GetGeometry()[0].Coordinates()[2]};
            DiscreteParticleConfigure<3>::TransformToClosestPeriodicCoordinates(my_coors, other_coors);
            other_to_me_vect[0] = my_coors[0] - other_coors[0];
            other_to_me_vect[1] = my_coors[1] - other_coors[1];
            other_to_me_vect[2] = my_coors[2] - other_coors[2];
        }

        double other_radius                  = ineighbour->GetInteractionRadius();
        double distance                      = DEM_MODULUS_3(other_to_me_vect);
        double radius_sum                    = GetInteractionRadius() + other_radius;
        double indentation                   = radius_sum - distance;

        r_current_max_indentation = (indentation > r_current_max_indentation) ? indentation : r_current_max_indentation;
    }
}

void SphericParticle::CalculateMaxBallToFaceIndentation(double& r_current_max_indentation)
{
    r_current_max_indentation = - std::numeric_limits<double>::max();

    std::vector<DEMWall*>& rNeighbours   = this->mNeighbourRigidFaces;

    for (unsigned int i = 0; i < rNeighbours.size(); i++) {

        double LocalCoordSystem[3][3]            = {{0.0}, {0.0}, {0.0}};
        array_1d<double, 3> wall_delta_disp_at_contact_point = ZeroVector(3);
        array_1d<double, 3> wall_velocity_at_contact_point = ZeroVector(3);
        double DistPToB = 0.0;
        int ContactType = -1;
        array_1d<double, 4>& Weight = this->mContactConditionWeights[i];

        ComputeConditionRelativeData(i,rNeighbours[i], LocalCoordSystem, DistPToB, Weight, wall_delta_disp_at_contact_point, wall_velocity_at_contact_point, ContactType);

        if(ContactType > 0){
            double indentation = GetInteractionRadius() - DistPToB;
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

void SphericParticle::CalculateMomentum(array_1d<double, 3>& r_momentum)
{
    const array_1d<double, 3>& vel = this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
    noalias(r_momentum) = GetMass() * vel;
}

void SphericParticle::CalculateLocalAngularMomentum(array_1d<double, 3>& r_angular_momentum)
{
    const array_1d<double, 3> ang_vel  = this->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
    const double moment_of_inertia     = this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA);
    noalias(r_angular_momentum) = moment_of_inertia * ang_vel;
}

void SphericParticle::ComputeNewNeighboursHistoricalData(boost::numeric::ublas::vector<int>& mTempNeighboursIds,
                                                         std::vector<array_1d<double, 3> >& mTempNeighbourElasticContactForces)
{
    std::vector<array_1d<double, 3> > mTempNeighbourElasticExtraContactForces;
    unsigned int new_size = mNeighbourElements.size();
    array_1d<double, 3> vector_of_zeros = ZeroVector(3);
    mTempNeighboursIds.resize(new_size);
    mTempNeighbourElasticContactForces.resize(new_size);
    mTempNeighbourElasticExtraContactForces.resize(new_size);

    boost::numeric::ublas::vector<int>& vector_of_ids_of_neighbours = GetValue(NEIGHBOUR_IDS);

    for (unsigned int i = 0; i < new_size; i++) {
        noalias(mTempNeighbourElasticContactForces[i]) = vector_of_zeros;
        noalias(mTempNeighbourElasticExtraContactForces[i]) = vector_of_zeros;

        if (mNeighbourElements[i] == NULL) { // This is required by the continuum sphere which reorders the neighbors
            mTempNeighboursIds[i] = -1;
            continue;
        }

        mTempNeighboursIds[i] = mNeighbourElements[i]->Id();

        for (unsigned int j = 0; j < vector_of_ids_of_neighbours.size(); j++) {
            if (int(mTempNeighboursIds[i]) == vector_of_ids_of_neighbours[j] && vector_of_ids_of_neighbours[j] != -1) {
                noalias(mTempNeighbourElasticContactForces[i]) = mNeighbourElasticContactForces[j];
                noalias(mTempNeighbourElasticExtraContactForces[i]) = mNeighbourElasticExtraContactForces[j]; //TODO: remove this from discontinuum!!
                break;
            }
        }
    }

    vector_of_ids_of_neighbours.swap(mTempNeighboursIds);
    mNeighbourElasticContactForces.swap(mTempNeighbourElasticContactForces);
    mNeighbourElasticExtraContactForces.swap(mTempNeighbourElasticExtraContactForces);
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

        noalias(temp_neighbours_elastic_contact_forces[i]) = vector_of_zeros;
        noalias(temp_neighbours_contact_forces[i]) = vector_of_zeros;

        if (rNeighbours[i] == NULL) { // This is required by the continuum sphere which reorders the neighbors
            temp_neighbours_ids[i] = -1;
            continue;
        }

        temp_neighbours_ids[i] = static_cast<int>(rNeighbours[i]->Id());

        for (unsigned int j = 0; j != mFemOldNeighbourIds.size(); j++) {
            if (static_cast<int>(temp_neighbours_ids[i]) == mFemOldNeighbourIds[j] && mFemOldNeighbourIds[j] != -1) {
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
    rMassMatrix(0,0) = GetMass();
}

void SphericParticle::EvaluateDeltaDisplacement(double RelDeltDisp[3],
                                                double RelVel[3],
                                                double LocalCoordSystem[3][3],
                                                double OldLocalCoordSystem[3][3],
                                                const array_1d<double, 3>& other_to_me_vect,
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

void SphericParticle::RelativeDisplacementAndVelocityOfContactPointDueToRotation(const double indentation,
                                                double RelDeltDisp[3],
                                                double RelVel[3],
                                                double LocalCoordSystem[3][3],
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
    const double other_young = p_neighbour->GetYoung();
    const double my_young = GetYoung();

    const double my_arm_length = GetInteractionRadius() - indentation * other_young / (other_young + my_young);
    const double other_arm_length  = other_radius - indentation * my_young / (other_young + my_young);

    my_arm_vector[0] = -LocalCoordSystem[2][0] * my_arm_length;
    my_arm_vector[1] = -LocalCoordSystem[2][1] * my_arm_length;
    my_arm_vector[2] = -LocalCoordSystem[2][2] * my_arm_length;

    GeometryFunctions::CrossProduct(my_ang_vel, my_arm_vector, my_vel_at_contact_point_due_to_rotation);

    other_arm_vector[0] = LocalCoordSystem[2][0] * other_arm_length;
    other_arm_vector[1] = LocalCoordSystem[2][1] * other_arm_length;
    other_arm_vector[2] = LocalCoordSystem[2][2] * other_arm_length;

    GeometryFunctions::CrossProduct(other_ang_vel, other_arm_vector, other_vel_at_contact_point_due_to_rotation);

    RelVel[0] += my_vel_at_contact_point_due_to_rotation[0] - other_vel_at_contact_point_due_to_rotation[0];
    RelVel[1] += my_vel_at_contact_point_due_to_rotation[1] - other_vel_at_contact_point_due_to_rotation[1];
    RelVel[2] += my_vel_at_contact_point_due_to_rotation[2] - other_vel_at_contact_point_due_to_rotation[2];

    GeometryFunctions::CrossProduct(my_delta_rotation,    my_arm_vector,    my_delta_disp_at_contact_point_due_to_rotation);
    GeometryFunctions::CrossProduct(other_delta_rotation, other_arm_vector, other_delta_disp_at_contact_point_due_to_rotation);


    // Contribution of the rotation
    RelDeltDisp[0] += my_delta_disp_at_contact_point_due_to_rotation[0] - other_delta_disp_at_contact_point_due_to_rotation[0];
    RelDeltDisp[1] += my_delta_disp_at_contact_point_due_to_rotation[1] - other_delta_disp_at_contact_point_due_to_rotation[1];
    RelDeltDisp[2] += my_delta_disp_at_contact_point_due_to_rotation[2] - other_delta_disp_at_contact_point_due_to_rotation[2];
}


void SphericParticle::RelativeDisplacementAndVelocityOfContactPointDueToRotationMatrix(double DeltDisp[3],
                                                double RelVel[3],
                                                double OldLocalCoordSystem[3][3],
                                                const double& other_radius,
                                                const double& dt,
                                                const array_1d<double, 3>& angular_vel,
                                                SphericParticle* p_neighbour)
{
        const array_1d<double, 3>& neigh_angular_vel = p_neighbour->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
        const double other_young = p_neighbour->GetYoung();
        const double my_young = GetYoung();


        array_1d<double, 3> temp_angular_vel;
        noalias(temp_angular_vel) = angular_vel;
        array_1d<double, 3> temp_neigh_angular_vel;
        noalias(temp_neigh_angular_vel) = neigh_angular_vel;
        DEM_MULTIPLY_BY_SCALAR_3(temp_angular_vel, dt); //THIS MACRO CONVERTS THE temp_angular_vel VARIABLE INTO AN ANGLE, DESPITE THE NAME
        DEM_MULTIPLY_BY_SCALAR_3(temp_neigh_angular_vel, dt); //THIS MACRO CONVERTS THE temp_angular_vel VARIABLE INTO AN ANGLE, DESPITE THE NAME

        const double my_rotated_angle = DEM_MODULUS_3(temp_angular_vel);
        const double other_rotated_angle = DEM_MODULUS_3(temp_neigh_angular_vel);

        array_1d<double, 3> other_to_me_vect;
        noalias(other_to_me_vect)  = this->GetGeometry()[0].Coordinates() - p_neighbour->GetGeometry()[0].Coordinates();
        const double distance            = DEM_MODULUS_3(other_to_me_vect);
        const double radius_sum          = GetInteractionRadius() + other_radius;
        const double indentation         = radius_sum - distance;

        const double arm = GetInteractionRadius() - indentation * other_young / (other_young + my_young);
        array_1d<double, 3> e1;
        DEM_COPY_SECOND_TO_FIRST_3(e1, OldLocalCoordSystem[2]);
        DEM_MULTIPLY_BY_SCALAR_3(e1, -arm);
        array_1d<double, 3> new_axes1 = e1;

        const double other_arm = other_radius - indentation * my_young / (other_young + my_young);
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

        GeometryFunctions::CrossProduct(angular_vel, radial_vector, vel);
        GeometryFunctions::CrossProduct(neigh_angular_vel, other_radial_vector, other_vel);

        RelVel[0] += vel[0] - other_vel[0];
        RelVel[1] += vel[1] - other_vel[1];
        RelVel[2] += vel[2] - other_vel[2];

        // Contribution of the rotation velocity


        DeltDisp[0] += (new_axes1[0] - new_axes2[0]) + (e2[0] - e1[0]);
        DeltDisp[1] += (new_axes1[1] - new_axes2[1]) + (e2[1] - e1[1]);
        DeltDisp[2] += (new_axes1[2] - new_axes2[2]) + (e2[2] - e1[2]);
}

void SphericParticle::ComputeMoments(double NormalLocalElasticContactForce,
                                     double Force[3],
                                     array_1d<double, 3>& rInitialRotaMoment,
                                     double LocalCoordSystem2[3],
                                     SphericParticle* p_neighbour,
                                     double indentation,
                                     bool wall)
{
    double arm_length = GetInteractionRadius() - indentation;

    if (!wall) {
        const double other_young = p_neighbour->GetYoung();
        arm_length = GetInteractionRadius() - indentation * other_young / (other_young + GetYoung());
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

        array_1d<double, 3>& rolling_resistance_moment = this->GetGeometry()[0].FastGetSolutionStepValue(ROLLING_RESISTANCE_MOMENT);
        
        double equiv_rolling_friction_coeff       = GetRollingFriction() * GetInteractionRadius();

        if (!wall) {
            double other_rolling_friction_coeff = p_neighbour->GetRollingFriction() * p_neighbour->GetInteractionRadius();
            equiv_rolling_friction_coeff = std::min(equiv_rolling_friction_coeff, other_rolling_friction_coeff);
        }

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

            if (DEM_MODULUS_3(CoordSystemMoment1) > 0.0 && DEM_MODULUS_3(CoordSystemMoment2) > 0.0) {
                GeometryFunctions::CrossProduct(CoordSystemMoment2, CoordSystemMoment1, MR);
                double AbsoluteNormalLocalElasticContactForce = fabs(NormalLocalElasticContactForce);
                DEM_MULTIPLY_BY_SCALAR_3(MR, AbsoluteNormalLocalElasticContactForce)                
            }
            
            double MR_now = DEM_MODULUS_3(MR) * equiv_rolling_friction_coeff;
            double MR_max = DEM_MODULUS_3(MaxRotaMoment);

            if (MR_max > MR_now) {
                mContactMoment[0] += MR[0] * equiv_rolling_friction_coeff;
                mContactMoment[1] += MR[1] * equiv_rolling_friction_coeff;
                mContactMoment[2] += MR[2] * equiv_rolling_friction_coeff;
                rolling_resistance_moment[0] += MR[0] * equiv_rolling_friction_coeff;
                rolling_resistance_moment[1] += MR[1] * equiv_rolling_friction_coeff;
                rolling_resistance_moment[2] += MR[2] * equiv_rolling_friction_coeff;
            }
            else {
                rolling_resistance_moment = - mContactMoment;      
                mContactMoment = - rInitialRotaMoment;
            }
        } // if (equiv_rolling_friction_coeff != 0.0)
    } // if (this->Is(DEMFlags::HAS_ROLLING_FRICTION) )

}

void SphericParticle::ComputeBallToBallContactForce(array_1d<double, 3>& r_elastic_force,
                                                    array_1d<double, 3>& r_contact_force,
                                                    array_1d<double, 3>& rInitialRotaMoment,
                                                    ProcessInfo& r_process_info,
                                                    const double dt,
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

        if (!DiscreteParticleConfigure<3>::GetDomainPeriodicity()){ // default infinite-domain case
            noalias(other_to_me_vect) = this->GetGeometry()[0].Coordinates() - ineighbour->GetGeometry()[0].Coordinates();
        }

        else { // periodic domain
            double my_coors[3] = {this->GetGeometry()[0].Coordinates()[0], this->GetGeometry()[0].Coordinates()[1], this->GetGeometry()[0].Coordinates()[2]};
            double other_coors[3] = {ineighbour->GetGeometry()[0].Coordinates()[0], ineighbour->GetGeometry()[0].Coordinates()[1], ineighbour->GetGeometry()[0].Coordinates()[2]};
            DiscreteParticleConfigure<3>::TransformToClosestPeriodicCoordinates(my_coors, other_coors);
            other_to_me_vect[0] = my_coors[0] - other_coors[0];
            other_to_me_vect[1] = my_coors[1] - other_coors[1];
            other_to_me_vect[2] = my_coors[2] - other_coors[2];
        }

        const double& other_radius = ineighbour->GetInteractionRadius();
        double distance            = DEM_MODULUS_3(other_to_me_vect);
        double radius_sum          = GetInteractionRadius() + other_radius;
        double indentation         = radius_sum - distance;

        DEM_SET_COMPONENTS_TO_ZERO_3(DeltDisp)
        DEM_SET_COMPONENTS_TO_ZERO_3(LocalDeltDisp)
        DEM_SET_COMPONENTS_TO_ZERO_3(RelVel)
        DEM_SET_COMPONENTS_TO_ZERO_3x3(LocalCoordSystem)
        DEM_SET_COMPONENTS_TO_ZERO_3x3(OldLocalCoordSystem)

        EvaluateDeltaDisplacement(DeltDisp, RelVel, LocalCoordSystem, OldLocalCoordSystem, other_to_me_vect, velocity, delta_displ, ineighbour, distance);

        if (this->Is(DEMFlags::HAS_ROTATION)) {
            RelativeDisplacementAndVelocityOfContactPointDueToRotation(indentation, DeltDisp, RelVel, LocalCoordSystem, other_radius, dt, ang_velocity, ineighbour);
        }

        RelativeDisplacementAndVelocityOfContactPointDueToOtherReasons(r_process_info, DeltDisp, RelVel, OldLocalCoordSystem, LocalCoordSystem, ineighbour);

        double LocalContactForce[3]             = {0.0};
        double GlobalContactForce[3]            = {0.0};
        double LocalElasticContactForce[3]      = {0.0};
        double LocalElasticExtraContactForce[3] = {0.0};        
        double GlobalElasticContactForce[3]     = {0.0};
        double GlobalElasticExtraContactForce[3] = {0.0};
        double TotalGlobalElasticContactForce[3] = {0.0};
        double ViscoDampingLocalContactForce[3] = {0.0};
        double cohesive_force                   =  0.0;

        if (indentation > 0.0) {
            double OldLocalElasticContactForce[3] = {0.0};
            RotateOldContactForces(OldLocalCoordSystem, LocalCoordSystem, mNeighbourElasticContactForces[i]);// still in global coordinates
            // Here we recover the old local forces projected in the new coordinates in the way they were in the old ones; Now they will be increased if necessary
            GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, mNeighbourElasticContactForces[i], OldLocalElasticContactForce);
            GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, DeltDisp, LocalDeltDisp);
            const double previous_indentation = indentation + LocalDeltDisp[2];
            double LocalRelVel[3] = {0.0};
            GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, RelVel, LocalRelVel);
            mDiscontinuumConstitutiveLaw->CalculateForces(r_process_info, OldLocalElasticContactForce, LocalElasticContactForce, LocalDeltDisp, LocalRelVel, indentation, previous_indentation, ViscoDampingLocalContactForce, cohesive_force, this, ineighbour, sliding);
        }

        array_1d<double, 3> other_ball_to_ball_forces(3,0.0);

        ComputeOtherBallToBallForces(other_ball_to_ball_forces);

        // Transforming to global forces and adding up
        AddUpForcesAndProject(OldLocalCoordSystem, LocalCoordSystem, LocalContactForce, LocalElasticContactForce, LocalElasticExtraContactForce, GlobalContactForce,
                              GlobalElasticContactForce, GlobalElasticExtraContactForce, TotalGlobalElasticContactForce, ViscoDampingLocalContactForce, cohesive_force, other_ball_to_ball_forces, r_elastic_force, r_contact_force, i, r_process_info);
        //TODO: make different AddUpForces for continuum and discontinuum (different arguments, different operations!)
        
        // ROTATION FORCES
        if (this->Is(DEMFlags::HAS_ROTATION) && !multi_stage_RHS) {
            if (this->Is(DEMFlags::HAS_ROLLING_FRICTION) && !multi_stage_RHS) {
                const double coeff_acc      = this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA) / dt;
                noalias(rInitialRotaMoment) = coeff_acc * ang_velocity; // the moment needed to stop the spin in one time step
            }

            ComputeMoments(LocalElasticContactForce[2], GlobalContactForce, rInitialRotaMoment, LocalCoordSystem[2], ineighbour, indentation);
        }

        if (this->Is(DEMFlags::HAS_STRESS_TENSOR)) {
                AddNeighbourContributionToStressTensor(GlobalElasticContactForce, LocalCoordSystem[2], distance, radius_sum, this);
        }

    }// for each neighbor

    KRATOS_CATCH("")
}// ComputeBallToBallContactForce

void SphericParticle::ComputeBallToRigidFaceContactForce(array_1d<double, 3>& r_elastic_force,
                                                         array_1d<double, 3>& r_contact_force,
                                                         array_1d<double, 3>& rInitialRotaMoment,
                                                         array_1d<double, 3>& rigid_element_force,
                                                         ProcessInfo& r_process_info,
                                                         double mTimeStep,
                                                         int search_control)
{
    KRATOS_TRY

    RenewData();

    std::vector<DEMWall*>& rNeighbours   = this->mNeighbourRigidFaces;
    array_1d<double, 3> vel              = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
    const array_1d<double,3>& AngularVel = GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);

    for (unsigned int i = 0; i < rNeighbours.size(); i++) {

        DEMWall* wall = rNeighbours[i];
        if(wall == NULL) continue;

        double LocalElasticContactForce[3]       = {0.0};
        double GlobalElasticContactForce[3]      = {0.0};
        double ViscoDampingLocalContactForce[3]  = {0.0};
        double cohesive_force                    =  0.0;
        double LocalCoordSystem[3][3]            = {{0.0}, {0.0}, {0.0}};
        array_1d<double, 3> wall_delta_disp_at_contact_point = ZeroVector(3);
        array_1d<double, 3> wall_velocity_at_contact_point = ZeroVector(3);
        bool sliding = false;

        double ini_delta = GetInitialDeltaWithFEM(i);
        double DistPToB = 0.0;

        int ContactType = -1;
        array_1d<double, 4>& Weight = this->mContactConditionWeights[i];

        ComputeConditionRelativeData(i,rNeighbours[i], LocalCoordSystem, DistPToB, Weight, wall_delta_disp_at_contact_point, wall_velocity_at_contact_point, ContactType);

        if (ContactType == 1 || ContactType == 2 || ContactType == 3) {

            double indentation = -(DistPToB - GetInteractionRadius()) - ini_delta;
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

            const double actual_distance_to_contact_point = GetInteractionRadius() - indentation;

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

            double Norm_SQ_Delta_Disp = DEM_INNER_PRODUCT_3(DeltDisp,DeltDisp);
            if(Norm_SQ_Delta_Disp > mBoundDeltaDispSq)
            {mBoundDeltaDispSq = Norm_SQ_Delta_Disp;}

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

            if (this->Is(DEMFlags::HAS_ROTATION)) {
                if (this->Is(DEMFlags::HAS_ROLLING_FRICTION)) {
                    const double coeff_acc      = this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA) / mTimeStep;
                    noalias(rInitialRotaMoment) = coeff_acc * AngularVel; // the moment needed to stop the spin in one time step
                }

              ComputeMoments(LocalElasticContactForce[2], GlobalContactForce, rInitialRotaMoment, LocalCoordSystem[2], this, indentation, true); //WARNING: sending itself as the neighbor!!
            }

            //WEAR
            if (wall->GetProperties()[COMPUTE_WEAR]) {
                const double area              = KRATOS_M_PI * GetInteractionRadius() * GetInteractionRadius();
                const double density           = GetDensity();
                const double inverse_of_volume = 1.0 / (4.0 * 0.333333333333333 * area * GetInteractionRadius());
                ComputeWear(LocalCoordSystem, vel, DeltVel, mTimeStep, density, sliding, inverse_of_volume, LocalElasticContactForce[2], wall);
            } //wall->GetProperties()[COMPUTE_WEAR] if

            if (this->Is(DEMFlags::HAS_STRESS_TENSOR)) {
                AddWallContributionToStressTensor(GlobalElasticContactForce, LocalCoordSystem[2], actual_distance_to_contact_point, 0.0);
            }
        } //ContactType if
    } //rNeighbours.size loop

    KRATOS_CATCH("")
}// ComputeBallToRigidFaceContactForce

void SphericParticle::RenewData()
{
  //To be redefined
}
void SphericParticle::ComputeConditionRelativeData(int rigid_neighbour_index,
                                                  DEMWall* const wall,
                                            double LocalCoordSystem[3][3],
                                            double& DistPToB,
                                            array_1d<double, 4>& Weight,
                                            array_1d<double, 3>& wall_delta_disp_at_contact_point,
                                            array_1d<double, 3>& wall_velocity_at_contact_point,
                                            int& ContactType)
{
    size_t FE_size = wall->GetGeometry().size();
    std::vector< array_1d <double,3> >Coord;
    Coord.resize(FE_size, array_1d<double,3>(3,0.0) );
    std::vector<double> TempWeight;
    TempWeight.resize(FE_size);

    double total_weight = 0.0;
    int points = 0;
    int inode1 = 0, inode2 = 0;

    for (unsigned int inode = 0; inode < FE_size; inode++) {

        if (Weight[inode] > 1.0e-12) {

            for (unsigned int j = 0; j < 3; j++)
            {
                Coord[inode][j] = wall->GetGeometry()[inode].Coordinates()[j];
            }
            total_weight = total_weight + Weight[inode];
            points++;
            if (points == 1) {inode1 = inode;}
            if (points == 2) {inode2 = inode;}
        }

        if (fabs(total_weight - 1.0) < 1.0e-12) {
            break;
        }
    }

    bool contact_exists = true;
    array_1d<double, 3>& node_coordinates = this->GetGeometry()[0].Coordinates();
    double node_coor[3] = {0.0};
    DEM_COPY_SECOND_TO_FIRST_3(node_coor, node_coordinates)

    const double radius = this->GetInteractionRadius();

    if (points == 3 || points == 4)
    {
        unsigned int dummy_current_edge_index;
        contact_exists = GeometryFunctions::FacetCheck(Coord, node_coor, radius, LocalCoordSystem, DistPToB, TempWeight, dummy_current_edge_index);
        ContactType = 1;
        Weight[0]=TempWeight[0];
        Weight[1]=TempWeight[1];
        Weight[2]=TempWeight[2];
        if (points == 4)
        {
            Weight[3] = TempWeight[3];
        }
        else
        {
            Weight[3] = 0.0;
        }
    }

    else if (points == 2) {

        double eta = 0.0;
        contact_exists = GeometryFunctions::EdgeCheck(Coord[inode1], Coord[inode2], node_coor, radius, LocalCoordSystem, DistPToB, eta);

        Weight[inode1] = 1-eta;
        Weight[inode2] = eta;
        ContactType = 2;

    }

    else if (points == 1) {
        contact_exists = GeometryFunctions::VertexCheck(Coord[inode1], node_coor, radius, LocalCoordSystem, DistPToB);
        Weight[inode1] = 1.0;
        ContactType = 3;
    }

    if (contact_exists == false) {ContactType = -1;}

    for (std::size_t inode = 0; inode < FE_size; inode++) {
        noalias(wall_velocity_at_contact_point) += wall->GetGeometry()[inode].FastGetSolutionStepValue(VELOCITY) * Weight[inode];

        array_1d<double, 3>  wall_delta_displacement = ZeroVector(3);
        wall->GetDeltaDisplacement(wall_delta_displacement, inode);
        noalias(wall_delta_disp_at_contact_point) += wall_delta_displacement* Weight[inode];

    }
}//ComputeConditionRelativeData

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

    //TODO: generalize for any wall (3 or 4 nodes)
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

void SphericParticle::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& r_process_info)
{
    KRATOS_TRY

    ElementalDofList.resize(0);

    for (unsigned int i = 0; i < GetGeometry().size(); i++) {
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

    mRadius = this->GetGeometry()[0].FastGetSolutionStepValue(RADIUS); //Just in case someone is overwriting the radius in Python
    mPartialRepresentativeVolume = 0.0;
    this->GetGeometry()[0].FastGetSolutionStepValue(REPRESENTATIVE_VOLUME) = CalculateVolume();
    double& elastic_energy = this->GetElasticEnergy();
    elastic_energy = 0.0;
    if (this->Is(DEMFlags::HAS_STRESS_TENSOR)) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                (*mStressTensor)(i,j) = 0.0;
            }
        }
    }

    KRATOS_CATCH("")
}

void SphericParticle::AddNeighbourContributionToStressTensor(const double Force[3],
                                                             const double other_to_me_vect[3],
                                                             const double distance,
                                                             const double radius_sum,
                                                             SphericParticle* element) {
    KRATOS_TRY

    double gap = distance - radius_sum;
    double real_distance = GetInteractionRadius() + 0.5 * gap;
    
    // Esto ya estaba comentado!!!
    //double& rRepresentative_Volume = this->GetGeometry()[0].FastGetSolutionStepValue(REPRESENTATIVE_VOLUME);
    //rRepresentative_Volume += 0.33333333333333 * (real_distance * contact_area);

    array_1d<double, 3> normal_vector_on_contact;
    normal_vector_on_contact[0] = -1 * other_to_me_vect[0]; //outwards
    normal_vector_on_contact[1] = -1 * other_to_me_vect[1]; //outwards
    normal_vector_on_contact[2] = -1 * other_to_me_vect[2]; //outwards

    array_1d<double, 3> x_centroid = real_distance * normal_vector_on_contact;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            (*mStressTensor)(i,j) += x_centroid[j] * Force[i]; //ref: Katalin Bagi 1995 Mean stress tensor
        }
    }
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

void SphericParticle::CorrectRepresentativeVolume(double& rRepresentative_Volume/*, bool& is_smaller_than_sphere*/) {

    KRATOS_TRY

    const double sphere_volume = CalculateVolume();

    if ((rRepresentative_Volume <= sphere_volume)) { //In case it gets 0.0 (discontinuum). Also sometimes the error can be too big. This puts some bound to the error for continuum.
        rRepresentative_Volume = sphere_volume;
        //is_smaller_than_sphere = true;
    }

    KRATOS_CATCH("")
}

void SphericParticle::FinalizeSolutionStep(ProcessInfo& r_process_info){

    KRATOS_TRY

    ComputeReactions();
	
    this->GetGeometry()[0].FastGetSolutionStepValue(REPRESENTATIVE_VOLUME) = mPartialRepresentativeVolume;
    double& rRepresentative_Volume = this->GetGeometry()[0].FastGetSolutionStepValue(REPRESENTATIVE_VOLUME);
    
    //bool is_smaller_than_sphere = false;
    CorrectRepresentativeVolume(rRepresentative_Volume/*, is_smaller_than_sphere*/);

    if (this->Is(DEMFlags::HAS_STRESS_TENSOR)) {

        //Divide Stress Tensor by the total volume:
        //const array_1d<double, 3>& reaction_force=this->GetGeometry()[0].FastGetSolutionStepValue(FORCE_REACTION);

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                (*mStressTensor)(i,j) /= rRepresentative_Volume;
            }
            //(*mStressTensor)(i,i) += GeometryFunctions::sign( (*mStressTensor)(i,i) ) * GetRadius() * fabs(reaction_force[i]) / rRepresentative_Volume;
        }
        /*if( this->Is(DEMFlags::HAS_ROTATION) ) { //THIS IS GIVING STABILITY PROBLEMS WHEN USING THE EXTRA TERMS FOR CONTINUUM
            const array_1d<double, 3>& reaction_moment=this->GetGeometry()[0].FastGetSolutionStepValue(MOMENT_REACTION);
            const double fabs_reaction_moment_modulus = fabs( DEM_MODULUS_3(reaction_moment) );
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    if(i!=j){
                        (*mStressTensor)(i,j) += GeometryFunctions::sign( (*mStressTensor)(i,j) ) * fabs_reaction_moment_modulus / rRepresentative_Volume;
                    }
                }
            }
        }*/

        SymmetrizeStressTensor();
    }
    KRATOS_CATCH("")
}

void SphericParticle::SymmetrizeStressTensor(){
    //The following operation symmetrizes the tensor. We will work with the symmetric stress tensor always, because the non-symmetric one is being filled while forces are being calculated
    for (int i = 0; i < 3; i++) {
        for (int j = i; j < 3; j++) {
            if(fabs((*mStressTensor)(i,j)) > fabs((*mStressTensor)(j,i))) {
                (*mSymmStressTensor)(i,j) = (*mSymmStressTensor)(j,i) = (*mStressTensor)(i,j);
            }
            else {
                (*mSymmStressTensor)(i,j) = (*mSymmStressTensor)(j,i) = (*mStressTensor)(j,i);
            }
        }
    }

    /*for (int i = 0; i < 3; i++) {
        for (int j = i; j < 3; j++) {
            (*mSymmStressTensor)(i,j) = (*mSymmStressTensor)(j,i) = 0.5 * ((*mStressTensor)(i,j) + (*mStressTensor)(j,i));
        }
    }*/
}

void SphericParticle::ComputeReactions(){
    KRATOS_TRY
    if (this->Is(DEMFlags::HAS_STRESS_TENSOR)) {
        Node<3>& node = GetGeometry()[0];
        array_1d<double, 3>& reaction_force=node.FastGetSolutionStepValue(FORCE_REACTION);
        array_1d<double, 3>& r_total_forces = node.FastGetSolutionStepValue(TOTAL_FORCES);
        reaction_force[0] = node.Is(DEMFlags::FIXED_VEL_X) * (-r_total_forces[0]);
        reaction_force[1] = node.Is(DEMFlags::FIXED_VEL_Y) * (-r_total_forces[1]);
        reaction_force[2] = node.Is(DEMFlags::FIXED_VEL_Z) * (-r_total_forces[2]);

        if( this->Is(DEMFlags::HAS_ROTATION) ) {
            array_1d<double, 3>& reaction_moment=this->GetGeometry()[0].FastGetSolutionStepValue(MOMENT_REACTION);
            array_1d<double, 3>& r_total_moment = this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT);
            reaction_moment[0] = node.Is(DEMFlags::FIXED_ANG_VEL_X) * (-r_total_moment[0]);
            reaction_moment[1] = node.Is(DEMFlags::FIXED_ANG_VEL_Y) * (-r_total_moment[1]);
            reaction_moment[2] = node.Is(DEMFlags::FIXED_ANG_VEL_Z) * (-r_total_moment[2]);
        }
    }
    KRATOS_CATCH("")
}

void SphericParticle::PrepareForPrinting(ProcessInfo& r_process_info){

    if (this->Is(DEMFlags::PRINT_STRESS_TENSOR)) {
        this->GetGeometry()[0].FastGetSolutionStepValue(DEM_STRESS_TENSOR) = (*mSymmStressTensor);
    }
}

void SphericParticle::ComputeAdditionalForces(array_1d<double, 3>& externally_applied_force,
                                              array_1d<double, 3>& externally_applied_moment,
                                              const ProcessInfo& r_process_info,
                                              const array_1d<double,3>& gravity)
{
    KRATOS_TRY
    noalias(externally_applied_force)  += ComputeWeight(gravity, r_process_info);
    noalias(externally_applied_force)  += this->GetGeometry()[0].FastGetSolutionStepValue(EXTERNAL_APPLIED_FORCE);
    noalias(externally_applied_moment) += this->GetGeometry()[0].FastGetSolutionStepValue(EXTERNAL_APPLIED_MOMENT);
    KRATOS_CATCH("")
}

array_1d<double,3> SphericParticle::ComputeWeight(const array_1d<double,3>& gravity, const ProcessInfo& r_process_info) {
    
    KRATOS_TRY
    
    return GetMass() * gravity;
    
    KRATOS_CATCH("")
}

void SphericParticle::AddUpForcesAndProject(double OldCoordSystem[3][3],
                                            double LocalCoordSystem[3][3],
                                            double LocalContactForce[3],
                                            double LocalElasticContactForce[3],
                                            double LocalElasticExtraContactForce[3],
                                            double GlobalContactForce[3],
                                            double GlobalElasticContactForce[3],
                                            double GlobalElasticExtraContactForce[3],
                                            double TotalGlobalElasticContactForce[3],
                                            double ViscoDampingLocalContactForce[3],
                                            const double cohesive_force,
                                            array_1d<double, 3>& other_ball_to_ball_forces,
                                            array_1d<double, 3>& r_elastic_force,
                                            array_1d<double, 3>& r_contact_force,
                                            const unsigned int i_neighbour_count,
                                            ProcessInfo& r_process_info)
{

    for (unsigned int index = 0; index < 3; index++) {
        LocalContactForce[index] = LocalElasticContactForce[index] + ViscoDampingLocalContactForce[index] + other_ball_to_ball_forces[index];
    }
    LocalContactForce[2] -= cohesive_force;

    DEM_ADD_SECOND_TO_FIRST(LocalElasticContactForce, other_ball_to_ball_forces);

    GeometryFunctions::VectorLocal2Global(LocalCoordSystem, LocalElasticContactForce, GlobalElasticContactForce);
    GeometryFunctions::VectorLocal2Global(LocalCoordSystem, LocalContactForce, GlobalContactForce);
    GeometryFunctions::VectorLocal2Global(LocalCoordSystem, LocalElasticExtraContactForce, GlobalElasticExtraContactForce);

    // Saving contact forces (We need to, since tangential elastic force is history-dependent)
    DEM_COPY_SECOND_TO_FIRST_3(mNeighbourElasticContactForces[i_neighbour_count], GlobalElasticContactForce)
    DEM_COPY_SECOND_TO_FIRST_3(mNeighbourElasticExtraContactForces[i_neighbour_count], GlobalElasticExtraContactForce)
                   
    TotalGlobalElasticContactForce[0] = GlobalElasticContactForce[0] + GlobalElasticExtraContactForce[0];
    TotalGlobalElasticContactForce[1] = GlobalElasticContactForce[1] + GlobalElasticExtraContactForce[1];
    TotalGlobalElasticContactForce[2] = GlobalElasticContactForce[2] + GlobalElasticExtraContactForce[2];
    DEM_ADD_SECOND_TO_FIRST(r_elastic_force, TotalGlobalElasticContactForce)
            
    double TotalGlobalContactForce[3];
    TotalGlobalContactForce[0] = GlobalContactForce[0] + GlobalElasticExtraContactForce[0];
    TotalGlobalContactForce[1] = GlobalContactForce[1] + GlobalElasticExtraContactForce[1];
    TotalGlobalContactForce[2] = GlobalContactForce[2] + GlobalElasticExtraContactForce[2];    
    DEM_ADD_SECOND_TO_FIRST(r_contact_force, TotalGlobalContactForce )
}

void SphericParticle::AddUpMomentsAndProject(double LocalCoordSystem[3][3],
                                             double LocalElasticRotationalMoment[3],
                                             double LocalViscoRotationalMoment[3]) {

    double LocalContactRotationalMoment[3] = {0.0};
    double GlobalContactRotationalMoment[3] = {0.0};

    for (unsigned int index = 0; index < 3; index++) {
        LocalContactRotationalMoment[index] = LocalElasticRotationalMoment[index] + LocalViscoRotationalMoment[index];
    }

    GeometryFunctions::VectorLocal2Global(LocalCoordSystem, LocalContactRotationalMoment, GlobalContactRotationalMoment);

    DEM_ADD_SECOND_TO_FIRST(mContactMoment, GlobalContactRotationalMoment)
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

    if (r_process_info[CRITICAL_TIME_OPTION])    this->Set(DEMFlags::HAS_CRITICAL_TIME, true);   //obsolete
    else                                         this->Set(DEMFlags::HAS_CRITICAL_TIME, false);

    if (r_process_info[COMPUTE_STRESS_TENSOR_OPTION]) this->Set(DEMFlags::HAS_STRESS_TENSOR, true);
    else                                              this->Set(DEMFlags::HAS_STRESS_TENSOR, false);

    if (r_process_info[PRINT_STRESS_TENSOR_OPTION]) this->Set(DEMFlags::PRINT_STRESS_TENSOR, true);
    else                                            this->Set(DEMFlags::PRINT_STRESS_TENSOR, false);
    
    if (this->Is(DEMFlags::HAS_STRESS_TENSOR)) {

        mStressTensor  = new Matrix(3,3);
        *mStressTensor = ZeroMatrix(3,3);

        mSymmStressTensor  = new Matrix(3,3);
        *mSymmStressTensor = ZeroMatrix(3,3);
    }
    else {

        mStressTensor     = NULL;
        mSymmStressTensor = NULL;
    }
}



double SphericParticle::CalculateLocalMaxPeriod(const bool has_mpi, const ProcessInfo& r_process_info) {
    KRATOS_TRY

    double max_sqr_period = 0.0;
    for (unsigned int i = 0; i < mNeighbourElements.size(); i++) {
         double sqr_period_discontinuum = mDiscontinuumConstitutiveLaw->LocalPeriod(i, this, mNeighbourElements[i]);
         if (sqr_period_discontinuum > max_sqr_period) { (max_sqr_period = sqr_period_discontinuum); }
    }

    return max_sqr_period;

    KRATOS_CATCH("")
}


void SphericParticle::ComputeOtherBallToBallForces(array_1d<double, 3>& other_ball_to_ball_forces) {}


double SphericParticle::GetInitialDeltaWithFEM(int index) //only available in continuum_particle
{
    return 0.0;
}

void SphericParticle::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& r_process_info)
{
    KRATOS_TRY

    //CRITICAL DELTA CALCULATION

    if (rVariable == DELTA_TIME) {
        double mass = GetMass();
        double coeff = r_process_info[NODAL_MASS_COEFF];

        if (coeff > 1.0) {
            KRATOS_THROW_ERROR(std::runtime_error, "The coefficient assigned for virtual mass is larger than one. Virtual_mass_coeff is ", coeff);
        }

        else if ((coeff == 1.0) && (r_process_info[VIRTUAL_MASS_OPTION])) {
            Output = 9.0E09;
        }

        else {

            if (r_process_info[VIRTUAL_MASS_OPTION]) {
                mass = mass / (1 - coeff);
            }

            double eq_mass = 0.5 * mass; //"mass" of the contact

            double kn = 0.0;
            double kt = 0.0;

            double ini_delta = 0.05 * GetInteractionRadius(); // Hertz needs an initial Delta, linear ignores it

            mDiscontinuumConstitutiveLaw->GetContactStiffness(this, this, ini_delta, kn, kt);

            //double K = KRATOS_M_PI * GetYoung() * GetRadius(); //M. Error, should be the same that the local definition.

            Output = 0.34 * sqrt(eq_mass / kn);

            if (this->Is(DEMFlags::HAS_ROTATION)) {
                //Output *= 0.5; //factor for critical time step when rotation is allowed.
            }
        }

        return;
    }

    if (rVariable == PARTICLE_TRANSLATIONAL_KINEMATIC_ENERGY) {

      const array_1d<double, 3>& vel = this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
      double square_of_celerity      = vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2];

      Output = 0.5 * (GetMass() * square_of_celerity);
      
      return;
    }
    
    if (rVariable == PARTICLE_ROTATIONAL_KINEMATIC_ENERGY) {

        const array_1d<double, 3> ang_vel = this->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
        const double moment_of_inertia    = this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA);
        double square_of_angular_celerity = ang_vel[0] * ang_vel[0] + ang_vel[1] * ang_vel[1] + ang_vel[2] * ang_vel[2];
        Output = 0.5 * moment_of_inertia * square_of_angular_celerity;
      
        return;
    }
    
    if (rVariable == PARTICLE_ELASTIC_ENERGY) {
        
        Output = GetElasticEnergy();
        
    }

    if (rVariable == PARTICLE_INELASTIC_FRICTIONAL_ENERGY) {
        
        Output = GetInelasticFrictionalEnergy();
        
    }
    
    if (rVariable == PARTICLE_INELASTIC_VISCODAMPING_ENERGY) {
        
        Output = GetInelasticViscodampingEnergy();
        
    }

    AdditionalCalculate(rVariable, Output, r_process_info);

    KRATOS_CATCH("")

} //Calculate

void SphericParticle::Calculate(const Variable<array_1d<double, 3> >& rVariable, array_1d<double, 3>& Output,
                                const ProcessInfo& r_process_info)
{
    if (rVariable == MOMENTUM) {
        CalculateMomentum(Output);
    }

    else if (rVariable == ANGULAR_MOMENTUM) {
        CalculateLocalAngularMomentum(Output);
    }
}

void SphericParticle::RotateOldContactForces(const double OldLocalCoordSystem[3][3], const double LocalCoordSystem[3][3], array_1d<double, 3>& mNeighbourElasticContactForces) {

    array_1d<double, 3> v1;
    array_1d<double, 3> v2;
    array_1d<double, 3> v3;
    array_1d<double, 3> mNeighbourElasticContactForcesFinal;

    v1[0] = OldLocalCoordSystem[2][0]; v1[1] = OldLocalCoordSystem[2][1]; v1[2] = OldLocalCoordSystem[2][2];
    v2[0] = LocalCoordSystem[2][0];    v2[1] = LocalCoordSystem[2][1];    v2[2] = LocalCoordSystem[2][2];

    GeometryFunctions::CrossProduct(v1, v2, v3);

    double v1_mod = GeometryFunctions::module(v1);
    double v2_mod = GeometryFunctions::module(v2);
    double v3_mod = GeometryFunctions::module(v3);

    double alpha = asin(v3_mod / (v1_mod * v2_mod));

    GeometryFunctions::normalize(v3);

    GeometryFunctions::RotateAVectorAGivenAngleAroundAUnitaryVector(mNeighbourElasticContactForces, v3, alpha, mNeighbourElasticContactForcesFinal);

    DEM_COPY_SECOND_TO_FIRST_3(mNeighbourElasticContactForces, mNeighbourElasticContactForcesFinal)
}

void SphericParticle::RelativeDisplacementAndVelocityOfContactPointDueToOtherReasons(const ProcessInfo& r_process_info,
                                                                                    double DeltDisp[3], //IN GLOBAL AXES
                                                                                    double RelVel[3], //IN GLOBAL AXES
                                                                                    double OldLocalCoordSystem[3][3],
                                                                                    double LocalCoordSystem[3][3],
                                                                                    SphericParticle* neighbour_iterator) {}

void SphericParticle::SendForcesToFEM(){};


void SphericParticle::Calculate(const Variable<Vector >& rVariable, Vector& Output, const ProcessInfo& r_process_info){}
void SphericParticle::Calculate(const Variable<Matrix >& rVariable, Matrix& Output, const ProcessInfo& r_process_info){}

void SphericParticle::AdditionalCalculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& r_process_info){}

int    SphericParticle::GetClusterId()                                                   { return mClusterId;      }
void   SphericParticle::SetClusterId(int givenId)                                        { mClusterId = givenId;   }
double SphericParticle::GetRadius()                                                      { return mRadius;         }
double SphericParticle::CalculateVolume()                                                { return 4.0 * KRATOS_M_PI_3 * mRadius * mRadius * mRadius;     }
void   SphericParticle::SetRadius(double radius)                                         { mRadius = radius;       }
void   SphericParticle::SetRadius()                                                      { mRadius = GetGeometry()[0].FastGetSolutionStepValue(RADIUS);       }
double SphericParticle::GetInteractionRadius()                                           { return mRadius;         }
void   SphericParticle::SetInteractionRadius(const double radius)                        { mRadius = radius; GetGeometry()[0].FastGetSolutionStepValue(RADIUS) = radius;}
double SphericParticle::GetSearchRadius()                                                { return mSearchRadius;   }
double SphericParticle::GetSearchRadiusWithFem()                                         { return mSearchRadiusWithFem;   }
void   SphericParticle::SetSearchRadius(const double radius)                             { mSearchRadius = radius; }
void   SphericParticle::SetSearchRadiusWithFem(const double radius)                      { mSearchRadiusWithFem = radius; }
double SphericParticle::GetMass()                                                        { return mRealMass;       }
void   SphericParticle::SetMass(double real_mass)                                        { mRealMass = real_mass;  GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) = real_mass;}
double SphericParticle::CalculateMomentOfInertia()                                       { return 0.4 * GetMass() * GetRadius() * GetRadius(); }

double SphericParticle::GetYoung()                                                       { return GetFastProperties()->GetYoung();                     }
double SphericParticle::GetRollingFriction()                                             { return GetFastProperties()->GetRollingFriction();           }
double SphericParticle::GetPoisson()                                                     { return GetFastProperties()->GetPoisson();                   }
double SphericParticle::GetTgOfFrictionAngle()                                           { return GetFastProperties()->GetTgOfFrictionAngle() ;        }
double SphericParticle::GetCoefficientOfRestitution()                                    { return GetFastProperties()->GetCoefficientOfRestitution();  }
double SphericParticle::GetLnOfRestitCoeff()                                             { return GetFastProperties()->GetLnOfRestitCoeff();           }
double SphericParticle::GetDensity()                                                     { return GetFastProperties()->GetDensity();                   }
int    SphericParticle::GetParticleMaterial()                                            { return GetFastProperties()->GetParticleMaterial();          }
double SphericParticle::GetParticleCohesion()                                            { return GetFastProperties()->GetParticleCohesion();          }
double SphericParticle::GetParticleKNormal()                                             { return GetFastProperties()->GetParticleKNormal();           }
double SphericParticle::GetParticleKTangential()                                         { return GetFastProperties()->GetParticleKTangential();       }

// Conical damage
double SphericParticle::GetParticleContactRadius()                                       { return GetFastProperties()->GetParticleContactRadius();     }
double SphericParticle::GetParticleMaxStress()                                           { return GetFastProperties()->GetParticleMaxStress();         }
double SphericParticle::GetParticleAlpha()                                               { return GetFastProperties()->GetParticleAlpha();             }
double SphericParticle::GetParticleGamma()                                               { return GetFastProperties()->GetParticleGamma();             }

array_1d<double, 3>& SphericParticle::GetForce()                                         { return GetGeometry()[0].FastGetSolutionStepValue(TOTAL_FORCES);}
double&              SphericParticle::GetElasticEnergy()                                 { return mElasticEnergy; }
double&              SphericParticle::GetInelasticFrictionalEnergy()                     { return mInelasticFrictionalEnergy; }
double&              SphericParticle::GetInelasticViscodampingEnergy()                   { return mInelasticViscodampingEnergy; }

void   SphericParticle::SetYoungFromProperties(double* young)                            { GetFastProperties()->SetYoungFromProperties( young);                             }
void   SphericParticle::SetRollingFrictionFromProperties(double* rolling_friction)       { GetFastProperties()->SetRollingFrictionFromProperties( rolling_friction);        }
void   SphericParticle::SetPoissonFromProperties(double* poisson)                        { GetFastProperties()->SetPoissonFromProperties( poisson);                         }
void   SphericParticle::SetTgOfFrictionAngleFromProperties(double* tg_of_friction_angle) { GetFastProperties()->SetTgOfFrictionAngleFromProperties( tg_of_friction_angle);  }
void   SphericParticle::SetCoefficientOfRestitutionFromProperties(double* coefficient_of_restitution)     { GetFastProperties()->SetCoefficientOfRestitutionFromProperties( coefficient_of_restitution);      }
void   SphericParticle::SetLnOfRestitCoeffFromProperties(double* ln_of_restit_coeff)     { GetFastProperties()->SetLnOfRestitCoeffFromProperties( ln_of_restit_coeff);      }
void   SphericParticle::SetDensityFromProperties(double* density)                        { GetFastProperties()->SetDensityFromProperties( density);                         }
void   SphericParticle::SetParticleMaterialFromProperties(int* particle_material)        { GetFastProperties()->SetParticleMaterialFromProperties( particle_material);      }
void   SphericParticle::SetParticleCohesionFromProperties(double* particle_cohesion)     { GetFastProperties()->SetParticleCohesionFromProperties( particle_cohesion);      }
void   SphericParticle::SetParticleKNormalFromProperties(double* particle_k_normal)      { GetFastProperties()->SetParticleKNormalFromProperties( particle_k_normal);       }
void   SphericParticle::SetParticleKTangentialFromProperties(double* particle_k_tangential) { GetFastProperties()->SetParticleKTangentialFromProperties( particle_k_tangential); }

// Conical damage
void   SphericParticle::SetParticleContactRadiusFromProperties(double* particle_contact_radius) { GetFastProperties()->SetParticleContactRadiusFromProperties( particle_contact_radius); }
void   SphericParticle::SetParticleMaxStressFromProperties(double* particle_max_stress)  { GetFastProperties()->SetParticleMaxStressFromProperties( particle_max_stress);   }
void   SphericParticle::SetParticleAlphaFromProperties(double* particle_alpha)           { GetFastProperties()->SetParticleAlphaFromProperties( particle_alpha);            }
void   SphericParticle::SetParticleGammaFromProperties(double* particle_gamma)           { GetFastProperties()->SetParticleGammaFromProperties( particle_gamma);            }

PropertiesProxy* SphericParticle::GetFastProperties()                                    { return mFastProperties;                                                          }
void   SphericParticle::SetFastProperties(PropertiesProxy* pProps)                       { mFastProperties = pProps;                                                        }
void   SphericParticle::SetFastProperties(std::vector<PropertiesProxy>& list_of_proxies) {
    for (unsigned int j = 0; j < list_of_proxies.size(); j++){
        if (list_of_proxies[j].GetId() == GetProperties().Id()) {
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
double SphericParticle::SlowGetDensity()                                                 { return GetProperties()[PARTICLE_DENSITY];                                        }
int    SphericParticle::SlowGetParticleMaterial()                                        { return GetProperties()[PARTICLE_MATERIAL];                                       }
double SphericParticle::SlowGetParticleCohesion()                                        { return GetProperties()[PARTICLE_COHESION];                                       }
double SphericParticle::GetBoundDeltaDispSq()                                            { return mBoundDeltaDispSq;   }

}  // namespace Kratos.
