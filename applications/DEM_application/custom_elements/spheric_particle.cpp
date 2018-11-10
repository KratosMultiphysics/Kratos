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
    mpTranslationalIntegrationScheme = NULL;
    mpRotationalIntegrationScheme = NULL;
    mpInlet = NULL;
}

SphericParticle::SphericParticle(IndexType NewId, GeometryType::Pointer pGeometry)
    : DiscreteElement(NewId, pGeometry), mRealMass(0){
    mRadius = 0;
    mRealMass = 0;
    mStressTensor = NULL;
    mSymmStressTensor = NULL;
    mpTranslationalIntegrationScheme = NULL;
    mpRotationalIntegrationScheme = NULL;
    mpInlet = NULL;
}

SphericParticle::SphericParticle(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : DiscreteElement(NewId, pGeometry, pProperties), mRealMass(0)
{
    mRadius = 0;
    mRealMass = 0;
    mStressTensor = NULL;
    mSymmStressTensor = NULL;
    mpTranslationalIntegrationScheme = NULL;
    mpRotationalIntegrationScheme = NULL;
    mpInlet = NULL;
}

SphericParticle::SphericParticle(IndexType NewId, NodesArrayType const& ThisNodes)
    : DiscreteElement(NewId, ThisNodes), mRealMass(0)
{
    mRadius = 0;
    mRealMass = 0;
    mStressTensor = NULL;
    mSymmStressTensor = NULL;
    mpTranslationalIntegrationScheme = NULL;
    mpRotationalIntegrationScheme = NULL;
    mpInlet = NULL;
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
    if (mpTranslationalIntegrationScheme!=NULL) {
        delete mpTranslationalIntegrationScheme;
    }
    if (mpRotationalIntegrationScheme!=NULL) {
        delete mpRotationalIntegrationScheme;
    }
}

SphericParticle& SphericParticle::operator=(const SphericParticle& rOther) {
    DiscreteElement::operator=(rOther);
    mElasticEnergy = rOther.mElasticEnergy;
    mInelasticFrictionalEnergy = rOther.mInelasticFrictionalEnergy;
    mInelasticViscodampingEnergy = rOther.mInelasticViscodampingEnergy;
    mNeighbourElements = rOther.mNeighbourElements;
    mContactingNeighbourIds = rOther.mContactingNeighbourIds;
    mContactingFaceNeighbourIds = rOther.mContactingFaceNeighbourIds;
    mNeighbourRigidFaces = rOther.mNeighbourRigidFaces;
    mNeighbourPotentialRigidFaces = rOther.mNeighbourPotentialRigidFaces;
    mContactConditionWeights = rOther.mContactConditionWeights;
    mContactConditionContactTypes = rOther.mContactConditionContactTypes;
    mConditionContactPoints = rOther.mConditionContactPoints;
    mNeighbourRigidFacesTotalContactForce = rOther.mNeighbourRigidFacesTotalContactForce;
    mNeighbourRigidFacesElasticContactForce = rOther.mNeighbourRigidFacesElasticContactForce;
    mNeighbourElasticContactForces = rOther.mNeighbourElasticContactForces;
    mNeighbourElasticExtraContactForces = rOther.mNeighbourElasticExtraContactForces;
    mContactMoment = rOther.mContactMoment;
    mPartialRepresentativeVolume = rOther.mPartialRepresentativeVolume; //TODO: to continuum!
    mFemOldNeighbourIds = rOther.mFemOldNeighbourIds;
    mRadius = rOther.mRadius;
    mSearchRadius = rOther.mSearchRadius;
    mRealMass = rOther.mRealMass;
    mClusterId = rOther.mClusterId;
    mGlobalDamping = rOther.mGlobalDamping;

    if(rOther.mStressTensor != NULL) {
        mStressTensor  = new Matrix(3,3);
        *mStressTensor = *rOther.mStressTensor;

        mSymmStressTensor  = new Matrix(3,3);
        *mSymmStressTensor = *rOther.mSymmStressTensor;
    }
    else {

        mStressTensor     = NULL;
        mSymmStressTensor = NULL;
    }

    mFastProperties = rOther.mFastProperties; //This might be unsafe

    DEMIntegrationScheme::Pointer& translational_integration_scheme = GetProperties()[DEM_TRANSLATIONAL_INTEGRATION_SCHEME_POINTER];
    DEMIntegrationScheme::Pointer& rotational_integration_scheme = GetProperties()[DEM_ROTATIONAL_INTEGRATION_SCHEME_POINTER];
    SetIntegrationScheme(translational_integration_scheme, rotational_integration_scheme);

    mDiscontinuumConstitutiveLaw = GetProperties()[DEM_DISCONTINUUM_CONSTITUTIVE_LAW_POINTER]->Clone();

    SetValue(WALL_POINT_CONDITION_POINTERS, std::vector<Condition*>());
    SetValue(WALL_POINT_CONDITION_ELASTIC_FORCES, std::vector<array_1d<double, 3> >());
    SetValue(WALL_POINT_CONDITION_TOTAL_FORCES, std::vector<array_1d<double, 3> >());

    return *this;
}

void SphericParticle::Initialize(const ProcessInfo& r_process_info)
{
    KRATOS_TRY

    SetValue(NEIGHBOUR_IDS, DenseVector<int>());

    MemberDeclarationFirstStep(r_process_info);

    NodeType& node = GetGeometry()[0];

    SetRadius(node.GetSolutionStepValue(RADIUS));
    SetMass(GetDensity() * CalculateVolume());

    if (this->IsNot(BLOCKED)) node.GetSolutionStepValue(PARTICLE_MATERIAL) = GetParticleMaterial();

    mClusterId = -1;

    if (this->Is(DEMFlags::HAS_ROTATION)) {
        node.GetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA) = CalculateMomentOfInertia();

        Quaternion<double  > Orientation = Quaternion<double>::Identity();
        node.GetSolutionStepValue(ORIENTATION) = Orientation;

        array_1d<double, 3> angular_momentum;
        CalculateLocalAngularMomentum(angular_momentum);
        noalias(node.GetSolutionStepValue(ANGULAR_MOMENTUM)) = angular_momentum;
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

    CreateDiscontinuumConstitutiveLaws(r_process_info);

    DEMIntegrationScheme::Pointer& translational_integration_scheme = GetProperties()[DEM_TRANSLATIONAL_INTEGRATION_SCHEME_POINTER];
    DEMIntegrationScheme::Pointer& rotational_integration_scheme = GetProperties()[DEM_ROTATIONAL_INTEGRATION_SCHEME_POINTER];
    SetIntegrationScheme(translational_integration_scheme, rotational_integration_scheme);

    SetValue(WALL_POINT_CONDITION_POINTERS, std::vector<Condition*>());
    SetValue(WALL_POINT_CONDITION_ELASTIC_FORCES, std::vector<array_1d<double, 3> >());
    SetValue(WALL_POINT_CONDITION_TOTAL_FORCES, std::vector<array_1d<double, 3> >());

    KRATOS_CATCH( "" )
}

void SphericParticle::SetIntegrationScheme(DEMIntegrationScheme::Pointer& translational_integration_scheme, DEMIntegrationScheme::Pointer& rotational_integration_scheme) {
    mpTranslationalIntegrationScheme = translational_integration_scheme->CloneRaw();
    mpRotationalIntegrationScheme = rotational_integration_scheme->CloneRaw();
}

void SphericParticle::CalculateRightHandSide(ProcessInfo& r_process_info, double dt, const array_1d<double,3>& gravity, int search_control)
{
    KRATOS_TRY

    // Creating a data buffer to store those variables that we want to reuse so that we can keep function parameter lists short

    SphericParticle::BufferPointerType p_buffer = CreateParticleDataBuffer(this); // all memory will be freed once this shared pointer goes out of scope
    ParticleDataBuffer& data_buffer = *p_buffer;
    data_buffer.SetBoundingBox(r_process_info[DOMAIN_IS_PERIODIC], r_process_info[DOMAIN_MIN_CORNER], r_process_info[DOMAIN_MAX_CORNER]);

    NodeType& this_node = GetGeometry()[0];

    data_buffer.mDt = dt;
    data_buffer.mMultiStageRHS = false;

    array_1d<double, 3> additional_forces(3, 0.0);
    array_1d<double, 3> additionally_applied_moment(3, 0.0);
    array_1d<double, 3>& elastic_force       = this_node.FastGetSolutionStepValue(ELASTIC_FORCES);
    array_1d<double, 3>& contact_force       = this_node.FastGetSolutionStepValue(CONTACT_FORCES);
    array_1d<double, 3>& rigid_element_force = this_node.FastGetSolutionStepValue(RIGID_ELEMENT_FORCE);

    mContactMoment.clear();
    elastic_force.clear();
    contact_force.clear();
    rigid_element_force.clear();

    InitializeForceComputation(r_process_info);

    double RollingResistance = 0.0;

    ComputeBallToBallContactForce(data_buffer, r_process_info, elastic_force, contact_force, RollingResistance);

    ComputeBallToRigidFaceContactForce(data_buffer, elastic_force, contact_force, RollingResistance, rigid_element_force, r_process_info, search_control);

    if (this->IsNot(DEMFlags::BELONGS_TO_A_CLUSTER)){
        ComputeAdditionalForces(additional_forces, additionally_applied_moment, r_process_info, gravity);
        #ifdef KRATOS_DEBUG
        DemDebugFunctions::CheckIfNan(additional_forces, "NAN in Additional Force in RHS of Ball");
        DemDebugFunctions::CheckIfNan(additionally_applied_moment, "NAN in Additional Torque in RHS of Ball");
        #endif
    }

    // ROLLING FRICTION
    if (this->Is(DEMFlags::HAS_ROTATION) && !data_buffer.mMultiStageRHS) {
        if (this->Is(DEMFlags::HAS_ROLLING_FRICTION) && !data_buffer.mMultiStageRHS) {
            array_1d<double, 3>& rolling_resistance_moment = this_node.FastGetSolutionStepValue(ROLLING_RESISTANCE_MOMENT);
            rolling_resistance_moment.clear();

            ComputeRollingFriction(rolling_resistance_moment, RollingResistance, data_buffer.mDt);
        }
    }

    array_1d<double,3>& total_forces = this_node.FastGetSolutionStepValue(TOTAL_FORCES);
    array_1d<double,3>& total_moment = this_node.FastGetSolutionStepValue(PARTICLE_MOMENT);

    total_forces[0] = contact_force[0] + additional_forces[0];
    total_forces[1] = contact_force[1] + additional_forces[1];
    total_forces[2] = contact_force[2] + additional_forces[2];

    total_moment[0] = mContactMoment[0] + additionally_applied_moment[0];
    total_moment[1] = mContactMoment[1] + additionally_applied_moment[1];
    total_moment[2] = mContactMoment[2] + additionally_applied_moment[2];

    ApplyGlobalDampingToContactForcesAndMoments(total_forces, total_moment);

    #ifdef KRATOS_DEBUG
    DemDebugFunctions::CheckIfNan(total_forces, "NAN in Total Forces in RHS of Ball");
    DemDebugFunctions::CheckIfNan(total_moment, "NAN in Total Torque in RHS of Ball");
    #endif

    FinalizeForceComputation(data_buffer);
    KRATOS_CATCH("")
}

void SphericParticle::InitializeForceComputation(ProcessInfo& r_process_info){}

void SphericParticle::FirstCalculateRightHandSide(ProcessInfo& r_process_info, double dt, int search_control){}

void SphericParticle::CollectCalculateRightHandSide(ProcessInfo& r_process_info){}

void SphericParticle::FinalCalculateRightHandSide(ProcessInfo& r_process_info, double dt, const array_1d<double,3>& gravity){}

void SphericParticle::CalculateMaxBallToBallIndentation(double& r_current_max_indentation, const ProcessInfo& r_process_info)
{
    r_current_max_indentation = - std::numeric_limits<double>::max();

    for (unsigned int i = 0; i < mNeighbourElements.size(); i++){
        SphericParticle* ineighbour = mNeighbourElements[i];

        array_1d<double, 3> other_to_me_vect;
        if (!r_process_info[DOMAIN_IS_PERIODIC]){ // default infinite-domain case
            noalias(other_to_me_vect) = this->GetGeometry()[0].Coordinates() - ineighbour->GetGeometry()[0].Coordinates();
        }

        else { // periodic domain
            double my_coors[3] = {this->GetGeometry()[0][0], this->GetGeometry()[0][1], this->GetGeometry()[0][2]};
            double other_coors[3] = {ineighbour->GetGeometry()[0][0], ineighbour->GetGeometry()[0][1], ineighbour->GetGeometry()[0][2]};

            TransformNeighbourCoorsToClosestInPeriodicDomain(r_process_info, my_coors, other_coors);
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

        //ComputeConditionRelativeData(i,rNeighbours[i], LocalCoordSystem, DistPToB, Weight, wall_delta_disp_at_contact_point, wall_velocity_at_contact_point, ContactType);

        rNeighbours[i]->ComputeConditionRelativeData(i, this, LocalCoordSystem, DistPToB, Weight, wall_delta_disp_at_contact_point, wall_velocity_at_contact_point, ContactType);

        if(ContactType > 0){
            double indentation = GetInteractionRadius() - DistPToB;
            r_current_max_indentation = (indentation > r_current_max_indentation) ? indentation : r_current_max_indentation;

        }

    } //for every rigidface neighbor
}

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

void SphericParticle::ComputeNewNeighboursHistoricalData(DenseVector<int>& mTempNeighboursIds,
                                                         std::vector<array_1d<double, 3> >& mTempNeighbourElasticContactForces)
{
    std::vector<array_1d<double, 3> > mTempNeighbourElasticExtraContactForces;
    unsigned int new_size = mNeighbourElements.size();
    array_1d<double, 3> vector_of_zeros = ZeroVector(3);
    mTempNeighboursIds.resize(new_size);
    mTempNeighbourElasticContactForces.resize(new_size);
    mTempNeighbourElasticExtraContactForces.resize(new_size);

    DenseVector<int>& vector_of_ids_of_neighbours = GetValue(NEIGHBOUR_IDS);

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

void SphericParticle::EvaluateDeltaDisplacement(ParticleDataBuffer & data_buffer,
                                                double RelDeltDisp[3],
                                                double RelVel[3],
                                                double LocalCoordSystem[3][3],
                                                double OldLocalCoordSystem[3][3],
                                                const array_1d<double, 3>& vel,
                                                const array_1d<double, 3>& delta_displ)
{
    // FORMING LOCAL COORDINATES

    // Notes: Since we will normally inherit the mesh from GiD, we respect the global system X,Y,Z [0],[1],[2]
    // In the local coordinates we will define the normal direction of the contact as the [2] component.
    // Compression is positive.
    GeometryFunctions::ComputeContactLocalCoordSystem(data_buffer.mOtherToMeVector, data_buffer.mDistance, LocalCoordSystem); //new Local Coordinate System (normalizes data_buffer.mOtherToMeVector)

    // FORMING OLD LOCAL COORDINATES
    array_1d<double, 3> old_coord_target;
    noalias(old_coord_target) = this->GetGeometry()[0].Coordinates() - delta_displ;

    const array_1d<double, 3>& other_delta_displ = data_buffer.mpOtherParticleNode->FastGetSolutionStepValue(DELTA_DISPLACEMENT);
    array_1d<double, 3> old_coord_neigh;
    noalias(old_coord_neigh) = data_buffer.mpOtherParticleNode->Coordinates() - other_delta_displ;

    if (data_buffer.mDomainIsPeriodic){ // the domain is periodic
        TransformNeighbourCoorsToClosestInPeriodicDomain(data_buffer, old_coord_target, old_coord_neigh);
    }

    array_1d<double, 3> old_other_to_me_vect;
    noalias(old_other_to_me_vect) = old_coord_target - old_coord_neigh;

    const double old_distance = DEM_MODULUS_3(old_other_to_me_vect);

    GeometryFunctions::ComputeContactLocalCoordSystem(old_other_to_me_vect, old_distance, OldLocalCoordSystem); //Old Local Coordinate System

    // VELOCITIES AND DISPLACEMENTS
    const array_1d<double, 3 >& other_vel = data_buffer.mpOtherParticleNode->FastGetSolutionStepValue(VELOCITY);

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
                                                const double LocalCoordSystem[3][3],
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
                                                const double OldLocalCoordSystem[3][3],
                                                const double& other_radius,
                                                const double& dt,
                                                const array_1d<double, 3>& angular_vel,
                                                SphericParticle* p_neighbour)
{
        const array_1d<double, 3>& neigh_angular_vel = p_neighbour->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
        const array_1d<double, 3>& my_delta_rotation = GetGeometry()[0].FastGetSolutionStepValue(DELTA_ROTATION);
        const array_1d<double, 3>& other_delta_rotation = p_neighbour->GetGeometry()[0].FastGetSolutionStepValue(DELTA_ROTATION);
        const double other_young = p_neighbour->GetYoung();
        const double my_young = GetYoung();

        const double my_rotated_angle = DEM_MODULUS_3(my_delta_rotation);
        const double other_rotated_angle = DEM_MODULUS_3(other_delta_rotation);

        const array_1d<double, 3>& coors = this->GetGeometry()[0].Coordinates();
        array_1d<double, 3> neigh_coors = p_neighbour->GetGeometry()[0].Coordinates();

        const array_1d<double, 3> other_to_me_vect = coors - neigh_coors;
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
            axis_1[0] = my_delta_rotation[0] / my_rotated_angle;
            axis_1[1] = my_delta_rotation[1] / my_rotated_angle;
            axis_1[2] = my_delta_rotation[2] / my_rotated_angle;

            GeometryFunctions::RotateAVectorAGivenAngleAroundAUnitaryVector(e1, axis_1, my_rotated_angle, new_axes1);

        } //if my_rotated_angle

        if (other_rotated_angle) {
            array_1d<double, 3> axis_2;
            axis_2[0] = other_delta_rotation[0] / other_rotated_angle;
            axis_2[1] = other_delta_rotation[1] / other_rotated_angle;
            axis_2[2] = other_delta_rotation[2] / other_rotated_angle;

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

void SphericParticle::RelativeDisplacementAndVelocityOfContactPointDueToRotationQuaternion(double RelDeltDisp[3],
                                                double RelVel[3],
                                                const double OldLocalCoordSystem[3][3],
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
    array_1d<double, 3> my_new_arm_vector;
    array_1d<double, 3> other_new_arm_vector;
    array_1d<double, 3> my_vel_at_contact_point_due_to_rotation;
    array_1d<double, 3> other_vel_at_contact_point_due_to_rotation;
    const double other_young = p_neighbour->GetYoung();
    const double my_young = GetYoung();

    const array_1d<double, 3>& coors = this->GetGeometry()[0].Coordinates();
    array_1d<double, 3> neigh_coors = p_neighbour->GetGeometry()[0].Coordinates();

    const array_1d<double, 3> other_to_me_vect = coors - neigh_coors;
    const double distance            = DEM_MODULUS_3(other_to_me_vect);
    const double radius_sum          = GetInteractionRadius() + other_radius;
    const double indentation         = radius_sum - distance;

    const double my_arm_length = GetInteractionRadius() - indentation * other_young / (other_young + my_young);
    const double other_arm_length  = other_radius - indentation * my_young / (other_young + my_young);

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

    Quaternion<double> MyDeltaOrientation = Quaternion<double>::Identity();
    Quaternion<double> OtherDeltaOrientation = Quaternion<double>::Identity();

    GeometryFunctions::OrientationFromRotationAngle(MyDeltaOrientation, my_delta_rotation);
    GeometryFunctions::OrientationFromRotationAngle(OtherDeltaOrientation, other_delta_rotation);

    MyDeltaOrientation.RotateVector3(my_arm_vector, my_new_arm_vector);
    OtherDeltaOrientation.RotateVector3(other_arm_vector, other_new_arm_vector);

    // Contribution of the rotation
    RelDeltDisp[0] += (my_new_arm_vector[0] - other_new_arm_vector[0]) + (other_arm_vector[0] - my_arm_vector[0]);
    RelDeltDisp[1] += (my_new_arm_vector[1] - other_new_arm_vector[1]) + (other_arm_vector[1] - my_arm_vector[1]);
    RelDeltDisp[2] += (my_new_arm_vector[2] - other_new_arm_vector[2]) + (other_arm_vector[2] - my_arm_vector[2]);
}

void SphericParticle::ComputeMoments(double NormalLocalContactForce,
                                     double Force[3],
                                     double& RollingResistance,
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

        double equiv_rolling_friction_coeff = 0.0;

        if (!wall) {
            const double my_rolling_friction_coeff    = GetRollingFriction() * GetRadius();
            const double other_rolling_friction_coeff = p_neighbour->GetRollingFriction() * p_neighbour->GetRadius();
            equiv_rolling_friction_coeff = std::min(my_rolling_friction_coeff, other_rolling_friction_coeff);
        }

        else if (wall) {
            equiv_rolling_friction_coeff = GetRollingFrictionWithWalls() * GetRadius();
        }

        if (equiv_rolling_friction_coeff != 0.0) {
            RollingResistance += fabs(NormalLocalContactForce) * equiv_rolling_friction_coeff;
        }
    }
}

void SphericParticle::ComputeRollingFriction(array_1d<double, 3>& rolling_resistance_moment, double& RollingResistance, double dt)
{
    const double coeff_acc                            = this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA) / dt;
    const array_1d<double, 3>& ang_velocity           = this->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
    const array_1d<double, 3> initial_rotation_moment = coeff_acc * ang_velocity; // the moment needed to stop the spin in one time step

    const double MaxRotaMoment[3] = {initial_rotation_moment[0] + mContactMoment[0], initial_rotation_moment[1] + mContactMoment[1], initial_rotation_moment[2] + mContactMoment[2]};
    double CoordSystemMoment[3]   = {0.0};

    double max_rota_moment_modulus_inv = 1.0 / DEM_MODULUS_3(MaxRotaMoment);
    CoordSystemMoment[0]         = MaxRotaMoment[0] * max_rota_moment_modulus_inv;
    CoordSystemMoment[1]         = MaxRotaMoment[1] * max_rota_moment_modulus_inv;
    CoordSystemMoment[2]         = MaxRotaMoment[2] * max_rota_moment_modulus_inv;

    const double MR_now = DEM_INNER_PRODUCT_3(CoordSystemMoment, CoordSystemMoment) * RollingResistance * RollingResistance;
    const double MR_max = DEM_INNER_PRODUCT_3(MaxRotaMoment, MaxRotaMoment);

    if (MR_max > MR_now) {
        mContactMoment[0] -= CoordSystemMoment[0] * RollingResistance;
        mContactMoment[1] -= CoordSystemMoment[1] * RollingResistance;
        mContactMoment[2] -= CoordSystemMoment[2] * RollingResistance;

        rolling_resistance_moment[0] -= CoordSystemMoment[0] * RollingResistance;
        rolling_resistance_moment[1] -= CoordSystemMoment[1] * RollingResistance;
        rolling_resistance_moment[2] -= CoordSystemMoment[2] * RollingResistance;
    }
    else {
        rolling_resistance_moment = - mContactMoment;
        mContactMoment = - initial_rotation_moment;
    }
}

void SphericParticle::ComputeBallToBallContactForce(SphericParticle::ParticleDataBuffer & data_buffer,
                                                    ProcessInfo& r_process_info,
                                                    array_1d<double, 3>& r_elastic_force,
                                                    array_1d<double, 3>& r_contact_force,
                                                    double& RollingResistance)
{
    KRATOS_TRY

    NodeType& this_node = this->GetGeometry()[0];
    DEM_COPY_SECOND_TO_FIRST_3(data_buffer.mMyCoors, this_node)

    //LOOP OVER NEIGHBORS:
    for (int i = 0; data_buffer.SetNextNeighbourOrExit(i); ++i){

        if (CalculateRelativePositionsOrSkipContact(data_buffer)) {
            DEM_SET_COMPONENTS_TO_ZERO_3x3(data_buffer.mLocalCoordSystem)
            DEM_SET_COMPONENTS_TO_ZERO_3x3(data_buffer.mOldLocalCoordSystem)
            double DeltDisp[3]                       = {0.0};
            double LocalDeltDisp[3]                  = {0.0};
            double RelVel[3]                         = {0.0};
            double LocalContactForce[3]              = {0.0};
            double GlobalContactForce[3]             = {0.0};
            double LocalElasticContactForce[3]       = {0.0};
            double LocalElasticExtraContactForce[3]  = {0.0};
            double GlobalElasticContactForce[3]      = {0.0};
            double GlobalElasticExtraContactForce[3] = {0.0};
            double TotalGlobalElasticContactForce[3] = {0.0};
            double ViscoDampingLocalContactForce[3]  = {0.0};
            double cohesive_force                    =  0.0;
            bool sliding = false;

            const array_1d<double, 3>& velocity     = this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
            const array_1d<double, 3>& delta_displ  = this->GetGeometry()[0].FastGetSolutionStepValue(DELTA_DISPLACEMENT);
            const array_1d<double, 3>& ang_velocity = this->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);

            EvaluateDeltaDisplacement(data_buffer, DeltDisp, RelVel, data_buffer.mLocalCoordSystem, data_buffer.mOldLocalCoordSystem, velocity, delta_displ);

            if (this->Is(DEMFlags::HAS_ROTATION)) {
                RelativeDisplacementAndVelocityOfContactPointDueToRotationQuaternion(DeltDisp, RelVel, data_buffer.mLocalCoordSystem, data_buffer.mOtherRadius, data_buffer.mDt, ang_velocity, data_buffer.mpOtherParticle);
            }

            RelativeDisplacementAndVelocityOfContactPointDueToOtherReasons(r_process_info, DeltDisp, RelVel, data_buffer.mOldLocalCoordSystem, data_buffer.mLocalCoordSystem, data_buffer.mpOtherParticle);


            EvaluateBallToBallForcesForPositiveIndentiations(data_buffer,
                                                             r_process_info,
                                                             LocalElasticContactForce,
                                                             DeltDisp,
                                                             LocalDeltDisp,
                                                             RelVel,
                                                             data_buffer.mIndentation,
                                                             ViscoDampingLocalContactForce,
                                                             cohesive_force,
                                                             data_buffer.mpOtherParticle,
                                                             sliding,
                                                             data_buffer.mLocalCoordSystem,
                                                             data_buffer.mOldLocalCoordSystem,
                                                             mNeighbourElasticContactForces[i]);


            array_1d<double, 3> other_ball_to_ball_forces(3, 0.0);
            ComputeOtherBallToBallForces(other_ball_to_ball_forces); //These forces can exist even with no indentation.

            // Transforming to global forces and adding up
            AddUpForcesAndProject(data_buffer.mOldLocalCoordSystem, data_buffer.mLocalCoordSystem, LocalContactForce, LocalElasticContactForce, LocalElasticExtraContactForce, GlobalContactForce,
                                  GlobalElasticContactForce, GlobalElasticExtraContactForce, TotalGlobalElasticContactForce, ViscoDampingLocalContactForce, cohesive_force, other_ball_to_ball_forces, r_elastic_force, r_contact_force, i, r_process_info);
            //TODO: make different AddUpForces for continuum and discontinuum (different arguments, different operations!)

            // ROTATION FORCES
            if (this->Is(DEMFlags::HAS_ROTATION) && !data_buffer.mMultiStageRHS) {
                ComputeMoments(LocalContactForce[2], GlobalContactForce, RollingResistance, data_buffer.mLocalCoordSystem[2], data_buffer.mpOtherParticle, data_buffer.mIndentation);
            }

            if (this->Is(DEMFlags::HAS_STRESS_TENSOR)) {
                AddNeighbourContributionToStressTensor(GlobalElasticContactForce, data_buffer.mLocalCoordSystem[2], data_buffer.mDistance, data_buffer.mRadiusSum, this);
            }

            if (r_process_info[IS_TIME_TO_PRINT] && r_process_info[CONTACT_MESH_OPTION] == 1) {
                unsigned int neighbour_iterator_id = data_buffer.mpOtherParticle->Id();
                if ((i < (int)mNeighbourElements.size()) && this->Id() < neighbour_iterator_id) {
                    CalculateOnContactElements(i, LocalContactForce);
                }
            }

            DEM_SET_COMPONENTS_TO_ZERO_3(DeltDisp)
            DEM_SET_COMPONENTS_TO_ZERO_3(LocalDeltDisp)
            DEM_SET_COMPONENTS_TO_ZERO_3(RelVel)
            DEM_SET_COMPONENTS_TO_ZERO_3x3(data_buffer.mLocalCoordSystem)
            DEM_SET_COMPONENTS_TO_ZERO_3x3(data_buffer.mOldLocalCoordSystem)

            #ifdef KRATOS_DEBUG
                DemDebugFunctions::CheckIfNan(GlobalContactForce, "NAN in Force in Ball to Ball contact");
                DemDebugFunctions::CheckIfNan(mContactMoment, "NAN in Torque in Ball to Ball contact");
            #endif
        }
    }// for each neighbor

    KRATOS_CATCH("")
}// ComputeBallToBallContactForce

void SphericParticle::EvaluateBallToBallForcesForPositiveIndentiations(SphericParticle::ParticleDataBuffer & data_buffer,
                                                                       const ProcessInfo& r_process_info,
                                                                       double LocalElasticContactForce[3],
                                                                       double DeltDisp[3],
                                                                       double LocalDeltDisp[3],
                                                                       double RelVel[3],
                                                                       const double indentation,
                                                                       double ViscoDampingLocalContactForce[3],
                                                                       double& cohesive_force,
                                                                       SphericParticle* p_neighbour_element,
                                                                       bool& sliding,
                                                                       double LocalCoordSystem[3][3],
                                                                       double OldLocalCoordSystem[3][3],
                                                                       array_1d<double, 3>& neighbour_elastic_contact_force)
{
    double OldLocalElasticContactForce[3] = {0.0};
    RotateOldContactForces(OldLocalCoordSystem, LocalCoordSystem, neighbour_elastic_contact_force);// still in global coordinates
    // Here we recover the old local forces projected in the new coordinates in the way they were in the old ones; Now they will be increased if necessary
    GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, neighbour_elastic_contact_force, OldLocalElasticContactForce);
    GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, DeltDisp, LocalDeltDisp);
    const double previous_indentation = indentation + LocalDeltDisp[2];
    data_buffer.mLocalRelVel[0] = 0.0;
    data_buffer.mLocalRelVel[1] = 0.0;
    data_buffer.mLocalRelVel[2] = 0.0;
    GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, RelVel, data_buffer.mLocalRelVel);
    mDiscontinuumConstitutiveLaw->CalculateForces(r_process_info, OldLocalElasticContactForce,
            LocalElasticContactForce, LocalDeltDisp, data_buffer.mLocalRelVel, indentation, previous_indentation,
            ViscoDampingLocalContactForce, cohesive_force, this, p_neighbour_element, sliding, LocalCoordSystem);
}

void SphericParticle::ComputeBallToRigidFaceContactForce(SphericParticle::ParticleDataBuffer & data_buffer,
                                                         array_1d<double, 3>& r_elastic_force,
                                                         array_1d<double, 3>& r_contact_force,
                                                         double& RollingResistance,
                                                         array_1d<double, 3>& rigid_element_force,
                                                         ProcessInfo& r_process_info,
                                                         int search_control)
{
    KRATOS_TRY

    RenewData();

    std::vector<DEMWall*>& rNeighbours   = this->mNeighbourRigidFaces;
    array_1d<double, 3> velocity         = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
    const array_1d<double,3>& AngularVel = GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);

    for (unsigned int i = 0; i < rNeighbours.size(); i++) {
        DEMWall* wall = rNeighbours[i];
        if(wall == NULL) continue;
        if(wall->IsPhantom()){
            wall->CheckSide(this);
            continue;
        }

        double LocalElasticContactForce[3]       = {0.0};
        double GlobalElasticContactForce[3]      = {0.0};
        double ViscoDampingLocalContactForce[3]  = {0.0};
        double cohesive_force                    =  0.0;
        DEM_SET_COMPONENTS_TO_ZERO_3x3(data_buffer.mLocalCoordSystem)
        DEM_SET_COMPONENTS_TO_ZERO_3x3(data_buffer.mOldLocalCoordSystem)
        array_1d<double, 3> wall_delta_disp_at_contact_point = ZeroVector(3);
        array_1d<double, 3> wall_velocity_at_contact_point = ZeroVector(3);
        bool sliding = false;

        double ini_delta = GetInitialDeltaWithFEM(i);
        double DistPToB = 0.0;

        int ContactType = -1;
        array_1d<double, 4>& Weight = this->mContactConditionWeights[i];

        rNeighbours[i]->ComputeConditionRelativeData(i, this, data_buffer.mLocalCoordSystem, DistPToB, Weight, wall_delta_disp_at_contact_point, wall_velocity_at_contact_point, ContactType);

        if (ContactType == 1 || ContactType == 2 || ContactType == 3) {

            double indentation = -(DistPToB - GetInteractionRadius()) - ini_delta;
            double DeltDisp[3] = {0.0};
            double DeltVel [3] = {0.0};

            DeltVel[0] = velocity[0] - wall_velocity_at_contact_point[0];
            DeltVel[1] = velocity[1] - wall_velocity_at_contact_point[1];
            DeltVel[2] = velocity[2] - wall_velocity_at_contact_point[2];

            // For translation movement delta displacement
            const array_1d<double, 3>& delta_displ  = this->GetGeometry()[0].FastGetSolutionStepValue(DELTA_DISPLACEMENT);
            DeltDisp[0] = delta_displ[0] - wall_delta_disp_at_contact_point[0];
            DeltDisp[1] = delta_displ[1] - wall_delta_disp_at_contact_point[1];
            DeltDisp[2] = delta_displ[2] - wall_delta_disp_at_contact_point[2];

            if (this->Is(DEMFlags::HAS_ROTATION)) {
                const array_1d<double,3>& delta_rotation = GetGeometry()[0].FastGetSolutionStepValue(DELTA_ROTATION);

                array_1d<double, 3> actual_arm_vector, new_arm_vector;
                actual_arm_vector[0] = -data_buffer.mLocalCoordSystem[2][0] * DistPToB;
                actual_arm_vector[1] = -data_buffer.mLocalCoordSystem[2][1] * DistPToB;
                actual_arm_vector[2] = -data_buffer.mLocalCoordSystem[2][2] * DistPToB;

                double tangential_vel[3]           = {0.0};
                double tangential_displacement_due_to_rotation[3]  = {0.0};
                GeometryFunctions::CrossProduct(AngularVel, actual_arm_vector, tangential_vel);

                Quaternion<double> DeltaOrientation = Quaternion<double>::Identity();
                GeometryFunctions::OrientationFromRotationAngle(DeltaOrientation, delta_rotation);

                DeltaOrientation.RotateVector3(actual_arm_vector, new_arm_vector);

                // Contribution of the rotation
                tangential_displacement_due_to_rotation[0] = (new_arm_vector[0] - actual_arm_vector[0]);
                tangential_displacement_due_to_rotation[1] = (new_arm_vector[1] - actual_arm_vector[1]);
                tangential_displacement_due_to_rotation[2] = (new_arm_vector[2] - actual_arm_vector[2]);

                DEM_ADD_SECOND_TO_FIRST(DeltVel, tangential_vel)
                DEM_ADD_SECOND_TO_FIRST(DeltDisp, tangential_displacement_due_to_rotation)
            }

            double LocalDeltDisp[3] = {0.0};
            GeometryFunctions::VectorGlobal2Local(data_buffer.mLocalCoordSystem, DeltDisp, LocalDeltDisp);

            double OldLocalElasticContactForce[3] = {0.0};

            GeometryFunctions::VectorGlobal2Local(data_buffer.mLocalCoordSystem, mNeighbourRigidFacesElasticContactForce[i], OldLocalElasticContactForce);
            const double previous_indentation = indentation + LocalDeltDisp[2];
            data_buffer.mLocalRelVel[0] = 0.0;
            data_buffer.mLocalRelVel[1] = 0.0;
            data_buffer.mLocalRelVel[2] = 0.0;

            if (indentation > 0.0) {

                GeometryFunctions::VectorGlobal2Local(data_buffer.mLocalCoordSystem, DeltVel, data_buffer.mLocalRelVel);
                mDiscontinuumConstitutiveLaw->CalculateForcesWithFEM(r_process_info,
                                                                    OldLocalElasticContactForce,
                                                                    LocalElasticContactForce,
                                                                    LocalDeltDisp,
                                                                    data_buffer.mLocalRelVel,
                                                                    indentation,
                                                                    previous_indentation,
                                                                    ViscoDampingLocalContactForce,
                                                                    cohesive_force,
                                                                    this,
                                                                    wall,
                                                                    sliding);
            }

            double LocalContactForce[3]  = {0.0};
            double GlobalContactForce[3] = {0.0};

            AddUpFEMForcesAndProject(data_buffer.mLocalCoordSystem, LocalContactForce, LocalElasticContactForce, GlobalContactForce,
                                     GlobalElasticContactForce, ViscoDampingLocalContactForce, cohesive_force, r_elastic_force,
                                     r_contact_force, mNeighbourRigidFacesElasticContactForce[i], mNeighbourRigidFacesTotalContactForce[i]);

            rigid_element_force[0] -= GlobalContactForce[0];
            rigid_element_force[1] -= GlobalContactForce[1];
            rigid_element_force[2] -= GlobalContactForce[2];

            if (this->Is(DEMFlags::HAS_ROTATION)) {
                ComputeMoments(LocalContactForce[2], GlobalContactForce, RollingResistance, data_buffer.mLocalCoordSystem[2], this, indentation, true); //WARNING: sending itself as the neighbor!!
            }

            //WEAR
            if (wall->GetProperties()[COMPUTE_WEAR]) {
                const double area              = Globals::Pi * GetInteractionRadius() * GetInteractionRadius();
                const double inverse_of_volume = 1.0 / (4.0 * 0.333333333333333 * area * GetInteractionRadius());
                ComputeWear(data_buffer.mLocalRelVel,
                            data_buffer.mDt,
                            sliding,
                            inverse_of_volume,
                            LocalElasticContactForce[2],
                            wall);
            } //wall->GetProperties()[COMPUTE_WEAR] if

            if (this->Is(DEMFlags::HAS_STRESS_TENSOR)) {
                AddWallContributionToStressTensor(GlobalElasticContactForce, data_buffer.mLocalCoordSystem[2], DistPToB, 0.0);
            }
        } //ContactType if
    } //rNeighbours.size loop

    auto& list_of_point_condition_pointers = this->GetValue(WALL_POINT_CONDITION_POINTERS);
    auto& neighbour_point_faces_elastic_contact_force = this->GetValue(WALL_POINT_CONDITION_ELASTIC_FORCES);
    auto& neighbour_point_faces_total_contact_force = this->GetValue(WALL_POINT_CONDITION_TOTAL_FORCES);

    for (unsigned int i=0; i<list_of_point_condition_pointers.size(); i++) {
        Condition* wall = list_of_point_condition_pointers[i];
        Node<3>& cond_node = wall->GetGeometry()[0];

        double RelVel[3]                         = {0.0};
        double LocalElasticContactForce[3]       = {0.0};
        double GlobalElasticContactForce[3]      = {0.0};
        double ViscoDampingLocalContactForce[3]  = {0.0};
        double cohesive_force                    =  0.0;
        DEM_SET_COMPONENTS_TO_ZERO_3x3(data_buffer.mLocalCoordSystem)
        DEM_SET_COMPONENTS_TO_ZERO_3x3(data_buffer.mOldLocalCoordSystem)

        array_1d<double, 3> wall_delta_disp_at_contact_point = cond_node.GetSolutionStepValue(DELTA_DISPLACEMENT);
        array_1d<double, 3> wall_velocity_at_contact_point = cond_node.GetSolutionStepValue(VELOCITY);
        bool sliding = false;
        double ini_delta = 0.0;

        array_1d<double, 3> cond_to_me_vect;
        noalias(cond_to_me_vect) = GetGeometry()[0].Coordinates() - cond_node.Coordinates();
        double DistPToB = DEM_MODULUS_3(cond_to_me_vect);;

        double indentation = -(DistPToB - GetInteractionRadius()) - ini_delta;
        double DeltDisp[3] = {0.0};
        double DeltVel [3] = {0.0};

        DeltVel[0] = velocity[0] - wall_velocity_at_contact_point[0];
        DeltVel[1] = velocity[1] - wall_velocity_at_contact_point[1];
        DeltVel[2] = velocity[2] - wall_velocity_at_contact_point[2];

        // For translation movement delta displacement
        const array_1d<double, 3>& delta_displ  = this->GetGeometry()[0].FastGetSolutionStepValue(DELTA_DISPLACEMENT);
        DeltDisp[0] = delta_displ[0] - wall_delta_disp_at_contact_point[0];
        DeltDisp[1] = delta_displ[1] - wall_delta_disp_at_contact_point[1];
        DeltDisp[2] = delta_displ[2] - wall_delta_disp_at_contact_point[2];

        data_buffer.mpOtherParticleNode = &cond_node;
        DEM_COPY_SECOND_TO_FIRST_3(data_buffer.mOtherToMeVector, cond_to_me_vect)
        data_buffer.mDistance = DistPToB;
        data_buffer.mDomainIsPeriodic = false;
        EvaluateDeltaDisplacement(data_buffer, DeltDisp, RelVel, data_buffer.mLocalCoordSystem, data_buffer.mOldLocalCoordSystem, velocity, delta_displ);

        if (this->Is(DEMFlags::HAS_ROTATION)) {
            const array_1d<double,3>& delta_rotation = GetGeometry()[0].FastGetSolutionStepValue(DELTA_ROTATION);

            array_1d<double, 3> actual_arm_vector;
            actual_arm_vector[0] = -data_buffer.mLocalCoordSystem[2][0] * DistPToB;
            actual_arm_vector[1] = -data_buffer.mLocalCoordSystem[2][1] * DistPToB;
            actual_arm_vector[2] = -data_buffer.mLocalCoordSystem[2][2] * DistPToB;

            double tangential_vel[3]           = {0.0};
            double tangential_displacement_due_to_rotation[3]  = {0.0};
            GeometryFunctions::CrossProduct(AngularVel, actual_arm_vector, tangential_vel);
            GeometryFunctions::CrossProduct(delta_rotation, actual_arm_vector, tangential_displacement_due_to_rotation);

            DEM_ADD_SECOND_TO_FIRST(DeltVel, tangential_vel)
            DEM_ADD_SECOND_TO_FIRST(DeltDisp, tangential_displacement_due_to_rotation)
        }

        double LocalDeltDisp[3] = {0.0};
        GeometryFunctions::VectorGlobal2Local(data_buffer.mLocalCoordSystem, DeltDisp, LocalDeltDisp);

        double OldLocalElasticContactForce[3] = {0.0};

        GeometryFunctions::VectorGlobal2Local(data_buffer.mLocalCoordSystem, neighbour_point_faces_elastic_contact_force[i], OldLocalElasticContactForce);
        const double previous_indentation = indentation + LocalDeltDisp[2];
        data_buffer.mLocalRelVel[0] = 0.0;
        data_buffer.mLocalRelVel[1] = 0.0;
        data_buffer.mLocalRelVel[2] = 0.0;

        if (indentation > 0.0) {

            GeometryFunctions::VectorGlobal2Local(data_buffer.mLocalCoordSystem, DeltVel, data_buffer.mLocalRelVel);
            mDiscontinuumConstitutiveLaw->CalculateForcesWithFEM(r_process_info,
                                                                OldLocalElasticContactForce,
                                                                LocalElasticContactForce,
                                                                LocalDeltDisp,
                                                                data_buffer.mLocalRelVel,
                                                                indentation,
                                                                previous_indentation,
                                                                ViscoDampingLocalContactForce,
                                                                cohesive_force,
                                                                this,
                                                                wall,
                                                                sliding);
        }

        double LocalContactForce[3]  = {0.0};
        double GlobalContactForce[3] = {0.0};

        AddUpFEMForcesAndProject(data_buffer.mLocalCoordSystem, LocalContactForce, LocalElasticContactForce, GlobalContactForce,
                                GlobalElasticContactForce, ViscoDampingLocalContactForce, cohesive_force, r_elastic_force,
                                r_contact_force,  neighbour_point_faces_elastic_contact_force[i], neighbour_point_faces_total_contact_force[i]);

        rigid_element_force[0] -= GlobalContactForce[0];
        rigid_element_force[1] -= GlobalContactForce[1];
        rigid_element_force[2] -= GlobalContactForce[2];

        if (this->Is(DEMFlags::HAS_ROTATION)) {
            ComputeMoments(LocalContactForce[2], GlobalContactForce, RollingResistance, data_buffer.mLocalCoordSystem[2], this, indentation, true); //WARNING: sending itself as the neighbor!!
        }

        if (this->Is(DEMFlags::HAS_STRESS_TENSOR)) {
            AddWallContributionToStressTensor(GlobalElasticContactForce, data_buffer.mLocalCoordSystem[2], DistPToB, 0.0);
        }
    }

    KRATOS_CATCH("")
}// ComputeBallToRigidFaceContactForce

void SphericParticle::RenewData()
{
  //To be redefined
}
void SphericParticle::ComputeConditionRelativeData(int rigid_neighbour_index,   // check for delete
                                                    DEMWall* const wall,
                                                    double LocalCoordSystem[3][3],
                                                    double& DistPToB,
                                                    array_1d<double, 4>& Weight,
                                                    array_1d<double, 3>& wall_delta_disp_at_contact_point,
                                                    array_1d<double, 3>& wall_velocity_at_contact_point,
                                                    int& ContactType)
{
    size_t FE_size = wall->GetGeometry().size();

    std::vector<double> TempWeight;
    TempWeight.resize(FE_size);

    double total_weight = 0.0;
    int points = 0;
    int inode1 = 0, inode2 = 0;

    for (unsigned int inode = 0; inode < FE_size; inode++) {

        if (Weight[inode] > 1.0e-12) {
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

    const double radius = this->GetInteractionRadius();

    if (points == 3 || points == 4)
    {
        unsigned int dummy_current_edge_index;
        contact_exists = GeometryFunctions::FacetCheck(wall->GetGeometry(), node_coordinates, radius, LocalCoordSystem, DistPToB, TempWeight, dummy_current_edge_index);
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
        contact_exists = GeometryFunctions::EdgeCheck(wall->GetGeometry()[inode1], wall->GetGeometry()[inode2], node_coordinates, radius, LocalCoordSystem, DistPToB, eta);

        Weight[inode1] = 1-eta;
        Weight[inode2] = eta;
        ContactType = 2;

    }

    else if (points == 1) {
        contact_exists = GeometryFunctions::VertexCheck(wall->GetGeometry()[inode1], node_coordinates, radius, LocalCoordSystem, DistPToB);
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
} //ComputeConditionRelativeData

void SphericParticle::ComputeWear(double LocalRelVel[3],
                                  double mTimeStep,
                                  bool sliding,
                                  double inverse_of_volume,
                                  double LocalElasticContactForce,
                                  DEMWall* wall) {

    array_1d<double, 3>& node_coor_array = this->GetGeometry()[0].Coordinates();
    double volume_wear = 0.0;
    double WallSeverityOfWear           = wall->GetProperties()[SEVERITY_OF_WEAR];
    double WallImpactSeverityOfWear     = wall->GetProperties()[IMPACT_WEAR_SEVERITY];
    double InverseOfWallBrinellHardness = 1.0 / (wall->GetProperties()[BRINELL_HARDNESS]);
    double Sliding_0 = LocalRelVel[0] * mTimeStep;
    double Sliding_1 = LocalRelVel[1] * mTimeStep;

    double impact_wear = WallImpactSeverityOfWear * InverseOfWallBrinellHardness * GetDensity() * mRadius * std::abs(LocalRelVel[2]);
    if (sliding) volume_wear = WallSeverityOfWear * InverseOfWallBrinellHardness * std::abs(LocalElasticContactForce) * sqrt(Sliding_0 * Sliding_0 + Sliding_1 * Sliding_1);

    double element_area = wall->GetGeometry().Area();

    if (element_area) {
        impact_wear /= element_area;
        volume_wear /= element_area;
    }

    //COMPUTING THE PROJECTED POINT
    array_1d<double, 3> inner_point = ZeroVector(3);
    array_1d<double, 3> relative_vector = wall->GetGeometry()[0].Coordinates() - node_coor_array; //We could have chosen [1] or [2], also.

    if (wall->GetGeometry().size()>2){
        array_1d<double, 3> normal_to_wall;

    wall->CalculateNormal(normal_to_wall);

    array_1d<double, 3> relative_vector = wall->GetGeometry()[0].Coordinates() - node_coor_array; //We could have chosen [1] or [2], also.

        double dot_prod = DEM_INNER_PRODUCT_3(relative_vector, normal_to_wall);

        DEM_MULTIPLY_BY_SCALAR_3(normal_to_wall, dot_prod);

        inner_point = node_coor_array + normal_to_wall;
    }
    else{
        // projection on a line element
        const double numerical_limit = std::numeric_limits<double>::epsilon();

        array_1d<double, 3> line_vector = wall->GetGeometry()[1].Coordinates()-wall->GetGeometry()[0].Coordinates();
        KRATOS_ERROR_IF(wall->GetGeometry().Length()<=numerical_limit) << "Line element has zero length" << std::endl;
        line_vector/=wall->GetGeometry().Length();

        DEM_COPY_SECOND_TO_FIRST_3(inner_point,line_vector);
        double dot_prod = DEM_INNER_PRODUCT_3(relative_vector, line_vector);
        DEM_MULTIPLY_BY_SCALAR_3(inner_point, dot_prod);
        DEM_ADD_SECOND_TO_FIRST(inner_point,wall->GetGeometry()[0].Coordinates());
    }



    array_1d<double, 3> point_local_coordinates;
    Vector shape_functions_coefs(3);
    wall->GetGeometry().PointLocalCoordinates(point_local_coordinates, inner_point);
    wall->GetGeometry().ShapeFunctionsValues(shape_functions_coefs, point_local_coordinates);

    if ((shape_functions_coefs[0] >= 0.0) && (shape_functions_coefs[1] >= 0.0) && (shape_functions_coefs[2] >= 0.0)) { // Which means the ball is inside this element
                                                                                                                       // This avoids negative shape functions that would give
                                                                                                                       // negative wear values
        wall->GetGeometry()[0].SetLock();
        wall->GetGeometry()[0].FastGetSolutionStepValue(NON_DIMENSIONAL_VOLUME_WEAR) += shape_functions_coefs[0] * volume_wear;
        wall->GetGeometry()[0].FastGetSolutionStepValue(IMPACT_WEAR) += shape_functions_coefs[0] * impact_wear;
        wall->GetGeometry()[0].UnSetLock();

        wall->GetGeometry()[1].SetLock();
        wall->GetGeometry()[1].FastGetSolutionStepValue(NON_DIMENSIONAL_VOLUME_WEAR) += shape_functions_coefs[1] * volume_wear;
        wall->GetGeometry()[1].FastGetSolutionStepValue(IMPACT_WEAR) += shape_functions_coefs[1] * impact_wear;
        wall->GetGeometry()[1].UnSetLock();

        wall->GetGeometry()[2].SetLock();
        wall->GetGeometry()[2].FastGetSolutionStepValue(NON_DIMENSIONAL_VOLUME_WEAR) += shape_functions_coefs[2] * volume_wear;
        wall->GetGeometry()[2].FastGetSolutionStepValue(IMPACT_WEAR) += shape_functions_coefs[2] * impact_wear;
        wall->GetGeometry()[2].UnSetLock();
    }
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

    array_1d<double, 3> normal_vector_on_contact;
    normal_vector_on_contact[0] = - other_to_me_vect[0]; //outwards
    normal_vector_on_contact[1] = - other_to_me_vect[1]; //outwards
    normal_vector_on_contact[2] = - other_to_me_vect[2]; //outwards

    array_1d<double, 3> x_centroid = real_distance * normal_vector_on_contact;

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
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

void SphericParticle::ComputeReactions() {

    KRATOS_TRY

    Node<3>& node = GetGeometry()[0];
    array_1d<double, 3>& reaction_force = node.FastGetSolutionStepValue(FORCE_REACTION);
    array_1d<double, 3>& r_total_forces = node.FastGetSolutionStepValue(TOTAL_FORCES);
    reaction_force[0] = node.Is(DEMFlags::FIXED_VEL_X) * (-r_total_forces[0]);
    reaction_force[1] = node.Is(DEMFlags::FIXED_VEL_Y) * (-r_total_forces[1]);
    reaction_force[2] = node.Is(DEMFlags::FIXED_VEL_Z) * (-r_total_forces[2]);

    if (this->Is(DEMFlags::HAS_ROTATION)) {
        array_1d<double, 3>& reaction_moment = this->GetGeometry()[0].FastGetSolutionStepValue(MOMENT_REACTION);
        array_1d<double, 3>& r_total_moment = this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT);
        reaction_moment[0] = node.Is(DEMFlags::FIXED_ANG_VEL_X) * (-r_total_moment[0]);
        reaction_moment[1] = node.Is(DEMFlags::FIXED_ANG_VEL_Y) * (-r_total_moment[1]);
        reaction_moment[2] = node.Is(DEMFlags::FIXED_ANG_VEL_Z) * (-r_total_moment[2]);
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

    if (this->Is(DEMFlags::CUMULATIVE_ZONE)) {
        const array_1d<double,3> gravity_force = ComputeWeight(gravity, r_process_info);
        const double gravity_force_magnitude = DEM_MODULUS_3(gravity_force);
        const array_1d<double, 3>& vel = this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
        const double vel_magnitude = DEM_MODULUS_3(vel);
        if (vel_magnitude != 0.0){
            const array_1d<double, 3> unitary_vel =  vel/vel_magnitude;
            const double inlet_damping_coefficient = 1e3;
            const array_1d<double, 3> damping_force = - inlet_damping_coefficient * GetMass() * vel_magnitude  * vel_magnitude * unitary_vel;
            const array_1d<double, 3> counter_force  = -5.0 * gravity_force_magnitude * unitary_vel;
            noalias(externally_applied_force)  += damping_force;
            noalias(externally_applied_force)  += counter_force;}
    } else {
        noalias(externally_applied_force)  += ComputeWeight(gravity, r_process_info);
        noalias(externally_applied_force)  += this->GetGeometry()[0].FastGetSolutionStepValue(EXTERNAL_APPLIED_FORCE);
        noalias(externally_applied_moment) += this->GetGeometry()[0].FastGetSolutionStepValue(EXTERNAL_APPLIED_MOMENT);
    }
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
                                               array_1d<double, 3>& elastic_force_backup,
                                               array_1d<double, 3>& total_force_backup) {
    for (unsigned int index = 0; index < 3; index++) {
        LocalContactForce[index] = LocalElasticContactForce[index] + ViscoDampingLocalContactForce[index];
    }
    LocalContactForce[2] -= cohesive_force;

    GeometryFunctions::VectorLocal2Global(LocalCoordSystem, LocalElasticContactForce, GlobalElasticContactForce);
    GeometryFunctions::VectorLocal2Global(LocalCoordSystem, LocalContactForce, GlobalContactForce);
    // Saving contact forces (We need to, since tangential elastic force is history-dependent)
    DEM_COPY_SECOND_TO_FIRST_3(elastic_force_backup,GlobalElasticContactForce)
    DEM_COPY_SECOND_TO_FIRST_3(total_force_backup,GlobalContactForce)
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

    mGlobalDamping = r_process_info[GLOBAL_DAMPING];
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

            //double K = Globals::Pi * GetYoung() * GetRadius(); //M. Error, should be the same that the local definition.

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

bool SphericParticle::CalculateRelativePositionsOrSkipContact(ParticleDataBuffer & data_buffer)
{
    const bool other_is_injecting_me = this->Is(NEW_ENTITY) && data_buffer.mpOtherParticle->Is(BLOCKED);
    const bool i_am_injecting_other = this->Is(BLOCKED) && data_buffer.mpOtherParticle->Is(NEW_ENTITY);
    const bool multistage_condition = data_buffer.mMultiStageRHS  &&  this->Id() > data_buffer.mpOtherParticle->Id();

    bool must_skip_contact_calculation = other_is_injecting_me || i_am_injecting_other || multistage_condition;

    if (must_skip_contact_calculation){
        return false;
    }

    NodeType& other_node = data_buffer.mpOtherParticle->GetGeometry()[0];
    DEM_COPY_SECOND_TO_FIRST_3(data_buffer.mOtherCoors, other_node)

    if (data_buffer.mDomainIsPeriodic){
        TransformNeighbourCoorsToClosestInPeriodicDomain(data_buffer);
    }

    data_buffer.mOtherToMeVector[0] = data_buffer.mMyCoors[0] - data_buffer.mOtherCoors[0];
    data_buffer.mOtherToMeVector[1] = data_buffer.mMyCoors[1] - data_buffer.mOtherCoors[1];
    data_buffer.mOtherToMeVector[2] = data_buffer.mMyCoors[2] - data_buffer.mOtherCoors[2];

    data_buffer.mDistance    = DEM_MODULUS_3(data_buffer.mOtherToMeVector);

    must_skip_contact_calculation = data_buffer.mDistance < std::numeric_limits<double>::epsilon();

    if (must_skip_contact_calculation){
        return false;
    }

    data_buffer.mOtherRadius = data_buffer.mpOtherParticle->GetInteractionRadius();
    data_buffer.mRadiusSum   = this->GetInteractionRadius() + data_buffer.mOtherRadius;
    data_buffer.mIndentation = data_buffer.mRadiusSum - data_buffer.mDistance;

    return data_buffer.mIndentation > 0.0;
}

void SphericParticle::RelativeDisplacementAndVelocityOfContactPointDueToOtherReasons(const ProcessInfo& r_process_info,
                                                                                    double DeltDisp[3], //IN GLOBAL AXES
                                                                                    double RelVel[3], //IN GLOBAL AXES
                                                                                    double OldLocalCoordSystem[3][3],
                                                                                    double LocalCoordSystem[3][3],
                                                                                    SphericParticle* neighbour_iterator) {}

void SphericParticle::SendForcesToFEM(){}

void SphericParticle::TransformNeighbourCoorsToClosestInPeriodicDomain(ParticleDataBuffer & data_buffer)
{
    const double periods[3] = {data_buffer.mDomainMax[0] - data_buffer.mDomainMin[0],
                               data_buffer.mDomainMax[1] - data_buffer.mDomainMin[1],
                               data_buffer.mDomainMax[2] - data_buffer.mDomainMin[2]};

    DiscreteParticleConfigure<3>::TransformToClosestPeriodicCoordinates(data_buffer.mMyCoors, data_buffer.mOtherCoors, periods);
}

void SphericParticle::TransformNeighbourCoorsToClosestInPeriodicDomain(ParticleDataBuffer & data_buffer,
                                                                       const array_1d<double, 3>& coors,
                                                                       array_1d<double, 3>& neighbour_coors)
{
    const double periods[3] = {data_buffer.mDomainMax[0] - data_buffer.mDomainMin[0],
                               data_buffer.mDomainMax[1] - data_buffer.mDomainMin[1],
                               data_buffer.mDomainMax[2] - data_buffer.mDomainMin[2]};

    DiscreteParticleConfigure<3>::TransformToClosestPeriodicCoordinates(coors, neighbour_coors, periods);
}

void SphericParticle::TransformNeighbourCoorsToClosestInPeriodicDomain(const ProcessInfo& r_process_info,
                                                                       const double coors[3],
                                                                       double neighbour_coors[3])
{
    const array_1d<double,3>& domain_min = r_process_info[DOMAIN_MIN_CORNER];
    const array_1d<double,3>& domain_max = r_process_info[DOMAIN_MAX_CORNER];

    const double periods[3] = {domain_min[0] - domain_max[0],
                               domain_min[1] - domain_max[1],
                               domain_min[2] - domain_max[2]};

    DiscreteParticleConfigure<3>::TransformToClosestPeriodicCoordinates(coors, neighbour_coors, periods);
}


void SphericParticle::Move(const double delta_t, const bool rotation_option, const double force_reduction_factor, const int StepFlag) {
    GetTranslationalIntegrationScheme().Move(GetGeometry()[0], delta_t, force_reduction_factor, StepFlag);
    if (rotation_option) {
        GetRotationalIntegrationScheme().Rotate(GetGeometry()[0], delta_t, force_reduction_factor, StepFlag);
    }
}

void SphericParticle::Calculate(const Variable<Vector >& rVariable, Vector& Output, const ProcessInfo& r_process_info){}
void SphericParticle::Calculate(const Variable<Matrix >& rVariable, Matrix& Output, const ProcessInfo& r_process_info){}

void SphericParticle::AdditionalCalculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& r_process_info){}

void SphericParticle::ApplyGlobalDampingToContactForcesAndMoments(array_1d<double,3>& total_forces, array_1d<double,3>& total_moment) {

        KRATOS_TRY

        const array_1d<double, 3> velocity =         this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
        const array_1d<double, 3> angular_velocity = this->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);

        if (this->GetGeometry()[0].IsNot(DEMFlags::FIXED_VEL_X)) {
            total_forces[0] *= (1.0 - mGlobalDamping * GeometryFunctions::sign(total_forces[0] * velocity[0]));
        }
        if (this->GetGeometry()[0].IsNot(DEMFlags::FIXED_VEL_Y)) {
            total_forces[1] *= (1.0 - mGlobalDamping * GeometryFunctions::sign(total_forces[1] * velocity[1]));
        }
        if (this->GetGeometry()[0].IsNot(DEMFlags::FIXED_VEL_Z)) {
            total_forces[2] *= (1.0 - mGlobalDamping * GeometryFunctions::sign(total_forces[2] * velocity[2]));
        }

        if (this->GetGeometry()[0].IsNot(DEMFlags::FIXED_ANG_VEL_X)) {
            total_moment[0] *= (1.0 - mGlobalDamping * GeometryFunctions::sign(total_moment[0] * angular_velocity[0]));
        }
        if (this->GetGeometry()[0].IsNot(DEMFlags::FIXED_ANG_VEL_Y)) {
            total_moment[1] *= (1.0 - mGlobalDamping * GeometryFunctions::sign(total_moment[1] * angular_velocity[1]));
        }
        if (this->GetGeometry()[0].IsNot(DEMFlags::FIXED_ANG_VEL_Z)) {
            total_moment[2] *= (1.0 - mGlobalDamping * GeometryFunctions::sign(total_moment[2] * angular_velocity[2]));
        }

        KRATOS_CATCH("")
    }

    void SphericParticle::CalculateOnContactElements(size_t i, double LocalContactForce[3]) {

        KRATOS_TRY

        if (!mBondElements.size()) return; // we skip this function if the vector of bonds hasn't been filled yet.
        ParticleContactElement* bond = mBondElements[i];
        if (bond == NULL) return; //This bond was never created (happens in some MPI cases, see CreateContactElements() in explicit_solve_continumm.h)

        bond->mLocalContactForce[0] = LocalContactForce[0];
        bond->mLocalContactForce[1] = LocalContactForce[1];
        bond->mLocalContactForce[2] = LocalContactForce[2];
        KRATOS_CATCH("")
    }

int    SphericParticle::GetClusterId()                                                           { return mClusterId;      }
void   SphericParticle::SetClusterId(int givenId)                                                { mClusterId = givenId;   }
double SphericParticle::GetRadius()                                                              { return mRadius;         }
double SphericParticle::CalculateVolume()                                                        { return 4.0 * Globals::Pi / 3.0 * mRadius * mRadius * mRadius;     }
void   SphericParticle::SetRadius(double radius)                                                 { mRadius = radius;       }
void   SphericParticle::SetRadius()                                                              { mRadius = GetGeometry()[0].FastGetSolutionStepValue(RADIUS);       }
double SphericParticle::GetInteractionRadius(const int radius_index)                             { return mRadius;         }
void   SphericParticle::SetInteractionRadius(const double radius, const int radius_index)        { mRadius = radius; GetGeometry()[0].FastGetSolutionStepValue(RADIUS) = radius;}
double SphericParticle::GetSearchRadius()                                                        { return mSearchRadius;   }
void   SphericParticle::SetSearchRadius(const double radius)                                     { mSearchRadius = radius; }
void SphericParticle::SetDefaultRadiiHierarchy(const double radius)
{
    SetRadius(radius);
    SetSearchRadius(radius);
}

double SphericParticle::GetMass()                                                                { return mRealMass;       }
void   SphericParticle::SetMass(double real_mass)                                                { mRealMass = real_mass;  GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) = real_mass;}
double SphericParticle::CalculateMomentOfInertia()                                               { return 0.4 * GetMass() * GetRadius() * GetRadius(); }

double SphericParticle::GetYoung()                                                               { return GetFastProperties()->GetYoung();                     }
double SphericParticle::GetRollingFriction()                                                     { return GetFastProperties()->GetRollingFriction();           }
double SphericParticle::GetRollingFrictionWithWalls()                                            { return GetFastProperties()->GetRollingFrictionWithWalls();  }
double SphericParticle::GetPoisson()                                                             { return GetFastProperties()->GetPoisson();                   }
double SphericParticle::GetTgOfFrictionAngle()                                                   { return GetFastProperties()->GetTgOfFrictionAngle() ;        }
double SphericParticle::GetCoefficientOfRestitution()                                            { return GetFastProperties()->GetCoefficientOfRestitution();  }
double SphericParticle::GetLnOfRestitCoeff()                                                     { return GetFastProperties()->GetLnOfRestitCoeff();           }
double SphericParticle::GetDensity()                                                             { return GetFastProperties()->GetDensity();                   }
int    SphericParticle::GetParticleMaterial()                                                    { return GetFastProperties()->GetParticleMaterial();          }
double SphericParticle::GetParticleCohesion()                                                    { return GetFastProperties()->GetParticleCohesion();          }
double SphericParticle::GetParticleKNormal()                                                     { return GetFastProperties()->GetParticleKNormal();           }
double SphericParticle::GetParticleKTangential()                                                 { return GetFastProperties()->GetParticleKTangential();       }
double SphericParticle::GetParticleContactRadius()                                               { return GetFastProperties()->GetParticleContactRadius();     }
double SphericParticle::GetParticleMaxStress()                                                   { return GetFastProperties()->GetParticleMaxStress();         }
double SphericParticle::GetParticleGamma()                                                       { return GetFastProperties()->GetParticleGamma();             }

array_1d<double, 3>& SphericParticle::GetForce()                                                 { return GetGeometry()[0].FastGetSolutionStepValue(TOTAL_FORCES);}
double&              SphericParticle::GetElasticEnergy()                                         { return mElasticEnergy; }
double&              SphericParticle::GetInelasticFrictionalEnergy()                             { return mInelasticFrictionalEnergy; }
double&              SphericParticle::GetInelasticViscodampingEnergy()                           { return mInelasticViscodampingEnergy; }

void   SphericParticle::SetYoungFromProperties(double* young)                                    { GetFastProperties()->SetYoungFromProperties( young);                             }
void   SphericParticle::SetRollingFrictionFromProperties(double* rolling_friction)               { GetFastProperties()->SetRollingFrictionFromProperties( rolling_friction);        }
void   SphericParticle::SetRollingFrictionWithWallsFromProperties(double* rolling_friction_with_walls) { GetFastProperties()->SetRollingFrictionWithWallsFromProperties( rolling_friction_with_walls); }
void   SphericParticle::SetPoissonFromProperties(double* poisson)                                { GetFastProperties()->SetPoissonFromProperties( poisson);                         }
void   SphericParticle::SetTgOfFrictionAngleFromProperties(double* tg_of_friction_angle)         { GetFastProperties()->SetTgOfFrictionAngleFromProperties( tg_of_friction_angle);  }
void   SphericParticle::SetCoefficientOfRestitutionFromProperties(double* coefficient_of_restitution) { GetFastProperties()->SetCoefficientOfRestitutionFromProperties( coefficient_of_restitution);      }
void   SphericParticle::SetLnOfRestitCoeffFromProperties(double* ln_of_restit_coeff)     { GetFastProperties()->SetLnOfRestitCoeffFromProperties( ln_of_restit_coeff);      }
void   SphericParticle::SetDensityFromProperties(double* density)                        { GetFastProperties()->SetDensityFromProperties( density);                         }
void   SphericParticle::SetParticleMaterialFromProperties(int* particle_material)        { GetFastProperties()->SetParticleMaterialFromProperties( particle_material);      }
void   SphericParticle::SetParticleCohesionFromProperties(double* particle_cohesion)     { GetFastProperties()->SetParticleCohesionFromProperties( particle_cohesion);      }
void   SphericParticle::SetParticleKNormalFromProperties(double* particle_k_normal)      { GetFastProperties()->SetParticleKNormalFromProperties( particle_k_normal);       }
void   SphericParticle::SetParticleKTangentialFromProperties(double* particle_k_tangential) { GetFastProperties()->SetParticleKTangentialFromProperties( particle_k_tangential); }
void   SphericParticle::SetParticleContactRadiusFromProperties(double* particle_contact_radius) { GetFastProperties()->SetParticleContactRadiusFromProperties( particle_contact_radius); }
void   SphericParticle::SetParticleMaxStressFromProperties(double* particle_max_stress)  { GetFastProperties()->SetParticleMaxStressFromProperties( particle_max_stress);   }
void   SphericParticle::SetParticleGammaFromProperties(double* particle_gamma)           { GetFastProperties()->SetParticleGammaFromProperties( particle_gamma);            }

DEMDiscontinuumConstitutiveLaw::Pointer SphericParticle::GetConstitutiveLawPointer(){return mDiscontinuumConstitutiveLaw;}


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
double SphericParticle::SlowGetRollingFriction()                                         { return GetProperties()[ROLLING_FRICTION];                                        }
double SphericParticle::SlowGetRollingFrictionWithWalls()                                { return GetProperties()[ROLLING_FRICTION_WITH_WALLS];                             }
double SphericParticle::SlowGetPoisson()                                                 { return GetProperties()[POISSON_RATIO];                                           }
double SphericParticle::SlowGetTgOfFrictionAngle()                                       { return GetProperties()[FRICTION];                                       }
double SphericParticle::SlowGetCoefficientOfRestitution()                                { return GetProperties()[COEFFICIENT_OF_RESTITUTION];                              }
double SphericParticle::SlowGetDensity()                                                 { return GetProperties()[PARTICLE_DENSITY];                                        }
int    SphericParticle::SlowGetParticleMaterial()                                        { return GetProperties()[PARTICLE_MATERIAL];                                       }
double SphericParticle::SlowGetParticleCohesion()                                        { return GetProperties()[PARTICLE_COHESION];                                       }

}  // namespace Kratos.
