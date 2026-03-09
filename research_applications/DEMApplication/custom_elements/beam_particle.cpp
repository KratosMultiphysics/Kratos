//
// Author: Joaquín Irazábal jirazabal@cimne.upc.edu
//

// System includes
#include <string>
#include <iostream>
#include <cmath>

#include <fstream>

// External includes

// Project includes
#include "beam_particle.h"

namespace Kratos {

    BeamParticle::BeamParticle() : SphericContinuumParticle() {}

    BeamParticle::BeamParticle(IndexType NewId, GeometryType::Pointer pGeometry) : SphericContinuumParticle(NewId, pGeometry) {}


    BeamParticle::BeamParticle(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : SphericContinuumParticle(NewId, pGeometry, pProperties) {}


    BeamParticle::BeamParticle(IndexType NewId, NodesArrayType const& ThisNodes)
    : SphericContinuumParticle(NewId, ThisNodes) {}


    BeamParticle::BeamParticle(Element::Pointer p_continuum_spheric_particle)
    {
        GeometryType::Pointer p_geom = p_continuum_spheric_particle->pGetGeometry();
        PropertiesType::Pointer pProperties = p_continuum_spheric_particle->pGetProperties();
        BeamParticle(p_continuum_spheric_particle->Id(), p_geom, pProperties);
    }

    Element::Pointer BeamParticle::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
    {
        GeometryType::Pointer p_geom = GetGeometry().Create(ThisNodes);

        return Element::Pointer(new BeamParticle(NewId, p_geom, pProperties));
    }

    void BeamParticle::CreateContinuumConstitutiveLaws() {

        mBeamConstitutiveLawArray.resize(mContinuumInitialNeighborsSize);

        for (unsigned int i = 0; i < mContinuumInitialNeighborsSize; i++) {
            Properties::Pointer properties_of_this_contact = GetProperties().pGetSubProperties(mNeighbourElements[i]->GetProperties().Id());
            mBeamConstitutiveLawArray[i] = (*properties_of_this_contact)[DEM_BEAM_CONSTITUTIVE_LAW_POINTER]-> Clone();
            SphericContinuumParticle* p_cont_neighbour_particle = dynamic_cast<SphericContinuumParticle*>(mNeighbourElements[i]);
            mBeamConstitutiveLawArray[i]->Initialize(this, p_cont_neighbour_particle, properties_of_this_contact);
        }
    }

    void BeamParticle::ContactAreaWeighting() {}

    void BeamParticle::CalculateMeanContactArea(const bool has_mpi, const ProcessInfo& r_process_info) {}

    double BeamParticle::CalculateMaxSearchDistance(const bool has_mpi, const ProcessInfo& r_process_info) { return 0.0; }

    void BeamParticle::Initialize(const ProcessInfo& r_process_info)
    {

        KRATOS_TRY

        SphericContinuumParticle::Initialize(r_process_info);

        double distance = GetProperties()[BEAM_PARTICLES_DISTANCE];
        auto& n0 = GetGeometry()[0];

        if (distance)
        {
            double contact_area = GetProperties()[CROSS_AREA];

            if (IsSkin()) distance *= 0.5;

            n0.FastGetSolutionStepValue(REPRESENTATIVE_VOLUME) = distance * contact_area;
            SetMass(GetDensity() * distance * contact_area);

            if (this->Is(DEMFlags::HAS_ROTATION)) {
                const double a = std::sqrt(12.0 * GetProperties()[BEAM_INERTIA_ROT_UNIT_LENGHT_Y] - 1.0);
                const double b = std::sqrt(12.0 * GetProperties()[BEAM_INERTIA_ROT_UNIT_LENGHT_Z] - 1.0);

                n0.FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA)[0] = GetProperties()[BEAM_INERTIA_ROT_UNIT_LENGHT_X] * GetDensity() * distance * contact_area;
                n0.FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA)[1] = 0.083333333 * (a * a + distance * distance) * GetDensity() * distance * contact_area;
                n0.FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA)[2] = 0.083333333 * (b * b + distance * distance) * GetDensity() * distance * contact_area;
            }
        }
        else
        {
            if (this->Is(DEMFlags::HAS_ROTATION)) {
                n0.FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA)[0] = n0.FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA);
                n0.FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA)[1] = n0.FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA);
                n0.FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA)[2] = n0.FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA);
            }
        }

        array_1d<double, 3> base_principal_moments_of_inertia = n0.FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA);

        Quaternion<double>& Orientation = n0.FastGetSolutionStepValue(ORIENTATION);
        Orientation.normalize();

        array_1d<double, 3> angular_velocity = n0.FastGetSolutionStepValue(ANGULAR_VELOCITY);
        array_1d<double, 3> angular_momentum;
        double LocalTensor[3][3];
        double GlobalTensor[3][3];
        GeometryFunctions::ConstructLocalTensor(base_principal_moments_of_inertia, LocalTensor);
        GeometryFunctions::QuaternionTensorLocal2Global(Orientation, LocalTensor, GlobalTensor);
        GeometryFunctions::ProductMatrix3X3Vector3X1(GlobalTensor, angular_velocity, angular_momentum);
        noalias(n0.FastGetSolutionStepValue(ANGULAR_MOMENTUM)) = angular_momentum;

        array_1d<double, 3> local_angular_velocity;
        GeometryFunctions::QuaternionVectorGlobal2Local(Orientation, angular_velocity, local_angular_velocity);
        noalias(n0.FastGetSolutionStepValue(LOCAL_ANGULAR_VELOCITY)) = local_angular_velocity;

        KRATOS_CATCH("")
    }

    void BeamParticle::InitializeSolutionStep(const ProcessInfo& r_process_info)
    {
        KRATOS_TRY

        mRadius = this->GetGeometry()[0].FastGetSolutionStepValue(RADIUS); //Just in case someone is overwriting the radius in Python
        mPartialRepresentativeVolume = 0.0;
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

    void BeamParticle::ComputeBallToBallContactForceAndMoment(SphericParticle::ParticleDataBuffer & data_buffer,
                                                     const ProcessInfo& r_process_info,
                                                     array_1d<double, 3>& rElasticForce,
                                                     array_1d<double, 3>& rContactForce)
    {

        KRATOS_TRY

        NodeType& this_node = this->GetGeometry()[0];
        DEM_COPY_SECOND_TO_FIRST_3(data_buffer.mMyCoors, this_node)

        const int time_steps = r_process_info[TIME_STEPS];

        const array_1d<double, 3>& vel         = this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
        const array_1d<double, 3>& delta_displ = this->GetGeometry()[0].FastGetSolutionStepValue(DELTA_DISPLACEMENT);
        const array_1d<double, 3>& ang_vel     = this->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);

        int NeighbourSize = mNeighbourElements.size();
        GetGeometry()[0].GetSolutionStepValue(NEIGHBOUR_SIZE) = NeighbourSize;

        for (int i = 0; data_buffer.SetNextNeighbourOrExit(i); ++i) {

            if (mNeighbourElements[i] == NULL) continue;
            if (this->Is(NEW_ENTITY) && mNeighbourElements[i]->Is(NEW_ENTITY)) continue;
            SphericContinuumParticle* neighbour_iterator = dynamic_cast<SphericContinuumParticle*>(mNeighbourElements[i]);
            data_buffer.mpOtherParticle = neighbour_iterator;

            unsigned int neighbour_iterator_id = data_buffer.mpOtherParticle->Id();

            noalias(data_buffer.mOtherToMeVector) = this->GetGeometry()[0].Coordinates() - data_buffer.mpOtherParticle->GetGeometry()[0].Coordinates();

            const double& other_radius = data_buffer.mpOtherParticle->GetRadius();

            data_buffer.mDistance = DEM_MODULUS_3(data_buffer.mOtherToMeVector);

            double radius_sum = GetRadius() + other_radius;
            double initial_delta = GetInitialDelta(i);
            double initial_dist = radius_sum - initial_delta;
            double indentation = initial_dist - data_buffer.mDistance;

            double myYoung = GetYoung();
            double myPoisson = GetPoisson();

            double kn_el = 0.0;
            double kt_el_0 = 0.0;
            double kt_el_1 = 0.0;
            double DeltDisp[3] = {0.0};
            double RelVel[3] = {0.0};
            DEM_SET_COMPONENTS_TO_ZERO_3x3(data_buffer.mLocalCoordSystem)
            DEM_SET_COMPONENTS_TO_ZERO_3x3(data_buffer.mOldLocalCoordSystem)
            bool sliding = false;

            double contact_tau = 0.0;
            double contact_sigma = 0.0;
            double failure_criterion_state = 0.0;
            double acumulated_damage = 0.0;

            // Getting neighbor properties
            double other_young = data_buffer.mpOtherParticle->GetYoung();
            double other_poisson = data_buffer.mpOtherParticle->GetPoisson();
            double equiv_poisson;
            if ((myPoisson + other_poisson) != 0.0) { equiv_poisson = 2.0 * myPoisson * other_poisson / (myPoisson + other_poisson); }
            else { equiv_poisson = 0.0; }

            double equiv_young = 2.0 * myYoung * other_young / (myYoung + other_young);
            double calculation_area = 0.0;
            const double equiv_shear = equiv_young / (2.0 * (1 + equiv_poisson));

            if (i < (int)mContinuumInitialNeighborsSize) {
                double area = this->GetProperties()[CROSS_AREA];
                double other_area = data_buffer.mpOtherParticle->GetProperties()[CROSS_AREA];
                calculation_area = std::max(area, other_area);
                mBeamConstitutiveLawArray[i]->CalculateElasticConstants(kn_el,
                                                                        kt_el_0,
                                                                        kt_el_1,
                                                                        initial_dist,
                                                                        equiv_young,
                                                                        equiv_poisson,
                                                                        calculation_area,
                                                                        this,
                                                                        neighbour_iterator,
                                                                        indentation);
            }

            EvaluateDeltaDisplacement(data_buffer,
                                      DeltDisp,
                                      RelVel,
                                      data_buffer.mLocalCoordSystem,
                                      data_buffer.mOldLocalCoordSystem,
                                      vel,
                                      delta_displ);

            if (this->Is(DEMFlags::HAS_ROTATION)) {
                RelativeDisplacementAndVelocityOfContactPointDueToRotationQuaternion(DeltDisp,
                                                                                     RelVel,
                                                                                     data_buffer.mOldLocalCoordSystem,
                                                                                     other_radius,
                                                                                     data_buffer.mDt,
                                                                                     ang_vel,
                                                                                     neighbour_iterator,
                                                                                     data_buffer);
            }

            RelativeDisplacementAndVelocityOfContactPointDueToOtherReasons(r_process_info,
                                                                           DeltDisp,
                                                                           RelVel,
                                                                           data_buffer.mOldLocalCoordSystem,
                                                                           data_buffer.mLocalCoordSystem,
                                                                           neighbour_iterator);

            double LocalDeltDisp[3] = {0.0};
            double LocalElasticContactForce[3] = {0.0}; // 0: first tangential, // 1: second tangential, // 2: normal force
            double LocalElasticExtraContactForce[3] = {0.0};
            double GlobalElasticContactForce[3] = {0.0};
            double GlobalElasticExtraContactForce[3] = {0.0};
            double TotalGlobalElasticContactForce[3] = {0.0};
            double OldLocalElasticContactForce[3] = {0.0};

            double cohesive_force =  0.0;

            GeometryFunctions::VectorGlobal2Local(data_buffer.mLocalCoordSystem, DeltDisp, LocalDeltDisp);

            RotateOldContactForces(data_buffer.mOldLocalCoordSystem, data_buffer.mLocalCoordSystem, mNeighbourElasticContactForces[i]);
            RotateOldContactForces(data_buffer.mOldLocalCoordSystem, data_buffer.mLocalCoordSystem, mNeighbourElasticExtraContactForces[i]);

            GeometryFunctions::VectorGlobal2Local(data_buffer.mLocalCoordSystem, mNeighbourElasticContactForces[i], OldLocalElasticContactForce);

            GlobalElasticContactForce[0] = mNeighbourElasticContactForces[i][0];
            GlobalElasticContactForce[1] = mNeighbourElasticContactForces[i][1];
            GlobalElasticContactForce[2] = mNeighbourElasticContactForces[i][2];

            GeometryFunctions::VectorGlobal2Local(data_buffer.mLocalCoordSystem, GlobalElasticContactForce, LocalElasticContactForce); //TODO: can we remove this? We should overwrite LocalElasticContactForce afterwards

            double ViscoDampingLocalContactForce[3] = {0.0};

            double LocalRelVel[3] = {0.0};
            GeometryFunctions::VectorGlobal2Local(data_buffer.mLocalCoordSystem, RelVel, LocalRelVel);

            if (i < (int)mContinuumInitialNeighborsSize) {

                double equiv_visco_damp_coeff_normal;
                double equiv_visco_damp_coeff_tangential_0;
                double equiv_visco_damp_coeff_tangential_1;

                mBeamConstitutiveLawArray[i]->CalculateForces(r_process_info,
                                                              OldLocalElasticContactForce,
                                                              LocalElasticContactForce,
                                                              LocalElasticExtraContactForce,
                                                              data_buffer.mLocalCoordSystem,
                                                              LocalDeltDisp,
                                                              kn_el,
                                                              kt_el_0,
                                                              kt_el_1,
                                                              contact_sigma,
                                                              contact_tau,
                                                              failure_criterion_state,
                                                              equiv_young,
                                                              equiv_shear,
                                                              indentation,
                                                              calculation_area,
                                                              acumulated_damage,
                                                              this,
                                                              neighbour_iterator,
                                                              i,
                                                              r_process_info[TIME_STEPS],
                                                              sliding,
                                                              equiv_visco_damp_coeff_normal,
                                                              equiv_visco_damp_coeff_tangential_0,
                                                              equiv_visco_damp_coeff_tangential_1,
                                                              LocalRelVel,
                                                              ViscoDampingLocalContactForce);

            } else if (indentation > 0.0) {
                const double previous_indentation = indentation + LocalDeltDisp[2];
                mDiscontinuumConstitutiveLaw = pCloneDiscontinuumConstitutiveLawWithNeighbour(data_buffer.mpOtherParticle);
                mDiscontinuumConstitutiveLaw->CalculateForces(r_process_info,
                                                              OldLocalElasticContactForce,
                                                              LocalElasticContactForce,
                                                              LocalDeltDisp,
                                                              LocalRelVel,
                                                              indentation,
                                                              previous_indentation,
                                                              ViscoDampingLocalContactForce,
                                                              cohesive_force,
                                                              this,
                                                              data_buffer.mpOtherParticle,
                                                              sliding,
                                                              data_buffer.mLocalCoordSystem);

            } else { //Not bonded and no idata_buffer.mpOtherParticlendentation
                LocalElasticContactForce[0] = 0.0; LocalElasticContactForce[1] = 0.0; LocalElasticContactForce[2] = 0.0;
                ViscoDampingLocalContactForce[0] = 0.0; ViscoDampingLocalContactForce[1] = 0.0; ViscoDampingLocalContactForce[2] = 0.0;
                cohesive_force= 0.0;
            }

            // Transforming to global forces and adding up
            double LocalContactForce[3] = {0.0};
            double GlobalContactForce[3] = {0.0};

            array_1d<double, 3> other_ball_to_ball_forces(3,0.0);
            ComputeOtherBallToBallForces(other_ball_to_ball_forces);

            AddUpForcesAndProject(data_buffer.mOldLocalCoordSystem,
                                  data_buffer.mLocalCoordSystem,
                                  LocalContactForce,
                                  LocalElasticContactForce,
                                  LocalElasticExtraContactForce,
                                  GlobalContactForce,
                                  GlobalElasticContactForce,
                                  GlobalElasticExtraContactForce,
                                  TotalGlobalElasticContactForce,
                                  ViscoDampingLocalContactForce,
                                  0.0,
                                  other_ball_to_ball_forces,
                                  rElasticForce,
                                  rContactForce,
                                  i,
                                  r_process_info); // 0.0 means null cohesive force

            if (this->Is(DEMFlags::HAS_ROTATION)) {

                double ElasticLocalRotationalMoment[3] = {0.0};
                double ViscoLocalRotationalMoment[3] = {0.0};

                ComputeMoments(LocalContactForce[2],
                               GlobalContactForce,
                               data_buffer.mLocalCoordSystem[2],
                               data_buffer.mpOtherParticle,
                               indentation,
                               i);

                if (i < (int)mContinuumInitialNeighborsSize) {

                    mBeamConstitutiveLawArray[i]->ComputeParticleRotationalMoments(this,
                                                                                   neighbour_iterator,
                                                                                   equiv_young,
                                                                                   data_buffer.mDistance,
                                                                                   calculation_area,
                                                                                   data_buffer.mLocalCoordSystem,
                                                                                   ElasticLocalRotationalMoment,
                                                                                   ViscoLocalRotationalMoment,
                                                                                   equiv_poisson,
                                                                                   indentation);
                }

                AddUpMomentsAndProject(data_buffer.mLocalCoordSystem, ElasticLocalRotationalMoment, ViscoLocalRotationalMoment);
            }

            if (r_process_info[CONTACT_MESH_OPTION] == 1 && (i < (int)mContinuumInitialNeighborsSize) && this->Id() < neighbour_iterator_id) {
                double total_local_elastic_contact_force[3] = {0.0};
                total_local_elastic_contact_force[0] = LocalElasticContactForce[0] + LocalElasticExtraContactForce[0];
                total_local_elastic_contact_force[1] = LocalElasticContactForce[1] + LocalElasticExtraContactForce[1];
                total_local_elastic_contact_force[2] = LocalElasticContactForce[2] + LocalElasticExtraContactForce[2];
                SphericContinuumParticle::CalculateOnContinuumContactElements(i,
                                                                              total_local_elastic_contact_force,
                                                                              contact_sigma,
                                                                              contact_tau,
                                                                              failure_criterion_state,
                                                                              acumulated_damage,
                                                                              time_steps,
                                                                              calculation_area,
                                                                              GlobalContactForce);
            }
        } // for each neighbor

        KRATOS_CATCH("")
    } //  ComputeBallToBallContactForceAndMoment

    void BeamParticle::Move(const double delta_t, const bool rotation_option, const double force_reduction_factor, const int StepFlag) {

        KRATOS_TRY

        GetTranslationalIntegrationScheme().Move(GetGeometry()[0], delta_t, force_reduction_factor, StepFlag);
        if (rotation_option) {
            GetRotationalIntegrationScheme().CalculateRotationalMotionOfRigidBodyElementNode(GetGeometry()[0], delta_t, force_reduction_factor, StepFlag);
        }

        KRATOS_CATCH("")
    }

    void BeamParticle::FinalizeSolutionStep(const ProcessInfo& r_process_info) {}

} // namespace Kratos
