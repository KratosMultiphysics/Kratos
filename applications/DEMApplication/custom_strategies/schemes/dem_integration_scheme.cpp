//
// Author: Miguel AngelCeligueta, maceli@cimne.upc.edu
//

#include "custom_utilities/GeometryFunctions.h"
#include "custom_elements/cluster3D.h"
#include "custom_elements/rigid_body_element.h"
#include "dem_integration_scheme.h"
#include "DEM_application_variables.h"


namespace Kratos {

    DEMIntegrationScheme::DEMIntegrationScheme(){}
    DEMIntegrationScheme::~DEMIntegrationScheme(){}

    void DEMIntegrationScheme::SetTranslationalIntegrationSchemeInProperties(Properties::Pointer pProp, bool verbose) const {
        //if(verbose) KRATOS_INFO("DEM") << "Assigning DEMDiscontinuumConstitutiveLaw to properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_TRANSLATIONAL_INTEGRATION_SCHEME_POINTER, this->CloneShared());
    }

    void DEMIntegrationScheme::SetRotationalIntegrationSchemeInProperties(Properties::Pointer pProp, bool verbose) const {
        //if(verbose) KRATOS_INFO("DEM") << "Assigning DEMDiscontinuumConstitutiveLaw to properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_ROTATIONAL_INTEGRATION_SCHEME_POINTER, this->CloneShared());
    }

    void DEMIntegrationScheme::Move(Node<3> & i, const double delta_t, const double force_reduction_factor, const int StepFlag) {
        if (i.Is(DEMFlags::BELONGS_TO_A_CLUSTER)) return;
        CalculateTranslationalMotionOfNode(i, delta_t, force_reduction_factor, StepFlag);
    }

    void DEMIntegrationScheme::Rotate(Node<3> & i, const double delta_t, const double moment_reduction_factor, const int StepFlag) {
        if (i.Is(DEMFlags::BELONGS_TO_A_CLUSTER)) return;
        CalculateRotationalMotionOfSphereNode(i, delta_t, moment_reduction_factor, StepFlag);
    }

    void DEMIntegrationScheme::MoveRigidBodyElement(RigidBodyElement3D* rigid_body_element, Node<3> & i, const double delta_t, const double force_reduction_factor, const int StepFlag) {
        CalculateTranslationalMotionOfNode(i, delta_t, force_reduction_factor, StepFlag);
        rigid_body_element->UpdateLinearDisplacementAndVelocityOfNodes();
    }

    void DEMIntegrationScheme::RotateRigidBodyElement(RigidBodyElement3D* rigid_body_element, Node<3> & i, const double delta_t, const double moment_reduction_factor, const int StepFlag) {
        CalculateRotationalMotionOfRigidBodyElementNode(i, delta_t, moment_reduction_factor, StepFlag);
        rigid_body_element->UpdateAngularDisplacementAndVelocityOfNodes();
    }

    void DEMIntegrationScheme::CalculateTranslationalMotionOfNode(Node<3> & i, const double delta_t, const double force_reduction_factor, const int StepFlag) {
        array_1d<double, 3 >& vel = i.FastGetSolutionStepValue(VELOCITY);
        array_1d<double, 3 >& displ = i.FastGetSolutionStepValue(DISPLACEMENT);
        array_1d<double, 3 >& delta_displ = i.FastGetSolutionStepValue(DELTA_DISPLACEMENT);
        array_1d<double, 3 >& coor = i.Coordinates();
        array_1d<double, 3 >& initial_coor = i.GetInitialPosition();
        array_1d<double, 3 >& force = i.FastGetSolutionStepValue(TOTAL_FORCES);

        #ifdef KRATOS_DEBUG
            DemDebugFunctions::CheckIfNan(force, "NAN in Force in Integration Scheme");
        #endif

        double mass = i.FastGetSolutionStepValue(NODAL_MASS);

        bool Fix_vel[3] = {false, false, false};

        Fix_vel[0] = i.Is(DEMFlags::FIXED_VEL_X);
        Fix_vel[1] = i.Is(DEMFlags::FIXED_VEL_Y);
        Fix_vel[2] = i.Is(DEMFlags::FIXED_VEL_Z);

        UpdateTranslationalVariables(StepFlag, i, coor, displ, delta_displ, vel, initial_coor, force, force_reduction_factor, mass, delta_t, Fix_vel);
    }

    void DEMIntegrationScheme::CalculateRotationalMotionOfSphereNode(Node<3> & i, const double delta_t, const double moment_reduction_factor, const int StepFlag) {

        double moment_of_inertia               = i.FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA);
        array_1d<double, 3 >& angular_velocity = i.FastGetSolutionStepValue(ANGULAR_VELOCITY);
        array_1d<double, 3 >& torque           = i.FastGetSolutionStepValue(PARTICLE_MOMENT);
        array_1d<double, 3 >& rotated_angle    = i.FastGetSolutionStepValue(PARTICLE_ROTATION_ANGLE);
        array_1d<double, 3 >& delta_rotation   = i.FastGetSolutionStepValue(DELTA_ROTATION);

        #ifdef KRATOS_DEBUG
        DemDebugFunctions::CheckIfNan(torque, "NAN in Torque in Integration Scheme");
        #endif

        bool Fix_Ang_vel[3] = {false, false, false};
        Fix_Ang_vel[0] = i.Is(DEMFlags::FIXED_ANG_VEL_X);
        Fix_Ang_vel[1] = i.Is(DEMFlags::FIXED_ANG_VEL_Y);
        Fix_Ang_vel[2] = i.Is(DEMFlags::FIXED_ANG_VEL_Z);

        CalculateNewRotationalVariablesOfSpheres(StepFlag, i, moment_of_inertia, angular_velocity, torque, moment_reduction_factor, rotated_angle, delta_rotation, delta_t, Fix_Ang_vel);
    }

    void DEMIntegrationScheme::CalculateRotationalMotionOfRigidBodyElementNode(Node<3> & i, const double delta_t, const double moment_reduction_factor, const int StepFlag) {

        array_1d<double, 3 >& moments_of_inertia = i.FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA);
        array_1d<double, 3 >& angular_velocity   = i.FastGetSolutionStepValue(ANGULAR_VELOCITY);
        array_1d<double, 3 >& torque             = i.FastGetSolutionStepValue(PARTICLE_MOMENT);
        array_1d<double, 3 >& rotated_angle      = i.FastGetSolutionStepValue(PARTICLE_ROTATION_ANGLE);
        array_1d<double, 3 >& delta_rotation     = i.FastGetSolutionStepValue(DELTA_ROTATION);
        Quaternion<double  >& Orientation        = i.FastGetSolutionStepValue(ORIENTATION);

        #ifdef KRATOS_DEBUG
        DemDebugFunctions::CheckIfNan(torque, "NAN in Torque in Integration Scheme");
        #endif

        bool Fix_Ang_vel[3] = {false, false, false};

        Fix_Ang_vel[0] = i.Is(DEMFlags::FIXED_ANG_VEL_X);
        Fix_Ang_vel[1] = i.Is(DEMFlags::FIXED_ANG_VEL_Y);
        Fix_Ang_vel[2] = i.Is(DEMFlags::FIXED_ANG_VEL_Z);

        CalculateNewRotationalVariablesOfRigidBodyElements(StepFlag, i, moments_of_inertia, angular_velocity, torque, moment_reduction_factor, rotated_angle, delta_rotation, Orientation, delta_t, Fix_Ang_vel);
    }

    void DEMIntegrationScheme::UpdateTranslationalVariables(
                int StepFlag,
                Node < 3 >& i,
                array_1d<double, 3 >& coor,
                array_1d<double, 3 >& displ,
                array_1d<double, 3 >& delta_displ,
                array_1d<double, 3 >& vel,
                const array_1d<double, 3 >& initial_coor,
                const array_1d<double, 3 >& force,
                const double force_reduction_factor,
                const double mass,
                const double delta_t,
                const bool Fix_vel[3])
    {
        KRATOS_THROW_ERROR(std::runtime_error, "This function (DEMIntegrationScheme::UpdateTranslationalVariables) shouldn't be accessed, use derived class instead", 0);
    }

    void DEMIntegrationScheme::CalculateNewRotationalVariablesOfSpheres(
                int StepFlag,
                Node < 3 >& i,
                const double moment_of_inertia,
                array_1d<double, 3 >& angular_velocity,
                array_1d<double, 3 >& torque,
                const double moment_reduction_factor,
                array_1d<double, 3 >& rotated_angle,
                array_1d<double, 3 >& delta_rotation,
                const double delta_t,
                const bool Fix_Ang_vel[3]) {
        KRATOS_THROW_ERROR(std::runtime_error, "This function (DEMIntegrationScheme::CalculateNewRotationalVariablesOfSpheres) shouldn't be accessed, use derived class instead", 0);
    }

    void DEMIntegrationScheme::CalculateNewRotationalVariablesOfRigidBodyElements(
                int StepFlag,
                Node < 3 >& i,
                const array_1d<double, 3 > moments_of_inertia,
                array_1d<double, 3 >& angular_velocity,
                array_1d<double, 3 >& torque,
                const double moment_reduction_factor,
                array_1d<double, 3 >& rotated_angle,
                array_1d<double, 3 >& delta_rotation,
                Quaternion<double  >& Orientation,
                const double delta_t,
                const bool Fix_Ang_vel[3]) {
        KRATOS_THROW_ERROR(std::runtime_error, "This function (DEMIntegrationScheme::CalculateNewRotationalVariablesOfRigidBodyElements) shouldn't be accessed, use derived class instead", 0);
    }

    void DEMIntegrationScheme::UpdateRotationalVariables(
                int StepFlag,
                Node < 3 >& i,
                array_1d<double, 3 >& rotated_angle,
                array_1d<double, 3 >& delta_rotation,
                array_1d<double, 3 >& angular_velocity,
                array_1d<double, 3 >& angular_acceleration,
                const double delta_t,
                const bool Fix_Ang_vel[3]) {
        KRATOS_THROW_ERROR(std::runtime_error, "This function (DEMIntegrationScheme::UpdateRotationalVariables) shouldn't be accessed, use derived class instead", 0);
    }

    void DEMIntegrationScheme::UpdateRotationalVariables(
                int StepFlag,
                Node < 3 >& i,
                const double& moment_of_inertia,
                array_1d<double, 3 >& rotated_angle,
                array_1d<double, 3 >& delta_rotation,
                Quaternion<double  >& Orientation,
                const array_1d<double, 3 >& angular_momentum,
                array_1d<double, 3 >& angular_velocity,
                const double delta_t,
                const bool Fix_Ang_vel[3]) {
        KRATOS_THROW_ERROR(std::runtime_error, "This function (DEMIntegrationScheme::UpdateRotationalVariables) shouldn't be accessed, use derived class instead", 0);
    }

    void DEMIntegrationScheme::UpdateRotationalVariables(
                int StepFlag,
                Node < 3 >& i,
                const array_1d<double, 3 >& moments_of_inertia,
                array_1d<double, 3 >& rotated_angle,
                array_1d<double, 3 >& delta_rotation,
                Quaternion<double  >& Orientation,
                const array_1d<double, 3 >& angular_momentum,
                array_1d<double, 3 >& angular_velocity,
                const double delta_t,
                const bool Fix_Ang_vel[3]) {
        KRATOS_THROW_ERROR(std::runtime_error, "This function (DEMIntegrationScheme::UpdateRotationalVariables) shouldn't be accessed, use derived class instead", 0);
    }

    void DEMIntegrationScheme::UpdateRotatedAngle(
                array_1d<double, 3 >& rotated_angle,
                array_1d<double, 3 >& delta_rotation,
                const array_1d<double, 3 >& angular_velocity,
                const double delta_t) {
        KRATOS_THROW_ERROR(std::runtime_error, "This function (DEMIntegrationScheme::UpdateRotatedAngle) shouldn't be accessed, use derived class instead", 0);
    }

    void DEMIntegrationScheme::UpdateAngularVelocity(
                const Quaternion<double>& Orientation,
                const double LocalTensorInv[3][3],
                const array_1d<double, 3>& angular_momentum,
                array_1d<double, 3>& angular_velocity) {
        KRATOS_THROW_ERROR(std::runtime_error, "This function (DEMIntegrationScheme::UpdateAngularVelocity) shouldn't be accessed, use derived class instead", 0);
    }

    void DEMIntegrationScheme::CalculateLocalAngularAcceleration(
                const double moment_of_inertia,
                const array_1d<double, 3 >& torque,
                const double moment_reduction_factor,
                array_1d<double, 3 >& angular_acceleration) {
        KRATOS_THROW_ERROR(std::runtime_error, "This function (DEMIntegrationScheme::CalculateLocalAngularAcceleration) shouldn't be accessed, use derived class instead", 0);
    }

    void DEMIntegrationScheme::CalculateLocalAngularAccelerationByEulerEquations(
                const array_1d<double, 3 >& local_angular_velocity,
                const array_1d<double, 3 >& moments_of_inertia,
                const array_1d<double, 3 >& local_torque,
                const double moment_reduction_factor,
                array_1d<double, 3 >& local_angular_acceleration) {
            KRATOS_THROW_ERROR(std::runtime_error, "This function (DEMIntegrationScheme::CalculateLocalAngularAccelerationByEulerEquations) shouldn't be accessed, use derived class instead", 0);
    }

    void DEMIntegrationScheme::CalculateAngularVelocityRK(
                const Quaternion<double  >& Orientation,
                const double& moment_of_inertia,
                const array_1d<double, 3 >& angular_momentum,
                array_1d<double, 3 >& angular_velocity,
                const double delta_t,
                const bool Fix_Ang_vel[3]) {
            KRATOS_THROW_ERROR(std::runtime_error, "This function (DEMIntegrationScheme::CalculateAngularVelocityRK) shouldn't be accessed, use derived class instead", 0);
    }

    void DEMIntegrationScheme::CalculateAngularVelocityRK(
                const Quaternion<double  >& Orientation,
                const array_1d<double, 3 >& moments_of_inertia,
                const array_1d<double, 3 >& angular_momentum,
                array_1d<double, 3 >& angular_velocity,
                const double delta_t,
                const bool Fix_Ang_vel[3]) {
            KRATOS_THROW_ERROR(std::runtime_error, "This function (DEMIntegrationScheme::CalculateAngularVelocityRK) shouldn't be accessed, use derived class instead", 0);
    }

    void DEMIntegrationScheme::QuaternionCalculateMidAngularVelocities(
                const Quaternion<double>& Orientation,
                const double LocalTensorInv[3][3],
                const array_1d<double, 3>& angular_momentum,
                const double dt,
                const array_1d<double, 3>& InitialAngularVel,
                array_1d<double, 3>& FinalAngularVel) {
        KRATOS_THROW_ERROR(std::runtime_error, "This function (DEMIntegrationScheme::QuaternionCalculateMidAngularVelocities) shouldn't be accessed, use derived class instead", 0);
    }
}
