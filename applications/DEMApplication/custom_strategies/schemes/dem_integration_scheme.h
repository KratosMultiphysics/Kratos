//
// Author: Miguel AngelCeligueta, maceli@cimne.upc.edu
//

#if !defined KRATOS_DEM_INTEGRATION_SCHEME_H_INCLUDED
#define KRATOS_DEM_INTEGRATION_SCHEME_H_INCLUDED

// Project includes
#include "includes/define.h"
#include "utilities/openmp_utils.h"
#include "includes/model_part.h"
#include "utilities/quaternion.h"

// System includes
#include <float.h>
#include <string>
#include <iostream>

namespace Kratos {

    class Cluster3D;
    class RigidBodyElement3D;

    class KRATOS_API(DEM_APPLICATION) DEMIntegrationScheme {
    public:

        typedef ModelPart::NodesContainerType NodesArrayType;
        KRATOS_CLASS_POINTER_DEFINITION(DEMIntegrationScheme);

        DEMIntegrationScheme();

        virtual ~DEMIntegrationScheme();

        virtual DEMIntegrationScheme* CloneRaw() const {
            DEMIntegrationScheme* cloned_scheme(new DEMIntegrationScheme(*this));
            return cloned_scheme;
        }

        virtual DEMIntegrationScheme::Pointer CloneShared() const {
            DEMIntegrationScheme::Pointer cloned_scheme(new DEMIntegrationScheme(*this));
            return cloned_scheme;
        }

        virtual void SetTranslationalIntegrationSchemeInProperties(Properties::Pointer pProp, bool verbose = true) const;
        virtual void SetRotationalIntegrationSchemeInProperties(Properties::Pointer pProp, bool verbose = true) const;

        virtual void Move(Node<3> & i, const double delta_t, const double force_reduction_factor, const int StepFlag);
        virtual void Rotate(Node<3> & i, const double delta_t, const double force_reduction_factor, const int StepFlag);
        virtual void MoveRigidBodyElement(RigidBodyElement3D* rigid_body_element, Node<3> & i, const double delta_t, const double force_reduction_factor, const int StepFlag);
        virtual void RotateRigidBodyElement(RigidBodyElement3D* rigid_body_element, Node<3> & i, const double delta_t, const double force_reduction_factor, const int StepFlag);

        virtual void UpdateTranslationalVariables(
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
                const bool Fix_vel[3]);

        virtual void CalculateTranslationalMotionOfNode(Node<3> & i, const double delta_t, const double force_reduction_factor, const int StepFlag);
        virtual void CalculateRotationalMotionOfSphereNode(Node<3> & i, const double delta_t, const double force_reduction_factor, const int StepFlag);
        virtual void CalculateRotationalMotionOfRigidBodyElementNode(Node<3> & i, const double delta_t, const double moment_reduction_factor, const int StepFlag);

        virtual void CalculateNewRotationalVariablesOfSpheres(
                int StepFlag,
                Node < 3 >& i,
                const double moment_of_inertia,
                array_1d<double, 3 >& angular_velocity,
                array_1d<double, 3 >& torque,
                const double moment_reduction_factor,
                array_1d<double, 3 >& rotated_angle,
                array_1d<double, 3 >& delta_rotation,
                const double delta_t,
                const bool Fix_Ang_vel[3]);

        virtual void CalculateNewRotationalVariablesOfRigidBodyElements(
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
                const bool Fix_Ang_vel[3]);

        virtual void UpdateRotationalVariables(
                int StepFlag,
                Node < 3 >& i,
                array_1d<double, 3 >& rotated_angle,
                array_1d<double, 3 >& delta_rotation,
                array_1d<double, 3 >& angular_velocity,
                array_1d<double, 3 >& angular_acceleration,
                const double delta_t,
                const bool Fix_Ang_vel[3]);

        virtual void UpdateRotationalVariables(
                int StepFlag,
                Node < 3 >& i,
                const double& moment_of_inertia,
                array_1d<double, 3 >& rotated_angle,
                array_1d<double, 3 >& delta_rotation,
                Quaternion<double  >& Orientation,
                const array_1d<double, 3 >& angular_momentum,
                array_1d<double, 3 >& angular_velocity,
                const double delta_t,
                const bool Fix_Ang_vel[3]);

        virtual void UpdateRotationalVariables(
                int StepFlag,
                Node < 3 >& i,
                const array_1d<double, 3 >& moments_of_inertia,
                array_1d<double, 3 >& rotated_angle,
                array_1d<double, 3 >& delta_rotation,
                Quaternion<double  >& Orientation,
                const array_1d<double, 3 >& angular_momentum,
                array_1d<double, 3 >& angular_velocity,
                const double delta_t,
                const bool Fix_Ang_vel[3]);

        virtual void UpdateRotatedAngle(
                array_1d<double, 3 >& rotated_angle,
                array_1d<double, 3 >& delta_rotation,
                const array_1d<double, 3 >& angular_velocity,
                const double delta_t);

        virtual void UpdateAngularVelocity(
                const Quaternion<double>& Orientation,
                const double LocalTensorInv[3][3],
                const array_1d<double, 3>& angular_momentum,
                array_1d<double, 3>& angular_velocity);

        virtual void CalculateLocalAngularAcceleration(
                const double moment_of_inertia,
                const array_1d<double, 3 >& torque,
                const double moment_reduction_factor,
                array_1d<double, 3 >& angular_acceleration);

        virtual void CalculateLocalAngularAccelerationByEulerEquations(
                const array_1d<double, 3 >& local_angular_velocity,
                const array_1d<double, 3 >& moments_of_inertia,
                const array_1d<double, 3 >& local_torque,
                const double moment_reduction_factor,
                array_1d<double, 3 >& local_angular_acceleration);

        virtual void CalculateAngularVelocityRK(
                const Quaternion<double  >& Orientation,
                const double& moment_of_inertia,
                const array_1d<double, 3 >& angular_momentum,
                array_1d<double, 3 >& angular_velocity,
                const double delta_t,
                const bool Fix_Ang_vel[3]);

        virtual void CalculateAngularVelocityRK(
                const Quaternion<double  >& Orientation,
                const array_1d<double, 3 >& moments_of_inertia,
                const array_1d<double, 3 >& angular_momentum,
                array_1d<double, 3 > & angular_velocity,
                const double delta_t,
                const bool Fix_Ang_vel[3]);

        virtual void QuaternionCalculateMidAngularVelocities(
                const Quaternion<double>& Orientation,
                const double LocalTensorInv[3][3],
                const array_1d<double, 3>& angular_momentum,
                const double dt,
                const array_1d<double, 3>& InitialAngularVel,
                array_1d<double, 3>& FinalAngularVel);

        virtual std::string Info() const {
            std::stringstream buffer;
            buffer << "DEMIntegrationScheme";
            return buffer.str();
        }

        /// Print information about this object.

        virtual void PrintInfo(std::ostream & rOStream) const {
            rOStream << "DEMIntegrationScheme";
        }

        /// Print object's data.

        virtual void PrintData(std::ostream & rOStream) const {
        }

        protected:

        private:

        //bool mRotationOption;

        DEMIntegrationScheme& operator=(DEMIntegrationScheme const& rOther) {
            return *this;
        }

        /// Copy constructor.

        DEMIntegrationScheme(DEMIntegrationScheme const& rOther) {
            *this = rOther;
        }

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const {
                    //rSerializer.save("MyMemberName",myMember);
        }

        virtual void load(Serializer& rSerializer) {
                    //rSerializer.load("MyMemberName",myMember);
        }

    }; // Class DEMIntegrationScheme

    /// input stream function

    inline std::istream& operator>>(std::istream& rIStream, DEMIntegrationScheme& rThis) {
        return rIStream;
    }

    /// output stream function
    inline std::ostream& operator<<(std::ostream& rOStream, const DEMIntegrationScheme& rThis) {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);
        return rOStream;
    }

} // namespace Kratos

#endif // KRATOS_DEM_INTEGRATION_SCHEME_H_INCLUDED defined

