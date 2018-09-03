//
// Author: Joaquin Irazabal, jirazabal@cimne.upc.edu
//

#if !defined(KRATOS_RUNGE_KUTTA_SCHEME_H_INCLUDED )
#define  KRATOS_RUNGE_KUTTA_SCHEME_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <cfloat>

// Project includes
#include "dem_integration_scheme.h"
#include "includes/define.h"
#include "utilities/openmp_utils.h"
#include "includes/model_part.h"
#include "custom_utilities/GeometryFunctions.h"
#include "utilities/quaternion.h"

namespace Kratos {

    class KRATOS_API(DEM_APPLICATION) RungeKuttaScheme : public DEMIntegrationScheme {
    public:

        typedef ModelPart::NodesContainerType NodesArrayType;

        /// Pointer definition of RungeKuttaScheme
        KRATOS_CLASS_POINTER_DEFINITION(RungeKuttaScheme);

        /// Default constructor.
        RungeKuttaScheme() {}

        /// Destructor.
        virtual ~RungeKuttaScheme() {}

        DEMIntegrationScheme* CloneRaw() const override {
            DEMIntegrationScheme* cloned_scheme(new RungeKuttaScheme(*this));
            return cloned_scheme;
        }

        DEMIntegrationScheme::Pointer CloneShared() const override {
            DEMIntegrationScheme::Pointer cloned_scheme(new RungeKuttaScheme(*this));
            return cloned_scheme;
        }

        void SetTranslationalIntegrationSchemeInProperties(Properties::Pointer pProp, bool verbose = true) const override;
        void SetRotationalIntegrationSchemeInProperties(Properties::Pointer pProp, bool verbose = true) const override;

        void UpdateTranslationalVariables(
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
                const bool Fix_vel[3]) override;

        void CalculateNewRotationalVariablesOfSpheres(
                int StepFlag,
                Node < 3 >& i,
                const double moment_of_inertia,
                array_1d<double, 3 >& angular_velocity,
                array_1d<double, 3 >& torque,
                const double moment_reduction_factor,
                array_1d<double, 3 >& rotated_angle,
                array_1d<double, 3 >& delta_rotation,
                const double delta_t,
                const bool Fix_Ang_vel[3]) override;

        void CalculateNewRotationalVariablesOfRigidBodyElements(
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
                const bool Fix_Ang_vel[3]) override;

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
                const bool Fix_Ang_vel[3]) override;

        void UpdateRotationalVariables(
                int StepFlag,
                Node < 3 >& i,
                const array_1d<double, 3 >& moments_of_inertia,
                array_1d<double, 3 >& rotated_angle,
                array_1d<double, 3 >& delta_rotation,
                Quaternion<double  >& Orientation,
                const array_1d<double, 3 >& angular_momentum,
                array_1d<double, 3 >& angular_velocity,
                const double delta_t,
                const bool Fix_Ang_vel[3]) override;

        void UpdateAngularVelocity(
                const Quaternion<double>& Orientation,
                const double LocalTensorInv[3][3],
                const array_1d<double, 3>& angular_momentum,
                array_1d<double, 3>& angular_velocity)  override;

        void CalculateLocalAngularAcceleration(
                const double moment_of_inertia,
                const array_1d<double, 3 >& torque,
                const double moment_reduction_factor,
                array_1d<double, 3 >& angular_acceleration) override;

        void CalculateLocalAngularAccelerationByEulerEquations(
                const array_1d<double, 3 >& local_angular_velocity,
                const array_1d<double, 3 >& moments_of_inertia,
                const array_1d<double, 3 >& local_torque,
                const double moment_reduction_factor,
                array_1d<double, 3 >& local_angular_acceleration) override;

        void CalculateAngularVelocityRK(
                const Quaternion<double  >& Orientation,
                const double& moment_of_inertia,
                const array_1d<double, 3 >& angular_momentum,
                array_1d<double, 3 >& angular_velocity,
                const double delta_t,
                const bool Fix_Ang_vel[3]) override;

        void CalculateAngularVelocityRK(
                const Quaternion<double  >& Orientation,
                const array_1d<double, 3 >& moments_of_inertia,
                const array_1d<double, 3 >& angular_momentum,
                array_1d<double, 3 > & angular_velocity,
                const double delta_t,
                const bool Fix_Ang_vel[3]) override;

        void QuaternionCalculateMidAngularVelocities(
                const Quaternion<double>& Orientation,
                const double LocalTensorInv[3][3],
                const array_1d<double, 3>& angular_momentum,
                const double dt,
                const array_1d<double, 3>& InitialAngularVel,
                array_1d<double, 3>& FinalAngularVel) override;

        /// Turn back information as a string.

        virtual std::string Info() const override{
            std::stringstream buffer;
            buffer << "RungeKuttaScheme";
            return buffer.str();
        }

        /// Print information about this object.

        virtual void PrintInfo(std::ostream& rOStream) const override{
            rOStream << "RungeKuttaScheme";
        }

        /// Print object's data.

        virtual void PrintData(std::ostream& rOStream) const override{
        }


    protected:


    private:

        /// Assignment operator.

        RungeKuttaScheme& operator=(RungeKuttaScheme const& rOther) {
            return *this;
        }

        /// Copy constructor.

        RungeKuttaScheme(RungeKuttaScheme const& rOther) {
            *this = rOther;
        }


        ///@}

    }; // Class RungeKuttaScheme


    inline std::istream& operator>>(std::istream& rIStream,
            RungeKuttaScheme& rThis) {
        return rIStream;
    }

    inline std::ostream& operator<<(std::ostream& rOStream,
            const RungeKuttaScheme& rThis) {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }

} // namespace Kratos.

#endif // KRATOS_RUNGE_KUTTA_SCHEME_H_INCLUDED  defined
