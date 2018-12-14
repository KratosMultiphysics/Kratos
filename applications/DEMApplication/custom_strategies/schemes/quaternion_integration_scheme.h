//
// Author: Joaquin Irazabal, jirazabal@cimne.upc.edu
//

#if !defined(KRATOS_QUATERNION_INTEGRATION_SCHEME_H_INCLUDED )
#define  KRATOS_QUATERNION_INTEGRATION_SCHEME_H_INCLUDED

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

    class KRATOS_API(DEM_APPLICATION) QuaternionIntegrationScheme : public DEMIntegrationScheme {
    public:

        typedef ModelPart::NodesContainerType NodesArrayType;

        /// Pointer definition of QuaternionIntegrationScheme
        KRATOS_CLASS_POINTER_DEFINITION(QuaternionIntegrationScheme);

        /// Default constructor.
        QuaternionIntegrationScheme() {}

        /// Destructor.
        virtual ~QuaternionIntegrationScheme() {}

        DEMIntegrationScheme* CloneRaw() const override {
            DEMIntegrationScheme* cloned_scheme(new QuaternionIntegrationScheme(*this));
            return cloned_scheme;
        }

        DEMIntegrationScheme::Pointer CloneShared() const override {
            DEMIntegrationScheme::Pointer cloned_scheme(new QuaternionIntegrationScheme(*this));
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

        void UpdateRotatedAngle(
                array_1d<double, 3 >& rotated_angle,
                array_1d<double, 3 >& delta_rotation,
                const array_1d<double, 3 >& angular_velocity,
                const double delta_t) override;

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

        /// Turn back information as a string.

        virtual std::string Info() const override{
            std::stringstream buffer;
            buffer << "QuaternionIntegrationScheme";
            return buffer.str();
        }

        /// Print information about this object.

        virtual void PrintInfo(std::ostream& rOStream) const override{
            rOStream << "QuaternionIntegrationScheme";
        }

        /// Print object's data.

        virtual void PrintData(std::ostream& rOStream) const override{
        }


    protected:


    private:

        /// Assignment operator.

        QuaternionIntegrationScheme& operator=(QuaternionIntegrationScheme const& rOther) {
            return *this;
        }

        /// Copy constructor.

        QuaternionIntegrationScheme(QuaternionIntegrationScheme const& rOther) {
            *this = rOther;
        }


        ///@}

    }; // Class QuaternionIntegrationScheme


    inline std::istream& operator>>(std::istream& rIStream,
            QuaternionIntegrationScheme& rThis) {
        return rIStream;
    }

    inline std::ostream& operator<<(std::ostream& rOStream,
            const QuaternionIntegrationScheme& rThis) {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }

} // namespace Kratos.

#endif // KRATOS_QUATERNION_INTEGRATION_SCHEME_H_INCLUDED  defined
