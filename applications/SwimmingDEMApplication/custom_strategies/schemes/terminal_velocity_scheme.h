//
// Author: Guillermo Casas gcasas@cimne.upc.edu
//

#if !defined(KRATOS_TERMINAL_VELOCITY_SCHEME_H_INCLUDED )
#define  KRATOS_TERMINAL_VELOCITY_SCHEME_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <cfloat>

// Project includes
#include "hybrid_bashforth_scheme.h"
#include "includes/define.h"
#include "utilities/openmp_utils.h"
#include "includes/model_part.h"
#include "../DEMApplication/custom_utilities/GeometryFunctions.h"
#include "utilities/quaternion.h"

namespace Kratos {

    class KRATOS_API(SWIMMING_DEM_APPLICATION) TerminalVelocityScheme : public HybridBashforthScheme {
    public:

        typedef ModelPart::NodesContainerType NodesArrayType;

        /// Pointer definition of TerminalVelocityScheme
        KRATOS_CLASS_POINTER_DEFINITION(TerminalVelocityScheme);

        /// Default constructor.
        TerminalVelocityScheme() {}

        /// Destructor.
        virtual ~TerminalVelocityScheme() {}

        DEMIntegrationScheme* CloneRaw() const override {
            DEMIntegrationScheme* cloned_scheme(new TerminalVelocityScheme(*this));
            return cloned_scheme;
        }

        DEMIntegrationScheme::Pointer CloneShared() const override {
            DEMIntegrationScheme::Pointer cloned_scheme(new TerminalVelocityScheme(*this));
            return cloned_scheme;
        }

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

        void UpdateRotationalVariables(
                int StepFlag,
                Node < 3 >& i,
                array_1d<double, 3 >& rotated_angle,
                array_1d<double, 3 >& delta_rotation,
                array_1d<double, 3 >& angular_velocity,
                array_1d<double, 3 >& angular_acceleration,
                const double delta_t,
                const bool Fix_Ang_vel[3]) override;

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

        virtual std::string Info() const override {
            std::stringstream buffer;
            buffer << "SymplecticEulerScheme";
            return buffer.str();
        }

        /// Print information about this object.

        virtual void PrintInfo(std::ostream& rOStream) const override {
            rOStream << "TerminalVelocityScheme";
        }

        /// Print object's data.

        virtual void PrintData(std::ostream& rOStream) const override {
        }


    protected:


    private:

    /// Assignment operator.

        TerminalVelocityScheme& operator=(TerminalVelocityScheme const& rOther) {
            return *this;
        }

        /// Copy constructor.

        TerminalVelocityScheme(TerminalVelocityScheme const& rOther) {
            *this = rOther;
        }

        ///@}

    }; // Class TerminalVelocityScheme

    inline std::istream& operator>>(std::istream& rIStream,
            TerminalVelocityScheme& rThis) {
        return rIStream;
    }

    inline std::ostream& operator<<(std::ostream& rOStream,
            const TerminalVelocityScheme& rThis) {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }

} // namespace Kratos.

#endif // KRATOS_TERMINAL_VELOCITY_SCHEME_H_INCLUDED  defined
