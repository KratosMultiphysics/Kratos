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
#include "includes/kratos_parameters.h"
#include "utilities/openmp_utils.h"
#include "includes/model_part.h"
#include "custom_utilities/GeometryFunctions.h"
#include "utilities/quaternion.h"

namespace Kratos {

    class KRATOS_API(SWIMMING_DEM_APPLICATION) TerminalVelocityScheme : public HybridBashforthScheme {
    public:

        typedef ModelPart::NodesContainerType NodesArrayType;

        /// Pointer definition of TerminalVelocityScheme
        KRATOS_CLASS_POINTER_DEFINITION(TerminalVelocityScheme);

        /// Default constructor: the scheme is NOT usable until the fluid viscosity and
        /// the gravity have been given (see the Parameters constructor); a move attempt
        /// with an unconfigured scheme raises an error instead of using hidden constants.
        TerminalVelocityScheme()
            : mDynamicViscosity(0.0), mGravity(ZeroVector(3)), mIsConfigured(false) {}

        /// Constructor from Parameters. Expected fields (both required, no hidden defaults):
        ///   "dynamic_viscosity" : fluid dynamic viscosity used in the Stokes drag 6 pi mu a
        ///   "gravity"           : gravity vector [gx, gy, gz]
        /// All quantities are taken as given, in whatever (dimensional or dimensionless)
        /// system of units the rest of the case uses.
        explicit TerminalVelocityScheme(Parameters rParameters);

        /// Destructor.
        virtual ~TerminalVelocityScheme() {}

        /// Copy: the physical parameters travel with the clone stored in the Properties
        /// (the base schemes keep their copy constructors private and carry no state
        /// that has to be copied, so the base is default-constructed).
        TerminalVelocityScheme(TerminalVelocityScheme const& rOther)
            : HybridBashforthScheme(),
              mDynamicViscosity(rOther.mDynamicViscosity),
              mGravity(rOther.mGravity),
              mIsConfigured(rOther.mIsConfigured) {}
        TerminalVelocityScheme& operator=(TerminalVelocityScheme const& rOther) {
            mDynamicViscosity = rOther.mDynamicViscosity;
            noalias(mGravity) = rOther.mGravity;
            mIsConfigured = rOther.mIsConfigured;
            return *this;
        }

        /// Explicit setters/getters (Python convenience and checks).
        void SetDynamicViscosity(const double Viscosity);
        void SetGravity(const array_1d<double, 3>& rGravity);
        double GetDynamicViscosity() const { return mDynamicViscosity; }
        const array_1d<double, 3>& GetGravity() const { return mGravity; }
        bool IsConfigured() const { return mIsConfigured; }

        static Parameters GetDefaultParameters()
        {
            return Parameters(R"({
                "dynamic_viscosity" : 0.0,
                "gravity"           : [0.0, 0.0, 0.0]
            })");
        }

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
                Node& i,
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
                Node& i,
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
            buffer << "TerminalVelocityScheme";
            return buffer.str();
        }

        /// Print information about this object.

        virtual void PrintInfo(std::ostream& rOStream) const override {
            rOStream << "TerminalVelocityScheme";
        }

        /// Print object's data.

        virtual void PrintData(std::ostream& rOStream) const override {
            rOStream << "dynamic_viscosity = " << mDynamicViscosity << ", gravity = " << mGravity
                     << (mIsConfigured ? "" : " (NOT CONFIGURED)");
        }


    protected:


    private:

        double mDynamicViscosity;      // fluid dynamic viscosity mu (Stokes drag 6 pi mu a)
        array_1d<double, 3> mGravity;  // gravity vector
        bool mIsConfigured;            // both values were given (Parameters ctor or setters)

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
