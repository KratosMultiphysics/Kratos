//
// Author: Miquel Santasusana msantasusana@cimne.upc.edu
//

#if !defined(KRATOS_TERMINAL_VELOCITY_SCHEME_H_INCLUDED )
#define  KRATOS_TERMINAL_VELOCITY_SCHEME_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <cfloat>

// Project includes

#include "../DEM_application/custom_strategies/schemes/dem_integration_scheme.h"
#include "includes/define.h"
#include "utilities/openmp_utils.h"
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"
#include "swimming_DEM_application.h"
#include "../DEM_application/DEM_application.h"



//#include "includes/define.h"
//#include "includes/kratos_application.h"
//#include "includes/variables.h"
//#include "includes/dem_variables.h"  //TODO: must be removed eventually
//#include "includes/cfd_variables.h"  //TODO: must be removed eventually
//#include "includes/legacy_structural_app_vars.h"  //TODO: must be removed eventually
//#include "custom_elements/monolithic_dem_coupled.h"
//#include "custom_elements/monolithic_dem_coupled_weak.h"
//#include "custom_elements/spheric_swimming_particle.h"
//#include "../DEM_application/custom_elements/spheric_particle.h"
//#include "../DEM_application/custom_elements/nanoparticle.h"

namespace Kratos {

    class TerminalVelocityScheme : public DEMIntegrationScheme
    {
    public:

        typedef ModelPart::NodesContainerType NodesArrayType;

        /// Pointer definition of TerminalVelocityScheme
        KRATOS_CLASS_POINTER_DEFINITION(TerminalVelocityScheme);

        /// Default constructor.
        TerminalVelocityScheme();

        /// Destructor.
        virtual ~TerminalVelocityScheme();

        void UpdateTranslationalVariables(
            const Node < 3 > & i,
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

        void UpdateRotationalVariables(
                const Node < 3 > & i,
                array_1d<double, 3 >& rotated_angle,
                array_1d<double, 3 >& delta_rotation,
                array_1d<double, 3 >& angular_velocity,
                array_1d<double, 3 >& angular_acceleration,
                const double delta_t,
                const bool Fix_Ang_vel[3]);

        void CalculateLocalAngularAcceleration(
                                const Node < 3 > & i,
                                const double moment_of_inertia,
                                const array_1d<double, 3 >& torque,
                                const double moment_reduction_factor,
                                array_1d<double, 3 >& angular_acceleration);

        void CalculateLocalAngularAccelerationByEulerEquations(
                                    const Node < 3 > & i,
                                    const array_1d<double, 3 >& local_angular_velocity,
                                    const array_1d<double, 3 >& moments_of_inertia,
                                    const array_1d<double, 3 >& local_torque,
                                    const double moment_reduction_factor,
                                    array_1d<double, 3 >& local_angular_acceleration);

        /// Turn back information as a string.

        virtual std::string Info() const {
            std::stringstream buffer;
            buffer << "TerminalVelocityScheme";
            return buffer.str();
        }

        /// Print information about this object.

        virtual void PrintInfo(std::ostream& rOStream) const {
            rOStream << "TerminalVelocityScheme";
        }

        /// Print object's data.

        virtual void PrintData(std::ostream& rOStream) const {
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
