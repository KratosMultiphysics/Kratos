// Author: Guillermo Casas (gcasas@cimne.upc.edu)
// Date: February 2019

#include "swimming_DEM_application.h"
#include "rotation_induced_lift_law.h"

namespace Kratos {

    void RotationInducedLiftLaw::Initialize(const ProcessInfo& r_process_info) {}

    std::string RotationInducedLiftLaw::GetTypeOfLaw() {
        std::string type_of_law = "Generic rotation-induced lift law";
        return type_of_law;
    }

    void RotationInducedLiftLaw::SetRotationInducedLiftLawInProperties(Properties::Pointer pProp) const {
    }

    RotationInducedLiftLaw::Pointer RotationInducedLiftLaw::Clone() const {
        RotationInducedLiftLaw::Pointer p_clone(new RotationInducedLiftLaw(*this));
        return p_clone;
    }

    double RotationInducedLiftLaw::ComputeParticleRotationReynoldsNumber(const double norm_of_slip_rot,
                                                                         const double particle_radius,
                                                                         const double fluid_kinematic_viscosity)
    {
        return 4 * SWIMMING_POW_2(particle_radius) * norm_of_slip_rot / fluid_kinematic_viscosity;
    }

    double RotationInducedLiftLaw::ComputeNondimensionalRotVelocity(const double norm_of_slip_vel,
                                                                    const double norm_of_slip_rot,
                                                                    const double particle_radius,
                                                                    const double fluid_kinematic_viscosity)
    {
        if (norm_of_slip_vel > 0.0){
            return 2.0 * particle_radius * norm_of_slip_rot / norm_of_slip_vel;
        }

        else {
            return 0.0;
        }
    }

} // namespace Kratos
