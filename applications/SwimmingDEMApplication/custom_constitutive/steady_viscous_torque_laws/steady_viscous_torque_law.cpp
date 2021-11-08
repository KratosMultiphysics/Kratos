// Author: Guillermo Casas (gcasas@cimne.upc.edu)
// Date: February 2019

#include "swimming_DEM_application.h"
#include "steady_viscous_torque_law.h"

namespace Kratos {

    void SteadyViscousTorqueLaw::Initialize(const ProcessInfo& r_process_info) {}

    std::string SteadyViscousTorqueLaw::GetTypeOfLaw() {
        std::string type_of_law = "Generic steady viscous torque law";
        return type_of_law;
    }

    void SteadyViscousTorqueLaw::SetSteadyViscousTorqueLawInProperties(Properties::Pointer pProp) const {
    }

    SteadyViscousTorqueLaw::Pointer SteadyViscousTorqueLaw::Clone() const {
        SteadyViscousTorqueLaw::Pointer p_clone(new SteadyViscousTorqueLaw(*this));
        return p_clone;
    }

    double SteadyViscousTorqueLaw::ComputeParticleRotationReynoldsNumber(const double norm_of_slip_rot,
                                                                         const double particle_radius,
                                                                         const double fluid_kinematic_viscosity)
    {
        return 4 * SWIMMING_POW_2(particle_radius) * norm_of_slip_rot / fluid_kinematic_viscosity;
    }

    double SteadyViscousTorqueLaw::ComputeNondimensionalRotVelocity(const double norm_of_slip_vel,
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
