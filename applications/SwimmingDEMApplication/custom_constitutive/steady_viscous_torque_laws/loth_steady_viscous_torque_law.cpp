// Author: Guillermo Casas (gcasas@cimne.upc.edu)
// Date: February 2019

#include "swimming_DEM_application.h"
#include "loth_steady_viscous_torque_law.h"

namespace Kratos {

    SteadyViscousTorqueLaw::Pointer LothSteadyViscousTorqueLaw::Clone() const {
        LothSteadyViscousTorqueLaw::Pointer p_clone(new LothSteadyViscousTorqueLaw(*this));
        return p_clone;
    }

    LothSteadyViscousTorqueLaw::LothSteadyViscousTorqueLaw(Parameters r_parameters)
    {

    }

    void LothSteadyViscousTorqueLaw::Initialize(const ProcessInfo& r_process_info){}

    std::string LothSteadyViscousTorqueLaw::GetTypeOfLaw() {
        std::string type_of_law = "Loth steady viscous torque law";
        return type_of_law;
    }

    void LothSteadyViscousTorqueLaw::ComputeMoment(Geometry<Node >& r_geometry,
                                                  const double reynolds_number,
                                                  double particle_radius,
                                                  double fluid_density,
                                                  double fluid_kinematic_viscosity,
                                                  array_1d<double, 3>& minus_slip_velocity,
                                                  array_1d<double, 3>& steady_viscous_torque,
                                                  const ProcessInfo& r_current_process_info)
    {
        Node& node = r_geometry[0];
        const array_1d<double, 3> minus_slip_rot = (0.5 * node.FastGetSolutionStepValue(FLUID_VORTICITY_PROJECTED)
                                              - node.FastGetSolutionStepValue(ANGULAR_VELOCITY));
        const double norm_of_slip_rot = SWIMMING_MODULUS_3(minus_slip_rot);

        if (norm_of_slip_rot){
            // First compute Rubinow and Keller's torque
            RubinowAndKellerTorqueLaw::ComputeMoment(r_geometry,
                                                    reynolds_number,
                                                    particle_radius,
                                                    fluid_density,
                                                    fluid_kinematic_viscosity,
                                                    minus_slip_velocity,
                                                    steady_viscous_torque,
                                                    r_current_process_info);

            // Then apply Loth's correction
            const double norm_of_slip_rot = SWIMMING_MODULUS_3(minus_slip_rot);
            const double rot_reynolds = this->ComputeParticleRotationReynoldsNumber(norm_of_slip_rot,
                                                                                    particle_radius,
                                                                                    fluid_kinematic_viscosity) / norm_of_slip_rot;

            steady_viscous_torque *= 1.0 + 5 / (64 * Globals::Pi) * std::pow(rot_reynolds, 0.6);
        }
    }

} // namespace Kratos
