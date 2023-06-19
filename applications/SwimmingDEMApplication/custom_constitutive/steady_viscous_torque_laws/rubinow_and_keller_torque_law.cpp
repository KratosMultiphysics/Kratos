// Author: Guillermo Casas (gcasas@cimne.upc.edu)
// Date: February 2019

#include "swimming_DEM_application.h"
#include "rubinow_and_keller_torque_law.h"

namespace Kratos {

    SteadyViscousTorqueLaw::Pointer RubinowAndKellerTorqueLaw::Clone() const {
        RubinowAndKellerTorqueLaw::Pointer p_clone(new RubinowAndKellerTorqueLaw(*this));
        return p_clone;
    }

    RubinowAndKellerTorqueLaw::RubinowAndKellerTorqueLaw(Parameters r_parameters)
    {

    }

    void RubinowAndKellerTorqueLaw::Initialize(const ProcessInfo& r_process_info){}

    std::string RubinowAndKellerTorqueLaw::GetTypeOfLaw() {
        std::string type_of_law = "Rubinow and Keller torque law";
        return type_of_law;
    }

    void RubinowAndKellerTorqueLaw::ComputeMoment(Geometry<Node >& r_geometry,
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
            array_1d<double, 3> minus_slip_rot_cross_slip_vel;
            SWIMMING_SET_TO_CROSS_OF_FIRST_TWO_3(minus_slip_rot, minus_slip_velocity, minus_slip_rot_cross_slip_vel)

            const double rot_reynolds = this->ComputeParticleRotationReynoldsNumber(norm_of_slip_rot,
                                                                                    particle_radius,
                                                                                    fluid_kinematic_viscosity) / norm_of_slip_rot;

            double rotational_coeff;

            if (rot_reynolds > 32){ // Rubinow and Keller, 1961 (Re_rot ~ 32 - 1000)
                rotational_coeff = 12.9 * std::sqrt(norm_of_slip_rot * rot_reynolds) + 128.4 / rot_reynolds;
            }

            else { // Rubinow and Keller, 1961 (Re_rot < 32)
                rotational_coeff = 64 * Globals::Pi / rot_reynolds;
            }

            noalias(steady_viscous_torque) = 0.5 * fluid_density * SWIMMING_POW_5(particle_radius) * rotational_coeff * minus_slip_rot;
        }

    }

} // namespace Kratos
