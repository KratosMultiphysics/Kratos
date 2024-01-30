// Author: Guillermo Casas (gcasas@cimne.upc.edu)
// Date: February 2019

#include "swimming_DEM_application.h"
#include "loth_rotation_induced_lift_law.h"

namespace Kratos {

    RotationInducedLiftLaw::Pointer LothRotationInducedLiftLaw::Clone() const {
        LothRotationInducedLiftLaw::Pointer p_clone(new LothRotationInducedLiftLaw(*this));
        return p_clone;
    }

    LothRotationInducedLiftLaw::LothRotationInducedLiftLaw(Parameters r_parameters)
    {

    }

    void LothRotationInducedLiftLaw::Initialize(const ProcessInfo& r_process_info){}

    std::string LothRotationInducedLiftLaw::GetTypeOfLaw() {
        std::string type_of_law = "Loth rotation induced lift law";
        return type_of_law;
    }

    void LothRotationInducedLiftLaw::ComputeForce(Geometry<Node >& r_geometry,
                                                  const double reynolds_number,
                                                  double particle_radius,
                                                  double fluid_density,
                                                  double fluid_kinematic_viscosity,
                                                  array_1d<double, 3>& minus_slip_velocity,
                                                  array_1d<double, 3>& rotation_induced_lift,
                                                  const ProcessInfo& r_current_process_info)
    {
        // First compute Rubinow and Keller's lift
        RubinowAndKellerLiftLaw::ComputeForce(r_geometry,
                                              reynolds_number,
                                              particle_radius,
                                              fluid_density,
                                              fluid_kinematic_viscosity,
                                              minus_slip_velocity,
                                              rotation_induced_lift,
                                              r_current_process_info);

        // Then apply Loth's correction
        Node& node = r_geometry[0];
        const array_1d<double, 3> minus_slip_rot = (0.5 * node.FastGetSolutionStepValue(FLUID_VORTICITY_PROJECTED)
                                              - node.FastGetSolutionStepValue(ANGULAR_VELOCITY));
        const double norm_of_slip_vel = SWIMMING_MODULUS_3(minus_slip_velocity);
        const double norm_of_slip_rot = SWIMMING_MODULUS_3(minus_slip_rot);
        const double nondimensional_minus_slip_rot_vel = this->ComputeNondimensionalRotVelocity(norm_of_slip_vel,
                                                                                          norm_of_slip_rot,
                                                                                          particle_radius,
                                                                                          fluid_kinematic_viscosity);
        const double coeff = 1 - (0.675 + 0.15 * (1 + std::tanh(0.28 * (nondimensional_minus_slip_rot_vel - 2)))) * std::tanh(0.18 * std::sqrt(reynolds_number));

        rotation_induced_lift *= coeff;
    }

} // namespace Kratos
