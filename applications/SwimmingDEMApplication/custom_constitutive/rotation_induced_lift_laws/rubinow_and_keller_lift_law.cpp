// Author: Guillermo Casas (gcasas@cimne.upc.edu)
// Date: February 2019

#include "swimming_DEM_application.h"
#include "rubinow_and_keller_lift_law.h"

namespace Kratos {

    RotationInducedLiftLaw::Pointer RubinowAndKellerLiftLaw::Clone() const {
        RubinowAndKellerLiftLaw::Pointer p_clone(new RubinowAndKellerLiftLaw(*this));
        return p_clone;
    }

    RubinowAndKellerLiftLaw::RubinowAndKellerLiftLaw(Parameters r_parameters)
    {

    }

    void RubinowAndKellerLiftLaw::Initialize(const ProcessInfo& r_process_info){}

    std::string RubinowAndKellerLiftLaw::GetTypeOfLaw() {
        std::string type_of_law = "Rubinow and Keller lift law";
        return type_of_law;
    }

    void RubinowAndKellerLiftLaw::ComputeForce(Geometry<Node >& r_geometry,
                                               const double reynolds_number,
                                               double particle_radius,
                                               double fluid_density,
                                               double fluid_kinematic_viscosity,
                                               array_1d<double, 3>& minus_slip_velocity,
                                               array_1d<double, 3>& rotation_induced_lift,
                                               const ProcessInfo& r_current_process_info)
    {
        Node& node = r_geometry[0];
        const array_1d<double, 3> minus_slip_rot = (0.5 * node.FastGetSolutionStepValue(FLUID_VORTICITY_PROJECTED)
                                              - node.FastGetSolutionStepValue(ANGULAR_VELOCITY));
        array_1d<double, 3> minus_slip_rot_cross_slip_vel;
        SWIMMING_SET_TO_CROSS_OF_FIRST_TWO_3(minus_slip_rot, minus_slip_velocity, minus_slip_rot_cross_slip_vel)

        // Rubinow and Keller, 1961 (Re_p < 0.1; nondimensional_minus_slip_rot_vel < 0.1)
        rotation_induced_lift = Globals::Pi * SWIMMING_POW_3(particle_radius) * fluid_density * minus_slip_rot_cross_slip_vel;
    }

} // namespace Kratos
