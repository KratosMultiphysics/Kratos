// Author: Guillermo Casas (gcasas@cimne.upc.edu)
// Date: February 2019

#include "swimming_DEM_application.h"
#include "oesterle_dinh_lift_law.h"

namespace Kratos {

    RotationInducedLiftLaw::Pointer OesterleAndDinhLiftLaw::Clone() const {
        OesterleAndDinhLiftLaw::Pointer p_clone(new OesterleAndDinhLiftLaw(*this));
        return p_clone;
    }

    OesterleAndDinhLiftLaw::OesterleAndDinhLiftLaw(Parameters& r_parameters)
    {

    }

    void OesterleAndDinhLiftLaw::Initialize(const ProcessInfo& r_process_info){}

    std::string OesterleAndDinhLiftLaw::GetTypeOfLaw() {
        std::string type_of_law = "Oesterle and Dinh lift law";
        return type_of_law;
    }

    void OesterleAndDinhLiftLaw::ComputeForce(Geometry<Node<3> >& r_geometry,
                                               const double reynolds_number,
                                               double particle_radius,
                                               double fluid_density,
                                               double fluid_kinematic_viscosity,
                                               array_1d<double, 3>& minus_slip_velocity,
                                               array_1d<double, 3>& rotation_induced_lift,
                                               const ProcessInfo& r_current_process_info)
    {
        Node<3>& node = r_geometry[0];
        const array_1d<double, 3> minus_slip_rot = (0.5 * node.FastGetSolutionStepValue(FLUID_VORTICITY_PROJECTED)
                                              - node.FastGetSolutionStepValue(ANGULAR_VELOCITY));
        array_1d<double, 3> minus_slip_rot_cross_slip_vel;
        SWIMMING_SET_TO_CROSS_OF_FIRST_TWO_3(minus_slip_rot, minus_slip_velocity, minus_slip_rot_cross_slip_vel)

        const double norm_of_slip_vel = SWIMMING_MODULUS_3(minus_slip_velocity);
        const double norm_of_slip_rot = SWIMMING_MODULUS_3(minus_slip_rot);
        const double rot_reynolds = this->ComputeParticleRotationReynoldsNumber(norm_of_slip_rot,
                                                                                particle_radius,
                                                                                fluid_kinematic_viscosity);

        if (std::abs(reynolds_number) < std::numeric_limits<double>::epsilon() || std::abs(rot_reynolds) < std::numeric_limits<double>::epsilon()){
            return;
        }

        else {
            const double lift_coeff = 0.45  + (rot_reynolds / reynolds_number - 0.45) * std::exp(- 0.05684 * std::pow(rot_reynolds, 0.4) * std::pow(reynolds_number, 0.3));
            noalias(rotation_induced_lift) = 0.5 *  fluid_density * Globals::Pi * SWIMMING_POW_2(particle_radius) * lift_coeff * norm_of_slip_vel * minus_slip_rot_cross_slip_vel / norm_of_slip_rot;
        }
    }

} // namespace Kratos
