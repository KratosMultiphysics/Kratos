// Author: Guillermo Casas (gcasas@cimne.upc.edu)
// Date: February 2019

#include "swimming_DEM_application.h"
#include "el_samni_lift_law.h"

namespace Kratos {

    VorticityInducedLiftLaw::Pointer ElSamniLiftLaw::Clone() const {
        ElSamniLiftLaw::Pointer p_clone(new ElSamniLiftLaw(*this));
        return p_clone;
    }

    ElSamniLiftLaw::ElSamniLiftLaw(Parameters r_parameters)
    {

    }

    void ElSamniLiftLaw::Initialize(const ProcessInfo& r_process_info){}

    std::string ElSamniLiftLaw::GetTypeOfLaw() {
        std::string type_of_law = "El Samni lift law";
        return type_of_law;
    }

    double ElSamniLiftLaw::ComputeElSamniLiftCoefficient(const double particle_radius,
                                                         const double fluid_density,
                                                         const double norm_of_slip_vel,
                                                         const double vorticity_norm,
                                                         const ProcessInfo& r_current_process_info)
    {
        if (vorticity_norm > 0.000000000001 && norm_of_slip_vel > 0.000000000001){
            const double yield_stress = 0.0; // we are considering a Bingham type fluid
            const double power_law_K = r_current_process_info[POWER_LAW_K];
            const double power_law_n = r_current_process_info[POWER_LAW_N];
            const double shear_rate_p = norm_of_slip_vel / particle_radius * (4.5 / power_law_n - 3.5); // graphic model by Unhlherr et al. (fit by Wallis, G.B. and Dobson, J.E., 1973)
            const double equivalent_viscosity = yield_stress / shear_rate_p + power_law_K * pow(shear_rate_p, power_law_n - 1);
            const double coeff = std::max(0.09 * norm_of_slip_vel, 5.82 * std::sqrt(0.5 * norm_of_slip_vel * equivalent_viscosity /  fluid_density));
            const double lift_coeff = 0.5 * Globals::Pi * SWIMMING_POW_2(particle_radius) *  fluid_density * coeff * norm_of_slip_vel / vorticity_norm;
            return(lift_coeff);
        }

        else {
            return 0.0;
        }
    }

    void ElSamniLiftLaw::ComputeForce(Geometry<Node >& r_geometry,
                                      const double reynolds_number,
                                      double particle_radius,
                                      double fluid_density,
                                      double fluid_kinematic_viscosity,
                                      array_1d<double, 3>& minus_slip_velocity,
                                      array_1d<double, 3>& vorticity_induced_lift,
                                      const ProcessInfo& r_current_process_info)
    {
        Node& node = r_geometry[0];
        const array_1d<double, 3>& vorticity = node.FastGetSolutionStepValue(FLUID_VORTICITY_PROJECTED);
        array_1d<double, 3> vort_cross_slip_vel;
        SWIMMING_SET_TO_CROSS_OF_FIRST_TWO_3(minus_slip_velocity, vorticity, vort_cross_slip_vel)
        const double norm_of_slip_vel = SWIMMING_MODULUS_3(minus_slip_velocity);
        const double vorticity_norm = SWIMMING_MODULUS_3(vorticity);

        const double lift_coeff = ComputeElSamniLiftCoefficient(particle_radius,
                                                                fluid_density,
                                                                norm_of_slip_vel,
                                                                vorticity_norm,
                                                                r_current_process_info);

        noalias(vorticity_induced_lift) = lift_coeff * vort_cross_slip_vel; // the direction is given by the vorticity x (- slip_vel) (Jackson, 2000), which is normalized here
    }

} // namespace Kratos
