// Author: Guillermo Casas (gcasas@cimne.upc.edu)
// Date: February 2019

#include "swimming_DEM_application.h"
#include "mei_lift_law.h"

namespace Kratos {

    VorticityInducedLiftLaw::Pointer MeiLiftLaw::Clone() const {
        MeiLiftLaw::Pointer p_clone(new MeiLiftLaw(*this));
        return p_clone;
    }

    MeiLiftLaw::MeiLiftLaw(Parameters r_parameters)
    {

    }

    void MeiLiftLaw::Initialize(const ProcessInfo& r_process_info){}

    std::string MeiLiftLaw::GetTypeOfLaw() {
        std::string type_of_law = "Mei lift law";
        return type_of_law;
    }

    double MeiLiftLaw::ComputeShearReynoldsNumber(const double particle_radius,
                                                  const double fluid_kinematic_viscosity,
                                                  const double norm_of_vorticity)
    {
        return 4 * norm_of_vorticity * particle_radius * particle_radius / fluid_kinematic_viscosity;
    }

    double MeiLiftLaw::ComputeMeiCorrectionOnSaffmanCoefficient(const double reynolds_number,
                                                                const double fluid_kinematic_viscosity,
                                                                const double particle_radius,
                                                                const double norm_of_vorticity)
    {
        const double reynolds_shear = ComputeShearReynoldsNumber(particle_radius,
                                                                 fluid_kinematic_viscosity,
                                                                 norm_of_vorticity);
        if (reynolds_number != 0.0 && reynolds_shear != 0.0 ){
            const double alpha = 0.5 * reynolds_shear / reynolds_number;
            double mei_over_saff;

            if (reynolds_number < 40){
                const double sqrt_alpha = std::sqrt(alpha);
                mei_over_saff = (1 - 0.3314 * sqrt_alpha) * std::exp(- 0.1 * reynolds_number) + 0.3314 * sqrt_alpha;
            }

            else {
                mei_over_saff = 0.0524 * std::sqrt(alpha * reynolds_number);
            }

            return mei_over_saff;
        }

        else {
            return 0.0;
        }
    }

    void MeiLiftLaw::ComputeForce(Geometry<Node >& r_geometry,
                                      const double reynolds_number,
                                      double particle_radius,
                                      double fluid_density,
                                      double fluid_kinematic_viscosity,
                                      array_1d<double, 3>& minus_slip_velocity,
                                      array_1d<double, 3>& vorticity_induced_lift,
                                      const ProcessInfo& r_current_process_info)
    {
        // First compute Saffman's lift
        SaffmanLiftLaw::ComputeForce(r_geometry,
                                     reynolds_number,
                                     particle_radius,
                                     fluid_density,
                                     fluid_kinematic_viscosity,
                                     minus_slip_velocity,
                                     vorticity_induced_lift,
                                     r_current_process_info);

        // Then apply Mei's correction
        Node& node = r_geometry[0];
        const array_1d<double, 3>& vorticity = node.FastGetSolutionStepValue(FLUID_VORTICITY_PROJECTED);
        array_1d<double, 3> vort_cross_slip_vel;
        SWIMMING_SET_TO_CROSS_OF_FIRST_TWO_3(minus_slip_velocity, vorticity, vort_cross_slip_vel)
        const double norm_of_vorticity = SWIMMING_MODULUS_3(vorticity);

        const double mei_over_saff = ComputeMeiCorrectionOnSaffmanCoefficient(reynolds_number,
                                                                              fluid_kinematic_viscosity,
                                                                              particle_radius,
                                                                              norm_of_vorticity);

        vorticity_induced_lift *= mei_over_saff;
    }

} // namespace Kratos
