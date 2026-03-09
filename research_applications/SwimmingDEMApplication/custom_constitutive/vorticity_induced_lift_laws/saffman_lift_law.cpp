// Author: Guillermo Casas (gcasas@cimne.upc.edu)
// Date: February 2019

#include "swimming_DEM_application.h"
#include "saffman_lift_law.h"

namespace Kratos {

    VorticityInducedLiftLaw::Pointer SaffmanLiftLaw::Clone() const {
        SaffmanLiftLaw::Pointer p_clone(new SaffmanLiftLaw(*this));
        return p_clone;
    }

    SaffmanLiftLaw::SaffmanLiftLaw(Parameters r_parameters)
    {

    }

    void SaffmanLiftLaw::Initialize(const ProcessInfo& r_process_info){}

    std::string SaffmanLiftLaw::GetTypeOfLaw() {
        std::string type_of_law = "Saffman lift law";
        return type_of_law;
    }

    double SaffmanLiftLaw::ComputeSaffmanLiftCoefficient(const double fluid_density,
                                                         const double fluid_kinematic_viscosity,
                                                         const double particle_radius,
                                                         const double norm_of_vorticity)
    {
        if (norm_of_vorticity != 0.0 ){
            return 6.46 * fluid_density * particle_radius * particle_radius * std::sqrt(fluid_kinematic_viscosity / norm_of_vorticity);
        }

        else {
            return 0.0;
        }
    }

    void SaffmanLiftLaw::ComputeForce(Geometry<Node >& r_geometry,
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
        const double norm_of_vorticity = SWIMMING_MODULUS_3(vorticity);
        const double lift_coeff = ComputeSaffmanLiftCoefficient(fluid_density, fluid_kinematic_viscosity, particle_radius, norm_of_vorticity);

        noalias(vorticity_induced_lift) = lift_coeff * vort_cross_slip_vel; // the direction is given by the vorticity x (- slip_vel) (Jackson, 2000), which is normalized here
    }

} // namespace Kratos
