// Author: Guillermo Casas (gcasas@cimne.upc.edu)
// Date: February 2019
#include "swimming_DEM_application.h"
#include "beetstra_drag_law.h"

namespace Kratos {

    DragLaw::Pointer BeetstraDragLaw::Clone() const {
        BeetstraDragLaw::Pointer p_clone(new BeetstraDragLaw(*this));
        return p_clone;
    }

    void BeetstraDragLaw::Initialize(const ProcessInfo& r_process_info) {}

    std::string BeetstraDragLaw::GetTypeOfLaw() {
        std::string type_of_law = "Beetstra drag law";
        return type_of_law;
    }

    void BeetstraDragLaw::ComputeForce(SphericParticle* p_particle,
                                       const double reynolds_number,
                                       double particle_radius,
                                       double fluid_density,
                                       double fluid_kinematic_viscosity,
                                       array_1d<double, 3>& minus_slip_velocity,
                                       array_1d<double, 3>& drag_force,
                                       const ProcessInfo& r_current_process_info)
    {

        if (reynolds_number < 1.0){
            return StokesDragLaw::ComputeForce(p_particle,
                                               reynolds_number,
                                               particle_radius,
                                               fluid_density,
                                               fluid_kinematic_viscosity,
                                               minus_slip_velocity,
                                               drag_force,
                                               r_current_process_info);
        }

        double drag_coeff;
        Geometry<Node >& r_geometry = p_particle->GetGeometry();
        double eps = r_geometry[0].FastGetSolutionStepValue(FLUID_FRACTION_PROJECTED);

        if (eps > 0.999){
            eps = 0.9;
        }

        const double eps_s = 1.0 - eps;
        const double mod_reynolds_number = reynolds_number * eps;

        double A = 180 + 18 * std::pow(eps, 4) / eps_s * (1 + 1.5 * std::sqrt(eps_s));
        double B = 0.31 * (1.0 / eps + 3 * eps_s * eps + 8.4 * std::pow(mod_reynolds_number, - 0.343)) / (1.0 + std::pow(10.0, 3 * eps_s) * std::pow(mod_reynolds_number, 2 * eps - 2.5));
        drag_coeff = Globals::Pi / 3.0 * fluid_kinematic_viscosity * fluid_density * particle_radius * (A * eps_s / eps + B * mod_reynolds_number);

        noalias(drag_force) = drag_coeff * minus_slip_velocity;
    }
} // namespace Kratos
