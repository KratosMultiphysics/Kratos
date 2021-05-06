// Author: Guillermo Casas (gcasas@cimne.upc.edu)
// Date: February 2019
#include "swimming_DEM_application.h"
#include "schiller_and_naumann_drag_law.h"

namespace Kratos {

    DragLaw::Pointer SchillerAndNaumannDragLaw::Clone() const {
        SchillerAndNaumannDragLaw::Pointer p_clone(new SchillerAndNaumannDragLaw(*this));
        return p_clone;
    }

    void SchillerAndNaumannDragLaw::Initialize(const ProcessInfo& r_process_info) {}

    std::string SchillerAndNaumannDragLaw::GetTypeOfLaw() {
        std::string type_of_law = "Schiller and Naumann drag law";
        return type_of_law;
    }

    void SchillerAndNaumannDragLaw::ComputeForce(SphericParticle* p_particle,
                                       const double reynolds_number,
                                       double particle_radius,
                                       double fluid_density,
                                       double fluid_kinematic_viscosity,
                                       array_1d<double, 3>& minus_slip_velocity,
                                       array_1d<double, 3>& drag_force,
                                       const ProcessInfo& r_current_process_info)
    {
        StokesDragLaw::ComputeForce(p_particle,
                                    reynolds_number,
                                    particle_radius,
                                    fluid_density,
                                    fluid_kinematic_viscosity,
                                    minus_slip_velocity,
                                    drag_force,
                                    r_current_process_info);

        if (reynolds_number < 1000){
            drag_force *= (1 + 0.15 * std::pow(reynolds_number, 0.687));
        }

        else {
            drag_force *= 0.01826 * reynolds_number;
        }
    }
} // namespace Kratos
