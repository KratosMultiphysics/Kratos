// Author: Guillermo Casas (gcasas@cimne.upc.edu)
// Date: February 2019
#include "swimming_DEM_application.h"
#include "schiller_and_naumann_drag_law.h"

namespace Kratos {

    DragLaw::Pointer SchillerAndNaumann::Clone() const {
        SchillerAndNaumann::Pointer p_clone(new SchillerAndNaumann(*this));
        return p_clone;
    }

    void SchillerAndNaumann::Initialize(const ProcessInfo& r_process_info) {}

    std::string SchillerAndNaumann::GetTypeOfLaw() {
        std::string type_of_law = "SchillerAndNaumann";
        return type_of_law;
    }

    void SchillerAndNaumann::ComputeForce(Geometry<Node<3> >& r_geometry,
                                       const double reynolds_number,
                                       double particle_radius,
                                       double fluid_density,
                                       double fluid_kinematic_viscosity,
                                       array_1d<double, 3>& slip_velocity,
                                       array_1d<double, 3>& drag_force,
                                       const ProcessInfo& r_current_process_info)
    {
        double drag_coeff  = 0.5 * Globals::Pi * SWIMMING_POW_2(particle_radius) * fluid_density * SWIMMING_MODULUS_3(slip_velocity);

        drag_coeff *= 24 / reynolds_number * (1 + 0.15 * pow(reynolds_number, 0.687));
        noalias(drag_force) = drag_coeff * slip_velocity;
    }
} // namespace Kratos
