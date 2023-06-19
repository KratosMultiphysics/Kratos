// Author: Guillermo Casas (gcasas@cimne.upc.edu)
// Date: February 2019
#include "swimming_DEM_application.h"
#include "shah_drag_law.h"

namespace Kratos {

    DragLaw::Pointer ShahDragLaw::Clone() const {
        ShahDragLaw::Pointer p_clone(new ShahDragLaw(*this));
        return p_clone;
    }

    void ShahDragLaw::Initialize(const ProcessInfo& r_process_info) {}

    std::string ShahDragLaw::GetTypeOfLaw() {
        std::string type_of_law = "ShahDragLaw";
        return type_of_law;
    }

    void ShahDragLaw::ComputeForce(SphericParticle* p_particle,
                                       const double reynolds_number,
                                       double particle_radius,
                                       double fluid_density,
                                       double fluid_kinematic_viscosity,
                                       array_1d<double, 3>& minus_slip_velocity,
                                       array_1d<double, 3>& drag_force,
                                       const ProcessInfo& r_current_process_info)
    {
    const double power_law_tol = 0.0001;
    const double K = r_current_process_info[POWER_LAW_K];
    const double n = r_current_process_info[POWER_LAW_N];
    const bool use_shahi_correction = false; //TODO: make specific law for this option

    if (std::abs(n) < power_law_tol || std::abs(K) < power_law_tol){
        std::cout << "WARNING: Shah's method is being used with Power Law data being zero (n = 0 or K = 0)!!" << std::endl << std::flush;
    }

    double A =   6.9148 * n * n - 24.838 * n + 22.642;
    double B = - 0.5067 * n * n + 1.3234 * n - 0.1744;

    if (use_shahi_correction){ // from 2016 Shahi (doi: 10.1016/j.minpro.2016.06.002)
        A = 1.5269 * A - 3.9375;
        B =  0.892 * B + 0.0326;
    }

    const double exponents_coeff = 1.0 / (2 - n);
    const double area = Globals::Pi * particle_radius * particle_radius;
    const double dimensional_coefficient = 0.5 * area * fluid_density * SWIMMING_MODULUS_3(minus_slip_velocity);

    const double drag_coeff = dimensional_coefficient * std::pow(A, exponents_coeff) * std::pow(reynolds_number, exponents_coeff * (2 * B - 2));
    noalias(drag_force) = drag_coeff * minus_slip_velocity;
    }
} // namespace Kratos
