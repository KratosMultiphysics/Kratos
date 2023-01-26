// Author: Guillermo Casas (gcasas@cimne.upc.edu)
// Date: February 2019
#include "swimming_DEM_application.h"
#include "haider_and_levenspiel_drag_law.h"

namespace Kratos {

    DragLaw::Pointer HaiderAndLevenspielDragLaw::Clone() const {
        HaiderAndLevenspielDragLaw::Pointer p_clone(new HaiderAndLevenspielDragLaw(*this));
        return p_clone;
    }

    void HaiderAndLevenspielDragLaw::Initialize(const ProcessInfo& r_process_info) {}

    std::string HaiderAndLevenspielDragLaw::GetTypeOfLaw() {
        std::string type_of_law = "HaiderAndLevenspielDragLaw";
        return type_of_law;
    }

    void HaiderAndLevenspielDragLaw::ComputeForce(SphericParticle* p_particle,
                                       const double reynolds_number,
                                       double particle_radius,
                                       double fluid_density,
                                       double fluid_kinematic_viscosity,
                                       array_1d<double, 3>& minus_slip_velocity,
                                       array_1d<double, 3>& drag_force,
                                       const ProcessInfo& r_current_process_info)
    {
        Geometry<Node<3> >& r_geometry = p_particle->GetGeometry();
        const double sphericity = r_geometry[0].FastGetSolutionStepValue(PARTICLE_SPHERICITY);
        double drag_coeff       = 0.5 * Globals::Pi * SWIMMING_POW_2(particle_radius) * fluid_density * SWIMMING_MODULUS_3(minus_slip_velocity);

        double A = exp(2.3288 - 6.4581 * sphericity + 2.4486 * sphericity * sphericity);
        double B = 0.0964 + 0.5565 * sphericity;
        double C = exp(4.905  - 13.8944 * sphericity + 18.4222 * sphericity * sphericity - 10.2599 * sphericity * sphericity * sphericity);
        double D = exp(1.4681 + 12.2584 * sphericity - 20.7322 * sphericity * sphericity + 15.8855 * sphericity * sphericity * sphericity);

        drag_coeff *= (24.0 * (1.0 + A * pow(reynolds_number, B))) / reynolds_number + C * reynolds_number / (reynolds_number + D);

        noalias(drag_force) = drag_coeff * minus_slip_velocity;
    }
} // namespace Kratos
