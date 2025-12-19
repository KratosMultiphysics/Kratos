// Author: Guillermo Casas (gcasas@cimne.upc.edu)
// Date: February 2019
#include "swimming_DEM_application.h"
#include "ganser_drag_law.h"

namespace Kratos {

    DragLaw::Pointer GanserDragLaw::Clone() const {
        GanserDragLaw::Pointer p_clone(new GanserDragLaw(*this));
        return p_clone;
    }

    void GanserDragLaw::Initialize(const ProcessInfo& r_process_info) {}

    std::string GanserDragLaw::GetTypeOfLaw() {
        std::string type_of_law = "GanserDragLaw";
        return type_of_law;
    }

    void GanserDragLaw::ComputeForce(SphericParticle* p_particle,
                                       const double reynolds_number,
                                       double particle_radius,
                                       double fluid_density,
                                       double fluid_kinematic_viscosity,
                                       array_1d<double, 3>& minus_slip_velocity,
                                       array_1d<double, 3>& drag_force,
                                       const ProcessInfo& r_current_process_info)
    {
        const int isometric_shape                = 1; // TEMPORARY!! yes (1) or no (0); should be given as data
        const double surface_area                = 4 * Globals::Pi * SWIMMING_POW_2(particle_radius); // TEMPORARY!! corresponding to a sphere; should be generalized b taking it as a parameter
        const double surface_area_circular_diam  = std::sqrt(4.0 * surface_area / Globals::Pi);
        Geometry<Node >& r_geometry = p_particle->GetGeometry();
        const double sphericity = r_geometry[0].FastGetSolutionStepValue(PARTICLE_SPHERICITY);
        double k_1;
        double k_2;

        // calculating ganser geometric parameters and equivalent Reynolds number
        if (isometric_shape){
            k_1 = 3 / (1 + 2 / std::sqrt(sphericity));
        }

        else {
            k_1 = 3 / (0.5 * surface_area_circular_diam / particle_radius + 2 / std::sqrt(sphericity));
        }

        k_2 = std::pow(10.0, 1.8148 * std::pow(- std::log10(sphericity), 0.5743));

        double equiv_reynolds = k_1 * k_2 * reynolds_number;

        // calculating nondimensional drag coefficient
        double drag_coeff =  k_2 * (24 * (1 + 0.1118 * std::pow((equiv_reynolds), 0.6567)) / (equiv_reynolds) + 0.4305 / (1 + 3305 / equiv_reynolds));

        // and then the dimensional drag coefficient
        drag_coeff *= 0.5 *  fluid_density * surface_area * SWIMMING_MODULUS_3(minus_slip_velocity);

        noalias(drag_force) = drag_coeff * minus_slip_velocity;
    }
} // namespace Kratos
