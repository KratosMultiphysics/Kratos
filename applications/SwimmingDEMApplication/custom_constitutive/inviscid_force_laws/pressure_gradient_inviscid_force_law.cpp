// Author: Joaquín González-Usúa (jgonzalez@cimne.upc.edu)
// Date: January 2024

#include "swimming_DEM_application.h"
#include "pressure_gradient_inviscid_force_law.h"

namespace Kratos {

    InviscidForceLaw::Pointer PressureGradientInviscidForceLaw::Clone() const {
        PressureGradientInviscidForceLaw::Pointer p_clone(new PressureGradientInviscidForceLaw(*this));
        return p_clone;
    }

    PressureGradientInviscidForceLaw::PressureGradientInviscidForceLaw(Parameters r_parameters)
    {}

    void PressureGradientInviscidForceLaw::Initialize(const ProcessInfo& r_process_info){}

    std::string PressureGradientInviscidForceLaw::GetTypeOfLaw() {
        std::string type_of_law = "Pressure Gradient inviscid force law";
        return type_of_law;
    }

    void PressureGradientInviscidForceLaw::ComputeForce(Geometry<Node >& r_geometry,
                                                const double fluid_density,
                                                const double displaced_volume,
                                                array_1d<double, 3>& virtual_mass_plus_undisturbed_flow_force,
                                                const ProcessInfo& r_current_process_info)
    {
        const array_1d<double, 3>& fluid_press_grad = r_geometry[0].FastGetSolutionStepValue(PRESSURE_GRAD_PROJECTED);
        // the particle acceleration is assumed to be treated implicitly through the added_mass

        noalias(virtual_mass_plus_undisturbed_flow_force) = - displaced_volume * fluid_press_grad;
    }

} // namespace Kratos
