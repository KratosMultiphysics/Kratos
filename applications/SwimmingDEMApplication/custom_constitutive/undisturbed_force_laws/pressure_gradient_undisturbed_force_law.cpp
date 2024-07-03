// Author: Joaquín González-Usúa (jgonzalez@cimne.upc.edu)
// Date: January 2024

#include "swimming_DEM_application.h"
#include "pressure_gradient_undisturbed_force_law.h"

namespace Kratos {

    UndisturbedForceLaw::Pointer PressureGradientUndisturbedForceLaw::Clone() const {
        PressureGradientUndisturbedForceLaw::Pointer p_clone(new PressureGradientUndisturbedForceLaw(*this));
        return p_clone;
    }

    PressureGradientUndisturbedForceLaw::PressureGradientUndisturbedForceLaw(Parameters r_parameters)
    {

    }

    void PressureGradientUndisturbedForceLaw::Initialize(const ProcessInfo& r_process_info){}

    std::string PressureGradientUndisturbedForceLaw::GetTypeOfLaw() {
        std::string type_of_law = "Pressure Gradient undisturbed force law";
        return type_of_law;
    }

    void PressureGradientUndisturbedForceLaw::ComputeForce(Geometry<Node >& r_geometry,
                                                const double fluid_density,
                                                const double displaced_volume,
                                                array_1d<double, 3>& undisturbed_flow_force,
                                                const ProcessInfo& r_current_process_info)
    {
        const array_1d<double, 3>& fluid_press_grad = r_geometry[0].FastGetSolutionStepValue(PRESSURE_GRAD_PROJECTED);
        // the particle acceleration is assumed to be treated implicitly through the added_mass
        noalias(undisturbed_flow_force) = -displaced_volume * fluid_press_grad;
    }

} // namespace Kratos
