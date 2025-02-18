// Author: Joaquin Gonzalez-Usua (jgonzalez@cimne.upc.edu)
// Date: January 2024
#include "swimming_DEM_application.h"
#include "undisturbed_force_law.h"

namespace Kratos {

    UndisturbedForceLaw::Pointer UndisturbedForceLaw::Clone() const {
        UndisturbedForceLaw::Pointer p_clone(new UndisturbedForceLaw(*this));
        return p_clone;
    }

    UndisturbedForceLaw::UndisturbedForceLaw(Parameters r_parameters)
    {
        Parameters default_parameters( R"(
            {
                "name":"UndisturbedForceLaw"
            }  )" );

        r_parameters.ValidateAndAssignDefaults(default_parameters);
    }

    void UndisturbedForceLaw::Initialize(const ProcessInfo& r_process_info){}

    std::string UndisturbedForceLaw::GetTypeOfLaw() {
        std::string type_of_law = "Undisturbed force law";
        return type_of_law;
    }

    void UndisturbedForceLaw::ComputeForce(Geometry<Node >& r_geometry,
                                                const double fluid_density,
                                                const double displaced_volume,
                                                array_1d<double, 3>& undisturbed_flow_force,
                                                const ProcessInfo& r_current_process_info)
    {
        const array_1d<double, 3>& fluid_acc = r_geometry[0].FastGetSolutionStepValue(FLUID_ACCEL_PROJECTED);
        // const double radius = r_geometry[0].FastGetSolutionStepValue(RADIUS);

        const double fluid_mass = displaced_volume * fluid_density;

        noalias(undisturbed_flow_force) = fluid_mass * fluid_acc;
    }

} // namespace Kratos