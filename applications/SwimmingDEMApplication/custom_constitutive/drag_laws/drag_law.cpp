// Author: Guillermo Casas (gcasas@cimne.upc.edu)
// Date: February 2019

#include "drag_law.h"

namespace Kratos {

    void DragLaw::Initialize(const ProcessInfo& r_process_info) {}

    std::string DragLaw::GetTypeOfLaw() {
        std::string type_of_law = "Generic drag law";
        return type_of_law;
    }


    void DragLaw::SetDragLawInProperties(Properties::Pointer pProp) const {
        pProp->SetValue(SDEM_DRAG_LAW_POINTER, this->Clone());
    }

    void DragLaw::ComputeForce(SphericParticle* p_particle,
                               const double reynolds_number,
                               double particle_radius,
                               double fluid_density,
                               double fluid_kinematic_viscosity,
                               array_1d<double, 3>& minus_slip_velocity,
                               array_1d<double, 3>& drag_force,
                               const ProcessInfo& r_current_process_info){}

    DragLaw::Pointer DragLaw::Clone() const {
        DragLaw::Pointer p_clone(new DragLaw(*this));
        return p_clone;
    }
} // namespace Kratos
