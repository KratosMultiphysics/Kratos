#include "hydrodynamic_interaction_law.h"

namespace Kratos {

    HydrodynamicInteractionLaw::HydrodynamicInteractionLaw() {
        auto mpDragLaw = make_unique<DragLaw>();
    }

    HydrodynamicInteractionLaw::HydrodynamicInteractionLaw(const HydrodynamicInteractionLaw &rHydrodynamicInteractionLaw) {
    }

    void HydrodynamicInteractionLaw::Initialize(const ProcessInfo& r_process_info) {
    }

    void HydrodynamicInteractionLaw::SetHydrodynamicInteractionLawInProperties(Properties::Pointer pProp, bool verbose) const {
        pProp->SetValue(SDEM_HYDRODYNAMIC_INTERACTION_LAW_POINTER, this->Clone());
    }

    std::string HydrodynamicInteractionLaw::GetTypeOfLaw() {
        std::string type_of_law = "HydrodynamicInteractionLaw";
        return type_of_law;
    }

    HydrodynamicInteractionLaw::Pointer HydrodynamicInteractionLaw::Clone() const {
        HydrodynamicInteractionLaw::Pointer p_clone(new HydrodynamicInteractionLaw(*this));
        return p_clone;
    }

    HydrodynamicInteractionLaw::~HydrodynamicInteractionLaw(){}


    void HydrodynamicInteractionLaw::ComputeDragForce(NodeType& node,
                                                      double particle_radius,
                                                      double fluid_density,
                                                      double fluid_kinematic_viscosity,
                                                      array_1d<double, 3>& slip_velocity,
                                                      array_1d<double, 3>& drag_force,
                                                      const ProcessInfo& r_current_process_info)
    {
        mpDragLaw->ComputeForce(node,
                                particle_radius,
                                fluid_density,
                                fluid_kinematic_viscosity,
                                slip_velocity,
                                drag_force,
                                r_current_process_info);

    }

} // Namespace Kratos
