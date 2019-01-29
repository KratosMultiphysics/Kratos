#include "hydrodynamic_interaction_law.h"

namespace Kratos {

    HydrodynamicInteractionLaw::HydrodynamicInteractionLaw() {

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

} // Namespace Kratos
