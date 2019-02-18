// Author: Guillermo Casas (gcasas@cimne.upc.edu)
// Date: February 2019

#include "rotation_induced_lift_law.h"

namespace Kratos {

    void RotationInducedLiftLaw::Initialize(const ProcessInfo& r_process_info) {}

    std::string RotationInducedLiftLaw::GetTypeOfLaw() {
        std::string type_of_law = "Generic rotation-induced lift law";
        return type_of_law;
    }

    void RotationInducedLiftLaw::SetRotationInducedLiftLawInProperties(Properties::Pointer pProp) const {
    }

    RotationInducedLiftLaw::Pointer RotationInducedLiftLaw::Clone() const {
        RotationInducedLiftLaw::Pointer p_clone(new RotationInducedLiftLaw(*this));
        return p_clone;
    }

} // namespace Kratos
