// Author: Guillermo Casas (gcasas@cimne.upc.edu)
// Date: February 2019

#include "vorticity_induced_lift_law.h"

namespace Kratos {

    void VorticityInducedLiftLaw::Initialize(const ProcessInfo& r_process_info) {}

    std::string VorticityInducedLiftLaw::GetTypeOfLaw() {
        std::string type_of_law = "Generic vorticity-induced lift law";
        return type_of_law;
    }

    void VorticityInducedLiftLaw::SetVorticityInducedLiftLawInProperties(Properties::Pointer pProp) const {
    }

    VorticityInducedLiftLaw::Pointer VorticityInducedLiftLaw::Clone() const {
        VorticityInducedLiftLaw::Pointer p_clone(new VorticityInducedLiftLaw(*this));
        return p_clone;
    }

} // namespace Kratos
