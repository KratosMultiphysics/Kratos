// Author: Guillermo Casas (gcasas@cimne.upc.edu)
// Date: February 2019

#include "buoyancy_law.h"

namespace Kratos {

    void BuoyancyLaw::Initialize(const ProcessInfo& r_process_info) {}

    std::string BuoyancyLaw::GetTypeOfLaw() {
        std::string type_of_law = "Generic buoyancy law";
        return type_of_law;
    }

    void BuoyancyLaw::SetBuoyancyLawInProperties(Properties::Pointer pProp) const {
        pProp->SetValue(SDEM_BUOYANCY_LAW_POINTER, this->Clone());
    }

    BuoyancyLaw::Pointer BuoyancyLaw::Clone() const {
        BuoyancyLaw::Pointer p_clone(new BuoyancyLaw(*this));
        return p_clone;
    }

} // namespace Kratos
