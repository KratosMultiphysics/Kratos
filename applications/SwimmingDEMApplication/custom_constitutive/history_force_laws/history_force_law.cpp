// Author: Guillermo Casas (gcasas@cimne.upc.edu)
// Date: February 2019

#include "history_force_law.h"

namespace Kratos {

    void HistoryForceLaw::Initialize(const ProcessInfo& r_process_info) {}

    std::string HistoryForceLaw::GetTypeOfLaw() {
        std::string type_of_law = "Generic history force law";
        return type_of_law;
    }

    void HistoryForceLaw::SetHistoryForceLawInProperties(Properties::Pointer pProp) const {
        pProp->SetValue(SDEM_HISTORY_FORCE_LAW_POINTER, this->Clone());
    }

    HistoryForceLaw::Pointer HistoryForceLaw::Clone() const {
        HistoryForceLaw::Pointer p_clone(new HistoryForceLaw(*this));
        return p_clone;
    }

} // namespace Kratos
