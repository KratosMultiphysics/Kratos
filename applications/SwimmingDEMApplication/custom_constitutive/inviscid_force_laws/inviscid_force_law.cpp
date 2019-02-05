// Author: Guillermo Casas (gcasas@cimne.upc.edu)
// Date: February 2019

#include "inviscid_force_law.h"

namespace Kratos {

    void InviscidForceLaw::Initialize(const ProcessInfo& r_process_info) {}

    std::string InviscidForceLaw::GetTypeOfLaw() {
        std::string type_of_law = "Generic Drag Law";
        return type_of_law;
    }

    void InviscidForceLaw::SetInviscidForceLawInProperties(Properties::Pointer pProp) const {
        pProp->SetValue(SDEM_INVISCID_FORCE_LAW_POINTER, this->Clone());
    }

    InviscidForceLaw::Pointer InviscidForceLaw::Clone() const {
        InviscidForceLaw::Pointer p_clone(new InviscidForceLaw(*this));
        return p_clone;
    }

} // namespace Kratos
