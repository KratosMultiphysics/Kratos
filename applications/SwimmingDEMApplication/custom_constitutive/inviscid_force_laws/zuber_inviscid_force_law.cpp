// Author: Guillermo Casas (gcasas@cimne.upc.edu)
// Date: February 2019

#include "swimming_DEM_application.h"
#include "zuber_inviscid_force_law.h"

namespace Kratos {

    InviscidForceLaw::Pointer ZuberInviscidForceLaw::Clone() const {
        ZuberInviscidForceLaw::Pointer p_clone(new ZuberInviscidForceLaw(*this));
        return p_clone;
    }

    ZuberInviscidForceLaw::ZuberInviscidForceLaw(Parameters r_parameters)
        : AutonHuntPrudhommeInviscidForceLaw::AutonHuntPrudhommeInviscidForceLaw(r_parameters)
    {}

    void ZuberInviscidForceLaw::Initialize(const ProcessInfo& r_process_info){}

    std::string ZuberInviscidForceLaw::GetTypeOfLaw() {
        std::string type_of_law = "Zuber inviscid force law";
        return type_of_law;
    }

    double ZuberInviscidForceLaw::GetVirtualMassCoefficient(Geometry<Node >& r_geometry,
                                                            const array_1d<double, 3>& minus_slip_acc)
    {
        Node& node = r_geometry[0];
        const double fluid_volume_fraction = node.FastGetSolutionStepValue(FLUID_FRACTION_PROJECTED);
        double standard_virtual_mass_coeff = AutonHuntPrudhommeInviscidForceLaw::GetVirtualMassCoefficient(r_geometry, minus_slip_acc);
        return standard_virtual_mass_coeff + 1.5 * (1 - fluid_volume_fraction);
    }

} // namespace Kratos
