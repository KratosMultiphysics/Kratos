// Author: Guillermo Casas (gcasas@cimne.upc.edu)
// Date: February 2019

#include "swimming_DEM_application.h"
#include "zuber_virtual_mass_force_law.h"

namespace Kratos {

    VirtualMassForceLaw::Pointer ZuberVirtualMassForceLaw::Clone() const {
        ZuberVirtualMassForceLaw::Pointer p_clone(new ZuberVirtualMassForceLaw(*this));
        return p_clone;
    }

    ZuberVirtualMassForceLaw::ZuberVirtualMassForceLaw(Parameters r_parameters)
        : AutonHuntPrudhommeVirtualMassForceLaw::AutonHuntPrudhommeVirtualMassForceLaw(r_parameters)
    {}

    void ZuberVirtualMassForceLaw::Initialize(const ProcessInfo& r_process_info){}

    std::string ZuberVirtualMassForceLaw::GetTypeOfLaw() {
        std::string type_of_law = "Zuber virtual mass force law";
        return type_of_law;
    }

    double ZuberVirtualMassForceLaw::GetVirtualMassCoefficient(Geometry<Node >& r_geometry,
                                                            const array_1d<double, 3>& minus_slip_acc)
    {
        Node& node = r_geometry[0];
        const double fluid_volume_fraction = node.FastGetSolutionStepValue(FLUID_FRACTION_PROJECTED);
        double standard_virtual_mass_coeff = AutonHuntPrudhommeVirtualMassForceLaw::GetVirtualMassCoefficient(r_geometry, minus_slip_acc);
        return standard_virtual_mass_coeff + 1.5 * (1 - fluid_volume_fraction);
    }

} // namespace Kratos
