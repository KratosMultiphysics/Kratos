// Author: Guillermo Casas (gcasas@cimne.upc.edu)
// Date: February 2019

#include "swimming_DEM_application.h"
#include "standard_buoyancy_law.h"

namespace Kratos {

    BuoyancyLaw::Pointer StandardBuoyancyLaw::Clone() const {
        StandardBuoyancyLaw::Pointer p_clone(new StandardBuoyancyLaw(*this));
        return p_clone;
    }

    void StandardBuoyancyLaw::Initialize(const ProcessInfo& r_process_info){}

    std::string StandardBuoyancyLaw::GetTypeOfLaw() {
        std::string type_of_law = "Standard Buoyancy Force Law";
        return type_of_law;
    }

    void StandardBuoyancyLaw::ComputeForce(Geometry<Node<3> >& r_geometry,
                                           const double fluid_density,
                                           const double displaced_volume,
                                           const array_1d<double, 3>& body_force,
                                           array_1d<double, 3>& buoyancy,
                                           const ProcessInfo& r_current_process_info)
    {
        const double fluid_mass = displaced_volume * fluid_density;
        noalias(buoyancy) = - body_force * fluid_mass;
    }

} // namespace Kratos
