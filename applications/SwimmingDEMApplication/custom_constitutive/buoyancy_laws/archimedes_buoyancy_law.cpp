// Author: Guillermo Casas (gcasas@cimne.upc.edu)
// Date: February 2019

#include "swimming_DEM_application.h"
#include "archimedes_buoyancy_law.h"

namespace Kratos {

    BuoyancyLaw::Pointer ArchimedesBuoyancyLaw::Clone() const {
        ArchimedesBuoyancyLaw::Pointer p_clone(new ArchimedesBuoyancyLaw(*this));
        return p_clone;
    }

    void ArchimedesBuoyancyLaw::Initialize(const ProcessInfo& r_process_info){}

    std::string ArchimedesBuoyancyLaw::GetTypeOfLaw() {
        std::string type_of_law = "ArchimedesBuoyancyLaw";
        return type_of_law;
    }

    void ArchimedesBuoyancyLaw::ComputeForce(Geometry<Node >& r_geometry,
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
