//
// Author: Salva Latorre    latorre@cimne.upc.edu
//

// System includes
#include <string>
#include <iostream>
#include <iomanip> // to improve std::cout precision

// Project includes
#include "ice_continuum_particle.h"

namespace Kratos {
    
    void IceContinuumParticle::Initialize(const ProcessInfo& r_process_info) {   
        SphericContinuumParticle::Initialize(r_process_info);
        double added_mass_coefficient = 1.0;
        SetMass(added_mass_coefficient * GetDensity() * CalculateVolume());
    }

    void IceContinuumParticle::MemberDeclarationFirstStep(const ProcessInfo& r_process_info) {
        SphericContinuumParticle::MemberDeclarationFirstStep(r_process_info);
    }

    double IceContinuumParticle::CalculateVolume()
    {
        //
        double representative_volume = 1.0; // Compute properly!!!
        //
        return representative_volume;
    }

    void IceContinuumParticle::SetInteractionRadius(const double radius){}

    double IceContinuumParticle::GetInteractionRadius()
    {
        return mSearchRadius;
    }
    void IceContinuumParticle::ComputeAdditionalForces(array_1d<double, 3>& externally_applied_force,
                                                  array_1d<double, 3>& externally_applied_moment,
                                                  const ProcessInfo& r_process_info,
                                                  const array_1d<double,3>& gravity)
    {
        KRATOS_TRY

        SphericParticle::ComputeAdditionalForces(externally_applied_force,
                                                         externally_applied_moment,
                                                         r_process_info,
                                                         gravity);

        double Z_coord = this->GetGeometry()[0].Coordinates()[2];
        double density_water = 1000.0;
        
        //Adding buoyancy:
        if (Z_coord < 1.0) noalias(externally_applied_force) -= CalculateVolume() * density_water * gravity;

        KRATOS_CATCH("")
    }
    
} // namespace Kratos
