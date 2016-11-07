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

        SetMass(GetDensity() * CalculateVolume());
    }

    void IceContinuumParticle::MemberDeclarationFirstStep(const ProcessInfo& r_process_info) {
        SphericContinuumParticle::MemberDeclarationFirstStep(r_process_info);
    }

    double IceContinuumParticle::CalculateVolume()
    {
        return 4.0 * KRATOS_M_PI_3 * mRadius * mRadius * mRadius;
    }

    void IceContinuumParticle::SetInteractionRadius(const double radius) {}

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

        SphericContinuumParticle::ComputeAdditionalForces(externally_applied_force,
                                                         externally_applied_moment,
                                                         r_process_info,
                                                         gravity);
        
        ComputeIceForces(this, externally_applied_force, gravity);

        KRATOS_CATCH("")
    }
    
    void IceContinuumParticle::ComputeIceForces(SphericContinuumParticle* continuum_particle,
                                                array_1d<double, 3>& externally_applied_force,
                                                const array_1d<double,3>& gravity) {
        KRATOS_TRY
        
        double Z_coord = continuum_particle->GetGeometry()[0].Coordinates()[2];
        double water_density = 1000.0; // It is assumed for the time being that we are dealing with water
        double radius = continuum_particle->GetRadius();
        array_1d<double, 3> velocity = continuum_particle->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
        double velocity_module = 1.0; //sqrt(velocity[0] * velocity[0] + velocity[1] * velocity[1] * velocity[2] * velocity[2]); // This would correspond to a quadratic drag law
        
        double water_level = 0.0; // We will always assume this value 
        double Cd = 2.0; // 3.0; // We can experiment a little bit with this
        
        // Adding buoyancy and drag. Now a linear (with the velocity) drag law is implemented
        if (Z_coord < water_level) {
            // Buoyancy:
            // noalias(externally_applied_force) -= this->GetGeometry()[0].FastGetSolutionStepValue(REPRESENTATIVE_VOLUME) * water_density * gravity;
            noalias(externally_applied_force) -= 8.0 * radius * radius * radius * water_density * gravity;
            if (continuum_particle->IsSkin()) {
                DEM_MULTIPLY_BY_SCALAR_3(velocity, 0.5 * Cd * velocity_module * water_density * 4.0 * radius * radius); // The correct area is still to be known
                // Drag force:
                noalias(externally_applied_force) -= velocity;              
            }
        }

        KRATOS_CATCH("")
    }
    
} // namespace Kratos
