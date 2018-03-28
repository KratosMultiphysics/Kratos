// Created by: Salva Latorre, latorre@cimne.upc.edu

// System includes
//#include <cmath> // for the sine in case we activate the wave equation

// Project includes
#include "ice_continuum_particle.h"

namespace Kratos {

    array_1d<double, 3> IceContinuumParticle::ComputeWeight(const array_1d<double, 3>& gravity, const ProcessInfo& r_process_info) {
    
        KRATOS_TRY

        array_1d<double, 3> total_weight = ZeroVector(3);
        double radius = this->GetRadius();
        array_1d<double, 3> velocity = this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
        double velocity_module = 1.0; //sqrt(velocity[0] * velocity[0] + velocity[1] * velocity[1] * velocity[2] * velocity[2]); // This would correspond to a quadratic (with the velocity) drag law
        double Cd = 2.0; // 3.0; // Drag coefficient. We can experiment a little bit with this number
        double water_level = 0.0; // We will always assume this value
        double density_of_fluid = 0.0; // At the moment we are surrounded by air
        
        // IN CASE OF A SINUSOIDAL WATER LEVEL IN XZ:
        // double amplitude = 1.0; double wave_celerity = 1.0; double wavelength = 8.0;
        // water_level = amplitude * sin ((2.0 * Globals::Pi * ((this->GetGeometry()[0].Coordinates()[0]) - wave_celerity * r_process_info[TIME])) / wavelength);
        
        // Adding buoyancy and drag. Now a linear (with the velocity) drag law is implemented
        if (this->GetGeometry()[0].Coordinates()[2] < water_level) { // Water level is measured in Z axis
            density_of_fluid = 1000.0; // If we are here, we are submerged into water
            if (this->IsSkin()) {
                DEM_MULTIPLY_BY_SCALAR_3(velocity, 0.5 * Cd * velocity_module * density_of_fluid * 4.0 * radius * radius); // The correct area is still to be known
                total_weight -= velocity; // Drag force             
            }
        }
        
        return total_weight += this->GetGeometry()[0].FastGetSolutionStepValue(REPRESENTATIVE_VOLUME) * gravity * (GetDensity() - density_of_fluid);

        KRATOS_CATCH("")
    }
} // namespace Kratos
