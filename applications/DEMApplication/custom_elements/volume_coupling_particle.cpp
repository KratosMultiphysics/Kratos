//
// Authors:
// Miguel Angel Celigueta maceli@cimne.upc.edu
// Salvador Latorre latorre@cimne.upc.edu
// Miquel Santasusana msantasusana@cimne.upc.edu
// Guillermo Casas gcasas@cimne.upc.edu
// Chengshun Shang cshang@cimne.upc.edu
//

// System includes
#include <string>
#include <iostream>
#include <cmath>

#include <fstream>

// External includes

// Project includes
#include "volume_coupling_particle.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "custom_utilities/discrete_particle_configure.h"
#include "custom_strategies/schemes/glued_to_wall_scheme.h"


namespace Kratos
{


void VolumeCouplingParticle::ComputeBallToRigidFaceContactForceAndMoment(SphericParticle::ParticleDataBuffer & data_buffer,
                                                        array_1d<double, 3>& r_elastic_force,
                                                        array_1d<double, 3>& r_contact_force,
                                                        array_1d<double, 3>& rigid_element_force,
                                                        const ProcessInfo& r_process_info)
{
    // Call the base class's function.
    SphericParticle::ParticleDataBuffer::ComputeBallToRigidFaceContactForceAndMoment(data_buffer, r_elastic_force, r_contact_force, rigid_element_force, r_process_info);
    
 
    double w = /* the value of w */;
    
    // Multiply each element of r_elastic_force and r_contact_force by w.
    for (int i = 0; i < 3; ++i)
    {
        r_elastic_force[i] *= w;
        r_contact_force[i] *= w;
        rigid_element_force[i] *= w;
    }

}

void VolumeCouplingParticle::ComputeBallToBallContactForceAndMoment(SphericParticle::ParticleDataBuffer & data_buffer,
                                                                    const ProcessInfo& r_process_info,
                                                                    array_1d<double, 3>& rElasticForce,
                                                                    array_1d<double, 3>& rContactForce)
{
    // Call the base class's function.
    SphericParticle::ParticleDataBuffer::ComputeBallToBallContactForceAndMoment(data_buffer, r_process_info, rElasticForce, rContactForce);
    
    
    double w = /* the value of w */;
    
    // Multiply each element of rElasticForce and rContactForce by w.
    for (int i = 0; i < 3; ++i)
    {
        rElasticForce[i] *= w;
        rContactForce[i] *= w;
    }
}

void VolumeCouplingParticle::ComputeAdditionalForces(array_1d<double, 3>& externally_applied_force, 
                                                     array_1d<double, 3>& externally_applied_moment, 
                                                     const ProcessInfo& r_process_info, 
                                                     const array_1d<double,3>& gravity)
{
    // Call the base class's function.
    SphericParticle::ComputeAdditionalForces(externally_applied_force, externally_applied_moment, r_process_info, gravity);
    
    // Assuming 'w' is known and is of type 'double'.
    double w = /* the value of w */;
    
    // Multiply each element of externally_applied_force, externally_applied_moment, and gravity by w.
    for (int i = 0; i < 3; ++i)
    {
        externally_applied_force[i] *= w;
        externally_applied_moment[i] *= w;
    }
}


}  // namespace Kratos.
