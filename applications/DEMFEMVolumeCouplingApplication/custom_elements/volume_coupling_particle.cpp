// Project includes
#include "volume_coupling_particle.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "custom_utilities/discrete_particle_configure.h"
#include "custom_strategies/schemes/glued_to_wall_scheme.h"
#include "../DEMFEM_volume_coupling_application.h"


namespace Kratos {



VolumeCouplingParticle::VolumeCouplingParticle( IndexType NewId, GeometryType::Pointer pGeometry )
    : SphericParticle( NewId, pGeometry )
{
    
}
VolumeCouplingParticle::VolumeCouplingParticle( ) // Default constructor needed for serialization
    : SphericParticle( )
{
    
}

// Function to compute contact forces and moments between ball and rigid face considering coupling weight.
void VolumeCouplingParticle::ComputeBallToRigidFaceContactForceAndMoment(SphericParticle::ParticleDataBuffer & data_buffer,
                                                        array_1d<double, 3>& r_elastic_force,
                                                        array_1d<double, 3>& r_contact_force,
                                                        array_1d<double, 3>& rigid_element_force,
                                                        const ProcessInfo& r_process_info)
{
    // Call the base class's function.
    SphericParticle::ComputeBallToRigidFaceContactForceAndMoment(data_buffer, r_elastic_force, r_contact_force, rigid_element_force, r_process_info);
    double particle_weight = this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_COUPLING_WEIGHT); // Coupling weight of the particle center

         r_elastic_force *= particle_weight; // Scale elastic force by the weight
         r_contact_force *= particle_weight; // Scale contact force by the  weight
         rigid_element_force *= particle_weight; // Scale rigid element force by the  weight


    // // Calculate the coordinates of the contact point
    // Node<3>& central_node = this->GetGeometry()[0];  // Central node of the geometry
    // double particle_radius = this->GetRadius(); // Radius of the particle
    // array_1d<double, 3> particle_centre_coordinates = central_node.Coordinates(); // Coordinates of the particle center
    // array_1d<double, 3> contact_point_coordinates = particle_centre_coordinates - data_buffer.mLocalCoordSystem[2] * particle_radius; // Contact point coordinates

    // // Get the weight at the particle center
    // double particle_weight = central_node.FastGetSolutionStepValue(PARTICLE_COUPLING_WEIGHT); // Coupling weight of the particle center

    // // Check if there are neighbour elements and if yes, use the first one
    // if (!data_buffer.mpThisParticle->mNeighbourElements.empty()) {
    //     SphericParticle& neighbor_particle = *data_buffer.mpThisParticle->mNeighbourElements[0]; // First neighbouring element

    //     // Calculate the slope using the particle center and the neighbor
    //     double neighbor_weight = neighbor_particle.GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_COUPLING_WEIGHT); // Neighbour's coupling weight
    //     double neighbor_y_coordinate = neighbor_particle.GetGeometry()[0].Coordinates()[1]; // Neighbour's y-coordinate

    //     double slope = (neighbor_weight - particle_weight) / (neighbor_y_coordinate - particle_centre_coordinates[1]); // Calculating the slope

    //     // Calculate the weight at the contact point, assuming linear variation with y-coordinate
    //     double contact_point_weight = particle_weight + slope * (contact_point_coordinates[1] - particle_centre_coordinates[1]); // Weight at the contact point

    //     r_elastic_force *= contact_point_weight; // Scale elastic force by the contact point weight
    //     r_contact_force *= contact_point_weight; // Scale contact force by the contact point weight
    //     rigid_element_force *= contact_point_weight; // Scale rigid element force by the contact point weight
    // }
}

// Function to compute forces for positive indentations between balls considering coupling weight.
void VolumeCouplingParticle::EvaluateBallToBallForcesForPositiveIndentiations(SphericParticle::ParticleDataBuffer & data_buffer,
                                                            const ProcessInfo& r_process_info,
                                                            double LocalElasticContactForce[3],
                                                            double DeltDisp[3],
                                                            double LocalDeltDisp[3],
                                                            double RelVel[3],
                                                            double indentation,
                                                            double ViscoDampingLocalContactForce[3],
                                                            double& cohesive_force,
                                                            SphericParticle* element2,
                                                            bool& sliding,
                                                            double LocalCoordSystem[3][3],
                                                            double OldLocalCoordSystem[3][3],
                                                            array_1d<double, 3>& neighbour_elastic_contact_force)
{
    // Call the base class's function.
    SphericParticle::EvaluateBallToBallForcesForPositiveIndentiations(data_buffer, r_process_info, LocalElasticContactForce,
                                                                        DeltDisp,LocalDeltDisp,RelVel,
                                                                        indentation,ViscoDampingLocalContactForce,
                                                                        cohesive_force,element2,sliding,LocalCoordSystem,
                                                                        OldLocalCoordSystem,neighbour_elastic_contact_force);
  


    const double other_young = element2->GetYoung();
    const double my_young = this->GetYoung();

    const double inverse_of_sum_of_youngs = 1.0 / (other_young + my_young);
    const double my_arm_length = this->GetRadius() - indentation * other_young * inverse_of_sum_of_youngs;
    const double other_arm_length  = element2->GetRadius() - indentation * my_young * inverse_of_sum_of_youngs;
    double interpolated_weight = (this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_COUPLING_WEIGHT)*other_arm_length+element2->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_COUPLING_WEIGHT)*my_arm_length)/(my_arm_length+other_arm_length);

    // Calculate the interpolated weight based on the current and neighboring particles, w'=w1+((3r1+r2-|y1-y2|)/2|y1-y2|)*(w2-w1) considering some indentation
    // double interpolated_weight = this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_COUPLING_WEIGHT)+
    // (( 3* this->GetGeometry()[0].FastGetSolutionStepValue(RADIUS) + element2->GetGeometry()[0].FastGetSolutionStepValue(RADIUS)-
    // std::abs(this->GetGeometry()[0][1]-element2->GetGeometry()[0][1]))/2*std::abs(this->GetGeometry()[0][1]-element2->GetGeometry()[0][1]))*(element2->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_COUPLING_WEIGHT)-this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_COUPLING_WEIGHT));

    // Scale forces by the interpolated weight
    for(int i=0; i<3; ++i)
    {
        LocalElasticContactForce[i] *= interpolated_weight; // Scale elastic contact force
        ViscoDampingLocalContactForce[i] *= interpolated_weight; // Scale visco-damping contact force
    }

    cohesive_force *=interpolated_weight; // Scale cohesive force

}

// Function to compute additional forces and moments on the particle considering coupling weight.
void VolumeCouplingParticle::ComputeAdditionalForces(array_1d<double, 3>& externally_applied_force, 
                                                     array_1d<double, 3>& externally_applied_moment, 
                                                     const ProcessInfo& r_process_info, 
                                                     const array_1d<double,3>& gravity)
{
    // Call the base class's function.
    SphericParticle::ComputeAdditionalForces(externally_applied_force, externally_applied_moment, r_process_info, gravity);
    
    
    double w = this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_COUPLING_WEIGHT); // Coupling weight of the particle center
    
    // Multiply each element of externally_applied_force and externally_applied_moment by w.
    for (int i = 0; i < 3; ++i)
    {
        externally_applied_force[i] *= w; // Scale externally applied force
        externally_applied_moment[i] *= w; // Scale externally applied moment
    }
}

}  // namespace Kratos.

// System includes
// #include <string>
// #include <iostream>
// #include <cmath>

// #include <fstream>

// External includes

// // Project includes
// #include "volume_coupling_particle.h"
// #include "custom_utilities/GeometryFunctions.h"
// #include "custom_utilities/AuxiliaryFunctions.h"
// #include "custom_utilities/discrete_particle_configure.h"
// #include "custom_strategies/schemes/glued_to_wall_scheme.h"


// namespace Kratos
// {



// void VolumeCouplingParticle::ComputeBallToRigidFaceContactForceAndMoment(SphericParticle::ParticleDataBuffer & data_buffer,
//                                                         array_1d<double, 3>& r_elastic_force,
//                                                         array_1d<double, 3>& r_contact_force,
//                                                         array_1d<double, 3>& rigid_element_force,
//                                                         const ProcessInfo& r_process_info)
// {
//     // Call the base class's function.
//     SphericParticle::ComputeBallToRigidFaceContactForceAndMoment(data_buffer, r_elastic_force, r_contact_force, rigid_element_force, r_process_info);

//     // Calculate the coordinates of the contact point
//     Node<3>& central_node = this->GetGeometry()[0];
//     double particle_radius = this->GetRadius();
//     array_1d<double, 3> particle_centre_coordinates = central_node.Coordinates();
//     array_1d<double, 3> contact_point_coordinates = particle_centre_coordinates - data_buffer.mLocalCoordSystem[2] * particle_radius;

//     // Get the weight at the particle center
//     double particle_weight = central_node.FastGetSolutionStepValue(PARTICLE_COUPLING_WEIGHT);

//     // Check if there are neighbour elements and if yes, use the first one
//     if (!data_buffer.mpThisParticle->mNeighbourElements.empty()) {
//         SphericParticle& neighbor_particle = *data_buffer.mpThisParticle->mNeighbourElements[0];

//         // Calculate the slope using the particle center and the neighbor
//         double neighbor_weight = neighbor_particle.GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_COUPLING_WEIGHT);
//         double neighbor_y_coordinate = neighbor_particle.GetGeometry()[0].Coordinates()[1];

//         double slope = (neighbor_weight - particle_weight) / (neighbor_y_coordinate - particle_centre_coordinates[1]);

//         // Calculate the weight at the contact point, assuming linear variation with y-coordinate
//         double contact_point_weight = particle_weight + slope * (contact_point_coordinates[1] - particle_centre_coordinates[1]);

//         r_elastic_force *= contact_point_weight;
//         r_contact_force *= contact_point_weight;
//         rigid_element_force *= contact_point_weight;
// }
// }


// virtual void EvaluateBallToBallForcesForPositiveIndentiations(SphericParticle::ParticleDataBuffer & data_buffer,
//                                                             const ProcessInfo& r_process_info,
//                                                             double LocalElasticContactForce[3],
//                                                             double DeltDisp[3],
//                                                             double LocalDeltDisp[3],
//                                                             double RelVel[3],
//                                                             double indentation,
//                                                             double ViscoDampingLocalContactForce[3],
//                                                             double& cohesive_force,
//                                                             SphericParticle* element2,
//                                                             bool& sliding,
//                                                             double LocalCoordSystem[3][3],
//                                                             double OldLocalCoordSystem[3][3],
//                                                             array_1d<double, 3>& neighbour_elastic_contact_force)
// {
//     // Call the base class's function.
//     SphericParticle::EvaluateBallToBallForcesForPositiveIndentiations(data_buffer, r_process_info, LocalElasticContactForce[3],
//                                                                         DeltDisp[3],LocalDeltDisp[3],RelVel[3],
//                                                                         indentation,ViscoDampingLocalContactForce[3],
//                                                                         cohesive_force,element2,sliding,LocalCoordSystem[3][3],
//                                                                         OldLocalCoordSystem[3][3],neighbour_elastic_contact_force);
    
  
//     double interpolated_weight = this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_COUPLING_WEIGHT)+
//     (( 3* this->GetGeometry()[0].FastGetSolutionStepValue(RADIUS) + element2->GetGeometry()[0].FastGetSolutionStepValue(RADIUS)-
//     std::abs(this->GetGeometry()[0][1]-element2->GetGeometry()[0][1]))/2*std::abs(this->GetGeometry()[0][1]-element2->GetGeometry()[0][1]))*(element2->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_COUPLING_WEIGHT)-this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_COUPLING_WEIGHT));

//     for(int i=0; i<3; ++i)
//     {
//         LocalElasticContactForce[i] *= interpolated_weight;
//         ViscoDampingLocalContactForce[i] *= interpolated_weight;
//     }

//     cohesive_force *=interpolated_weight;

// }

// void VolumeCouplingParticle::ComputeAdditionalForces(array_1d<double, 3>& externally_applied_force, 
//                                                      array_1d<double, 3>& externally_applied_moment, 
//                                                      const ProcessInfo& r_process_info, 
//                                                      const array_1d<double,3>& gravity)
// {
//     // Call the base class's function.
//     SphericParticle::ComputeAdditionalForces(externally_applied_force, externally_applied_moment, r_process_info, gravity);
    
//     // Assuming 'w' is known and is of type 'double'.
//     double w = /* the value of w */;
    
//     // Multiply each element of externally_applied_force, externally_applied_moment, and gravity by w.
//     for (int i = 0; i < 3; ++i)
//     {
//         externally_applied_force[i] *= w;
//         externally_applied_moment[i] *= w;
//     }
// }


// }  // namespace Kratos.

