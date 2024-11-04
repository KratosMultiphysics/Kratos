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
VolumeCouplingParticle::VolumeCouplingParticle(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : SphericParticle(NewId, pGeometry, pProperties)
{

}

Element::Pointer VolumeCouplingParticle::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new VolumeCouplingParticle(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

// Function to compute contact forces and moments between ball and rigid face considering coupling weight.
void VolumeCouplingParticle::ComputeBallToRigidFaceContactForceAndMoment(
    SphericParticle::ParticleDataBuffer & data_buffer,
    array_1d<double, 3>& r_elastic_force,
    array_1d<double, 3>& r_contact_force,
    array_1d<double, 3>& rigid_element_force,
    const ProcessInfo& r_process_info)
{
    // Call the base class's function.
    SphericParticle::ComputeBallToRigidFaceContactForceAndMoment(data_buffer, r_elastic_force, r_contact_force, rigid_element_force, r_process_info);

        // double particle_weight = this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_COUPLING_WEIGHT); // Coupling weight of the particle center
        // // to check if particle coupling weight is non zer0,
        //  if (particle_weight!=0)
        //  {
        //   r_elastic_force *= particle_weight; // Scale elastic force by the weight
        //     if (std::sqrt(r_elastic_force[0] * r_elastic_force[0] + 
        //                 r_elastic_force[1] * r_elastic_force[1] + 
        //                 r_elastic_force[2] * r_elastic_force[2]) != 0) 
        //     {
        //         // std::cout << r_elastic_force;
        //     }

        //     r_contact_force *= particle_weight; // Scale contact force by the weight
        //     if (std::sqrt(r_contact_force[0] * r_contact_force[0] + 
        //                 r_contact_force[1] * r_contact_force[1] + 
        //                 r_contact_force[2] * r_contact_force[2]) != 0) 
        //     {
        //         // std::cout << r_contact_force;
        //     }

        //     rigid_element_force *= particle_weight; // Scale rigid element force by the weight
        //     if (std::sqrt(rigid_element_force[0] * rigid_element_force[0] + 
        //                 rigid_element_force[1] * rigid_element_force[1] + 
        //                 rigid_element_force[2] * rigid_element_force[2]) != 0) 
        //     {
        //         // std::cout << rigid_element_force;
        //     }
        //  }
}

double VolumeCouplingParticle::GetMass ()
  { 
    double Pparticle_weight = this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_COUPLING_WEIGHT); 
    // std::cout<<"particle ID: "<<this->Id()<<" mass before = "<<mRealMass<<" mass after = "<<mRealMass*Pparticle_weight<<" particle weight = "<<Pparticle_weight<<std::endl;
    //to check if particle id is 1 or not
    // if (this->Id()==1)
    // {
    //     Pparticle_weight=1;
    // }
    // else
    // {
    //     Pparticle_weight=0.5;
    // }
    return mRealMass*Pparticle_weight;     
  }
void VolumeCouplingParticle::Initialize(const ProcessInfo& r_process_info)
{
    
    // Initialize the particle coupling weight
    this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_COUPLING_WEIGHT) = 1.0;

    // Call the base class's function.
    SphericParticle::Initialize(r_process_info);

    
}
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
    double particle_weight = this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_COUPLING_WEIGHT);
    double other_particle_weight = element2->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_COUPLING_WEIGHT);
    // unweight the forces
    if(particle_weight!=1 && other_particle_weight!=1)
    {

        const double other_young = element2->GetYoung();
        const double my_young = this->GetYoung();
        // KRATOS_WATCH(r_process_info[TIME])
        // KRATOS_WATCH(this->Id())
        // KRATOS_WATCH(element2->Id())
        const double inverse_of_sum_of_youngs = 1.0 / (other_young + my_young);
        const double my_arm_length = this->GetRadius() - indentation * other_young * inverse_of_sum_of_youngs;
        const double other_arm_length  = element2->GetRadius() - indentation * my_young * inverse_of_sum_of_youngs;
        double interpolated_weight = (particle_weight*other_arm_length+other_particle_weight*my_arm_length)/(my_arm_length+other_arm_length);

        // Scale forces by the interpolated weight
        // KRATOS_WATCH(r_process_info[TIME])
        //std::cout<<"Before interpolation Ball_contact particle Id "<<this->Id()<<" neighbour id "<<data_buffer.mpOtherParticle->Id()<<" "<<LocalElasticContactForce[0]<<" "<<LocalElasticContactForce[1]<<" "<<LocalElasticContactForce[2]<<std::endl<<std::flush;
        // KRATOS_WATCH(interpolated_weight)
        // KRATOS_WATCH(particle_weight)
        for(int i=0; i<3; ++i)
        {
            LocalElasticContactForce[i] = LocalElasticContactForce[i]/(interpolated_weight); // unScale elastic contact force
            ViscoDampingLocalContactForce[i] = ViscoDampingLocalContactForce[i]/(interpolated_weight); // unScale visco-damping contact force
        }

        cohesive_force =cohesive_force/(interpolated_weight); // unScale cohesive force
    }

    SphericParticle::EvaluateBallToBallForcesForPositiveIndentiations(data_buffer, r_process_info, LocalElasticContactForce,
                                                                        DeltDisp,LocalDeltDisp,RelVel,
                                                                        indentation,ViscoDampingLocalContactForce,
                                                                        cohesive_force,element2,sliding,LocalCoordSystem,
                                                                        OldLocalCoordSystem,neighbour_elastic_contact_force);
  


        // to rule out weighting contact forces with particle outsides of the hybrid domain
    if(particle_weight!=1 && other_particle_weight!=1)
    {

        const double other_young = element2->GetYoung();
        const double my_young = this->GetYoung();
        // KRATOS_WATCH(r_process_info[TIME])
        // KRATOS_WATCH(this->Id())
        // KRATOS_WATCH(element2->Id())
        const double inverse_of_sum_of_youngs = 1.0 / (other_young + my_young);
        const double my_arm_length = this->GetRadius() - indentation * other_young * inverse_of_sum_of_youngs;
        const double other_arm_length  = element2->GetRadius() - indentation * my_young * inverse_of_sum_of_youngs;
        double interpolated_weight = (particle_weight*other_arm_length+other_particle_weight*my_arm_length)/(my_arm_length+other_arm_length);

        // Scale forces by the interpolated weight
        // KRATOS_WATCH(r_process_info[TIME])

        //std::cout<<"Before interpolation Ball_contact particle Id "<<this->Id()<<" neighbour id "<<data_buffer.mpOtherParticle->Id()<<" "<<LocalElasticContactForce[0]<<" "<<LocalElasticContactForce[1]<<" "<<LocalElasticContactForce[2]<<std::endl<<std::flush;
        // KRATOS_WATCH(interpolated_weight)
        // KRATOS_WATCH(particle_weight)
        for(int i=0; i<3; ++i)
        {
            LocalElasticContactForce[i] *= (interpolated_weight); // Scale elastic contact force
            ViscoDampingLocalContactForce[i] *= (interpolated_weight); // Scale visco-damping contact force
        }

        cohesive_force *=(interpolated_weight); // Scale cohesive force
    }


    //std::cout<<"After interpolation Ball_contact particle Id "<<this->Id()<<" neighbour id "<<data_buffer.mpOtherParticle->Id()<<" "<<LocalElasticContactForce[0]<<" "<<LocalElasticContactForce[1]<<" "<<LocalElasticContactForce[2]<<std::endl<<std::flush;
       
    // std::cout<<"LocalElasticContactForce "<<LocalElasticContactForce[0]<<" " << LocalElasticContactForce[1]<<" " << LocalElasticContactForce[2]<<std::endl;
}

void VolumeCouplingParticle::EvaluateBallToRigidFaceForcesForPositiveIndentations(SphericParticle::ParticleDataBuffer &data_buffer,
                                                                   const int rigid_neighbour_index,
                                                                   const double DeltVel[3],
                                                                   const ProcessInfo& r_process_info,
                                                                   double OldLocalElasticContactForce[3],
                                                                   double LocalElasticContactForce[3],
                                                                   double LocalDeltDisp[3],
                                                                   const double indentation,
                                                                   const double  previous_indentation,
                                                                   double ViscoDampingLocalContactForce[3],
                                                                   double& cohesive_force,
                                                                   Condition* const wall,
                                                                   bool& sliding)
{


    double particle_weight = this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_COUPLING_WEIGHT);

    for(int i=0; i<3; ++i)
        {
            LocalElasticContactForce[i] = LocalElasticContactForce[i]/(particle_weight); // unScale elastic contact force
            ViscoDampingLocalContactForce[i] = ViscoDampingLocalContactForce[i]/(particle_weight); // unScale visco-damping contact force
        }

    cohesive_force =cohesive_force/(particle_weight); // unScale cohesive force



    // Call the base class's function.
    SphericParticle::EvaluateBallToRigidFaceForcesForPositiveIndentations(data_buffer, rigid_neighbour_index, DeltVel,
                                                                        r_process_info,OldLocalElasticContactForce,LocalElasticContactForce,
                                                                        LocalDeltDisp,indentation,
                                                                        previous_indentation,ViscoDampingLocalContactForce,cohesive_force,wall,
                                                                        sliding);
  

    
    // print particle weights
    // std::cout<< " Particle weights"<<particle_weight<<" \n";
    // std::cout<<"TESTING########################################## EvaluateBallToRigidFaceForcesForPositiveIndentations\n";
    // to rule out weighting contact forces with particle outsides of the hybrid domain
    for(int i=0; i<3; ++i)
        {
            LocalElasticContactForce[i] *= (particle_weight); // Scale elastic contact force
            ViscoDampingLocalContactForce[i] *= (particle_weight); // Scale visco-damping contact force
        }

    cohesive_force *=(particle_weight); // Scale cohesive force

    //std::cout<<"After interpolation Ball_contact particle Id "<<this->Id()<<" neighbour id "<<data_buffer.mpOtherParticle->Id()<<" "<<LocalElasticContactForce[0]<<" "<<LocalElasticContactForce[1]<<" "<<LocalElasticContactForce[2]<<std::endl<<std::flush;
       
    // std::cout<<"LocalElasticContactForce "<<LocalElasticContactForce[0]<<" " << LocalElasticContactForce[1]<<" " << LocalElasticContactForce[2]<<std::endl;
}

// Function to compute additional forces and moments on the particle considering coupling weight.
// void VolumeCouplingParticle::ComputeAdditionalForces(array_1d<double, 3>& externally_applied_force, 
//                                                      array_1d<double, 3>& externally_applied_moment, 
//                                                      const ProcessInfo& r_process_info, 
//                                                      const array_1d<double,3>& gravity)
// {
//     // statement for setting EXTERNAL_APPLIED_FORCE = DEMFEM_VOLUME_COUPLING_FORCE


//     // Call the base class's function.
//     SphericParticle::ComputeAdditionalForces(externally_applied_force, externally_applied_moment, r_process_info, gravity);
 
//     // double particle_weight = this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_COUPLING_WEIGHT); // Coupling weight of the particle center

//     //     // to check if particle coupling weight is non zer0 and to check multiple calls of this function ,
//     // if (particle_weight!=0)
//     // {
//     // // Multiply each element of externally_applied_force and externally_applied_moment by weight.
//     //     for (int i = 0; i < 3; ++i)
//     //     {
//     //         externally_applied_force[i] *= particle_weight; // Scale externally applied force
//     //         externally_applied_moment[i] *= particle_weight; // Scale externally applied moment
//     //     }
//     //    externally_applied_force += this->GetGeometry()[0].FastGetSolutionStepValue(DEMFEM_VOLUME_COUPLING_FORCE);  //adding the coupling forces exreted on particles
//     //     //KRATOS_WATCH(this->Id());
//         //KRATOS_WATCH(*externally_applied_force);
//         //std::cout<<"particle_id"<<this->Id()<<std::endl;
//         //std::cout<<"particle coupling force"<<this->GetGeometry()[0].FastGetSolutionStepValue(DEMFEM_VOLUME_COUPLING_FORCE)<<std::endl;
//         //std::cout<<" externally_applied_force"<<externally_applied_force<<std::endl;
//    // }

// }

}  // namespace Kratos.




// void VolumeCouplingParticle::EvaluateBallToRigidFaceForcesForPositiveIndentations(SphericParticle::ParticleDataBuffer &data_buffer,
//                                                                    const int rigid_neighbour_index,
//                                                                    double DeltVel[3],
//                                                                    const ProcessInfo& r_process_info,
//                                                                    double OldLocalElasticContactForce[3],
//                                                                    double LocalElasticContactForce[3],
//                                                                    double LocalDeltDisp[3],
//                                                                    const double indentation,
//                                                                    const double  previous_indentation,
//                                                                    double ViscoDampingLocalContactForce[3],
//                                                                    double& cohesive_force,
//                                                                    Condition* const wall,
//                                                                    bool& sliding)
//  {

//      // Call the base class's function.
//     SphericParticle::EvaluateBallToRigidFaceForcesForPositiveIndentations(data_buffer,
//                                                                     rigid_neighbour_index,
//                                                                     DeltVel,
//                                                                     r_process_info,
//                                                                     OldLocalElasticContactForce,
//                                                                     LocalElasticContactForce,
//                                                                     LocalDeltDisp,
//                                                                     indentation,
//                                                                     previous_indentation,
//                                                                     ViscoDampingLocalContactForce,
//                                                                     cohesive_force,
//                                                                     wall,
//                                                                     sliding);

// double weight_at_contact_point =0;
// array_1d<double, 4>& Weight = this->mContactConditionWeights[rigid_neighbour_index];

// for (IndexType i = 0; i < wall->GetGeometry().size(); ++i)
//  {
//      weight_at_contact_point += Weight[i]*wall->GetGeometry()[i].GetSolutionStepValue(NODAL_COUPLING_WEIGHT);
//  }

//  for(int i=0; i<3; ++i)
//  {
//     LocalElasticContactForce[i] *=  weight_at_contact_point; // Scale elastic contact force
//     ViscoDampingLocalContactForce[i] *=  weight_at_contact_point; // Scale visco-damping contact force
//  }                                        
//  cohesive_force *= weight_at_contact_point; // Scale cohesive force
//  }

// Function to compute forces for positive indentations between balls considering coupling weight.



// // Function to compute contact forces and moments between ball and rigid face considering coupling weight.
// void VolumeCouplingParticle::ComputeBallToRigidFaceContactForceAndMoment(
//     SphericParticle::ParticleDataBuffer & data_buffer,
//     array_1d<double, 3>& r_elastic_force,
//     array_1d<double, 3>& r_contact_force,
//     array_1d<double, 3>& rigid_element_force,
//     const ProcessInfo& r_process_info)
// {
//     // Call the base class's function.
//     SphericParticle::ComputeBallToRigidFaceContactForceAndMoment(data_buffer, r_elastic_force, r_contact_force, rigid_element_force, r_process_info);

//     // Retrieve the list of point condition pointers
//     auto& list_of_point_condition_pointers = this->GetValue(WALL_POINT_CONDITION_POINTERS);
//     auto& neighbour_point_faces_elastic_contact_force = this->GetValue(WALL_POINT_CONDITION_ELASTIC_FORCES);
//     auto& neighbour_point_faces_total_contact_force = this->GetValue(WALL_POINT_CONDITION_TOTAL_FORCES);

//     for (unsigned int i = 0; i < list_of_point_condition_pointers.size(); i++) {
//         Condition* wall = list_of_point_condition_pointers[i];

//         // Retrieve the shape functions
//         const Matrix& r_N = wall->GetGeometry().ShapeFunctionsValues();
//         //wall->GetGEometry()

//         // Initialize the weight at the contact point
//         double weight_at_contact_point = 0.0;

//         // Interpolate the weight at the contact point using the shape functions
//         for (IndexType i = 0; i < wall->GetGeometry().size(); ++i) {
//             double nodal_weight = wall->GetGeometry()[i].GetValue(NODAL_COUPLING_WEIGHT);
//             weight_at_contact_point += r_N(0, i) * nodal_weight;
//         }

//         // Scale the contact forces using the weight at the contact point
//         array_1d<double, 3>& LocalContactForce = neighbour_point_faces_elastic_contact_force[i];
//         array_1d<double, 3>& GlobalContactForce = neighbour_point_faces_total_contact_force[i];
//         LocalContactForce *= weight_at_contact_point;
//         GlobalContactForce *= weight_at_contact_point;
//     }
// }





    //double particle_weight = this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_COUPLING_WEIGHT); // Coupling weight of the particle center
        //  r_elastic_force *= particle_weight; // Scale elastic force by the weight
        //  r_contact_force *= particle_weight; // Scale contact force by the  weight
        //  rigid_element_force *= particle_weight; // Scale rigid element force by the  weight
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

