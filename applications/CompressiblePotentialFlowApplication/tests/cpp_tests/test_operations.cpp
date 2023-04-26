//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Marco Antonio Zuñiga Perez
//
//

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "fluid_dynamics_application_variables.h"
#include "compressible_potential_flow_application_variables.h"
#include "custom_operations/potential_to_compressible_navier_stokes_operation.h"

namespace Kratos::Testing {

/**
* @brief This method creates a potential_model_part and a "fake" compressible_model_part. 
* The potential nodal values are setting and the nodal velocities are passed with 
* the PotentialToCompressibleNavierStokesOperation to the compressible_model_part.
* The values are check with the explicit local calculate for each node.
*/
KRATOS_TEST_CASE_IN_SUITE(PotentialToCompressibleNavierStokesOperation, CompressiblePotentialApplicationFastSuite)
{
    Model this_model;
    
    // Create potential_model_part
    auto& r_model_part = this_model.CreateModelPart("Main", 3);
    model_part.GetProcessInfo()[DOMAIN_SIZE] = 2;
    
    // Set potential_model_part properties
    BoundedVector<double, 3> free_stream_velocity = ZeroVector(3);
    free_stream_velocity(0) = 10.0;      
    const double gamma = 1.4;
    const double sound_velocity = 340;
    const double free_stream_density = 1.0;
    const double free_stream_mach = 0.2;
    const double reference_temperature = 273;
    model_part.GetProcessInfo()[FREE_STREAM_VELOCITY] = free_stream_velocity;
    model_part.GetProcessInfo()[HEAT_CAPACITY_RATIO] = gamma;
    model_part.GetProcessInfo()[SOUND_VELOCITY] = sound_velocity;
    model_part.GetProcessInfo()[FREE_STREAM_DENSITY] = free_stream_density;
    model_part.GetProcessInfo()[FREE_STREAM_MACH] = free_stream_mach;
    const double specific_heat = (std::pow(sound_velocity,2) / (gamma * reference_temperature)) / (gamma - 1.0);
    
    // Variables addition
    model_part.AddNodalSolutionStepVariable(VELOCITY_POTENTIAL);
    model_part.AddNodalSolutionStepVariable(AUXILIARY_VELOCITY_POTENTIAL);
    
    // Set the element properties
    auto p_elem_prop = model_part.CreateNewProperties(0);
   
    // Geometry creation
    model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    model_part.CreateNewNode(3, 1.0, 1.0, 0.0);
    std::vector<ModelPart::IndexType> elemNodes_1{ 1, 2, 3 };
    auto p_element = model_part.CreateNewElement("IncompressiblePotentialFlowElement2D3N", 1, elemNodes_1, pElemProp);
    r_element.GetValue(WAKE) = 1;
    Vector wake_elemental_distances(3);
    wake_elemental_distances(0) = 1.0;
    wake_elemental_distances(1) = 1.0;
    wake_elemental_distances(2) = -1.0;
    r_element.GetValue(WAKE_ELEMENTAL_DISTANCES) = wake_elemental_distances;
    
    // Define the nodal values
    Vector potential(3);
    potential(0) = 1.0;
    potential(1) = 2.0;
    potential(2) = 3.0;
    for (unsigned int i = 0; i < 3; i++){
      if (wake_elemental_distances(i) > 0.0)
        r_element.GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = potential(i);
      else
        r_element.GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) = potential(i);
    }
    for (unsigned int i = 0; i < 3; i++){
      if (wake_elemental_distances(i) < 0.0)
        r_element.GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = potential(i)+1;
      else
        r_element.GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) = potential(i)+1;
    }      
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // Create fake compressible_model_part
    ModelPart& compressible_model_part = this_model.CreateModelPart("Compressible", 3);
    compressible_model_part.GetProcessInfo()[DOMAIN_SIZE] = 2;
    
    // Variables addition
    compressible_model_part.AddNodalSolutionStepVariable(DENSITY);
    compressible_model_part.AddNodalSolutionStepVariable(MOMENTUM);
    compressible_model_part.AddNodalSolutionStepVariable(TOTAL_ENERGY);
    
    // Set the element properties
    auto p_compressible_elem_prop = compressible_model_part.CreateNewProperties(0);
    
    // Geometry creation
    compressible_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    compressible_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    compressible_model_part.CreateNewNode(3, 1.0, 1.0, 0.0);
    std::vector<ModelPart::IndexType> elemNodes_2{ 1, 2, 3 };
    compressible_model_part.CreateNewElement("Element2D3N", 1, elemNodes_2, compressible_pElemProp);
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // Operation parameters
    Parameters operation_parameters = Parameters(R"(
      {
      "origin_model_part"       : "Main",
      "destination_model_part"  : "Compressible",
      "reference_temperature"   : 273,
      "compute_nodal_velocities": true
      })" );
    
    // Execute PotentialToCompressibleNavierStokesOperation
    PotentialToCompressibleNavierStokesOperation PotentialToCompressibleNavierStokesOperation(this_model, operation_parameters);
    PotentialToCompressibleNavierStokesOperation.Execute();
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   
    for(unsigned int index = 0; index < model_part.NumberOfNodes(); ++index)
    {  
      auto it_node = compressible_model_part.NodesBegin() + index;
      const auto it_potential_node = model_part.NodesBegin() + index;

      // Getting potential flow velocities
      const auto &velocity = potential_node->GetValue(VELOCITY);

      // Check potential velocities and nodal areas
      auto nodal_area = potential_node->GetValue(NODAL_AREA);
      KRATOS_CHECK_NEAR(nodal_area, 0.166667, 1e-6);
      KRATOS_CHECK_NEAR(velocity[0], 1, 1e-6);
      KRATOS_CHECK_NEAR(velocity[1], 1, 1e-6);
      KRATOS_CHECK_NEAR(velocity[2], 0, 1e-6);

      // Explicit calculate of local nodal values
      const double velocity_norm_2 = velocity[0] * velocity[0] + velocity[1] * velocity[1];
      const double velocity_norm = std::sqrt(velocity_norm_2);
      const double mach = velocity_norm / sound_velocity;
      const double num = 1.0 + 0.5 * (gamma - 1.0) * std::pow(free_stream_mach,2);
      const double den = 1.0 + 0.5 * (gamma - 1.0) * std::pow(mach,2);
      const double density = free_stream_density * std::pow((num / den),(1.0 / (gamma - 1.0)));
      const double energy = specific_heat * reference_temperature + 0.5 * velocity_norm_2;

      // Check nodal values
      auto nodal_density = r_node->FastGetSolutionStepValue(DENSITY);
      KRATOS_CHECK_NEAR(nodal_density, density, 1e-6);
      auto& r_nodal_momentum = r_node->FastGetSolutionStepValue(MOMENTUM);
      KRATOS_CHECK_NEAR(nodal_momentum[0], density * velocity[0], 1e-6);
      KRATOS_CHECK_NEAR(nodal_momentum[1], density * velocity[1], 1e-6);
      auto nodal_total_energy = r_node->FastGetSolutionStepValue(TOTAL_ENERGY);
      KRATOS_CHECK_NEAR(nodal_total_energy, density * energy, 1e-6);
    }
}
} // namespace Kratos::Testing
