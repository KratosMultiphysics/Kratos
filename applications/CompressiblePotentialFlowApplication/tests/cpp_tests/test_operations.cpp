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
#include "tests/cpp_tests/compressible_potential_flow_fast_suite.h"
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
    r_model_part.GetProcessInfo()[DOMAIN_SIZE] = 2;
    
    // Set potential_model_part properties
    BoundedVector<double, 3> free_stream_velocity = ZeroVector(3);
    free_stream_velocity(0) = 10.0;      
    const double heat_capacity_ratio = 1.4;
    const double sound_velocity = 340;
    const double reference_temperature = 273;    
    const double free_stream_density = 1.225;
    const double free_stream_mach = 0.2;
    const double specific_heat = (std::pow(sound_velocity,2) / (heat_capacity_ratio * reference_temperature)) / (heat_capacity_ratio - 1.0);
    r_model_part.GetProcessInfo()[FREE_STREAM_VELOCITY] = free_stream_velocity;
    r_model_part.GetProcessInfo()[HEAT_CAPACITY_RATIO] = heat_capacity_ratio;
    r_model_part.GetProcessInfo()[SOUND_VELOCITY] = sound_velocity;
    r_model_part.GetProcessInfo()[FREE_STREAM_DENSITY] = free_stream_density;
    r_model_part.GetProcessInfo()[FREE_STREAM_MACH] = free_stream_mach;
    
    // Variables addition
    r_model_part.AddNodalSolutionStepVariable(VELOCITY_POTENTIAL);
    r_model_part.AddNodalSolutionStepVariable(AUXILIARY_VELOCITY_POTENTIAL);
    
    // Set the element properties
    auto p_elem_prop = r_model_part.CreateNewProperties(0);
   
    // Geometry creation
    r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    r_model_part.CreateNewNode(3, 1.0, 1.0, 0.0);
    std::vector<ModelPart::IndexType> elemNodes_1{ 1, 2, 3 };
    r_model_part.CreateNewElement("IncompressiblePotentialFlowElement2D3N", 1, elemNodes_1, p_elem_prop);
    
    auto& r_element = r_model_part.GetElement(1);

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
    ModelPart& r_compressible_model_part = this_model.CreateModelPart("Compressible", 3);
    r_compressible_model_part.GetProcessInfo()[DOMAIN_SIZE] = 2;
    
    // Variables addition
    r_compressible_model_part.AddNodalSolutionStepVariable(DENSITY);
    r_compressible_model_part.AddNodalSolutionStepVariable(MOMENTUM);
    r_compressible_model_part.AddNodalSolutionStepVariable(TOTAL_ENERGY);
    
    // Set the element properties
    auto p_compressible_elem_prop = r_compressible_model_part.CreateNewProperties(0);
    
    // Geometry creation
    r_compressible_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_compressible_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    r_compressible_model_part.CreateNewNode(3, 1.0, 1.0, 0.0);
    std::vector<ModelPart::IndexType> elemNodes_2{ 1, 2, 3 };
    r_compressible_model_part.CreateNewElement("Element2D3N", 1, elemNodes_2, p_compressible_elem_prop);
    
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
   
    for(unsigned int index = 0; index < r_model_part.NumberOfNodes(); ++index)
    {  
      auto it_node = r_compressible_model_part.NodesBegin() + index;
      const auto it_potential_node = r_model_part.NodesBegin() + index;

      // Getting potential flow velocities
      const auto &r_velocity = it_potential_node->GetValue(VELOCITY);

      // Check potential velocities and nodal areas
      auto nodal_area = it_potential_node->GetValue(NODAL_AREA);
      KRATOS_EXPECT_NEAR(nodal_area, 0.166666667, 1e-9);

      KRATOS_EXPECT_NEAR(r_velocity[0], 1, 1e-9);
      KRATOS_EXPECT_NEAR(r_velocity[1], 1, 1e-9);
      KRATOS_EXPECT_NEAR(r_velocity[2], 0, 1e-9);

      // Explicit calculate of local nodal values
      const double velocity_norm_2 = r_velocity[0] * r_velocity[0] + r_velocity[1] * r_velocity[1] + r_velocity[2] * r_velocity[2];
      const double velocity_norm = std::sqrt(velocity_norm_2);
      const double mach = velocity_norm / sound_velocity;
      const double num = 1.0 + 0.5 * (heat_capacity_ratio - 1.0) * std::pow(free_stream_mach,2);
      const double den = 1.0 + 0.5 * (heat_capacity_ratio - 1.0) * std::pow(mach,2);
      const double density = free_stream_density * std::pow((num / den),(1.0 / (heat_capacity_ratio - 1.0)));
      const double energy = specific_heat * reference_temperature + 0.5 * velocity_norm_2;  

      // Check nodal values
      auto r_nodal_density = it_node->FastGetSolutionStepValue(DENSITY);
      KRATOS_EXPECT_NEAR(r_nodal_density, density, 1e-9);

      auto& r_nodal_momentum = it_node->FastGetSolutionStepValue(MOMENTUM);
      KRATOS_EXPECT_NEAR(r_nodal_momentum[0], density * r_velocity[0], 1e-9);
      KRATOS_EXPECT_NEAR(r_nodal_momentum[1], density * r_velocity[1], 1e-9);
      KRATOS_EXPECT_NEAR(r_nodal_momentum[2], density * r_velocity[2], 1e-9);

      auto r_nodal_total_energy = it_node->FastGetSolutionStepValue(TOTAL_ENERGY);
      KRATOS_EXPECT_NEAR(r_nodal_total_energy, density * energy, 1e-9);
    }
}
} // namespace Kratos::Testing
