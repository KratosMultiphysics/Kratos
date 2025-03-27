//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolo' Antonelli
//  Main authors:    Andrea Gorgi

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "custom_processes/assign_iga_external_conditions_process.h"
#include "includes/kratos_parameters.h"
#include "iga_application_variables.h"

namespace Kratos::Testing
{

// Tests the SnakeSbmUProcess with a square outer geometry
KRATOS_TEST_CASE_IN_SUITE(AssignIgaConditionProcessTest, KratosIgaFastSuite)
{
    
    Model model;
    ModelPart& iga_model_part = model.CreateModelPart("IgaModelPart");
    
    iga_model_part.CreateNewProperties(0);

    auto& convection_diffusion_domain_sub_model_part = iga_model_part.CreateSubModelPart("ConvectionDiffusionDomain");
    auto& sbm_outer_sub_model_part = iga_model_part.CreateSubModelPart("SBM_Support_outer");

    convection_diffusion_domain_sub_model_part.CreateNewNode(1, 1.0, 2.0, 0.0);
    convection_diffusion_domain_sub_model_part.CreateNewNode(2, 1.0, 2.0, 0.0); 
    convection_diffusion_domain_sub_model_part.CreateNewElement("LaplacianIGAElement", 1, {{1,2}}, 0);
    
    sbm_outer_sub_model_part.CreateNewCondition("SbmLaplacianConditionDirichlet", 1, {{1,2}}, 0);

    Kratos::Parameters assign_iga_parameters(R"(
        {
            "echo_level": 0,
            "analysis_model_part_name": "IgaModelPart",
            "element_condition_list": [
                {
                    "iga_model_part": "ConvectionDiffusionDomain",
                    "variables": [
                    {
                        "variable_name": "HEAT_FLUX",
                        "value": ["0.0", "0.0", "0.0"]
                    }]
                },
                {
                    "iga_model_part": "SBM_Support_outer",
                    "variables": [
                    {
                        "variable_name": "TEMPERATURE",
                        "value": 0.0
                    }]
                }
            ]
        }
    )");

    AssignIgaExternalConditionsProcess assign_iga_condition_process(model, assign_iga_parameters);    
    assign_iga_condition_process.Execute() ;
    
    
    // // Ensure the number of nodes matches expectation
    // KRATOS_EXPECT_NEAR(surrogate_sub_model_part_outer.NumberOfConditions(), expected_coordinates.size(), tolerance);

}

}
