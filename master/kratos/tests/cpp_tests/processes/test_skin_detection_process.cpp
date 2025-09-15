//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo
//                   Vicente Mataix Ferrandiz
//

// Project includes
#include "includes/process_info.h"
#include "testing/testing.h"
#include "containers/model.h"
#include "includes/model_part.h"
#include "processes/skin_detection_process.h"

namespace Kratos::Testing
{

void CreateSimpleGeometry(ModelPart& rModelPart)
{
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
    rModelPart.CreateNewNode(4, 0.0, 1.0, 0.0);

    Properties::Pointer p_elem_prop = rModelPart.CreateNewProperties(0);
    Element::Pointer p_elem_1 = rModelPart.CreateNewElement("Element2D3N", 1, {{1,2,3}}, p_elem_prop);
    Element::Pointer p_elem_2 = rModelPart.CreateNewElement("Element2D3N", 2, {{1,3,4}}, p_elem_prop);
}

/**
* Checks the correct work of the skin detection process
* Active elements
*/
KRATOS_TEST_CASE_IN_SUITE(SkinDetectionProcess1, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("test_model_part",2);

    CreateSimpleGeometry(r_model_part);

    Parameters default_parameters = Parameters(R"(
    {
        "name_auxiliar_model_part"              : "SkinModelPart",
        "name_auxiliar_condition"               : "Condition",
        "list_model_parts_to_assign_conditions" : [],
        "echo_level"                            : 0
    })");

    SkinDetectionProcess<2> skin_process = SkinDetectionProcess<2>(
                                           r_model_part, default_parameters);

    // We generate in several iterations to see if it crashes
    for (int i = 0; i < 2; i++) {
        skin_process.Execute();
        KRATOS_EXPECT_EQ(r_model_part.GetSubModelPart("SkinModelPart").NumberOfConditions(), 4);
    }

    // Now we remove one element
    r_model_part.Elements()[2].Set(TO_ERASE);
    r_model_part.RemoveElementsFromAllLevels(TO_ERASE);

    // We execute again
    skin_process.Execute();
    KRATOS_EXPECT_EQ(r_model_part.GetSubModelPart("SkinModelPart").NumberOfConditions(), 3);
}

/**
* Checks the correct work of the skin detection process
* Inactive elements
*/
KRATOS_TEST_CASE_IN_SUITE(SkinDetectionProcess2, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("test_model_part",2);

    CreateSimpleGeometry(r_model_part);

    Parameters default_parameters = Parameters(R"(
    {
        "name_auxiliar_model_part"              : "SkinModelPart",
        "name_auxiliar_condition"               : "Condition",
        "list_model_parts_to_assign_conditions" : [],
        "echo_level"                            : 0
    })");

    SkinDetectionProcess<2> skin_process = SkinDetectionProcess<2>(
                                           r_model_part, default_parameters);
    // Now we inactive one element
    r_model_part.Elements()[2].Set(ACTIVE, false);

    // We execute again
    skin_process.Execute();
    KRATOS_EXPECT_EQ(r_model_part.GetSubModelPart("SkinModelPart").NumberOfConditions(), 3);
}
} // namespace Kratos::Testing
