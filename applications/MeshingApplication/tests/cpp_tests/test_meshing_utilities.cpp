// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                       license: MeshingApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "tests/test_utilities/cpp_tests_utilities.h"

// Application includes
#include "custom_utilities/meshing_utilities.h"
#include "tests/cpp_tests/meshing_fast_suite.h"

namespace Kratos::Testing 
{

typedef Geometry<Node> GeometryType;

void CreateDummy2DNoModelPartPropertiesModelPart(ModelPart& rModelPart)
{
    Properties::Pointer p_elem_prop = Kratos::make_shared<Properties>(0);

    // First we create the nodes
    rModelPart.CreateNewNode(1, 0.0 , 0.0 , 0.0);
    rModelPart.CreateNewNode(2, 1.0 , 0.0 , 0.0);
    rModelPart.CreateNewNode(3, 1.0 , 1.0 , 0.0);
    rModelPart.CreateNewNode(4, 0.0 , 1.0 , 0.0);
    rModelPart.CreateNewNode(5, 2.0 , 0.0 , 0.0);
    rModelPart.CreateNewNode(6, 2.0 , 1.0 , 0.0);

    // Now we create the elements
    rModelPart.CreateNewElement("Element2D3N", 1, {{1,2,3}}, p_elem_prop);
    rModelPart.CreateNewElement("Element2D3N", 2, {{1,3,4}}, p_elem_prop);
    rModelPart.CreateNewElement("Element2D3N", 3, {{2,5,3}}, p_elem_prop);
    rModelPart.CreateNewElement("Element2D3N", 4, {{5,6,3}}, p_elem_prop);
}

/**
* Checks the correct work of the BlockThresholdSizeElements
* Test triangle
*/
KRATOS_TEST_CASE_IN_SUITE(BlockThresholdSizeElements2D, KratosMeshingApplicationFastSuite)
{
    Model this_model;
    ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);
    ProcessInfo& r_current_process_info = r_model_part.GetProcessInfo();
    r_current_process_info[DOMAIN_SIZE] = 2;

    CppTestsUtilities::Create2DGeometry(r_model_part, "Element2D3N");

    Parameters parameters = Parameters(R"(
    {
        "minimal_size" : 2.0,
        "maximal_size" : 10.0
    })" );

    MeshingUtilities::BlockThresholdSizeElements(r_model_part, parameters);

    for (auto& r_element: r_model_part.Elements()) {
        KRATOS_EXPECT_TRUE(r_element.Is(BLOCKED));
        // KRATOS_WATCH(r_element.GetValue(ELEMENT_H))
    }
}

/**
 * Checks the correct work of the BlockThresholdSizeElements
 * Test tetrahedra
*/
KRATOS_TEST_CASE_IN_SUITE(BlockThresholdSizeElements3D, KratosMeshingApplicationFastSuite)
{
    Model this_model;
    ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);
    ProcessInfo& r_current_process_info = r_model_part.GetProcessInfo();
    r_current_process_info[DOMAIN_SIZE] = 3;

    CppTestsUtilities::Create3DGeometry(r_model_part, "Element3D4N");

    Parameters parameters = Parameters(R"(
    {
        "minimal_size" : 2.0,
        "maximal_size" : 10.0
    })" );

    MeshingUtilities::BlockThresholdSizeElements(r_model_part, parameters);

    for (auto& r_element: r_model_part.Elements()) {
        KRATOS_EXPECT_TRUE(r_element.Is(BLOCKED));
        // KRATOS_WATCH(r_element.GetValue(ELEMENT_H))
    }
}

} // namespace Kratos::Testing
