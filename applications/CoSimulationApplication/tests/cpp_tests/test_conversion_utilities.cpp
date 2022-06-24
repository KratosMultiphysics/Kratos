// KRATOS  / ___|___/ ___|(_)_ __ ___  _   _| | __ _| |_(_) ___  _ ___
//        | |   / _ \___ \| | '_ ` _ \| | | | |/ _` | __| |/ _ \| '_  |
//        | |__| (_) |__) | | | | | | | |_| | | (_| | |_| | (_) | | | |
//         \____\___/____/|_|_| |_| |_|\__,_|_|\__,_|\__|_|\___/|_| |_|
//
//  License:		 BSD License
//                   license: CoSimulationApplication/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//                   Ashish Darekar
//

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "containers/model.h"
#include "testing/testing.h"
#include "custom_utilities/conversion_utilities.h"

// Utilities
#include "utilities/cpp_tests_utilities.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(ElementalDataToNodalDataForce, KratosCosimulationFastSuite)
{
    // Creating a Model
    Model model;
    auto &test_model_part = model.CreateModelPart("TestModelPart");
    test_model_part.SetBufferSize(1);
    test_model_part.AddNodalSolutionStepVariable(FORCE);

    // Creating Nodes
    test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    test_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    test_model_part.CreateNewNode(3, 1.0, 1.0, 0.0);
    test_model_part.CreateNewNode(4, 0.0, 1.0, 0.0);
    test_model_part.CreateNewNode(5, 2.0, 0.0, 0.0);
    test_model_part.CreateNewNode(6, 2.0, 1.0, 0.0);

    // Creating Elements
    auto p_props = test_model_part.CreateNewProperties(0);
    std::vector<ModelPart::IndexType> elem_nodes_1{1, 2, 3, 4};
    test_model_part.CreateNewElement("Element2D4N", 1, elem_nodes_1, p_props);
    std::vector<ModelPart::IndexType> elem_nodes_2{2, 5, 6, 3};
    test_model_part.CreateNewElement("Element2D4N", 2, elem_nodes_2, p_props);

    // Assign values at elements
    array_1d<double, 3> force_value {12.0, 8.0, 0.0};
    for(auto& r_elem : test_model_part.Elements()){
        r_elem.SetValue(FORCE, force_value);
    }

    // ConversionUtilities::ConvertElementalDataToNodalData(test_model_part);
    ConversionUtilities::ConvertElementalDataToNodalData(test_model_part, FORCE, FORCE);

    // Converted Values at Nodes
    Vector nodal_force_values(18);
    int counter = 0;

    for (auto& r_node : test_model_part.Nodes()){
        auto nodal_force = r_node.FastGetSolutionStepValue(FORCE);
        for(int i = 0 ; i< 3; i++){
            nodal_force_values[counter++] = nodal_force[i];
        }
    }

    // Expected Values at Nodes
    array_1d<double, 18> expected_values {3, 2 ,0, 6, 4, 0, 6, 4 ,0, 3, 2 ,0, 3, 2 ,0, 3, 2 ,0};

    KRATOS_CHECK_VECTOR_NEAR(nodal_force_values, expected_values, 1.0e-12)
}

} // namespace Testing
} // namespace Kratos
