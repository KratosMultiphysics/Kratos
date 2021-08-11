// KRATOS  / ___|___/ ___|(_)_ __ ___  _   _| | __ _| |_(_) ___  _ ___
//        | |   / _ \___ \| | '_ ` _ \| | | | |/ _` | __| |/ _ \| '_  |
//        | |__| (_) |__) | | | | | | | |_| | | (_| | |_| | (_) | | | |
//         \____\___/____/|_|_| |_| |_|\__,_|_|\__,_|\__|_|\___/|_| |_|
//
//  License:		 BSD License
//                   license: CoSimulationApplication/license.txt
//
//  Main authors:    Philipp Bucher
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
    Model model;
    auto &test_model_part = model.CreateModelPart("TestModelPart");
    test_model_part.SetBufferSize(1);
    test_model_part.AddNodalSolutionStepVariable(FORCE);

    //Node creation
    test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    test_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    test_model_part.CreateNewNode(3, 1.0, 1.0, 0.0);
    test_model_part.CreateNewNode(4, 0.0, 1.0, 0.0);
    test_model_part.CreateNewNode(5, 2.0, 0.0, 0.0);
    test_model_part.CreateNewNode(6, 2.0, 1.0, 0.0);

    std::vector<ModelPart::IndexType> elem_nodes_1{1, 2, 3, 4};
    test_model_part.CreateNewElement("Element2D4N", 1, elem_nodes_1, test_model_part.pGetProperties(0));
    std::vector<ModelPart::IndexType> elem_nodes_2{2, 5, 6, 3};
    test_model_part.CreateNewElement("Element2D4N", 2, elem_nodes_2, test_model_part.pGetProperties(0));

    const auto it_elem_begin = test_model_part.ElementsBegin();
    auto it_elem = it_elem_begin + 0;
    array_1d<double, 3> forces {12.0, 8.0, 0.0};
    it_elem->SetValue(FORCE, forces);

    it_elem = it_elem_begin + 1;
    array_1d<double, 3> forces_2 {16.0, 4.0, 0.0};
    it_elem->SetValue(FORCE, forces_2);

    ConversionUtilities::ConvertElementalDataToNodalData(test_model_part);
    Vector nodal_force_values(18);
    int counter = 0;

    for (auto& r_node : test_model_part.Nodes()){
        auto nodal_force = r_node.FastGetSolutionStepValue(FORCE);
        for(int i = 0 ; i< 3; i++)
        {
            nodal_force_values[counter] = nodal_force[i];
            counter ++;
        }
    }

    array_1d<double, 18> expected_values {3, 2 ,0, 7, 3, 0, 7, 3 ,0, 3, 2 ,0, 4, 1 ,0, 4, 1 ,0};

    KRATOS_CHECK_VECTOR_NEAR(nodal_force_values, expected_values, 1.0e-4)

}


} // namespace Testing
} // namespace Kratos
