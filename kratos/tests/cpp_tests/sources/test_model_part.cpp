//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Philipp Bucher
//                   Vicente Mataix Ferrandiz
//

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "includes/model_part.h"
#include "utilities/auxiliar_model_part_utilities.h"
#include "utilities/cpp_tests_utilities.h"

namespace Kratos {
  namespace Testing {

    typedef Node<3> NodeType;

    void GenerateGenericModelPart(ModelPart& rModelPart)
    {
        Properties::Pointer p_elem_prop = rModelPart.CreateNewProperties(0);

        CppTestsUtilities::Create2DGeometry(rModelPart, "Element2D3N");

        rModelPart.CreateNewCondition("LineCondition2D2N", 1, {{1,2}}, p_elem_prop);
        rModelPart.CreateNewCondition("LineCondition2D2N", 2, {{1,4}}, p_elem_prop);
        rModelPart.CreateNewCondition("LineCondition2D2N", 3, {{2,5}}, p_elem_prop);
        rModelPart.CreateNewCondition("LineCondition2D2N", 4, {{5,6}}, p_elem_prop);

    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartSubModelPartsIterator, KratosCoreFastSuite)
    {
        Model current_model;

        ModelPart& r_model_part = current_model.CreateModelPart("Main");

        r_model_part.CreateSubModelPart("Inlet1");
        r_model_part.CreateSubModelPart("Inlet2");
        r_model_part.CreateSubModelPart("Outlet");
        r_model_part.CreateSubModelPart("AnotherOutlet");

        std::size_t id = 1;
        for(auto i_SubModelPart = r_model_part.SubModelPartsBegin() ; i_SubModelPart != r_model_part.SubModelPartsEnd() ; i_SubModelPart++){
            i_SubModelPart->CreateNewNode(id++, 0.00,0.00,0.00);
        }

        KRATOS_CHECK_EQUAL(r_model_part.NumberOfNodes(), 4);
        KRATOS_CHECK_EQUAL(r_model_part.GetSubModelPart("Inlet1").NumberOfNodes(), 1);
        KRATOS_CHECK_EQUAL(r_model_part.GetSubModelPart("Outlet").NumberOfNodes(), 1);
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartAddNodalSolutionStepVariable, KratosCoreFastSuite)
    {
        Model current_model;

        ModelPart& r_model_part = current_model.CreateModelPart("Main");

        r_model_part.AddNodalSolutionStepVariable(VELOCITY);

        r_model_part.CreateNewNode(123, 0.00,0.00,0.00);

        KRATOS_CHECK_EXCEPTION_IS_THROWN(r_model_part.AddNodalSolutionStepVariable(PRESSURE),
            "Error: Attempting to add the variable \"PRESSURE\" to the model part with name \"Main\" which is not empty");

        r_model_part.AddNodalSolutionStepVariable(VELOCITY); // Adding the same Variable twice is fine bcs it wont do anything
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartHasNodalSolutionStepVariable, KratosCoreFastSuite)
    {
        Model current_model;
        ModelPart& r_model_part = current_model.CreateModelPart("Main");


        r_model_part.AddNodalSolutionStepVariable(VELOCITY);

        KRATOS_CHECK(r_model_part.HasNodalSolutionStepVariable(VELOCITY));
        KRATOS_CHECK_IS_FALSE(r_model_part.HasNodalSolutionStepVariable(PRESSURE));
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartFullName, KratosCoreFastSuite)
    {
        Model current_model;
        ModelPart& r_model_part = current_model.CreateModelPart("Main");
        ModelPart& r_sub_model_part = r_model_part.CreateSubModelPart("SubModelPart");
        ModelPart& r_sub_sub_model_part = r_sub_model_part.CreateSubModelPart("SubSubModelPart");

        KRATOS_CHECK_STRING_EQUAL(r_model_part.FullName(), "Main");
        KRATOS_CHECK_STRING_EQUAL(r_sub_model_part.FullName(), "Main.SubModelPart");
        KRATOS_CHECK_STRING_EQUAL(r_sub_sub_model_part.FullName(), "Main.SubModelPart.SubSubModelPart");
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartEmptyName, KratosCoreFastSuite)
    {
        Model current_model;

        // Constructor with name
        KRATOS_CHECK_EXCEPTION_IS_THROWN(current_model.CreateModelPart(""),
            "Error: Please don't use empty names (\"\") when creating a ModelPart");
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartNameContainingPoint, KratosCoreFastSuite)
    {
        Model current_model;

        // Constructor with name
        KRATOS_CHECK_EXCEPTION_IS_THROWN(current_model.CreateModelPart("name.other"),
            "Error: Please don't use names containing (\".\") when creating a ModelPart (used in \"name.other\")");
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartRemoveElements, KratosCoreFastSuite)
    {
        Model current_model;

        ModelPart& r_model_part = current_model.CreateModelPart("Main");

        // Fill model part
        GenerateGenericModelPart(r_model_part);

        // Set flags
        r_model_part.pGetElement(1)->Set(TO_ERASE, true);
        r_model_part.pGetElement(2)->Set(TO_ERASE, true);

        // Call method
        r_model_part.RemoveElements(TO_ERASE);

        // Check results
        KRATOS_CHECK(r_model_part.NumberOfNodes() == 6);
        KRATOS_CHECK(r_model_part.NumberOfElements() == 2);
        KRATOS_CHECK(r_model_part.NumberOfConditions() == 4);
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartRemoveElementsAndBelongings, KratosCoreFastSuite)
    {
        Model current_model;

        ModelPart& r_model_part = current_model.CreateModelPart("Main");

        // Fill model part
        GenerateGenericModelPart(r_model_part);

        // Set flags
        r_model_part.pGetElement(1)->Set(TO_ERASE, true);
        r_model_part.pGetElement(2)->Set(TO_ERASE, true);

        // Call method
        auto aux_util = AuxiliarModelPartUtilities(r_model_part);
        aux_util.RemoveElementsAndBelongings(TO_ERASE);

        // Check results
        KRATOS_CHECK(r_model_part.NumberOfNodes() == 4);
        KRATOS_CHECK(r_model_part.NumberOfElements() == 2);
        KRATOS_CHECK(r_model_part.NumberOfConditions() == 2);
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartRemoveConditions, KratosCoreFastSuite)
    {
        Model current_model;

        ModelPart& r_model_part = current_model.CreateModelPart("Main");

        // Fill model part
        GenerateGenericModelPart(r_model_part);

        // Set flags
        r_model_part.pGetCondition(1)->Set(TO_ERASE, true);
        r_model_part.pGetCondition(2)->Set(TO_ERASE, true);

        // Call method
        r_model_part.RemoveConditions(TO_ERASE);

        // Check results
        KRATOS_CHECK(r_model_part.NumberOfNodes() == 6);
        KRATOS_CHECK(r_model_part.NumberOfConditions() == 2);
        KRATOS_CHECK(r_model_part.NumberOfElements() == 4);
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartRemoveConditionsAndBelongings, KratosCoreFastSuite)
    {
        Model current_model;

        ModelPart& r_model_part = current_model.CreateModelPart("Main");

        // Fill model part
        GenerateGenericModelPart(r_model_part);

        // Set flags
        r_model_part.pGetCondition(1)->Set(TO_ERASE, true);
        r_model_part.pGetCondition(2)->Set(TO_ERASE, true);

        // Call method
        auto aux_util = AuxiliarModelPartUtilities(r_model_part);
        aux_util.RemoveConditionsAndBelongings(TO_ERASE);

        // Check results
        KRATOS_CHECK(r_model_part.NumberOfNodes() == 3);
        KRATOS_CHECK(r_model_part.NumberOfConditions() == 2);
        KRATOS_CHECK(r_model_part.NumberOfElements() == 2);
    }
  }  // namespace Testing.
}  // namespace Kratos.
