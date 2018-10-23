//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Philipp Bucher
//

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "includes/model_part.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(ModelPartSubModelPartsIterator, KratosCoreFastSuite)
{
    Model current_model;

    ModelPart& model_part = current_model.CreateModelPart("Main");

    model_part.CreateSubModelPart("Inlet1");
    model_part.CreateSubModelPart("Inlet2");
    model_part.CreateSubModelPart("Outlet");
    model_part.CreateSubModelPart("AnotherOutlet");

    std::size_t id = 1;
    for(auto i_SubModelPart = model_part.SubModelPartsBegin() ; i_SubModelPart != model_part.SubModelPartsEnd() ; i_SubModelPart++){
        i_SubModelPart->CreateNewNode(id++, 0.00,0.00,0.00);
    }

    KRATOS_CHECK_EQUAL(model_part.NumberOfNodes(), 4);
    KRATOS_CHECK_EQUAL(model_part.GetSubModelPart("Inlet1").NumberOfNodes(), 1);
    KRATOS_CHECK_EQUAL(model_part.GetSubModelPart("Outlet").NumberOfNodes(), 1);
}

KRATOS_TEST_CASE_IN_SUITE(ModelPartAddNodalSolutionStepVariable, KratosCoreFastSuite)
{
    Model current_model;

    ModelPart& model_part = current_model.CreateModelPart("Main");

    model_part.AddNodalSolutionStepVariable(VELOCITY);

    model_part.CreateNewNode(123, 0.00,0.00,0.00);

    KRATOS_CHECK_EXCEPTION_IS_THROWN(model_part.AddNodalSolutionStepVariable(PRESSURE),
        "Error: Attempting to add the variable \"PRESSURE\" to the model part with name \"Main\" which is not empty");

    model_part.AddNodalSolutionStepVariable(VELOCITY); // Adding the same Variable twice is fine bcs it wont do anything
}

KRATOS_TEST_CASE_IN_SUITE(ModelPartHasNodalSolutionStepVariable, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& model_part = current_model.CreateModelPart("Main");

    model_part.AddNodalSolutionStepVariable(VELOCITY);

    KRATOS_CHECK(model_part.HasNodalSolutionStepVariable(VELOCITY));
    KRATOS_CHECK_IS_FALSE(model_part.HasNodalSolutionStepVariable(PRESSURE));
}

KRATOS_TEST_CASE_IN_SUITE(ModelPartEmptyName, KratosCoreFastSuite)
{
    Model current_model;

    // Constructor with name
    KRATOS_CHECK_EXCEPTION_IS_THROWN(current_model.CreateModelPart(""),
        "Error: Please don't use empty names (\"\") when creating a ModelPart");

    // Constructor with name and bufferSize
    KRATOS_CHECK_EXCEPTION_IS_THROWN(current_model.CreateModelPart("", 2),
        "Error: Please don't use empty names (\"\") when creating a ModelPart");
}

KRATOS_TEST_CASE_IN_SUITE(ModelPartNameContainingPoint, KratosCoreFastSuite)
{
    Model current_model;

    // Constructor with name
    KRATOS_CHECK_EXCEPTION_IS_THROWN(current_model.CreateModelPart("name.other"),
        "Error: Please don't use names containing (\".\") when creating a ModelPart");

    // Constructor with name and bufferSize
    KRATOS_CHECK_EXCEPTION_IS_THROWN(current_model.CreateModelPart("name.other", 2),
        "Error: Please don't use names containing (\".\") when creating a ModelPart");
}
}
}  // namespace Kratos.
