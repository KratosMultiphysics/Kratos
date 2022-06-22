//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//
//

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "includes/model_part.h"
#include "includes/variables.h"

// Application includes
#include "custom_utilities/indirect_variable.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(IndirectVariable, FluidDynamicsApplicationFastSuite) {
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("TestPart");
    r_model_part.SetBufferSize(2);
    r_model_part.AddNodalSolutionStepVariable(STEP);

    // non const check
    auto& r_node_1 = *r_model_part.CreateNewNode(1, 1.0, 0.0, 0.0);

    // const check
    const auto& r_const_node_1 = r_node_1;

    r_node_1.FastGetSolutionStepValue(STEP) = 1;
    r_node_1.FastGetSolutionStepValue(STEP, 1) = -1;

    IndirectVariable<int> indirect_step(STEP);
    IndirectVariable<int> indirect_none;

    indirect_step(r_node_1, 0) += 10;
    indirect_step(r_node_1, 1) += 10;
    indirect_none(r_node_1, 0) += 10;
    indirect_none(r_node_1, 1) += 10;

    KRATOS_CHECK_EQUAL(indirect_step(r_node_1, 0), 11);
    KRATOS_CHECK_EQUAL(indirect_step(r_node_1, 1), 9);
    KRATOS_CHECK_EQUAL(indirect_none(r_node_1, 0), 0);
    KRATOS_CHECK_EQUAL(indirect_none(r_node_1, 1), 0);

    KRATOS_CHECK_EQUAL(indirect_step(r_const_node_1, 0), 11);
    KRATOS_CHECK_EQUAL(indirect_step(r_const_node_1, 1), 9);
    KRATOS_CHECK_EQUAL(indirect_none(r_const_node_1, 0), 0);
    KRATOS_CHECK_EQUAL(indirect_none(r_const_node_1, 1), 0);

    indirect_step(r_node_1) += 10;
    indirect_none(r_node_1) += 10;

    KRATOS_CHECK_EQUAL(indirect_step(r_node_1), 21);
    KRATOS_CHECK_EQUAL(indirect_step(r_node_1, 1), 9);
    KRATOS_CHECK_EQUAL(indirect_none(r_node_1), 0);
    KRATOS_CHECK_EQUAL(indirect_none(r_node_1, 1), 0);

    KRATOS_CHECK_EQUAL(indirect_step(r_const_node_1), 21);
    KRATOS_CHECK_EQUAL(indirect_step(r_const_node_1, 1), 9);
    KRATOS_CHECK_EQUAL(indirect_none(r_const_node_1), 0);
    KRATOS_CHECK_EQUAL(indirect_none(r_const_node_1, 1), 0);

    KRATOS_CHECK_EQUAL(indirect_none(r_node_1, 1) * indirect_step(r_node_1), 0);

    KRATOS_CHECK_EQUAL(r_node_1.FastGetSolutionStepValue(STEP), 21);
    KRATOS_CHECK_EQUAL(r_node_1.FastGetSolutionStepValue(STEP, 1), 9);

}

}
}
