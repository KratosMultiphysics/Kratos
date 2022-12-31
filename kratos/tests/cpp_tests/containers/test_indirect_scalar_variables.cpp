//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//
//

// system includes
#include <vector>

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "containers/indirect_scalar_variable.h"

// Application includes

namespace Kratos {

namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(IndirectScalarVariable, KratosCoreFastSuite) {
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("TestPart");
    r_model_part.SetBufferSize(2);
    r_model_part.AddNodalSolutionStepVariable(PRESSURE);

    // non const check
    auto& r_node_1 = *r_model_part.CreateNewNode(1, 1.0, 0.0, 0.0);

    // const check
    const auto& r_const_node_1 = r_node_1;

    // creating indirect variables
    IndirectScalarVariable indirect_pressure(PRESSURE);
    IndirectScalarVariable indirect_zero;

    indirect_pressure(r_node_1, 0) = 1.0;
    indirect_pressure(r_node_1, 1) = -1.0;
    indirect_zero(r_node_1, 0) = 1.0;
    indirect_zero(r_node_1, 1) = -1.0;

    indirect_pressure(r_node_1, 0) += 10.0;
    indirect_pressure(r_node_1, 1) += 10.0;
    indirect_zero(r_node_1, 0) += 10.0;
    indirect_zero(r_node_1, 1) += 10.0;

    indirect_pressure(r_node_1, 0) += indirect_zero(r_node_1, 0);
    indirect_pressure(r_node_1, 1) += indirect_zero(r_node_1, 0);
    indirect_zero(r_node_1, 0) += indirect_zero(r_node_1, 0);
    indirect_zero(r_node_1, 1) += indirect_zero(r_node_1, 0);

    indirect_pressure(r_node_1, 0) += indirect_pressure(r_node_1, 0);
    indirect_pressure(r_node_1, 1) += indirect_pressure(r_node_1, 0);
    indirect_zero(r_node_1, 0) += indirect_pressure(r_node_1, 0);
    indirect_zero(r_node_1, 1) += indirect_pressure(r_node_1, 0);

    KRATOS_CHECK_EQUAL(indirect_pressure(r_node_1, 0), 22.0);
    KRATOS_CHECK_EQUAL(indirect_pressure(r_node_1, 1), 31.0);
    KRATOS_CHECK_EQUAL(indirect_zero(r_node_1, 0), 0.0);
    KRATOS_CHECK_EQUAL(indirect_zero(r_node_1, 1), 0.0);

    KRATOS_CHECK_EQUAL(indirect_pressure(r_const_node_1), 22.0);
    KRATOS_CHECK_EQUAL(indirect_pressure(r_const_node_1, 1), 31.0);
    KRATOS_CHECK_EQUAL(indirect_zero(r_const_node_1), 0.0);
    KRATOS_CHECK_EQUAL(indirect_zero(r_const_node_1, 1), 0.0);

    KRATOS_CHECK_EQUAL(indirect_zero(r_node_1, 1) * indirect_pressure(r_node_1), 0.0);
    KRATOS_CHECK_EQUAL(indirect_pressure(r_node_1, 1) * indirect_pressure(r_node_1), 682.0);

    KRATOS_CHECK_EQUAL(r_node_1.FastGetSolutionStepValue(PRESSURE), 22.0);
    KRATOS_CHECK_EQUAL(r_node_1.FastGetSolutionStepValue(PRESSURE, 1), 31.0);
}

}
}
