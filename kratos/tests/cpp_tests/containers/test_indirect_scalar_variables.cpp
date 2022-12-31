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
#include "utilities/parallel_utilities.h"

// Application includes

namespace Kratos {

namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(IndirectScalarVariable, KratosCoreFastSuite) {
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("TestPart");
    r_model_part.SetBufferSize(2);
    r_model_part.AddNodalSolutionStepVariable(PRESSURE);

    // non const check
    for (IndexType i = 1; i <= 100; ++i){
        r_model_part.CreateNewNode(i, 1.0, 0.0, 0.0);
    }

    // creating indirect variables
    IndirectScalarVariable indirect_pressure(PRESSURE);
    IndirectScalarVariable indirect_zero;

    block_for_each(r_model_part.Nodes(), [&](auto& rNode){
        indirect_pressure(rNode, 0) = 1.0;
        indirect_pressure(rNode, 1) = -1.0;
        indirect_zero(rNode, 0) = 1.0;
        indirect_zero(rNode, 1) = -1.0;

        indirect_pressure(rNode, 0) += 10.0;
        indirect_pressure(rNode, 1) += 10.0;
        indirect_zero(rNode, 0) += 10.0;
        indirect_zero(rNode, 1) += 10.0;

        indirect_pressure(rNode, 0) += indirect_zero(rNode, 0);
        indirect_pressure(rNode, 1) += indirect_zero(rNode, 0);
        indirect_zero(rNode, 0) += indirect_zero(rNode, 0);
        indirect_zero(rNode, 1) += indirect_zero(rNode, 0);

        indirect_pressure(rNode, 0) += indirect_pressure(rNode, 0);
        indirect_pressure(rNode, 1) += indirect_pressure(rNode, 0);
        indirect_zero(rNode, 0) += indirect_pressure(rNode, 0);
        indirect_zero(rNode, 1) += indirect_pressure(rNode, 0);
    });

    block_for_each(r_model_part.Nodes(), [&](auto& rNode) {
        KRATOS_CHECK_EQUAL(indirect_pressure(rNode, 0), 22.0);
        KRATOS_CHECK_EQUAL(indirect_pressure(rNode, 1), 31.0);
        KRATOS_CHECK_EQUAL(indirect_zero(rNode, 0), 0.0);
        KRATOS_CHECK_EQUAL(indirect_zero(rNode, 1), 0.0);
    });

    block_for_each(r_model_part.Nodes(), [&](const auto& rNode){
        KRATOS_CHECK_EQUAL(indirect_pressure(rNode), 22.0);
        KRATOS_CHECK_EQUAL(indirect_pressure(rNode, 1), 31.0);
        KRATOS_CHECK_EQUAL(indirect_zero(rNode), 0.0);
        KRATOS_CHECK_EQUAL(indirect_zero(rNode, 1), 0.0);
    });

    block_for_each(r_model_part.Nodes(), [&](const auto& rNode){
        KRATOS_CHECK_EQUAL(indirect_zero(rNode, 1) * indirect_pressure(rNode), 0.0);
        KRATOS_CHECK_EQUAL(indirect_pressure(rNode, 1) * indirect_pressure(rNode), 682.0);
        KRATOS_CHECK_EQUAL(rNode.FastGetSolutionStepValue(PRESSURE), 22.0);
        KRATOS_CHECK_EQUAL(rNode.FastGetSolutionStepValue(PRESSURE, 1), 31.0);
    });
}

}
}
