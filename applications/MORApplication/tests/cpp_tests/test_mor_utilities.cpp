//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Quirin Aumann
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
// #include "includes/ublas_interface.h"
#include "includes/ublas_complex_interface.h"
#include "mor_application_variables.h"
#include "custom_utilities/complex_dof_updater.hpp"
#include "containers/model.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(ComplexDofUpdaterTest, KratosMORFastSuite)
{
    using complex = std::complex<double>;

    // create model part
    Model current_model;

    // Generate a model part with the previous
    ModelPart& model_part = current_model.CreateModelPart("test_model_part");
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(REAL_DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(IMAG_DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(PRESSURE);
    model_part.AddNodalSolutionStepVariable(REAL_PRESSURE);
    model_part.AddNodalSolutionStepVariable(IMAG_PRESSURE);

    // Fill the model part geometry data
    auto p_node_1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    auto p_node_2 = model_part.CreateNewNode(2, 1.0, 1.0, 1.0);

    // Assign dofs and equation ids
    auto dof_1 = p_node_1->pAddDof(DISPLACEMENT_X);
    dof_1->SetEquationId(0);
    auto dof_2 = p_node_1->pAddDof(DISPLACEMENT_Y);
    dof_2->SetEquationId(1);
    auto dof_3 = p_node_1->pAddDof(DISPLACEMENT_Z);
    dof_3->SetEquationId(2);
    auto dof_4 = p_node_1->pAddDof(PRESSURE);
    dof_4->SetEquationId(3);
    auto dof_5 = p_node_2->pAddDof(PRESSURE);
    dof_5->SetEquationId(4);

    // Generate test data
    ComplexVector X = ComplexVector(5);
    for( size_t i=0; i<X.size(); ++i)
    {
        X(i) = complex(i, -1*(1+(int)i));
    }

    // Assign dof values
    ComplexDofUpdater::AssignDofs<ComplexVector>(model_part, X);

    // Test
    KRATOS_CHECK_DOUBLE_EQUAL(p_node_1->FastGetSolutionStepValue(DISPLACEMENT_Y), std::abs(X(1)));
    KRATOS_CHECK_DOUBLE_EQUAL(p_node_1->FastGetSolutionStepValue(REAL_DISPLACEMENT_Y), std::real(X(1)));
    KRATOS_CHECK_DOUBLE_EQUAL(p_node_1->FastGetSolutionStepValue(IMAG_DISPLACEMENT_Y), std::imag(X(1)));
    KRATOS_CHECK_DOUBLE_EQUAL(p_node_1->FastGetSolutionStepValue(IMAG_DISPLACEMENT_Z), std::imag(X(2)));
    KRATOS_CHECK_DOUBLE_EQUAL(p_node_2->FastGetSolutionStepValue(PRESSURE), std::abs(X(4)));
    KRATOS_CHECK_DOUBLE_EQUAL(p_node_2->FastGetSolutionStepValue(REAL_PRESSURE), std::real(X(4)));
    KRATOS_CHECK_DOUBLE_EQUAL(p_node_2->FastGetSolutionStepValue(IMAG_PRESSURE), std::imag(X(4)));
}

} // namespace Testing
} // namespace Kratos
