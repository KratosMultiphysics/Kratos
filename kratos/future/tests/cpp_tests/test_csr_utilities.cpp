//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

// System includes

// External includes

// Project includes
#include "containers/nd_data.h"
#include "testing/testing.h"
#include "future/containers/define_linear_algebra_serial.h"
#include "future/containers/linear_system_container.h"
#include "future/solving_strategies/builders/block_builder.h"
#include "future/utilities/csr_utilities.h"
#include "utilities/dof_utilities/dof_array_utilities.h"
#include "test_utilities/solving_strategies_test_utilities.h"

namespace Kratos::Testing
{

namespace
{
// USE THE ONES FROM THE FUTURE!!!!!
}

KRATOS_TEST_CASE_IN_SUITE(CsrUtilitiesGetCsrEquationIdIndices, KratosCoreFastSuite)
{
    // Set up the test model part
    Model test_model;
    auto& r_test_model_part = test_model.CreateModelPart("TestModelPart");
    const std::size_t num_elems_x = 5;
    const std::size_t num_elems_y = 3;
    const double elem_size_x = 1.0;
    const double elem_size_y = 1.0;
    SolvingStrategiesTestUtilities::SetUpTestModelPart2D(num_elems_x, num_elems_y, elem_size_x, elem_size_y, r_test_model_part);

    // Set the DOF array
    ModelPart::DofsArrayType dof_set;
    DofArrayUtilities::SetUpDofArray(r_model_part, dof_set);

    // // Block builder set up to allocate the linear system matrix
    // using BuilderType = BlockBuilder<SerialLinearAlgebraTraits::MatrixType, SerialLinearAlgebraTraits::VectorType, SerialLinearAlgebraTraits::GraphType>;
    // auto p_builder = Kratos::make_shared<BuilderType>(r_test_model_part);
    // //FIXME: avoid builders

    // Output NdData
    NDDData<std::size_t> nd_data(DenseVector<unsigned int>(0));

    // Call the utility
    // CsrUtilities::GetCsrEquationIdIndices(container, dummy_matrix, nd_data);

    // // Verify results
    // const auto& shape = nd_data.Shape();
    // KRATOS_EXPECT_EQ(shape.size(), 3);
    // KRATOS_EXPECT_EQ(shape[0], 2); // 2 entities
    // KRATOS_EXPECT_EQ(shape[1], 3); // 3 DoFs per entity
    // KRATOS_EXPECT_EQ(shape[2], 3); // 3 DoFs per entity
}

} // namespace Kratos::Testing
