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
#include "containers/csr_matrix.h"
#include "containers/nd_data.h"
#include "containers/sparse_contiguous_row_graph.h"
#include "testing/testing.h"
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
    DofArrayUtilities::SetUpDofArray(r_test_model_part, dof_set);

    // Build the sparse matrix graph
    SparseContiguousRowGraph<std::size_t> sparse_matrix_graph(dof_set.size());
    IndexPartition<std::size_t>(r_test_model_part.NumberOfElements()).for_each([&](std::size_t Index) {
        Element::EquationIdVectorType eq_ids;
        auto it_elem = r_test_model_part.ElementsBegin() + Index;
        it_elem->EquationIdVector(eq_ids, r_test_model_part.GetProcessInfo());
        sparse_matrix_graph.AddEntries(eq_ids);
    });

    // Allocate the CSR matrix
    CsrMatrix<> csr_matrix(sparse_matrix_graph);

    // Call the utility
    DenseVector<unsigned int> aux_size(0);
    NDData<unsigned int> eq_ids_csr_indices(aux_size);
    CsrUtilities::GetEquationIdCsrIndices(r_test_model_part.Elements(), r_test_model_part.GetProcessInfo(), csr_matrix, eq_ids_csr_indices);

    // // Verify results
    // const auto& shape = nd_data.Shape();
    // KRATOS_EXPECT_EQ(shape.size(), 3);
    // KRATOS_EXPECT_EQ(shape[0], 2); // 2 entities
    // KRATOS_EXPECT_EQ(shape[1], 3); // 3 DoFs per entity
    // KRATOS_EXPECT_EQ(shape[2], 3); // 3 DoFs per entity
}

} // namespace Kratos::Testing
