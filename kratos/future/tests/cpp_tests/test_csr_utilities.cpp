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
    const DenseVector<unsigned int> aux_size(0);
    NDData<int> eq_ids_csr_indices(aux_size);
    KRATOS_WATCH("A")
    CsrUtilities::GetEquationIdCsrIndices(r_test_model_part.Elements(), r_test_model_part.GetProcessInfo(), csr_matrix, eq_ids_csr_indices);
    KRATOS_WATCH("B")

    // Verify results
    const auto& r_shape = eq_ids_csr_indices.Shape();
    const auto& r_data = eq_ids_csr_indices.ViewData();
    KRATOS_EXPECT_EQ(r_shape.size(), 3);
    KRATOS_EXPECT_EQ(r_shape[0], 2);
    KRATOS_EXPECT_EQ(r_shape[1], 3);
    KRATOS_EXPECT_EQ(r_shape[2], 3);
    KRATOS_EXPECT_EQ(r_data[0], 0);
    KRATOS_EXPECT_EQ(r_data[1], 1);
    KRATOS_EXPECT_EQ(r_data[2], 2);
    KRATOS_EXPECT_EQ(r_data[3], 3);
    KRATOS_EXPECT_EQ(r_data[4], 4);
    KRATOS_EXPECT_EQ(r_data[5], 5);
}

KRATOS_TEST_CASE_IN_SUITE(CsrUtilitiesAssemble, KratosCoreFastSuite)
{
    // Set up the indices array
    DenseVector<unsigned int> eq_ids_shape(3);
    eq_ids_shape[0] = 3;
    eq_ids_shape[1] = 2;
    eq_ids_shape[2] = 2;
    int eq_ids_csr_indices_data[12];
    eq_ids_csr_indices_data[0] = 0; eq_ids_csr_indices_data[1] = 1;
    eq_ids_csr_indices_data[2] = 2; eq_ids_csr_indices_data[3] = 3;
    eq_ids_csr_indices_data[4] = 3; eq_ids_csr_indices_data[5] = 4;
    eq_ids_csr_indices_data[6] = 5; eq_ids_csr_indices_data[7] = 6;
    eq_ids_csr_indices_data[8] = 6; eq_ids_csr_indices_data[9] = 7;
    eq_ids_csr_indices_data[10] = 8; eq_ids_csr_indices_data[11] = 9;
    NDData<int> eq_ids_csr_indices(eq_ids_csr_indices_data, eq_ids_shape);

    // Set up the contributions array
    DenseVector<unsigned int> contributions_shape(3);
    contributions_shape[0] = 3;
    contributions_shape[1] = 2;
    contributions_shape[2] = 2;
    double contributions_data[12];
    contributions_data[0] = 1.0; contributions_data[1] = -1.0;
    contributions_data[2] = -1.0; contributions_data[3] = 1.0;
    contributions_data[4] = 2.0; contributions_data[5] = -2.0;
    contributions_data[6] = -2.0; contributions_data[7] = 2.0;
    contributions_data[8] = 1.0; contributions_data[9] = -1.0;
    contributions_data[10] = -1.0; contributions_data[11] = 1.0;
    NDData<double> contributions(contributions_data, contributions_shape);

    // Set up and initialize the CSR matrix
    SparseContiguousRowGraph<std::size_t> sparse_matrix_graph(4);
    sparse_matrix_graph.AddEntry(0,0); sparse_matrix_graph.AddEntry(0,1);
    sparse_matrix_graph.AddEntry(1,0); sparse_matrix_graph.AddEntry(1,1); sparse_matrix_graph.AddEntry(1,2);
    sparse_matrix_graph.AddEntry(2,1); sparse_matrix_graph.AddEntry(2,2); sparse_matrix_graph.AddEntry(2,3);
    sparse_matrix_graph.AddEntry(3,2); sparse_matrix_graph.AddEntry(3,3);
    CsrMatrix<> csr_matrix(sparse_matrix_graph);
    csr_matrix.SetValue(0.0);

    // Call the utility
    CsrUtilities::Assemble(contributions, eq_ids_csr_indices, csr_matrix);

    // Verify results
}

} // namespace Kratos::Testing
