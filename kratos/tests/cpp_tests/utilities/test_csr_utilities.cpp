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
#include "testing/testing.h"
#include "utilities/csr_utilities.h"

namespace Kratos::Testing
{

namespace
{
// USE THE ONES FROM THE FUTURE!!!!!
}

KRATOS_TEST_CASE_IN_SUITE(CsrUtilitiesGetCsrEquationIdIndices, KratosCoreFastSuite)
{
    // Create a container of mock entities
    std::vector<MockEntity> container;
    container.emplace_back(3); // Entity with 3 DoFs
    container.emplace_back(3); // Entity with 3 DoFs

    // Create a dummy CsrMatrix (type doesn't matter for this test as usage is not yet implemented in utility)
    int dummy_matrix = 0;

    // Output NdData
    NdData<std::size_t> nd_data(DenseVector<unsigned int>(0));

    // Call the utility
    CsrUtilities::GetCsrEquationIdIndices(container, dummy_matrix, nd_data);

    // Verify results
    const auto& shape = nd_data.Shape();
    KRATOS_EXPECT_EQ(shape.size(), 3);
    KRATOS_EXPECT_EQ(shape[0], 2); // 2 entities
    KRATOS_EXPECT_EQ(shape[1], 3); // 3 DoFs per entity
    KRATOS_EXPECT_EQ(shape[2], 3); // 3 DoFs per entity
}

} // namespace Kratos::Testing
