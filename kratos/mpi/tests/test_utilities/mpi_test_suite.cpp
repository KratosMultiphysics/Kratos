//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos A. Roig
//
//

// System includes

// External includes
#include "mpi.h"

// Project includes
#include "testing/testing.h"

// Parallel Extension
#include "mpi/tests/test_utilities/mpi_test_suite.h"

namespace Kratos::Testing
{

void KratosMPICoreFastSuite::TearDown() {
    int has_failure = ::testing::Test::HasFailure();

    // Synchronize the failre status
    MPI_Allreduce(MPI_IN_PLACE, &has_failure, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    // If any of the processes has issued a failure, fail the whole tests
    if (has_failure) {
        KRATOS_FAIL();
    }
}

} // namespace Kratos::Testing