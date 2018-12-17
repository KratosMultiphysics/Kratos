//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//
//

#include "includes/kratos_flags.h"
#include "mpi/includes/mpi_data_communicator.h"

#include "testing/testing.h"

namespace Kratos {

namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommunicatorFlagsAndAll, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_self_communicator(MPI_COMM_SELF);
    //const int rank = mpi_self_communicator.Rank();
    //const int size = mpi_self_communicator.Size();

    Kratos::Flags test_flag;
    test_flag.Set(STRUCTURE, true);

    Kratos::Flags out = mpi_self_communicator.AndAll(test_flag, STRUCTURE);

    KRATOS_CHECK_EQUAL(test_flag.Is(STRUCTURE), true);
    KRATOS_CHECK_EQUAL(test_flag.IsDefined(PERIODIC), false);
}

}

}