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

#include "includes/data_communicator.h"
#include "mpi/includes/mpi_data_communicator.h"

#include "testing/testing.h"

namespace Kratos {

namespace Testing {

// MPI Communicator splitting /////////////////////////////////////////////////

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorSplit, KratosMPICoreFastSuite)
{
    const DataCommunicator& r_comm = DataCommunicator::GetDefault();
    const int global_rank = r_comm.Rank();
    const int global_size = r_comm.Size();

    for (int i = 1; i < global_size; i++)
    {
        int color = global_rank < i ? 0 : 1;
        int key = global_rank < i ? global_rank : global_size - global_rank;
        std::stringstream name;
        name << "split_communicator_step_" << i;

        const DataCommunicator& r_split_comm = MPIDataCommunicator::SplitDataCommunicator(r_comm, color, key, name.str());

        int expected_size = global_rank < i ? i : global_size - i;
        int expected_rank = global_rank < i ? key : key - 1;
        // MPI_Comm_split assigns rank by order of increasing key.
        // Here keys in the global_rank >= i side of the split range between 1 and (global_size-i)
        KRATOS_CHECK_EQUAL(r_split_comm.Size(), expected_size);
        KRATOS_CHECK_EQUAL(r_split_comm.Rank(), expected_rank);
    }
}

}
}