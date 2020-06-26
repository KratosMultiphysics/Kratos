//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//
//

#include "mpi.h"
#include "includes/parallel_environment.h"
#include "utilities/communication_coloring_utilities.h"

#include "testing/testing.h"

namespace Kratos {

namespace Testing {

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIColoringUtilities_ComputeRecvList, KratosMPICoreFastSuite)
{
    DataCommunicator& r_default_comm = ParallelEnvironment::GetDefaultDataCommunicator();

    const int world_size = r_default_comm.Size();
    const int current_rank = r_default_comm.Rank();

    KRATOS_SKIP_TEST_IF_NOT(world_size == 4) << "This test is designed for 4 MPI ranks." << std::endl;

    // send lists
    std::vector< std::vector< int > > send_list(4);
    send_list[0] = {1,3};
    send_list[1] = {0,2,3};
    send_list[2]; //does not send to anyone!!
    send_list[3] = {0};

    // //expected_recv_list;
    std::vector< std::vector< int > > expected_recv_list(4);
    expected_recv_list[0] = {1,3};
    expected_recv_list[1] = {0};
    expected_recv_list[2] = {1};
    expected_recv_list[3] = {0,1};

    auto recv_list = MPIColoringUtilities::ComputeRecvList(send_list[current_rank], r_default_comm);

    for(unsigned int j=0; j<recv_list.size(); ++j)
    {
        KRATOS_CHECK_EQUAL(recv_list[j], expected_recv_list[current_rank][j]);
    }
};

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIColoringUtilities_ComputeCommunicationScheduling, KratosMPICoreFastSuite)
{
    DataCommunicator& r_default_comm = ParallelEnvironment::GetDefaultDataCommunicator();

    const int world_size = r_default_comm.Size();
    const int current_rank = r_default_comm.Rank();

    KRATOS_SKIP_TEST_IF_NOT(world_size == 4) << "This test is designed for 4 MPI ranks." << std::endl;

    // send lists
    std::vector< std::vector< int > > send_list(4);
    send_list[0] = {1,3};
    send_list[1] = {0,2,3};
    send_list[2]; //does not send to anyone!!
    send_list[3] = {0};

    // //expected_recv_list;
    std::vector< std::vector< int > > expected_colors(4);
    expected_colors[0] = {1,3,-1};
    expected_colors[1] = {0,2,3};
    expected_colors[2] = {-1,1,-1};
    expected_colors[3] = {-1,0,1};

    auto colors = MPIColoringUtilities::ComputeCommunicationScheduling(send_list[current_rank], r_default_comm);

    for(unsigned int j=0; j<colors.size(); ++j)
    {
        KRATOS_CHECK_EQUAL(colors[j], expected_colors[current_rank][j]);
    }
};


}
}