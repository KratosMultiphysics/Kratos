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
#include "mpi/mpi_environment.h"
#include "mpi/utilities/mpi_coloring_utilities.h"
#include "includes/parallel_environment.h"

#include "testing/testing.h"

namespace Kratos {

namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(MPIColoringUtilities_ComputeRecvList, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);

    const int world_size = mpi_world_communicator.Size();
    const int current_rank = mpi_world_communicator.Rank();

    if(world_size == 4) //only implemented for the case of 4 mpi ranks
    {
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

        auto recv_list = MPIColoringUtilities::ComputeRecvList(send_list[current_rank], mpi_world_communicator);

        for(unsigned int j=0; j<recv_list.size(); ++j)
        {
           KRATOS_CHECK_EQUAL(recv_list[j], expected_recv_list[current_rank][j]);
        }
    }
};

KRATOS_TEST_CASE_IN_SUITE(MPIColoringUtilities_ComputeCommunicationScheduling, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);

    const int world_size = mpi_world_communicator.Size();
    const int current_rank = mpi_world_communicator.Rank();

    if(world_size == 4) //only implemented for the case of 4 mpi ranks
    {
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

        auto colors = MPIColoringUtilities::ComputeCommunicationScheduling(send_list[current_rank], mpi_world_communicator);

        for(unsigned int j=0; j<colors.size(); ++j)
        {
           KRATOS_CHECK_EQUAL(colors[j], expected_colors[current_rank][j]);
        }
    }
};


}
}