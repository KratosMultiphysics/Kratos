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

#include "includes/parallel_environment.h"
#include "mpi/includes/mpi_communicator.h"
#include "testing/testing.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(MPICommunicatorCreation, KratosMPICoreFastSuite)
{
    const DataCommunicator& r_world = ParallelEnvironment::GetDataCommunicator("World");

    /*const MPICommunicator constructed_communicator = MPICommunicator(r_world);

    KRATOS_CHECK_EQUAL(constructed_communicator.IsDistributed(), true);
    KRATOS_CHECK_EQUAL(constructed_communicator.MyPID(), r_world.Rank());
    KRATOS_CHECK_EQUAL(constructed_communicator.TotalProcesses(), r_world.Size());

    Communicator::Pointer p_created_communicator = constructed_communicator.Create(r_world);

    KRATOS_CHECK_EQUAL(p_created_communicator->IsDistributed(), true);
    KRATOS_CHECK_EQUAL(p_created_communicator->MyPID(), r_world.Rank());
    KRATOS_CHECK_EQUAL(p_created_communicator->TotalProcesses(), r_world.Size());*/
}

}
}