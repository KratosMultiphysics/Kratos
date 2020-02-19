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

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorFlagsAndAll, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int rank = mpi_world_communicator.Rank();
    const int size = mpi_world_communicator.Size();

    Kratos::Flags test_flag;
    test_flag.Set(STRUCTURE, rank == 0);
    test_flag.Set(INLET, rank == 0);

    Kratos::Flags synchronized_flag = mpi_world_communicator.AndReduceAll(test_flag, STRUCTURE);

    KRATOS_CHECK_EQUAL(synchronized_flag.Is(STRUCTURE), (size == 1)); // true for single-rank runs, false for multiple ranks.
    KRATOS_CHECK_EQUAL(synchronized_flag.Is(INLET), (rank == 0)); // This value does not participate in the synchronization
    KRATOS_CHECK_EQUAL(synchronized_flag.IsDefined(PERIODIC), false);
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorFlagsAndAllUnset, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int rank = mpi_world_communicator.Rank();
    const int size = mpi_world_communicator.Size();

    Kratos::Flags test_flag;
    if (rank != size - 1) // last rank does not define a value
    {
        test_flag.Set(STRUCTURE, true);
    }
    test_flag.Set(INLET, rank == 0);

    Kratos::Flags synchronized_flag = mpi_world_communicator.AndReduceAll(test_flag, STRUCTURE);

    if (size > 1)
    {
        KRATOS_CHECK_EQUAL(synchronized_flag.Is(STRUCTURE), false); // all ranks (including last) are set to false (since one rank was unset).
    }
    else
    {
        KRATOS_CHECK_EQUAL(synchronized_flag.IsDefined(STRUCTURE), false); // only one rank, which did not set it.
    }
    KRATOS_CHECK_EQUAL(synchronized_flag.Is(INLET), (rank == 0)); // This value does not participate in the synchronization
    KRATOS_CHECK_EQUAL(synchronized_flag.IsDefined(PERIODIC), false);
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorFlagsOrAll, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int rank = mpi_world_communicator.Rank();

    Kratos::Flags test_flag;
    test_flag.Set(STRUCTURE, rank == 0);
    test_flag.Set(INLET, rank == 0);

    Kratos::Flags synchronized_flag = mpi_world_communicator.OrReduceAll(test_flag, STRUCTURE);

    KRATOS_CHECK_EQUAL(synchronized_flag.Is(STRUCTURE), true);
    KRATOS_CHECK_EQUAL(synchronized_flag.Is(INLET), (rank == 0)); // This value does not participate in the synchronization
    KRATOS_CHECK_EQUAL(synchronized_flag.IsDefined(PERIODIC), false);
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorFlagsOrAllUnset, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int rank = mpi_world_communicator.Rank();
    const int size = mpi_world_communicator.Size();

    Kratos::Flags test_flag;
    if (rank != size - 1) // last rank does not define a value
    {
        test_flag.Set(STRUCTURE, true);
    }
    test_flag.Set(INLET, rank == 0);

    Kratos::Flags synchronized_flag = mpi_world_communicator.OrReduceAll(test_flag, STRUCTURE);

    if (size > 1)
    {
        KRATOS_CHECK_EQUAL(synchronized_flag.Is(STRUCTURE), true); // all ranks (including last) are set.
    }
    else
    {
        KRATOS_CHECK_EQUAL(synchronized_flag.IsDefined(STRUCTURE), false); // only one rank, which did not set it.
    }
    KRATOS_CHECK_EQUAL(synchronized_flag.Is(INLET), (rank == 0)); // This value does not participate in the synchronization
    KRATOS_CHECK_EQUAL(synchronized_flag.IsDefined(PERIODIC), false);
}

// Synchronization using ALL_DEFINED as mask //////////////////////////////////

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorFlagsReduceWithAllDefinedAsMask, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int rank = mpi_world_communicator.Rank();
    const int size = mpi_world_communicator.Size();

    Kratos::Flags test_flag;
    test_flag.Set(STRUCTURE, rank == 0);
    test_flag.Set(INLET, rank == 0);

    // Using ALL_DEFINED as argument should synchronize all defined flags
    Kratos::Flags synchronized_and_flag = mpi_world_communicator.AndReduceAll(test_flag, ALL_DEFINED);

    KRATOS_CHECK_EQUAL(synchronized_and_flag.Is(STRUCTURE), (size == 1)); // true for single-rank runs, false for multiple ranks.
    KRATOS_CHECK_EQUAL(synchronized_and_flag.Is(INLET), (size == 1));
    KRATOS_CHECK_EQUAL(synchronized_and_flag.IsDefined(PERIODIC), false);

    Kratos::Flags synchronized_or_flag = mpi_world_communicator.OrReduceAll(test_flag, ALL_DEFINED);

    KRATOS_CHECK_EQUAL(synchronized_or_flag.Is(STRUCTURE), true);
    KRATOS_CHECK_EQUAL(synchronized_or_flag.Is(INLET), true);
    KRATOS_CHECK_EQUAL(synchronized_or_flag.IsDefined(PERIODIC), false);
}

// Flags And //////////////////////////////////////////////////////////////////

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorFlagsAndOperations, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_rank = mpi_world_communicator.Rank();
    const int world_size = mpi_world_communicator.Size();
    constexpr int root = 0;

    Kratos::Flags flags;
    //       both true | both false | opposite sets | first true   | first false | second true | second false
    if (world_rank == root) {
        flags = ACTIVE | RIGID.AsFalse()  | STRUCTURE    | MPI_BOUNDARY | PERIODIC.AsFalse();
    }
    else {
        flags = ACTIVE | RIGID.AsFalse()  | STRUCTURE.AsFalse() |                              INLET       | OUTLET.AsFalse();
    }

    // Setting an extra flag, not involved in communication
    if (world_rank == root)
        flags.Set(CONTACT, true);

    Flags output = mpi_world_communicator.AndReduce(flags, ACTIVE | RIGID | STRUCTURE | MPI_BOUNDARY | PERIODIC | INLET | OUTLET | ISOLATED , root);

    if (world_size > 1 && world_rank == root) {
        // true (defined) & true (defined) = true (defined)
        KRATOS_CHECK_EQUAL(output.IsDefined(ACTIVE), true);
        KRATOS_CHECK_EQUAL(output.Is(ACTIVE), true);
        // false (defined) & false (defined) = false (defined)
        KRATOS_CHECK_EQUAL(output.IsDefined(RIGID), true);
        KRATOS_CHECK_EQUAL(output.Is(RIGID), false);
        // true (defined) & false (defined) = false (defined)
        KRATOS_CHECK_EQUAL(output.IsDefined(STRUCTURE), true);
        KRATOS_CHECK_EQUAL(output.Is(STRUCTURE), false);
        // true (defined) & (undefined) = false (defined) (undefined defaults to false)
        KRATOS_CHECK_EQUAL(output.IsDefined(MPI_BOUNDARY), true);
        KRATOS_CHECK_EQUAL(output.Is(MPI_BOUNDARY), false);
        // false (defined) & (undefined) = false (defined) (undefined defaults to false)
        KRATOS_CHECK_EQUAL(output.IsDefined(PERIODIC), true);
        KRATOS_CHECK_EQUAL(output.Is(PERIODIC), false);
        // (undefined) & true (defined) = false (defined) (undefined defaults to false)
        KRATOS_CHECK_EQUAL(output.IsDefined(INLET), true);
        KRATOS_CHECK_EQUAL(output.Is(INLET), false);
        // (undefied) & false (defined) = false (defined) (undefined defaults to false)
        KRATOS_CHECK_EQUAL(output.IsDefined(OUTLET), true);
        KRATOS_CHECK_EQUAL(output.Is(OUTLET), false);
        // (undefined) & (undefined) = (undefined)
        KRATOS_CHECK_EQUAL(output.IsDefined(ISOLATED), false);
        KRATOS_CHECK_EQUAL(output.Is(ISOLATED), false);

        // set, but not involved in communication
        KRATOS_CHECK_EQUAL(output.IsDefined(CONTACT), world_rank == root);
        KRATOS_CHECK_EQUAL(output.Is(CONTACT), world_rank == root);
        // not set and not involved in communication
        KRATOS_CHECK_EQUAL(output.IsDefined(TO_ERASE), false);
    }
    else {
        KRATOS_CHECK_EQUAL(output, flags);
    }
}

// Flags Or ///////////////////////////////////////////////////////////////////

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorFlagsOrOperations, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_rank = mpi_world_communicator.Rank();
    const int world_size = mpi_world_communicator.Size();
    constexpr int root = 0;

    Kratos::Flags flags;
    //       both true | both false | opposite sets | first true   | first false | second true | second false
    if (world_rank == root) {
        flags = ACTIVE | RIGID.AsFalse()  | STRUCTURE     | MPI_BOUNDARY | PERIODIC.AsFalse();
    }
    else {
        flags = ACTIVE | RIGID.AsFalse()  | STRUCTURE.AsFalse() |                              INLET       | OUTLET.AsFalse();
    }

    // Setting an extra flag, not involved in communication
    if (world_rank == root)
        flags.Set(CONTACT, true);

    Flags output = mpi_world_communicator.OrReduce(flags, ACTIVE | RIGID | STRUCTURE | MPI_BOUNDARY | PERIODIC | INLET | OUTLET | ISOLATED , root);

    if (world_size > 1 && world_rank == root) {
        // true (defined) | true (defined) = true (defined)
        KRATOS_CHECK_EQUAL(output.IsDefined(ACTIVE), true);
        KRATOS_CHECK_EQUAL(output.Is(ACTIVE), true);
        // false (defined) | false (defined) = false (defined)
        KRATOS_CHECK_EQUAL(output.IsDefined(RIGID), true);
        KRATOS_CHECK_EQUAL(output.Is(RIGID), false);
        // true (defined) | false (defined) = true (defined)
        KRATOS_CHECK_EQUAL(output.IsDefined(STRUCTURE), true);
        KRATOS_CHECK_EQUAL(output.Is(STRUCTURE), true);
        // true (defined) | (undefined) = true (defined) (undefined defaults to false)
        KRATOS_CHECK_EQUAL(output.IsDefined(MPI_BOUNDARY), true);
        KRATOS_CHECK_EQUAL(output.Is(MPI_BOUNDARY), true);
        // false (defined) | (undefined) = false (defined) (undefined defaults to false)
        KRATOS_CHECK_EQUAL(output.IsDefined(PERIODIC), true);
        KRATOS_CHECK_EQUAL(output.Is(PERIODIC), false);
        // (undefined) | true (defined) = true (defined) (undefined defaults to false)
        KRATOS_CHECK_EQUAL(output.IsDefined(INLET), true);
        KRATOS_CHECK_EQUAL(output.Is(INLET), true);
        // (undefied) | false (defined) = false (defined) (undefined defaults to false)
        KRATOS_CHECK_EQUAL(output.IsDefined(OUTLET), true);
        KRATOS_CHECK_EQUAL(output.Is(OUTLET), false);
        // (undefined) | (undefined) = (undefined)
        KRATOS_CHECK_EQUAL(output.IsDefined(ISOLATED), false);
        KRATOS_CHECK_EQUAL(output.Is(ISOLATED), false);

        // set, but not involved in communication
        KRATOS_CHECK_EQUAL(output.IsDefined(CONTACT), world_rank == root);
        KRATOS_CHECK_EQUAL(output.Is(CONTACT), world_rank == root);
        // not set and not involved in communication
        KRATOS_CHECK_EQUAL(output.IsDefined(TO_ERASE), false);
    }
    else {
        KRATOS_CHECK_EQUAL(output, flags);
    }
}


// Flags AndAll ///////////////////////////////////////////////////////////////

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorFlagsAndAllOperations, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_rank = mpi_world_communicator.Rank();
    const int world_size = mpi_world_communicator.Size();
    constexpr int root = 0;

    Kratos::Flags flags;
    //       both true | both false | opposite sets | first true   | first false | second true | second false
    if (world_rank == root) {
        flags = ACTIVE | RIGID.AsFalse()  | STRUCTURE     | MPI_BOUNDARY | PERIODIC.AsFalse();
    }
    else {
        flags = ACTIVE | RIGID.AsFalse()  | STRUCTURE.AsFalse() |                              INLET       | OUTLET.AsFalse();
    }

    // Setting an extra flag, not involved in communication
    if (world_rank == root)
        flags.Set(CONTACT, true);

    Flags output = mpi_world_communicator.AndReduceAll(flags, ACTIVE | RIGID | STRUCTURE | MPI_BOUNDARY | PERIODIC | INLET | OUTLET | ISOLATED);

    if (world_size > 1) {
        // true (defined) & true (defined) = true (defined)
        KRATOS_CHECK_EQUAL(output.IsDefined(ACTIVE), true);
        KRATOS_CHECK_EQUAL(output.Is(ACTIVE), true);
        // false (defined) & false (defined) = false (defined)
        KRATOS_CHECK_EQUAL(output.IsDefined(RIGID), true);
        KRATOS_CHECK_EQUAL(output.Is(RIGID), false);
        // true (defined) & false (defined) = false (defined)
        KRATOS_CHECK_EQUAL(output.IsDefined(STRUCTURE), true);
        KRATOS_CHECK_EQUAL(output.Is(STRUCTURE), false);
        // true (defined) & (undefined) = false (defined) (undefined defaults to false)
        KRATOS_CHECK_EQUAL(output.IsDefined(MPI_BOUNDARY), true);
        KRATOS_CHECK_EQUAL(output.Is(MPI_BOUNDARY), false);
        // false (defined) & (undefined) = false (defined) (undefined defaults to false)
        KRATOS_CHECK_EQUAL(output.IsDefined(PERIODIC), true);
        KRATOS_CHECK_EQUAL(output.Is(PERIODIC), false);
        // (undefined) & true (defined) = false (defined) (undefined defaults to false)
        KRATOS_CHECK_EQUAL(output.IsDefined(INLET), true);
        KRATOS_CHECK_EQUAL(output.Is(INLET), false);
        // (undefied) & false (defined) = false (defined) (undefined defaults to false)
        KRATOS_CHECK_EQUAL(output.IsDefined(OUTLET), true);
        KRATOS_CHECK_EQUAL(output.Is(OUTLET), false);
        // (undefined) & (undefined) = (undefined)
        KRATOS_CHECK_EQUAL(output.IsDefined(ISOLATED), false);
        KRATOS_CHECK_EQUAL(output.Is(ISOLATED), false);

        // set, but not involved in communication
        KRATOS_CHECK_EQUAL(output.IsDefined(CONTACT), world_rank == root);
        KRATOS_CHECK_EQUAL(output.Is(CONTACT), world_rank == root);
        // not set and not involved in communication
        KRATOS_CHECK_EQUAL(output.IsDefined(TO_ERASE), false);
    }
    else {
        KRATOS_CHECK_EQUAL(output, flags);
    }
}

// Flags OrAll ////////////////////////////////////////////////////////////////

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorFlagsOrAllOperations, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_rank = mpi_world_communicator.Rank();
    const int world_size = mpi_world_communicator.Size();
    constexpr int root = 0;

    Kratos::Flags flags;
    //       both true | both false | opposite sets | first true   | first false | second true | second false
    if (world_rank == root) {
        flags = ACTIVE | RIGID.AsFalse()  | STRUCTURE     | MPI_BOUNDARY | PERIODIC.AsFalse();
    }
    else {
        flags = ACTIVE | RIGID.AsFalse()  | STRUCTURE.AsFalse() |                              INLET       | OUTLET.AsFalse();
    }

    // Setting an extra flag, not involved in communication
    if (world_rank == root)
        flags.Set(CONTACT, true);

    Flags output = mpi_world_communicator.OrReduceAll(flags, ACTIVE | RIGID | STRUCTURE | MPI_BOUNDARY | PERIODIC | INLET | OUTLET | ISOLATED);

    if (world_size > 1) {
        // true (defined) | true (defined) = true (defined)
        KRATOS_CHECK_EQUAL(output.IsDefined(ACTIVE), true);
        KRATOS_CHECK_EQUAL(output.Is(ACTIVE), true);
        // false (defined) | false (defined) = false (defined)
        KRATOS_CHECK_EQUAL(output.IsDefined(RIGID), true);
        KRATOS_CHECK_EQUAL(output.Is(RIGID), false);
        // true (defined) | false (defined) = true (defined)
        KRATOS_CHECK_EQUAL(output.IsDefined(STRUCTURE), true);
        KRATOS_CHECK_EQUAL(output.Is(STRUCTURE), true);
        // true (defined) | (undefined) = true (defined) (undefined defaults to false)
        KRATOS_CHECK_EQUAL(output.IsDefined(MPI_BOUNDARY), true);
        KRATOS_CHECK_EQUAL(output.Is(MPI_BOUNDARY), true);
        // false (defined) | (undefined) = false (defined) (undefined defaults to false)
        KRATOS_CHECK_EQUAL(output.IsDefined(PERIODIC), true);
        KRATOS_CHECK_EQUAL(output.Is(PERIODIC), false);
        // (undefined) | true (defined) = true (defined) (undefined defaults to false)
        KRATOS_CHECK_EQUAL(output.IsDefined(INLET), true);
        KRATOS_CHECK_EQUAL(output.Is(INLET), true);
        // (undefied) | false (defined) = false (defined) (undefined defaults to false)
        KRATOS_CHECK_EQUAL(output.IsDefined(OUTLET), true);
        KRATOS_CHECK_EQUAL(output.Is(OUTLET), false);
        // (undefined) | (undefined) = (undefined)
        KRATOS_CHECK_EQUAL(output.IsDefined(ISOLATED), false);
        KRATOS_CHECK_EQUAL(output.Is(ISOLATED), false);

        // set, but not involved in communication
        KRATOS_CHECK_EQUAL(output.IsDefined(CONTACT), world_rank == root);
        KRATOS_CHECK_EQUAL(output.Is(CONTACT), world_rank == root);
        // not set and not involved in communication
        KRATOS_CHECK_EQUAL(output.IsDefined(TO_ERASE), false);
    }
    else {
        KRATOS_CHECK_EQUAL(output, flags);
    }
}


}

}