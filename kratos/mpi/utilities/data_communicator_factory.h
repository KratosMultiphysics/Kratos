//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Jordi Cotela
//

#ifndef KRATOS_DATA_COMMUNICATOR_FACTORY_INCLUDED
#define KRATOS_DATA_COMMUNICATOR_FACTORY_INCLUDED

#include <string>

#include "includes/data_communicator.h"

namespace Kratos
{

/// Common tools to define new MPI DataCommunicators.
/** Note that new DataCommunicators created by the functions in the
 *  DataCommunicatorFactory namespace are automatically registered
 *  in the ParallelEnvironment, that takes care of managing their
 *  life cycle. In this way, we ensure that the underlying MPI_Comm
 *  object is always properly freed before calling MPI_Finalize.
 *
 *  Note that some of the functions can create a communicator that is
 *  undefined in some processes, because that process does not participate
 *  in the communication. In MPI terms, the returned DataCommunicator is
 *  a wrapper for MPI_COMM_NULL. Since MPI calls are undefined for MPI_COMM_NULL,
 *  communicators created using those functions should be used carefully
 *  checking DataCommunicator::IsDefinedOnThisRank or
 *  DataCommunicator::IsNullOnThisRank as appropriate. All MPI calls,
 *  including Rank and Size, will result in errors for null communicators.
 */
namespace DataCommunicatorFactory
{

/// Create a new MPIDataCommunicator as a duplicate of an existing one.
/** This function is a wrapper for MPI_Comm_dup. Note that, if a serial
 *  DataCommunicator is passed as an argument, the returned DataCommunicator
 *  will be a wrapper for a duplicate of MPI_COMM_SELF.
 *  @param rOriginalCommunicator The DataCommunicator to be copied.
 *  @param rNewCommunicatorName The name to register the new DataCommunicator with.
 *  @return A reference to the duplicate DataCommunicator.
 */
const DataCommunicator& DuplicateAndRegister(
    const DataCommunicator& rOriginalCommunicator,
    const std::string& rNewCommunicatorName);

/// Create a new MPIDataCommunicator by splitting an existing one.
/** This function is a wrapper for MPI_Comm_split. The resulting communicator
 *  divides all processes in the original one into disjoint groups, one for
 *  each distinct value of Color.
 *
 *  A process may provide MPI_UNDEFINED as color. In this case, the returned
 *  communicator will be MPI_COMM_NULL for that process (the process will not
 *  participate in communication).
 *
 *  The Key argument is used to determine the relative ordering of processes
 *  within the same color group: ranks will be assigned in order of increasing
 *  Key. If ranks with the same Key and Color are found, they will be sorted
 *  according to increasing rank in the original communicator. In practice,
 *  assigning Key=0 to all processes will assign ranks respecting the rank
 *  ordering in the parent communicator.
 *
 *  @param rOriginalCommunicator The DataCommunicator to be split.
 *  @param Color The tag identifying which group the process will belong to.
 *  @param Key A value indicating the relative ordering (rank) of the process
 *  within the new communicator.
 *  @param rNewCommunicatorName The name to register the new DataCommunicator with.
 *  @return A reference to the split DataCommunicator.
 */
const DataCommunicator& SplitAndRegister(
    const DataCommunicator& rOriginalCommunicator,
    int Color,
    int Key,
    const std::string& rNewCommunicatorName);

/// Create a new MPIDataCommunicator connecting the provided ranks.
/** This function defines and registers in ParallelEnvironment a new
 *  communicator involving the selected ranks of rOriginalCommunicator.
 *
 *  Note that this returns a wrapper for MPI_COMM_NULL in all other ranks.
 *
 *  @param rOriginalCommunicator Reference DataCommunicator.
 *  @param rRanks list of ranks of rOriginalCommunicator to be used in
 *  the new communicator.
 *  @param rNewCommunicatorName The name to register the new DataCommunicator with.
 *  @return A reference to the new DataCommunicator.
 */
const DataCommunicator& CreateFromRanksAndRegister(
    const DataCommunicator& rOriginalCommunicator,
    const std::vector<int>& rRanks,
    const std::string& rNewCommunicatorName);

/// Create a new MPIDataCommunicator as the union of the given ones.
/** This function creates and registers in ParallelEnvironment a new
 *  MPIDataCommunicator involving all ranks that participate on either
 *  of the arguments. A parent DataCommunicator encompassing all ranks
 *  in either of the two arguments (and possibly more) is required.
 *  The operation is collective in all processes in rParentDataCommunicator.
 *  For processes involved in rParentDataCommunicator but not on either of
 *  the first two arguments, the new communicator will be a wrapper to
 *  MPI_COMM_NULL.
 *  @param rFirstDataCommunicator first argument of the union.
 *  @param rSecondDataCommunicator second argument of the union.
 *  @param rParentDataCommunicator auxiliary DataCommunicator involving
 *  all ranks on either of the first two arguments (a wrapper to
 *  MPI_COMM_WORLD is always a valid argument here).
 *  @param rNewCommunicatorName The name to register the new DataCommunicator with.
 *  @return A reference to the new DataCommunicator.
 */
const DataCommunicator& CreateUnionAndRegister(
    const DataCommunicator& rFirstDataCommunicator,
    const DataCommunicator& rSecondDataCommunicator,
    const DataCommunicator& rParentDataCommunicator,
    const std::string& rNewCommunicatorName);

/// Create a new MPIDataCommunicator as the intersection of the given ones.
/** This function creates and registers in ParallelEnvironment a new
 *  MPIDataCommunicator involving all ranks that participate on both
 *  of the arguments. A parent DataCommunicator encompassing all ranks
 *  in either of the two arguments (and possibly more) is required.
 *  The operation is collective in all processes in rParentDataCommunicator.
 *  For processes involved in rParentDataCommunicator but not on both of
 *  the two arguments, the new communicator will be a wrapper to MPI_COMM_NULL.
 *
 *  @param rFirstDataCommunicator first argument of the intersection.
 *  @param rSecondDataCommunicator second argument of the intersection.
 *  @param rParentDataCommunicator auxiliary DataCommunicator involving
 *  all ranks on either of the first two arguments (a wrapper to
 *  MPI_COMM_WORLD is always a valid argument here).
 *  @param rNewCommunicatorName The name to register the new DataCommunicator with.
 *  @return A reference to the new DataCommunicator.
 */
const DataCommunicator& CreateIntersectionAndRegister(
    const DataCommunicator& rFirstDataCommunicator,
    const DataCommunicator& rSecondDataCommunicator,
    const DataCommunicator& rParentDataCommunicator,
    const std::string& rNewCommunicatorName);

}

}

#endif // KRATOS_DATA_COMMUNICATOR_FACTORY_INCLUDED