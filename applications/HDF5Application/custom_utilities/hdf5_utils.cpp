#include "custom_utilities/hdf5_utils.h"

#include "utilities/openmp_utils.h"
#ifdef KRATOS_USING_MPI
#include "mpi.h"
#endif

namespace Kratos
{
namespace HDF5
{
namespace Detail
{
void DivideNodes(NodesContainerType const& rNodes,
                 std::vector<NodeType*>& rLocalNodes,
                 std::vector<NodeType*>& rGhostNodes)
{
    KRATOS_TRY;

    bool is_partitioned = false;
    // Check if we are simulating a partitioned problem.
#ifdef KRATOS_USING_MPI
    int mpi_is_initialized, num_proc, my_pid, ierr;
    ierr = MPI_Initialized(&mpi_is_initialized);
    KRATOS_ERROR_IF(ierr != MPI_SUCCESS) << "MPI_Initialized failed." << std::endl;
    if (mpi_is_initialized)
    {
        ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
        KRATOS_ERROR_IF(ierr != MPI_SUCCESS) << "MPI_Comm_size failed." << std::endl;
        is_partitioned = (num_proc > 1);
    }
    else
        num_proc = 1;

    if (is_partitioned)
    {
        ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_pid);
        KRATOS_ERROR_IF(ierr != MPI_SUCCESS) << "MPI_Comm_rank failed." << std::endl;
    }
    else
        my_pid = 0;
#endif

    if (is_partitioned)
    {
        rLocalNodes.reserve(rNodes.size());
        rGhostNodes.reserve(0.1*rNodes.size());
        for (auto it = rNodes.begin(); it != rNodes.end(); ++it)
        {
            if (it->FastGetSolutionStepValue(PARTITION_INDEX) == my_pid)
                rLocalNodes.push_back(&(*it));
            else
                rGhostNodes.push_back(&(*it));
        }
    }
    else
    {
        rLocalNodes.resize(rNodes.size());
        rGhostNodes.resize(0);
        const int num_threads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector partition;
        OpenMPUtils::DivideInPartitions(rNodes.size(), num_threads, partition);
#pragma omp parallel
        {
            const int thread_id = OpenMPUtils::ThisThread();
            NodesContainerType::const_iterator it = rNodes.begin() + partition[thread_id];
            for (auto i = partition[thread_id]; i < partition[thread_id + 1]; ++i)
            {
                rLocalNodes[i] = &(*it);
                ++it;
            }
        }
    }

    KRATOS_CATCH("");
}

void GetLocalNodes(NodesContainerType const& rNodes, std::vector<NodeType*>& rLocalNodes)
{
    KRATOS_TRY;

    std::vector<NodeType*> ghost_nodes;
    DivideNodes(rNodes, rLocalNodes, ghost_nodes);

    KRATOS_CATCH("");
}

} // namespace Detail.
} // namespace HDF5.
} // namespace Kratos.
