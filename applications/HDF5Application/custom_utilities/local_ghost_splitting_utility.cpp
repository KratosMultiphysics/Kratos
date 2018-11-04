#include "custom_utilities/local_ghost_splitting_utility.h"

#include "utilities/openmp_utils.h"
#ifdef KRATOS_USING_MPI
#include "mpi.h"
#endif

namespace Kratos
{

void SplitNodesIntoLocalAndGhost(ModelPart::NodesContainerType const& rNodes,
                                 std::vector<ModelPart::NodeType*>& rLocalNodes,
                                 std::vector<ModelPart::NodeType*>& rGhostNodes)
{
    KRATOS_TRY;

    rLocalNodes.clear();
    rGhostNodes.clear();
    if (rNodes.size() == 0)
        return;

#ifdef KRATOS_USING_MPI
    int mpi_is_initialized;
    KRATOS_ERROR_IF(MPI_Initialized(&mpi_is_initialized) != MPI_SUCCESS) << "MPI_Initialized failed." << std::endl;
    int my_pid{0}, num_proc{1};
    if (mpi_is_initialized)
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &my_pid);
        MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    }
    const bool is_partitioned = (num_proc > 1);
#else /* KRATOS_USING_MPI */
    const int my_pid = 0;
    const bool is_partitioned = false;
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
#pragma omp parallel for
        for (int i = 0; i < static_cast<int>(rNodes.size()); ++i)
        {
            auto it = rNodes.begin() + i;
            rLocalNodes[i] = &(*it);
            ++it;
        }
    }

    KRATOS_CATCH("");
}

void GetLocalNodes(ModelPart::NodesContainerType const& rNodes,
                   std::vector<ModelPart::NodeType*>& rLocalNodes)
{
    KRATOS_TRY;

    std::vector<ModelPart::NodeType*> ghost_nodes;
    SplitNodesIntoLocalAndGhost(rNodes, rLocalNodes, ghost_nodes);

    KRATOS_CATCH("");
}

} // namespace Kratos.