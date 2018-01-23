#include "custom_utilities/hdf5_points_data.h"

#include "utilities/openmp_utils.h"

namespace Kratos
{
namespace HDF5
{
namespace Internals
{
void PointsData::ReadData(File& rFile, std::string Path, unsigned StartIndex, unsigned BlockSize)
{
    KRATOS_TRY;
    Clear();
    rFile.ReadDataSet(Path + "/Ids", mIds, StartIndex, BlockSize);
    rFile.ReadDataSet(Path + "/Coordinates", mCoords, StartIndex, BlockSize);
    KRATOS_CATCH("");
}

void PointsData::WriteData(File& rFile, std::string Path)
{
    KRATOS_TRY;
    rFile.WriteDataSet(Path + "/Ids", mIds);
    rFile.WriteDataSet(Path + "/Coordinates", mCoords);
    KRATOS_CATCH("");
}

void PointsData::CreateNodes(NodesContainerType& rNodes)
{
    KRATOS_TRY;
    const unsigned num_new_nodes = mIds.size();
    rNodes.reserve(rNodes.size() + num_new_nodes);
    for (unsigned i = 0; i < num_new_nodes; ++i)
    {
        const array_1d<double, 3>& r_coord = mCoords[i];
        NodeType::Pointer p_node = boost::make_shared<NodeType>(
            mIds[i], r_coord[0], r_coord[1], r_coord[2]);
        rNodes.push_back(p_node);
    }
    KRATOS_CATCH("");
}

void PointsData::SetData(NodesContainerType const& rNodes)
{
    KRATOS_TRY;

    const unsigned num_nodes = rNodes.size();
    mIds.resize(num_nodes, false);
    mCoords.resize(num_nodes, false);

    const int num_threads = OpenMPUtils::GetNumThreads();
    OpenMPUtils::PartitionVector partition;
    OpenMPUtils::DivideInPartitions(num_nodes, num_threads, partition);
#pragma omp parallel
    {
        const int thread_id = OpenMPUtils::ThisThread();
        NodesContainerType::const_iterator it = rNodes.begin() + partition[thread_id];
        for (auto i = partition[thread_id]; i < partition[thread_id + 1]; ++i)
        {
            const auto& r_node = *it;
            mIds[i] = r_node.Id();
            mCoords[i] = r_node.Coordinates();
            ++it;
        }
    }
    
    KRATOS_CATCH("");
}

void PointsData::Clear()
{
    mIds.clear();
    mCoords.clear();
}

} // namespace Internals.
} // namespace HDF5.
} // namespace Kratos.
