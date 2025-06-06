//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  license: HDF5Application/license.txt
//
//  Main author:    Michael Andre, https://github.com/msandre
//                  Suneth Warnakulasuriya
//

// System includes

// Project includes

// Application includes
#include "custom_io/hdf5_file.h"
#include "custom_utilities/container_io_utils.h"
#include "custom_utilities/hdf5_data_set_partition_utility.h"

// Include base h
#include "custom_io/hdf5_points_data.h"

namespace Kratos
{
namespace HDF5
{
namespace Internals
{

template<class TContainerDataIO>
PointsData<TContainerDataIO>::PointsData(
    const std::string& rPrefix,
    File::Pointer pFile)
    : mpFile(pFile),
      mPrefix(rPrefix)
{
}

template<class TContainerDataIO>
Parameters PointsData<TContainerDataIO>::Read(
    typename TContainerDataIO::ContainerType& rContainer,
    const TContainerDataIO& rContainerDataIO)
{
    KRATOS_TRY;

    IndexType start_index, block_size;
    std::tie(start_index, block_size) = StartIndexAndBlockSize(*mpFile, mPrefix);

    Vector<int> ids;
    Vector<array_1d<double, 3>> coords;

    mpFile->ReadDataSet(mPrefix + "/Ids", ids, start_index, block_size);
    mpFile->ReadDataSet(mPrefix + "/Coordinates", coords, start_index, block_size);
    auto attributes = mpFile->ReadAttribute(mPrefix);

    const unsigned num_new_nodes = ids.size();
    rContainer.reserve(rContainer.size() + num_new_nodes);

    for (unsigned i = 0; i < num_new_nodes; ++i) {
        rContainerDataIO.AddPoint(rContainer, ids[i], coords[i]);
    }

    return attributes;

    KRATOS_CATCH("");
}

template<class TContainerDataIO>
void PointsData<TContainerDataIO>::Write(
    const typename TContainerDataIO::ContainerType& rContainer,
    const TContainerDataIO& rContainerDataIO,
    const Parameters Attributes)
{
    KRATOS_TRY;

    const unsigned num_nodes = rContainer.size();

    Vector<int> ids(num_nodes);
    Vector<array_1d<double, 3>> coords(num_nodes);

    IndexPartition<IndexType>(num_nodes).for_each([&rContainer, &rContainerDataIO, &ids, &coords](const auto Index) {
        const auto& r_point = *(rContainer.begin() + Index);
        rContainerDataIO.GetData(ids[Index], coords[Index], r_point);
    });

    WriteInfo info;
    mpFile->WriteDataSet(mPrefix + "/Ids", ids, info);
    mpFile->WriteDataSet(mPrefix + "/Coordinates", coords, info);
    mpFile->WriteAttribute(mPrefix, Attributes);

    WritePartitionTable(*mpFile, mPrefix, info);

    KRATOS_CATCH("");
}

// template instantiations
template class KRATOS_API(HDF5_APPLICATION) PointsData<Internals::NodesIO>;
template class KRATOS_API(HDF5_APPLICATION) PointsData<Internals::VertexIO>;

} // namespace Internals.
} // namespace HDF5.
} // namespace Kratos.
