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

#pragma once

// System includes
#include <string>

// External includes

// Project includes
#include "includes/define.h"

// Application includes
#include "custom_io/hdf5_file.h"
#include "hdf5_application_define.h"
#include "custom_utilities/hdf5_data_set_partition_utility.h"

namespace Kratos
{
namespace HDF5
{
namespace Internals
{
///@addtogroup HDF5Application
///@{
///@name Kratos Classes
///@{

/// A class representing points in a mesh.
/**
 * @tparam TContainerDataIO A IO class which have the @p ContainerType defined in public scope
 *                          with @p SetValue, @p GetValue methods implemented.
*/
template<class TContainerDataIO>
class PointsData
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(PointsData);

    ///@}
    ///@name Life Cycle
    ///@{

    PointsData(
        const std::string& rPrefix,
        File::Pointer pFile) 
        : mpFile(pFile),
          mPrefix(rPrefix) 
    {

    }

    ///@}
    ///@name Operations
    ///@{

    Parameters Read(
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

    void Write(
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

    ///@}

protected:
    ///@name Protected member Variables
    ///@{

    File::Pointer mpFile;

    std::string mPrefix;

    ///@}
}; // class PointsData

///@} // Kratos Classes
///@} addtogroup
} // namespace Internals.
} // namespace HDF5.
} // namespace Kratos.