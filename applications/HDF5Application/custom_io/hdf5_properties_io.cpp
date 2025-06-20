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
//

// System includes
#include <sstream>

// Project includes
#include "includes/kratos_components.h"
#include "includes/serializer.h"

// Application includes
#include "custom_io/hdf5_data_value_container_io.h"
#include "custom_io/hdf5_file.h"
#include "custom_utilities/hdf5_data_set_partition_utility.h"

// Include base h
#include "custom_io/hdf5_properties_io.h"

namespace Kratos
{
namespace HDF5
{
namespace Internals
{

void ReadProperties(
    File& rFile,
    const std::string& rPrefix,
    PropertiesContainerType& rProperties)
{
    KRATOS_TRY;

    KRATOS_ERROR_IF_NOT(rProperties.empty()) << "The properties container is not empty.";

    const auto [start_index, block_size] = StartIndexAndBlockSize(rFile, rPrefix + "/Properties/");

    std::string serialized_data;
    rFile.ReadDataSet(rPrefix + "/Properties/ListOfProperties", serialized_data, start_index, block_size);

    StreamSerializer serializer(serialized_data);
    serializer.load("properties", rProperties);

    KRATOS_CATCH("");
}

void WriteProperties(
    File& rFile,
    const std::string& rPrefix,
    const PropertiesContainerType & rProperties)
{
    KRATOS_TRY;

    rFile.AddPath(rPrefix + "/Properties");

    StreamSerializer serializer;
    serializer.save("properties", rProperties);

    WriteInfo info;
    rFile.WriteDataSet(rPrefix + "/Properties/ListOfProperties", serializer.GetStringRepresentation(), info);
    WritePartitionTable(rFile, rPrefix + "/Properties/", info);

    KRATOS_CATCH("");
}

} // namespace Internals.
} // namespace HDF5.
} // namespace Kratos.
