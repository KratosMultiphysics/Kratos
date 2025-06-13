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

    IndexType start_index, block_size;
    std::tie(start_index, block_size) = StartIndexAndBlockSize(rFile, rPrefix + "/Properties/");

    Matrix<char> properties_data;
    rFile.ReadDataSet(rPrefix + "/Properties/ListOfProperties", properties_data, start_index, block_size);

    std::string current_str_data;
    current_str_data.resize(properties_data.size2());
    for (IndexType i = 0; i < properties_data.size1(); ++i) {
        const Vector<char> current_vec_data = row(properties_data, i);
        std::copy(current_vec_data.begin(), current_vec_data.end(), current_str_data.begin());

        // now look for the delimiter
        const auto pos = current_str_data.find(Delimiter);
        KRATOS_ERROR_IF(pos == std::string::npos) << "The delimiter is not found for properties:" << current_str_data;

        std::string serialized_data;
        serialized_data.resize(pos);
        std::copy(current_str_data.begin(), current_str_data.begin() + pos, serialized_data.begin());
        StreamSerializer serializer(serialized_data);
        auto p_props = Kratos::make_shared<Properties>();
        serializer.load("properties", *p_props);
        rProperties.insert(rProperties.end(), p_props);
    }

    KRATOS_CATCH("");
}

void WriteProperties(
    File& rFile,
    const std::string& rPrefix,
    const PropertiesContainerType & rProperties)
{
    KRATOS_TRY;

    std::vector<std::string> serialized_properties;
    IndexType max_length = 0;
    for (const PropertiesType& r_properties : rProperties) {
        StreamSerializer serializer;
        serializer.save("properties", r_properties);
        const auto& serialized_string = serializer.GetStringRepresentation();
        serialized_properties.push_back(serialized_string);
        max_length = std::max(max_length, serialized_string.size());
    }

    // now get the max string length by MPI communication
    // added an additional char to the length to mark the termination char
    max_length = rFile.GetDataCommunicator().MaxAll(max_length) + Delimiter.size();

    // now fill in the serialized data
    Matrix<char> properties_data;
    properties_data.resize(serialized_properties.size(), max_length);
    for (IndexType i = 0; i < serialized_properties.size(); ++i) {
        std::copy(serialized_properties[i].begin(), serialized_properties[i].end(), row(properties_data, i).begin());
        // mark the end of the char array
        std::copy(Delimiter.begin(), Delimiter.end(),  row(properties_data, i).begin() + serialized_properties[i].size());
    }

    rFile.AddPath(rPrefix + "/Properties");
    WriteInfo info;
    rFile.WriteDataSet(rPrefix + "/Properties/ListOfProperties", properties_data, info);
    WritePartitionTable(rFile, rPrefix + "/Properties/", info);

    KRATOS_CATCH("");
}

} // namespace Internals.
} // namespace HDF5.
} // namespace Kratos.
