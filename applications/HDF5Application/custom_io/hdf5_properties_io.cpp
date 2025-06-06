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

// Application includes
#include "custom_io/hdf5_data_value_container_io.h"
#include "custom_io/hdf5_file.h"

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

    Vector<int> prop_ids;
    if (rFile.HasAttribute(rPrefix + "/Properties", "Ids")) {
        rFile.ReadAttribute(rPrefix + "/Properties", "Ids", prop_ids);
    }

    for (auto pid : prop_ids) {
        std::stringstream pstream;
        pstream << rPrefix << "/Properties/(" << pid << ")";
        std::string path = pstream.str();

        PropertiesType::ContainerType& r_data = rProperties[pid].Data();
        Internals::ReadDataValueContainer(rFile, path, r_data);
    }

    KRATOS_CATCH("");
}

void WriteProperties(
    File& rFile,
    const std::string& rPrefix,
    const PropertiesType& rProperties)
{
    KRATOS_TRY;

    std::stringstream pstream;
    pstream << rPrefix << "/Properties/(" << rProperties.Id() << ")";
    std::string path = pstream.str();
    rFile.AddPath(path);

    const PropertiesType::ContainerType& r_data = rProperties.Data();
    Internals::WriteDataValueContainer(rFile, path, r_data);

    KRATOS_CATCH("");
}

void WriteProperties(
    File& rFile,
    const std::string& rPrefix,
    const PropertiesContainerType & rProperties)
{
    KRATOS_TRY;

    Vector<int> prop_ids(rProperties.size());
    unsigned i = 0;
    for (const PropertiesType& r_properties : rProperties) {
        prop_ids[i++] = r_properties.Id();
        WriteProperties(rFile, rPrefix, r_properties);
    }
    rFile.AddPath(rPrefix + "/Properties");
    if (rProperties.size() > 0) { // H5Awrite fails for empty container (sub model parts)
        rFile.WriteAttribute(rPrefix + "/Properties", "Ids", prop_ids);
    }

    KRATOS_CATCH("");
}

} // namespace Internals.
} // namespace HDF5.
} // namespace Kratos.
