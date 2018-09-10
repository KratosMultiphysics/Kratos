#include "custom_io/hdf5_properties_io.h"

#include <sstream>
#include "includes/kratos_components.h"
#include "custom_io/hdf5_data_value_container_io.h"
#include "custom_io/hdf5_file.h"

namespace Kratos
{
namespace HDF5
{
namespace Internals
{

void ReadProperties(File& rFile, std::string const& rPrefix, PropertiesContainerType& rProperties)
{
    KRATOS_TRY;

    Vector<int> prop_ids;
    rFile.ReadAttribute(rPrefix + "/Properties", "Ids", prop_ids);

    for (auto pid : prop_ids)
    {
        std::stringstream pstream;
        pstream << rPrefix << "/Properties/(" << pid << ")";
        std::string path = pstream.str();

        PropertiesType::ContainerType& r_data = rProperties[pid].Data();
        Internals::ReadDataValueContainer(rFile, path, r_data);
    }

    KRATOS_CATCH("");
}

void WriteProperties(File& rFile, std::string const& rPrefix, PropertiesType const& rProperties)
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

void WriteProperties(File& rFile, std::string const& rPrefix, PropertiesContainerType const& rProperties)
{
    KRATOS_TRY;

    Vector<int> prop_ids(rProperties.size());
    unsigned i = 0;
    for (const PropertiesType& r_properties : rProperties)
    {
        prop_ids[i++] = r_properties.Id();
        WriteProperties(rFile, rPrefix, r_properties);
    }
    rFile.AddPath(rPrefix + "/Properties");
    rFile.WriteAttribute(rPrefix + "/Properties", "Ids", prop_ids);

    KRATOS_CATCH("");
}

} // namespace Internals.
} // namespace HDF5.
} // namespace Kratos.
