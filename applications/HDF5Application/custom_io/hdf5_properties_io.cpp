#include "custom_io/hdf5_properties_io.h"

#include <sstream>
#include "includes/kratos_components.h"
#include "custom_io/hdf5_data_value_container_io.h"

namespace Kratos
{
namespace HDF5
{
namespace Internals
{

PropertiesIO::PropertiesIO(std::string Prefix, File::Pointer pFile)
    : mPrefix(Prefix), mpFile(pFile)
{
}

void PropertiesIO::ReadProperties(PropertiesContainerType& rProperties)
{
    KRATOS_TRY;

    Vector<int> prop_ids;
    mpFile->ReadAttribute(mPrefix + "/Properties", "Ids", prop_ids);

    for (unsigned i = 0; i < prop_ids.size(); ++i)
    {
        std::stringstream pstream;
        pstream << mPrefix << "/Properties/(" << prop_ids[i] << ")";
        std::string path = pstream.str();

        DataValueContainerIO data_io(path, mpFile);
        PropertiesType::ContainerType& r_data = rProperties[prop_ids[i]].Data();
        data_io.ReadDataValueContainer(r_data);
    }

    KRATOS_CATCH("");
}

void PropertiesIO::WriteProperties(PropertiesType const& rProperties)
{
    KRATOS_TRY;

    std::stringstream pstream;
    pstream << mPrefix << "/Properties/(" << rProperties.Id() << ")";
    std::string path = pstream.str();
    mpFile->AddPath(path);

    DataValueContainerIO data_io(path, mpFile);
    const PropertiesType::ContainerType& r_data = rProperties.Data();
    data_io.WriteDataValueContainer(r_data);

    KRATOS_CATCH("");
}

void PropertiesIO::WriteProperties(PropertiesContainerType const& rProperties)
{
    KRATOS_TRY;

    Vector<int> prop_ids(rProperties.size());
    unsigned i = 0;
    for (const PropertiesType& r_properties : rProperties)
    {
        prop_ids[i++] = r_properties.Id();
        WriteProperties(r_properties);
    }
    mpFile->AddPath(mPrefix + "/Properties");
    mpFile->WriteAttribute(mPrefix + "/Properties", "Ids", prop_ids);

    KRATOS_CATCH("");
}

} // namespace Internals.
} // namespace HDF5.
} // namespace Kratos.
