#include "custom_io/hdf5_properties_io.h"

namespace Kratos
{
namespace HDF5
{
PropertiesIO::PropertiesIO(std::string Prefix, File::Pointer pFile)
    : mPrefix(Prefix), mpFile(pFile)
{
}

void PropertiesIO::ReadProperties(PropertiesContainerType& rThisProperties)
{
}

void PropertiesIO::WriteProperties(Properties const& rThisProperties)
{
}

void PropertiesIO::WriteProperties(PropertiesContainerType const& rThisProperties)
{
}

} // namespace HDF5.
} // namespace Kratos.
