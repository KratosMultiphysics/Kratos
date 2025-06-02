//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: HDF5Application/license.txt
//
//  Main author:     Máté Kelemen
//

// Project includes
#include "utilities/parallel_utilities.h"

// Application includes

// Include base h
#include "hdf5_vertex_container_io.h"

namespace Kratos
{
namespace HDF5
{

namespace Internals
{
std::string AddMissingAndGetPrefix(Parameters Settings)
{
    Parameters default_params(R"(
        {
            "prefix": ""
        })");

    Settings.AddMissingParameters(default_params);
    return Settings["prefix"].GetString();
}
}

VertexContainerCoordinateIO::VertexContainerCoordinateIO(
    Parameters Settings,
    File::Pointer pFile)
    : BaseType(Internals::AddMissingAndGetPrefix(Settings), pFile)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(Settings["prefix"].GetString().back() == '/')
        << "The prefix for vertex coordinates assumed to be a group hence no need to have an ending \"\\\" [ prefix = \""
        << Settings["prefix"].GetString() << "\" ].\n";

    KRATOS_CATCH("");
}

void VertexContainerCoordinateIO::Write(
    const Detail::VertexContainerType& rVertices,
    Parameters Attributes)
{
    BaseType::Write(rVertices, Internals::VertexIO{}, Attributes);
}

template<class TVertexDataIOType>
VertexContainerVariableIO<TVertexDataIOType>::VertexContainerVariableIO(
    Parameters Settings,
    File::Pointer pFile)
    : BaseType(Settings, pFile)
{
}

template<class TVertexDataIOType>
void VertexContainerVariableIO<TVertexDataIOType>::Write(
    const Detail::VertexContainerType& rVertices,
    const TVertexDataIOType& rVertexDataIO,
    const Parameters Atttributes)
{
    BaseType::Write(rVertices, rVertexDataIO, Atttributes);
}

// template instantiations
template class KRATOS_API(HDF5_APPLICATION) VertexContainerVariableIO<Internals::VertexHistoricalValueIO>;
template class KRATOS_API(HDF5_APPLICATION) VertexContainerVariableIO<Internals::VertexNonHistoricalValueIO>;

} // namespace HDF5
} // namespace Kratos