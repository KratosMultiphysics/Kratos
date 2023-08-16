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
}

void VertexContainerCoordinateIO::Write(
    const Detail::VertexContainerType& rVertices,
    Parameters Attributes)
{
    BaseType::Write(rVertices, Internals::VerticesIO{});
    mpFile->WriteAttribute(mPrefix, Attributes);
}

VertexContainerVariableIO::VertexContainerVariableIO(
    Parameters Settings,
    File::Pointer pFile)
    : BaseType(Settings, pFile)
{
}

void VertexContainerVariableIO::Write(
    const Detail::VertexContainerType& rVertices,
    Parameters Atttributes)
{
    BaseType::Write(rVertices, Internals::VertexValueIO{}, Atttributes);
}

} // namespace HDF5
} // namespace Kratos