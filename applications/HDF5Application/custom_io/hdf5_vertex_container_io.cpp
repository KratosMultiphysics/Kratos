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

VertexContainerCoordinateIO::VertexContainerCoordinateIO(
    Parameters Settings,
    File::Pointer pFile)
    : mpFile(pFile)
{
    KRATOS_TRY

    Parameters default_params(R"(
        {
            "prefix"           : "",
            "write_vertex_ids" : false
        })");

    Settings.AddMissingParameters(default_params);
    mWriteIDs = Settings["write_vertex_ids"].GetBool();
    mPathPrefix = Settings["prefix"].GetString();

    KRATOS_CATCH("");
}

void VertexContainerCoordinateIO::Write(
    const Detail::VertexContainerType& rVertices,
    Parameters Attributes)
{
    if (mWriteIDs) {
        WriteWithIDs(rVertices, Attributes);
    } else {
        WriteWithoutIDs(rVertices, Attributes);
    }
}


void VertexContainerCoordinateIO::WriteWithIDs(
    const Detail::VertexContainerType& rVertices,
    Parameters Attributes)
{
    KRATOS_TRY

    const std::string path = mPathPrefix + "/POSITION";

    // Collect vertex coordinates into a buffer
    Matrix<double> buffer;
    buffer.resize(rVertices.size(), 4,  false);

    IndexPartition<IndexType>(rVertices.size()).for_each([&rVertices, &buffer](const auto Index) {
        const auto& rVertex = *(rVertices.begin() + Index);

        for (std::size_t i_component=0; i_component < 3; ++i_component) {
            buffer(Index, i_component) = rVertex[i_component];
        }

        buffer(Index, 3) = rVertex.GetID();
    });

    // Write buffer to the requested path
    WriteInfo write_info;
    mpFile->WriteDataSet(path, buffer, write_info);
    mpFile->WriteAttribute(path, Attributes);

    KRATOS_CATCH("");
}

void VertexContainerCoordinateIO::WriteWithoutIDs(
    const Detail::VertexContainerType& rVertices,
    Parameters Attributes)
{
    KRATOS_TRY

    const std::string path = mPathPrefix + "/POSITION";

    // Collect vertex coordinates into a buffer
    Matrix<double> buffer;
    buffer.resize(rVertices.size(), 3,  false);

    IndexPartition<IndexType>(rVertices.size()).for_each([&rVertices, &buffer](const auto Index) {
        const auto& rVertex = *(rVertices.begin() + Index);

        for (std::size_t i_component=0; i_component < 3; ++i_component) {
            buffer(Index, i_component) = rVertex[i_component];
        }
    });

    // Write buffer to the requested path
    WriteInfo write_info;
    mpFile->WriteDataSet(path, buffer, write_info);
    mpFile->WriteAttribute(path, Attributes);

    KRATOS_CATCH("");
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