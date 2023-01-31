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

// Internal includes
#include "hdf5_vertex_container_io.h"


namespace Kratos
{
namespace HDF5
{


VertexContainerIO::VertexContainerIO(Parameters parameters, File::Pointer pFile)
        : VertexContainerIO::BaseType(
            VertexContainerIO::FormatParameters(parameters),
            pFile,
            VertexContainerIO::GetPath(parameters))
{
    KRATOS_TRY

    parameters.AddMissingParameters(VertexContainerIO::GetDefaultParameters());

    KRATOS_CATCH("");
}


Parameters VertexContainerIO::GetDefaultParameters()
{
    return Parameters(R"({
        "prefix"            : "/",
        "list_of_variables" : []
    })");
}


Parameters VertexContainerIO::FormatParameters(Parameters parameters)
{
    KRATOS_TRY

    Parameters output = parameters.Clone();
    output.AddMissingParameters(VertexContainerIO::GetDefaultParameters());
    output["prefix"].SetString("");
    return output;

    KRATOS_CATCH("");
}


std::string VertexContainerIO::GetPath(Parameters parameters)
{
    KRATOS_TRY

    parameters.AddMissingParameters(VertexContainerIO::GetDefaultParameters());
    return parameters["prefix"].GetString();

    KRATOS_CATCH("");
}


VertexContainerCoordinateIO::VertexContainerCoordinateIO(Parameters parameters,
                                                         File::Pointer pFile)
    : VertexContainerIO(VertexContainerCoordinateIO::FormatParameters(parameters), pFile)
{
    KRATOS_TRY

    parameters.AddMissingParameters(VertexContainerCoordinateIO::GetDefaultParameters());
    mWriteIDs = parameters["write_vertex_ids"].GetBool();

    KRATOS_CATCH("");
}


void VertexContainerCoordinateIO::Write(const Detail::VertexContainerType& rVertices)
{
    if (mWriteIDs)
        WriteWithIDs(rVertices);
    else
        WriteWithoutIDs(rVertices);
}


void VertexContainerCoordinateIO::WriteWithIDs(const Detail::VertexContainerType& rVertices)
{
    KRATOS_TRY

    const std::string path = mComponentPath + "/POSITION";

    // Collect vertex coordinates into a buffer
    Matrix<double> buffer;
    buffer.resize(rVertices.size(),
                  4,
                  false);

    for (std::size_t i_vertex=0; i_vertex<rVertices.size(); ++i_vertex) {
        const auto& rVertex = *(rVertices.begin() + i_vertex);

        for (std::size_t i_component=0; i_component<3; ++i_component) {
            buffer(i_vertex, i_component) = rVertex[i_component];
        }

        buffer(i_vertex, 3) = rVertex.GetID();
    }

    // Write buffer to the requested path
    WriteInfo write_info;
    mpFile->WriteDataSet(path,
                         buffer,
                         write_info);

    KRATOS_CATCH("");
}


void VertexContainerCoordinateIO::WriteWithoutIDs(const Detail::VertexContainerType& rVertices)
{
    KRATOS_TRY

    const std::string path = mComponentPath + "/POSITION";

    // Collect vertex coordinates into a buffer
    Matrix<double> buffer;
    buffer.resize(rVertices.size(),
                  3,
                  false);

    for (std::size_t i_vertex=0; i_vertex<rVertices.size(); ++i_vertex) {
        const auto& rVertex = *(rVertices.begin() + i_vertex);

        for (std::size_t i_component=0; i_component<3; ++i_component) {
            buffer(i_vertex, i_component) = rVertex[i_component];
        }
    }

    // Write buffer to the requested path
    WriteInfo write_info;
    mpFile->WriteDataSet(path,
                         buffer,
                         write_info);

    KRATOS_CATCH("");
}


Parameters VertexContainerCoordinateIO::GetDefaultParameters()
{
    return Parameters(R"({
        "prefix"            : "/POSITION",
        "write_vertex_ids"  : false
    })");
}


Parameters VertexContainerCoordinateIO::FormatParameters(Parameters parameters)
{
    KRATOS_TRY

    parameters.AddMissingParameters(VertexContainerCoordinateIO::GetDefaultParameters());

    Parameters output = parameters.Clone();
    output.RemoveValue("write_vertex_ids");

    return output;

    KRATOS_CATCH("");
}


void VertexContainerVariableIO::Write(const Detail::VertexContainerType& rVertices)
{
    KRATOS_TRY

    this->WriteContainerComponents(rVertices);

    KRATOS_CATCH("");
}


Parameters VertexContainerVariableIO::GetDefaultParameters()
{
    return Parameters(R"({
        "prefix"            : "/",
        "list_of_variables" : []
    })");
}


} // namespace HDF5
} // namespace Kratos