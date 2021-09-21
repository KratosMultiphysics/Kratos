//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: HDF5Application/license.txt
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
            VertexContainerIO::GetComponentPath(parameters))
{
    KRATOS_TRY

    // Parse input parameters
    parameters.ValidateAndAssignDefaults(VertexContainerIO::GetDefaultParameters());
    
    const std::string prefix = parameters["group_prefix"].GetString(); 
    mCoordinatesPath = prefix + parameters["coordinates_path"].GetString();

    // Check/create prefix path
    // TODO: this doesn't take care of multi-level prefixes
    if (!mpFile->IsGroup(prefix)) {
        KRATOS_ERROR_IF(mpFile->HasPath(prefix)) << "'prefix' points to an existing dataset!";
        mpFile->CreateGroup(prefix);
    }

    KRATOS_CATCH("");
}


void VertexContainerIO::WriteCoordinates(const Detail::VertexContainerType& rVertices)
{
    KRATOS_TRY

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
    mpFile->WriteDataSet(mCoordinatesPath,
                         buffer,
                         write_info);

    KRATOS_CATCH("");
}


void VertexContainerIO::WriteCoordinatesAndIDs(const Detail::VertexContainerType& rVertices)
{
        KRATOS_TRY

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
    mpFile->WriteDataSet(mCoordinatesPath,
                         buffer,
                         write_info);

    KRATOS_CATCH("");
}


void VertexContainerIO::WriteVariables(const Detail::VertexContainerType& rVertices)
{
    KRATOS_TRY

    // Check/create variables path
    // TODO: this doesn't take care of multi-level prefixes
    if (!mpFile->IsGroup(mComponentPath)) {
        KRATOS_ERROR_IF(mpFile->HasPath(mComponentPath)) << "'variables_path' points to an existing dataset!";
        mpFile->CreateGroup(mComponentPath);
    }

    this->WriteContainerComponents(rVertices);

    KRATOS_CATCH("");
}


Parameters VertexContainerIO::GetDefaultParameters()
{
    return Parameters(R"({
        "group_prefix" : "/",
        "coordinates_path" : "/coordinates",
        "variables_path" : "/variables",
        "list_of_variables" : []
    })");
}


Parameters VertexContainerIO::FormatParameters(Parameters parameters)
{
    KRATOS_TRY

    Parameters output = parameters.Clone();
    output.ValidateAndAssignDefaults(VertexContainerIO::GetDefaultParameters());
    output.RemoveValue("coordinates_path");
    output.RemoveValue("variables_path");

    output.RemoveValue("group_prefix");
    output.AddValue("prefix", parameters["group_prefix"].Clone());

    return output;

    KRATOS_CATCH("");
}


std::string VertexContainerIO::GetComponentPath(Parameters parameters)
{
    KRATOS_TRY

    parameters.ValidateAndAssignDefaults(VertexContainerIO::GetDefaultParameters());
    return parameters["variables_path"].GetString();

    KRATOS_CATCH("");
}


} // namespace HDF5
} // namespace Kratos