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

#pragma once

// Application includes
#include "custom_io/hdf5_container_component_io.h"
#include "custom_utilities/container_io_utils.h"
#include "custom_utilities/vertex.h"

namespace Kratos
{
namespace HDF5
{

/** IO class for writing the coordinates (and IDs) of a set of vertices
 *
 *  @details Vertex coordinates are written to an (Nx3) matrix to
 *  <prefix>/POSITION. If "write_vertex_ids" is true, @ref{Vertex} IDs are
 *  written to an additional column (making the dataset an (Nx4) matrix).
 *  IDs are necessary for line outputs with MPI runs to identify the order
 *  of points along the line.
 *
 *  Default parameters:
 *  {
 *      "prefix"           : "/",
 *      "write_vertex_ids" : false
 *  }
 */
class KRATOS_API(HDF5_APPLICATION) VertexContainerCoordinateIO
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(VertexContainerCoordinateIO);

    ///@}
    ///@name Life Cycle
    ///@{

    VertexContainerCoordinateIO(
        Parameters Settings,
        File::Pointer pFile);

    ///@}
    ///@name Operations
    ///@{

    void Write(
        const Detail::VertexContainerType& rVertices,
        Parameters Attributes);

    ///@}

private:
    ///@name Private Operations
    ///@{

    void WriteWithoutIDs(
        const Detail::VertexContainerType& rVertices,
        Parameters Attributes);

    void WriteWithIDs(
        const Detail::VertexContainerType& rVertices,
        Parameters Attributes);

    ///@}
    ///@name Private Member Variables
    ///@{

    File::Pointer mpFile;

    std::string mPathPrefix;

    bool mWriteIDs;

    ///@}
}; // class VertexContainerCoordinateIO

/** IO class for writing variables of a set of vertices
 *
 *  Variables are written to (NxM) datasets at 'prefix'/variable.Name().
 *  N:  number of vertices
 *  M:  number of variable components
 *
 *  Default parameters:
 *  {
 *      "prefix"            : "/",
 *      "list_of_variables" : []
 *  }
 */
class KRATOS_API(HDF5_APPLICATION) VertexContainerVariableIO: protected ContainerComponentIO<
                                                                                Detail::VertexContainerType,
                                                                                Internals::VertexValueIO,
                                                                                Variable<int>,
                                                                                Variable<double>,
                                                                                Variable<array_1d<double, 3>>,
                                                                                Variable<array_1d<double, 4>>,
                                                                                Variable<array_1d<double, 6>>,
                                                                                Variable<array_1d<double, 9>>,
                                                                                Variable<Kratos::Vector>,
                                                                                Variable<Kratos::Matrix>>
{
public:
    ///@name Type Definitions
    ///@{

    using BaseType = ContainerComponentIO<
                                Detail::VertexContainerType,
                                Internals::VertexValueIO,
                                Variable<int>,
                                Variable<double>,
                                Variable<array_1d<double, 3>>,
                                Variable<array_1d<double, 4>>,
                                Variable<array_1d<double, 6>>,
                                Variable<array_1d<double, 9>>,
                                Variable<Kratos::Vector>,
                                Variable<Kratos::Matrix>>;

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(VertexContainerVariableIO);

    ///@}
    ///@name Life Cycle
    ///@{

    VertexContainerVariableIO(
        Parameters Settings,
        File::Pointer pFile);

    ///@}
    ///@name Operations
    ///@{

    void Write(
        const Detail::VertexContainerType& rVertices,
        Parameters Attributes);

    ///@}

}; // class VertexContainerVariableIO

} // namespace HDF5
} // namespace Kratos