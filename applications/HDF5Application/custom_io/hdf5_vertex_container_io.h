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

#ifndef KRATOS_HDF5APPLICATION_VERTEX_CONTAINER_IO_H
#define KRATOS_HDF5APPLICATION_VERTEX_CONTAINER_IO_H

// Application includes
#include "custom_io/hdf5_container_component_io.h"
#include "custom_utilities/vertex.h"


namespace Kratos
{
namespace HDF5
{


/** Interface for writing @ref{Vertex} data
 *
 *  Default parameters:
 *  {
 *      "prefix" : "/"
 *  }
 * 
 *  @note Derived classes take care of writing coordinates and variables of the vertices
 *  separately (see @ref{VertexContainerCoordinateIO} and @ref{VertexContainerVariableIO} respectively).
 */
class KRATOS_API(HDF5_APPLICATION) VertexContainerIO : protected ContainerComponentIO<Detail::VertexContainerType,
                                                         Detail::Vertex,
                                                         Variable<array_1d<double, 3>>,
                                                         Variable<double>,
                                                         Variable<int>,
                                                         Variable<Vector<double>>,
                                                         Variable<Matrix<double>>>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(VertexContainerIO);

    ///@}
    ///@name Life Cycle
    ///@{

    VertexContainerIO(Parameters parameters, File::Pointer pFile);

    virtual ~VertexContainerIO() {}

    ///@}
    ///@name Operations
    ///@{

    virtual void Write(const Detail::VertexContainerType& rVertices) = 0;

    ///@}
private:
    using BaseType = ContainerComponentIO<Detail::VertexContainerType,
                                          Detail::Vertex,
                                          Variable<array_1d<double, 3>>,
                                          Variable<double>,
                                          Variable<int>,
                                          Variable<Vector<double>>,
                                          Variable<Matrix<double>>>;

    static Parameters GetDefaultParameters();

    static Parameters FormatParameters(Parameters parameters);

    static std::string GetPath(Parameters parameters);

}; // class VertexContainerIO




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
class KRATOS_API(HDF5_APPLICATION) VertexContainerCoordinateIO final : public VertexContainerIO
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(VertexContainerCoordinateIO);

    ///@}
    ///@name Life Cycle
    ///@{

    VertexContainerCoordinateIO(Parameters parameters, File::Pointer pFile);

    ///@}
    ///@name Operations
    ///@{

    virtual void Write(const Detail::VertexContainerType& rVertices) override;

    ///@}

private:
    ///@name Private Operations
    ///@{

    void WriteWithoutIDs(const Detail::VertexContainerType& rVertices);

    void WriteWithIDs(const Detail::VertexContainerType& rVertices);

    // Format input parameters for initializing the base class
    static Parameters FormatParameters(Parameters parameters);

    static Parameters GetDefaultParameters();

    ///@}
    ///@name Private Member Variables
    ///@{

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
class KRATOS_API(HDF5_APPLICATION) VertexContainerVariableIO final : public VertexContainerIO
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(VertexContainerVariableIO);

    ///@}
    ///@name Life Cycle
    ///@{

    using VertexContainerIO::VertexContainerIO;

    ///@}
    ///@name Operations
    ///@{

    virtual void Write(const Detail::VertexContainerType& rVertices) override;

    ///@}

private:
    static Parameters GetDefaultParameters();
}; // class VertexContainerVariableIO




} // namespace HDF5
} // namespace Kratos


#endif