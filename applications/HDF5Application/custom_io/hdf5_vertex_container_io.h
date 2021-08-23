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


class VertexContainerIO : private ContainerComponentIO<Detail::VertexContainerType,
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

    ///@}
    ///@name Operations
    ///@{

    void WriteCoordinates(const Detail::VertexContainerType& rVertices);

    void WriteVariables(const Detail::VertexContainerType& rVertices);

    static Parameters GetDefaultParameters();

    ///@}
private:
    using BaseType = ContainerComponentIO<Detail::VertexContainerType,
                                          Detail::Vertex,
                                          Variable<array_1d<double, 3>>,
                                          Variable<double>,
                                          Variable<int>,
                                          Variable<Vector<double>>,
                                          Variable<Matrix<double>>>;

    ///@name Member Variables
    ///@{

    std::string mCoordinatesPath;

    ///@}
    ///@name Private Operations
    ///@{

    // Format the input parameters for use in the base class
    static Parameters FormatParameters(Parameters parameters);

    // Get 'component path' from parameters required by the base class
    static std::string GetComponentPath(Parameters parameters);

    ///@}
};


} // namespace HDF5
} // namespace Kratos


#endif