//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  license: HDF5Application/license.txt
//
//  Main author:    Michael Andre, https://github.com/msandre
//                  Suneth Warnakulasuriya
//

#pragma once

// System includes
#include <string>

// External includes

// Project includes
#include "includes/define.h"

// Application includes
#include "custom_io/hdf5_file.h"
#include "hdf5_application_define.h"

namespace Kratos
{
namespace HDF5
{
namespace Internals
{
///@addtogroup HDF5Application
///@{
///@name Kratos Classes
///@{

/// A class representing points in a mesh.
/**
 * @tparam TContainerDataIO A IO class which have the @p ContainerType defined in public scope
 *                          with @p SetValue, @p GetValue methods implemented.
*/
template<class TContainerDataIO>
class PointsData
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(PointsData);

    ///@}
    ///@name Life Cycle
    ///@{

    PointsData(
        const std::string& rPrefix,
        File::Pointer mpFile);

    ///@}
    ///@name Operations
    ///@{

    Parameters Read(
        typename TContainerDataIO::ContainerType& rContainer,
        const TContainerDataIO& rContainerDataIO);

    void Write(
        const typename TContainerDataIO::ContainerType& rContainer,
        const TContainerDataIO& rContainerDataIO,
        const Parameters Attributes);

    ///@}

protected:
    ///@name Protected member Variables
    ///@{

    File::Pointer mpFile;

    std::string mPrefix;

    ///@}
}; // class PointsData

///@} // Kratos Classes
///@} addtogroup
} // namespace Internals.
} // namespace HDF5.
} // namespace Kratos.