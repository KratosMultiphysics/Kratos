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

class File;
struct WriteInfo;

namespace Internals
{
///@addtogroup HDF5Application
///@{
///@name Kratos Classes
///@{

/// Represents connectivities information of a single element or condition type in a mesh.
/**
 * Acts as the intermediary between the HDF5 file and the Kratos elements and conditions.
 */
template<class TContainerType>
class KRATOS_API(HDF5_APPLICATION) ConnectivitiesData
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = unsigned int;

    using EntityType = typename TContainerType::value_type;

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(ConnectivitiesData);

    ///@}
    ///@name Life Cycle
    ///@{

    ConnectivitiesData(
        const std::string& rPrefix,
        File::Pointer pFile);

    ///@}
    ///@name Operations
    ///@{

    /// Read data from a file.
    /**
     * Ensures valid element or condition data is read from the given path on
     * return. Previously stored data is replaced.
     */
    void Read(
        const std::string& rEntityName,
        NodesContainerType& rNodes,
        PropertiesContainerType& rProperties,
        TContainerType& rEntities);

    /// Write data to a file.
    void Write(
        const TContainerType& rEntities,
        const bool WriteProperties = true);

    ///@}

private:
    ///@name Member Variables
    ///@{

    File::Pointer mpFile;

    std::string mPrefix;

    ///@}
}; // class ConnectivitiesData

///@} // Kratos Classes
///@} addtogroup
} // namespace Internals.
} // namespace HDF5.
} // namespace Kratos.