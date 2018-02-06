//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: HDF5Application/license.txt
//
//  Main author:    Michael Andre, https://github.com/msandre
//

#if !defined(KRATOS_HDF5_DATA_VALUE_CONTAINER_IO_H_INCLUDED)
#define KRATOS_HDF5_DATA_VALUE_CONTAINER_IO_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "containers/data_value_container.h"

// Application includes
#include "hdf5_application_define.h"
#include "custom_io/hdf5_file.h"

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

/// A class for IO of a data value container in HDF5.
class DataValueContainerIO
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(DataValueContainerIO);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    DataValueContainerIO(std::string Prefix, File::Pointer pFile);

    ///@}
    ///@name Operations
    ///@{
    void ReadDataValueContainer(DataValueContainer& rData);

    void WriteDataValueContainer(DataValueContainer const& rData);
    ///@}

private:
    ///@name Member Variables
    ///@{
    std::string mPrefix;
    File::Pointer mpFile;
    ///@}
};

///@} // Kratos Classes
///@} addtogroup
} // namespace Internals.
} // namespace HDF5.
} // namespace Kratos.

#endif // KRATOS_HDF5_DATA_VALUE_CONTAINER_IO_H_INCLUDED defined
