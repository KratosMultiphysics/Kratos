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
//                  Suneth Warnakulasuriya, https://github.com/sunethwarna
//

#if !defined(KRATOS_HDF5_CONTAINER_COMPONENT_IO_H_INCLUDED)
#define KRATOS_HDF5_CONTAINER_COMPONENT_IO_H_INCLUDED

// System includes
#include <string>
#include <vector>

// External includes

// Project includes
#include "includes/define.h"

// Application includes
#include "custom_io/hdf5_file.h"
#include "hdf5_application_define.h"

namespace Kratos
{
class Parameters;

namespace HDF5
{
///@addtogroup HDF5Application
///@{

///@name Kratos Classes
///@{

/// A class for IO of element data in HDF5.
template <typename TContainerType, typename TContainerItemType, typename... TComponents>
class ContainerComponentIO
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(ContainerComponentIO);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    ContainerComponentIO(Parameters Settings, File::Pointer pFile, const std::string& rComponentPath);

    ///@}
    ///@name Operations
    ///@{

    ///@}

protected:
    ///@name Protected Operations
    ///@{

    void WriteContainerComponents(TContainerType const& rContainerItems);

    void ReadContainerComponents(TContainerType& rContainerItems, Communicator& rCommunicator);

    ///@}

private:
    ///@name Member Variables
    ///@{
    File::Pointer mpFile;
    std::string mComponentPath;
    std::vector<std::string> mComponentNames;
    ///@}
    ///@name Private Operations
    ///@{

    template <typename... Targs>
    void WriteRegisteredComponent(const std::string& rComponentName, Targs&... args);

    template <typename... Targs>
    void ReadRegisteredComponent(const std::string& rComponentName, Targs&... args);    

    ///@}

}; // class ContainerComponentIO.

///@} // Kratos Classes
///@} addtogroup
} // namespace HDF5.
} // namespace Kratos.

#endif // KRATOS_HDF5_CONTAINER_COMPONENT_IO_H_INCLUDED defined
