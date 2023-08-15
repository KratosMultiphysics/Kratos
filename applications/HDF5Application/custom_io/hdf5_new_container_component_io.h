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
//                  Suneth Warnakulasuriya, https://github.com/sunethwarna
//

#pragma once

// System includes
#include <string>
#include <vector>
#include <map>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/communicator.h"

// Application includes
#include "custom_io/hdf5_file.h"

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
template <class TContainerType, class TContainerDataIO, class... TComponents>
class KRATOS_API(HDF5_APPLICATION) NewContainerComponentIO
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(NewContainerComponentIO);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    NewContainerComponentIO(
        Parameters Settings,
        File::Pointer pFile);

    ///@}
    ///@name Operations
    ///@{

    void Write(
        const TContainerType& rContainer,
        Parameters Attributes);

    std::map<std::string, Parameters> Read(
        TContainerType& rContainer,
        Communicator& rCommunicator);

    ///@}

private:
    ///@name Private memeber variables
    ///@{

    File::Pointer mpFile;

    std::vector<std::string> mComponentNames;

    std::string mComponentPath;

    ///@}
    ///@name Private Operations
    ///@{

    template<class TComponentType>
    bool WriteComponentData(
        const std::string& rComponentName,
        const TContainerType& rContainer,
        Parameters Attributes,
        WriteInfo& rInfo);

    template<class TComponentType>
    bool ReadComponentData(
        const std::string& rComponentName,
        TContainerType& rContainer,
        Communicator& rCommunicator,
        std::map<std::string, Parameters>& rAttributesMap,
        const IndexType StartIndex,
        const IndexType BlockSize);

    ///@}

}; // class NewContainerComponentIO.


///@} // Kratos Classes
///@} addtogroup
} // namespace HDF5.
} // namespace Kratos.
