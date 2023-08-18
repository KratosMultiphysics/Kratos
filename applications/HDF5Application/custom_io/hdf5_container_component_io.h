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
#include "includes/model_part.h"

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
class KRATOS_API(HDF5_APPLICATION) ContainerComponentIO
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(ContainerComponentIO);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    ContainerComponentIO(
        Parameters Settings,
        File::Pointer pFile);

    ///@}
    ///@name Operations
    ///@{

    void Write(
        const ModelPart& rModelPart,
        const TContainerDataIO& rContainerDataIO,
        const Parameters Attributes);

    std::map<std::string, Parameters> Read(
        ModelPart& rModelPart,
        const TContainerDataIO& rContainerDataIO);

    void Write(
        const TContainerType& rLocalContainer,
        const TContainerDataIO& rContainerDataIO,
        const Parameters Attributes);

    std::map<std::string, Parameters> Read(
        TContainerType& rLocalContainer,
        const TContainerDataIO& rContainerDataIO,
        Communicator& rCommunicator);

    std::map<std::string, Parameters> ReadAttributes();

    ///@}

protected:
    ///@name Protected member variables
    ///@{

    File::Pointer mpFile;

    std::vector<std::string> mComponentNames;

    std::string mComponentPath;

    ///@}

private:
    ///@name Private Operations
    ///@{

    template<class TComponentType>
    bool WriteComponentData(
        const std::string& rComponentName,
        const TContainerDataIO& rContainerDataIO,
        const TContainerType& rLocalContainer,
        Parameters Attributes,
        WriteInfo& rInfo);

    template<class TComponentType>
    bool ReadComponentData(
        const std::string& rComponentName,
        const TContainerDataIO& rContainerDataIO,
        TContainerType& rLocalContainer,
        Communicator& rCommunicator,
        std::map<std::string, Parameters>& rAttributesMap,
        const IndexType StartIndex,
        const IndexType BlockSize);

    ///@}

}; // class ContainerComponentIO.


///@} // Kratos Classes
///@} addtogroup
} // namespace HDF5.
} // namespace Kratos.
