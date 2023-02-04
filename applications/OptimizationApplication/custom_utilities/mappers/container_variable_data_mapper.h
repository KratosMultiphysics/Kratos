//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: HDF5Application/license.txt
//
//  Main author:     Suneth Warnakulasuriya
//


#pragma once

// System includes
#include <string>

// Project includes
#include "includes/define.h"

// Application includes
#include "custom_utilities/container_variable_data_holder/container_variable_data_holder_base.h"

namespace Kratos {

///@name Kratos Classes
///@{

template<class TContainerType>
class KRATOS_API(OPTIMIZATION_APPLICATION) ContainerVariableDataMapper
{
public:
    ///@name Type definitions
    ///@{

    using ContainerVariableDataHolderType = ContainerVariableDataHolderBase<TContainerType>;

    /// Pointer definition of ContainerData
    KRATOS_CLASS_POINTER_DEFINITION(ContainerVariableDataMapper);

    ///@}
    ///@name Life cycle
    ///@{

    ContainerVariableDataMapper() = default;

    virtual ~ContainerVariableDataMapper() = default;

    ///@}
    ///@name Public operations

    virtual void Update()
    {
        KRATOS_ERROR << "Calling ContainerVariableDataMapper::Update. This "
                        "needs to be implemented in the derrived class.";
    }

    virtual void Map(
        const ContainerVariableDataHolderType& rOriginDataContainer,
        ContainerVariableDataHolderType& rDestinationDataContainer) const
    {
        KRATOS_ERROR << "Calling ContainerVariableDataMapper::Map. This needs to be "
                        "implemented in the derrived class.";
    }

    virtual void InverseMap(
        ContainerVariableDataHolderType& rOriginDataContainer,
        const ContainerVariableDataHolderType& rDestinationDataContainer) const
    {
        KRATOS_ERROR << "Calling ContainerVariableDataMapper::Map. This needs to be "
                        "implemented in the derrived class.";
    }

    virtual std::string Info() const
    {
        KRATOS_ERROR << "Calling ContainerVariableDataMapper::Info. This needs to be "
                        "implemented in the derrived class.";
        return "";
    }

    ///@}
};

///@}
///@name Input and output
///@{

/// output stream function
template<class TContainerType>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const ContainerVariableDataMapper<TContainerType>& rThis)
{
    return rOStream << rThis.Info();
}

///@}
} // namespace Kratos