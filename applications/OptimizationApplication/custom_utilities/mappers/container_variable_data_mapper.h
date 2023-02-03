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
        ContainerVariableDataHolderType& rOutputDataContainer,
        const ContainerVariableDataHolderType& rInputDataContainer)
    {
        KRATOS_ERROR << "Calling ContainerVariableDataMapper::Map. This needs to be "
                        "implemented in the derrived class.";
    }

    virtual void InverseMap(
        ContainerVariableDataHolderType& rOutputDataContainer,
        const ContainerVariableDataHolderType& rInputDataContainer)
    {
        KRATOS_ERROR << "Calling ContainerVariableDataMapper::Map. This needs to be "
                        "implemented in the derrived class.";
    }

    ///@}
};

///@}
} // namespace Kratos