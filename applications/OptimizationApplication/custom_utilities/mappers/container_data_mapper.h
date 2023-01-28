//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: HDF5Application/license.txt
//
//  Main author:     Suneth Warnakulasuriya
//

#pragma once

// System includes

// Project includes
#include "includes/define.h"

// Application includes
#include "custom_utilities/container_data.h"

namespace Kratos {

///@name Kratos Classes
///@{

class KRATOS_API(OPTIMIZATION_APPLICATION) ContainerDataMapper
{
public:
    ///@name Type definitions
    ///@{

    /// Pointer definition of ContainerData
    KRATOS_CLASS_POINTER_DEFINITION(ContainerDataMapper);

    ///@}
    ///@name Life cycle
    ///@{

    ContainerDataMapper() = default;

    virtual ~ContainerDataMapper() = default;

    ///@}
    ///@name Public operations

    virtual void Update()
    {
        KRATOS_ERROR << "Calling ContainerDataMapper::Update. This "
                        "needs to be implemented in the derrived class.";
    }

    virtual void Map(
        ContainerData& rOutputDataContainer,
        const ContainerData& rInputDataContainer)
    {
        KRATOS_ERROR << "Calling ContainerDataMapper::Map. This needs to be "
                        "implemented in the derrived class.";
    }

    virtual void InverseMap(
        ContainerData& rOutputDataContainer,
        const ContainerData& rInputDataContainer)
    {
        KRATOS_ERROR << "Calling ContainerDataMapper::Map. This needs to be "
                        "implemented in the derrived class.";
    }

    ///@}
};

///@}
} // namespace Kratos