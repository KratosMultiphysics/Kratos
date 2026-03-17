//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         SystemIdentificationApplication/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/model_part.h"

// Application includes

namespace Kratos {
///@name Kratos Classes
///@{

class KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) ControlUtils
{
public:
    ///@name Public static operations
    ///@{

    template<class TContainerType>
    static void AssignEquivalentProperties(
        TContainerType& rSourceContainer,
        TContainerType& rDestinationContainer);

    ///@}
};

///@} // Kratos Classes

} /* namespace Kratos.*/