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
#include "custom_utilities/container_variable_data_holder.h"

namespace Kratos
{

///@name Kratos Classes
///@{

class KRATOS_API(OPTIMIZATION_APPLICATION) ContainerVariableDataHolderUtils
{
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    ///@}
    ///@name Static operations
    ///@{

    static double NormInf(const ContainerVariableDataHolderBase& rContainer);

    static double EntityMaxNormL2(const ContainerVariableDataHolderBase& rContainer);

    static double InnerProduct(
        const ContainerVariableDataHolderBase& rContainer1,
        const ContainerVariableDataHolderBase& rContainer2);

    ///@}
};

///@}
}