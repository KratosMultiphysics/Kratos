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

    template<class TContainerType>
    static double NormInf(const ContainerVariableDataHolderBase<TContainerType>& rContainer);

    template<class TContainerType>
    static double EntityMaxNormL2(const ContainerVariableDataHolderBase<TContainerType>& rContainer);

    template<class TContainerType>
    static double InnerProduct(
        const ContainerVariableDataHolderBase<TContainerType>& rContainer1,
        const ContainerVariableDataHolderBase<TContainerType>& rContainer2);

    ///@}
};

///@}
}