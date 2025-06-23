//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#pragma once

// System includes
#include <type_traits>

// External includes

// Project includes
#include "utilities/data_type_traits.h"

namespace Kratos {

///@name Kratos Classes
///@{

class TensorAdaptorUtils {
public:
    ///@name Public static methods
    ///@{

    template<class TContainerType, class TContainerIO>
    static void GetShape(
        std::vector<int>& rShape,
        const TContainerType& rContainer,
        const TContainerIO& rContainerIO)
    {
        using return_type = typename TContainerIO::ReturnType;

        rShape.resize(DataTypeTraits<return_type>::Dimension + 1);

        return_type dummy_value{};
        if (!rContainer.empty()) {
            rContainerIO.GetValue(dummy_value, rContainer.front());
        }

        DataTypeTraits<return_type>::Shape(dummy_value, rShape.data() + 1, rShape.data() + rShape.size());
        rShape[0] = rContainer.size();
    }

    template<class TIndexType>
    static TIndexType GetFlatLength(
        TIndexType const * pShapeBegin,
        TIndexType const * pShapeEnd)
    {
        return std::accumulate(pShapeBegin, pShapeEnd, 1, std::multiplies<TIndexType>{});
    }

    ///@}
};

///@}

} // namespace Kratos