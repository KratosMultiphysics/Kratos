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

class KRATOS_API(KRATOS_CORE) TensorAdaptorUtils {
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

    template<template<class> class TContainerIOType, class TDataType, class TContainerPointerType, class... TArgs>
    static typename TContainerIOType<TDataType>::Pointer const * InitializeAndGetContainerIO(
        DenseVector<int>& rContainerShape,
        TContainerPointerType pContainer,
        const Variable<TDataType>& rVariable,
        TArgs&&... rArgs)
    {
        using return_type = typename TContainerIOType<TDataType>::ReturnType;

        // creating the io
        auto p_container_io = Kratos::make_shared<TContainerIOType<TDataType>>(rVariable, rArgs...);

        if (rUserGivenShape.size() != 1 || rUserGivenShape[0] != -1) {
            KRATOS_ERROR_IF_NOT(DataTypeTraits<return_type>::Dimension + 1 == rUserGivenShape.size())
                << "Dimensions mismatch for " << pVariable->Name()
                << " [ Required dimensions by variable = " << DataTypeTraits<return_type>::Dimension
                << ", shape = " << rUserGivenShape << " ].\n";
        }

        TensorAdaptorUtils::GetShape(rContainerShape, *pContainer, *p_container_io);
        if constexpr (!DataTypeTraits<TDataType>::IsDynamic) {
            // now we know the shape exactly because, rContainerShape should have the
            // exact sizes, hence we can check whether user has provided the shape
            // correctly. These types are double, array3, array4, ...

            if (rUserGivenShape.size() != 1 || rUserGivenShape[0] != -1) {
                for (IndexType i = 0; i < rUserGivenShape.size(); ++i) {
                    KRATOS_ERROR_IF(rUserGivenShape[i] != -1 && rUserGivenShape[i] != rContainerShape[i])
                        << "Shape mismatch for " << pVariable->Name()
                        << " [ Required variable shape = " << rContainerShape
                        << ", shape = " << rUserGivenShape << " ].\n";
                }
            }
        }
        else {
            // here we may not know the exact sizes. so we use the user given sizes
            for (IndexType i = 0; i < rUserGivenShape.size(); ++i) {
                if (rUserGivenShape[i] != -1) {
                    rContainerShape[i] = rUserGivenShape[i];
                }
            }
        }

        return p_container_io;
    }

    ///@}
};

///@}

} // namespace Kratos