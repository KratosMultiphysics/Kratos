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
#include <string>
#include <sstream>

// External includes

// Project includes
#include "includes/define.h"
#include "utilities/data_type_traits.h"
#include "utilities/container_io_utils.h"

namespace Kratos {

///@name Kratos Classes
///@{

class KRATOS_API(KRATOS_CORE) TensorAdaptorUtils
{
public:
    ///@name Public static operations
    ///@{

    template<class TContainerType, class TIOType>
    static std::string Info(
        const TContainerType& rContainer,
        const TIOType& rIO,
        const DenseVector<unsigned int>& rShape)
    {
        std::stringstream msg;
        msg << "TensorAdaptor: " << rIO.Info() << " with " << rContainer.size()
            << " " << ModelPart::Container<TContainerType>::GetEntityName() << "(s) having shape "
            << rShape << " ].";
        return msg.str();
    }

    template<class TContainerType, class TIOType, class TPrimitiveDataType>
    static void InitializeData(
        const TContainerType& rContainer,
        const TIOType& rIO,
        const DenseVector<unsigned int>& rTensorAdaptorShape,
        DenseVector<TPrimitiveDataType>& rData)
    {
        using return_type = typename TIOType::ReturnType;

        KRATOS_ERROR_IF_NOT(rTensorAdaptorShape[0] == rContainer.size())
            << "First dimension of the initialized tensor adaptor mismatch with the container size [ "
            << Info(rContainer, rIO, rTensorAdaptorShape) << " ].\n";

        KRATOS_ERROR_IF_NOT(DataTypeTraits<return_type>::IsValidShape(rTensorAdaptorShape.data().begin() + 1, rTensorAdaptorShape.data().begin() + rTensorAdaptorShape.size()))
            << "Invalid data shape provided. [ data shape provided = " << rTensorAdaptorShape
            << ", max possible sizes in each dimension  = "
            << DataTypeTraits<return_type>::Shape(return_type{})
            << ", " << Info(rContainer, rIO, rTensorAdaptorShape) << " ].\n";

        rData.resize(DataTypeTraits<return_type>::Size(rTensorAdaptorShape.begin() + 1, rTensorAdaptorShape.end()) * rContainer.size(), false);
    }

    template<class TContainerType, class TIOType, class TPrimitiveDataType>
    static void CollectData(
        const TContainerType& rContainer,
        const TIOType& rIO,
        const DenseVector<unsigned int>& rTensorAdaptorShape,
        DenseVector<TPrimitiveDataType>& rData)
    {
        KRATOS_TRY

        // sanity checks
        KRATOS_ERROR_IF_NOT(rTensorAdaptorShape[0] == rContainer.size())
            << "First dimension of the initialized tensor adaptor mismatch with the container size [ "
            << Info(rContainer, rIO, rTensorAdaptorShape) << " ].\n";

        KRATOS_ERROR_IF_NOT(rData.size() == rContainer.size() * DataTypeTraits<typename TIOType::ReturnType>::Size(rTensorAdaptorShape.begin() + 1, rTensorAdaptorShape.end()))
            << "Data container size mismatch [ mData.size() = "
            << rData.size() << ", shape size = " << rContainer.size() * DataTypeTraits<typename TIOType::ReturnType>::Size(rTensorAdaptorShape.begin() + 1, rTensorAdaptorShape.end())
            << ", " << Info(rContainer, rIO, rTensorAdaptorShape) << " ].\n";

        if constexpr(TIOType::template IsAllowedContainer<TContainerType>) {
            CopyToContiguousArray(rContainer, rIO, rData.data().begin(), rTensorAdaptorShape.begin(), rTensorAdaptorShape.end());
        } else {
            KRATOS_ERROR << "The IOType " << rIO.Info() << " does not support containers with "
                         << ModelPart::Container<TContainerType>::GetEntityName() << ".\n";
        }

        KRATOS_CATCH("");
    }

    template<class TContainerType, class TIOType, class TPrimitiveDataType>
    static void StoreData(
        TContainerType& rContainer,
        const TIOType& rIO,
        const DenseVector<unsigned int>& rTensorAdaptorShape,
        const DenseVector<TPrimitiveDataType>& rData)
    {
        KRATOS_TRY

        // sanity checks
        KRATOS_ERROR_IF_NOT(rTensorAdaptorShape[0] == rContainer.size())
            << "First dimension of the initialized tensor adaptor mismatch with the container size [ "
            << Info(rContainer, rIO, rTensorAdaptorShape) << " ].\n";

        KRATOS_ERROR_IF_NOT(rData.size() == rContainer.size() * DataTypeTraits<typename TIOType::ReturnType>::Size(rTensorAdaptorShape.begin() + 1, rTensorAdaptorShape.end()))
            << "Data container size mismatch [ mData.size() = "
            << rData.size() << ", shape size = " << rContainer.size() * DataTypeTraits<typename TIOType::ReturnType>::Size(rTensorAdaptorShape.begin() + 1, rTensorAdaptorShape.end())
            << ", " << Info(rContainer, rIO, rTensorAdaptorShape) << " ].\n";

        const std::vector<unsigned int> shape(rTensorAdaptorShape.begin() + 1, rTensorAdaptorShape.end());

        if constexpr(TIOType::template IsAllowedContainer<TContainerType>) {
            CopyFromContiguousDataArray(rContainer, rIO, rData.data().begin(), rTensorAdaptorShape.begin(), rTensorAdaptorShape.end());
        } else {
            KRATOS_ERROR << "The IOType " << rIO.Info() << " does not support containers with "
                         << ModelPart::Container<TContainerType>::GetEntityName() << ".\n";
        }

        KRATOS_CATCH("");
    }

    ///@}
};

/// @}
} // namespace Kratos