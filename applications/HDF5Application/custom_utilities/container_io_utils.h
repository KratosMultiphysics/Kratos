//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  license: HDF5Application/license.txt
//
//  Main author:    Suneth Warnakulasuriya
//

#pragma once

// System includes
#include <vector>
#include <variant>
#include <type_traits>

// Project includes
#include "containers/flags.h"
#include "includes/model_part.h"
#include "utilities/data_type_traits.h"

// Application includes
#include "custom_utilities/data_type_utilities.h"
#include "custom_utilities/vertex.h"

namespace Kratos
{
namespace HDF5
{
namespace Internals
{

class FlagIO
{
public:
    ///@name Type definitions
    ///@{

    template<class TDataType>
    using TLSType = char;

    template<class TDataType>
    using ComponentDataType = char;

    ///@}
    ///@name Public operations
    ///@{

    template<class TEntityType>
    inline char GetValue(
        const TEntityType& rEntity,
        const Flags& rFlag,
        char& rTLS) const
    {
        return rEntity.Is(rFlag);
    }

    template<class TEntityType>
    inline void SetValue(
        TEntityType& rEntity,
        const Flags& rFlag,
        const char rValue) const
    {
        rEntity.Set(rFlag, static_cast<bool>(rValue));
    }

    ///@}
};

class HistoricalIO
{
public:
    ///@name Type definitions
    ///@{

    template<class TDataType>
    using TLSType = char;

    template<class TDataType>
    using ComponentDataType = TDataType;

    ///@}
    ///@name Life cycle
    ///@{

    HistoricalIO(const unsigned int StepIndex) : mStepIndex(StepIndex) {}

    ///@}
    ///@name Public operations
    ///@{

    template<class TDataType>
    inline const TDataType& GetValue(
        const ModelPart::NodeType& rNode,
        const Variable<TDataType>& rVariable,
        char& rTLS) const
    {
        return rNode.FastGetSolutionStepValue(rVariable, mStepIndex);
    }

    template<class TDataType>
    inline void SetValue(
        ModelPart::NodeType& rNode,
        const Variable<TDataType>& rVariable,
        const TDataType& rValue) const
    {
        rNode.FastGetSolutionStepValue(rVariable, mStepIndex) = rValue;
    }

    ///@}

private:
    ///@name Private member variables
    ///@{

    const unsigned int mStepIndex;

    ///@}
};

class NonHistoricalIO
{
public:
    ///@name Type definitions
    ///@{

    template<class TDataType>
    using TLSType = char;

    template<class TDataType>
    using ComponentDataType = TDataType;

    ///@}
    ///@name Public operations
    ///@{

    template<class TEntityType, class TDataType>
    inline const TDataType& GetValue(
        const TEntityType& rEntity,
        const Variable<TDataType>& rVariable,
        char& rTLS) const
    {
        return rEntity.GetValue(rVariable);
    }

    template<class TEntityType, class TDataType>
    inline void SetValue(
        TEntityType& rEntity,
        const Variable<TDataType>& rVariable,
        const TDataType& rValue) const
    {
        rEntity.SetValue(rVariable, rValue);
    }

    ///@}
};

class BossakIO
{
public:
    ///@name Type definitions
    ///@{

    template<class TDataType>
    using TLSType = char;

    template<class TDataType>
    using ComponentDataType = TDataType;

    ///@}
    ///@name Life cycle
    ///@{

    BossakIO(const double AlphaBossak) : mAlphaBossak(AlphaBossak) {}

    ///@}
    ///@name Public operations
    ///@{

    template<class TDataType>
    inline TDataType GetValue(
        const ModelPart::NodeType& rNode,
        const Variable<TDataType>& rVariable,
        char& rTLS) const
    {
        if (rVariable == ACCELERATION) {
            return (1.0 - mAlphaBossak) * rNode.FastGetSolutionStepValue(rVariable, 0) + mAlphaBossak * rNode.FastGetSolutionStepValue(rVariable, 1);
        } else {
            return rNode.FastGetSolutionStepValue(rVariable);
        }
    }

    template<class TDataType>
    inline void SetValue(
        ModelPart::NodeType& rNode,
        const Variable<TDataType>& rVariable,
        const TDataType& rValue) const
    {
        rNode.FastGetSolutionStepValue(rVariable) =  rValue;
    }

    ///@}

private:
    ///@name Private member variables
    ///@{

    const double mAlphaBossak;

    ///@}
};

class VertexValueIO
{
public:
    ///@name Type definitions
    ///@{

    template<class TDataType>
    using TLSType = char;

    template<class TDataType>
    using ComponentDataType = TDataType;

    ///@}
    ///@name Public operations
    ///@{

    template<class TDataType>
    inline TDataType GetValue(
        const Detail::Vertex& rVertex,
        const Variable<TDataType>& rVariable,
        char& rTLS) const
    {
        return rVertex.GetValue(rVariable);
    }

    template<class TDataType>
    inline void SetValue(
        Detail::Vertex& rEntity,
        const Variable<TDataType>& rVariable,
        const TDataType& rValue) const
    {
        KRATOS_ERROR << "VertexValueIO does not support setting values from vertices";
    }

    ///@}
};

class GaussPointValueIO
{
public:
    ///@name Calss definitions
    ///@{

    template<class TDataType>
    struct TLSStruct {
        std::vector<TDataType> mList;
        Vector<typename DataTypeTraits<TDataType>::PrimitiveType> mVector;
    };

    ///@}
    ///@name Type definitions
    ///@{


    template<class TDataType>
    using TLSType = TLSStruct<TDataType>;

    template<class TDataType>
    using ComponentDataType = Vector<typename DataTypeTraits<TDataType>::PrimitiveType>;

    ///@}
    ///@name Life cycle
    ///@{

    GaussPointValueIO(const ProcessInfo& rProcessInfo) : mrProcessInfo(rProcessInfo) {}

    ///@}
    ///@name Public operations
    ///@{

    template<class TEntityType, class TDataType>
    inline const ComponentDataType<TDataType>& GetValue(
        const TEntityType& rEntity,
        const Variable<TDataType>& rVariable,
        TLSType<TDataType>& rTLS) const
    {
        // Since the array_1d<double, 4> and array_1d<double 9> interfaces are not present in
        // elements and conditions.
        if constexpr(!std::is_same_v<TDataType, array_1d<double, 4>> && !std::is_same_v<TDataType, array_1d<double, 9>>) {

            using list_traits = DataTypeTraits<std::vector<TDataType>>;

            using vector_traits = DataTypeTraits<Vector<typename list_traits::PrimitiveType>>;

            const_cast<TEntityType&>(rEntity).CalculateOnIntegrationPoints(rVariable, rTLS.mList, mrProcessInfo);

            const auto size = list_traits::Size(rTLS.mList);
            if (rTLS.mVector.size() != size) {
                rTLS.mVector.resize(size, false);
            }

            list_traits::CopyToContiguousData(vector_traits::GetContiguousData(rTLS.mVector), rTLS.mList);

            return rTLS.mVector;
        } else {
            KRATOS_ERROR << "The gauss point interface in elements/conditions for Variable<array_1d<double, 4>> and Variable<array_1d<double, 9>> are not defined.";
            return rTLS.mVector;
        }
    }

    template<class TEntityType, class TDataType>
    inline void SetValue(
        TEntityType& rEntity,
        const Variable<TDataType>& rVariable,
        const ComponentDataType<TDataType>& rValue) const
    {
        KRATOS_ERROR << "GaussPointValueIO does not support setting values from vertices";
    }

    ///@}

private:
    ///@name Private member variables
    ///@{

    const ProcessInfo& mrProcessInfo;

    ///@}
};

template<class TContainerIOType>
std::string GetContainerIOType()
{
    if constexpr(std::is_same_v<TContainerIOType, HistoricalIO>) {
        return "HISTORICAL";
    } else if constexpr(std::is_same_v<TContainerIOType, NonHistoricalIO>) {
        return "NONHISTORICAL";
    } else if constexpr(std::is_same_v<TContainerIOType, FlagIO>) {
        return "FLAGS";
    } else if constexpr(std::is_same_v<TContainerIOType, BossakIO>) {
        return "HISTORICAL_BOSSAK";
    } else if constexpr(std::is_same_v<TContainerIOType, VertexValueIO>) {
        return "INTERPOLATED";
    } else if constexpr(std::is_same_v<TContainerIOType, GaussPointValueIO>) {
        return "GAUSS_POINT";
    } else {
        static_assert(!std::is_same_v<TContainerIOType, TContainerIOType>, "Unsupported container io type.");
    }
}

template<class TContainerType>
std::string GetContainerType()
{
    if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
        return "NODES";
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>) {
        return "CONDITIONS";
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ElementsContainerType>) {
        return "ELEMENTS";
    } else if constexpr(std::is_same_v<TContainerType, Detail::VertexContainerType>) {
        return "VERTICES";
    } else {
        static_assert(!std::is_same_v<TContainerType, TContainerType>, "Unsupported container type.");
    }
}

template<class TContainerType>
auto& GetLocalContainer(ModelPart& rModelPart)
{
    if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
        return rModelPart.GetCommunicator().LocalMesh().Nodes();
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>) {
        return rModelPart.GetCommunicator().LocalMesh().Conditions();
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ElementsContainerType>) {
        return rModelPart.GetCommunicator().LocalMesh().Elements();
    } else {
        static_assert(!std::is_same_v<TContainerType, TContainerType>, "Unsupported container type.");
    }
}

template<class TContainerType>
const auto& GetLocalContainer(const ModelPart& rModelPart)
{
    if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
        return rModelPart.GetCommunicator().LocalMesh().Nodes();
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>) {
        return rModelPart.GetCommunicator().LocalMesh().Conditions();
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ElementsContainerType>) {
        return rModelPart.GetCommunicator().LocalMesh().Elements();
    } else {
        static_assert(!std::is_same_v<TContainerType, TContainerType>, "Unsupported container type.");
    }
}

template<class TContainerType, class TComponentType, class TContainerDataIO, class TDataType>
void CopyToContiguousDataArray(
    const TContainerType& rContainer,
    const TComponentType& rComponent,
    const TContainerDataIO& rContainerDataIO,
    TDataType* pBegin)
{
    KRATOS_TRY

    using component_data_type = typename Internals::template ComponentTraits<TComponentType>::ValueType;

    using value_type = typename TContainerDataIO::template ComponentDataType<component_data_type>;

    using value_type_traits = DataTypeTraits<value_type>;

    static_assert(value_type_traits::IsContiguous, "Only contiguous data types are supported.");

    if (rContainer.empty()) {
        // do nothing if the container is empty.
        return;
    }

    // get the first item for sizing.
    typename TContainerDataIO::template TLSType<component_data_type> tls;
    const auto& initial_value = rContainerDataIO.GetValue(rContainer.front(), rComponent, tls);

    // get the stride from the first element to support dynamic types.
    const auto stride = value_type_traits::Size(initial_value);

    IndexPartition<unsigned int>(rContainer.size()).for_each(tls, [&rContainer, &rComponent, &rContainerDataIO, pBegin, stride](const auto Index, auto& rTLS) {
        const auto& value = rContainerDataIO.GetValue(*(rContainer.begin() + Index), rComponent, rTLS);
        TDataType const* p_value_begin = value_type_traits::GetContiguousData(value);
        auto p_array_start = pBegin + Index * stride;
        std::copy(p_value_begin, p_value_begin + stride, p_array_start);
    });

    KRATOS_CATCH("");
}

template<class TContainerType, class TComponentType, class TContainerDataIO, class TDataType>
void CopyFromContiguousDataArray(
    TContainerType& rContainer,
    const TComponentType& rComponent,
    const TContainerDataIO& rContainerDataIO,
    TDataType const* pBegin,
    const std::vector<unsigned int>& rShape)
{
    KRATOS_TRY

    using value_type = typename TContainerDataIO::template ComponentDataType<typename Internals::template ComponentTraits<TComponentType>::ValueType>;

    using value_type_traits = DataTypeTraits<value_type>;

    static_assert(value_type_traits::IsContiguous, "Only contiguous data types are supported.");

    if (rContainer.empty()) {
        // do nothing if the container is empty.
        return;
    }

    value_type tls_prototype;
    value_type_traits::Reshape(tls_prototype, rShape);

    const auto stride = value_type_traits::Size(tls_prototype);

    IndexPartition<unsigned int>(rContainer.size()).for_each(tls_prototype, [&rContainer, &rComponent, &rContainerDataIO, pBegin, stride](const auto Index, auto& rTLS) {
        TDataType * p_value_begin = value_type_traits::GetContiguousData(rTLS);
        TDataType const * p_array_start = pBegin + Index * stride;
        std::copy(p_array_start, p_array_start + stride, p_value_begin);
        rContainerDataIO.SetValue(*(rContainer.begin() + Index), rComponent, rTLS);
    });

    KRATOS_CATCH("");
}

} // namespace Internals
} // namespace HDF5
} // namespace Kratos