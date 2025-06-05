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
    using ComponentType = char;

    constexpr static bool RequireAvailabilityCheck = true;

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

    template<class TEntityType>
    inline bool HasValue(
        const TEntityType& rEntity,
        const Flags& rFlag) const
    {
        return rEntity.IsDefined(rFlag);
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
    using ComponentType = TDataType;

    constexpr static bool RequireAvailabilityCheck = false;

    ///@}
    ///@name Life cycle
    ///@{

    HistoricalIO(const unsigned int StepIndex) : mStepIndex(StepIndex) {}

    ///@}
    ///@name Public operations
    ///@{

    template<class TDataType>
    inline const TDataType& GetValue(
        const Node& rNode,
        const Variable<TDataType>& rVariable,
        char& rTLS) const
    {
        return rNode.FastGetSolutionStepValue(rVariable, mStepIndex);
    }

    template<class TDataType>
    inline void SetValue(
        Node& rNode,
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
    using ComponentType = TDataType;

    constexpr static bool RequireAvailabilityCheck = true;

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

    template<class TEntityType, class TDataType>
    inline bool HasValue(
        const TEntityType& rEntity,
        const Variable<TDataType>& rVariable) const
    {
        return rEntity.Has(rVariable);
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
    using ComponentType = TDataType;

    constexpr static bool RequireAvailabilityCheck = false;

    ///@}
    ///@name Life cycle
    ///@{

    BossakIO(const double AlphaBossak) : mAlphaBossak(AlphaBossak) {}

    ///@}
    ///@name Public operations
    ///@{

    template<class TDataType>
    inline TDataType GetValue(
        const Node& rNode,
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
        Node& rNode,
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

class VertexHistoricalValueIO
{
public:
    ///@name Type definitions
    ///@{

    template<class TDataType>
    using TLSType = char;

    template<class TDataType>
    using ComponentType = TDataType;

    constexpr static bool RequireAvailabilityCheck = false;

    ///@}
    ///@name Life cycle
    ///@{

    VertexHistoricalValueIO(const unsigned int StepIndex) : mHistoricalIO(StepIndex) {}

    ///@}
    ///@name Public operations
    ///@{

    template<class TDataType>
    inline TDataType GetValue(
        const Detail::Vertex& rVertex,
        const Variable<TDataType>& rVariable,
        char& rTLS) const
    {
        return rVertex.GetValue(rVariable, mHistoricalIO);
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

private:
    ///@name Private member variables
    ///@{

    const HistoricalIO mHistoricalIO;

    ///@}
};

class VertexNonHistoricalValueIO
{
public:
    ///@name Type definitions
    ///@{

    template<class TDataType>
    using TLSType = char;

    template<class TDataType>
    using ComponentType = TDataType;

    constexpr static bool RequireAvailabilityCheck = false;

    ///@}
    ///@name Public operations
    ///@{

    template<class TDataType>
    inline TDataType GetValue(
        const Detail::Vertex& rVertex,
        const Variable<TDataType>& rVariable,
        char& rTLS) const
    {
        return rVertex.GetValue(rVariable, mNonHistoricalIO);
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

private:
    ///@name Private member variables
    ///@{

    const NonHistoricalIO mNonHistoricalIO{};

    ///@}
};

class GaussPointIO
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
    using ComponentType = Vector<typename DataTypeTraits<TDataType>::PrimitiveType>;

    constexpr static bool RequireAvailabilityCheck = false;

    ///@}
    ///@name Life cycle
    ///@{

    GaussPointIO(const ProcessInfo& rProcessInfo) : mrProcessInfo(rProcessInfo) {}

    ///@}
    ///@name Public operations
    ///@{

    template<class TEntityType, class TDataType>
    inline const ComponentType<TDataType>& GetValue(
        const TEntityType& rEntity,
        const Variable<TDataType>& rVariable,
        TLSType<TDataType>& rTLS) const
    {
        using list_traits = DataTypeTraits<std::vector<TDataType>>;

        using vector_traits = DataTypeTraits<Vector<typename list_traits::PrimitiveType>>;

        const_cast<TEntityType&>(rEntity).CalculateOnIntegrationPoints(rVariable, rTLS.mList, mrProcessInfo);

        const auto size = list_traits::Size(rTLS.mList);
        if (rTLS.mVector.size() != size) {
            rTLS.mVector.resize(size, false);
        }

        list_traits::CopyToContiguousData(vector_traits::GetContiguousData(rTLS.mVector), rTLS.mList);

        return rTLS.mVector;
    }

    template<class TEntityType, class TDataType>
    inline void SetValue(
        TEntityType& rEntity,
        const Variable<TDataType>& rVariable,
        const ComponentType<TDataType>& rValue) const
    {
        KRATOS_ERROR << "GaussPointIO does not support setting values from vertices";
    }

    ///@}

private:
    ///@name Private member variables
    ///@{

    const ProcessInfo& mrProcessInfo;

    ///@}
};

class NodesIO
{
public:
    ///@name Type definitions
    ///@{

    using IndexType = int;

    using ContainerType = ModelPart::NodesContainerType;

    ///@}
    ///@name Public operations
    ///@{

    inline void GetData(
        IndexType& rId,
        array_1d<double, 3>& rCoordinates,
        const ModelPart::NodeType& rNode) const
    {
        rId = rNode.Id();
        noalias(rCoordinates) = rNode.GetInitialPosition().Coordinates();
    }

    inline void AddPoint(
        ModelPart::NodesContainerType& rNodes,
        const IndexType Id,
        const array_1d<double, 3>& rCoordinates) const
    {
        auto p_node = Kratos::make_intrusive<ModelPart::NodeType>(Id, rCoordinates[0], rCoordinates[1], rCoordinates[2]);
        rNodes.push_back(p_node);
    }

    ///@}
};

class VertexIO
{
public:
    ///@name Type definitions
    ///@{

    using IndexType = int;

    using ContainerType = Detail::VertexContainerType;

    ///@}
    ///@name Public operations
    ///@{

    inline void GetData(
        IndexType& rId,
        array_1d<double, 3>& rCoordinates,
        const Detail::Vertex& rVertex) const
    {
        rId = rVertex.GetID();
        noalias(rCoordinates) = rVertex.Coordinates();
    }

    inline void AddPoint(
        Detail::VertexContainerType& rVertices,
        const IndexType Id,
        const array_1d<double, 3>& rCoordinates) const
    {
        KRATOS_ERROR << "Reading vertices is not allowed.";
    }

    ///@}
};

template<class TContainerIOType>
std::string GetContainerIOName()
{
    if constexpr(std::is_same_v<TContainerIOType, HistoricalIO>) {
        return "HISTORICAL";
    } else if constexpr(std::is_same_v<TContainerIOType, NonHistoricalIO>) {
        return "NONHISTORICAL";
    } else if constexpr(std::is_same_v<TContainerIOType, FlagIO>) {
        return "FLAGS";
    } else if constexpr(std::is_same_v<TContainerIOType, BossakIO>) {
        return "HISTORICAL_BOSSAK";
    } else if constexpr(std::is_same_v<TContainerIOType, VertexHistoricalValueIO>) {
        return "INTERPOLATED_HISTORICAL";
    } else if constexpr(std::is_same_v<TContainerIOType, VertexNonHistoricalValueIO>) {
        return "INTERPOLATED_NONHISTORICAL";
    } else if constexpr(std::is_same_v<TContainerIOType, GaussPointIO>) {
        return "GAUSS_POINT";
    } else {
        static_assert(!std::is_same_v<TContainerIOType, TContainerIOType>, "Unsupported container io type.");
    }
}

template<class TContainerType>
std::string GetContainerName()
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
void CopyToContiguousArray(
    const TContainerType& rContainer,
    const TComponentType& rComponent,
    const TContainerDataIO& rContainerDataIO,
    TDataType* pBegin,
    const IndexType Size)
{
    KRATOS_TRY

    using component_type = typename Internals::template ComponentTraits<TComponentType>::ValueType;

    using value_type = typename TContainerDataIO::template ComponentType<component_type>;

    using value_type_traits = DataTypeTraits<value_type>;

    static_assert(value_type_traits::IsContiguous, "Only contiguous data types are supported.");

    if (rContainer.empty()) {
        // do nothing if the container is empty.
        return;
    }

    // get the first item for sizing.
    typename TContainerDataIO::template TLSType<component_type> tls;
    const auto& initial_value = rContainerDataIO.GetValue(rContainer.front(), rComponent, tls);

    // get the stride from the first element to support dynamic types.
    const auto stride = value_type_traits::Size(initial_value);

    KRATOS_ERROR_IF_NOT(Size == rContainer.size() * stride)
        << "The contiguous array size mismatch with data in the container [ "
        << "Contiguous array size = " << Size << ", number of entities = "
        << rContainer.size() << ", data stride = " << stride << " ].";

    IndexPartition<unsigned int>(rContainer.size()).for_each(tls, [&rContainer, &rComponent, &rContainerDataIO, pBegin, stride](const auto Index, auto& rTLS) {
        const auto& value = rContainerDataIO.GetValue(*(rContainer.begin() + Index), rComponent, rTLS);
        TDataType const* p_value_begin = value_type_traits::GetContiguousData(value);
        auto p_subrange_begin = pBegin + Index * stride;
        std::copy(p_value_begin, p_value_begin + stride, p_subrange_begin);
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

    using value_type = typename TContainerDataIO::template ComponentType<typename Internals::template ComponentTraits<TComponentType>::ValueType>;

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
        TDataType const * p_subrange_begin = pBegin + Index * stride;
        std::copy(p_subrange_begin, p_subrange_begin + stride, p_value_begin);
        rContainerDataIO.SetValue(*(rContainer.begin() + Index), rComponent, rTLS);
    });

    KRATOS_CATCH("");
}

} // namespace Internals
} // namespace HDF5
} // namespace Kratos