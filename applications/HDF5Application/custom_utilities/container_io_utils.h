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

/**
 * @brief This enum is used to indicate whether a given variable in a list of entities are always available,
 *        always not available, or not sure about the availability.
 *          The entity containers which this will be used are:
 *              Nodes, Conditions, Elements
 *
 *        CONSISTENTLY_AVAILABLE: This refers the given variable is available in all of the entities. (Eg. Nodal solution step variable, or nodal data value container variable
 *                                 having all the nodes the same variable)
 *        CONSISTENTLY_UNAVAILABLE: This refers non of the entities has the given variable. (Eg. nodal data value container variable, where non of the nodes has the variable)
 *        INCONCLUSIVE: This refers some of the entities has the value defined, but not all.
 */
enum DataAvailabilityStatesList
{
    CONSISTENTLY_AVAILABLE = 0,
    CONSISTENTLY_UNAVAILABLE  = 1,
    INCONCLUSIVE = 2
};

/**
 * @class FlagIO
 * @brief Utility class for reading and writing flag values on entities.
 *
 * This class provides methods to get, set, and check the presence of flag values
 * on entities. Flags are represented as unsigned char values, where 255 indicates
 * an undefined flag state.
 *
 * @tparam TEntityType The type of entity on which flags are operated.
 *
 * @section TypeDefinitions Type Definitions
 * - TLSType: Alias for unsigned char, used for thread-local storage.
 * - ComponentType: Alias for unsigned char, used for flag components.
 * - DataAvailability: Indicates the data availability state (INCONCLUSIVE).
 *
 * @section PublicOperations Public Operations
 * - GetValue: Returns the value of the specified flag for the given entity.
 * - SetValue: Sets or resets the specified flag for the given entity based on the value.
 * - HasValue: Checks if the specified flag is defined for the given entity.
 */
class FlagIO
{
public:
    ///@name Type definitions
    ///@{

    template<class TDataType>
    using TLSType = unsigned char;

    template<class TDataType>
    using ComponentType = unsigned char;

    constexpr static DataAvailabilityStatesList DataAvailability = DataAvailabilityStatesList::INCONCLUSIVE;

    ///@}
    ///@name Public operations
    ///@{

    template<class TEntityType>
    unsigned char GetValue(
        const TEntityType& rEntity,
        const Flags& rFlag,
        unsigned char& rTLS) const
    {
        return (rEntity.IsDefined(rFlag) ? rEntity.Is(rFlag) : 255);
    }

    template<class TEntityType>
    void SetValue(
        TEntityType& rEntity,
        const Flags& rFlag,
        const unsigned char rValue) const
    {
        if (rValue == 255) {
            rEntity.Reset(rFlag);
        } else {
            rEntity.Set(rFlag, static_cast<bool>(rValue));
        }
    }

    template<class TEntityType>
    bool HasValue(
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

    constexpr static DataAvailabilityStatesList DataAvailability = DataAvailabilityStatesList::CONSISTENTLY_AVAILABLE;

    ///@}
    ///@name Life cycle
    ///@{

    HistoricalIO(const unsigned int StepIndex) : mStepIndex(StepIndex) {}

    ///@}
    ///@name Public operations
    ///@{

    template<class TDataType>
    const TDataType& GetValue(
        const Node& rNode,
        const Variable<TDataType>& rVariable,
        char& rTLS) const
    {
        return rNode.FastGetSolutionStepValue(rVariable, mStepIndex);
    }

    template<class TDataType>
    void SetValue(
        Node& rNode,
        const Variable<TDataType>& rVariable,
        const TDataType& rValue) const
    {
        rNode.FastGetSolutionStepValue(rVariable, mStepIndex) = rValue;
    }

    template<class TDataType>
    bool HasValue(
        const Node& rNode,
        const Variable<TDataType>& rVariable) const
    {
        return rNode.SolutionStepsDataHas(rVariable);
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

    constexpr static DataAvailabilityStatesList DataAvailability = DataAvailabilityStatesList::INCONCLUSIVE;

    ///@}
    ///@name Public operations
    ///@{

    template<class TEntityType, class TDataType>
    const TDataType& GetValue(
        const TEntityType& rEntity,
        const Variable<TDataType>& rVariable,
        char& rTLS) const
    {
        return rEntity.GetValue(rVariable);
    }

    template<class TEntityType, class TDataType>
    void SetValue(
        TEntityType& rEntity,
        const Variable<TDataType>& rVariable,
        const TDataType& rValue) const
    {
        rEntity.SetValue(rVariable, rValue);
    }

    template<class TEntityType, class TDataType>
    bool HasValue(
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

    constexpr static DataAvailabilityStatesList DataAvailability = DataAvailabilityStatesList::CONSISTENTLY_AVAILABLE;

    ///@}
    ///@name Life cycle
    ///@{

    BossakIO(const double AlphaBossak) : mAlphaBossak(AlphaBossak) {}

    ///@}
    ///@name Public operations
    ///@{

    template<class TDataType>
    TDataType GetValue(
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
    void SetValue(
        Node& rNode,
        const Variable<TDataType>& rVariable,
        const TDataType& rValue) const
    {
        rNode.FastGetSolutionStepValue(rVariable) =  rValue;
    }

    template<class TDataType>
    bool HasValue(
        const Node& rNode,
        const Variable<TDataType>& rVariable) const
    {
        return rNode.SolutionStepsDataHas(rVariable);
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

    constexpr static DataAvailabilityStatesList DataAvailability = DataAvailabilityStatesList::CONSISTENTLY_UNAVAILABLE ;

    ///@}
    ///@name Life cycle
    ///@{

    VertexHistoricalValueIO(const unsigned int StepIndex) : mHistoricalIO(StepIndex) {}

    ///@}
    ///@name Public operations
    ///@{

    template<class TDataType>
    TDataType GetValue(
        const Detail::Vertex& rVertex,
        const Variable<TDataType>& rVariable,
        char& rTLS) const
    {
        return rVertex.GetValue(rVariable, mHistoricalIO);
    }

    template<class TDataType>
    void SetValue(
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

    constexpr static DataAvailabilityStatesList DataAvailability = DataAvailabilityStatesList::CONSISTENTLY_UNAVAILABLE ;

    ///@}
    ///@name Public operations
    ///@{

    template<class TDataType>
    TDataType GetValue(
        const Detail::Vertex& rVertex,
        const Variable<TDataType>& rVariable,
        char& rTLS) const
    {
        return rVertex.GetValue(rVariable, mNonHistoricalIO);
    }

    template<class TDataType>
    void SetValue(
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

    constexpr static DataAvailabilityStatesList DataAvailability = DataAvailabilityStatesList::CONSISTENTLY_UNAVAILABLE ;

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
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::PropertiesContainerType>){
        return "PROPERTIES";
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
    const IndexType Size,
    const IndexType Stride,
    const Vector<bool>& rAvailability)
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

    // size check when the rAvailability is used.
    if constexpr(TContainerDataIO::DataAvailability == Internals::DataAvailabilityStatesList::INCONCLUSIVE) {
        KRATOS_ERROR_IF_NOT(rAvailability.size() == rContainer.size())
            << "The container size and the component availability size mismatch [ rContainer.size() = "
            << rContainer.size() << ", availability vector size = " << rAvailability.size() << " ].\n";
    }

    KRATOS_ERROR_IF_NOT(Size == rContainer.size() * Stride)
        << "The contiguous array size mismatch with data in the container [ "
        << "Contiguous array size = " << Size << ", number of entities = "
        << rContainer.size() << ", data stride = " << Stride << " ].";

    typename TContainerDataIO::template TLSType<component_type> tls{};
    IndexPartition<unsigned int>(rContainer.size()).for_each(tls, [&rAvailability, &rContainer, &rComponent, &rContainerDataIO, pBegin, Stride](const auto Index, auto& rTLS) {
        (void)(&rAvailability); // <== suppress unused capture errors (GCC disregards uses where constexpr if conditions do not apply)
        TDataType const* p_value_begin;
        if constexpr(TContainerDataIO::DataAvailability == Internals::DataAvailabilityStatesList::INCONCLUSIVE) {
            auto p_subrange_begin = pBegin + Index * Stride;

            if (rAvailability[Index]) {
                // the component is available in the entity. Then get the value.
                const auto& value = rContainerDataIO.GetValue(*(rContainer.begin() + Index), rComponent, rTLS);
                p_value_begin = value_type_traits::GetContiguousData(value);
                std::copy(p_value_begin, p_value_begin + Stride, p_subrange_begin);
            } else {
                // the component is not available. Put the default values of the primitive data type
                std::fill_n(p_subrange_begin, Stride, TDataType{});
            }

        } else {
            // the component is always available in the entity or it is always computed even if it is not available. Then get the value.
            const auto& value = rContainerDataIO.GetValue(*(rContainer.begin() + Index), rComponent, rTLS);
            p_value_begin = value_type_traits::GetContiguousData(value);
            auto p_subrange_begin = pBegin + Index * Stride;
            std::copy(p_value_begin, p_value_begin + Stride, p_subrange_begin);
        }
    });

    KRATOS_CATCH("");
}

template<class TContainerType, class TComponentType, class TContainerDataIO, class TDataType>
void CopyFromContiguousDataArray(
    TContainerType& rContainer,
    const TComponentType& rComponent,
    const TContainerDataIO& rContainerDataIO,
    TDataType const* pBegin,
    const std::vector<unsigned int>& rShape,
    const Vector<int>& rAvailability)
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

    IndexPartition<unsigned int>(rContainer.size()).for_each(tls_prototype, [&rAvailability, &rContainer, &rComponent, &rContainerDataIO, pBegin, stride](const auto Index, auto& rTLS) {
        // <== suppress unused capture errors (GCC disregards uses where constexpr if conditions do not apply)
        (void)(&rAvailability);
        (void)(&rContainer);
        (void)(&rComponent);
        (void)(&rContainerDataIO);
        (void)(&pBegin);
        (void)(&stride);

        if constexpr(TContainerDataIO::DataAvailability == Internals::DataAvailabilityStatesList::INCONCLUSIVE) {
            if (rAvailability[Index]) {
                TDataType * p_value_begin = value_type_traits::GetContiguousData(rTLS);
                TDataType const * p_subrange_begin = pBegin + Index * stride;
                std::copy(p_subrange_begin, p_subrange_begin + stride, p_value_begin);
                rContainerDataIO.SetValue(*(rContainer.begin() + Index), rComponent, rTLS);
            }
        } else if constexpr(TContainerDataIO::DataAvailability == Internals::DataAvailabilityStatesList::CONSISTENTLY_AVAILABLE) {
            TDataType * p_value_begin = value_type_traits::GetContiguousData(rTLS);
            TDataType const * p_subrange_begin = pBegin + Index * stride;
            std::copy(p_subrange_begin, p_subrange_begin + stride, p_value_begin);
            rContainerDataIO.SetValue(*(rContainer.begin() + Index), rComponent, rTLS);
        } else {
            // do nothing
        }
    });

    KRATOS_CATCH("");
}

} // namespace Internals
} // namespace HDF5
} // namespace Kratos