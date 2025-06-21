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

// System includes
#include <variant>
#include <sstream>

// External includes

// Project includes
#include "includes/define.h"
#include "containers/variable.h"
#include "includes/model_part.h"
#include "utilities/data_type_traits.h"
#include "utilities/container_io_utils.h"
#include "utilities/parallel_utilities.h"

// Include base h
#include "entity_variable_tensor_adaptor.h"

namespace Kratos {

namespace Detail {

template <class TContainerType, class TContainerIOType>
struct DoubleDataValueContainerIO
{
    using PrimitiveDataType = double;

    using ContainerType = TContainerType;

    using VariableType = std::variant<
                                const Variable<double>*,
                                const Variable<array_1d<double, 3>>*,
                                const Variable<array_1d<double, 4>>*,
                                const Variable<array_1d<double, 6>>*,
                                const Variable<array_1d<double, 9>>*,
                                const Variable<Vector>*,
                                const Variable<Matrix>*
                                >;

    using ContainerIOType = TContainerIOType;

    static std::string TypeInfo()
    {
        if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
            return "NodeNonHistoricalVariableTensorAdaptor";
        } else if constexpr(std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>) {
            return "ConditionVariableTensorAdaptor";
        } else if constexpr(std::is_same_v<TContainerType, ModelPart::ElementsContainerType>) {
            return "ElementVariableTensorAdaptor";
        } else if constexpr(std::is_same_v<TContainerType, ModelPart::GeometryContainerType::GeometriesMapType>) {
            return "GeometryVariableTensorAdaptor";
        } else if constexpr(std::is_same_v<TContainerType, ModelPart::PropertiesContainerType>) {
            return "PropertyVariableTensorAdaptor";
        } else {
            static_assert(std::is_same_v<TContainerType, TContainerType>, "Unsupported type");
            return "";
        }
    }
};
} // namespace Detail

template <class TIOType, class... TArgs>
EntityVariableTensorAdaptor<TIOType, TArgs...>::EntityVariableTensorAdaptor(
    typename ContainerType::Pointer pContainer,
    typename TIOType::VariableType pVariable,
    TArgs&&... rArgs)
    : BaseType(pContainer),
      mpVariable(pVariable),
      mContainerIO(rArgs...)
{
    // fix the shape correctly
    std::visit([this](const auto pVariable) {
        using variable_data_type = typename std::remove_cv_t<std::decay_t<decltype(*pVariable)>>::Type;
        using type_trait = DataTypeTraits<variable_data_type>;

        if constexpr(!type_trait::IsDynamic) {
            // these are for the static type such as double, array_1d<double, 3>, ...
            const auto& r_shape = type_trait::Shape(variable_data_type{});
            this->mShape.resize(r_shape.size() + 1);
            std::copy(r_shape.begin(), r_shape.end(), this->mShape.data().begin() + 1);
        } else {
            // these are for the dynamic types such as Vector, Matrix
            if (this->mpContainer->empty()) {
                this->mShape.resize(1);
            } else {
                typename TIOType::ContainerIOType::TLSType<variable_data_type> tls;
                const auto& value = this->mContainerIO.GetValue(*(this->mpContainer->begin()), *pVariable, tls);
                const auto& r_shape = type_trait::Shape(value);
                this->mShape.resize(r_shape.size() + 1);
                std::copy(r_shape.begin(), r_shape.end(), this->mShape.data().begin() + 1);
            }
        }

        this->mShape[0] = this->mpContainer->size();
    }, mpVariable);

    // now resize
    this->mData.resize(std::accumulate(this->mShape.begin(), this->mShape.end(), 1, std::multiplies<int>{}));
}

template <class TIOType, class... TArgs>
void EntityVariableTensorAdaptor<TIOType, TArgs...>::CollectData()
{
    KRATOS_TRY

    // sanity checks
    KRATOS_ERROR_IF_NOT(this->mShape[0] == static_cast<int>(this->mpContainer->size()))
        << "First dimension of the initialized tensor adaptor mismatch with the container size [ "
        << "Tensor adapter shape = " << this->mShape << ", container size = " << this->mpContainer->size() << " ].\n";

    std::visit([this](auto pVariable) {
        CopyToContiguousArray(*(this->mpContainer), *pVariable, this->mContainerIO, this->mData.data().begin(), this->mData.size());
    }, this->mpVariable);

    KRATOS_CATCH("");
}

template <class TIOType, class... TArgs>
void EntityVariableTensorAdaptor<TIOType, TArgs...>::StoreData()
{
    KRATOS_TRY

    // sanity checks
    KRATOS_ERROR_IF_NOT(this->mShape[0] == static_cast<int>(this->mpContainer->size()))
        << "First dimension of the initialized tensor adaptor mismatch with the container size [ "
        << "Tensor adapter shape = " << this->mShape << ", container size = " << this->mpContainer->size() << " ].\n";

    std::visit([this](auto pVariable) {
        std::vector<unsigned int> shape;
        shape.resize(this->mShape.size() - 1);
        std::copy(this->mShape.begin(), this->mShape.end(), shape.begin());
        CopyFromContiguousDataArray(*(this->mpContainer), *pVariable, this->mContainerIO, this->mData.data().begin(), shape);
    }, this->mpVariable);

    KRATOS_CATCH("");
}

template <class TIOType, class... TArgs>
std::string EntityVariableTensorAdaptor<TIOType, TArgs...>::Info() const
{
    return std::visit([this](auto pVariable) {
        std::stringstream msg;
        msg << TIOType::TypeInfo() << " - " << pVariable->Name() << " with " << this->mpContainer->size() << " entities.";
        return msg.str();
    }, mpVariable);
}

// template instantiations

// instantiations for the data value container
template class KRATOS_API(KRATOS_CORE) EntityVariableTensorAdaptor<Detail::DoubleDataValueContainerIO<ModelPart::NodesContainerType, NonHistoricalIO>>;
template class KRATOS_API(KRATOS_CORE) EntityVariableTensorAdaptor<Detail::DoubleDataValueContainerIO<ModelPart::ConditionsContainerType, NonHistoricalIO>>;
template class KRATOS_API(KRATOS_CORE) EntityVariableTensorAdaptor<Detail::DoubleDataValueContainerIO<ModelPart::ElementsContainerType, NonHistoricalIO>>;
template class KRATOS_API(KRATOS_CORE) EntityVariableTensorAdaptor<Detail::DoubleDataValueContainerIO<ModelPart::PropertiesContainerType, NonHistoricalIO>>;
template class KRATOS_API(KRATOS_CORE) EntityVariableTensorAdaptor<Detail::DoubleDataValueContainerIO<ModelPart::GeometryContainerType::GeometriesMapType, NonHistoricalIO>>;

// instantiations for the historical data value container
template class KRATOS_API(KRATOS_CORE) EntityVariableTensorAdaptor<Detail::DoubleDataValueContainerIO<ModelPart::NodesContainerType, HistoricalIO>, const unsigned int>;


} // namespace Kratos