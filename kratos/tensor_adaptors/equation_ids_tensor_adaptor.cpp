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

// External includes

// Project includes
#include "utilities/data_type_traits.h"
#include "utilities/parallel_utilities.h"
#include "tensor_adaptor_utils.h"

// Include base h
#include "equation_ids_tensor_adaptor.h"

namespace Kratos {

namespace EquationIdsTensorAdaptorHelpers {

template <class TContainerType>
void CollectData(
    Kratos::span<int> Span,
    const unsigned int Stride,
    const TContainerType& rContainer,
    const ProcessInfo& rProcessInfo)
{
    if constexpr(IsInList<TContainerType, ModelPart::ConditionsContainerType, ModelPart::ElementsContainerType>::value) {
        IndexPartition<IndexType>(rContainer.size()).for_each(std::vector<IndexType>{}, [Span, &rContainer, &rProcessInfo, Stride](const auto Index, auto& rTLS){
            const auto& r_entity = *(rContainer.begin() + Index);
            r_entity.EquationIdVector(rTLS, rProcessInfo);

            KRATOS_DEBUG_ERROR_IF_NOT(rTLS.size() == Stride)
                << "Non-uniform equation id vectors are found in container having " << rContainer.size()
                << " " << ModelPart::Container<TContainerType>::GetEntityName() << "(s)"
                << " [ Required size of equation ids = " << Stride << ", found EquationIds = "
                << rTLS << " ].\n";

            std::copy(rTLS.begin(), rTLS.end(), Span.data() + Index * Stride);
        });
    } else {
        KRATOS_ERROR << "EquationIdsTensorAdaptor only works with ModelPart::ConditionsContainerType or ModelPart::ElementsContainerType.\n";
    }
}

template<class TContainerType>
std::string InfoImpl(
    const DenseVector<unsigned int>& rShape,
    const TContainerType& rContainer)
{
    std::stringstream msg;
    msg << "EquationIdsTensorAdaptor: " << " with " << rContainer.size()
        << " " << ModelPart::Container<TContainerType>::GetEntityName() << "(s) having shape "
        << rShape << " ].";
    return msg.str();
}

template<class TContainerPointerType, class TSpanType>
TensorAdaptor<int>::Pointer Clone(
    TContainerPointerType pContainer,
    const ProcessInfo::Pointer pProcessInfo,
    const TSpanType& rDataSpan)
{
    if constexpr(IsInList<std::remove_cv_t<std::decay_t<decltype(*pContainer)>>, ModelPart::ConditionsContainerType, ModelPart::ElementsContainerType>::value) {
        auto p_tensor_adaptor = Kratos::make_intrusive<EquationIdsTensorAdaptor>(pContainer, pProcessInfo);
        IndexPartition<IndexType>(p_tensor_adaptor->Size()).for_each([p_tensor_adaptor, &rDataSpan](const auto Index) {
            p_tensor_adaptor->ViewData()[Index] = rDataSpan[Index];
        });
        return p_tensor_adaptor;
    } else {
        KRATOS_ERROR << "EquationIdsTensorAdaptor only works with ModelPart::ConditionsContainerType or ModelPart::ElementsContainerType.\n";
        return nullptr;
    }
}

} // namespace EquationIdsTensorAdaptorHelpers

template<class TContainerPointerType>
EquationIdsTensorAdaptor::EquationIdsTensorAdaptor(
    TContainerPointerType pContainer,
    ProcessInfo::Pointer pProcessInfo)
    : mpProcessInfo(pProcessInfo),
      mpContainer(pContainer)
{
    const auto& r_getter = [pProcessInfo](auto& rEquationIde, const auto& rEntity){
        rEntity.EquationIdVector(rEquationIde, *pProcessInfo);
    };

    this->SetShape(TensorAdaptorUtils::GetTensorShape<std::vector<IndexType>>(*pContainer, r_getter));
}

EquationIdsTensorAdaptor::BaseType::Pointer EquationIdsTensorAdaptor::Clone() const
{
    return std::visit([this](auto pContainer) {
        return EquationIdsTensorAdaptorHelpers::Clone(pContainer, this->mpProcessInfo, this->ViewData());
    }, mpContainer);
}

void EquationIdsTensorAdaptor::CollectData()
{
    std::visit([this](auto pContainer) {
        EquationIdsTensorAdaptorHelpers::CollectData(this->ViewData(), this->Shape()[1], *pContainer, *(this->mpProcessInfo));
    }, mpContainer);
}

void EquationIdsTensorAdaptor::StoreData()
{
    KRATOS_ERROR << "Equation ids storing is not allowed.";
}

EquationIdsTensorAdaptor::ContainerPointerType EquationIdsTensorAdaptor::GetContainer() const
{
    return mpContainer;
}

std::string EquationIdsTensorAdaptor::Info() const
{
    return std::visit([this](auto pContainer) {
        return TensorAdaptorUtils::Info("EquationIdsTensorAdaptor ", this->Shape(), *pContainer);
    }, mpContainer);
}

// template instantiations
template KRATOS_API(KRATOS_CORE) EquationIdsTensorAdaptor::EquationIdsTensorAdaptor(ModelPart::ConditionsContainerType::Pointer, ProcessInfo::Pointer);
template KRATOS_API(KRATOS_CORE) EquationIdsTensorAdaptor::EquationIdsTensorAdaptor(ModelPart::ElementsContainerType::Pointer, ProcessInfo::Pointer);

} // namespace Kratos