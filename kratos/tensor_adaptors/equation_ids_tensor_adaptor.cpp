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

namespace {

template <class TContainerType>
void EquationIdsTensorAdaptorCollectData(
    Kratos::span<int> Span,
    const unsigned int Stride,
    const TContainerType& rContainer,
    const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY

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

    KRATOS_CATCH("");
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

} // namespace EquationIdsTensorAdaptorHelpers

template<class TContainerPointerType>
EquationIdsTensorAdaptor::EquationIdsTensorAdaptor(
    TContainerPointerType pContainer,
    ProcessInfo::Pointer pProcessInfo)
    : mpProcessInfo(pProcessInfo)
{
    this->mpContainer = pContainer;

    const auto& r_getter = [pProcessInfo](auto& rEquationIde, const auto& rEntity){
        rEntity.EquationIdVector(rEquationIde, *pProcessInfo);
    };

    this->mpStorage = Kratos::make_shared<Storage>(
        TensorAdaptorUtils::GetTensorShape<std::vector<IndexType>>(*pContainer, r_getter));
}

EquationIdsTensorAdaptor::EquationIdsTensorAdaptor(
    const TensorAdaptor& rOther,
    ProcessInfo::Pointer pProcessInfo,
    const bool Copy)
    : BaseType(rOther, Copy),
      mpProcessInfo(pProcessInfo)
{
    KRATOS_TRY

    if (!HoldsAlternative<ModelPart::ConditionsContainerType::Pointer,
                          ModelPart::ElementsContainerType::Pointer>::Evaluate(this->GetContainer())) {
        KRATOS_ERROR << "EquationIdsTensorAdaptor can only be used with tensor data having condition or element containers "
                     << "[ tensor adaptor = " << rOther << " ].\n";
    }

    KRATOS_CATCH("");
}

void EquationIdsTensorAdaptor::CollectData()
{
    KRATOS_TRY

    std::visit([this](auto pContainer) {
        using container_type = BareType<decltype(*pContainer)>;

        if constexpr(IsInList<container_type, ModelPart::ConditionsContainerType, ModelPart::ElementsContainerType>) {
            KRATOS_ERROR_IF_NOT(this->Shape()[0] == pContainer->size())
                << "First dimension mismatch [ Container size = " << pContainer->size()
                << ", shape = " << this->Shape() << " ].\n";

            EquationIdsTensorAdaptorCollectData(this->ViewData(), this->Shape()[1], *pContainer, *(this->mpProcessInfo));
        }
    }, this->GetContainer());

    KRATOS_CATCH("");
}

void EquationIdsTensorAdaptor::StoreData()
{
    KRATOS_ERROR << "Equation ids storing is not allowed.";
}

std::string EquationIdsTensorAdaptor::Info() const
{
    std::stringstream info;
    info << "EquationIdsTensorAdaptor:";
    info << " " << BaseType::Info();
    return info.str();
}

// template instantiations
template KRATOS_API(KRATOS_CORE) EquationIdsTensorAdaptor::EquationIdsTensorAdaptor(ModelPart::ConditionsContainerType::Pointer, ProcessInfo::Pointer);
template KRATOS_API(KRATOS_CORE) EquationIdsTensorAdaptor::EquationIdsTensorAdaptor(ModelPart::ElementsContainerType::Pointer, ProcessInfo::Pointer);

} // namespace Kratos