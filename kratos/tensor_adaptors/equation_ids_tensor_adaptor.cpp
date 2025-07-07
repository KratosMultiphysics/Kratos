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

// Include base h
#include "equation_ids_tensor_adaptor.h"

namespace Kratos {

namespace EquationIdsTensorAdaptorHelpers {

template <class TContainerType>
void CollectDataImpl(
    DenseVector<int>& rData,
    const unsigned int Stride,
    const TContainerType& rContainer,
    const ProcessInfo& rProcessInfo)
{
    if constexpr(IsInList<TContainerType, ModelPart::ConditionsContainerType, ModelPart::ElementsContainerType>::value) {
        IndexPartition<IndexType>(rContainer.size()).for_each(std::vector<IndexType>{}, [&rData, &rContainer, &rProcessInfo, Stride](const auto Index, auto& rTLS){
            const auto& r_entity = *(rContainer.begin() + Index);
            r_entity.EquationIdVector(rTLS, rProcessInfo);

            KRATOS_ERROR_IF_NOT(rTLS.size() == Stride)
                << "Non-uniform equation id vectors are found in container having " << rContainer.size()
                << " " << ModelPart::Container<TContainerType>::GetEntityName() << "(s)"
                << " [ Required size of equation ids = " << Stride << ", found EquationIds = "
                << rTLS << " ].\n";

            std::copy(rTLS.begin(), rTLS.end(), rData.data().begin() + Index * Stride);
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

} // namespace EquationIdsTensorAdaptorHelpers

template<class TContainerPointerType>
EquationIdsTensorAdaptor::EquationIdsTensorAdaptor(
    TContainerPointerType pContainer,
    ProcessInfo::Pointer pProcessInfo)
    : mpProcessInfo(pProcessInfo),
      mpContainer(pContainer)
{
    KRATOS_ERROR_IF(pContainer->empty())
        << "Cannot create an EquationIdsTensorAdaptor with an empty container.\n";

    std::vector<IndexType> equation_ids;
    pContainer->front().EquationIdVector(equation_ids, *mpProcessInfo);

    // setting the tensor shape
    this->mShape.resize(2);
    this->mShape[0] = pContainer->size();  // first dimension is number of entities
    this->mShape[1] = equation_ids.size(); // second dimension is number of equation ids per entity

    // correctly size the underlying container
    this->mData.resize(this->Size(), false);
}


void EquationIdsTensorAdaptor::CollectData()
{
    std::visit([this](auto pContainer) {
        EquationIdsTensorAdaptorHelpers::CollectDataImpl(this->mData, this->mShape[1], *pContainer, *(this->mpProcessInfo));
    }, mpContainer);
}

void EquationIdsTensorAdaptor::StoreData()
{
    KRATOS_ERROR << "Equation ids storing is not allowed.";
}

EquationIdsTensorAdaptor::ContainerType EquationIdsTensorAdaptor::GetContainer() const
{
    return mpContainer;
}

std::string EquationIdsTensorAdaptor::Info() const
{
    return std::visit([this](auto pContainer) {
        return EquationIdsTensorAdaptorHelpers::InfoImpl(this->mShape, *pContainer);
    }, mpContainer);
}

// template instantiations
template EquationIdsTensorAdaptor::EquationIdsTensorAdaptor(ModelPart::ConditionsContainerType::Pointer, ProcessInfo::Pointer);
template EquationIdsTensorAdaptor::EquationIdsTensorAdaptor(ModelPart::ElementsContainerType::Pointer, ProcessInfo::Pointer);

} // namespace Kratos