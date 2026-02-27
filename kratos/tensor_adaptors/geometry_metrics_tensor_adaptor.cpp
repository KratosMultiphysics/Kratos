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
#include <unordered_map>

// External includes

// Project includes
#include "utilities/atomic_utilities.h"
#include "utilities/data_type_traits.h"

// Include base h
#include "geometry_metrics_tensor_adaptor.h"

namespace Kratos {

namespace GeometryMetricsTensorAdaptorHelperUtils {

DenseVector<unsigned int> GetShape(
    const unsigned int NumberOfEntities,
    GeometryMetricsTensorAdaptor::DatumType Datum)
{
    switch (Datum) {
        case GeometryMetricsTensorAdaptor::DatumType::DomainSize:
            return DenseVector<unsigned int>(1, NumberOfEntities);
    }

    return DenseVector<unsigned int>(0);
}

template<class TEntityType>
const ModelPart::GeometryType& GetGeometry(const TEntityType& rEntity)
{
    if constexpr(std::is_same_v<TEntityType, ModelPart::GeometryType>) {
        return rEntity;
    } else {
        return rEntity.GetGeometry();
    }
}

template<class TContainerType>
void FillDomainSize(
    Kratos::span<double> DataSpan,
    const TContainerType& rContainer)
{
    IndexPartition<IndexType>(rContainer.size()).for_each([&rContainer, &DataSpan](const auto Index) {
        DataSpan[Index] = GetGeometry(*(rContainer.begin() + Index)).DomainSize();
    });
}

template<class TContainerType>
void FillData(
    Kratos::span<double> DataSpan,
    const GeometryMetricsTensorAdaptor::DatumType Datum,
    const TContainerType& rContainer)
{
    switch (Datum) {
        case GeometryMetricsTensorAdaptor::DatumType::DomainSize:
            FillDomainSize(DataSpan, rContainer);
            break;
    }
}

} // namespace GeometryMetricsTensorAdaptorHelperUtils

template<class TContainerPointerType>
GeometryMetricsTensorAdaptor::GeometryMetricsTensorAdaptor(
    TContainerPointerType pContainer,
    const DatumType Datum)
    : mDatum(Datum)
{
    this->mpContainer = pContainer;
    this->mpStorage = Kratos::make_shared<Storage>(DenseVector<unsigned int>(1, pContainer->size()));
}

GeometryMetricsTensorAdaptor::GeometryMetricsTensorAdaptor(
    const TensorAdaptor& rOther,
    const DatumType Datum,
    const bool Copy)
    : BaseType(rOther, Copy),
      mDatum(Datum)
{
    KRATOS_TRY

    if (!HoldsAlternative<ModelPart::GeometryContainerType::Pointer, ModelPart::ConditionsContainerType::Pointer,
                          ModelPart::ElementsContainerType::Pointer>::Evaluate(this->GetContainer())) {
        KRATOS_ERROR << "GeometryMetricsTensorAdaptor can only be used with tensor data "
                        "having geometry, condition or element containers.";
    }

    const auto& current_shape = this->Shape();

    KRATOS_ERROR_IF(current_shape.size() < 1) << "Tensor data's first dimension should represent number of entities [ tensor adaptor = " << *this << " ].\n";

    const auto& required_shape = GeometryMetricsTensorAdaptorHelperUtils::GetShape(current_shape[0], Datum);

    for (IndexType i = 0; i < required_shape.size(); ++i) {
        KRATOS_ERROR_IF_NOT(current_shape[i] == required_shape[i])
            << "Incompatible tensor data shape [ required tensor data shape = " << required_shape
            << ", current tensor data shape = " << current_shape << ", tensor adaptor = " << *this << " ].\n";
    }

    KRATOS_CATCH("");
}

TensorAdaptor<double>::Pointer GeometryMetricsTensorAdaptor::Clone() const
{
    return Kratos::make_shared<GeometryMetricsTensorAdaptor>(*this);
}

void GeometryMetricsTensorAdaptor::CollectData()
{
    KRATOS_TRY

    std::visit([this](auto pContainer) {
        using container_type = BareType<decltype(*pContainer)>;
        if constexpr(IsInList<container_type, ModelPart::GeometryContainerType, ModelPart::ConditionsContainerType, ModelPart::ElementsContainerType>) {
            GeometryMetricsTensorAdaptorHelperUtils::FillData(this->ViewData(), this->mDatum, *pContainer);
        }

    }, this->GetContainer());

    KRATOS_CATCH("");
}

void GeometryMetricsTensorAdaptor::StoreData()
{
    KRATOS_ERROR << "StoreData method is not implemented.";
}

std::string GeometryMetricsTensorAdaptor::Info() const
{
    std::stringstream info;
    // info << "GeometryMetricsTensorAdaptor:";
    // std::visit([&info](auto pContainer) {
    //     info << " number of  " << ModelPart::Container<BareType<decltype(*pContainer)>>::GetEntityName() << " = " << pContainer->size();
    // }, this->mpEntityContainer);
    // info << ", " << BaseType::Info();
    return info.str();
}

template KRATOS_API(KRATOS_CORE) GeometryMetricsTensorAdaptor::GeometryMetricsTensorAdaptor(ModelPart::GeometryContainerType::Pointer, GeometryMetricsTensorAdaptor::DatumType);
template KRATOS_API(KRATOS_CORE) GeometryMetricsTensorAdaptor::GeometryMetricsTensorAdaptor(ModelPart::ConditionsContainerType::Pointer, GeometryMetricsTensorAdaptor::DatumType);
template KRATOS_API(KRATOS_CORE) GeometryMetricsTensorAdaptor::GeometryMetricsTensorAdaptor(ModelPart::ElementsContainerType::Pointer, GeometryMetricsTensorAdaptor::DatumType);

} // namespace Kratos