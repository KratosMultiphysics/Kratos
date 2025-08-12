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

// Include base h
#include "variable_tensor_adaptor.h"

namespace Kratos {

template<class TContainerPointerType>
VariableTensorAdaptor::VariableTensorAdaptor(
    TContainerPointerType pContainer,
    VariablePointerType pVariable)
    : mpVariable(pVariable)
{
    std::visit([this, pContainer](auto pVariable) {
        this->mpStorage = Kratos::make_intrusive<Storage>(
                                pContainer, TensorAdaptorUtils::GetTensorShape(
                                                *pContainer, *pVariable,
                                                [pVariable](auto& rValue, const auto& rEntity) {
                                                    rValue = rEntity.GetValue(*pVariable);
                                                }));
    }, mpVariable);
}

template<class TContainerPointerType>
VariableTensorAdaptor::VariableTensorAdaptor(
    TContainerPointerType pContainer,
    VariablePointerType pVariable,
    const std::vector<unsigned int>& rDataShape)
    : mpVariable(pVariable)
{
    std::visit([&rDataShape, this, pContainer](auto pVariable) {
        this->mpStorage = Kratos::make_intrusive<Storage>(
                                pContainer, TensorAdaptorUtils::GetTensorShape(
                                *pContainer, *pVariable, rDataShape.data(),
                                rDataShape.data() + rDataShape.size()));
    }, mpVariable);
}

VariableTensorAdaptor::VariableTensorAdaptor(
    const TensorAdaptor& rOther,
    VariablePointerType pVariable,
    const bool Copy)
    : BaseType(rOther, Copy),
      mpVariable(pVariable)
{
    KRATOS_TRY

    if (!HoldsAlternative<ModelPart::NodesContainerType::Pointer,
                          ModelPart::ConditionsContainerType::Pointer,
                          ModelPart::ElementsContainerType::Pointer,
                          ModelPart::PropertiesContainerType::Pointer,
                          ModelPart::MasterSlaveConstraintContainerType::Pointer,
                          ModelPart::GeometryContainerType::Pointer>::Evaluate(this->GetContainer())) {
        KRATOS_ERROR << "VariableTensorAdaptor can only be used with tensor data having nodal, condition, element, properties, master-slave or geometry containers "
                     << "[ tensor adaptor = " << rOther << " ].\n";
    }

    // now check whether the given storage is compatible with the variable.
    std::visit([this, &rOther](auto pVariable) {
        using data_type = typename BareType<decltype(*pVariable)>::Type;
        const auto& r_data_shape = this->mpStorage->DataShape();
        KRATOS_ERROR_IF_NOT(DataTypeTraits<data_type>::IsValidShape(r_data_shape.data().begin(), r_data_shape.data().begin() + r_data_shape.size()))
            << "The data storage within the tensor data is not compatible with the " << pVariable->Name()
            << " [ origin data shape = " << rOther.DataShape() << ", tensor adaptor = " << *this << " ].\n";

    }, mpVariable);

    KRATOS_CATCH("");
}

void VariableTensorAdaptor::Check() const
{
    KRATOS_TRY

    BaseType::Check();

    std::visit([](auto pContainer, auto pVariable) {
        using container_type = BareType<decltype(*pContainer)>;

        if constexpr(IsInList<container_type,
                              ModelPart::NodesContainerType,
                              ModelPart::ConditionsContainerType,
                              ModelPart::ElementsContainerType,
                              ModelPart::PropertiesContainerType,
                              ModelPart::MasterSlaveConstraintContainerType,
                              ModelPart::GeometryContainerType>) {

            block_for_each(*pContainer, [&pVariable](const auto& rEntity) {
                KRATOS_ERROR_IF_NOT(rEntity.Has(*pVariable))
                    << "The entity with id = " << rEntity.Id() << " does not have the variable " << pVariable->Name() << ".\n";
            });
        }
    }, this->mpStorage->GetContainer(), mpVariable);

    KRATOS_CATCH("");
}

void VariableTensorAdaptor::CollectData()
{
    std::visit([this](auto pContainer, auto pVariable) {
        using container_type = BareType<decltype(*pContainer)>;

        if constexpr(IsInList<container_type,
                              ModelPart::NodesContainerType,
                              ModelPart::ConditionsContainerType,
                              ModelPart::ElementsContainerType,
                              ModelPart::PropertiesContainerType,
                              ModelPart::MasterSlaveConstraintContainerType,
                              ModelPart::GeometryContainerType>) {

            using variable_type = BareType<decltype(*pVariable)>;
            using data_type = typename variable_type::Type;

            const auto& r_tensor_shape = this->Shape();

            KRATOS_ERROR_IF_NOT(r_tensor_shape[0] == pContainer->size())
                << "Underlying container of the tensor data has changed size [ tensor data = "
                << this->mpStorage->Info() << ", container size = " << pContainer->size() << " ].\n";

            ContainerIOUtils::CopyToContiguousArray<data_type>(
                *pContainer, this->ViewData(), r_tensor_shape.data().begin(),
                r_tensor_shape.data().begin() + r_tensor_shape.size(),
                [pVariable](auto& rValue, const auto& rEntity) {
                    rValue = rEntity.GetValue(*pVariable);
                });
        }
    }, this->mpStorage->GetContainer(), mpVariable);
}

void VariableTensorAdaptor::StoreData()
{
    std::visit([this](auto pContainer, auto pVariable) {
        using container_type = BareType<decltype(*pContainer)>;

        if constexpr(IsInList<container_type,
                              ModelPart::NodesContainerType,
                              ModelPart::ConditionsContainerType,
                              ModelPart::ElementsContainerType,
                              ModelPart::PropertiesContainerType,
                              ModelPart::MasterSlaveConstraintContainerType,
                              ModelPart::GeometryContainerType>) {

            using variable_type = BareType<decltype(*pVariable)>;
            using data_type = typename variable_type::Type;

            const auto& r_tensor_shape = this->Shape();

            KRATOS_ERROR_IF_NOT(r_tensor_shape[0] == pContainer->size())
                << "Underlying container of the tensor data has changed size [ tensor data = "
                << this->mpStorage->Info() << ", container size = " << pContainer->size() << " ].\n";

            if constexpr(DataTypeTraits<data_type>::IsDynamic) {
                // this zero value may be different from the Variable::Zero()
                // in the case where the variable type is Variable<Vector> or Variable<Matrix>.
                // Because, in these dynamic data type variables, Variable::Zero will create
                // a zero sized vector or matrix, which is useless in the StoreData method.
                // The following method creates correctly sized zero Vector or Matrix
                // accordingly to be assigned to the entities, if they don't have the specified
                // variable in their DataValueContainer.
                const auto& zero = TensorAdaptorUtils::GetZeroValue(*pVariable, this->DataShape());

                // -----------------------------------------------------------------------------------------
                // TODO: Following block of code is not required if the discussion in PR #13685 is accepted.
                //       This is costly, which can be heavily improved with the above mentioned PR.
                block_for_each(*pContainer, [&pVariable, &zero](auto& rEntity) {
                    // following code will need to do two lookups, one for Has, and one for SetValue.
                    // The PR #13685 have a mechanism to do both with one lookup, halving the cost.
                    if (!rEntity.Has(*pVariable)) {
                        // this variable only needs to be set, if it is not present because,
                        // in a case where this tensor is created in the following scenario
                        //
                        // a = VariableTensorAdaptor(nodes, VELOCITY, data_shape = [2]) //
                        // a.data = numpy.ndarray
                        // a.StoreData() # here i only want to write to VELOCITY_X and VELOCITY_Y, not touching VELOCITY_Z
                        // in order to do that, i must only set the values in the entities which are not having the variable.
                        // This can be easily avoided if as suggested in the PR #13685 we have the Emplace method.
                        rEntity.SetValue(*pVariable, zero);
                    }
                });
                // -----------------------------------------------------------------------------------------

                ContainerIOUtils::CopyFromContiguousDataArray<data_type>(
                    *pContainer, this->ViewData(), r_tensor_shape.data().begin(),
                    r_tensor_shape.data().begin() + r_tensor_shape.size(),
                    [pVariable](auto& rEntity) -> auto& {
                        // -------------------------------------------------------------
                        // with the PR #13685, we can reduce the 3 lookups per entity
                        // to 1 lookup.
                        return rEntity.GetValue(*pVariable);
                        // -------------------------------------------------------------
                        // TODO: Needs the approval of PR #13685
                        // only dynamic data_type types require the zero
                        // to initialize the uninitialized variables, because
                        // they need to be correctly sized.
                        //      return rEntity.Emplace(*pVariable, zero);
                        // -------------------------------------------------------------
                    });
            } else {
                // -----------------------------------------------------------------------------------------
                // TODO: Following block of code is not required if the discussion in PR #13685 is accepted.
                //       This is costly, which can be heavily improved with the above mentioned PR.
                block_for_each(*pContainer, [&pVariable](auto& rEntity) {
                    // following code will need to do two lookups, one for Has, and one for SetValue.
                    // The PR #13685 have a mechanism to do both with one lookup, halving the cost.
                    if (!rEntity.Has(*pVariable)) {
                        // this variable only needs to be set, if it is not present because,
                        // in a case where this tensor is created in the following scenario
                        //
                        // a = VariableTensorAdaptor(nodes, VELOCITY, data_shape = [2]) //
                        // a.data = numpy.ndarray
                        // a.StoreData() # here i only want to write to VELOCITY_X and VELOCITY_Y, not touching VELOCITY_Z
                        // in order to do that, i must only set the values in the entities which are not having the variable.
                        // This can be easily avoided if as suggested in the PR #13685 we have the Emplace method.
                        rEntity.SetValue(*pVariable, pVariable->Zero());
                    }
                });
                // -----------------------------------------------------------------------------------------

                ContainerIOUtils::CopyFromContiguousDataArray<data_type>(
                    *pContainer, this->ViewData(), r_tensor_shape.data().begin(),
                    r_tensor_shape.data().begin() + r_tensor_shape.size(),
                    [pVariable](auto& rEntity) -> auto& {
                        // -------------------------------------------------------------
                        // with the PR #13685, we can reduce the 3 lookups per entity
                        // to 1 lookup.
                        return rEntity.GetValue(*pVariable);
                        // -------------------------------------------------------------
                        // TODO: Needs the approval of PR #13685
                        // only dynamic data_type types require the zero
                        // to initialize the uninitialized variables, because
                        // they need to be correctly sized.
                        //      return rEntity.Emplace(*pVariable);
                        // -------------------------------------------------------------
                    });
            }
        }
    }, this->mpStorage->GetContainer(), mpVariable);
}

std::string VariableTensorAdaptor::Info() const
{
    std::stringstream info;
    info << "VariableTensorAdaptor:";
    std::visit([&info](auto pVariable) {
        info << " Variable = " << pVariable->Name();
    }, this->mpVariable);
    info << ", " << this->mpStorage->Info();
    return info.str();
}

template KRATOS_API(KRATOS_CORE) VariableTensorAdaptor::VariableTensorAdaptor(ModelPart::NodesContainerType::Pointer, VariablePointerType);
template KRATOS_API(KRATOS_CORE) VariableTensorAdaptor::VariableTensorAdaptor(ModelPart::NodesContainerType::Pointer, VariablePointerType, const std::vector<unsigned int>&);
template KRATOS_API(KRATOS_CORE) VariableTensorAdaptor::VariableTensorAdaptor(ModelPart::ConditionsContainerType::Pointer, VariablePointerType);
template KRATOS_API(KRATOS_CORE) VariableTensorAdaptor::VariableTensorAdaptor(ModelPart::ConditionsContainerType::Pointer, VariablePointerType, const std::vector<unsigned int>&);
template KRATOS_API(KRATOS_CORE) VariableTensorAdaptor::VariableTensorAdaptor(ModelPart::ElementsContainerType::Pointer, VariablePointerType);
template KRATOS_API(KRATOS_CORE) VariableTensorAdaptor::VariableTensorAdaptor(ModelPart::ElementsContainerType::Pointer, VariablePointerType, const std::vector<unsigned int>&);
template KRATOS_API(KRATOS_CORE) VariableTensorAdaptor::VariableTensorAdaptor(ModelPart::PropertiesContainerType::Pointer, VariablePointerType);
template KRATOS_API(KRATOS_CORE) VariableTensorAdaptor::VariableTensorAdaptor(ModelPart::PropertiesContainerType::Pointer, VariablePointerType, const std::vector<unsigned int>&);
template KRATOS_API(KRATOS_CORE) VariableTensorAdaptor::VariableTensorAdaptor(ModelPart::GeometryContainerType::Pointer, VariablePointerType);
template KRATOS_API(KRATOS_CORE) VariableTensorAdaptor::VariableTensorAdaptor(ModelPart::GeometryContainerType::Pointer, VariablePointerType, const std::vector<unsigned int>&);
template KRATOS_API(KRATOS_CORE) VariableTensorAdaptor::VariableTensorAdaptor(ModelPart::MasterSlaveConstraintContainerType::Pointer, VariablePointerType);
template KRATOS_API(KRATOS_CORE) VariableTensorAdaptor::VariableTensorAdaptor(ModelPart::MasterSlaveConstraintContainerType::Pointer, VariablePointerType, const std::vector<unsigned int>&);

} // namespace Kratos