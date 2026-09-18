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
#include "containers/nd_data.h"
#include "includes/kratos_components.h"
#include "utilities/data_type_traits.h"
#include "utilities/parallel_utilities.h"
#include "processes/find_global_nodal_entity_neighbours_process.h"

// Include base h
#include "tensor_adaptor_utils.h"

namespace Kratos {

void TensorAdaptorUtils::SaveVariablePointer(
    Serializer& rSerializer,
    const std::string& rTag,
    const VariablePointerType& rVariable)
{
    rSerializer.save(rTag + "Index", rVariable.index());
    std::visit([&rSerializer, &rTag](auto pVariable) {
        rSerializer.save(rTag + "Name", pVariable->Name());
    }, rVariable);
}

void TensorAdaptorUtils::LoadVariablePointer(
    Serializer& rSerializer,
    const std::string& rTag,
    VariablePointerType& rVariable)
{
    std::size_t index;
    rSerializer.load(rTag + "Index", index);
    std::string name;
    rSerializer.load(rTag + "Name", name);
    switch (index) {
        case 0: rVariable = &KratosComponents<Variable<double>>::Get(name); break;
        case 1: rVariable = &KratosComponents<Variable<array_1d<double, 3>>>::Get(name); break;
        case 2: rVariable = &KratosComponents<Variable<array_1d<double, 4>>>::Get(name); break;
        case 3: rVariable = &KratosComponents<Variable<array_1d<double, 6>>>::Get(name); break;
        case 4: rVariable = &KratosComponents<Variable<array_1d<double, 9>>>::Get(name); break;
        case 5: rVariable = &KratosComponents<Variable<Vector>>::Get(name); break;
        case 6: rVariable = &KratosComponents<Variable<Matrix>>::Get(name); break;
        default: KRATOS_ERROR << "Unknown tensor adaptor variable variant index: " << index << std::endl;
    }
}

template<class TContainerType>
TensorAdaptor<int>::Pointer TensorAdaptorUtils::CreateNodalNeighboursCountTensorAdaptor(ModelPart& rModelPart)
{
    KRATOS_TRY

    // find the nodal neighbours for the TContainerType
    FindGlobalNodalEntityNeighboursProcess<TContainerType>(rModelPart).Execute();

    auto p_nd_data = Kratos::make_shared<NDData<int>>(DenseVector<unsigned int>(1, rModelPart.NumberOfNodes()));
    auto p_result = Kratos::make_shared<TensorAdaptor<int>>(rModelPart.pNodes(), p_nd_data, false);
    auto data_view = p_nd_data->ViewData();

    IndexPartition<IndexType>(p_nd_data->Size()).for_each([&data_view, &rModelPart](const auto Index) {
        int& value = *(data_view.begin() + Index);
        if constexpr(std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>) {
            value = (rModelPart.NodesBegin() + Index)->GetValue(NEIGHBOUR_CONDITIONS).size();
        } else if constexpr(std::is_same_v<TContainerType, ModelPart::ElementsContainerType>) {
            value = (rModelPart.NodesBegin() + Index)->GetValue(NEIGHBOUR_ELEMENTS).size();
        } else {
            value = 0.0;
            static_assert(sizeof(TContainerType) == 0, "Only supports conditions and element types");
        }
    });

    return p_result;

    KRATOS_CATCH("");
}

// template instantiations
template KRATOS_API(KRATOS_CORE) TensorAdaptor<int>::Pointer TensorAdaptorUtils::CreateNodalNeighboursCountTensorAdaptor<ModelPart::ConditionsContainerType>(ModelPart&);
template KRATOS_API(KRATOS_CORE) TensorAdaptor<int>::Pointer TensorAdaptorUtils::CreateNodalNeighboursCountTensorAdaptor<ModelPart::ElementsContainerType>(ModelPart&);

} // namespace Kratos