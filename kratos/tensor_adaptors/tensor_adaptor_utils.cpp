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
#include "utilities/data_type_traits.h"
#include "utilities/parallel_utilities.h"
#include "processes/find_global_nodal_entity_neighbours_process.h"

// Include base h
#include "tensor_adaptor_utils.h"

namespace Kratos {


template<class TContainerType>
TensorAdaptor<int>::Pointer TensorAdaptorUtils::GetNodalNeighboursCountTensorAdaptor(ModelPart& rModelPart)
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
            static_assert(sizeof(TContainerType) == 0, "Only supports conditions and element types");
        }
    });

    return p_result;

    KRATOS_CATCH("");
}

// template instantiations
template KRATOS_API(KRATOS_CORE) TensorAdaptor<int>::Pointer TensorAdaptorUtils::GetNodalNeighboursCountTensorAdaptor<ModelPart::ConditionsContainerType>(ModelPart&);
template KRATOS_API(KRATOS_CORE) TensorAdaptor<int>::Pointer TensorAdaptorUtils::GetNodalNeighboursCountTensorAdaptor<ModelPart::ElementsContainerType>(ModelPart&);

} // namespace Kratos