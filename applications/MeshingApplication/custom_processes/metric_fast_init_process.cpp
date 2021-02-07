// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                       license: MeshingApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "utilities/parallel_utilities.h"
#include "custom_processes/metric_fast_init_process.h"

namespace Kratos
{
template<SizeType TDim>
void MetricFastInit<TDim>::Execute()
{
    KRATOS_TRY;

    const TensorArrayType zero_array(3 * (TDim - 1), 0.0);

    // We iterate over the nodes
    NodesArrayType& r_nodes_array = mrThisModelPart.Nodes();
    const auto it_node_begin = r_nodes_array.begin();

    // Tensor variable definition
    const Variable<TensorArrayType>& r_tensor_variable = KratosComponents<Variable<TensorArrayType>>::Get("METRIC_TENSOR_"+std::to_string(TDim)+"D");

    // Iterate over nodes
    IndexPartition<std::size_t>(r_nodes_array.size()).for_each(
        [&it_node_begin,&r_tensor_variable,&zero_array](std::size_t i_node)
        { (it_node_begin + i_node)->SetValue(r_tensor_variable, zero_array); }
    );

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template class MetricFastInit<2>;
template class MetricFastInit<3>;

}
