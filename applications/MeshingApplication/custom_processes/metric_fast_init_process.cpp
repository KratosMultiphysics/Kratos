// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                       license: MeshingApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
//

// System includes

// External includes

// Project includes
#include "custom_processes/metric_fast_init_process.h"

namespace Kratos
{
template<unsigned int TDim>
void MetricFastInit<TDim>::Execute()
{
    KRATOS_TRY;

    constexpr unsigned int size = TDim == 2  ? 3: 6;

    const array_1d<double, size> zerovector(size, 0.0);

    // We iterate over the node
    NodesArrayType& nodes_array = mrThisModelPart.Nodes();
    int num_nodes = mrThisModelPart.NumberOfNodes();

    #pragma omp parallel for firstprivate(zerovector)
    for(int i = 0; i < num_nodes; ++i) 
    {
        auto it_node = nodes_array.begin() + i;

        // The metric
        it_node->SetValue(MMG_METRIC, zerovector);
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template class MetricFastInit<2>;
template class MetricFastInit<3>;

}
