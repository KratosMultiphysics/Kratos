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

    // We iterate over the nodes
    NodesArrayType& nodes_array = mrThisModelPart.Nodes();

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) 
        (nodes_array.begin() + i)->SetValue(MMG_METRIC, ZeroVector(size));

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template class MetricFastInit<2>;
template class MetricFastInit<3>;

}
