//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//                   Ashish Darekar
//

// System includes

// External includes

// Project includes
#include "utilities/parallel_utilities.h"
#include "utilities/variable_utils.h"
#include "utilities/atomic_utilities.h"
#include "conversion_utilities.h"

namespace Kratos {

void ConversionUtilities::ConvertElementalDataToNodalData(
    ModelPart& rModelPart,
    const Variable<array_1d<double,3> >& rElementalVariable,
    const Variable<array_1d<double,3> >& rNodalVariable )
{
    // prepare nodal variable
    VariableUtils().SetHistoricalVariableToZero(rNodalVariable, rModelPart.Nodes());

    block_for_each(rModelPart.Elements(), [&](Element& rElement){
        const array_1d<double, 3>& elem_rVariable =  rElement.GetValue(rElementalVariable);

        const std::size_t num_nodes = rElement.GetGeometry().PointsNumber();

        for (auto& r_node : rElement.GetGeometry().Points()){
            AtomicAddVector( r_node.FastGetSolutionStepValue(rNodalVariable), (elem_rVariable / static_cast<double>(num_nodes)) );
        }
    });

    rModelPart.GetCommunicator().AssembleCurrentData(rNodalVariable);
}

}  // namespace Kratos.
