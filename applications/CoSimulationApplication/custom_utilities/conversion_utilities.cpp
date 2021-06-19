//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Author1 Fullname
//                   Author2 Fullname
//


// System includes


// External includes


// Project includes
#include "includes/variables.h"
#include "utilities/parallel_utilities.h"
#include "utilities/variable_utilities.h"
#include "conversion_utilities.h"


namespace Kratos {


static void ConversionUtilities::ConvertPressureToForces()
{
    // initialize Forces
    // ...

    VariableUtils().SetHistoricalVariableToZero(
        FORCE,
        ModelPart.Nodes());

    block_for_each(ModelPart.Elements(), (Element& rElement)[]{
        const array_1d<double, 3>& elem_force = rElement.GetGeometry().Area() * rElement.GetValue(PRESSURE);

        const std::size_t num_nodes = rElement.GetGeometry().PointsNumber();

        for (auto& r_node : rElement.GetGeometry().Points()){
            AtomicAdd(noalias(r_node.FastGetSolutionStepValue(FORCE)) += elem_force/static_cast<double>(num_nodes));
        }
    });

    ModelPart.GetCommunicator().AssembleCurrentData(FORCE);
}


}  // namespace Kratos.


