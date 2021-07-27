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
#include "utilities/variable_utils.h"
#include "utilities/atomic_utilities.h"
#include "conversion_utilities.h"

namespace Kratos {


void ConversionUtilities::ConvertPressureToForces(ModelPart& model_part_interface)
{
    // initialize Forces
    /* VariableUtils().SetHistoricalVariableToZero(FORCE, model_part_interface.Nodes());

    block_for_each(model_part_interface.Elements(),[&](Element& rElement){
        const array_1d<double, 3>& elem_force = rElement.GetGeometry().Area() * rElement.GetValue(PRESSURE);

        const std::size_t num_nodes = rElement.GetGeometry().PointsNumber();

        for (auto& r_node : rElement.GetGeometry().Points()){
            AtomicAdd( noalias(r_node.FastGetSolutionStepValue(FORCE)), (elem_force/static_cast<double>(num_nodes)) );
        }
    });

    model_part_interface.GetCommunicator().AssembleCurrentData(FORCE); */
}

void ConversionUtilities::ConvertElementalDataToNodalData(ModelPart& model_part_interface)
{
    // initialize Forces
    VariableUtils().SetHistoricalVariableToZero(FORCE, model_part_interface.Nodes());

    block_for_each(model_part_interface.Elements(), [&](Element& rElement){
        const array_1d<double, 3>& elem_force =  rElement.GetValue(FORCE);

        const std::size_t num_nodes = rElement.GetGeometry().PointsNumber();

        for (auto& r_node : rElement.GetGeometry().Points()){
            //AtomicAdd( noalias(r_node.FastGetSolutionStepValue(FORCE)), (elem_force/static_cast<double>(num_nodes)) );
            r_node.FastGetSolutionStepValue(FORCE) += (elem_force/static_cast<double>(num_nodes)) ;
        }
    });

    auto total_elemental_force = VariableUtils().SumElementVectorVariable(FORCE, model_part_interface);
    array_1d<double, 3> total_nodal_forces = VariableUtils().SumHistoricalVariable<array_1d<double, 3> >(FORCE, model_part_interface);

    std::cout<< "Total Force values on the given interface (Elemental Calculation)= " << total_elemental_force << std::endl;
    std::cout<< "Total Force values on the given interface (Nodal Calculation)= " << total_nodal_forces << std::endl;

    //model_part_interface.GetCommunicator().AssembleCurrentData(FORCE);
}


}  // namespace Kratos.


