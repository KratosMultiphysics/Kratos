//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//

// System includes

// External includes
#include "includes/define.h"
#include "includes/condition.h"
#include "includes/model_part.h"
#include "custom_utilities/FEMDEM_coupling_utilities.h"
#include "dem_structures_coupling_application_variables.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDim>
FEMDEMCouplingUtilities<TDim>::FEMDEMCouplingUtilities()
{
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDim>
void FEMDEMCouplingUtilities<TDim>::SaveStructuralSolution(
        ModelPart& rStructureModelPart
    ) 
{

    KRATOS_TRY

    const int number_of_nodes = static_cast<int>(rStructureModelPart.Nodes().size());
    ModelPart::NodesContainerType::iterator node_begin = rStructureModelPart.NodesBegin();

    #pragma omp parallel for
    for (int i = 0; i < number_of_nodes; i++) {

        ModelPart::NodesContainerType::iterator it_node = node_begin + i;

        array_1d<double,3>& r_current_velocity = it_node->FastGetSolutionStepValue(BACKUP_LAST_STRUCTURAL_VELOCITY);
        noalias(r_current_velocity) = it_node->FastGetSolutionStepValue(VELOCITY);

        array_1d<double,3>& r_current_displacement = it_node->FastGetSolutionStepValue(BACKUP_LAST_STRUCTURAL_DISPLACEMENT);
        noalias(r_current_displacement) = it_node->FastGetSolutionStepValue(DISPLACEMENT);

        array_1d<double,3>& r_smoothed_velocity = it_node->FastGetSolutionStepValue(SMOOTHED_STRUCTURAL_VELOCITY);
        noalias(r_smoothed_velocity) = 1.0 / 3.0 * (it_node->FastGetSolutionStepValue(VELOCITY) + 2.0 * it_node->FastGetSolutionStepValue(SMOOTHED_STRUCTURAL_VELOCITY, 1));
    }
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDim>
void FEMDEMCouplingUtilities<TDim>::InterpolateStructuralSolution(
    ModelPart& rStructureModelPart, 
    const double FemDeltaTime, 
    const double FemTime, 
    const double DemDeltaTime, 
    const double DemTime
    )
{
    KRATOS_TRY

    const double previous_FemTime = FemTime - FemDeltaTime;
    const double time_factor = (DemTime - previous_FemTime) / FemDeltaTime;
    const double previous_time_factor = (DemTime - DemDeltaTime - previous_FemTime) / FemDeltaTime;

    const int number_of_nodes = static_cast<int>(rStructureModelPart.Nodes().size());
    ModelPart::NodesContainerType::iterator node_begin = rStructureModelPart.NodesBegin();

    #pragma omp parallel for
    for (int i = 0; i < number_of_nodes; i++) {

        ModelPart::NodesContainerType::iterator it_node = node_begin + i;

        noalias(it_node->Coordinates()) = it_node->GetInitialPosition().Coordinates() + it_node->FastGetSolutionStepValue(DISPLACEMENT, 1) + (it_node->FastGetSolutionStepValue(BACKUP_LAST_STRUCTURAL_DISPLACEMENT) - it_node->FastGetSolutionStepValue(DISPLACEMENT, 1)) * time_factor;

        array_1d<double,3> &r_velocity = it_node->FastGetSolutionStepValue(VELOCITY);
        const array_1d<double,3> &previous_velocity = it_node->FastGetSolutionStepValue(SMOOTHED_STRUCTURAL_VELOCITY, 1);
        noalias(r_velocity) = previous_velocity + (it_node->FastGetSolutionStepValue(SMOOTHED_STRUCTURAL_VELOCITY) - previous_velocity) * time_factor;

        array_1d<double,3>& r_displacement = it_node->FastGetSolutionStepValue(DISPLACEMENT);
        noalias(r_displacement) = it_node->Coordinates() - it_node->GetInitialPosition().Coordinates();

        array_1d<double,3> previous_coordinates;
        noalias(previous_coordinates) = it_node->GetInitialPosition().Coordinates() + it_node->FastGetSolutionStepValue(DISPLACEMENT, 1) + (it_node->FastGetSolutionStepValue(BACKUP_LAST_STRUCTURAL_DISPLACEMENT) - it_node->FastGetSolutionStepValue(DISPLACEMENT, 1)) * previous_time_factor;

        array_1d<double,3>& delta_displacement = it_node->FastGetSolutionStepValue(DELTA_DISPLACEMENT);
        noalias(delta_displacement) = it_node->Coordinates() - previous_coordinates;
    }
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDim>
void FEMDEMCouplingUtilities<TDim>::RestoreStructuralSolution(
    ModelPart& rStructureModelPart
    )
{
    KRATOS_TRY

    const int number_of_nodes = static_cast<int>(rStructureModelPart.Nodes().size());
    ModelPart::NodesContainerType::iterator node_begin = rStructureModelPart.NodesBegin();

    #pragma omp parallel for
    for (int i = 0; i < number_of_nodes; i++) {

        ModelPart::NodesContainerType::iterator it_node = node_begin + i;

        array_1d<double,3>& r_velocity = it_node->FastGetSolutionStepValue(VELOCITY);
        noalias(r_velocity) = it_node->FastGetSolutionStepValue(BACKUP_LAST_STRUCTURAL_VELOCITY);

        array_1d<double,3>& r_displacement = it_node->FastGetSolutionStepValue(DISPLACEMENT);
        noalias(r_displacement) = it_node->FastGetSolutionStepValue(BACKUP_LAST_STRUCTURAL_DISPLACEMENT);

        noalias(it_node->Coordinates()) = it_node->GetInitialPosition().Coordinates() + it_node->FastGetSolutionStepValue(DISPLACEMENT);
    }

    KRATOS_CATCH("") 
}

/***********************************************************************************/
/***********************************************************************************/

template class FEMDEMCouplingUtilities<3>;
template class FEMDEMCouplingUtilities<2>;

} // namespace Kratos