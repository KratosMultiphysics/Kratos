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
#include "fem_to_dem_application_variables.h"
#include "DEM_application_variables.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

void FEMDEMCouplingUtilities::SaveStructuralSolution(
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

void FEMDEMCouplingUtilities::InterpolateStructuralSolution(
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

void FEMDEMCouplingUtilities::RestoreStructuralSolution(
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

void FEMDEMCouplingUtilities::AddExplicitImpulses(
    ModelPart& rStructureModelPart,
    const double DEMTimeStep
    )
{
    auto& r_sub_model_conditions = rStructureModelPart.GetSubModelPart("ContactForcesDEMConditions");
    const auto it_cond_begin = r_sub_model_conditions.ConditionsBegin();
    auto& r_process_info = rStructureModelPart.GetProcessInfo();
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(r_sub_model_conditions.Conditions().size()); i++) {
        auto it_cond = it_cond_begin + i;
        auto& r_geometry = it_cond->GetGeometry();
        auto& r_node = r_geometry[0];

        if (r_node.GetValue(IS_DEM)) {
            auto p_spheric_particle_associated = r_node.GetValue(DEM_PARTICLE_POINTER);
            array_1d<double,3>& r_explicit_impulse_node = r_node.FastGetSolutionStepValue(CONTACT_IMPULSE);
            array_1d<double,3>& r_explicit_impulse_DEM = (p_spheric_particle_associated->GetGeometry()[0]).FastGetSolutionStepValue(CONTACT_IMPULSE);
            array_1d<double,3> dem_forces;

            if (!r_process_info[DEMFEM_CONTACT]) {
                dem_forces = (p_spheric_particle_associated->GetGeometry()[0]).FastGetSolutionStepValue(TOTAL_FORCES);
                r_explicit_impulse_DEM += dem_forces * DEMTimeStep;
            } else { // In the DE-FE contact the force is stored at the FEM nodes
                auto& r_dem_forces_ball = (p_spheric_particle_associated->GetGeometry()[0]).FastGetSolutionStepValue(TOTAL_FORCES);
                auto& r_dem_forces_wall = r_node.FastGetSolutionStepValue(TOTAL_FORCES);
                r_explicit_impulse_DEM  += r_dem_forces_ball * DEMTimeStep;
                r_explicit_impulse_node += r_dem_forces_wall * DEMTimeStep;
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void FEMDEMCouplingUtilities::ComputeAndTranferAveragedContactTotalForces(
    ModelPart& rStructureModelPart,
    const double FEMtimeStep
    )
{
    auto& r_sub_model_conditions = rStructureModelPart.GetSubModelPart("ContactForcesDEMConditions");
    const auto it_cond_begin = r_sub_model_conditions.ConditionsBegin();
    auto& r_process_info = rStructureModelPart.GetProcessInfo();
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(r_sub_model_conditions.Conditions().size()); i++) {
        auto it_cond = it_cond_begin + i;
        auto& r_geometry = it_cond->GetGeometry();
        auto& r_node = r_geometry[0];

        if (r_node.GetValue(IS_DEM)) {
            array_1d<double, 3> dem_forces;
            auto p_spheric_particle_associated = r_node.GetValue(DEM_PARTICLE_POINTER);
            array_1d<double,3> explicit_impulse_node = r_node.FastGetSolutionStepValue(CONTACT_IMPULSE);
            array_1d<double,3> explicit_impulse_DEM = (p_spheric_particle_associated->GetGeometry()[0]).FastGetSolutionStepValue(CONTACT_IMPULSE);
            if (!r_process_info[DEMFEM_CONTACT]) {
                dem_forces = (p_spheric_particle_associated->GetGeometry()[0]).FastGetSolutionStepValue(CONTACT_IMPULSE) / FEMtimeStep;
            } else { // In the DE-FE contact the force is stored at the FEM nodes
                array_1d<double,3>& r_dem_forces_ball = (p_spheric_particle_associated->GetGeometry()[0]).FastGetSolutionStepValue(TOTAL_FORCES);
                array_1d<double,3>& r_dem_forces_wall = r_node.FastGetSolutionStepValue(TOTAL_FORCES);
                r_dem_forces_ball = explicit_impulse_DEM / FEMtimeStep;
                r_dem_forces_wall = explicit_impulse_node / FEMtimeStep;
                dem_forces = r_dem_forces_ball + r_dem_forces_wall;
            }
            it_cond->SetValue(FORCE_LOAD, dem_forces);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

} // namespace Kratos