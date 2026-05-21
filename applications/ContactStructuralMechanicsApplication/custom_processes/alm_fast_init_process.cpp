// KRATOS    ______            __             __  _____ __                  __                   __
//          / ____/___  ____  / /_____ ______/ /_/ ___// /________  _______/ /___  ___________ _/ /
//         / /   / __ \/ __ \/ __/ __ `/ ___/ __/\__ \/ __/ ___/ / / / ___/ __/ / / / ___/ __ `/ / 
//        / /___/ /_/ / / / / /_/ /_/ / /__/ /_ ___/ / /_/ /  / /_/ / /__/ /_/ /_/ / /  / /_/ / /  
//        \____/\____/_/ /_/\__/\__,_/\___/\__//____/\__/_/   \__,_/\___/\__/\__,_/_/   \__,_/_/  MECHANICS
//
//  License:         BSD License
//                   license: ContactStructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "utilities/parallel_utilities.h"
#include "contact_structural_mechanics_application_variables.h"
#include "custom_processes/alm_fast_init_process.h"
#include "utilities/atomic_utilities.h"

namespace Kratos
{
void ALMFastInit::Execute()
{
    KRATOS_TRY;

    // First we reorder the conditions ids (may methods and utilities assume that conditions are ordered)
    ConditionsArrayType& r_root_conditions_array = mrThisModelPart.GetRootModelPart().Conditions();
    const auto it_root_cond_begin = r_root_conditions_array.begin();

    IndexPartition<std::size_t>(r_root_conditions_array.size()).for_each([&it_root_cond_begin](std::size_t i) {
        auto it_cond = it_root_cond_begin + i;
        it_cond->SetId(i + 1);
    });

    // We differentiate between frictional or frictionless
    const bool is_frictional = mrThisModelPart.Is(SLIP);

    // We initialize the penalty parameter
    const double epsilon = mrThisModelPart.GetProcessInfo()[INITIAL_PENALTY];

    // Auxiliary zero array
    const array_1d<double, 3> zero_array = ZeroVector(3);

    // We iterate over the nodes
    NodesArrayType& r_nodes_array = mrThisModelPart.Nodes();
    block_for_each(r_nodes_array, [&](Node& rNode) {
        const bool is_slave = rNode.IsDefined(SLAVE) ? rNode.Is(SLAVE) : true;
        if (is_slave) {
            // Weighted values
            rNode.FastGetSolutionStepValue(WEIGHTED_GAP) = 0.0;
            if (is_frictional) {
                noalias(rNode.FastGetSolutionStepValue(WEIGHTED_SLIP)) = zero_array;
            }

            // Penalty parameter
            rNode.SetValue(INITIAL_PENALTY, epsilon);

            // Auxiliary values
            rNode.SetValue(DYNAMIC_FACTOR, 1.0);
            rNode.SetValue(AUGMENTED_NORMAL_CONTACT_PRESSURE, 0.0);
            if (is_frictional) {
                rNode.SetValue(AUGMENTED_TANGENT_CONTACT_PRESSURE, zero_array);
            }
        }
    });

    // Now we iterate over the conditions
    ConditionsArrayType& r_conditions_array = mrThisModelPart.Conditions();
    block_for_each(r_conditions_array,[&](Condition& rCond) {
        rCond.SetValue(NORMAL, zero_array); // The normal and tangents vectors
    });

    if (is_frictional) {
        // We initialize the frictional coefficient. The evolution of the frictional coefficient it is supposed to be controled by a law
        block_for_each(r_nodes_array, [&](Node& rNode) {
            rNode.SetValue(FRICTION_COEFFICIENT, 0.0);
            rNode.SetValue(NODAL_AREA, 0.0);
        });

        // Iterate over submodelparts
        for (auto& r_sub_model_part : mrThisModelPart.SubModelParts()) {
            // Now we iterate over the conditions
            ConditionsArrayType& r_contact_conditions_array = r_sub_model_part.Conditions();
            block_for_each(r_contact_conditions_array,[&](Condition& rCond) {
                auto p_prop = rCond.pGetProperties();
                auto& r_geom = rCond.GetGeometry();

                for (auto& r_node : r_geom) {
                    double& r_nodal_area = r_node.GetValue(NODAL_AREA);
                    AtomicAdd(r_nodal_area, 1.0);
                }

                if (p_prop->Has(FRICTION_COEFFICIENT)) {
                    const double friction_coefficient = p_prop->GetValue(FRICTION_COEFFICIENT);
                    for (auto& r_node : r_geom) {
                        double& r_friction_coefficient = r_node.GetValue(FRICTION_COEFFICIENT);
                        AtomicAdd(r_friction_coefficient, friction_coefficient);
                    }
                } else {
                    KRATOS_WARNING("ALMFastInit") << "WARNING:: Friction coefficient not defined, zero will be considered" << std::endl;
                }
            });
        }

        block_for_each(r_nodes_array, [&](Node& rNode) {
            double& friction_coefficient = rNode.GetValue(FRICTION_COEFFICIENT);
            friction_coefficient /= rNode.GetValue(NODAL_AREA);
        });
    }

    KRATOS_CATCH("");
} // class ALMFastInit
} // namespace Kratos
