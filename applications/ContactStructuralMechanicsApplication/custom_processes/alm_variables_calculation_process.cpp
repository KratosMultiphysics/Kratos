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
#include "utilities/reduction_utilities.h"
#include "contact_structural_mechanics_application_variables.h"
#include "custom_processes/alm_variables_calculation_process.h"

namespace Kratos
{
void ALMVariablesCalculationProcess::Execute()
{
    KRATOS_TRY

    /* We compute the penalty factor */

    // Auxiliary values
    const double epsilon = std::numeric_limits<double>::epsilon();
    const double tolerance = 1.0e-12;

    // We initialize the mean values
    double mean_young_modulus_slave  = 0.0;
    double mean_nodal_h_slave        = 0.0;
    double mean_young_modulus_master = 0.0;
    double mean_nodal_h_master       = 0.0;

    // We initialize the total areas and volumes
    double total_volume_slave  = 0.0;
    double total_area_slave    = 0.0;
    double total_volume_master = 0.0;
    double total_area_master   = 0.0;

    // Now we iterate over the conditions to calculate the nodal area
    using EightReduction = CombinedReduction<SumReduction<double>, SumReduction<double>, SumReduction<double>, SumReduction<double>,SumReduction<double>, SumReduction<double>, SumReduction<double>, SumReduction<double>>;
    std::tie(total_volume_slave, total_area_slave, mean_young_modulus_slave, mean_nodal_h_slave, total_volume_master, total_area_master, mean_young_modulus_master, mean_nodal_h_master) = block_for_each<EightReduction>(mrThisModelPart.Conditions(), [&](Condition& rCond) {
        // We get the condition geometry
        const GeometryType& r_this_geometry = rCond.GetGeometry();

        // We get the values from the condition
        const Properties& r_properties = rCond.GetProperties();

        const double young_modulus = r_properties.Has(YOUNG_MODULUS) ? r_properties[YOUNG_MODULUS] : 0.0;
        KRATOS_WARNING_IF("ALMVariablesCalculationProcess", young_modulus < epsilon) << "Not assigned Young modulus. Using zero" << std::endl;
        const double element_volume = r_this_geometry.DomainSize();
        KRATOS_WARNING_IF("ALMVariablesCalculationProcess", element_volume < epsilon) << "Null element volume. Please check!!!" << std::endl;

        // We get the values from the condition
        const double condition_area = r_this_geometry.DomainSize();
        const double nodal_condition_area = condition_area/static_cast<double>(r_this_geometry.size());

        // Slave conditions
        if (rCond.Is(SLAVE)) {
            double aux_sum = 0.0;
            for (auto& r_node : r_this_geometry) {
                aux_sum += r_node.FastGetSolutionStepValue(mrNodalLengthVariable) * nodal_condition_area;
            }
            return std::make_tuple(element_volume,condition_area,young_modulus * element_volume,aux_sum,0.0,0.0,0.0,0.0);
        } else if (rCond.Is(MASTER)) { // Master conditions
            double aux_sum = 0.0;
            for (auto& r_node : r_this_geometry) {
                aux_sum += r_node.FastGetSolutionStepValue(mrNodalLengthVariable) * nodal_condition_area;
            }
            return std::make_tuple(0.0,0.0,0.0,0.0,element_volume,condition_area,young_modulus * element_volume,aux_sum);
        } else { // Not defined
            double aux_sum_slave = 0.0;
            double aux_sum_master = 0.0;
            for (auto& r_node : r_this_geometry) {
                aux_sum_slave += r_node.FastGetSolutionStepValue(mrNodalLengthVariable) * nodal_condition_area;
                aux_sum_master += r_node.FastGetSolutionStepValue(mrNodalLengthVariable) * nodal_condition_area;
            }
            return std::make_tuple(element_volume,condition_area,young_modulus * element_volume,aux_sum_slave,element_volume,condition_area,young_modulus * element_volume,aux_sum_master);
        }
        return std::make_tuple(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
    });

    // Now we divide between the total areas and volumes
    mean_nodal_h_slave        /= (total_area_slave   + tolerance);
    mean_young_modulus_slave  /= (total_volume_slave + tolerance);

    mean_nodal_h_master       /= (total_area_master   + tolerance);
    mean_young_modulus_master /= (total_volume_master + tolerance);

    // Prepare the string stream
    std::stringstream buffer;

    // Finally we compute the penalty factor
    if (mComputePenalty) {
        const double penalty_parameter_slave  = (mean_young_modulus_slave > epsilon) ? mFactorStiffness * mean_young_modulus_slave/(mean_nodal_h_slave + tolerance) : mFactorStiffness * mean_young_modulus_master/(mean_nodal_h_master + tolerance);
        const double penalty_parameter_master = mFactorStiffness * mean_young_modulus_master/(mean_nodal_h_master + tolerance);

        // Fill the buffer
        if (mComputeScaleFactor) buffer << "The potential parameters for penalty and scale factor are: \n";
        buffer << "PENALTY SLAVE:\t" << penalty_parameter_slave << "\tPENALTY MASTER:\t" << penalty_parameter_master;

        // NOTE: > or <? , we are supposed to take the smallest of the values (less stiff)
        const double final_penalty = (penalty_parameter_slave < penalty_parameter_master) ? penalty_parameter_slave : penalty_parameter_master;
        
        // Setting value
        mrThisModelPart.GetProcessInfo().SetValue(INITIAL_PENALTY, final_penalty);
    }

    // Finally we compute the scale factor
    if (mComputeScaleFactor) {
        const double scale_factor_slave  = mPenaltyScale * mFactorStiffness * mean_young_modulus_slave/(mean_nodal_h_slave + tolerance);
        const double scale_factor_master = mPenaltyScale * mFactorStiffness * mean_young_modulus_master/(mean_nodal_h_master + tolerance);

        // Fill the buffer
        if (!mComputePenalty) buffer << "The potential parameters for scale factor are: \n";
        buffer << "\nSCALE FACTOR SLAVE:\t" << scale_factor_slave << "\tSCALE FACTOR MASTER:\t" << scale_factor_master;

        // NOTE: > or <? , we are supposed to take the smallest of the values (less stiff)
        const double final_scale_factor = (scale_factor_slave < scale_factor_master) ? scale_factor_slave : scale_factor_master;

        // Setting value
        mrThisModelPart.GetProcessInfo().SetValue(SCALE_FACTOR, final_scale_factor);
    }

    KRATOS_INFO("ALM Values") << buffer.str() << std::endl;

    KRATOS_CATCH("")
} // class ALMVariablesCalculationProcess
} // namespace Kratos
