// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "contact_structural_mechanics_application_variables.h"
#include "custom_processes/alm_variables_calculation_process.h"

namespace Kratos
{
void ALMVariablesCalculationProcess::Execute()
{
    KRATOS_TRY

    /* We compute the penalty factor */

    // Auxiliar values
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
    ConditionsArrayType& r_conditions_array = mrThisModelPart.Conditions();

    // First condition iterator
    const auto it_cond_begin = r_conditions_array.begin();

    #pragma omp parallel for reduction(+:total_volume_slave, total_area_slave, mean_young_modulus_slave, mean_nodal_h_slave, total_volume_master, total_area_master, mean_young_modulus_master, mean_nodal_h_master)
    for(int i = 0; i < static_cast<int>(r_conditions_array.size()); ++i) {
        auto it_cond = it_cond_begin + i;

        // We get the condition geometry
        const GeometryType& r_this_geometry = it_cond->GetGeometry();

        // We get the values from the condition
        const Properties& r_properties = it_cond->GetProperties();

        const double young_modulus = r_properties.Has(YOUNG_MODULUS) ? r_properties[YOUNG_MODULUS] : 0.0;
        KRATOS_WARNING_IF("ALMVariablesCalculationProcess", young_modulus < epsilon) << "Not assigned Young modulus. Using zero" << std::endl;
        const double element_volume = r_this_geometry.DomainSize();
        KRATOS_WARNING_IF("ALMVariablesCalculationProcess", element_volume < epsilon) << "Null element volume. Please check!!!" << std::endl;

        // We get the values from the condition
        const double condition_area = r_this_geometry.DomainSize();
        const double nodal_condition_area = condition_area/static_cast<double>(r_this_geometry.size());

        // Slave conditions
        if (it_cond->Is(SLAVE)) {
            total_volume_slave += element_volume;
            total_area_slave += condition_area;
            mean_young_modulus_slave += young_modulus * element_volume;

            for (auto& r_node : r_this_geometry) {
                mean_nodal_h_slave += r_node.FastGetSolutionStepValue(mrNodalLengthVariable) * nodal_condition_area;
            }
        } else if (it_cond->Is(MASTER)) { // Master conditions
            total_volume_master += element_volume;
            total_area_master += condition_area;
            mean_young_modulus_master += young_modulus * element_volume;

            for (auto& r_node : r_this_geometry) {
                mean_nodal_h_master += r_node.FastGetSolutionStepValue(mrNodalLengthVariable) * nodal_condition_area;
            }
        } else { // Not defined
            total_volume_slave += element_volume;
            total_area_slave += condition_area;
            mean_young_modulus_slave += young_modulus * element_volume;

            total_volume_master += element_volume;
            total_area_master += condition_area;
            mean_young_modulus_master += young_modulus * element_volume;

            for (auto& r_node : r_this_geometry) {
                mean_nodal_h_slave += r_node.FastGetSolutionStepValue(mrNodalLengthVariable) * nodal_condition_area;
                mean_nodal_h_master += r_node.FastGetSolutionStepValue(mrNodalLengthVariable) * nodal_condition_area;
            }
        }
    }

    // Now we divide between the total areas and volumes
    mean_nodal_h_slave        /= (total_area_slave   + tolerance);
    mean_young_modulus_slave  /= (total_volume_slave + tolerance);

    mean_nodal_h_master       /= (total_area_master   + tolerance);
    mean_young_modulus_master /= (total_volume_master + tolerance);

    // Finally we compute the penalty factor
    const double penalty_parameter_slave  = (mean_young_modulus_slave > epsilon) ? mFactorStiffness * mean_young_modulus_slave/(mean_nodal_h_slave + tolerance) : mFactorStiffness * mean_young_modulus_master/(mean_nodal_h_master + tolerance);
    const double scale_factor_slave       = mPenaltyScale * mFactorStiffness * mean_young_modulus_slave/(mean_nodal_h_slave + tolerance);
    const double penalty_parameter_master = mFactorStiffness * mean_young_modulus_master/(mean_nodal_h_master + tolerance);
    const double scale_factor_master      = mPenaltyScale * mFactorStiffness * mean_young_modulus_master/(mean_nodal_h_master + tolerance);

    KRATOS_INFO("ALM Values") << "The potential parameters for penalty and scale factor are: \n" << "PENALTY SLAVE:\t" << penalty_parameter_slave << "\tPENALTY MASTER:\t" << penalty_parameter_master << "\nSCALE FACTOR SLAVE:\t" << scale_factor_slave << "\tSCALE FACTOR MASTER:\t" << scale_factor_master << std::endl;

    // NOTE: > or <? , we are supposed to take the smallest of the values (less stiff)
    const double final_penalty = (penalty_parameter_slave < penalty_parameter_master) ? penalty_parameter_slave : penalty_parameter_master;
    const double final_scale_factor = (scale_factor_slave < scale_factor_master) ? scale_factor_slave : scale_factor_master;
    mrThisModelPart.GetProcessInfo().SetValue(INITIAL_PENALTY, final_penalty);
    mrThisModelPart.GetProcessInfo().SetValue(SCALE_FACTOR, final_scale_factor);

    KRATOS_CATCH("")
} // class ALMVariablesCalculationProcess
} // namespace Kratos
