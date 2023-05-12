// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Sergio Jimenez Reyes
//

#include "constitutive_laws_application_variables.h"
#include "custom_processes/advance_in_time_high_cycle_fatigue_process.h"
#include "utilities/parallel_utilities.h"
#include "utilities/atomic_utilities.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

// Constructor
	AdvanceInTimeHighCycleFatigueProcess::AdvanceInTimeHighCycleFatigueProcess(
		ModelPart &rModelPart,
		Parameters ThisParameters) : mrModelPart(rModelPart), mThisParameters(ThisParameters)
{
}

/***********************************************************************************/
/***********************************************************************************/

void AdvanceInTimeHighCycleFatigueProcess::Execute()
{

    auto& process_info = mrModelPart.GetProcessInfo();
    bool cycle_found = false;
    std::vector<double> damage;
    process_info[ADVANCE_STRATEGY_APPLIED] = false;

    if (!process_info[DAMAGE_ACTIVATION]) {


        KRATOS_ERROR_IF(mrModelPart.NumberOfElements() == 0) << "The number of elements in the domain is zero. The process can not be applied."<< std::endl;

        for (auto& r_elem : mrModelPart.Elements()) {
            unsigned int number_of_ip = r_elem.GetGeometry().IntegrationPoints(r_elem.GetIntegrationMethod()).size();
            r_elem.CalculateOnIntegrationPoints(DAMAGE, damage, process_info);
            for (unsigned int i = 0; i < number_of_ip; i++) {
                    if (damage[i] > 0.0) {
                        process_info[DAMAGE_ACTIVATION] = true;
                        break;
                    }
            }
        }
    }

    this->CyclePeriodPerIntegrationPoint(cycle_found);  //This method detects if a cycle has finished somewhere in the model and
                                                        //computes the time period of the cycle that has just finished.

    if (cycle_found) {  //If a cycle has finished then it is possible to apply the advancing strategy
        bool advancing_strategy = false;
        this->StableConditionForAdvancingStrategy(advancing_strategy, process_info[DAMAGE_ACTIVATION]);  //Check if the conditions are optimal to apply the advancing strategy in
                                                                        //terms of max stress and reversion factor variation.
        if (advancing_strategy) {
            double increment = 0.0;
            if (!process_info[DAMAGE_ACTIVATION]) {
                this->TimeIncrement(increment);
                if (increment > 0.0) {
                    this->TimeAndCyclesUpdate(increment);
                    process_info[ADVANCE_STRATEGY_APPLIED] = true;
                }
            } else {
                increment = mThisParameters["fatigue"]["advancing_strategy_damage"].GetDouble();
                this->TimeAndCyclesUpdate(increment);
                process_info[ADVANCE_STRATEGY_APPLIED] = true;
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void AdvanceInTimeHighCycleFatigueProcess::CyclePeriodPerIntegrationPoint(bool& rCycleFound)
{
    auto& process_info = mrModelPart.GetProcessInfo();
    std::vector<bool> cycle_identifier;
    std::vector<double> previous_cycle_time;    //time when the previous cycle finished. It is used to obtain the new period for the current cycle
    std::vector<double> period;
    double time = process_info[TIME];

    KRATOS_ERROR_IF(mrModelPart.NumberOfElements() == 0) << "The number of elements in the domain is zero. The process can not be applied."<< std::endl;

    for (auto& r_elem : mrModelPart.Elements()) {   //This loop is done for all the integration points even if the cycle has not changed
        unsigned int number_of_ip = r_elem.GetGeometry().IntegrationPoints(r_elem.GetIntegrationMethod()).size();
        r_elem.CalculateOnIntegrationPoints(CYCLE_INDICATOR, cycle_identifier, process_info);
        r_elem.CalculateOnIntegrationPoints(PREVIOUS_CYCLE, previous_cycle_time, process_info);
        r_elem.CalculateOnIntegrationPoints(CYCLE_PERIOD, period, process_info);

        // Check if the HCF CL is being used. Otherwise there is no reason for entering into the advance in time strategy
        std::vector<ConstitutiveLaw::Pointer> constitutive_law_vector(number_of_ip);
        r_elem.CalculateOnIntegrationPoints(CONSTITUTIVE_LAW,constitutive_law_vector,process_info);

        const bool is_fatigue = constitutive_law_vector[0]->Has(CYCLE_INDICATOR);

        if (is_fatigue) {
            for (unsigned int i = 0; i < number_of_ip; i++) {
                    if (cycle_identifier[i]) {
                        period[i] = time - previous_cycle_time[i];
                        previous_cycle_time[i] = time;
                        rCycleFound = true;
                    }
            }
            r_elem.SetValuesOnIntegrationPoints(PREVIOUS_CYCLE, previous_cycle_time, process_info);
            r_elem.SetValuesOnIntegrationPoints(CYCLE_PERIOD, period, process_info);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void AdvanceInTimeHighCycleFatigueProcess::StableConditionForAdvancingStrategy(bool& rAdvancingStrategy, bool DamageIndicator)
{
    rAdvancingStrategy = false;
    auto& process_info = mrModelPart.GetProcessInfo();
    std::vector<double> max_stress_rel_error;
    std::vector<double> rev_factor_rel_error;
    std::vector<double> s_th;
    std::vector<double> max_stress;

    double acumulated_max_stress_rel_error = 0.0;
    double acumulated_rev_factor_rel_error = 0.0;
    bool fatigue_in_course = false;

    KRATOS_ERROR_IF(mrModelPart.NumberOfElements() == 0) << "The number of elements in the domain is zero. The process can not be applied."<< std::endl;

    for (auto& r_elem : mrModelPart.Elements()) {   //This loop is done for all the integration points even if the cycle has not changed
                                                    //This loop is done for all the integration points even if the cycle has not changed
                                                    //in order to guarantee the stable condition in the WHOLE model (when Smax > Sth)
        unsigned int number_of_ip = r_elem.GetGeometry().IntegrationPoints(r_elem.GetIntegrationMethod()).size();
        r_elem.CalculateOnIntegrationPoints(MAX_STRESS_RELATIVE_ERROR, max_stress_rel_error, process_info);
        r_elem.CalculateOnIntegrationPoints(REVERSION_FACTOR_RELATIVE_ERROR, rev_factor_rel_error, process_info);
        r_elem.CalculateOnIntegrationPoints(THRESHOLD_STRESS, s_th, process_info);
        r_elem.CalculateOnIntegrationPoints(MAX_STRESS, max_stress, process_info);

        for (unsigned int i = 0; i < number_of_ip; i++) {
            if (max_stress[i] > s_th[i]) {
                fatigue_in_course = true;
                acumulated_max_stress_rel_error += max_stress_rel_error[i];
                acumulated_rev_factor_rel_error += rev_factor_rel_error[i];
            }
        }
    }
    if ((acumulated_max_stress_rel_error < 1e-4 && acumulated_rev_factor_rel_error < 1e-4 && fatigue_in_course) || (DamageIndicator && acumulated_max_stress_rel_error < 1e-3 && acumulated_rev_factor_rel_error < 1e-3 && fatigue_in_course)) {
        rAdvancingStrategy = true;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void AdvanceInTimeHighCycleFatigueProcess::TimeIncrement(double& rIncrement)
{
    auto& process_info = mrModelPart.GetProcessInfo();
    std::vector<double>  cycles_to_failure_element;
    std::vector<int>  local_number_of_cycles;
    std::vector<double> period;
    std::vector<double> s_th;
    std::vector<double> max_stress;
    double min_time_increment;
    double time = process_info[TIME];
    bool current_constraints_process_list_detected = false;

    double user_avancing_time = mThisParameters["fatigue"]["advancing_strategy_time"].GetDouble();
    double user_avancing_cycles = mThisParameters["fatigue"]["advancing_strategy_cycles"].GetDouble();
    std::vector<std::string> constraints_list = mThisParameters["fatigue"]["constraints_process_list"].GetStringArray();
    double model_part_final_time = mThisParameters["problem_data"]["end_time"].GetDouble();
    for (unsigned int i = 0; i < constraints_list.size(); i++) {
        for (unsigned int j = 0; j < mThisParameters["processes"]["constraints_process_list"].size(); j++) {
            std::string model_part_name = mThisParameters["processes"]["constraints_process_list"][j]["Parameters"]["model_part_name"].GetString();
            double model_part_end_time = mThisParameters["processes"]["constraints_process_list"][j]["Parameters"]["interval"][1].GetDouble();
            if (constraints_list[i] == model_part_name && time <= model_part_end_time && !current_constraints_process_list_detected) {
                model_part_final_time = model_part_end_time;
                current_constraints_process_list_detected = true;
            }
        }
        if (current_constraints_process_list_detected) {
            break;
        }
    }

    double model_part_advancing_time = model_part_final_time - time;
    min_time_increment = std::min(user_avancing_time, model_part_advancing_time);

    // KRATOS_ERROR_IF(mrModelPart.NumberOfElements() == 0) << "The number of elements in the domain is zero. The process can not be applied."<< std::endl;

    for (auto& r_elem : mrModelPart.Elements()) {
        unsigned int number_of_ip = r_elem.GetGeometry().IntegrationPoints(r_elem.GetIntegrationMethod()).size();
        r_elem.CalculateOnIntegrationPoints(CYCLES_TO_FAILURE, cycles_to_failure_element, process_info);
        r_elem.CalculateOnIntegrationPoints(LOCAL_NUMBER_OF_CYCLES, local_number_of_cycles, process_info);
        r_elem.CalculateOnIntegrationPoints(CYCLE_PERIOD, period, process_info);
        r_elem.CalculateOnIntegrationPoints(THRESHOLD_STRESS, s_th, process_info);
        r_elem.CalculateOnIntegrationPoints(MAX_STRESS, max_stress, process_info);
        for (unsigned int i = 0; i < number_of_ip; i++) {
            if (max_stress[i] > s_th[i]) {  //This is used to guarantee that only IP in fatigue conditions are taken into account
                double Nf_conversion_to_time = (cycles_to_failure_element[i] - local_number_of_cycles[i]) * period[i];
                double user_avancing_cycles_conversion_to_time = user_avancing_cycles * period[i];
                if (Nf_conversion_to_time < min_time_increment) {
                    min_time_increment = Nf_conversion_to_time;
                }
                if (user_avancing_cycles_conversion_to_time < min_time_increment) {
                    min_time_increment = user_avancing_cycles_conversion_to_time;

                }
            }
        }
    }
	rIncrement = min_time_increment;
}

/***********************************************************************************/
/***********************************************************************************/

void AdvanceInTimeHighCycleFatigueProcess::TimeAndCyclesUpdate(const double Increment)
{
    auto& r_process_info = mrModelPart.GetProcessInfo();

    // KRATOS_ERROR_IF(mrModelPart.NumberOfElements() == 0) << "The number of elements in the domain is zero. The process can not be applied."<< std::endl;

    for (auto& r_elem : mrModelPart.Elements()) {
        std::vector<bool> cycle_identifier;
        std::vector<int>  local_number_of_cycles;
        std::vector<int>  global_number_of_cycles;
        std::vector<double> period;
        double time_increment = 0.0;
        std::vector<double> previous_cycle_time;    //time when the previous cycle finished. It is used to obtain the new period for the current cycle

        unsigned int number_of_ip = r_elem.GetGeometry().IntegrationPoints(r_elem.GetIntegrationMethod()).size();
        r_elem.CalculateOnIntegrationPoints(CYCLE_INDICATOR, cycle_identifier, r_process_info);
        r_elem.CalculateOnIntegrationPoints(LOCAL_NUMBER_OF_CYCLES, local_number_of_cycles, r_process_info);
        r_elem.CalculateOnIntegrationPoints(NUMBER_OF_CYCLES, global_number_of_cycles, r_process_info);
        r_elem.CalculateOnIntegrationPoints(CYCLE_PERIOD, period, r_process_info);
        r_elem.CalculateOnIntegrationPoints(PREVIOUS_CYCLE, previous_cycle_time, r_process_info);

        // Check if the HCF CL is being used. Otherwise there is no reason for entering into the advance in time strategy
        std::vector<ConstitutiveLaw::Pointer> constitutive_law_vector(number_of_ip);
        r_elem.CalculateOnIntegrationPoints(CONSTITUTIVE_LAW,constitutive_law_vector,r_process_info);

        const bool is_fatigue = constitutive_law_vector[0]->Has(CYCLE_INDICATOR);

        if (is_fatigue) {
            for (unsigned int i = 0; i < number_of_ip; i++) {
                unsigned int local_cycles_increment;
                if (period[i] == 0.0) {
                    local_cycles_increment = 0;
                } else {
                    local_cycles_increment = std::trunc(Increment / period[i]);
                    time_increment = std::trunc(Increment / period[i]) * period[i];
                    previous_cycle_time[i] += time_increment;
                }
                local_number_of_cycles[i] += local_cycles_increment;
                global_number_of_cycles[i] += local_cycles_increment;
            }
            r_elem.SetValuesOnIntegrationPoints(LOCAL_NUMBER_OF_CYCLES, local_number_of_cycles, r_process_info);
            r_elem.SetValuesOnIntegrationPoints(NUMBER_OF_CYCLES, global_number_of_cycles, r_process_info);
            r_elem.SetValuesOnIntegrationPoints(PREVIOUS_CYCLE, previous_cycle_time, r_process_info);
            // #pragma omp critical
            r_process_info[TIME_INCREMENT] = time_increment;
        }
    }
}

} // namespace Kratos