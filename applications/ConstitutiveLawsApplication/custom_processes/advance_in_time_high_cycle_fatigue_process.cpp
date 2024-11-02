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
    // std::vector<double> previous_cycle_damage;
    std::vector<double> cycles_to_failure_element;
    std::vector<int> local_number_of_cycles;
    
    // process_info[NO_LINEARITY_ACTIVATION] = false;
    process_info[ADVANCE_STRATEGY_APPLIED] = false;

    this->CyclicLoad();  //This method checks if a cyclic load is being applied.

    KRATOS_ERROR_IF(mrModelPart.NumberOfElements() == 0) << "The number of elements in the domain is zero. The process can not be applied."<< std::endl;

    for (auto& r_elem : mrModelPart.Elements()) {
        unsigned int number_of_ip = r_elem.GetGeometry().IntegrationPoints(r_elem.GetIntegrationMethod()).size();
        r_elem.CalculateOnIntegrationPoints(DAMAGE, damage, process_info);
        // r_elem.CalculateOnIntegrationPoints(PREVIOUS_CYCLE_DAMAGE, previous_cycle_damage, process_info);
        r_elem.CalculateOnIntegrationPoints(CYCLES_TO_FAILURE, cycles_to_failure_element, process_info);
        r_elem.CalculateOnIntegrationPoints(LOCAL_NUMBER_OF_CYCLES, local_number_of_cycles, process_info);

        std::vector<ConstitutiveLaw::Pointer> constitutive_law_vector(number_of_ip);
        r_elem.CalculateOnIntegrationPoints(CONSTITUTIVE_LAW,constitutive_law_vector,process_info);

        const bool is_fatigue = constitutive_law_vector[0]->Has(CYCLE_INDICATOR);

        if (is_fatigue) {    
            for (unsigned int i = 0; i < number_of_ip; i++) {
                
                // double delta_damage = damage[i] - previous_cycle_damage[i];
                double delta_number_of_cycles = cycles_to_failure_element[i] - local_number_of_cycles[i];               

                if ((damage[i] > 0.0) || (delta_number_of_cycles <= 1.0)) {
                    process_info[NO_LINEARITY_ACTIVATION] = true;
                    break;
                }
            }
        }

        if (process_info[NO_LINEARITY_ACTIVATION]){
            break;
        }         
    }

    this->CyclePeriodPerIntegrationPoint(cycle_found);  //This method detects if a cycle has finished somewhere in the model and
                                                        //computes the time period of the cycle that has just finished.

    if (cycle_found) {  //If a cycle has finished then it is possible to apply the advancing strategy
   
        bool advancing_strategy = false;
        this->StableConditionForAdvancingStrategy(advancing_strategy, process_info[NO_LINEARITY_ACTIVATION]);  //Check if the conditions are optimal to apply the advancing strategy in                                                                                                                    
        double increment = 0.0;
        if (advancing_strategy) {
            this->TimeIncrement(increment);
            if (increment > 0.0) {
                this->TimeAndCyclesUpdate(increment);
                process_info[ADVANCE_STRATEGY_APPLIED] = true;
            }
        }
        // process_info[NO_LINEARITY_ACTIVATION] = false;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void AdvanceInTimeHighCycleFatigueProcess::CyclicLoad()
{
    auto& process_info = mrModelPart.GetProcessInfo();
    double time = process_info[TIME];
    double delta_time;
    if (mThisParameters["solver_settings"]["time_stepping"].Has("time_step")) {
        delta_time = mThisParameters["solver_settings"]["time_stepping"]["time_step"].GetDouble();
    } else if (mThisParameters["solver_settings"]["time_stepping"].Has("time_step_table")) {
        const Matrix time_step_table = mThisParameters["solver_settings"]["time_stepping"]["time_step_table"].GetMatrix();
        KRATOS_ERROR << "Advance in time process not prepared yet for time_step_table!" << std::endl;
    } else {
        KRATOS_ERROR << "Time stepping not defined!" << std::endl;
    }

    //Monotonic (false) or cyclic (true) load. If a monotonic and a cyclic load coexist, cyclic type will be considered
    // bool current_load_type = false;
    bool new_model_part = false;
    // bool break_condition = false;

    const bool has_constraints_list = mThisParameters["fatigue"].Has("constraints_process_list");
    if (!has_constraints_list) {
        KRATOS_ERROR << "Using the advance in time strategy without using the constraints_process_list" << std::endl;
    }

    if (has_constraints_list) {
        
        const int load_increments_per_cycle = mThisParameters["fatigue"]["load_increments_per_cycle"].GetInt();
        process_info[LOAD_INCREMENTS_PER_CYCLE] = load_increments_per_cycle;
              
        std::vector<std::string> constraints_list = mThisParameters["fatigue"]["constraints_process_list"].GetStringArray();
        //Loop on the cyclic constraints list
        for (unsigned int i = 0; i < constraints_list.size(); i++) {
            for (unsigned int j = 0; j < mThisParameters["processes"]["constraints_process_list"].size(); j++) {
                std::string model_part_name = mThisParameters["processes"]["constraints_process_list"][j]["Parameters"]["model_part_name"].GetString();
                double model_part_start_time = mThisParameters["processes"]["constraints_process_list"][j]["Parameters"]["interval"][0].GetDouble();
                double model_part_end_time = mThisParameters["processes"]["constraints_process_list"][j]["Parameters"]["interval"][1].GetDouble();
                if (constraints_list[i] == model_part_name && time >= model_part_start_time && time <= model_part_end_time) {
                    // current_load_type = true;
                    new_model_part = false; //This is done just in case a new monotonic load coexists with a cyclic load.
                                            //Then the thresholds used as reference should not be updated. This needs to be checked.

                    //Checking if this is the first step of a new model part
                    double model_part_start_time = mThisParameters["processes"]["constraints_process_list"][j]["Parameters"]["interval"][0].GetDouble();
                    if (time - delta_time <= model_part_start_time) {
                        new_model_part = true;
                        double new_model_part_start_time = mThisParameters["processes"]["constraints_process_list"][j]["Parameters"]["interval"][0].GetDouble();;
                        process_info[NEW_MODEL_PART_START_TIME] = new_model_part_start_time;
                    }
                }
            }
            // if (current_load_type) {
            //     break;
            // }
        }
    }

    // process_info[CURRENT_LOAD_TYPE] = current_load_type;
    process_info[NEW_MODEL_PART] = new_model_part;
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
    // std::vector<double> s_th;
    // std::vector<double> max_stress;

    double accumulated_max_stress_rel_error = 0.0;
    double accumulated_rev_factor_rel_error = 0.0;
    // bool fatigue_in_course = false;
    const double max_stress_relative_tolerance = mThisParameters["fatigue"]["max_stress_relative_tolerance"].GetDouble();
    const double rev_factor_relative_tolerance = mThisParameters["fatigue"]["rev_factor_relative_tolerance"].GetDouble();
    const double max_stress_relative_tolerance_damage = mThisParameters["fatigue"]["max_stress_relative_tolerance_damage"].GetDouble();
    const double rev_factor_relative_tolerance_damage = mThisParameters["fatigue"]["rev_factor_relative_tolerance_damage"].GetDouble();

    KRATOS_ERROR_IF(mrModelPart.NumberOfElements() == 0) << "The number of elements in the domain is zero. The process can not be applied."<< std::endl;

    for (auto& r_elem : mrModelPart.Elements()) {   //This loop is done for all the integration points even if the cycle has not changed
                                                    //This loop is done for all the integration points even if the cycle has not changed
                                                    //in order to guarantee the stable condition in the WHOLE model (when Smax > Sth)
        unsigned int number_of_ip = r_elem.GetGeometry().IntegrationPoints(r_elem.GetIntegrationMethod()).size();
        r_elem.CalculateOnIntegrationPoints(MAX_STRESS_RELATIVE_ERROR, max_stress_rel_error, process_info);
        r_elem.CalculateOnIntegrationPoints(REVERSION_FACTOR_RELATIVE_ERROR, rev_factor_rel_error, process_info);
        // r_elem.CalculateOnIntegrationPoints(THRESHOLD_STRESS, s_th, process_info);
        // r_elem.CalculateOnIntegrationPoints(MAX_STRESS, max_stress, process_info);

        std::vector<ConstitutiveLaw::Pointer> constitutive_law_vector(number_of_ip);
        r_elem.CalculateOnIntegrationPoints(CONSTITUTIVE_LAW,constitutive_law_vector,process_info);
                
        const bool is_fatigue = constitutive_law_vector[0]->Has(CYCLE_INDICATOR);
            
        if (is_fatigue) {  
            for (unsigned int i = 0; i < number_of_ip; i++) {
                accumulated_max_stress_rel_error += max_stress_rel_error[i];
                accumulated_rev_factor_rel_error += rev_factor_rel_error[i];
            }
        }
    }

    if ((!DamageIndicator && accumulated_max_stress_rel_error < max_stress_relative_tolerance && accumulated_rev_factor_rel_error < rev_factor_relative_tolerance) || (DamageIndicator && accumulated_max_stress_rel_error < max_stress_relative_tolerance_damage && accumulated_rev_factor_rel_error < rev_factor_relative_tolerance_damage)) {
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
    // std::vector<double> s_th;
    // std::vector<double> max_stress;
    double min_time_increment;
    double time = process_info[TIME];
    bool current_constraints_process_list_detected = false;

    double user_advancing_time = mThisParameters["fatigue"]["advancing_strategy_time"].GetDouble();
    double user_advancing_cycles = mThisParameters["fatigue"]["advancing_strategy_cycles"].GetDouble();
    double user_advancing_time_damage = mThisParameters["fatigue"]["advancing_strategy_damage"].GetDouble();
    double user_advancing_cycles_damage = mThisParameters["fatigue"]["advancing_strategy_cycles_damage"].GetDouble();
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
    if (!process_info[NO_LINEARITY_ACTIVATION]) {
        min_time_increment = std::min(user_advancing_time, model_part_advancing_time);
    }else{
        min_time_increment = std::min(user_advancing_time_damage, model_part_advancing_time);
    }

    // KRATOS_ERROR_IF(mrModelPart.NumberOfElements() == 0) << "The number of elements in the domain is zero. The process can not be applied."<< std::endl;

        if (!process_info[NO_LINEARITY_ACTIVATION]) {
            for (auto &r_elem : mrModelPart.Elements()){
                unsigned int number_of_ip = r_elem.GetGeometry().IntegrationPoints(r_elem.GetIntegrationMethod()).size();
                r_elem.CalculateOnIntegrationPoints(CYCLES_TO_FAILURE, cycles_to_failure_element, process_info);
                r_elem.CalculateOnIntegrationPoints(LOCAL_NUMBER_OF_CYCLES, local_number_of_cycles, process_info);
                r_elem.CalculateOnIntegrationPoints(CYCLE_PERIOD, period, process_info);
                // r_elem.CalculateOnIntegrationPoints(THRESHOLD_STRESS, s_th, process_info);
                // r_elem.CalculateOnIntegrationPoints(MAX_STRESS, max_stress, process_info);

                std::vector<ConstitutiveLaw::Pointer> constitutive_law_vector(number_of_ip);
                r_elem.CalculateOnIntegrationPoints(CONSTITUTIVE_LAW,constitutive_law_vector,process_info);
                        
                const bool is_fatigue = constitutive_law_vector[0]->Has(CYCLE_INDICATOR);

                if(is_fatigue){
                    for (unsigned int i = 0; i < number_of_ip; i++)
                    {
                        // if (max_stress[i] > s_th[i])
                        // { // This is used to guarantee that only IP in fatigue conditions are taken into account
                        double Nf_conversion_to_time = (cycles_to_failure_element[i] - local_number_of_cycles[i]) * period[i];
                        double user_advancing_cycles_conversion_to_time = user_advancing_cycles * period[i];
                        if (Nf_conversion_to_time < min_time_increment)
                        {   
                            min_time_increment = Nf_conversion_to_time;
                        }
                        if (user_advancing_cycles_conversion_to_time < min_time_increment)
                        {   
                            min_time_increment = user_advancing_cycles_conversion_to_time;
                        }
                    }
                }
            }
        }else{
            for (auto &r_elem : mrModelPart.Elements()){
                unsigned int number_of_ip = r_elem.GetGeometry().IntegrationPoints(r_elem.GetIntegrationMethod()).size();
                r_elem.CalculateOnIntegrationPoints(CYCLE_PERIOD, period, process_info);
                // r_elem.CalculateOnIntegrationPoints(THRESHOLD_STRESS, s_th, process_info);
                // r_elem.CalculateOnIntegrationPoints(MAX_STRESS, max_stress, process_info);

                std::vector<ConstitutiveLaw::Pointer> constitutive_law_vector(number_of_ip);
                r_elem.CalculateOnIntegrationPoints(CONSTITUTIVE_LAW,constitutive_law_vector,process_info);
                        
                const bool is_fatigue = constitutive_law_vector[0]->Has(CYCLE_INDICATOR);

                if(is_fatigue){
                    for (unsigned int i = 0; i < number_of_ip; i++)
                    {
                        // if (max_stress[i] > s_th[i])
                        // { // This is used to guarantee that only IP in fatigue conditions are taken into account
                        double user_advancing_cycles_conversion_to_time = user_advancing_cycles_damage * period[i];
                        if (user_advancing_cycles_conversion_to_time < min_time_increment)
                        {
                            min_time_increment = user_advancing_cycles_conversion_to_time;
                        }
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
        double time_increment;
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