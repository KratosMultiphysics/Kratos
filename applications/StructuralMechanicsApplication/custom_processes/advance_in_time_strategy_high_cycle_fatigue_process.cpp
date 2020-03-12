//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Sergio Jimenez Reyes
//

#include "custom_processes/advance_in_time_strategy_high_cycle_fatigue_process.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

// Constructor
	AdvanceInTimeStrategyHighCycleFatigueProcess::AdvanceInTimeStrategyHighCycleFatigueProcess(
		ModelPart &rModelPart,
		Parameters ThisParameters) : mrModelPart(rModelPart), mThisParameters(ThisParameters)
{
}

/***********************************************************************************/
/***********************************************************************************/

void AdvanceInTimeStrategyHighCycleFatigueProcess::Execute()
{
    KRATOS_WATCH("uno")
    auto& process_info = mrModelPart.GetProcessInfo();
    bool cycle_found = false;
    std::vector<bool> cycle_identificator;
    std::vector<int> cycles_after_advance_strategy;
    std::vector<double> damage;
    bool damage_indicator = false;
    process_info[ADVANCE_STRATEGY_APPLIED] = false;
    bool cycles_from_last_advance = false;

    for (auto& r_elem : mrModelPart.Elements()) {
        r_elem.GetValueOnIntegrationPoints(DAMAGE, damage, process_info);
        r_elem.GetValueOnIntegrationPoints(CYCLES_AFTER_ADVANCE_STRATEGY, cycles_after_advance_strategy, process_info);
		const int number_of_ip = r_elem.GetGeometry().IntegrationPoints(r_elem.GetIntegrationMethod()).size();
        for (unsigned int i = 0; i < number_of_ip; i++){
                if (damage[i] > 0.0){
                    damage_indicator = true;
                    process_info[DAMAGE_ACTIVATION] = true;
                }
                if (cycles_after_advance_strategy[i] > 1){
                    cycles_from_last_advance = true;
                }
        }
    }

    this->CyclePeriodPerIntegrationPoint(cycle_found);  //This method detects if a cycle has finished somewhere in the model and 
                                                        //computes the time period of the cycle that has just finished.
    
    if (cycle_found) {  //If a cycle has finished then it is possible to apply the advancing strategy
        bool advancing_strategy = false;
        this->StableConditionForAdvancingStrategy(advancing_strategy, damage_indicator);  //Check if the conditions are optimal to apply the advancing strategy in 
                                                                        //terms of max stress and reversion factor variation.  
        double increment = 0.0;
        // if (advancing_strategy & !damage_indicator & cycles_from_last_advance) {
        if (advancing_strategy) {
            if (!damage_indicator) {
                this->TimeIncrement(increment);
                if(increment > 0){
                    this->TimeAndCyclesUpdate(increment);
                    process_info[ADVANCE_STRATEGY_APPLIED] = true;
                }
            // } else {
            //     increment = 1000.0*12.0;
            //         this->TimeAndCyclesUpdate(increment);   
            //         process_info[ADVANCE_STRATEGY_APPLIED] = true;             
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void AdvanceInTimeStrategyHighCycleFatigueProcess::CyclePeriodPerIntegrationPoint(bool& rCycleFound) 
{
    auto& process_info = mrModelPart.GetProcessInfo();
    std::vector<bool> cycle_identificator;
    std::vector<double> previous_cycle_time;    //time when the previous cycle finished. It is used to obtain the new period for the current cycle
    std::vector<double> period;
    double time = process_info[TIME];
    
    for (auto& r_elem : mrModelPart.Elements()) {

        r_elem.GetValueOnIntegrationPoints(CYCLE_INDICATOR, cycle_identificator, process_info);
        r_elem.GetValueOnIntegrationPoints(PREVIOUS_CYCLE, previous_cycle_time, process_info);
        r_elem.GetValueOnIntegrationPoints(CYCLE_PERIOD, period, process_info);

		const int number_of_ip = r_elem.GetGeometry().IntegrationPoints(r_elem.GetIntegrationMethod()).size();
        
        // Check if the HCF CL is being used. Otherwise there is no reason for entering into the advance in time strategy
        std::vector<ConstitutiveLaw::Pointer> constitutive_law_vector(number_of_ip);
        r_elem.GetValueOnIntegrationPoints(CONSTITUTIVE_LAW,constitutive_law_vector,process_info);

        const bool is_fatigue = constitutive_law_vector[0]->Has(CYCLE_INDICATOR);
        
        if (is_fatigue){
            for (unsigned int i = 0; i < number_of_ip; i++){
                    if (cycle_identificator[i]){
                        period[i] = time - previous_cycle_time[i];
                        previous_cycle_time[i] = time;
                        rCycleFound = true;
                    }
            }
            r_elem.SetValueOnIntegrationPoints(PREVIOUS_CYCLE, previous_cycle_time, process_info);
            r_elem.SetValueOnIntegrationPoints(CYCLE_PERIOD, period, process_info);
        }    
    }
}

/***********************************************************************************/
/***********************************************************************************/

void AdvanceInTimeStrategyHighCycleFatigueProcess::StableConditionForAdvancingStrategy(bool& rAdvancingStrategy, bool DamageIndicator) 
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
    for (auto& r_elem : mrModelPart.Elements()) {   //Este loop se hace por todos los elementos y PI independientemente que haya habido o no cambio 
                                                    //de ciclo, porque se debe garantizar condición de estabilidad en TODO el modelo (siempre que 
                                                    //Smax > Sth)
        
        r_elem.Id();
        r_elem.GetValueOnIntegrationPoints(MAX_STRESS_RELATIVE_ERROR, max_stress_rel_error, process_info);
        r_elem.GetValueOnIntegrationPoints(REVERSION_FACTOR_RELATIVE_ERROR, rev_factor_rel_error, process_info);
        r_elem.GetValueOnIntegrationPoints(THRESHOLD_STRESS, s_th, process_info);
        r_elem.GetValueOnIntegrationPoints(MAX_STRESS, max_stress, process_info);        

        // for (unsigned int i=0; i < max_stress_rel_error.size(); i++) {
        const int number_of_ip = r_elem.GetGeometry().IntegrationPoints(r_elem.GetIntegrationMethod()).size();
        for (unsigned int i = 0; i < number_of_ip; i++){
            if (max_stress[i] > s_th[i]) {
                fatigue_in_course = true;
                acumulated_max_stress_rel_error += max_stress_rel_error[i];
                acumulated_rev_factor_rel_error += rev_factor_rel_error[i];
            }
        }
    }
    KRATOS_WATCH(acumulated_max_stress_rel_error)
    KRATOS_WATCH(acumulated_rev_factor_rel_error)
    if ((acumulated_max_stress_rel_error < 1e-4 && acumulated_rev_factor_rel_error < 1e-4 && fatigue_in_course) || (DamageIndicator && acumulated_max_stress_rel_error < 1e-2 && acumulated_rev_factor_rel_error < 1e-2 && fatigue_in_course)) {
        rAdvancingStrategy = true;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void AdvanceInTimeStrategyHighCycleFatigueProcess::TimeIncrement(double& rIncrement) 
{
    auto& process_info = mrModelPart.GetProcessInfo();
    std::vector<double>  cycles_to_failure_element;
    std::vector<int>  local_number_of_cycles;
    std::vector<double> period;
    std::vector<double> s_th;
    std::vector<double> max_stress;
    double min_time_increment;
    double time = process_info[TIME];
    // const double period_json = mThisParameters["fatigue"]["period"].GetDouble();
    bool current_constraints_process_list_detected = false;

    double user_avancing_time = mThisParameters["fatigue"]["advancing_strategy_time"].GetDouble();
    double user_avancing_cycles = mThisParameters["fatigue"]["advancing_strategy_cycles"].GetDouble();
    std::vector<std::string> constraints_list = mThisParameters["fatigue"]["constraints_process_list"].GetStringArray();
    
    double model_part_final_time = mThisParameters["problem_data"]["end_time"].GetDouble();
    for (unsigned int i = 0; i < constraints_list.size(); i++){
        for (unsigned int j = 0; j < mThisParameters["processes"]["constraints_process_list"].size(); j++){
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
    for (auto& r_elem : mrModelPart.Elements()) {
        r_elem.GetValueOnIntegrationPoints(CYCLES_TO_FAILURE, cycles_to_failure_element, process_info); 
        r_elem.GetValueOnIntegrationPoints(LOCAL_NUMBER_OF_CYCLES, local_number_of_cycles, process_info); 
        r_elem.GetValueOnIntegrationPoints(CYCLE_PERIOD, period, process_info);
        r_elem.GetValueOnIntegrationPoints(THRESHOLD_STRESS, s_th, process_info);
        r_elem.GetValueOnIntegrationPoints(MAX_STRESS, max_stress, process_info);
		const int number_of_ip = r_elem.GetGeometry().IntegrationPoints(r_elem.GetIntegrationMethod()).size();
        for (unsigned int i = 0; i < number_of_ip; i++){
            if (max_stress[i] > s_th[i]) {  //This is used to guarantee that only IP in fatigue conditions are taken into account
                double Nf_conversion_to_time = (cycles_to_failure_element[i] - local_number_of_cycles[i]) * period[i];
                // double Nf_conversion_to_time = (cycles_to_failure_element[i] - local_number_of_cycles[i]) * period_json;
                double user_avancing_cycles_conversion_to_time = user_avancing_cycles * period[i];
                // double user_avancing_cycles_conversion_to_time = user_avancing_cycles * period_json;
                if (Nf_conversion_to_time < min_time_increment){
                    min_time_increment = Nf_conversion_to_time;
                }
                if (user_avancing_cycles_conversion_to_time < min_time_increment){
                    min_time_increment = user_avancing_cycles_conversion_to_time;
                }
            }
        }
    }
	rIncrement = min_time_increment;
}

/***********************************************************************************/
/***********************************************************************************/

void AdvanceInTimeStrategyHighCycleFatigueProcess::TimeAndCyclesUpdate(double Increment) 
{
    auto& r_process_info = mrModelPart.GetProcessInfo();
    std::vector<int>  local_number_of_cycles;
    std::vector<int>  global_number_of_cycles;
    std::vector<double> period;
    double time_increment;
    std::vector<double> previous_cycle_time;    //time when the previous cycle finished. It is used to obtain the new period for the current cycle

    for (auto& r_elem : mrModelPart.Elements()) {
        r_elem.GetValueOnIntegrationPoints(LOCAL_NUMBER_OF_CYCLES, local_number_of_cycles, r_process_info); 
        r_elem.GetValueOnIntegrationPoints(NUMBER_OF_CYCLES, global_number_of_cycles, r_process_info); 
        r_elem.GetValueOnIntegrationPoints(CYCLE_PERIOD, period, r_process_info);
        r_elem.GetValueOnIntegrationPoints(PREVIOUS_CYCLE, previous_cycle_time, r_process_info);

		const int number_of_ip = r_elem.GetGeometry().IntegrationPoints(r_elem.GetIntegrationMethod()).size();
        for (unsigned int i = 0; i < number_of_ip; i++){
            unsigned int local_cycles_increment;
            if (period[i] == 0.0) {
                local_cycles_increment = 0;
            } else {
                local_cycles_increment = std::trunc(Increment / period[i]);
                time_increment = std::trunc(Increment / period[i]) * period[i];
                previous_cycle_time[i] += time_increment;
            }
            // unsigned int local_cycles_increment = std::trunc(Increment / period_json);
            local_number_of_cycles[i] += local_cycles_increment;
            global_number_of_cycles[i] += local_cycles_increment;
        }
        r_elem.SetValueOnIntegrationPoints(LOCAL_NUMBER_OF_CYCLES, local_number_of_cycles, r_process_info);
        r_elem.SetValueOnIntegrationPoints(NUMBER_OF_CYCLES, global_number_of_cycles, r_process_info);
        r_elem.SetValueOnIntegrationPoints(PREVIOUS_CYCLE, previous_cycle_time, r_process_info);
    }
    r_process_info[TIME_INCREMENT] = time_increment - 0.0;

}



} // namespace Kratos
