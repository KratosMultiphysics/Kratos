// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: kratos/license.txt
//
//  Main authors:    Luis Antonio Goncalves Junior
//                   Alejandro Cornejo
//

// System includes

// External includes



// Project includes
#include "includes/model_part.h"
#include "custom_processes/set_automated_initial_damage_process.h"
#include "utilities/parallel_utilities.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "constitutive_laws_application_variables.h"

namespace Kratos
{
SetAutomatedInitialDamageProcess::SetAutomatedInitialDamageProcess(
    ModelPart& rThisModelPart,
    Parameters ThisParameters
    ):mrThisModelPart(rThisModelPart),
      mThisParameters(ThisParameters)
{
    mThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());
}

/***********************************************************************************/
/***********************************************************************************/

void SetAutomatedInitialDamageProcess::ExecuteInitializeSolutionStep()
{
    
    if (mrThisModelPart.GetProcessInfo()[STEP] == 1){

        KRATOS_TRY

        const array_1d<double, 3> hole_generatrix_axis = mThisParameters["hole_generatrix_axis"].GetVector();
        KRATOS_ERROR_IF(MathUtils<double>::Norm3(hole_generatrix_axis) < machine_tolerance) << "The hole generatrix axis has norm zero" << std::endl;
        
        const array_1d<double, 3> hole_generatrix_point = mThisParameters["hole_generatrix_point"].GetVector();
        auto& process_info = mrThisModelPart.GetProcessInfo();  

        array_1d<double, 3> normalized_generatrix_vector = hole_generatrix_axis;
        ConstitutiveLawUtilities<3>::CheckAndNormalizeVector<array_1d<double,3>>(normalized_generatrix_vector);

        const double hole_radius_offset = mThisParameters["hole_radius_offset"].GetDouble();

        const int table_id = mThisParameters["table_id"].GetInt();

        block_for_each(mrThisModelPart.Elements(), [&](auto& rElement) {

            const array_1d<double, 3>& r_element_centroid = rElement.GetGeometry().Center();
            
            array_1d<double, 3> relative_position_vector;
            relative_position_vector = r_element_centroid - hole_generatrix_point;

            const double vector_scaler = MathUtils<double>::Dot3(relative_position_vector, normalized_generatrix_vector);

            const array_1d<double, 3> intersection_point = hole_generatrix_point + vector_scaler * normalized_generatrix_vector;

            const array_1d<double, 3> radial_position_vector = r_element_centroid - intersection_point;

            double centroid_relative_distance = MathUtils<double>::Norm3(radial_position_vector) - hole_radius_offset;

            if (centroid_relative_distance < 0.0){
                if (std::abs(centroid_relative_distance) <= tolerance) {
                    centroid_relative_distance = 0.0;
                } else {
                    KRATOS_ERROR << "The relative centroid distance may not be negative. Check the hole radius offset and the thickness of element " << rElement.Id() << std::endl;
                }
            }

            double initial_damage_value = mrThisModelPart.GetTable(table_id).GetValue(centroid_relative_distance);

            if (initial_damage_value < 0.0 ){
                initial_damage_value = 0.0;

            }  else if (initial_damage_value >= 1.0 ){
                initial_damage_value = 0.999;
            }

            const unsigned int number_of_integration_points = rElement.GetGeometry().IntegrationPoints(rElement.GetIntegrationMethod()).size();

            std::vector<double> r_threshold(number_of_integration_points);           
            std::vector<double> initial_damage(number_of_integration_points);       

            rElement.CalculateOnIntegrationPoints(THRESHOLD, r_threshold, process_info);

            for ( IndexType point_number = 0; point_number < number_of_integration_points; ++point_number ) {
                initial_damage[point_number] = initial_damage_value;
                r_threshold[point_number] *= (1-initial_damage_value);
            }

            rElement.SetValuesOnIntegrationPoints(DAMAGE, initial_damage, process_info);
            rElement.SetValuesOnIntegrationPoints(THRESHOLD, r_threshold, process_info);
            
        });
        
        KRATOS_CATCH("")

    }
}

/***********************************************************************************/
/***********************************************************************************/

const Parameters SetAutomatedInitialDamageProcess::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"(
    {
        "hole_generatrix_axis"     : [0.0,0.0,1.0],
        "hole_generatrix_point"    : [0.0,0.0,0.0],
        "hole_radius_offset"       : 0.0,
        "table_id"                 : 0
    })");

    return default_parameters;
}

} // namespace Kratos.