// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Luis Antônio Gonçalves Junior
//                   Alejandro Cornejo
//

#include "includes/model_part.h"
#include "custom_processes/set_automated_initial_variable_process.h"
#include "utilities/parallel_utilities.h"
#include "utilities/math_utils.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
SetAutomatedInitialVariableProcess::SetAutomatedInitialVariableProcess(
    ModelPart& rThisModelPart,
    Parameters ThisParameters
    ):mrThisModelPart(rThisModelPart),
      mThisParameters(ThisParameters)
{
    mThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());
}

/***********************************************************************************/
/***********************************************************************************/

void SetAutomatedInitialVariableProcess::ExecuteInitialize()
{

    KRATOS_TRY

    array_1d<double, 3> hole_generatrix_axis = mThisParameters["hole_generatrix_axis"].GetVector();
    KRATOS_ERROR_IF(MathUtils<double>::Norm3(hole_generatrix_axis) < std::numeric_limits<double>::epsilon()) << "The hole generatrix axis has norm zero" << std::endl;
    
    array_1d<double, 3> hole_generatrix_point = mThisParameters["hole_generatrix_point"].GetVector();  

    array_1d<double, 3> normalized_generatrix_vector;   
    normalized_generatrix_vector[0] = (hole_generatrix_axis[0] * hole_generatrix_axis[0]) / std::sqrt(hole_generatrix_axis[0] * hole_generatrix_axis[0] + hole_generatrix_axis[1] * hole_generatrix_axis[1] + hole_generatrix_axis[2] * hole_generatrix_axis[2]);
    normalized_generatrix_vector[1] = (hole_generatrix_axis[1] * hole_generatrix_axis[1]) / std::sqrt(hole_generatrix_axis[0] * hole_generatrix_axis[0] + hole_generatrix_axis[1] * hole_generatrix_axis[1] + hole_generatrix_axis[2] * hole_generatrix_axis[2]);
    normalized_generatrix_vector[2] = (hole_generatrix_axis[2] * hole_generatrix_axis[2]) / std::sqrt(hole_generatrix_axis[0] * hole_generatrix_axis[0] + hole_generatrix_axis[1] * hole_generatrix_axis[1] + hole_generatrix_axis[2] * hole_generatrix_axis[2]);

    const double hole_radius_offset = mThisParameters["hole_radius_offset"].GetDouble();
    const double tol = 1.0e-6;

    block_for_each(mrThisModelPart.Elements(), [&](Element &rElement) {

        array_1d<double, 3> element_centroid = rElement.GetGeometry().Center();
        
        array_1d<double, 3> relative_position_vector;
        relative_position_vector[0] = element_centroid[0] - hole_generatrix_point[0];
        relative_position_vector[1] = element_centroid[1] - hole_generatrix_point[1];
        relative_position_vector[2] = element_centroid[2] - hole_generatrix_point[2];

        double vector_scaler = relative_position_vector[0] * normalized_generatrix_vector[0] + relative_position_vector[1] * normalized_generatrix_vector[1] + relative_position_vector[2] * normalized_generatrix_vector[2];
        
        array_1d<double, 3> intersection_point = hole_generatrix_point + vector_scaler * normalized_generatrix_vector;

        array_1d<double, 3> radial_position_vector = element_centroid - intersection_point;

        double centroid_relative_distance = std::sqrt(radial_position_vector[0] * radial_position_vector[0] + radial_position_vector[1] * radial_position_vector[1] + radial_position_vector[2] * radial_position_vector[2]) - hole_radius_offset;

        if (centroid_relative_distance < 0.0 && abs(centroid_relative_distance) <= tol) {
            centroid_relative_distance = 0.0;
        }
        else if(centroid_relative_distance < 0.0 && abs(centroid_relative_distance) >= tol){
            int ElemId = rElement.Id();
            KRATOS_ERROR << "Thickness of element " << ElemId << " is too small." << std::endl;
        }

        int TableFirstId = mThisParameters["initial_variable_table"]["table_id"].GetInt() - mrThisModelPart.Tables().size()+1;
        
        int count = 0;
        array_1d<double, 6> initial_variable_vector = ZeroVector(); 

        for (IndexType TableId = TableFirstId; TableId < TableFirstId + mrThisModelPart.Tables().size(); ++TableId) {
                initial_variable_vector[count] = mrThisModelPart.GetTable(TableId).GetValue(centroid_relative_distance);  
                count += 1;
        }
    
        const Variable<Vector>& r_variable = KratosComponents<Variable<Vector>>::Get(mThisParameters["variable_name"].GetString());
        rElement.SetValue(r_variable, initial_variable_vector);

    });
    
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

const Parameters SetAutomatedInitialVariableProcess::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"(
    {
        "help"                     : "This automates the application of initial strain/stress variables using csv tables",
        "model_part_name"          : "please_specify_model_part_name",
        "variable_name"            : "SPECIFY_VARIABLE_NAME",
        "hole_generatrix_axis"     : [0.0,0.0,1.0],
        "hole_generatrix_point"    : [0.0,0.0,0.0],
        "hole_radius_offset"       : 0.0,
        "last_layer"               : false,
        "initial_variable_table"     : {
                    "name"             : "csv_table",
                    "filename"         : "sample.csv",
                    "delimiter"        : ",",
                    "skiprows"         : 1,
                    "first_column_id"  : 0,
                    "second_column_id" : 1,
                    "table_id"         : -1,
                    "na_replace"       : 0.0
                }
    })");

    return default_parameters;
}

} // namespace Kratos.