// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Luis Antonio Goncalves Junior
//                   Alejandro Cornejo
//

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "custom_processes/set_automated_initial_variable_process_cylindrical.h"
#include "utilities/parallel_utilities.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
SetAutomatedInitialVariableProcessCylindrical::SetAutomatedInitialVariableProcessCylindrical(
    ModelPart& rThisModelPart,
    Parameters ThisParameters
    ):mrThisModelPart(rThisModelPart),
      mThisParameters(ThisParameters)
{
    mThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());
}

/***********************************************************************************/
/***********************************************************************************/

void SetAutomatedInitialVariableProcessCylindrical::ExecuteInitialize()
{
    
    KRATOS_TRY

    const array_1d<double, 3> hole_generatrix_axis = mThisParameters["hole_generatrix_axis"].GetVector();
    KRATOS_ERROR_IF(MathUtils<double>::Norm3(hole_generatrix_axis) < machine_tolerance) << "The hole generatrix axis has norm zero" << std::endl;
    
    const array_1d<double, 3> hole_generatrix_point = mThisParameters["hole_generatrix_point"].GetVector();  

    array_1d<double, 3> normalized_generatrix_vector = hole_generatrix_axis;
    ConstitutiveLawUtilities<3>::CheckAndNormalizeVector<array_1d<double,3>>(normalized_generatrix_vector);

    const double hole_radius_offset = mThisParameters["hole_radius_offset"].GetDouble();

    Vector table_id_vector = mThisParameters["table_id_vector"].GetVector();
    const Variable<Vector>& r_variable = KratosComponents<Variable<Vector>>::Get(mThisParameters["variable_name"].GetString());

    block_for_each(mrThisModelPart.Elements(), [&](auto& rElement) {

        const array_1d<double, 3>& r_element_centroid = rElement.GetGeometry().Center();
        
        array_1d<double, 3> relative_position_vector;
        relative_position_vector = r_element_centroid - hole_generatrix_point;

        double vector_scaler = MathUtils<double>::Dot3(relative_position_vector, normalized_generatrix_vector);
        
        const array_1d<double, 3> intersection_point = hole_generatrix_point + vector_scaler * normalized_generatrix_vector;

        const array_1d<double, 3> radial_position_vector = r_element_centroid - intersection_point;

        double centroid_relative_distance = MathUtils<double>::Norm3(radial_position_vector) - hole_radius_offset;

        if (centroid_relative_distance < tolerance){
            // centroid_relative_distance = 2.05E-06; //CP800
            centroid_relative_distance = 3.25E-06; //CP980
        }

        // if (centroid_relative_distance < tolerance){
        //     centroid_relative_distance = 2.4878275E-06;
        // }
        
        // if (centroid_relative_distance < 0.0){
        //     if (std::abs(centroid_relative_distance) <= tolerance) {
        //         centroid_relative_distance = 0.0;
        //     } else {
        //         KRATOS_ERROR << "The relative centroid distance may not be negative. Check the hole radius offset and the thickness of element " << rElement.Id() << std::endl;
        //     }
        // }
        
        array_1d<double, 6> initial_variable_vector;
        noalias(initial_variable_vector) = ZeroVector(6); 
 
        int table_id;
        int varible_id;

        for (IndexType count = 0; count < table_id_vector.size(); ++count) {
                table_id = table_id_vector[count];
                varible_id = table_id % 10;
                initial_variable_vector[varible_id] = mrThisModelPart.GetTable(table_id).GetValue(centroid_relative_distance);  
        }
    
        rElement.SetValue(r_variable, initial_variable_vector);

    });
    
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

const Parameters SetAutomatedInitialVariableProcessCylindrical::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"(
    {
        "variable_name"            : "SPECIFY_VARIABLE_NAME",
        "hole_generatrix_axis"     : [0.0,0.0,1.0],
        "hole_generatrix_point"    : [0.0,0.0,0.0],
        "hole_radius_offset"       : 0.0,
        "table_id_vector"          : [10,11,12,13,14,15]
    })");

    return default_parameters;
}

} // namespace Kratos.