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
#include "custom_processes/set_automated_initial_variable_process_cartesian.h"
#include "utilities/parallel_utilities.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
SetAutomatedInitialVariableProcessCartesian::SetAutomatedInitialVariableProcessCartesian(
    ModelPart& rThisModelPart,
    Parameters ThisParameters
    ):mrThisModelPart(rThisModelPart),
      mThisParameters(ThisParameters)
{
    mThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());
}

/***********************************************************************************/
/***********************************************************************************/

void SetAutomatedInitialVariableProcessCartesian::ExecuteInitialize()
{
    
    KRATOS_TRY

    const array_1d<double, 3> surface_normal_vector = mThisParameters["surface_normal_vector"].GetVector();
    KRATOS_ERROR_IF(MathUtils<double>::Norm3(surface_normal_vector) < machine_tolerance) << "The surface normal vector has norm zero" << std::endl;
    
    const array_1d<double, 3> surface_reference_point = mThisParameters["surface_reference_point"].GetVector();

    array_1d<double, 3> normalized_normal_vector = surface_normal_vector;
    ConstitutiveLawUtilities<3>::CheckAndNormalizeVector<array_1d<double,3>>(normalized_normal_vector);

    double surface_reference_projection = MathUtils<double>::Dot3(surface_reference_point, normalized_normal_vector);

    Vector table_id_vector = mThisParameters["table_id_vector"].GetVector();
    const Variable<Vector>& r_variable = KratosComponents<Variable<Vector>>::Get(mThisParameters["variable_name"].GetString());

    block_for_each(mrThisModelPart.Elements(), [&](auto& rElement) {

        const array_1d<double, 3>& r_element_centroid = rElement.GetGeometry().Center();
        double r_element_centroid_projection = MathUtils<double>::Dot3(r_element_centroid, normalized_normal_vector);

        double centroid_relative_distance = surface_reference_projection - r_element_centroid_projection;

        if (centroid_relative_distance < 0.0){
           KRATOS_ERROR << "The relative centroid distance may not be negative. Check if the surface normal vector is an outward vector" << std::endl;
        }

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

const Parameters SetAutomatedInitialVariableProcessCartesian::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"(
    {
        "variable_name"            : "SPECIFY_VARIABLE_NAME",
        "surface_normal_vector"    : [0.0,0.0,1.0],
        "surface_reference_point"  : [0.0,0.0,0.0],
        "table_id_vector"          : [10,11,12,13,14,15]
    })");

    return default_parameters;
}

} // namespace Kratos.