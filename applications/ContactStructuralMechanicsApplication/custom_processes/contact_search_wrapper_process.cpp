// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "custom_processes/contact_search_wrapper_process.h"
#include "custom_processes/simple_contact_search_process.h"
#include "custom_processes/advanced_contact_search_process.h"

namespace Kratos
{
ContactSearchWrapperProcess::ContactSearchWrapperProcess(
    ModelPart& rMainModelPart,
    Parameters ThisParameters
    )
{
    // The default parameters
    Parameters default_parameters = GetDefaultParameters();
    ThisParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

    // The dimensions
    const SizeType dimension = rMainModelPart.GetProcessInfo()[DOMAIN_SIZE];
    const bool predefined_master_slave = ThisParameters["predefined_master_slave"].GetBool();
    SizeType size_1 = 0;
    SizeType size_2 = 0;
    for (auto& r_cond : rMainModelPart.Conditions()) {

        if (predefined_master_slave) {
            if (r_cond.Is(MASTER))
                size_1 = r_cond.GetGeometry().size();
            if (r_cond.Is(SLAVE))
                size_2 = r_cond.GetGeometry().size();
        } else {
            size_1 = r_cond.GetGeometry().size();
            size_2 = size_1;
        }

        // Once both defined we break
        if (size_1 > 0 && size_2 > 0)
            break;
    }

    // Simple search bool
    const bool simple_search = ThisParameters["simple_search"].GetBool();
    ThisParameters.RemoveValue("simple_search"); // Removing to avoid problems

    // Creating the mapper
    if (dimension == 2) {
        // 2D
        if (simple_search) {
            mpContactProcess = Kratos::make_shared<SimpleContactSearchProcess<2, 2>>(rMainModelPart, ThisParameters);
        } else {
            mpContactProcess = Kratos::make_shared<AdvancedContactSearchProcess<2, 2>>(rMainModelPart, ThisParameters);
        }
    } else {
        // 3D
        if (size_1 == 3 && size_2 == 3) {
            if (simple_search) {
                mpContactProcess = Kratos::make_shared<SimpleContactSearchProcess<3, 3>>(rMainModelPart, ThisParameters);
            } else {
                mpContactProcess = Kratos::make_shared<AdvancedContactSearchProcess<3, 3>>(rMainModelPart, ThisParameters);
            }
        } else if (size_1 == 4 && size_2 == 4) {
            if (simple_search) {
                mpContactProcess = Kratos::make_shared<SimpleContactSearchProcess<3, 4>>(rMainModelPart, ThisParameters);
            } else {
                mpContactProcess = Kratos::make_shared<AdvancedContactSearchProcess<3, 4>>(rMainModelPart, ThisParameters);
            }
        } else if (size_1 == 4 && size_2 == 3) {
            if (simple_search) {
                mpContactProcess = Kratos::make_shared<SimpleContactSearchProcess<3, 3, 4>>(rMainModelPart, ThisParameters);
            } else {
                mpContactProcess = Kratos::make_shared<AdvancedContactSearchProcess<3, 3, 4>>(rMainModelPart, ThisParameters);
            }
        } else if (size_1 == 3 && size_2 == 4) {
            if (simple_search) {
                mpContactProcess = Kratos::make_shared<SimpleContactSearchProcess<3, 4, 3>>(rMainModelPart, ThisParameters);
            } else {
                mpContactProcess = Kratos::make_shared<AdvancedContactSearchProcess<3, 4, 3>>(rMainModelPart, ThisParameters);
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

Parameters ContactSearchWrapperProcess::GetDefaultParameters()
{
    Parameters default_parameters = Parameters(R"(
    {
        "simple_search"                        : false,
        "allocation_size"                      : 1000,
        "bucket_size"                          : 4,
        "search_factor"                        : 3.5,
        "type_search"                          : "InRadius",
        "check_gap"                            : "MappingCheck",
        "condition_name"                       : "",
        "final_string"                         : "",
        "inverted_search"                      : false,
        "dynamic_search"                       : false,
        "static_check_movement"                : false,
        "predefined_master_slave"              : true,
        "id_name"                              : "",
        "normal_orientation_threshold"         : 1.0e-1,
        "consider_gap_threshold"               : false,
        "predict_correct_lagrange_multiplier"  : false,
        "pure_slip"                            : false,
        "debug_mode"                           : false,
        "octree_search_parameters" : {
            "bounding_box_factor"             : 0.1,
            "debug_obb"                       : false,
            "OBB_intersection_type"           : "SeparatingAxisTheorem",
            "build_from_bounding_box"         : true,
            "lower_bounding_box_coefficient"  : 0.0,
            "higher_bounding_box_coefficient" : 1.0
            }
    })" );

    return default_parameters;
}

}  // namespace Kratos.
