//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

// System includes
//#include <string>
//#include <iostream>

// External includes

// Project includes
// #include "includes/define.h"
#include "distance_modification_process.h"

// #include "processes/process.h"
// #include "includes/cfd_variables.h"
// #include "processes/find_nodal_h_process.h"
// #include "utilities/openmp_utils.h"
// #include "includes/kratos_parameters.h"

// Application includes


namespace Kratos
{

/* Public functions *******************************************************/
DistanceModificationProcess::DistanceModificationProcess(
    ModelPart& rModelPart, 
    const bool CheckAtEachStep, 
    const bool NegElemDeactivation,
    const bool RecoverOriginalDistance)
    : Process(), mrModelPart(rModelPart) {

    mFactorCoeff = 2.0;
    mCheckAtEachStep = CheckAtEachStep;
    mNegElemDeactivation = NegElemDeactivation;
    mRecoverOriginalDistance = RecoverOriginalDistance;
}

DistanceModificationProcess::DistanceModificationProcess(
    ModelPart& rModelPart,
    Parameters& rParameters)
    : Process(), mrModelPart(rModelPart) {

    Parameters default_parameters( R"(
    {
        "mesh_id"                                : 0,
        "model_part_name"                        : "default_model_part_name",
        "distance_factor"                        : 2.0,
        "check_at_each_time_step"                : false,
        "deactivate_full_negative_elements"      : true,
        "recover_original_distance_at_each_step" : false
    }  )" );

    rParameters.ValidateAndAssignDefaults(default_parameters);

    mFactorCoeff = rParameters["distance_factor"].GetDouble();
    mCheckAtEachStep = rParameters["check_at_each_time_step"].GetBool();
    mNegElemDeactivation = rParameters["deactivate_full_negative_elements"].GetBool();
    mRecoverOriginalDistance = rParameters["recover_original_distance_at_each_step"].GetBool();
}

/* Protected functions ****************************************************/


/* Private functions ****************************************************/

};  // namespace Kratos.

