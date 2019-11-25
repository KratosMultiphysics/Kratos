//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Saransh Saxena
//

// System includes

// External includes

// Project includes
#include "custom_processes/simple_error_calculator_process.h"
#include "utilities/variable_utils.h"

namespace Kratos
{
template<SizeType TDim>
SimpleErrorCalculatorProcess<TDim>::SimpleErrorCalculatorProcess(ModelPart& rThisModelPart, Parameters ThisParameters):mrThisModelPart(rThisModelPart)
{
    //WIP
    Parameters default_parameters = Parameters(R"(
    {
        "minimal_size"                        : 0.01,
        "maximal_size"                        : 10.0,
        "refinement_strategy"                 : "Simple_Error_Calculator",
        "reference_variable_name"             : "ERROR_RATIO",
        "mean_distribution_strategy":
        {
            "target_refinement_coefficient"       : 0.9,
            "refinement_bound"                    : 2.0,
            "reference_norm_name"                 : "VELOCITY_H1_SEMINORM"
        },
        "maximum_strategy":
        {
            "target_refinement_coefficient"       : 0.1,
            "refinement_coefficient"              : 2.0
        },
        "global_tolerance_strategy":
        {
            "global_tolerance"                 : 0.1
        },
        "echo_level"                          : 0
    })"
    );

    ThisParameters.RecursivelyValidateAndAssignDefaults(default_parameters);