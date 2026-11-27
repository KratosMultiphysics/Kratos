//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

// System includes
// External includes

// Project includes
#include "snake_gap_sbm_process.h"
#include "integration/integration_point_utilities.h"
#include "iga_application_variables.h"
#include "includes/global_pointer_variables.h"

namespace Kratos
{
namespace
{
} // unnamed namespace

SnakeGapSbmProcess::SnakeGapSbmProcess(
    Model& rModel, Parameters ThisParameters) : 
    SnakeSbmProcess(
        rModel,
        // Reject user-defined lambdas before the base constructor adds their default values.
        (ThisParameters.Has("lambda_inner") || ThisParameters.Has("lambda_outer"))
            ? KRATOS_ERROR << "[SnakeGapSbmProcess]: \"lambda_inner\" and \"lambda_outer\" "
                           << "must not be provided in ProjectParameters. GapSBM assigns these "
                           << "values internally." << std::endl
            : ThisParameters)
{
    
    KRATOS_ERROR_IF_NOT(ThisParameters.Has("gap_element_name")) << "::[SnakeGapSbmProcess]::" 
                    << "Missing \"gap_element_name\" section." << std::endl;
    KRATOS_ERROR_IF_NOT(ThisParameters.Has("gap_interface_condition_name")) << "::[SnakeGapSbmProcess]::" 
                    << "Missing \"gap_interface_condition_name\" section." << std::endl;

    ThisParameters.AddMissingParameters(this->GetDefaultParameters());
    
    mpGapElementsSubModelPart = &(mpIgaModelPart->CreateSubModelPart("GapElements"));
    mpGapInterfaceSubModelPart = &(mpIgaModelPart->CreateSubModelPart("GapInterfaces"));
    mGapElementName = ThisParameters["gap_element_name"].GetString();
    mGapInterfaceConditionName = ThisParameters["gap_interface_condition_name"].GetString();
    mGapSbmType = ThisParameters["gap_sbm_type"].GetString(); 

    if (mGapSbmType != "default" && mGapSbmType != "interpolation" && mGapSbmType != "sbm") {
        KRATOS_ERROR << "::[SnakeGapSbmProcess]::"
                     << "The gap_sbm_type \"" << mGapSbmType << "\" is not supported. Available options are: "
                     << "default, interpolation, sbm." << std::endl;
    }

    if (ThisParameters.Has("gap_approximation_order"))
        mGapApproximationOrder = ThisParameters["gap_approximation_order"].GetInt();
    if (ThisParameters.Has("gap_relative_tolerance_for_subdivisions"))
        mGapRelativeToleranceForSubdivisions = ThisParameters["gap_relative_tolerance_for_subdivisions"].GetDouble();
    if (ThisParameters.Has("number_of_interpolation_levels"))
        mNumberOfInterpolationLevels = ThisParameters["number_of_interpolation_levels"].GetInt();
    if (ThisParameters.Has("polynomial_order"))
        mGapInterpolationOrder = ThisParameters["polynomial_order"].GetVector()[0];

    mLambdaInner = 0.0;
    mLambdaOuter = 1.0;
    mInternalDivisions = ThisParameters["number_internal_divisions"].GetInt();

}

void SnakeGapSbmProcess::CreateSbmExtendedGeometries()
{
    CreateSbmExtendedGeometries2D();
}

const Parameters SnakeGapSbmProcess::GetDefaultParameters() const
{
    return Parameters(R"(
    {
        "echo_level": 1,
        "model_part_name" : "IgaModelPart",
        "lower_point_xyz": [0.0, 0.0, 0.0],
        "upper_point_xyz": [1.0, 1.0, 0.0],
        "polynomial_order" : [2, 2],
        "number_of_knot_spans" : [10, 10],
        "gap_relative_tolerance_for_subdivisions": 0.1,
        "number_of_interpolation_levels": 3
    })");
}

const Parameters SnakeGapSbmProcess::GetValidParameters() const
{
    return Parameters(R"(
    {
        "echo_level": 0,
        "lower_point_xyz": [-0.5, -0.5,0.0],
        "upper_point_xyz": [0.5,0.5,0.0],
        "polynomial_order" : [2, 2],
        "number_of_knot_spans" : [7, 7],
        "number_of_inner_loops": 0,
        "number_initial_points_if_importing_nurbs": 1000,
        "number_internal_divisions": 0,
        "gap_relative_tolerance_for_subdivisions": 0.1,
        "number_of_interpolation_levels": 3,
        "gap_sbm_type": "default",
        "skin_model_part_outer_initial_name": "initial_skin_model_part_out",    
        "skin_model_part_inner_initial_name": "initial_skin_model_part_in",           
        "skin_model_part_name": "skin_model_part",
        "gap_element_name": "CutSbmSolidElement",
        "gap_interface_condition_name": "CutSbmSolidInterfaceCondition"
    })");
}

}  // namespace Kratos.
