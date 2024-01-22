// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Wijtze Pieter Kikstra
//

// System includes

// External includes

// Project includes
#include "extrapolate_and_smooth_results_process.hpp"

#include "geo_mechanics_application_variables.h"


namespace Kratos
{

    ExtrapolateAndSmoothResultsProcess::ExtrapolateAndSmoothResultsProcess(ModelPart& rModelPart,
        const Parameters& rSettings)
        : Process(),
        mrModelPart(rModelPart),
        mParameters(rSettings)
    {
        // function type: python, cpp, input
        const Parameters default_parameters(R"(
        {
            "help"              : "This process sets a parameter field on a model part, where each element can have different material properties.",
            "model_part_name"   : "please_specify_model_part_name",
            "variable_name"     : "CUSTOM",
            "func_type"         : "input",               
            "function"          : "0",
            "dataset"           : "dummy",
            "dataset_file_name" : "dummy",
            "vector_variable_indices"      : []
        }  )"
        );

        //mParameters.RecursivelyValidateAndAssignDefaults(default_parameters);
    }

void ExtrapolateAndSmoothResultsProcess::ExecuteBeforeOutputStep()
{
    if (!mrModelPart.GetProcessInfo().Has(NODAL_SMOOTHING) ||
        !mrModelPart.GetProcessInfo()[NODAL_SMOOTHING]) return;
        // loop over variables from process parameters ( currently CAUCHY_STRESS, DAMAGE_VARIABLE, JOINT_WIDTH, JOINT_DAMAGE )
        // From parameters find i.p. variable name to extrapolate
        // check for its existence.
        // from i.p. variable get its size and then claim memory and fill with zeros
        const unsigned int dim = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];
        // bij modellen met structurele elementen en interfaces erin gaat dit mank.
        const auto stress_tensor_size =
                dim == N_DIM_3D ? STRESS_TENSOR_SIZE_3D : STRESS_TENSOR_SIZE_2D;

        // Clear nodal variables
        Matrix zero_nodal_stress_tensor = zero_matrix(stress_tensor_size, stress_tensor_size);
        block_for_each(mrModelPart.Nodes(), [&zero_nodal_stress_tensor](Node &rNode) {
            rNode.FastGetSolutionStepValue(NODAL_AREA)                 = 0.0;
            rNode.FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR) = zero_nodal_stress_tensor;
        });
/*
        // Walk the elements and extrapolate values to nodes
        block_for_each(mrModelPart.Elements(), [](Element& rElement)
        {
            // Walk integration points and retrieve i.p. result
            // extrapolate those to element nodes
            // add element nodal values to system nodal vector
        }
        // Compute smoothed nodal variables
        block_for_each(mrModelPart.Nodes(), [](Node& rNode)
        {
            if (const double& nodal_area = rNode.FastGetSolutionStepValue(NODAL_AREA);
                    nodal_area > 1.0e-20)
            {
                const double inv_nodal_area = 1.0 / nodal_area;
                rNode.FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR) *= inv_nodal_area;
            }
        });

    }
    */
}

} // namespace Kratos.