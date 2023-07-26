#ifndef GRADUAL_VARIABLE_INTERPOLATION_UTILITY_H
#define GRADUAL_VARIABLE_INTERPOLATION_UTILITY_H

// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                       license: MeshingApplication/license.txt
//
//  Main authors:    Sebastian Ares de Parga Regalado

// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes

// Project includes
#include "processes/process.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "custom_processes/nodal_values_interpolation_process.h"
#include "includes/parallel_environment.h"
#include "utilities/parallel_utilities.h"


namespace Kratos {
class GradualVariableInterpolationUtility {
public:
    KRATOS_CLASS_POINTER_DEFINITION(GradualVariableInterpolationUtility);

    static void InitializeInterpolationAndConstraints(
        ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart,
        std::vector<Variable<double>>& rInterpolationVariablesList,
        double AlphaRampUpIncrement,
        int DomainSize,
        bool ConstrainVariables) 
        {
            // Set variables field to origin model part
            block_for_each(rOriginModelPart.Nodes(), [&](Node& node) {
                for(auto& variable : rInterpolationVariablesList) {
                    const double value = node.FastGetSolutionStepValue(variable);
                    node.SetValue(variable, AlphaRampUpIncrement * value);
                }
            });

            // Interpolate variables to destination model part
            if (DomainSize == 2) {
                NodalValuesInterpolationProcess<2> interpolation(rOriginModelPart, rDestinationModelPart);
                interpolation.Execute();
            } else if (DomainSize == 3) {
                NodalValuesInterpolationProcess<3> interpolation(rOriginModelPart, rDestinationModelPart);
                interpolation.Execute();
            } else {
                KRATOS_ERROR << "Invalid DomainSize: " << DomainSize << ". DomainSize should be 2 or 3." << std::endl;
            }

            // Using the parallel utilities
            if(ConstrainVariables)
            {
                // Apply constraints to variables
                block_for_each(rDestinationModelPart.Nodes(), [&](Node& rNode){
                    for(auto& variable : rInterpolationVariablesList){
                        auto& r_value = rNode.GetValue(variable);
                        rNode.FastGetSolutionStepValue(variable) = r_value;
                        rNode.Fix(variable);
                    }
                });
            }
            else
            {
                // Simply set the variables
                block_for_each(rDestinationModelPart.Nodes(), [&](Node& rNode){
                    for(auto& variable : rInterpolationVariablesList){
                        auto& r_value = rNode.GetValue(variable);
                        rNode.FastGetSolutionStepValue(variable) = r_value;
                    }
                });
            }
        }
    
    static void UpdateSolutionStepVariables(
        ModelPart& rDestinationModelPart,
        std::vector<Variable<double>>& rInterpolationVariablesList,
        double& rAlpha,
        double& rOldAlpha,
        bool ConstrainVariables) 
        {
            if(ConstrainVariables)
            {
                block_for_each(rDestinationModelPart.Nodes(), [&](Node& rNode){
                    for(auto& variable : rInterpolationVariablesList){
                        double value = rNode.FastGetSolutionStepValue(variable);
                        rNode.FastGetSolutionStepValue(variable) = rAlpha * value / rOldAlpha;
                        rNode.Fix(variable);
                    }
                });
            }
            else
            {
                block_for_each(rDestinationModelPart.Nodes(), [&](Node& rNode){
                    for(auto& variable : rInterpolationVariablesList){
                        double value = rNode.FastGetSolutionStepValue(variable);
                        rNode.FastGetSolutionStepValue(variable) = rAlpha * value / rOldAlpha;
                    }
                });
            }
        }
};

}  // namespace Kratos.

#endif /* GRADUAL_VARIABLE_INTERPOLATION_UTILITY_H defined */


