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
            for(auto& r_variable : rInterpolationVariablesList){
                block_for_each(rOriginModelPart.Nodes(), [&](Node& rNode) {
                    const double value = rNode.FastGetSolutionStepValue(r_variable);
                    rNode.SetValue(r_variable, AlphaRampUpIncrement * value);
                });
            }

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

            for(auto& r_variable : rInterpolationVariablesList){
                for(auto& rNode: rDestinationModelPart.Nodes()) {
                    rNode.AddDof(r_variable);
                }
            }

            // Using the parallel utilities
            if(ConstrainVariables)
            {
                // Apply constraints to variables
                for(auto& r_variable : rInterpolationVariablesList){
                    block_for_each(rDestinationModelPart.Nodes(), [&](Node& rNode){
                        const double value = rNode.GetValue(r_variable);
                        rNode.FastGetSolutionStepValue(r_variable) = value;
                        rNode.Fix(r_variable);
                    });
                }
            }
            else
            {
                // Simply set the variables
                for(auto& r_variable : rInterpolationVariablesList){
                    block_for_each(rDestinationModelPart.Nodes(), [&](Node& rNode){
                        const double value = rNode.GetValue(r_variable);
                        rNode.FastGetSolutionStepValue(r_variable) = value;
                    });
                }

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
                for(auto& r_variable : rInterpolationVariablesList){
                    block_for_each(rDestinationModelPart.Nodes(), [&](Node& rNode){
                        const double value = rNode.FastGetSolutionStepValue(r_variable);
                        rNode.FastGetSolutionStepValue(r_variable) = rAlpha * value / rOldAlpha;
                        rNode.Fix(r_variable);
                    });
                }
            }
            else
            {
                for(auto& variable : rInterpolationVariablesList){
                    block_for_each(rDestinationModelPart.Nodes(), [&](Node& rNode){
                        const double value = rNode.FastGetSolutionStepValue(variable);
                        rNode.FastGetSolutionStepValue(variable) = rAlpha * value / rOldAlpha;
                    });
                }
            }
        }
};

}  // namespace Kratos.

#endif /* GRADUAL_VARIABLE_INTERPOLATION_UTILITY_H defined */


