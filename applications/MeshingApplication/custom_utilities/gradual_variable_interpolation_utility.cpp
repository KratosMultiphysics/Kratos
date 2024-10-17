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
#include "custom_utilities/gradual_variable_interpolation_utility.h"
#include "includes/kratos_components.h"

namespace Kratos {

void GradualVariableInterpolationUtility::InitializeInterpolationAndConstraints(
    ModelPart& rOriginModelPart,
    ModelPart& rDestinationModelPart,
    std::vector<std::string>& rInterpolationVariablesList,
    double AlphaRampUpIncrement,
    int DomainSize,
    bool ConstrainVariables) 
{
    {
        // Set variables field to origin model part
        for(auto& r_variable_name : rInterpolationVariablesList){
            auto& r_variable = KratosComponents<Variable<double>>::Get(r_variable_name);
            for(auto& rNode: rOriginModelPart.Nodes()) {
                const double value = rNode.FastGetSolutionStepValue(r_variable);
                rNode.SetValue(r_variable, AlphaRampUpIncrement * value);
            }
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

        for(auto& r_variable_name : rInterpolationVariablesList){
            auto& r_variable = KratosComponents<Variable<double>>::Get(r_variable_name);
            for(auto& rNode: rDestinationModelPart.Nodes()) {
                rNode.AddDof(r_variable);
            }
        }

        // Using the parallel utilities
        if(ConstrainVariables)
        {
            // Apply constraints to variables
            for(auto& r_variable_name : rInterpolationVariablesList){
                auto& r_variable = KratosComponents<Variable<double>>::Get(r_variable_name);
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
            for(auto& r_variable_name : rInterpolationVariablesList){
                auto& r_variable = KratosComponents<Variable<double>>::Get(r_variable_name);
                block_for_each(rDestinationModelPart.Nodes(), [&](Node& rNode){
                    const double value = rNode.GetValue(r_variable);
                    rNode.FastGetSolutionStepValue(r_variable) = value;
                });
            }
        }
    }
}

void GradualVariableInterpolationUtility::UpdateSolutionStepVariables(
    ModelPart& rDestinationModelPart,
    std::vector<std::string>& rInterpolationVariablesList,
    double& rAlpha,
    double& rOldAlpha,
    bool ConstrainVariables) 
{
    {
        if(ConstrainVariables)
        {
            for(auto& r_variable_name : rInterpolationVariablesList){
                auto& r_variable = KratosComponents<Variable<double>>::Get(r_variable_name);
                block_for_each(rDestinationModelPart.Nodes(), [&](Node& rNode){
                    const double value = rNode.FastGetSolutionStepValue(r_variable);
                    rNode.FastGetSolutionStepValue(r_variable) = rAlpha * value / rOldAlpha;
                    rNode.Fix(r_variable);
                });
            }
        }
        else
        {
            for(auto& r_variable_name : rInterpolationVariablesList){
                auto& r_variable = KratosComponents<Variable<double>>::Get(r_variable_name);
                block_for_each(rDestinationModelPart.Nodes(), [&](Node& rNode){
                    const double value = rNode.FastGetSolutionStepValue(r_variable);
                    rNode.FastGetSolutionStepValue(r_variable) = rAlpha * value / rOldAlpha;
                });
            }
        }
    }
}

}  // namespace Kratos.
