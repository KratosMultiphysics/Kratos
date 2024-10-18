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

#pragma once

// Project includes
#include "processes/process.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "custom_processes/nodal_values_interpolation_process.h"
#include "includes/parallel_environment.h"
#include "utilities/parallel_utilities.h"

namespace Kratos {

class KRATOS_API(MESHING_APPLICATION) GradualVariableInterpolationUtility {
public:
    KRATOS_CLASS_POINTER_DEFINITION(GradualVariableInterpolationUtility);

    /**
     * @brief Initializes interpolation and applies constraints to the given variables
     * @param rOriginModelPart The original model part
     * @param rDestinationModelPart The destination model part
     * @param rInterpolationVariablesList List of variables to be interpolated
     * @param AlphaRampUpIncrement Alpha increment for the ramp-up
     * @param DomainSize Size of the domain (2D or 3D)
     * @param ConstrainVariables Flag to apply constraints to the variables
     */
    static void InitializeInterpolationAndConstraints(
        ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart,
        std::vector<std::string>& rInterpolationVariablesList,
        double AlphaRampUpIncrement,
        int DomainSize,
        bool ConstrainVariables);

    /**
     * @brief Updates the solution step variables
     * @param rDestinationModelPart The destination model part
     * @param rInterpolationVariablesList List of variables to be interpolated
     * @param rAlpha Alpha value for updating the variables
     * @param rOldAlpha Previous alpha value
     * @param ConstrainVariables Flag to apply constraints to the variables
     */
    static void UpdateSolutionStepVariables(
        ModelPart& rDestinationModelPart,
        std::vector<std::string>& rInterpolationVariablesList,
        double& rAlpha,
        double& rOldAlpha,
        bool ConstrainVariables);
};

}  // namespace Kratos.

