//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:   Juan I. Camarotti

// Project includes
#include "calculate_L2_norm_error_process.h"

namespace Kratos
{

    CalculateL2NormErrorProcess::CalculateL2NormErrorProcess(
        Model& rModel, Parameters ThisParameters) : mpModel(&rModel), mParameters(ThisParameters)
    {
        ThisParameters.ValidateAndAssignDefaults(this->GetDefaultParameters());
    }

    void CalculateL2NormErrorProcess::Execute(){
        // Read the compute error model part name
        std::string model_part_name = mParameters["compute_error_model_part_name"].GetString();
        mpComputeErrorModelPart = &mpModel->GetModelPart(model_part_name);
        // Read the unknown variable for the problem
        mUnknownVariable = mParameters["unknown_variable_name"].GetString();
        // Read the analytical solution
        std::string analytical_solution = mParameters["analytical_solution"].GetString();
        // Create the function to be evaluated
        mpAnalyticalSolution = Kratos::make_unique<GenericFunctionUtility>(analytical_solution);
    
        double L2_norm_error = 0.0;

        const Variable<double>& UNKNOWN_VARIABLE = KratosComponents<Variable<double>>::Get(mUnknownVariable);

        for (auto element_it = mpComputeErrorModelPart->ElementsBegin(); element_it != mpComputeErrorModelPart->ElementsEnd(); element_it++){
            const auto p_gauss_point_geometry = element_it->pGetGeometry();
            const auto gauss_point_position = p_gauss_point_geometry->Center();

            // Get the integration weight  
            const double weight = p_gauss_point_geometry->IntegrationPoints()[0].Weight();  

            // Determinant of jacobian
            Vector det_jacobian;
            p_gauss_point_geometry->DeterminantOfJacobian(det_jacobian);

            // Compute the value of the analytical solution at the gauss point
            double current_time = mpComputeErrorModelPart->GetProcessInfo()[TIME];
            const double value_analytical_solution = mpAnalyticalSolution->CallFunction(gauss_point_position.X(), gauss_point_position.Y(), gauss_point_position.Z(), current_time);

            // Compute the numerical solution at the gauss point
            double value_numerical_solution = 0.0;

            // Shape functions evaluated at the GP
            const CompressedMatrix& N = p_gauss_point_geometry->ShapeFunctionsValues();
            const GeometryType::ShapeFunctionsGradientsType& r_DN_De = p_gauss_point_geometry->ShapeFunctionsLocalGradients(element_it->GetIntegrationMethod());
            auto N_row = row(N, 0);

            double gradient_x = 0.0;
            double gradient_y = 0.0;
            // Iterate over the nodes in the geometry and retrieve the numerical solution and the numerical gradient
            for (std::size_t i = 0; i < p_gauss_point_geometry->size(); ++i) {
                auto& node = (*p_gauss_point_geometry)[i]; 
                double nodal_value = node.FastGetSolutionStepValue(UNKNOWN_VARIABLE, current_time); 
                value_numerical_solution += N_row[i] * nodal_value; 
                gradient_x += r_DN_De[0](i, 0) * nodal_value; 
                gradient_y += r_DN_De[0](i, 1) * nodal_value; 
            }


            // Compute the local error
            double analytical_gradient = std::sqrt(std::pow(gradient_x, 2) + std::pow(gradient_y, 2));
            double error = std::pow(value_analytical_solution - value_numerical_solution, 2);
            
            // Contribution to the L2 error 
            L2_norm_error += error * det_jacobian[0] * weight;
        }

        L2_norm_error = std::sqrt(L2_norm_error);

        KRATOS_WATCH(L2_norm_error)
        
    }
} // End namespace Kratos
