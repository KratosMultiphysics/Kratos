//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Joaquin Gonzalez-Usua
//
//

//These formulas are derived from "Derivatives.py" script

// System includes

// External includes
#include <omp.h>
#include <vector>
#include <iostream>
// Project includes
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"
#include "includes/kratos_parameters.h"
#include "includes/process_info.h"

// Application includes
#include "swimming_DEM_application.h"
#include "swimming_dem_application_variables.h"
#include "custom_utilities/error_norm_calculator_utility.h"

// Other applications includes
#include "fluid_dynamics_application_variables.h"
#include "fluid_dynamics_application.h"

namespace Kratos
{

void ErrorNormCalculator::ComputeDofsErrors(ModelPart& r_model_part)
{
    block_for_each(r_model_part.Nodes(), [&](Node& r_node)
    {
        auto& r_vectorial_error = r_node.FastGetSolutionStepValue(VECTORIAL_ERROR);
        r_vectorial_error = r_node.FastGetSolutionStepValue(VELOCITY) - r_node.FastGetSolutionStepValue(EXACT_VELOCITY);

        auto& r_scalar_error = r_node.FastGetSolutionStepValue(SCALAR_ERROR);
        r_scalar_error = r_node.FastGetSolutionStepValue(PRESSURE) - r_node.FastGetSolutionStepValue(EXACT_PRESSURE);
    });
}

double ErrorNormCalculator::GetL2VectorErrorNorm(ModelPart& r_model_part, const Variable<array_1d<double,3>>& rVariable, const Variable<array_1d<double,3>>& rExact_Variable)
{
    double total_area = 0.0, result = 0.0;
    ProcessInfo process_info = r_model_part.GetProcessInfo();
    const unsigned int dim = process_info[DOMAIN_SIZE];
    const int number_of_elements = r_model_part.NumberOfElements();

    // Compute L2 error = sum_e sum_g sum_d W_g * (u_d(x_g) - u_d,exact(x_g))^2,
    // where u_d and u_d,exact are the piecewise functions of the variable and the exact variable, respectively
    for (int e = 0; e < number_of_elements; e++)
    {
        ModelPart::ElementsContainerType::iterator rElement = r_model_part.ElementsBegin() + e;

        const GeometryType& r_geometry = rElement->GetGeometry();
        unsigned int num_nodes = r_geometry.size();
        const GeometryData::IntegrationMethod integration_method = rElement->GetIntegrationMethod();
        const auto& integration_points = r_geometry.IntegrationPoints(integration_method);
        const auto& r_number_integration_points = r_geometry.IntegrationPointsNumber(integration_method);
        Matrix NContainer = r_geometry.ShapeFunctionsValues(integration_method);

        ShapeFunctionDerivativesArrayType shape_derivatives;
        Vector DetJ;
        r_geometry.ShapeFunctionsIntegrationPointsGradients(shape_derivatives, DetJ, integration_method);

        double variables_diff_gauss_point = 0.0;

        for (unsigned int g = 0; g < r_number_integration_points; g++){
            double weight = DetJ[g] * integration_points[g].Weight();
            for (unsigned int d = 0; d < dim; d++)
            {
                // Interpolate the difference of the values of the variables at the integration points
                for (unsigned int a = 0; a < num_nodes; a++)
                {
                    array_1d<double, 3> var_nodal_values = r_geometry[a].FastGetSolutionStepValue(rVariable);
                    array_1d<double, 3> exact_var_nodal_values = r_geometry[a].FastGetSolutionStepValue(rExact_Variable);
                    variables_diff_gauss_point += (var_nodal_values[d] - exact_var_nodal_values[d]) * NContainer(g, a);
                }

                // Integrate the difference between exact and computed variable
                result += weight * std::pow(variables_diff_gauss_point, 2.0);
                variables_diff_gauss_point = 0.0;
            }
        }
        total_area += r_geometry.Area();
    }
    return std::sqrt(result / total_area);
}

double ErrorNormCalculator::GetL2ScalarErrorNorm(ModelPart& r_model_part, const Variable<double>& rVariable, const Variable<double>& rExact_Variable)
{
    double total_area = 0.0, result = 0.0;
    ProcessInfo process_info = r_model_part.GetProcessInfo();
    const int number_of_elements = r_model_part.NumberOfElements();

    for (int e = 0; e < number_of_elements; e++){
        ModelPart::ElementsContainerType::iterator rElement = r_model_part.ElementsBegin() + e;
        std::vector<double> exact_scalar, computed_scalar, scalar_error;

        const GeometryType& r_geometry = rElement->GetGeometry();
        const GeometryData::IntegrationMethod integration_method = rElement->GetIntegrationMethod();
        const auto& integration_points = r_geometry.IntegrationPoints(integration_method);
        const auto& r_number_integration_points = r_geometry.IntegrationPointsNumber(integration_method);

        total_area += r_geometry.Area();

        ShapeFunctionDerivativesArrayType shape_derivatives;
        Vector DetJ;
        Vector gauss_weights;
        r_geometry.ShapeFunctionsIntegrationPointsGradients(shape_derivatives,DetJ,integration_method);

        rElement->CalculateOnIntegrationPoints(rVariable,computed_scalar,process_info);

        rElement->CalculateOnIntegrationPoints(rExact_Variable,exact_scalar,process_info);

        std::transform(computed_scalar.begin(), computed_scalar.end(), exact_scalar.begin(), std::back_inserter(scalar_error),
                      [](double computed_scalar, double exact_scalar)
        { return std::pow(computed_scalar - exact_scalar,2.0);});


        if (gauss_weights.size() != r_number_integration_points)
            gauss_weights.resize(r_number_integration_points,false);

        for (unsigned int g = 0; g < r_number_integration_points; g++){
            gauss_weights[g] = DetJ[g] * integration_points[g].Weight();
            result += gauss_weights[g] * scalar_error[g];
        }
    }

    return std::sqrt(fabs(result/total_area));
}

double ErrorNormCalculator::GetH1ScalarErrorSemiNorm(ModelPart& r_model_part)
{
    double total_area = 0.0, result = 0.0;
    ProcessInfo process_info = r_model_part.GetProcessInfo();
    const unsigned int dim = process_info[DOMAIN_SIZE];

    //double squared_L2_scalar_norm = std::pow(this->GetL2ScalarErrorNorm(r_model_part),2.0);
    const int number_of_elements = r_model_part.NumberOfElements();

    for (int e = 0; e < number_of_elements; e++){
        ModelPart::ElementsContainerType::iterator rElement = r_model_part.ElementsBegin() + e;
        const GeometryType& r_geometry = rElement->GetGeometry();
        const GeometryData::IntegrationMethod integration_method = rElement->GetIntegrationMethod();
        const auto& integration_points = r_geometry.IntegrationPoints(integration_method);
        const auto& r_number_integration_points = r_geometry.IntegrationPointsNumber(integration_method);
        std::vector<array_1d<double,3>> exact_scalar_gradient;
        std::vector<array_1d<double,3>> computed_scalar_gradient;
        std::vector<double> scalar_gradient_error(r_number_integration_points,0.0);

        total_area += r_geometry.Area();

        ShapeFunctionDerivativesArrayType shape_derivatives;
        Vector DetJ;
        Vector gauss_weights;
        r_geometry.ShapeFunctionsIntegrationPointsGradients(shape_derivatives,DetJ,integration_method);

        rElement->CalculateOnIntegrationPoints(PRESSURE_GRADIENT,computed_scalar_gradient,process_info);

        rElement->CalculateOnIntegrationPoints(EXACT_PRESSURE_GRADIENT,exact_scalar_gradient,process_info);

        if (gauss_weights.size() != r_number_integration_points)
            gauss_weights.resize(r_number_integration_points,false);

        for (unsigned int g = 0; g < r_number_integration_points; g++){
            gauss_weights[g] = DetJ[g] * integration_points[g].Weight();
            for (unsigned int d = 0; d < dim; d++){
                scalar_gradient_error[g] += std::pow(exact_scalar_gradient[g][d] - computed_scalar_gradient[g][d], 2.0);
            }
            result += gauss_weights[g] * scalar_gradient_error[g];
        }
    }

    double squared_L2_gradient_scalar_norm = result/total_area;

    return  std::sqrt(squared_L2_gradient_scalar_norm);
}

double ErrorNormCalculator::GetH1VectorErrorSemiNorm(ModelPart& r_model_part)
{
    double total_area = 0.0, result = 0.0;
    ProcessInfo process_info = r_model_part.GetProcessInfo();
    const unsigned int dim = process_info[DOMAIN_SIZE];

    const int number_of_elements = r_model_part.NumberOfElements();

    for (int e = 0; e < number_of_elements; e++){
        ModelPart::ElementsContainerType::iterator rElement = r_model_part.ElementsBegin() + e;
        const GeometryType& r_geometry = rElement->GetGeometry();
        const GeometryData::IntegrationMethod integration_method = rElement->GetIntegrationMethod();
        const auto& integration_points = r_geometry.IntegrationPoints(integration_method);
        const auto& r_number_integration_points = r_geometry.IntegrationPointsNumber(integration_method);
        std::vector<Matrix> exact_vector_gradient;
        std::vector<Matrix> computed_vector_gradient;
        std::vector<double> vector_gradient_error(r_number_integration_points,0.0);

        total_area += r_geometry.Area();

        ShapeFunctionDerivativesArrayType shape_derivatives;
        Vector DetJ;
        Vector gauss_weights;
        r_geometry.ShapeFunctionsIntegrationPointsGradients(shape_derivatives,DetJ,integration_method);

        rElement->CalculateOnIntegrationPoints(VELOCITY_GRADIENT,computed_vector_gradient,process_info);

        rElement->CalculateOnIntegrationPoints(EXACT_VELOCITY_GRADIENT,exact_vector_gradient,process_info);

        if (gauss_weights.size() != r_number_integration_points)
            gauss_weights.resize(r_number_integration_points,false);

        for (unsigned int g = 0; g < r_number_integration_points; g++){
            gauss_weights[g] = DetJ[g] * integration_points[g].Weight();
            for (unsigned int d = 0; d < dim; d++){
                for (unsigned int e = 0; e < dim; e++)
                    vector_gradient_error[g] += std::pow(exact_vector_gradient[g](d,e) - computed_vector_gradient[g](d,e), 2.0);
            }
            result += gauss_weights[g] * vector_gradient_error[g];
        }
    }

    double squared_L2_gradient_vector_norm = result/total_area;

    return  std::sqrt(squared_L2_gradient_vector_norm);
}

/* Protected functions ****************************************************/



/* Private functions ****************************************************/

};  // namespace Kratos.