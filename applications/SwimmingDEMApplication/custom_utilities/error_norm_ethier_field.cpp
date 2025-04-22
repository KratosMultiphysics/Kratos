#include "error_norm_ethier_field.h"

namespace Kratos
{
    double ErrorNormEthierFieldCalculator::getL2NormFluidAccelWithoutRecoveryUsingGaussInterpolatedValues(ModelPart &r_model_part)
    {
        ProcessInfo process_info = r_model_part.GetProcessInfo();
        const unsigned int dim = process_info[DOMAIN_SIZE];
        const int number_of_elements = r_model_part.NumberOfElements();

        // Compute L2 error = sum_e sum_g sum_d W_g * (u_d(x_g) - u_d,exact(x_g))^2,
        // where u_d and u_d,exact are the piecewise functions of the variable and the exact variable, respectively
        double l2_norm_factor = 0.0, result = 0.0;
        for (int e = 0; e < number_of_elements; e++)
        {
            ModelPart::ElementsContainerType::iterator rElement = r_model_part.ElementsBegin() + e;

            const GeometryType &r_geometry = rElement->GetGeometry();
            unsigned int num_nodes = r_geometry.size();
            const GeometryData::IntegrationMethod integration_method = rElement->GetIntegrationMethod();
            const auto &integration_points = r_geometry.IntegrationPoints(integration_method);
            const auto &r_number_integration_points = r_geometry.IntegrationPointsNumber(integration_method);
            Matrix NContainer = r_geometry.ShapeFunctionsValues(integration_method);

            ShapeFunctionDerivativesArrayType shape_derivatives;
            Vector DetJ;
            r_geometry.ShapeFunctionsIntegrationPointsGradients(shape_derivatives, DetJ, integration_method);

            // Velocities and exact var values evaluated at gauss points
            Matrix vel_gauss_points = ZeroMatrix(r_number_integration_points, dim), exact_var_values_at_gauss_points = ZeroMatrix(r_number_integration_points, dim);
            for (unsigned g = 0; g < r_number_integration_points; g++)
            {
                for (unsigned n = 0; n < num_nodes; n++)
                {
                    array_1d<double, 3> u_nodal = r_geometry[n].FastGetSolutionStepValue(VELOCITY);
                    // array_1d<double, 3> exact_var_node_values = r_geometry[n].FastGetSolutionStepValue(EXACT_MATERIAL_ACCELERATION);
                    array_1d<double, 3> exact_var_node_values;
                    mFlowField.CalculateMaterialAcceleration(0.0, r_geometry[n].Coordinates(), exact_var_node_values, omp_get_thread_num());
                    for (unsigned i = 0; i < dim; i++)
                    {
                        vel_gauss_points(g, i) += u_nodal[i] * NContainer(g, n);
                        exact_var_values_at_gauss_points(g, i) += exact_var_node_values[i] * NContainer(g, n);
                    }
                }
            }

            // Fluid acceleration at gauss points
            Matrix fluid_accel_gauss_points = ZeroMatrix(r_number_integration_points, dim);
            for (unsigned g = 0; g < r_number_integration_points; g++)
            {
                for (unsigned n = 0; n < num_nodes; n++)
                {
                    array_1d<double, 3> u_nodal = r_geometry[n].FastGetSolutionStepValue(VELOCITY);
                    for (unsigned i = 0; i < dim; i++)
                    {
                        for (unsigned j = 0; j < dim; j++)
                        {
                            fluid_accel_gauss_points(g, i) += vel_gauss_points(g, j) * (shape_derivatives[g](n, j) * u_nodal[i]);
                        }
                    }
                }
            }

            // std::cout << "fluid_accel_gauss_points = \n" << fluid_accel_gauss_points << std::endl;
            // std::cout << "exact_var_values_at_gauss_points = \n" << exact_var_values_at_gauss_points << std::endl << std::endl;

            // Compute the error of the material acceleration computed using the derivatives of the shape function
            for (unsigned int g = 0; g < r_number_integration_points; g++)
            {
                double weight = DetJ[g] * integration_points[g].Weight();
                for (unsigned int d = 0; d < dim; d++)
                {
                    // Integrate the difference between exact and computed variable
                    result += weight * std::pow(fluid_accel_gauss_points(g, d) - exact_var_values_at_gauss_points(g, d), 2.0);
                    if(mNormalizeResult)
                        l2_norm_factor += weight * exact_var_values_at_gauss_points(g, d) * exact_var_values_at_gauss_points(g, d);
                }
            }
            if(!mNormalizeResult)
                l2_norm_factor += r_geometry.Area();
        }
        std::cout << "L2 norm factor: " << l2_norm_factor << ", result = " << result << std::endl;
        return std::sqrt(result / l2_norm_factor);
    }

    double ErrorNormEthierFieldCalculator::getL2NormFluidAccelWithRecoveryUsingGaussInterpolatedValues(ModelPart &r_model_part)
    {
        ProcessInfo process_info = r_model_part.GetProcessInfo();
        const unsigned int dim = process_info[DOMAIN_SIZE];
        const int number_of_elements = r_model_part.NumberOfElements();

        // Compute L2 error = sum_e sum_g sum_d W_g * (u_d(x_g) - u_d,exact(x_g))^2,
        // where u_d and u_d,exact are the piecewise functions of the variable and the exact variable, respectively
        double l2_norm_factor = 0.0, result = 0.0;
        for (int e = 0; e < number_of_elements; e++)
        {
            ModelPart::ElementsContainerType::iterator rElement = r_model_part.ElementsBegin() + e;

            const GeometryType &r_geometry = rElement->GetGeometry();
            unsigned int num_nodes = r_geometry.size();
            const GeometryData::IntegrationMethod integration_method = rElement->GetIntegrationMethod();
            const auto &integration_points = r_geometry.IntegrationPoints(integration_method);
            const auto &r_number_integration_points = r_geometry.IntegrationPointsNumber(integration_method);
            Matrix NContainer = r_geometry.ShapeFunctionsValues(integration_method);

            ShapeFunctionDerivativesArrayType shape_derivatives;
            Vector DetJ;
            r_geometry.ShapeFunctionsIntegrationPointsGradients(shape_derivatives, DetJ, integration_method);

            // Material acceleration at gauss points for the recovered field
            Matrix fluid_accel_gauss_points = ZeroMatrix(r_number_integration_points, dim), exact_fluid_accel_gauss_points = ZeroMatrix(r_number_integration_points, dim);
            for (unsigned g = 0; g < r_number_integration_points; g++)
            {
                for (unsigned n = 0; n < num_nodes; n++)
                {
                    array_1d<double, 3> fluid_accel_nodal = r_geometry[n].FastGetSolutionStepValue(MATERIAL_ACCELERATION);
                    array_1d<double, 3> exact_fluid_accel_nodal;
                    mFlowField.CalculateMaterialAcceleration(0.0, r_geometry[n].Coordinates(), exact_fluid_accel_nodal, omp_get_thread_num());
                    for (unsigned i = 0; i < dim; i++)
                    {
                        fluid_accel_gauss_points(g, i) += fluid_accel_nodal[i] * NContainer(g, n);
                        exact_fluid_accel_gauss_points(g, i) += exact_fluid_accel_nodal[i] * NContainer(g, n);
                    }
                }
            }

            // std::cout << "fluid_accel_gauss_points = \n" << fluid_accel_gauss_points << std::endl;
            // std::cout << "exact_var_values_at_gauss_points = \n" << exact_var_values_at_gauss_points << std::endl << std::endl;

            // Compute the error of the material acceleration computed using the derivatives of the shape function
            for (unsigned int g = 0; g < r_number_integration_points; g++)
            {
                double weight = DetJ[g] * integration_points[g].Weight();
                for (unsigned int d = 0; d < dim; d++)
                {
                    // Integrate the difference between exact and computed variable
                    result += weight * std::pow(fluid_accel_gauss_points(g, d) - exact_fluid_accel_gauss_points(g, d), 2.0);
                    if(mNormalizeResult)
                        l2_norm_factor += weight * exact_fluid_accel_gauss_points(g, d) * exact_fluid_accel_gauss_points(g, d);
                }
            }
            if(!mNormalizeResult)
                l2_norm_factor += r_geometry.Area();
        }
        std::cout << "L2 norm factor: " << l2_norm_factor << ", result = " << result << std::endl;
        return std::sqrt(result / l2_norm_factor);
    }

    double ErrorNormEthierFieldCalculator::getL2NormFluidAccelWithoutRecoveryUsingGaussExactValues(ModelPart &r_model_part)
    {
        ProcessInfo process_info = r_model_part.GetProcessInfo();
        const unsigned int dim = process_info[DOMAIN_SIZE];
        const int number_of_elements = r_model_part.NumberOfElements();

        // Compute L2 error = sum_e sum_g sum_d W_g * (u_d(x_g) - u_d,exact(x_g))^2,
        // where u_d and u_d,exact are the piecewise functions of the variable and the exact variable, respectively
        double l2_norm_factor = 0.0, result = 0.0;
        for (int e = 0; e < number_of_elements; e++)
        {
            ModelPart::ElementsContainerType::iterator rElement = r_model_part.ElementsBegin() + e;

            const GeometryType &r_geometry = rElement->GetGeometry();
            unsigned int num_nodes = r_geometry.size();
            const GeometryData::IntegrationMethod integration_method = rElement->GetIntegrationMethod();
            const auto &integration_points = r_geometry.IntegrationPoints(integration_method);
            const auto &r_number_integration_points = r_geometry.IntegrationPointsNumber(integration_method);
            Matrix NContainer = r_geometry.ShapeFunctionsValues(integration_method);

            ShapeFunctionDerivativesArrayType shape_derivatives;
            Vector DetJ;
            r_geometry.ShapeFunctionsIntegrationPointsGradients(shape_derivatives, DetJ, integration_method);

            // Velocities and exact var values evaluated at gauss points
            Matrix vel_gauss_points = ZeroMatrix(r_number_integration_points, dim);
            for (unsigned g = 0; g < r_number_integration_points; g++)
            {
                for (unsigned n = 0; n < num_nodes; n++)
                {
                    array_1d<double, 3> u_nodal = r_geometry[n].FastGetSolutionStepValue(VELOCITY);
                    for (unsigned i = 0; i < dim; i++)
                    {
                        vel_gauss_points(g, i) += u_nodal[i] * NContainer(g, n);
                    }
                }
            }

            // Fluid acceleration at gauss points
            Matrix fluid_accel_gauss_points = ZeroMatrix(r_number_integration_points, dim);
            for (unsigned g = 0; g < r_number_integration_points; g++)
            {
                for (unsigned n = 0; n < num_nodes; n++)
                {
                    array_1d<double, 3> u_nodal = r_geometry[n].FastGetSolutionStepValue(VELOCITY);
                    for (unsigned i = 0; i < dim; i++)
                    {
                        for (unsigned j = 0; j < dim; j++)
                        {
                            fluid_accel_gauss_points(g, i) += vel_gauss_points(g, j) * (shape_derivatives[g](n, j) * u_nodal[i]);
                        }
                    }
                }
            }

            // std::cout << "fluid_accel_gauss_points = \n" << fluid_accel_gauss_points << std::endl;
            // std::cout << "exact_var_values_at_gauss_points = \n" << exact_var_values_at_gauss_points << std::endl << std::endl;

            // Gauss global coordinates
            std::vector<array_1d<double, 3>> gauss_points_coordinates(r_number_integration_points, ZeroVector(3));
            for (unsigned g = 0; g < r_number_integration_points; g++)
            {
                // gauss_points_coordinates[g] = ZeroVector(3);
                for (unsigned n = 0; n < r_geometry.size(); n++)
                {
                    Vector node_coords = r_geometry[n].Coordinates();
                    for (unsigned d = 0; d < dim; d++)
                    {
                        gauss_points_coordinates[g][d] += node_coords[d] * NContainer(g, n);
                    }
                }
            }

            // Compute the error of the material acceleration computed using the derivatives of the shape function
            for (unsigned int g = 0; g < r_number_integration_points; g++)
            {
                array_1d<double, 3> exact_material_acceleration_gauss_point;
                // CalculateMaterialAcceleration(gauss_points_coordinates[g], exact_material_acceleration_gauss_point);
                mFlowField.CalculateMaterialAcceleration(0.0, gauss_points_coordinates[g], exact_material_acceleration_gauss_point, omp_get_thread_num());

                double weight = DetJ[g] * integration_points[g].Weight();
                for (unsigned int d = 0; d < dim; d++)
                {
                    // Integrate the difference between exact and computed variable
                    result += weight * std::pow(fluid_accel_gauss_points(g, d) - exact_material_acceleration_gauss_point[d], 2.0);
                    if(mNormalizeResult)
                        l2_norm_factor += weight * exact_material_acceleration_gauss_point[d] * exact_material_acceleration_gauss_point[d];
                }
            }
            if(!mNormalizeResult)
                l2_norm_factor += r_geometry.Area();
        }
        std::cout << "L2 norm factor: " << l2_norm_factor << ", result = " << result << std::endl;
        return std::sqrt(result / l2_norm_factor);
    }

    double ErrorNormEthierFieldCalculator::getL2NormFluidAccelWithRecoveryUsingGaussExactValues(ModelPart &r_model_part)
    {
        ProcessInfo process_info = r_model_part.GetProcessInfo();
        const unsigned int dim = process_info[DOMAIN_SIZE];
        const int number_of_elements = r_model_part.NumberOfElements();

        // Compute L2 error = sum_e sum_g sum_d W_g * (u_d(x_g) - u_d,exact(x_g))^2,
        // where u_d and u_d,exact are the piecewise functions of the variable and the exact variable, respectively
        double l2_norm_factor = 0.0, result = 0.0;
        for (int e = 0; e < number_of_elements; e++)
        {
            ModelPart::ElementsContainerType::iterator rElement = r_model_part.ElementsBegin() + e;

            const GeometryType &r_geometry = rElement->GetGeometry();
            unsigned int num_nodes = r_geometry.size();
            const GeometryData::IntegrationMethod integration_method = rElement->GetIntegrationMethod();
            const auto &integration_points = r_geometry.IntegrationPoints(integration_method);
            const auto &r_number_integration_points = r_geometry.IntegrationPointsNumber(integration_method);
            Matrix NContainer = r_geometry.ShapeFunctionsValues(integration_method);

            ShapeFunctionDerivativesArrayType shape_derivatives;
            Vector DetJ;
            r_geometry.ShapeFunctionsIntegrationPointsGradients(shape_derivatives, DetJ, integration_method);

            // Fluid acceleration at gauss points
            Matrix fluid_accel_gauss_points = ZeroMatrix(r_number_integration_points, dim);
            for (unsigned g = 0; g < r_number_integration_points; g++)
            {
                for (unsigned n = 0; n < num_nodes; n++)
                {
                    array_1d<double, 3> fluid_accel_nodal = r_geometry[n].FastGetSolutionStepValue(MATERIAL_ACCELERATION);
                    for (unsigned j = 0; j < dim; j++)
                    {
                        fluid_accel_gauss_points(g, j) += fluid_accel_nodal[j] * NContainer(g, n);
                    }
                }
            }

            // std::cout << "fluid_accel_gauss_points = \n" << fluid_accel_gauss_points << std::endl;
            // std::cout << "exact_var_values_at_gauss_points = \n" << exact_var_values_at_gauss_points << std::endl << std::endl;

            // Gauss global coordinates
            std::vector<array_1d<double, 3>> gauss_points_coordinates(r_number_integration_points, ZeroVector(3));
            for (unsigned g = 0; g < r_number_integration_points; g++)
            {
                // gauss_points_coordinates[g] = ZeroVector(3);
                for (unsigned n = 0; n < r_geometry.size(); n++)
                {
                    Vector node_coords = r_geometry[n].Coordinates();
                    for (unsigned d = 0; d < dim; d++)
                    {
                        gauss_points_coordinates[g][d] += node_coords[d] * NContainer(g, n);
                    }
                }
            }

            // Compute the error of the material acceleration computed using the derivatives of the shape function
            for (unsigned int g = 0; g < r_number_integration_points; g++)
            {
                double weight = DetJ[g] * integration_points[g].Weight();

                // Exact fluid's acceleration at gauss points
                array_1d<double, 3> exact_material_acceleration_gauss_point;
                // CalculateMaterialAcceleration(gauss_points_coordinates[g], exact_material_acceleration_gauss_point);
                mFlowField.CalculateMaterialAcceleration(0.0, gauss_points_coordinates[g], exact_material_acceleration_gauss_point, omp_get_thread_num());

                // std::cout << "Gauss point xg = " << gauss_points_coordinates[g] << " = " << gauss_points_coordinates[g] / mMajorRadius << ":" << std::endl;
                // std::cout << "  - rho   = " << getDistanceToCenter(gauss_points_coordinates[g]) / mMinorRadius << std::endl;
                // std::cout << "  - d     = " << std::sqrt(gauss_points_coordinates[g][0] * gauss_points_coordinates[g][0] + gauss_points_coordinates[g][1] * gauss_points_coordinates[g][1]) / mMajorRadius << std::endl;
                // std::cout << "  - u(xg) = " << getVelocity(gauss_points_coordinates[g]) / mCenterVelocity << std::endl;
                // std::cout << "  - Dtu   = " << exact_material_acceleration_gauss_point << std::endl;
                // std::cout << "  - x_c   = " << r_geometry.Center() << std::endl;
                // for (size_t i = 0; i < r_geometry.PointsNumber(); i++)
                // {
                //     std::cout << "  - x_n   = " << r_geometry.Points()[i].Coordinates() << std::endl;
                // }

                for (unsigned int d = 0; d < dim; d++)
                {
                    // Integrate the difference between exact and computed variable
                    result += weight * std::pow(fluid_accel_gauss_points(g, d) - exact_material_acceleration_gauss_point[d], 2.0);
                    if(mNormalizeResult)
                        l2_norm_factor += weight * exact_material_acceleration_gauss_point[d] * exact_material_acceleration_gauss_point[d];
                }
            }
            if(!mNormalizeResult)
                l2_norm_factor += r_geometry.Area();
        }
        std::cout << "L2 norm factor: " << l2_norm_factor << ", result = " << result << std::endl;
        return std::sqrt(result / l2_norm_factor);
    }
}