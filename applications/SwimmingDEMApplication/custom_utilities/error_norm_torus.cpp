//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Aniol Sala
//
//

// These formulas are derived from "Derivatives.py" script

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
#include "custom_utilities/error_norm_torus.h"

// Other applications includes
#include "fluid_dynamics_application_variables.h"
#include "fluid_dynamics_application.h"

namespace Kratos
{
    double ErrorNormTorusCalculator::getDistanceToCenter(const array_1d<double, 3> &coor)
    {
        double major_dist = std::sqrt(coor[0] * coor[0] + coor[1] * coor[1]);
        double minor_dist = major_dist - mMajorRadius;
        double distance = std::sqrt(minor_dist * minor_dist + coor[2] * coor[2]);

        return distance;
    }

    double ErrorNormTorusCalculator::getVelocityModule(const array_1d<double, 3> &coor)
    {
        double rho = getDistanceToCenter(coor) / mMinorRadius;
        return mCenterVelocity * (1 - rho * rho);
    }

    array_1d<double, 3> ErrorNormTorusCalculator::getVelocity(const array_1d<double, 3> &coor)
    {
        double velocity_module = getVelocityModule(coor);
        double major_dist_2 = coor[0] * coor[0] + coor[1] * coor[1];
        double major_dist = std::sqrt(major_dist_2);
        double sin_theta = coor[0] / major_dist, cos_theta = coor[1] / major_dist;

        array_1d<double, 3> vel;
        vel[0] = velocity_module * cos_theta;
        vel[1] = velocity_module * (-1.0 * sin_theta);
        vel[2] = 0.0;

        return vel;
    }

    void ErrorNormTorusCalculator::CalculateMaterialAcceleration(const array_1d<double, 3> &coor, array_1d<double, 3> &accel)
    {
        double velocity_module = getVelocityModule(coor);
        double major_dist_2 = coor[0] * coor[0] + coor[1] * coor[1];

        accel[0] = -(coor[0] / major_dist_2) * velocity_module * velocity_module;
        accel[1] = -(coor[1] / major_dist_2) * velocity_module * velocity_module;
        accel[2] = 0.0;
    }

    double ErrorNormTorusCalculator::getL2NormFluidAccelWithoutRecoveryUsingGaussInterpolatedValues(ModelPart &r_model_part)
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
                    array_1d<double, 3> exact_var_node_values = r_geometry[n].FastGetSolutionStepValue(EXACT_MATERIAL_ACCELERATION);
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
                }
            }
            total_area += r_geometry.Area();
        }
        return std::sqrt(result / total_area);
    }

    double ErrorNormTorusCalculator::getL2NormFluidAccelWithRecoveryUsingGaussInterpolatedValues(ModelPart &r_model_part)
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
                    array_1d<double, 3> exact_fluid_accel_nodal = r_geometry[n].FastGetSolutionStepValue(EXACT_MATERIAL_ACCELERATION);
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
                }
            }
            total_area += r_geometry.Area();
        }
        return std::sqrt(result / total_area);
    }

    double ErrorNormTorusCalculator::getL2NormFluidAccelWithoutRecoveryUsingGaussExactValues(ModelPart &r_model_part)
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
                CalculateMaterialAcceleration(gauss_points_coordinates[g], exact_material_acceleration_gauss_point);

                double weight = DetJ[g] * integration_points[g].Weight();
                for (unsigned int d = 0; d < dim; d++)
                {
                    // Integrate the difference between exact and computed variable
                    result += weight * std::pow(fluid_accel_gauss_points(g, d) - exact_material_acceleration_gauss_point[d], 2.0);
                }
            }
            total_area += r_geometry.Area();
        }
        return std::sqrt(result / total_area);
    }

    double ErrorNormTorusCalculator::getL2NormFluidAccelWithRecoveryUsingGaussExactValues(ModelPart &r_model_part)
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
                CalculateMaterialAcceleration(gauss_points_coordinates[g], exact_material_acceleration_gauss_point);

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
                }
            }
            total_area += r_geometry.Area();
        }
        return std::sqrt(result / total_area);
    }
}
