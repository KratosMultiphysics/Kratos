// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Gennady Markelov
//

#pragma once

// Project includes

// Application includes
#include "custom_elements/stress_state_policy.h"
#include "custom_retention/retention_law.h"
#include "custom_utilities/element_utilities.hpp"
#include "custom_utilities/transport_equation_utilities.hpp"
#include "geo_mechanics_application_variables.h"
#include "geometries/geometry.h"
#include "includes/element.h"
#include "includes/smart_pointers.h"
#include "includes/variables.h"

namespace Kratos
{

class GeoEquationOfMotionUtilities
{
public:
    static Matrix CalculateMassMatrix(const Geometry<Node>&                 rGeom,
                                      const GeometryData::IntegrationMethod IntegrationMethod,
                                      const Vector&                         rSolidDensities,
                                      const Vector& rIntegrationCoefficients)
    {
        const SizeType dimension          = rGeom.WorkingSpaceDimension();
        const SizeType number_U_nodes     = rGeom.PointsNumber();
        const SizeType block_element_size = number_U_nodes * dimension;
        const Geometry<Node>::IntegrationPointsArrayType& integration_points =
            rGeom.IntegrationPoints(IntegrationMethod);
        const Matrix& Nu_container = rGeom.ShapeFunctionsValues(IntegrationMethod);

        Matrix Nu                 = ZeroMatrix(dimension, number_U_nodes * dimension);
        Matrix aux_density_matrix = ZeroMatrix(dimension, number_U_nodes * dimension);
        Matrix density_matrix     = ZeroMatrix(dimension, dimension);
        Matrix mass_matrix        = ZeroMatrix(block_element_size, block_element_size);

        for (unsigned int g_point = 0; g_point < integration_points.size(); ++g_point) {
            GeoElementUtilities::AssembleDensityMatrix(density_matrix, rSolidDensities(g_point));
            GeoElementUtilities::CalculateNuMatrix(dimension, number_U_nodes, Nu, Nu_container, g_point);
            noalias(aux_density_matrix) = prod(density_matrix, Nu);
            mass_matrix += prod(trans(Nu), aux_density_matrix) * rIntegrationCoefficients(g_point);
        }
        return mass_matrix;
    }

    static Vector CalculateIntegrationCoefficientsInitialConfiguration(const Geometry<Node>& rGeom,
                                                                       const GeometryData::IntegrationMethod IntegrationMethod,
                                                                       const StressStatePolicy& rStressStatePolicy)
    {
        const Geometry<Node>::IntegrationPointsArrayType& integration_points =
            rGeom.IntegrationPoints(IntegrationMethod);
        const std::size_t number_G_points = integration_points.size();

        Vector integration_coefficient_initial_configuration(number_G_points);
        Matrix J0;
        Matrix inv_J0;
        for (unsigned int g_point = 0; g_point < number_G_points; ++g_point) {
            GeometryUtils::JacobianOnInitialConfiguration(rGeom, integration_points[g_point], J0);
            const Matrix& dN_De = rGeom.ShapeFunctionsLocalGradients(IntegrationMethod)[g_point];
            double        det_J_initial_configuration;
            MathUtils<double>::InvertMatrix(J0, inv_J0, det_J_initial_configuration);
            Matrix dNu_dX_initial_configuration;
            GeometryUtils::ShapeFunctionsGradients(dN_De, inv_J0, dNu_dX_initial_configuration);

            integration_coefficient_initial_configuration(g_point) =
                rStressStatePolicy.CalculateIntegrationCoefficient(
                    integration_points[g_point], det_J_initial_configuration, rGeom);
        }
        return integration_coefficient_initial_configuration;
    }

}; /* Class GeoTransportEquationUtilities*/
} /* namespace Kratos.*/
