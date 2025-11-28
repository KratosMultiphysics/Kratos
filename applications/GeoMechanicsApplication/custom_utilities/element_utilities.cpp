// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//

#include "element_utilities.hpp"

#include "custom_elements/interface_element.h"

#include <cstddef>

namespace Kratos
{

void GeoElementUtilities::FillArray1dOutput(array_1d<double, 3>& rOutputValue, const array_1d<double, 2>& ComputedValue)
{
    rOutputValue[0] = ComputedValue[0];
    rOutputValue[1] = ComputedValue[1];
    rOutputValue[2] = 0.0;
}

void GeoElementUtilities::FillArray1dOutput(array_1d<double, 3>& rOutputValue, const array_1d<double, 3>& ComputedValue)
{
    rOutputValue[0] = ComputedValue[0];
    rOutputValue[1] = ComputedValue[1];
    rOutputValue[2] = ComputedValue[2];
}

void GeoElementUtilities::FillPermeabilityMatrix(BoundedMatrix<double, 1, 1>&   rPermeabilityMatrix,
                                                 const Element::PropertiesType& Prop)
{
    // 1D
    if (Prop[RETENTION_LAW] == "PressureFilterLaw") {
        const auto equivalent_radius_square = Prop[CROSS_AREA] / Globals::Pi;
        rPermeabilityMatrix(0, 0)           = equivalent_radius_square * 0.125;
    } else {
        rPermeabilityMatrix(0, 0) = Prop[PERMEABILITY_XX];
    }
}

void GeoElementUtilities::FillPermeabilityMatrix(BoundedMatrix<double, 2, 2>&   rPermeabilityMatrix,
                                                 const Element::PropertiesType& Prop)
{
    // 2D
    rPermeabilityMatrix(0, 0) = Prop[PERMEABILITY_XX];
    rPermeabilityMatrix(1, 1) = Prop[PERMEABILITY_YY];

    rPermeabilityMatrix(0, 1) = Prop[PERMEABILITY_XY];
    rPermeabilityMatrix(1, 0) = rPermeabilityMatrix(0, 1);
}

void GeoElementUtilities::FillPermeabilityMatrix(BoundedMatrix<double, 3, 3>&   rPermeabilityMatrix,
                                                 const Element::PropertiesType& Prop)
{
    // 3D
    rPermeabilityMatrix(0, 0) = Prop[PERMEABILITY_XX];
    rPermeabilityMatrix(1, 1) = Prop[PERMEABILITY_YY];
    rPermeabilityMatrix(2, 2) = Prop[PERMEABILITY_ZZ];

    rPermeabilityMatrix(0, 1) = Prop[PERMEABILITY_XY];
    rPermeabilityMatrix(1, 0) = rPermeabilityMatrix(0, 1);

    rPermeabilityMatrix(1, 2) = Prop[PERMEABILITY_YZ];
    rPermeabilityMatrix(2, 1) = rPermeabilityMatrix(1, 2);

    rPermeabilityMatrix(2, 0) = Prop[PERMEABILITY_ZX];
    rPermeabilityMatrix(0, 2) = rPermeabilityMatrix(2, 0);
}

Matrix GeoElementUtilities::FillPermeabilityMatrix(const Element::PropertiesType& Prop, std::size_t Dimension)
{
    switch (Dimension) {
    case 1: {
        BoundedMatrix<double, 1, 1> result;
        FillPermeabilityMatrix(result, Prop);
        return result;
    }
    case 2: {
        BoundedMatrix<double, 2, 2> result;
        FillPermeabilityMatrix(result, Prop);
        return result;
    }
    case 3: {
        BoundedMatrix<double, 3, 3> result;
        FillPermeabilityMatrix(result, Prop);
        return result;
    }
    default:
        KRATOS_ERROR << "Dimension " << Dimension << " is not supported" << std::endl;
    }
}

double GeoElementUtilities::CalculateRadius(const Vector& rN, const GeometryType& rGeometry)
{
    auto radius = 0.0;

    for (unsigned int node = 0; node < rGeometry.size(); ++node) {
        // Displacement from the reference to the current configuration
        const array_1d<double, 3>& current_position = rGeometry[node].Coordinates();
        radius += current_position[0] * rN[node];
    }

    return radius;
}

double GeoElementUtilities::CalculateAxisymmetricCircumference(const Vector& rN, const GeometryType& rGeometry)
{
    return 2.0 * Globals::Pi * CalculateRadius(rN, rGeometry);
}

Vector GeoElementUtilities::CalculateNodalHydraulicHeadFromWaterPressures(const GeometryType& rGeom,
                                                                          const Properties&   rProp)
{
    constexpr auto numerical_limit = std::numeric_limits<double>::epsilon();
    // Defining necessary variables

    Vector nodal_hydraulic_heads(rGeom.PointsNumber());
    for (unsigned int node = 0; node < rGeom.PointsNumber(); ++node) {
        array_1d<double, 3> node_volume_acceleration;
        noalias(node_volume_acceleration) = rGeom[node].FastGetSolutionStepValue(VOLUME_ACCELERATION, 0);
        const auto g = norm_2(node_volume_acceleration);
        if (g > numerical_limit) {
            const auto fluid_weight = g * rProp[DENSITY_WATER];

            array_1d<double, 3> node_coordinates;
            noalias(node_coordinates) = rGeom[node].Coordinates();
            array_1d<double, 3> node_volume_acceleration_unit_vector;
            noalias(node_volume_acceleration_unit_vector) = node_volume_acceleration / g;

            const auto water_pressure = rGeom[node].FastGetSolutionStepValue(WATER_PRESSURE);
            nodal_hydraulic_heads[node] = -inner_prod(node_coordinates, node_volume_acceleration_unit_vector) -
                                          PORE_PRESSURE_SIGN_FACTOR * water_pressure / fluid_weight;
        } else {
            nodal_hydraulic_heads[node] = 0.0;
        }
    }
    return nodal_hydraulic_heads;
}

Geo::IntegrationPointVectorType GeoElementUtilities::GetIntegrationPointsOf(const Element& rElement)
{
    auto p_interface_element = dynamic_cast<const InterfaceElement*>(&rElement);
    return p_interface_element
               ? p_interface_element->GetIntegrationPoints()
               : rElement.GetGeometry().IntegrationPoints(rElement.GetIntegrationMethod());
}

std::vector<Vector> GeoElementUtilities::EvaluateShapeFunctionsAtIntegrationPoints(
    const Geo::IntegrationPointVectorType& rIntegrationPoints, const Geometry<Node>& rGeometry)
{
    auto evaluate_shape_function_values = [&rGeometry](const auto& rIntegrationPoint) {
        auto result = Vector{};
        rGeometry.ShapeFunctionsValues(result, rIntegrationPoint);
        return result;
    };

    auto result = std::vector<Vector>{};
    result.reserve(rIntegrationPoints.size());
    std::ranges::transform(rIntegrationPoints, std::back_inserter(result), evaluate_shape_function_values);

    return result;
}

} // namespace Kratos