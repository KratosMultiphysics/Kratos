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
        const double equivalent_radius_square = Prop[CROSS_AREA] / Globals::Pi;
        rPermeabilityMatrix(0, 0)             = equivalent_radius_square * 0.125;
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

void GeoElementUtilities::InvertMatrix2(BoundedMatrix<double, 2, 2>&       rInvertedMatrix,
                                        const BoundedMatrix<double, 2, 2>& InputMatrix,
                                        double&                            InputMatrixDet)
{
    KRATOS_TRY

    const double numerical_limit = std::numeric_limits<double>::epsilon();

    InputMatrixDet = InputMatrix(0, 0) * InputMatrix(1, 1) - InputMatrix(0, 1) * InputMatrix(1, 0);

    if (InputMatrixDet < numerical_limit) {
        KRATOS_ERROR << "determinant zero or negative" << std::endl;
    }

    rInvertedMatrix(0, 0) = InputMatrix(1, 1) / InputMatrixDet;
    rInvertedMatrix(0, 1) = -InputMatrix(0, 1) / InputMatrixDet;
    rInvertedMatrix(1, 0) = -InputMatrix(1, 0) / InputMatrixDet;
    rInvertedMatrix(1, 1) = InputMatrix(0, 0) / InputMatrixDet;

    KRATOS_CATCH("")
}

void GeoElementUtilities::InvertMatrix2(BoundedMatrix<double, 2, 2>&       rInvertedMatrix,
                                        const BoundedMatrix<double, 2, 2>& InputMatrix)
{
    KRATOS_TRY

    double InputMatrixDet;

    InvertMatrix2(rInvertedMatrix, InputMatrix, InputMatrixDet);

    KRATOS_CATCH("")
}

void GeoElementUtilities::CalculateNewtonCotesLocalShapeFunctionsGradients(BoundedMatrix<double, 2, 2>& DN_DeContainer)
{
    // Line 2-noded
    //  nodes:
    //  0------------1
    const unsigned int        NumNodes = 2;
    const std::vector<double> Xi{-1.0, 1.0};

    noalias(DN_DeContainer) = ZeroMatrix(NumNodes, NumNodes);

    for (unsigned int integrationPoint = 0; integrationPoint < NumNodes; ++integrationPoint) {
        DN_DeContainer(integrationPoint, 0) = -0.5;
        DN_DeContainer(integrationPoint, 1) = 0.5;
    }
}

void GeoElementUtilities::CalculateNewtonCotesLocalShapeFunctionsGradients(BoundedMatrix<double, 3, 3>& DN_DeContainer)
{
    // Line 3-noded
    //  nodes:
    //  0------2------1
    const unsigned int        NumNodes = 3;
    const std::vector<double> Xi{-1.0, 0.0, 1.0};

    noalias(DN_DeContainer) = ZeroMatrix(NumNodes, NumNodes);

    for (unsigned int integrationPoint = 0; integrationPoint < NumNodes; ++integrationPoint) {
        DN_DeContainer(integrationPoint, 0) = Xi[integrationPoint] - 0.5;
        DN_DeContainer(integrationPoint, 1) = Xi[integrationPoint] + 0.5;
        DN_DeContainer(integrationPoint, 2) = -Xi[integrationPoint] * 2.0;
    }
}

void GeoElementUtilities::CalculateNewtonCotesShapeFunctions(BoundedMatrix<double, 3, 3>& NContainer)
{
    // Line 3-noded
    //  nodes:
    //  0------2------1
    const unsigned int NumNodes = 3;

    for (unsigned int integrationPoint = 0; integrationPoint < NumNodes; ++integrationPoint) {
        for (unsigned int node = 0; node < NumNodes; ++node) {
            NContainer(integrationPoint, node) = (integrationPoint == node ? 1.0 : 0.0);
        }
    }
}

void GeoElementUtilities::CalculateEquallyDistributedPointsLineShapeFunctions3N(Matrix& NContainer)
{
    // Line 3-noded
    //  nodes:
    //  0------2------1
    const unsigned int        NumNodes = 3;
    const std::vector<double> Xi{-1.0 / 3.0, 1.0 / 3.0};

    if (NContainer.size1() != Xi.size() || NContainer.size2() != NumNodes)
        NContainer.resize(Xi.size(), NumNodes, false);

    for (unsigned int integrationPoint = 0; integrationPoint < Xi.size(); ++integrationPoint) {
        const double& X                 = Xi[integrationPoint];
        NContainer(integrationPoint, 0) = -0.5 * (1.0 - X) * X;
        NContainer(integrationPoint, 1) = 0.5 * (1.0 + X) * X;
        NContainer(integrationPoint, 2) = (1.0 + X) * (1.0 - X);
    }
}

void GeoElementUtilities::CalculateEquallyDistributedPointsLineGradientShapeFunctions3N(GeometryData::ShapeFunctionsGradientsType& DN_DeContainer)
{
    // Line 3-noded
    //  nodes:
    //  0------2------1
    const unsigned int        NumNodes = 3;
    const std::vector<double> Xi{-1.0 / 3.0, 1.0 / 3.0};

    if (DN_DeContainer.size() != Xi.size()) DN_DeContainer.resize(Xi.size());

    for (unsigned int integrationPoint = 0; integrationPoint < Xi.size(); ++integrationPoint) {
        if (DN_DeContainer[integrationPoint].size1() != NumNodes ||
            DN_DeContainer[integrationPoint].size2() != 1)
            DN_DeContainer[integrationPoint].resize(NumNodes, 1);

        DN_DeContainer[integrationPoint](0, 0) = Xi[integrationPoint] - 0.5;
        DN_DeContainer[integrationPoint](1, 0) = Xi[integrationPoint] + 0.5;
        DN_DeContainer[integrationPoint](2, 0) = -Xi[integrationPoint] * 2.0;
    }
}

double GeoElementUtilities::CalculateRadius(const Vector& N, const GeometryType& Geom)
{
    double Radius = 0.0;

    for (unsigned int iNode = 0; iNode < Geom.size(); ++iNode) {
        // Displacement from the reference to the current configuration
        const array_1d<double, 3>& CurrentPosition = Geom[iNode].Coordinates();
        Radius += CurrentPosition[0] * N[iNode];
    }

    return Radius;
}

double GeoElementUtilities::CalculateAxisymmetricCircumference(const Vector& N, const GeometryType& Geom)
{
    const double Radius        = CalculateRadius(N, Geom);
    const double Circumference = 2.0 * Globals::Pi * Radius;
    return Circumference;
}

Vector GeoElementUtilities::CalculateNodalHydraulicHeadFromWaterPressures(const GeometryType& rGeom,
                                                                          const Properties&   rProp)
{
    const auto NumericalLimit = std::numeric_limits<double>::epsilon();
    // Defining necessary variables

    Vector nodal_hydraulic_heads(rGeom.PointsNumber());
    for (unsigned int node = 0; node < rGeom.PointsNumber(); ++node) {
        array_1d<double, 3> node_volume_acceleration;
        noalias(node_volume_acceleration) = rGeom[node].FastGetSolutionStepValue(VOLUME_ACCELERATION, 0);
        const double g = norm_2(node_volume_acceleration);
        if (g > NumericalLimit) {
            const double FluidWeight = g * rProp[DENSITY_WATER];

            array_1d<double, 3> NodeCoordinates;
            noalias(NodeCoordinates) = rGeom[node].Coordinates();
            array_1d<double, 3> NodeVolumeAccelerationUnitVector;
            noalias(NodeVolumeAccelerationUnitVector) = node_volume_acceleration / g;

            const double WaterPressure = rGeom[node].FastGetSolutionStepValue(WATER_PRESSURE);
            nodal_hydraulic_heads[node] = -inner_prod(NodeCoordinates, NodeVolumeAccelerationUnitVector) -
                                          PORE_PRESSURE_SIGN_FACTOR * WaterPressure / FluidWeight;
        } else {
            nodal_hydraulic_heads[node] = 0.0;
        }
    }
    return nodal_hydraulic_heads;
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
    std::transform(rIntegrationPoints.begin(), rIntegrationPoints.end(), std::back_inserter(result),
                   evaluate_shape_function_values);

    return result;
}

} // namespace Kratos