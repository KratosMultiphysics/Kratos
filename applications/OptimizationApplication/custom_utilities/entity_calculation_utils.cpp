//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main author:     Reza Najian Asl
//                   Suneth Warnakulasuriya
//

// System includes

// Project includes

#include "geometries/tetrahedra_3d_4.h"
#include "geometries/pyramid_3d_5.h"

// Application includes

// Include base h
#include "entity_calculation_utils.h"

namespace Kratos {

Geometry<Node>::Pointer EntityCalculationUtils::CreateSolidGeometry(const GeometryType& rSurfaceGeometry)
{
    // add the surface nodes
    GeometryType::PointsArrayType nodes(rSurfaceGeometry.ptr_begin(), rSurfaceGeometry.ptr_end());

    // add the dummy additional point for the extrusion
    nodes.push_back(Kratos::make_intrusive<Node>(0, 0.0, 0.0, 0.0));

    // create the new solid geometry from the surface geometry
    switch (rSurfaceGeometry.GetGeometryType())
    {
        case GeometryData::KratosGeometryType::Kratos_Triangle3D3:
            return Kratos::make_shared<Tetrahedra3D4<Node>>(nodes);
        case GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4:
            return Kratos::make_shared<Pyramid3D5<Node>>(nodes);
        default:
            KRATOS_ERROR << "Unsupported geometry types. "
                            "CalculateSurfaceElementShapeDerivatives only "
                            "supports surface geometries of following types:"
                            "\n\tTriagle3D3"
                            "\n\tQuadrilateral3D4.\n";
    }

    return nullptr;
}

void EntityCalculationUtils::CalculateElementGaussPointData(
    Vector& rGaussWeights,
    Matrix& rShapeFunctionValues,
    const GeometryType& rGeometry,
    const IntegrationMethodType& rIntegrationMethod)
{
    rShapeFunctionValues = rGeometry.ShapeFunctionsValues(rIntegrationMethod);

    const auto& integration_points = rGeometry.IntegrationPoints(rIntegrationMethod);
    const IndexType number_of_integration_points = integration_points.size();

    if (rGaussWeights.size() != number_of_integration_points) {
        rGaussWeights.resize(number_of_integration_points, false);
    }

    for (IndexType g = 0; g < number_of_integration_points; ++g) {
        const double detJ0 = rGeometry.DeterminantOfJacobian(g, rIntegrationMethod);
        rGaussWeights[g] = detJ0 * integration_points[g].Weight();
    }
}

template<class TMatrixType>
void EntityCalculationUtils::CalculateSurfaceElementShapeDerivatives(
    TMatrixType& rOutput,
    GeometryType& rSolidGeometry,
    const GeometryType& rSurfaceGeometry,
    const IntegrationMethodType& rIntegrationMethod,
    const IndexType PointNumber)
{
    // now calculate the out of surface point
    const auto& integration_points = rSurfaceGeometry.IntegrationPoints(rIntegrationMethod);
    const auto& r_surface_normal = rSurfaceGeometry.UnitNormal(PointNumber);
    const auto length = rSurfaceGeometry.Length();

    const auto& surf_gp_local_pt = integration_points[PointNumber].Coordinates();
    Point surf_gp_global_pt;
    rSurfaceGeometry.GlobalCoordinates(surf_gp_global_pt, surf_gp_local_pt);

    // move the last point the solid geometry properly
    rSolidGeometry.back().Coordinates() = surf_gp_local_pt + r_surface_normal * length;

    Point elem_surf_gp_local_pt;
    rSolidGeometry.PointLocalCoordinates(elem_surf_gp_local_pt, surf_gp_global_pt);

    Matrix dN_de;
    rSolidGeometry.ShapeFunctionsLocalGradients(dN_de, elem_surf_gp_local_pt);
    Matrix inv_j0;
    rSolidGeometry.InverseOfJacobian(inv_j0, elem_surf_gp_local_pt);

    const Matrix& solid_dN_dX = prod(dN_de, inv_j0);

    if (rOutput.size1() != rSurfaceGeometry.size()) {
        rOutput.resize(rSurfaceGeometry.size(), 3, false);
    }

    std::copy(&solid_dN_dX.data()[0], &solid_dN_dX.data()[0] + rSurfaceGeometry.size() * 3, &rOutput.data()[0]);
}

// template instantiations
template void EntityCalculationUtils::CalculateSurfaceElementShapeDerivatives(Matrix&, GeometryType&, const GeometryType&, const IntegrationMethodType&, const IndexType);

} // namespace Kratos