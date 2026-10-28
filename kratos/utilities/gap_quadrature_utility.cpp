//     |  /           |
//     ' /   __| _` | __|  _ \   __|
//     . \  |   (   | |   (   |\__ `
//    _|\_\_|  \__,_|\__|\___/ ____/
//                    Multi-Physics
//
//   License:         BSD License
//                    Kratos default license: kratos/license.txt
//
//   Main author:     Nicolò Antonelli
//

#include "utilities/gap_quadrature_utility.h"
#include "geometries/geometry_shape_function_container.h"
#include "geometries/quadrature_point_geometry.h"

namespace Kratos {
Geometry<Node>::Pointer GapQuadratureUtility::Create(
    Geometry<Node>& rParent,
    const Geometry<Node>::IntegrationPointsArrayType& rPhysicalPoints)
{
    KRATOS_TRY
    using QPGeometry = QuadraturePointGeometry<Node, 2, 2>;
    using Container = GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>;
    Geometry<Node>::IntegrationPointsArrayType points;
    Matrix values(rPhysicalPoints.size(), rParent.size());
    Geometry<Node>::ShapeFunctionsGradientsType gradients(rPhysicalPoints.size());
    for (std::size_t g = 0; g < rPhysicalPoints.size(); ++g) {
        const auto& r_point = rPhysicalPoints[g];
        Geometry<Node>::CoordinatesArrayType local = ZeroVector(3);
        rParent.PointLocalCoordinates(local, Point(r_point.X(), r_point.Y(), 0.0));
        const double det_j = rParent.DeterminantOfJacobian(local);
        KRATOS_ERROR_IF(det_j <= 0.0 || r_point.Weight() <= 0.0)
            << "Gap quadrature requires positive parent Jacobians and weights." << std::endl;
        Vector n;
        rParent.ShapeFunctionsValues(n, local);
        for (std::size_t i = 0; i < rParent.size(); ++i) values(g, i) = n[i];
        rParent.ShapeFunctionsLocalGradients(gradients[g], local);
        // FluidElement multiplies the stored weight by the parent Jacobian.
        points.emplace_back(local[0], local[1], local[2], r_point.Weight() / det_j);
    }
    QPGeometry::IntegrationPointsContainerType all_points;
    QPGeometry::ShapeFunctionsValuesContainerType all_values;
    QPGeometry::ShapeFunctionsLocalGradientsContainerType all_gradients;
    for (const auto method : {GeometryData::IntegrationMethod::GI_GAUSS_1,
                              GeometryData::IntegrationMethod::GI_GAUSS_2}) {
        all_points[static_cast<int>(method)] = points;
        all_values[static_cast<int>(method)] = values;
        all_gradients[static_cast<int>(method)] = gradients;
    }
    Container container(GeometryData::IntegrationMethod::GI_GAUSS_2,
                        all_points, all_values, all_gradients);
    return Kratos::make_shared<QPGeometry>(rParent.Points(), container, &rParent);
    KRATOS_CATCH("")
}
} // namespace Kratos
