//    Kratos Multiphysics
//  License: BSD License (kratos/license.txt)

#include "testing/testing.h"
#include "containers/model.h"
#include "custom_elements/isogeometric_beam_element.h"
#include "geometries/nurbs_curve_geometry.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(IsogeometricBeamReferenceFrameUsesInitialCoordinates, KratosIgaFastSuite)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("beam");
    Geometry<Node>::PointsArrayType points;
    for (std::size_t i = 0; i < 3; ++i) {
        points.push_back(r_model_part.CreateNewNode(i + 1, static_cast<double>(i), 0.0, 0.0));
    }
    Vector knots(4);
    knots[0] = knots[1] = 0.0;
    knots[2] = knots[3] = 1.0;
    NurbsCurveGeometry<3, PointerVector<Node>> curve(points, 2, knots);
    Geometry<Node>::IntegrationPointsArrayType integration_points;
    integration_points.emplace_back(0.37, 0.0, 0.0, 1.0);
    Geometry<Node>::GeometriesArrayType quadrature_points;
    auto integration_info = curve.GetDefaultIntegrationInfo();
    curve.CreateQuadraturePointGeometries(quadrature_points, 3, integration_points, integration_info);
    auto p_properties = r_model_part.CreateNewProperties(1);
    array_1d<double, 3> tangent = ZeroVector(3);
    tangent[0] = 1.0;
    array_1d<double, 3> normal = ZeroVector(3);
    normal[2] = 1.0;
    p_properties->SetValue(T_0, tangent);
    p_properties->SetValue(N_0, normal);
    Matrix orientation = ZeroMatrix(2, 4);
    orientation(1, 0) = 1.0;
    orientation(0, 3) = orientation(1, 3) = 1.0;
    p_properties->SetValue(LOCAL_AXIS_ORIENTATION, orientation);
    IsogeometricBeamElement element(1, quadrature_points(0), p_properties);
    // The query must not require material initialization or solution-step data.
    r_model_part.GetNode(3).Y() = 10.0;
    std::vector<Matrix> frames;
    element.CalculateOnIntegrationPoints(LOCAL_AXES_MATRIX, frames, r_model_part.GetProcessInfo());
    Matrix expected = ZeroMatrix(3, 3);
    expected(0, 0) = 1.0;
    expected(1, 2) = 1.0;
    expected(2, 1) = -1.0;
    KRATOS_EXPECT_EQ(frames.size(), 1);
    KRATOS_EXPECT_MATRIX_NEAR(frames[0], expected, 1e-12);
    KRATOS_EXPECT_NEAR(r_model_part.GetNode(3).Y(), 10.0, 1e-12);
}

} // namespace Kratos::Testing
