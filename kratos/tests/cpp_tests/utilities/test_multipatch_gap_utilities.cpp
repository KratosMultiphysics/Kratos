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

#include <array>
#include <vector>
#include "testing/testing.h"
#include "containers/model.h"
#include "geometries/quadrilateral_2d_4.h"
#include "processes/structured_mesh_generator_process.h"
#include "utilities/multipatch_gap_utilities.h"

namespace Kratos::Testing {
namespace {
void Mesh(ModelPart& rPart, double X, double Y, double Size, double Angle,
          int Divisions, int Start)
{
    const double c = std::cos(Angle), s = std::sin(Angle);
    Geometry<Node>::PointsArrayType points;
    for (const auto& p:std::vector<std::array<double, 2>>{{0, 0}, {0, 1}, {1, 1}, {1, 0}}) {
        const double x = Size*(p[0]-0.5), y = Size*(p[1]-0.5);
        points.push_back(Kratos::make_intrusive<Node>(points.size()+1, X+Size/2+c*x-s*y, Y+Size/2+s*x+c*y, 0.0));
    }
    Quadrilateral2D4<Node> domain(points);
    Parameters settings(R"({"element_name":"Element2D4N", "number_of_divisions":[1,1,1],
        "create_skin_sub_model_part":false,"start_node_id":1,"start_element_id":1})");
    settings["number_of_divisions"][0].SetInt(Divisions);
    settings["number_of_divisions"][1].SetInt(Divisions);
    settings["start_node_id"].SetInt(Start);
    settings["start_element_id"].SetInt(Start);
    StructuredMeshGeneratorProcess(domain, rPart, settings).Execute();
}
}
KRATOS_TEST_CASE_IN_SUITE(FemMultipatchGapCoverageAndQ1Extension, KratosCoreFastSuite)
{
    for (const double angle:{0.0, 0.4014257279586958, 0.7853981633974483}) {
        Model model;
        auto& root = model.CreateModelPart("Fluid");
        auto& bg = root.CreateSubModelPart("Background");
        auto& inner = root.CreateSubModelPart("Inner");
        Mesh(bg, 0.0, 0.0, 1.0, 0.0, 12, 1);
        Mesh(inner, 0.33, 0.37, 0.32, angle, 5, 1000);
        const auto nodes = root.NumberOfNodes();
        const auto owners = FemMultipatchGapUtilities::UpdateActivation(bg, {&inner}, 1e-10);
        const auto data = FemMultipatchGapUtilities::Build(bg, inner, 1e-10);
        KRATOS_EXPECT_GT(data.Cells.size(), 0);
        double active_area = 0.0, gap_area = 0.0, perimeter = 0.0;
        for (const auto& e:bg.Elements()) if (e.Is(ACTIVE)) active_area+=e.GetGeometry().Area();
        for (const auto& segment:data.Physical) perimeter+=norm_2(segment.End-segment.Start);
        KRATOS_EXPECT_NEAR(perimeter, 4*0.32, 1e-10);
        for (const auto& cell:data.Cells) {
            const auto& g = *cell.pQuadrature;
            const auto method = GeometryData::IntegrationMethod::GI_GAUSS_2;
            const auto& ips = g.IntegrationPoints(method);
            const auto& n = g.ShapeFunctionsValues(method);
            Geometry<Node>::ShapeFunctionsGradientsType gradients;
            Vector determinants;
            g.ShapeFunctionsIntegrationPointsGradients(gradients, determinants, method);
            for (std::size_t k = 0; k<ips.size(); ++k) {
                Geometry<Node>::CoordinatesArrayType p = ZeroVector(3);
                g.GlobalCoordinates(p, ips[k]);
                gap_area+=ips[k].Weight()*determinants[k];
                double one = 0, x = 0, y = 0, xy = 0, gx = 0, gy = 0, xyx = 0, xyy = 0;
                for (std::size_t i = 0; i<4; ++i) {
                    one+=n(k, i); x+=n(k, i)*g[i].X(); y+=n(k, i)*g[i].Y();
                    xy+=n(k, i)*g[i].X()*g[i].Y();
                    gx+=gradients[k](i, 0); gy+=gradients[k](i, 1);
                    xyx+=gradients[k](i, 0)*g[i].X()*g[i].Y();
                    xyy+=gradients[k](i, 1)*g[i].X()*g[i].Y();
                }
                KRATOS_EXPECT_NEAR(one, 1.0, 1e-12);
                KRATOS_EXPECT_NEAR(x, p[0], 1e-12); KRATOS_EXPECT_NEAR(y, p[1], 1e-12);
                KRATOS_EXPECT_NEAR(xy, p[0]*p[1], 1e-12);
                KRATOS_EXPECT_NEAR(gx, 0.0, 1e-12); KRATOS_EXPECT_NEAR(gy, 0.0, 1e-12);
                KRATOS_EXPECT_NEAR(xyx, p[1], 1e-12); KRATOS_EXPECT_NEAR(xyy, p[0], 1e-12);
            }
        }
        KRATOS_EXPECT_NEAR(active_area+gap_area+0.32*0.32, 1.0, 1e-10);
        KRATOS_EXPECT_EQ(nodes, root.NumberOfNodes());
        KRATOS_EXPECT_EQ(bg.NumberOfElements(), 144);
        KRATOS_EXPECT_GT(owners.size(), 0);
    }
}
KRATOS_TEST_CASE_IN_SUITE(FemMultipatchGapFittedLimitAndInvalidOwnership, KratosCoreFastSuite)
{
    Model model;
    auto& root = model.CreateModelPart("Fluid");
    auto& bg = root.CreateSubModelPart("Background");
    auto& inner = root.CreateSubModelPart("Inner");
    Mesh(bg, 0, 0, 1, 0, 8, 1); Mesh(inner, 0.25, 0.25, 0.5, 0, 3, 1000);
    FemMultipatchGapUtilities::UpdateActivation(bg, {&inner}, 1e-10);
    const auto data = FemMultipatchGapUtilities::Build(bg, inner, 1e-10);
    KRATOS_EXPECT_EQ(data.Cells.size(), 0);
    double length = 0.0;
    for (const auto& s:data.Physical) length+=norm_2(s.End-s.Start);
    KRATOS_EXPECT_NEAR(length, 2.0, 1e-10);
    auto& other = root.CreateSubModelPart("Other");
    Mesh(other, 0.4, 0.4, 0.2, 0, 2, 2000);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        FemMultipatchGapUtilities::UpdateActivation(bg, {&inner, &other}, 1e-10),"overlap");
}
} // namespace Kratos::Testing
