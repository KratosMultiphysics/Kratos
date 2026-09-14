//    Kratos Multiphysics
//
//  License: BSD License
//           Kratos default license: kratos/license.txt

#include "mapping_fast_suite.h"
#include "custom_mappers/iga_beam_mapper.h"
#include "factories/mapper_factory.h"
#include "mappers/mapper_define.h"

namespace Kratos::Testing
{
namespace
{
using MapperType = IgaBeamMapper<MapperDefinitions::SparseSpaceType, MapperDefinitions::DenseSpaceType>;
using FactoryType = MapperFactory<MapperDefinitions::SparseSpaceType, MapperDefinitions::DenseSpaceType>;

// Test double for the element-side reference frame protocol. The actual beam
// implementation is tested in IgaApplication.
class ReferenceFrameElement : public Element
{
public:
    using Element::Element;
    Element::Pointer Create(IndexType Id, GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<ReferenceFrameElement>(Id, pGeometry, pProperties);
    }

    void CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable,
        std::vector<Matrix>& rOutput, const ProcessInfo& rProcessInfo) override
    {
        array_1d<double, 3> tangent = ZeroVector(3);
        const auto& r_dn = GetGeometry().ShapeFunctionDerivatives(1, 0);
        for (std::size_t i = 0; i < GetGeometry().size(); ++i) {
            tangent += r_dn(i, 0) * GetGeometry()[i].GetInitialPosition();
        }
        tangent /= norm_2(tangent);
        rOutput.assign(1, ZeroMatrix(3, 3));
        for (std::size_t i = 0; i < 3; ++i) rOutput[0](0, i) = tangent[i];
        rOutput[0](1, 0) = -tangent[1];
        rOutput[0](1, 1) = tangent[0];
        rOutput[0](2, 2) = 1.0;
    }
};

MapperType::NurbsCurveType::Pointer AddBeamCurve(ModelPart& rBeam, const bool Curved = false)
{
    Geometry<Node>::PointsArrayType points;
    const auto first_node_id = rBeam.NumberOfNodes() + 1;
    for (std::size_t i = 0; i < 3; ++i) {
        points.push_back(rBeam.CreateNewNode(first_node_id + i,
            Curved ? (i == 2 ? 0.0 : 1.0) : static_cast<double>(i),
            Curved && i > 0 ? 1.0 : 0.0, 0.0));
    }
    Vector knots(4);
    knots[0] = 0.0;
    knots[1] = 0.0;
    knots[2] = 1.0;
    knots[3] = 1.0;
    Vector weights(3, 1.0);
    if (Curved) weights[1] = std::sqrt(0.5);
    auto p_curve = Kratos::make_shared<MapperType::NurbsCurveType>(points, 2, knots, weights);
    Geometry<Node>::GeometriesArrayType quadrature_points;
    auto integration_info = p_curve->GetDefaultIntegrationInfo();
    Geometry<Node>::IntegrationPointsArrayType integration_points;
    p_curve->CreateIntegrationPoints(integration_points, integration_info);
    p_curve->CreateQuadraturePointGeometries(quadrature_points, 2, integration_points, integration_info);
    for (std::size_t i = 0; i < quadrature_points.size(); ++i) {
        rBeam.AddElement(Kratos::make_intrusive<ReferenceFrameElement>(
            rBeam.NumberOfElements() + 1, quadrature_points(i)));
    }
    return p_curve;
}

void AddSurfaceNode(ModelPart& rSurface)
{
    rSurface.AddNodalSolutionStepVariable(DISPLACEMENT);
    rSurface.CreateNewNode(1, 1.0, 1.0, 0.0);
}
} // namespace

KRATOS_TEST_CASE_IN_SUITE(IgaBeamMapperCachesUndeformedAttachments, KratosMappingApplicationSerialTestSuite)
{
    Model model;
    auto& r_beam = model.CreateModelPart("beam");
    auto& r_surface = model.CreateModelPart("surface");
    r_beam.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_beam.AddNodalSolutionStepVariable(ROTATION);
    const auto p_curve = AddBeamCurve(r_beam);
    AddSurfaceNode(r_surface);
    for (auto& r_node : r_beam.Nodes()) r_node.Y() += 10.0;
    r_surface.GetNode(1).Y() += 20.0;

    MapperType mapper(r_beam, r_surface, Parameters("{}"));
    const auto& r_attachment = mapper.GetReferenceAttachments().at(0);
    KRATOS_EXPECT_NEAR(r_attachment.Parameter, 0.5, 1e-10);
    KRATOS_EXPECT_NEAR(r_attachment.NormalOffset, 1.0, 1e-10);
    KRATOS_EXPECT_NEAR(r_attachment.BinormalOffset, 0.0, 1e-10);
    KRATOS_EXPECT_EQ(r_attachment.ControlPoints.size(), 3);
    KRATOS_EXPECT_EQ(r_attachment.ControlPoints[0].get(), r_beam.pGetNode(1).get());
    KRATOS_EXPECT_NEAR(r_attachment.ShapeFunctions[0], 0.25, 1e-10);
    KRATOS_EXPECT_NEAR(r_attachment.ShapeFunctions[1], 0.5, 1e-10);
    KRATOS_EXPECT_NEAR(r_attachment.ShapeFunctionDerivatives(0, 0), -1.0, 1e-10);
    KRATOS_EXPECT_NEAR(r_attachment.ShapeFunctionDerivatives(2, 0), 1.0, 1e-10);
    KRATOS_EXPECT_NEAR(r_attachment.ReferencePosition[1], 1.0, 1e-10);
    KRATOS_EXPECT_NEAR(r_beam.GetNode(1).Y(), 10.0, 1e-10);
    KRATOS_EXPECT_NEAR(r_surface.GetNode(1).Y(), 21.0, 1e-10);
}

KRATOS_TEST_CASE_IN_SUITE(IgaBeamMapperProjectsRationalCurve, KratosMappingApplicationSerialTestSuite)
{
    Model model;
    auto& r_beam = model.CreateModelPart("beam");
    auto& r_surface = model.CreateModelPart("surface");
    r_beam.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_beam.AddNodalSolutionStepVariable(ROTATION);
    const auto p_curve = AddBeamCurve(r_beam, true);
    r_surface.AddNodalSolutionStepVariable(DISPLACEMENT);
    array_1d<double, 3> local = ZeroVector(3);
    local[0] = 0.37;
    array_1d<double, 3> position = ZeroVector(3);
    p_curve->GlobalCoordinates(position, local);
    r_surface.CreateNewNode(1, 1.2 * position[0], 1.2 * position[1], 0.3);

    MapperType mapper(r_beam, r_surface, Parameters("{}"));
    const auto& r_attachment = mapper.GetReferenceAttachments().at(0);
    KRATOS_EXPECT_NEAR(r_attachment.Parameter, local[0], 1e-8);
    KRATOS_EXPECT_VECTOR_NEAR(r_attachment.CenterlinePosition, position, 1e-8);
    KRATOS_EXPECT_NEAR(r_attachment.NormalOffset, -0.2, 1e-8);
    KRATOS_EXPECT_NEAR(r_attachment.BinormalOffset, 0.3, 1e-8);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        MapperType(r_beam, r_surface, Parameters(R"({"projection_max_iterations":1})")),
        "minimum-distance projection failed");
}

KRATOS_TEST_CASE_IN_SUITE(IgaBeamMapperRejectsTangentialEndpointOffset, KratosMappingApplicationSerialTestSuite)
{
    Model model;
    auto& r_beam = model.CreateModelPart("beam");
    auto& r_surface = model.CreateModelPart("surface");
    r_beam.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_beam.AddNodalSolutionStepVariable(ROTATION);
    const auto p_curve = AddBeamCurve(r_beam);
    r_surface.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_surface.CreateNewNode(1, 3.0, 1.0, 0.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        MapperType(r_beam, r_surface, Parameters("{}")), "has a tangential offset");
}

KRATOS_TEST_CASE_IN_SUITE(IgaBeamMapperSelectsClosestKnotSpan, KratosMappingApplicationSerialTestSuite)
{
    Model model;
    auto& r_beam = model.CreateModelPart("beam");
    auto& r_surface = model.CreateModelPart("surface");
    r_beam.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_beam.AddNodalSolutionStepVariable(ROTATION);
    Geometry<Node>::PointsArrayType points;
    const std::vector<std::array<double, 3>> coordinates = {
        {0.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {2.0, 1.0, 0.0},
        {2.0, 2.0, 0.0}, {0.0, 2.0, 0.0}};
    for (std::size_t i = 0; i < coordinates.size(); ++i) {
        points.push_back(r_beam.CreateNewNode(10 * (i + 1),
            coordinates[i][0], coordinates[i][1], coordinates[i][2]));
    }
    Vector knots(6);
    knots[0] = knots[1] = 0.0;
    knots[2] = knots[3] = 0.5;
    knots[4] = knots[5] = 1.0;
    MapperType::NurbsCurveType curve(points, 2, knots);
    Geometry<Node>::GeometriesArrayType quadrature_points;
    Geometry<Node>::IntegrationPointsArrayType integration_points;
    auto integration_info = curve.GetDefaultIntegrationInfo();
    curve.CreateIntegrationPoints(integration_points, integration_info);
    curve.CreateQuadraturePointGeometries(quadrature_points, 3, integration_points, integration_info);
    for (std::size_t i = 0; i < quadrature_points.size(); ++i) {
        r_beam.AddElement(Kratos::make_intrusive<ReferenceFrameElement>(i + 1, quadrature_points(i)));
    }
    array_1d<double, 3> local = ZeroVector(3);
    local[0] = 0.85;
    std::vector<array_1d<double, 3>> derivatives;
    curve.GlobalSpaceDerivatives(derivatives, local, 1);
    const array_1d<double, 3> tangent = derivatives[1] / norm_2(derivatives[1]);
    r_surface.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_surface.CreateNewNode(1, derivatives[0][0] + 0.1 * tangent[1],
        derivatives[0][1] - 0.1 * tangent[0], 0.2);
    MapperType mapper(r_beam, r_surface, Parameters("{}"));
    const auto& r_attachment = mapper.GetReferenceAttachments().at(0);
    KRATOS_EXPECT_NEAR(r_attachment.Parameter, 0.85, 1e-8);
    KRATOS_EXPECT_EQ(r_attachment.ControlPoints[0]->Id(), 30);
    KRATOS_EXPECT_EQ(r_attachment.ControlPoints[2]->Id(), 50);
    KRATOS_EXPECT_NEAR(r_attachment.NormalOffset, -0.1, 1e-8);
    KRATOS_EXPECT_NEAR(r_attachment.BinormalOffset, 0.2, 1e-8);
}


KRATOS_TEST_CASE_IN_SUITE(IgaBeamMapperCreatesFromFactory, KratosMappingApplicationSerialTestSuite)
{
    Model model;
    auto& r_beam = model.CreateModelPart("beam");
    auto& r_surface = model.CreateModelPart("surface");
    r_beam.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_beam.AddNodalSolutionStepVariable(ROTATION);
    const auto p_curve = AddBeamCurve(r_beam);
    AddSurfaceNode(r_surface);

    auto p_mapper = FactoryType::CreateMapper(
        r_beam, r_surface, Parameters(R"({"mapper_type":"iga_beam_mapper"})"));
    KRATOS_EXPECT_EQ(p_mapper->Info(), "IgaBeamMapper");
    KRATOS_EXPECT_EQ(&p_mapper->GetInterfaceModelPartOrigin(), &r_beam);
    KRATOS_EXPECT_EQ(&p_mapper->GetInterfaceModelPartDestination(), &r_surface);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_mapper->Map(VELOCITY, DISPLACEMENT, Flags()), "maps DISPLACEMENT to DISPLACEMENT");
    p_mapper->Map(DISPLACEMENT, DISPLACEMENT, Flags());
    KRATOS_EXPECT_VECTOR_NEAR(r_surface.GetNode(1).FastGetSolutionStepValue(DISPLACEMENT), ZeroVector(3), 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(IgaBeamMapperRequiresTwistData, KratosMappingApplicationSerialTestSuite)
{
    Model model;
    auto& r_beam = model.CreateModelPart("beam");
    auto& r_surface = model.CreateModelPart("surface");
    r_beam.AddNodalSolutionStepVariable(DISPLACEMENT);
    const auto p_curve = AddBeamCurve(r_beam);
    AddSurfaceNode(r_surface);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        MapperType(r_beam, r_surface, Parameters("{}")), "missing historical ROTATION_X");
}

KRATOS_TEST_CASE_IN_SUITE(IgaBeamMapperRejectsMultipleCurves, KratosMappingApplicationSerialTestSuite)
{
    Model model;
    auto& r_beam = model.CreateModelPart("beam");
    auto& r_surface = model.CreateModelPart("surface");
    r_beam.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_beam.AddNodalSolutionStepVariable(ROTATION);
    const auto p_first_curve = AddBeamCurve(r_beam);
    const auto p_second_curve = AddBeamCurve(r_beam);
    AddSurfaceNode(r_surface);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        MapperType(r_beam, r_surface, Parameters("{}")), "exactly one parent NURBS curve");
}

KRATOS_TEST_CASE_IN_SUITE(IgaBeamMapperRequiresSurfaceDisplacement, KratosMappingApplicationSerialTestSuite)
{
    Model model;
    auto& r_beam = model.CreateModelPart("beam");
    auto& r_surface = model.CreateModelPart("surface");
    r_beam.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_beam.AddNodalSolutionStepVariable(ROTATION);
    const auto p_curve = AddBeamCurve(r_beam);
    r_surface.CreateNewNode(1, 1.0, 1.0, 0.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        MapperType(r_beam, r_surface, Parameters("{}")), "surface node 1 is missing historical DISPLACEMENT");
}

} // namespace Kratos::Testing
