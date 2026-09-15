//    Kratos Multiphysics
//  License: BSD License (kratos/license.txt)

#include "testing/testing.h"
#include "containers/model.h"
#include "custom_elements/isogeometric_beam_element.h"
#include "custom_constitutive/bernoulli_beam_elastic_constitutive_law.h"
#include "geometries/nurbs_curve_geometry.h"

namespace Kratos::Testing
{
namespace
{
struct BeamMassFixture
{
    Model CurrentModel;
    ModelPart& Beam;
    Geometry<Node>::GeometriesArrayType Quadrature;
    Properties::Pointer pProperties;
    NurbsCurveGeometry<3, PointerVector<Node>>::Pointer Curve;

    BeamMassFixture(const bool Curved = false, const bool DynamicData = true) : Beam(CurrentModel.CreateModelPart("beam"))
    {
        Beam.SetBufferSize(2);
        if (DynamicData) {
            Beam.AddNodalSolutionStepVariable(VELOCITY);
            Beam.AddNodalSolutionStepVariable(ACCELERATION);
            Beam.AddNodalSolutionStepVariable(ANGULAR_VELOCITY);
            Beam.AddNodalSolutionStepVariable(ANGULAR_ACCELERATION);
        }
        Beam.AddNodalSolutionStepVariable(DISPLACEMENT);
        Beam.AddNodalSolutionStepVariable(ROTATION);
        Geometry<Node>::PointsArrayType points;
        for (std::size_t i = 0; i < 3; ++i) {
            points.push_back(Beam.CreateNewNode(i + 1, static_cast<double>(i), 0.0, 0.0));
        }
        Vector knots(4);
        knots[0] = knots[1] = 0.0;
        knots[2] = knots[3] = 1.0;
        if (Curved) {
            points[0].X0() = points[0].X() = 2.0;
            points[1].X0() = points[1].X() = 2.0;
            points[1].Y0() = points[1].Y() = 2.0;
            points[2].X0() = points[2].X() = 0.0;
            points[2].Y0() = points[2].Y() = 2.0;
            Vector weights(3);
            weights[0] = weights[2] = 1.0;
            weights[1] = std::sqrt(0.5);
            Curve = Kratos::make_shared<NurbsCurveGeometry<3, PointerVector<Node>>>(points, 2, knots, weights);
        } else {
            Curve = Kratos::make_shared<NurbsCurveGeometry<3, PointerVector<Node>>>(points, 2, knots);
        }
        auto integration_info = Curve->GetDefaultIntegrationInfo();
        Geometry<Node>::IntegrationPointsArrayType integration_points;
        Curve->CreateIntegrationPoints(integration_points, integration_info);
        Curve->CreateQuadraturePointGeometries(Quadrature, 3, integration_points, integration_info);
        pProperties = Beam.CreateNewProperties(1);
        pProperties->SetValue(DENSITY, 3.0);
        pProperties->SetValue(CROSS_AREA, 2.0);
        pProperties->SetValue(I_N, 0.4);
        pProperties->SetValue(I_V, 0.7);
        pProperties->SetValue(I_T, 99.0); // Deliberately different from the polar area moment.
        Vector tangent = ZeroVector(3);
        tangent[Curved ? 1 : 0] = 1.0;
        Vector normal = ZeroVector(3);
        normal[2] = 1.0;
        pProperties->SetValue(T_0, tangent);
        pProperties->SetValue(N_0, normal);
        Matrix orientation = ZeroMatrix(2, 4);
        orientation(1, 0) = 1.0;
        orientation(0, 3) = orientation(1, 3) = 1.0;
        pProperties->SetValue(LOCAL_AXIS_ORIENTATION, orientation);
    }

    Matrix Mass(const std::string& rFormulation)
    {
        pProperties->SetValue(BEAM_MASS_FORMULATION, rFormulation);
        Matrix total = ZeroMatrix(12, 12);
        for (std::size_t p = 0; p < Quadrature.size(); ++p) {
            IsogeometricBeamElement element(p + 1, Quadrature(p), pProperties);
            Matrix mass;
            element.CalculateMassMatrix(mass, Beam.GetProcessInfo());
            total += mass;
        }
        return total;
    }
};
}

KRATOS_TEST_CASE_IN_SUITE(IsogeometricBeamSimplifiedMass, KratosIgaFastSuite)
{
    BeamMassFixture fixture;
    const Matrix mass = fixture.Mass("simplified");
    // Exact quadratic Bernstein integral on a beam of length two.
    const double coefficients[3][3] = {{6, 3, 1}, {3, 4, 3}, {1, 3, 6}};
    for (std::size_t i = 0; i < 12; ++i) {
        for (std::size_t j = 0; j < 12; ++j) {
            const double expected = i % 4 == j % 4
                ? 6.0 * coefficients[i/4][j/4] / 30.0 * (i % 4 == 3 ? 1.1 : 2.0) : 0.0;
            KRATOS_EXPECT_NEAR(mass(i, j), expected, 1e-12);
        }
    }
    for (auto& r_node : fixture.Beam.Nodes()) {
        r_node.FastGetSolutionStepValue(DISPLACEMENT_Y) = r_node.Id();
        r_node.FastGetSolutionStepValue(ROTATION_X) = 1.5 * r_node.Id();
        r_node.X() += 100.0;
    }
    KRATOS_EXPECT_MATRIX_NEAR(fixture.Mass("simplified"), mass, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(IsogeometricBeamFullMassReferenceEnergy, KratosIgaFastSuite)
{
    BeamMassFixture fixture;
    const Matrix mass = fixture.Mass("full");
    const Matrix simplified = fixture.Mass("simplified");
    const Matrix transpose = trans(mass);
    KRATOS_EXPECT_MATRIX_NEAR(mass, transpose, 1e-12);
    // Translation along each global direction must give rho*A*L.
    for (std::size_t k = 0; k < 3; ++k) {
        Vector velocity = ZeroVector(12);
        for (std::size_t i = 0; i < 3; ++i) velocity[4*i+k] = 1.0;
        KRATOS_EXPECT_NEAR(inner_prod(velocity, prod(mass, velocity)), 12.0, 1e-12);
    }
    Vector twist = ZeroVector(12);
    Vector bending = ZeroVector(12);
    for (std::size_t i = 0; i < 3; ++i) {
        twist[4*i+3] = 1.0;
        bending[4*i+1] = static_cast<double>(i); // unit rotation rate about Z
    }
    KRATOS_EXPECT_NEAR(inner_prod(twist, prod(mass, twist)), 6.6, 1e-12);
    // Integral rho*A*x^2 dx = 16; section rotation adds rho*I_N*L = 2.4.
    KRATOS_EXPECT_NEAR(inner_prod(bending, prod(simplified, bending)), 16.0, 1e-12);
    KRATOS_EXPECT_NEAR(inner_prod(bending, prod(mass, bending)), 18.4, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(IsogeometricBeamFullMassFiniteTwist, KratosIgaFastSuite)
{
    BeamMassFixture fixture;
    const Matrix reference = fixture.Mass("full");
    for (auto& r_node : fixture.Beam.Nodes()) r_node.FastGetSolutionStepValue(ROTATION_X) = 0.5 * std::acos(-1.0);
    const Matrix twisted = fixture.Mass("full");
    Vector bending = ZeroVector(12);
    for (std::size_t i = 0; i < 3; ++i) bending[4*i+1] = static_cast<double>(i);
    // A quarter-turn swaps the principal rotary inertias.
    KRATOS_EXPECT_NEAR(inner_prod(bending, prod(twisted, bending)), 16.0 + 6.0 * 0.7, 1e-12);
    KRATOS_EXPECT_GT(norm_frobenius(twisted - reference), 0.1);
    for (auto& r_node : fixture.Beam.Nodes()) r_node.Y() += 100.0;
    KRATOS_EXPECT_MATRIX_NEAR(fixture.Mass("full"), twisted, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(IsogeometricBeamFullMassCurvedKineticEnergy, KratosIgaFastSuite)
{
    BeamMassFixture fixture(true);
    Vector state = ZeroVector(12);
    for (std::size_t i = 0; i < 3; ++i) {
        auto& r_node = fixture.Beam.GetNode(i + 1);
        for (std::size_t k = 0; k < 3; ++k) {
            state[4*i+k] = 0.1 * (i + 1) * (k + 1);
            r_node.FastGetSolutionStepValue(DISPLACEMENT)[k] = state[4*i+k];
        }
        state[4*i+3] = 0.8 * (i + 1);
        r_node.FastGetSolutionStepValue(ROTATION_X) = state[4*i+3];
    }
    Matrix expected = ZeroMatrix(12, 12);
    const double epsilon = 1e-6;
    for (std::size_t p = 0; p < fixture.Quadrature.size(); ++p) {
        IsogeometricBeamElement element(p+1, fixture.Quadrature(p), fixture.pProperties);
        auto& r_geometry = element.GetGeometry();
        std::vector<Matrix> frames;
        element.CalculateOnIntegrationPoints(LOCAL_AXES_MATRIX, frames, fixture.Beam.GetProcessInfo());
        const auto& r_shape = r_geometry.ShapeFunctionsValues();
        const auto& r_derivatives = r_geometry.ShapeFunctionDerivatives(1, 0);
        // Independent minimal-rotation identity for a perpendicular director:
        // b = B - (B.t)/(1+T.t)*(T+t). Differentiate numerically.
        const auto evaluate = [&](const Vector& rState) {
            Matrix result = ZeroMatrix(3, 3); // centerline, normal, binormal
            array_1d<double, 3> derivative = ZeroVector(3);
            double phi = 0.0;
            for (std::size_t i = 0; i < 3; ++i) {
                for (std::size_t k = 0; k < 3; ++k) {
                    const double x = r_geometry[i].GetInitialPosition()[k] + rState[4*i+k];
                    result(0, k) += r_shape(0, i) * x;
                    derivative[k] += r_derivatives(i, 0) * x;
                }
                phi += r_shape(0, i) * rState[4*i+3];
            }
            const array_1d<double, 3> t = derivative / norm_2(derivative);
            const array_1d<double, 3> T = row(frames[0], 0);
            for (std::size_t d = 1; d < 3; ++d) {
                const array_1d<double, 3> B = row(frames[0], d);
                const array_1d<double, 3> b = B - (inner_prod(B, t) / (1.0 + inner_prod(T, t))) * (T + t);
                const array_1d<double, 3> director = std::cos(phi) * b
                    + std::sin(phi) * MathUtils<double>::CrossProduct(t, b);
                row(result, d) = director;
            }
            return result;
        };
        Matrix derivatives = ZeroMatrix(9, 12);
        for (std::size_t j = 0; j < 12; ++j) {
            Vector plus = state;
            Vector minus = state;
            plus[j] += epsilon;
            minus[j] -= epsilon;
            const Matrix diff = (evaluate(plus) - evaluate(minus)) / (2.0 * epsilon);
            for (std::size_t d = 0; d < 3; ++d) {
                for (std::size_t k = 0; k < 3; ++k) derivatives(3*d+k, j) = diff(d, k);
            }
        }
        array_1d<double, 3> reference_derivative = ZeroVector(3);
        for (std::size_t i = 0; i < 3; ++i) reference_derivative += r_derivatives(i, 0) * r_geometry[i].GetInitialPosition();
        const double weight = 3.0 * norm_2(reference_derivative) * r_geometry.IntegrationPoints()[0].Weight();
        const double section_weights[3] = {2.0, 0.7, 0.4};
        for (std::size_t i = 0; i < 12; ++i) {
            for (std::size_t j = 0; j < 12; ++j) {
                for (std::size_t k = 0; k < 9; ++k) {
                    expected(i, j) += weight * section_weights[k/3] * derivatives(k, i) * derivatives(k, j);
                }
            }
        }
    }
    const Matrix mass = fixture.Mass("full");
    KRATOS_EXPECT_MATRIX_NEAR(mass, expected, 1e-8);
    const Matrix transpose = trans(mass);
    KRATOS_EXPECT_MATRIX_NEAR(mass, transpose, 1e-12);
    KRATOS_EXPECT_GT(inner_prod(state, prod(mass, state)), 0.0);
}

KRATOS_TEST_CASE_IN_SUITE(IsogeometricBeamDerivativeVectorsAndDamping, KratosIgaFastSuite)
{
    BeamMassFixture fixture;
    IsogeometricBeamElement element(1, fixture.Quadrature(0), fixture.pProperties);
    for (auto& r_node : fixture.Beam.Nodes()) {
        for (int step = 0; step < 2; ++step) {
            for (std::size_t k = 0; k < 3; ++k) {
                const double value = 100.0 * step + 10.0 * r_node.Id() + k;
                r_node.FastGetSolutionStepValue(VELOCITY, step)[k] = value;
                r_node.FastGetSolutionStepValue(ACCELERATION, step)[k] = -value;
                r_node.FastGetSolutionStepValue(ANGULAR_VELOCITY, step)[k] = value + 3.0;
                r_node.FastGetSolutionStepValue(ANGULAR_ACCELERATION, step)[k] = -value - 3.0;
            }
        }
    }
    Vector first(1, 999.0);
    Vector second(20, 999.0);
    for (int step = 0; step < 2; ++step) {
        element.GetFirstDerivativesVector(first, step);
        element.GetSecondDerivativesVector(second, step);
        KRATOS_EXPECT_EQ(first.size(), 12);
        KRATOS_EXPECT_EQ(second.size(), 12);
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t k = 0; k < 4; ++k) {
                const double expected = 100.0 * step + 10.0 * (i + 1) + k;
                KRATOS_EXPECT_NEAR(first[4*i+k], expected, 1e-12);
                KRATOS_EXPECT_NEAR(second[4*i+k], -expected, 1e-12);
            }
        }
    }
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(element.GetFirstDerivativesVector(first, -1), "Invalid solution step");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(element.GetSecondDerivativesVector(second, 2), "Invalid solution step");
    BeamMassFixture missing_data(false, false);
    IsogeometricBeamElement missing_element(1, missing_data.Quadrature(0), missing_data.pProperties);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(missing_element.GetFirstDerivativesVector(first), "requires historical VELOCITY");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(missing_element.GetSecondDerivativesVector(second), "requires historical ACCELERATION");
    // Damping does not require a mass formulation or material initialization.
    Matrix damping(2, 7, 123.0);
    element.CalculateDampingMatrix(damping, fixture.Beam.GetProcessInfo());
    const Matrix zero = ZeroMatrix(12, 12);
    KRATOS_EXPECT_MATRIX_NEAR(damping, zero, 1e-12);
    damping(0, 0) = 99.0;
    element.CalculateDampingMatrix(damping, fixture.Beam.GetProcessInfo());
    KRATOS_EXPECT_MATRIX_NEAR(damping, zero, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(IsogeometricBeamRayleighDamping, KratosIgaFastSuite)
{
    BeamMassFixture fixture;
    auto& r_properties = *fixture.pProperties;
    auto& r_info = fixture.Beam.GetProcessInfo();
    r_properties.SetValue(YOUNG_MODULUS, 1000.0);
    r_properties.SetValue(POISSON_RATIO, 0.25);
    r_properties.SetValue(CONSTITUTIVE_LAW, Kratos::make_shared<BernoulliBeamElasticConstitutiveLaw>());
    IsogeometricBeamElement element(1, fixture.Quadrature(0), fixture.pProperties);
    element.Initialize(r_info);
    Matrix stiffness;
    element.CalculateLeftHandSide(stiffness, r_info);
    KRATOS_EXPECT_GT(norm_frobenius(stiffness), 0.0);
    // Stiffness-only damping must not require the mass formulation property.
    r_info.SetValue(RAYLEIGH_BETA, 0.02);
    Matrix damping;
    element.CalculateDampingMatrix(damping, r_info);
    Matrix expected = 0.02 * stiffness;
    KRATOS_EXPECT_MATRIX_NEAR(damping, expected, 1e-10);
    for (const std::string formulation : {"simplified", "full"}) {
        r_properties.SetValue(BEAM_MASS_FORMULATION, formulation);
        Matrix mass;
        element.CalculateMassMatrix(mass, r_info);
        r_info.SetValue(RAYLEIGH_ALPHA, 0.3);
        element.CalculateDampingMatrix(damping, r_info);
        expected = 0.3 * mass + 0.02 * stiffness;
        KRATOS_EXPECT_MATRIX_NEAR(damping, expected, 1e-10);
        // Explicit property values override ProcessInfo, including zero.
        r_properties.SetValue(RAYLEIGH_ALPHA, 0.5);
        r_properties.SetValue(RAYLEIGH_BETA, 0.0);
        element.CalculateDampingMatrix(damping, r_info);
        expected = 0.5 * mass;
        KRATOS_EXPECT_MATRIX_NEAR(damping, expected, 1e-10);
        r_properties.Erase(RAYLEIGH_ALPHA);
        r_properties.Erase(RAYLEIGH_BETA);
    }
    r_properties.SetValue(RAYLEIGH_ALPHA, -1.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(element.CalculateDampingMatrix(damping, r_info), "must be finite and non-negative");
}

KRATOS_TEST_CASE_IN_SUITE(IsogeometricBeamMassValidation, KratosIgaFastSuite)
{
    BeamMassFixture fixture;
    IsogeometricBeamElement element(1, fixture.Quadrature(0), fixture.pProperties);
    Matrix mass;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(element.CalculateMassMatrix(mass, fixture.Beam.GetProcessInfo()), "Specify BEAM_MASS_FORMULATION");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(fixture.Mass("invalid"), "Unknown BEAM_MASS_FORMULATION");
    fixture.pProperties->SetValue(DENSITY, -1.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(fixture.Mass("simplified"), "must be finite and positive");
    fixture.pProperties->SetValue(DENSITY, 3.0);
    for (auto& r_node : fixture.Beam.Nodes()) r_node.FastGetSolutionStepValue(DISPLACEMENT_X) = -2.0 * r_node.X0();
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(fixture.Mass("full"), "opposite reference and current tangents");
}

} // namespace Kratos::Testing
