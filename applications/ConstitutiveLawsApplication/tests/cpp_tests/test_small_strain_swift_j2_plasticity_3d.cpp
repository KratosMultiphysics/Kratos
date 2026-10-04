// License: BSD License
// Kratos default license: kratos/license.txt

#include <cmath>
#include <limits>

#include "tests/cpp_tests/constitutive_laws_fast_suite.h"
#include "containers/model.h"
#include "geometries/tetrahedra_3d_4.h"
#include "includes/stream_serializer.h"
#include "utilities/read_materials_utility.h"
#include "constitutive_laws_application_variables.h"
#include "custom_constitutive/small_strains/plasticity/small_strain_swift_j2_plasticity_3d.h"

namespace Kratos::Testing
{
namespace
{

struct SwiftMaterialPoint
{
    Model CurrentModel;
    ModelPart& rModelPart;
    Properties Material{1};
    Geometry<Node>::Pointer pGeometry;
    Vector ShapeFunctions = ScalarVector(4, 0.25);
    Vector Strain = ZeroVector(6);
    Vector Stress = ZeroVector(6);
    Matrix Tangent = ZeroMatrix(6, 6);
    Matrix F = IdentityMatrix(3);
    ConstitutiveLaw::Parameters Values;
    SmallStrainSwiftJ2Plasticity3D Law;

    SwiftMaterialPoint() : rModelPart(CurrentModel.CreateModelPart("Main"))
    {
        pGeometry = Kratos::make_shared<Tetrahedra3D4<Node>>(
            rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0),
            rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0),
            rModelPart.CreateNewNode(3, 0.0, 1.0, 0.0),
            rModelPart.CreateNewNode(4, 0.0, 0.0, 1.0));
        Material.SetValue(YOUNG_MODULUS, 210000.0); // MPa
        Material.SetValue(POISSON_RATIO, 0.3);
        Material.SetValue(SWIFT_COEFFICIENT, 1000.0); // MPa
        Material.SetValue(SWIFT_INITIAL_STRAIN, 0.01);
        Material.SetValue(SWIFT_HARDENING_EXPONENT, 0.2);
        Values.SetElementGeometry(*pGeometry);
        Values.SetMaterialProperties(Material);
        Values.SetProcessInfo(rModelPart.GetProcessInfo());
        Values.SetShapeFunctionsValues(ShapeFunctions);
        Values.SetStrainVector(Strain);
        Values.SetStressVector(Stress);
        Values.SetConstitutiveMatrix(Tangent);
        Values.SetDeformationGradientF(F);
        Values.SetDeterminantF(1.0);
        Values.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
        Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
        Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
        Law.InitializeMaterial(Material, *pGeometry, ShapeFunctions);
    }

    double P()
    {
        double p = 0.0;
        return Law.GetValue(ACCUMULATED_PLASTIC_STRAIN, p);
    }

    Vector PlasticStrain()
    {
        Vector strain;
        return Law.GetValue(PLASTIC_STRAIN_VECTOR, strain);
    }

    void Evaluate() { Law.CalculateMaterialResponseCauchy(Values); }
    void Commit() { Law.FinalizeMaterialResponseCauchy(Values); }

    // Each iteration is a trial from the same committed state. Only the caller commits.
    void SolveUniaxialStress(const double AxialStrain)
    {
        Strain[0] = AxialStrain;
        for (unsigned int iteration = 0; iteration < 30; ++iteration) {
            Evaluate();
            if (std::hypot(Stress[1], Stress[2]) < 1.0e-8) {
                return;
            }
            const double det = Tangent(1, 1) * Tangent(2, 2) - Tangent(1, 2) * Tangent(2, 1);
            const double dy = (Tangent(2, 2) * Stress[1] - Tangent(1, 2) * Stress[2]) / det;
            const double dz = (Tangent(1, 1) * Stress[2] - Tangent(2, 1) * Stress[1]) / det;
            Strain[1] -= dy;
            Strain[2] -= dz;
        }
        KRATOS_ERROR << "Uniaxial-stress test could not equilibrate transverse stresses";
    }
};

double VonMises(const Vector& rStress)
{
    return std::sqrt(0.5 * (std::pow(rStress[0] - rStress[1], 2) +
        std::pow(rStress[1] - rStress[2], 2) + std::pow(rStress[2] - rStress[0], 2)) +
        3.0 * (rStress[3] * rStress[3] + rStress[4] * rStress[4] + rStress[5] * rStress[5]));
}

Matrix ElasticReference(const Properties& rMaterial)
{
    const double e = rMaterial[YOUNG_MODULUS];
    const double nu = rMaterial[POISSON_RATIO];
    const double mu = e / (2.0 * (1.0 + nu));
    const double lambda = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    Matrix matrix = ZeroMatrix(6, 6);
    for (unsigned int i = 0; i < 3; ++i) {
        for (unsigned int j = 0; j < 3; ++j) {
            matrix(i, j) = lambda + (i == j ? 2.0 * mu : 0.0);
        }
        matrix(i + 3, i + 3) = mu;
    }
    return matrix;
}

// Independent uniaxial solution: eps_x = sigma/E + p. Bisection does not
// use the constitutive implementation, its multiplier or its Newton derivative.
double UniaxialReferenceP(const Properties& rMaterial, const double AxialStrain)
{
    const double e = rMaterial[YOUNG_MODULUS];
    const double k = rMaterial[SWIFT_COEFFICIENT];
    const double eps0 = rMaterial[SWIFT_INITIAL_STRAIN];
    const double n = rMaterial[SWIFT_HARDENING_EXPONENT];
    if (e * AxialStrain <= k * std::pow(eps0, n)) {
        return 0.0;
    }
    double lower = 0.0;
    double upper = AxialStrain;
    for (unsigned int i = 0; i < 100; ++i) {
        const double p = 0.5 * (lower + upper);
        if (e * (AxialStrain - p) > k * std::pow(eps0 + p, n)) {
            lower = p;
        } else {
            upper = p;
        }
    }
    return 0.5 * (lower + upper);
}

void CheckTangent(SwiftMaterialPoint& rPoint)
{
    rPoint.Evaluate();
    const Matrix analytical = rPoint.Tangent;
    const Vector strain = rPoint.Strain;
    const double committed_p = rPoint.P();
    const Vector committed_plastic_strain = rPoint.PlasticStrain();
    Matrix numerical = ZeroMatrix(6, 6);
    // h is below the distance to the yield switch, including the just-plastic
    // case. A 3e-7*E absolute tolerance allows residual/h and roundoff error
    // without hiding relative errors in the dominant stiffness terms.
    constexpr double h = 1.0e-8;
    rPoint.Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
    for (unsigned int j = 0; j < 6; ++j) {
        rPoint.Strain = strain;
        rPoint.Strain[j] += h;
        rPoint.Evaluate();
        const Vector plus = rPoint.Stress;
        rPoint.Strain[j] -= 2.0 * h;
        rPoint.Evaluate();
        for (unsigned int i = 0; i < 6; ++i) {
            numerical(i, j) = (plus[i] - rPoint.Stress[i]) / (2.0 * h);
            KRATOS_EXPECT_NEAR(analytical(i, j), analytical(j, i), 1.0e-10);
        }
    }
    KRATOS_EXPECT_MATRIX_NEAR(analytical, numerical, 3.0e-7 * rPoint.Material[YOUNG_MODULUS]);
    KRATOS_EXPECT_EQ(rPoint.P(), committed_p);
    KRATOS_EXPECT_VECTOR_NEAR(rPoint.PlasticStrain(), committed_plastic_strain, 0.0);
    rPoint.Strain = strain;
    rPoint.Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
}

void CheckContinuation(SwiftMaterialPoint& rPoint, ConstitutiveLaw& rOther)
{
    double other_p = 0.0;
    Vector other_plastic;
    KRATOS_EXPECT_NEAR(rOther.GetValue(ACCUMULATED_PLASTIC_STRAIN, other_p), rPoint.P(), 1.0e-15);
    KRATOS_EXPECT_VECTOR_NEAR(rOther.GetValue(PLASTIC_STRAIN_VECTOR, other_plastic), rPoint.PlasticStrain(), 1.0e-15);
    rPoint.Strain[0] += 0.005;
    rPoint.Strain[4] += 0.007;
    rPoint.Evaluate();
    const Vector expected_stress = rPoint.Stress;
    const Matrix expected_tangent = rPoint.Tangent;
    rOther.CalculateMaterialResponseCauchy(rPoint.Values);
    KRATOS_EXPECT_VECTOR_NEAR(rPoint.Stress, expected_stress, 1.0e-10);
    KRATOS_EXPECT_MATRIX_NEAR(rPoint.Tangent, expected_tangent, 1.0e-10);
    rPoint.Commit();
    rOther.FinalizeMaterialResponseCauchy(rPoint.Values);
    KRATOS_EXPECT_NEAR(rOther.GetValue(ACCUMULATED_PLASTIC_STRAIN, other_p), rPoint.P(), 1.0e-15);
    KRATOS_EXPECT_VECTOR_NEAR(rOther.GetValue(PLASTIC_STRAIN_VECTOR, other_plastic), rPoint.PlasticStrain(), 1.0e-15);
}

} // namespace

KRATOS_TEST_CASE_IN_SUITE(SwiftJ2FlowStressEquation, KratosConstitutiveLawsFastSuite)
{
    SwiftMaterialPoint point;
    for (const double n : {0.2, 1.0, 2.0}) {
        point.Material[SWIFT_HARDENING_EXPONENT] = n;
        for (const double p : {0.0, 0.001, 0.02, 0.1, 0.5}) {
            point.Law.SetValue(ACCUMULATED_PLASTIC_STRAIN, p, point.rModelPart.GetProcessInfo());
            double flow_stress = 0.0;
            point.Law.CalculateValue(point.Values, YIELD_STRESS, flow_stress);
            KRATOS_EXPECT_NEAR(flow_stress, 1000.0 * std::pow(0.01 + p, n), 1.0e-11);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(SwiftJ2InitialYieldStress, KratosConstitutiveLawsFastSuite)
{
    SwiftMaterialPoint point;
    double flow_stress = 0.0;
    point.Law.CalculateValue(point.Values, YIELD_STRESS, flow_stress);
    // Independent numerical value of 1000*(0.01)^0.2 in MPa.
    KRATOS_EXPECT_NEAR(flow_stress, 398.1071705534972, 1.0e-10);
    KRATOS_EXPECT_EQ(point.P(), 0.0);
}

KRATOS_TEST_CASE_IN_SUITE(SwiftJ2ElasticResponse, KratosConstitutiveLawsFastSuite)
{
    SwiftMaterialPoint point;
    for (unsigned int i = 0; i < 6; ++i) {
        point.Strain[i] = (static_cast<double>(i) - 2.0) * 1.0e-6;
    }
    point.Evaluate();
    const Matrix elastic = ElasticReference(point.Material);
    const Vector expected = prod(elastic, point.Strain);
    KRATOS_EXPECT_VECTOR_NEAR(point.Stress, expected, 1.0e-12);
    KRATOS_EXPECT_MATRIX_NEAR(point.Tangent, elastic, 1.0e-10);
    point.Commit();
    KRATOS_EXPECT_EQ(point.P(), 0.0);
}

KRATOS_TEST_CASE_IN_SUITE(SwiftJ2ElasticPlasticTransitionAndPureShear, KratosConstitutiveLawsFastSuite)
{
    SwiftMaterialPoint point;
    const double mu = 210000.0 / 2.6;
    const double yield = 1000.0 * std::pow(0.01, 0.2);
    for (const double factor : {0.5, 0.999, 1.001, 1.1, 2.0, 5.0}) {
        point.Strain[3] = factor * yield / (std::sqrt(3.0) * mu);
        const double previous_p = point.P();
        point.Evaluate();
        KRATOS_EXPECT_EQ(point.P(), previous_p);
        point.Commit();
        const double q = VonMises(point.Stress);
        const double flow = 1000.0 * std::pow(0.01 + point.P(), 0.2);
        if (factor < 1.0) {
            KRATOS_EXPECT_EQ(point.P(), 0.0);
            KRATOS_EXPECT_NEAR(q, factor * yield, 1.0e-9);
        } else {
            KRATOS_EXPECT_GT(point.P(), previous_p);
            KRATOS_EXPECT_NEAR(q, flow, 1.0e-8);
            KRATOS_EXPECT_NEAR(point.P(), point.PlasticStrain()[3] / std::sqrt(3.0), 1.0e-14);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(SwiftJ2PureShearIndependentAnalyticalReference, KratosConstitutiveLawsFastSuite)
{
    SwiftMaterialPoint point;
    const double young = point.Material[YOUNG_MODULUS];
    const double nu = point.Material[POISSON_RATIO];
    const double k = point.Material[SWIFT_COEFFICIENT];
    const double eps0 = point.Material[SWIFT_INITIAL_STRAIN];
    const double n = point.Material[SWIFT_HARDENING_EXPONENT];
    const double shear_modulus = young / (2.0 * (1.0 + nu));
    const double sqrt_three = std::sqrt(3.0);
    const double gamma = 0.02; // Well beyond initial yield at gamma = 0.002845728...

    // Pure shear has trace(sigma) = 0 and s:s = 2*tau^2, hence
    // sigma_eq = sqrt(3/2 * 2*tau^2) = sqrt(3)*abs(tau).
    // Associated flow has only xy and yx plastic strain components. Since
    // gamma_p = 2*epsilon_p_xy, epsilon_p_dot:epsilon_p_dot = gamma_p_dot^2/2.
    // Thus p_dot = abs(gamma_p_dot)/sqrt(3). Starting from virgin material,
    // monotonic pure shear gives p = abs(gamma_p)/sqrt(3); for positive shear,
    // gamma_p = sqrt(3)*p. Combining tau = G*(gamma-gamma_p) with Swift yield
    // gives F(p) = sqrt(3)*G*(gamma-sqrt(3)*p) - K*(eps0+p)^n = 0.
    // This reference uses bisection in p, without the law's return mapping.
    const auto reference_residual = [&](const double P) {
        return sqrt_three * shear_modulus * (gamma - sqrt_three * P)
            - k * std::pow(eps0 + P, n);
    };
    double lower = 0.0;
    double upper = gamma / sqrt_three;
    // F is strictly decreasing; these endpoints bracket its positive root.
    KRATOS_EXPECT_GT(reference_residual(lower), 0.0);
    KRATOS_EXPECT_LT(reference_residual(upper), 0.0);
    for (unsigned int iteration = 0; iteration < 100; ++iteration) {
        const double p = lower + 0.5 * (upper - lower);
        if (reference_residual(p) > 0.0) {
            lower = p;
        } else {
            upper = p;
        }
    }
    const double p_ref = lower + 0.5 * (upper - lower);
    const double gamma_p_ref = sqrt_three * p_ref;
    const double tau_ref = shear_modulus * (gamma - gamma_p_ref);
    const double sigma_eq_ref = sqrt_three * std::abs(tau_ref);
    const double sigma_y_ref = k * std::pow(eps0 + p_ref, n);
    // Match the existing pure-shear stress/yield and uniaxial strain tolerances.
    constexpr double stress_tolerance = 1.0e-8; // MPa
    constexpr double strain_tolerance = 1.0e-12;
    KRATOS_EXPECT_GT(p_ref, 0.0);
    KRATOS_EXPECT_NEAR(sigma_eq_ref, sigma_y_ref, stress_tolerance);

    const double initial_p = point.P();
    const Vector initial_plastic_strain = point.PlasticStrain();
    point.Strain[3] = gamma; // Only engineering xy shear is prescribed.
    point.Evaluate();
    KRATOS_EXPECT_NEAR(point.Stress[3], tau_ref, stress_tolerance);
    KRATOS_EXPECT_EQ(point.P(), initial_p);
    KRATOS_EXPECT_VECTOR_NEAR(point.PlasticStrain(), initial_plastic_strain, 0.0);

    point.Commit();
    const double p = point.P();
    const Vector plastic_strain = point.PlasticStrain();
    KRATOS_EXPECT_NEAR(point.Stress[3], tau_ref, stress_tolerance);
    KRATOS_EXPECT_NEAR(p, p_ref, strain_tolerance);
    KRATOS_EXPECT_NEAR(plastic_strain[3], gamma_p_ref, strain_tolerance);
    KRATOS_EXPECT_NEAR(plastic_strain[3], sqrt_three * p, strain_tolerance);
    KRATOS_EXPECT_NEAR(p, std::abs(plastic_strain[3]) / sqrt_three, strain_tolerance);
    KRATOS_EXPECT_NEAR(sqrt_three * std::abs(point.Stress[3]), sigma_eq_ref, stress_tolerance);
    KRATOS_EXPECT_NEAR(sqrt_three * std::abs(point.Stress[3]),
        k * std::pow(eps0 + p, n), stress_tolerance);
    for (unsigned int i = 0; i < 6; ++i) {
        if (i != 3) {
            KRATOS_EXPECT_NEAR(point.Stress[i], 0.0, stress_tolerance);
            KRATOS_EXPECT_NEAR(plastic_strain[i], 0.0, strain_tolerance);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(SwiftJ2MonotonicUniaxialStress, KratosConstitutiveLawsFastSuite)
{
    for (const double n : {0.2, 1.0, 2.0}) {
        SwiftMaterialPoint point;
        point.Material[SWIFT_HARDENING_EXPONENT] = n;
        for (const double axial : {0.0001, 0.001, 0.003, 0.006, 0.01, 0.02}) {
            point.SolveUniaxialStress(axial);
            point.Commit();
            const double p = UniaxialReferenceP(point.Material, axial);
            KRATOS_EXPECT_NEAR(point.Stress[0], 210000.0 * (axial - p), 1.0e-7);
            KRATOS_EXPECT_NEAR(point.Stress[1], 0.0, 1.0e-8);
            KRATOS_EXPECT_NEAR(point.Stress[2], 0.0, 1.0e-8);
            KRATOS_EXPECT_NEAR(point.P(), p, 1.0e-12);
            const Vector plastic = point.PlasticStrain();
            KRATOS_EXPECT_NEAR(plastic[0], p, 1.0e-12);
            KRATOS_EXPECT_NEAR(plastic[1], -0.5 * p, 1.0e-12);
            KRATOS_EXPECT_NEAR(plastic[2], -0.5 * p, 1.0e-12);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(SwiftJ2Unloading, KratosConstitutiveLawsFastSuite)
{
    SwiftMaterialPoint point;
    point.SolveUniaxialStress(0.02);
    point.Commit();
    const double stress = point.Stress[0];
    const double p = point.P();
    const Vector plastic = point.PlasticStrain();
    const double decrement = 0.5 * stress / 210000.0;
    point.SolveUniaxialStress(0.02 - decrement);
    point.Commit();
    KRATOS_EXPECT_NEAR(stress - point.Stress[0], 210000.0 * decrement, 1.0e-8);
    KRATOS_EXPECT_EQ(point.P(), p);
    KRATOS_EXPECT_VECTOR_NEAR(point.PlasticStrain(), plastic, 0.0);
    point.SolveUniaxialStress(p);
    point.Commit();
    KRATOS_EXPECT_NEAR(point.Stress[0], 0.0, 1.0e-8);
    KRATOS_EXPECT_EQ(point.P(), p);
}

KRATOS_TEST_CASE_IN_SUITE(SwiftJ2UniaxialCompressionSymmetry, KratosConstitutiveLawsFastSuite)
{
    SwiftMaterialPoint tension;
    SwiftMaterialPoint compression;
    for (const double axial : {0.001, 0.003, 0.01, 0.02}) {
        tension.SolveUniaxialStress(axial);
        compression.SolveUniaxialStress(-axial);
        tension.Commit();
        compression.Commit();
        const double p = UniaxialReferenceP(compression.Material, axial);
        KRATOS_EXPECT_NEAR(compression.Stress[0], -210000.0 * (axial - p), 1.0e-7);
        KRATOS_EXPECT_NEAR(compression.Stress[1], 0.0, 1.0e-8);
        KRATOS_EXPECT_NEAR(compression.Stress[2], 0.0, 1.0e-8);
        KRATOS_EXPECT_NEAR(compression.P(), p, 1.0e-12);
        const Vector opposite = -tension.PlasticStrain();
        KRATOS_EXPECT_VECTOR_NEAR(compression.PlasticStrain(), opposite, 1.0e-14);
        KRATOS_EXPECT_MATRIX_NEAR(compression.Tangent, tension.Tangent, 1.0e-8);
    }
}

KRATOS_TEST_CASE_IN_SUITE(SwiftJ2StronglyNonlinearHardening, KratosConstitutiveLawsFastSuite)
{
    SwiftMaterialPoint point;
    point.Material[SWIFT_INITIAL_STRAIN] = 1.0e-12;
    point.Material[SWIFT_HARDENING_EXPONENT] = 0.1;
    point.SolveUniaxialStress(0.02);
    point.Commit();
    const double p = UniaxialReferenceP(point.Material, 0.02);
    KRATOS_EXPECT_NEAR(point.P(), p, 1.0e-12);
    KRATOS_EXPECT_NEAR(point.Stress[0], 210000.0 * (0.02 - p), 1.0e-7);
}

KRATOS_TEST_CASE_IN_SUITE(SwiftJ2AccumulatedPlasticStrainNonproportionalPath, KratosConstitutiveLawsFastSuite)
{
    SwiftMaterialPoint point;
    double accumulated = 0.0;
    const double path[][6] = {
        {0.01, -0.004, -0.006, 0.0, 0.0, 0.0},
        {0.005, -0.002, -0.003, 0.02, 0.0, 0.0},
        {-0.01, 0.006, 0.004, -0.01, 0.012, 0.0},
        {0.0, 0.0, 0.0, 0.0, -0.015, 0.02}};
    for (const auto& r_step : path) {
        const Vector previous = point.PlasticStrain();
        for (unsigned int i = 0; i < 6; ++i) {
            point.Strain[i] = r_step[i];
        }
        point.Commit();
        const Vector increment = point.PlasticStrain() - previous;
        double squared_norm = 0.0;
        for (unsigned int i = 0; i < 6; ++i) {
            squared_norm += increment[i] * increment[i] * (i < 3 ? 1.0 : 0.5);
        }
        accumulated += std::sqrt(2.0 / 3.0 * squared_norm);
        KRATOS_EXPECT_NEAR(point.P(), accumulated, 1.0e-14);
        const Vector plastic = point.PlasticStrain();
        KRATOS_EXPECT_NEAR(plastic[0] + plastic[1] + plastic[2], 0.0, 1.0e-14);
    }
    const Vector plastic = point.PlasticStrain();
    double squared_norm = 0.0;
    for (unsigned int i = 0; i < 6; ++i) {
        squared_norm += plastic[i] * plastic[i] * (i < 3 ? 1.0 : 0.5);
    }
    KRATOS_EXPECT_GT(point.P(), std::sqrt(2.0 / 3.0 * squared_norm));
}

KRATOS_TEST_CASE_IN_SUITE(SwiftJ2TangentElastic, KratosConstitutiveLawsFastSuite)
{
    SwiftMaterialPoint point;
    point.Strain[0] = 1.0e-5;
    point.Strain[4] = 2.0e-5;
    CheckTangent(point);
}

KRATOS_TEST_CASE_IN_SUITE(SwiftJ2TangentJustPlastic, KratosConstitutiveLawsFastSuite)
{
    SwiftMaterialPoint point;
    point.Strain[3] = 1.001 * 1000.0 * std::pow(0.01, 0.2) / (std::sqrt(3.0) * 210000.0 / 2.6);
    CheckTangent(point);
}

KRATOS_TEST_CASE_IN_SUITE(SwiftJ2TangentPlasticWithHistory, KratosConstitutiveLawsFastSuite)
{
    for (const double n : {0.2, 1.0, 2.0}) {
        SwiftMaterialPoint point;
        point.Material[SWIFT_HARDENING_EXPONENT] = n;
        point.Strain[0] = 0.004;
        point.Commit();
        const double strain[] = {0.008, 0.0005, -0.002, 0.012, 0.004, -0.005};
        for (unsigned int i = 0; i < 6; ++i) {
            point.Strain[i] = strain[i];
        }
        CheckTangent(point);
    }
}

KRATOS_TEST_CASE_IN_SUITE(SwiftJ2ClonePreservesHistory, KratosConstitutiveLawsFastSuite)
{
    SwiftMaterialPoint point;
    point.SolveUniaxialStress(0.02);
    point.Commit();
    auto p_clone = point.Law.Clone();
    CheckContinuation(point, *p_clone);
    point.Strain[0] += 0.01;
    point.Commit();
    double clone_p = 0.0;
    KRATOS_EXPECT_GT(point.P(), p_clone->GetValue(ACCUMULATED_PLASTIC_STRAIN, clone_p));
}

KRATOS_TEST_CASE_IN_SUITE(SwiftJ2SerializationPreservesHistory, KratosConstitutiveLawsFastSuite)
{
    SwiftMaterialPoint point;
    point.SolveUniaxialStress(0.02);
    point.Commit();
    ConstitutiveLaw::Pointer p_original = point.Law.Clone();
    StreamSerializer serializer;
    serializer.save("Law", p_original);
    ConstitutiveLaw::Pointer p_loaded;
    serializer.load("Law", p_loaded);
    KRATOS_EXPECT_EQ(p_loaded->Info(), "SmallStrainSwiftJ2Plasticity3DLaw");
    CheckContinuation(point, *p_loaded);
}

KRATOS_TEST_CASE_IN_SUITE(SwiftJ2MaterialRegistration, KratosConstitutiveLawsFastSuite)
{
    SwiftMaterialPoint point;
    Parameters materials(R"({"properties": [{"model_part_name": "Main", "properties_id": 2,
        "Material": {"constitutive_law": {"name": "SmallStrainSwiftJ2Plasticity3DLaw"},
        "Variables": {"YOUNG_MODULUS": 210000.0, "POISSON_RATIO": 0.3,
        "SWIFT_COEFFICIENT": 1000.0, "SWIFT_INITIAL_STRAIN": 0.01,
        "SWIFT_HARDENING_EXPONENT": 0.2}, "Tables": {}}}]})");
    ReadMaterialsUtility reader(point.CurrentModel);
    reader.ReadMaterials(materials);
    auto& r_material = point.rModelPart.GetProperties(2);
    auto p_law = r_material[CONSTITUTIVE_LAW];
    KRATOS_EXPECT_EQ(p_law->Info(), "SmallStrainSwiftJ2Plasticity3DLaw");
    KRATOS_EXPECT_EQ(p_law->WorkingSpaceDimension(), 3);
    KRATOS_EXPECT_EQ(p_law->GetStrainSize(), 6);
    KRATOS_EXPECT_EQ(p_law->Check(r_material, *point.pGeometry, point.rModelPart.GetProcessInfo()), 0);
    ConstitutiveLaw::Features features;
    p_law->GetLawFeatures(features);
    KRATOS_EXPECT_TRUE(features.mOptions.Is(ConstitutiveLaw::INFINITESIMAL_STRAINS));
    KRATOS_EXPECT_TRUE(p_law->Has(ACCUMULATED_PLASTIC_STRAIN));
    KRATOS_EXPECT_TRUE(p_law->Has(PLASTIC_STRAIN_VECTOR));
}

KRATOS_TEST_CASE_IN_SUITE(SwiftJ2InvalidParameters, KratosConstitutiveLawsFastSuite)
{
    SwiftMaterialPoint point;
    for (const auto* p_variable : {&YOUNG_MODULUS, &SWIFT_COEFFICIENT,
            &SWIFT_INITIAL_STRAIN, &SWIFT_HARDENING_EXPONENT}) {
        const double original = point.Material[*p_variable];
        for (const double invalid : {0.0, -1.0, std::numeric_limits<double>::infinity(),
                std::numeric_limits<double>::quiet_NaN()}) {
            point.Material[*p_variable] = invalid;
            KRATOS_EXPECT_EXCEPTION_IS_THROWN(point.Law.Check(point.Material, *point.pGeometry,
                point.rModelPart.GetProcessInfo()), "must be finite and positive");
        }
        point.Material[*p_variable] = original;
    }
    for (const double invalid : {-1.0, 0.5, std::numeric_limits<double>::quiet_NaN()}) {
        point.Material[POISSON_RATIO] = invalid;
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(point.Law.Check(point.Material, *point.pGeometry,
            point.rModelPart.GetProcessInfo()), "POISSON_RATIO must satisfy");
    }
    point.Material[POISSON_RATIO] = 0.3;
    point.Material[SWIFT_HARDENING_EXPONENT] = 1000.0;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(point.Evaluate(), "not representable");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(point.Law.SetValue(ACCUMULATED_PLASTIC_STRAIN,
        -0.1, point.rModelPart.GetProcessInfo()), "p must be finite and non-negative");
}

KRATOS_TEST_CASE_IN_SUITE(SwiftJ2OutputFlagsAndInfinitesimalStrain, KratosConstitutiveLawsFastSuite)
{
    SwiftMaterialPoint point;
    point.Values.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
    point.F(0, 0) += 1.0e-5;
    point.F(0, 1) = 2.0e-5;
    point.F(1, 0) = 3.0e-5;
    point.Evaluate();
    KRATOS_EXPECT_NEAR(point.Strain[0], 1.0e-5, 1.0e-15);
    KRATOS_EXPECT_NEAR(point.Strain[3], 5.0e-5, 1.0e-15);
    const Vector expected = prod(ElasticReference(point.Material), point.Strain);
    KRATOS_EXPECT_VECTOR_NEAR(point.Stress, expected, 1.0e-10);
    point.Values.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    point.Strain = ZeroVector(6);
    point.Strain[0] = point.Strain[1] = point.Strain[2] = 0.001;
    point.Commit();
    KRATOS_EXPECT_EQ(point.P(), 0.0); // Hydrostatic: zero deviator, no unsafe normalization.
    point.Strain[3] = 0.02;
    point.Stress = ScalarVector(6, -123.0);
    point.Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS, false);
    point.Evaluate();
    KRATOS_EXPECT_VECTOR_NEAR(point.Stress, ScalarVector(6, -123.0), 0.0);
    KRATOS_EXPECT_EQ(point.P(), 0.0);
    const Matrix tangent = point.Tangent;
    point.Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
    point.Commit();
    KRATOS_EXPECT_GT(point.P(), 0.0);
    KRATOS_EXPECT_MATRIX_NEAR(point.Tangent, tangent, 0.0);
}

} // namespace Kratos::Testing
