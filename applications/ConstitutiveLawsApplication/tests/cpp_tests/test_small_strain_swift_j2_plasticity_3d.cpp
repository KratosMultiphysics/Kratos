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
    Model mModel;
    ModelPart& mrModelPart;
    Properties mMaterial{1};
    Geometry<Node>::Pointer mpGeometry;
    Vector mShapeFunctions = ScalarVector(4, 0.25);
    Vector mStrain = ZeroVector(6);
    Vector mStress = ZeroVector(6);
    Matrix mTangent = ZeroMatrix(6, 6);
    Matrix mF = IdentityMatrix(3);
    ConstitutiveLaw::Parameters mValues;
    SmallStrainSwiftJ2Plasticity3D mLaw;

    SwiftMaterialPoint() : mrModelPart(mModel.CreateModelPart("Main"))
    {
        mpGeometry = Kratos::make_shared<Tetrahedra3D4<Node>>(
            mrModelPart.CreateNewNode(1, 0.0, 0.0, 0.0),
            mrModelPart.CreateNewNode(2, 1.0, 0.0, 0.0),
            mrModelPart.CreateNewNode(3, 0.0, 1.0, 0.0),
            mrModelPart.CreateNewNode(4, 0.0, 0.0, 1.0));
        mMaterial.SetValue(YOUNG_MODULUS, 210000.0); // MPa
        mMaterial.SetValue(POISSON_RATIO, 0.3);
        mMaterial.SetValue(SWIFT_COEFFICIENT, 1000.0); // MPa
        mMaterial.SetValue(SWIFT_INITIAL_STRAIN, 0.01);
        mMaterial.SetValue(SWIFT_HARDENING_EXPONENT, 0.2);
        mValues.SetElementGeometry(*mpGeometry);
        mValues.SetMaterialProperties(mMaterial);
        mValues.SetProcessInfo(mrModelPart.GetProcessInfo());
        mValues.SetShapeFunctionsValues(mShapeFunctions);
        mValues.SetStrainVector(mStrain);
        mValues.SetStressVector(mStress);
        mValues.SetConstitutiveMatrix(mTangent);
        mValues.SetDeformationGradientF(mF);
        mValues.SetDeterminantF(1.0);
        mValues.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
        mValues.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
        mValues.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
        mLaw.InitializeMaterial(mMaterial, *mpGeometry, mShapeFunctions);
    }

    double P()
    {
        double p = 0.0;
        return mLaw.GetValue(ACCUMULATED_PLASTIC_STRAIN, p);
    }

    Vector PlasticStrain()
    {
        Vector strain;
        return mLaw.GetValue(PLASTIC_STRAIN_VECTOR, strain);
    }

    void Evaluate() { mLaw.CalculateMaterialResponseCauchy(mValues); }
    void Commit() { mLaw.FinalizeMaterialResponseCauchy(mValues); }

    // Each iteration is a trial from the same committed state. Only the caller commits.
    void SolveUniaxialStress(const double AxialStrain)
    {
        mStrain[0] = AxialStrain;
        for (unsigned int iteration = 0; iteration < 30; ++iteration) {
            Evaluate();
            if (std::hypot(mStress[1], mStress[2]) < 1.0e-8) {
                return;
            }
            const double det = mTangent(1, 1) * mTangent(2, 2) - mTangent(1, 2) * mTangent(2, 1);
            const double dy = (mTangent(2, 2) * mStress[1] - mTangent(1, 2) * mStress[2]) / det;
            const double dz = (mTangent(1, 1) * mStress[2] - mTangent(2, 1) * mStress[1]) / det;
            mStrain[1] -= dy;
            mStrain[2] -= dz;
        }
        KRATOS_ERROR << "Uniaxial-stress test could not equilibrate transverse stresses";
    }
};

double von_mises(const Vector& rStress)
{
    return std::sqrt(0.5 * (std::pow(rStress[0] - rStress[1], 2) +
        std::pow(rStress[1] - rStress[2], 2) + std::pow(rStress[2] - rStress[0], 2)) +
        3.0 * (rStress[3] * rStress[3] + rStress[4] * rStress[4] + rStress[5] * rStress[5]));
}

Matrix elastic_reference(const Properties& rMaterial)
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
double uniaxial_reference_p(const Properties& rMaterial, const double AxialStrain)
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

void check_tangent(SwiftMaterialPoint& rPoint)
{
    rPoint.Evaluate();
    const Matrix analytical = rPoint.mTangent;
    const Vector strain = rPoint.mStrain;
    const double committed_p = rPoint.P();
    const Vector committed_plastic_strain = rPoint.PlasticStrain();
    Matrix numerical = ZeroMatrix(6, 6);
    // h is below the distance to the yield switch, including the just-plastic
    // case. A 3e-7*E absolute tolerance allows residual/h and roundoff error
    // without hiding relative errors in the dominant stiffness terms.
    constexpr double h = 1.0e-8;
    rPoint.mValues.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
    for (unsigned int j = 0; j < 6; ++j) {
        rPoint.mStrain = strain;
        rPoint.mStrain[j] += h;
        rPoint.Evaluate();
        const Vector plus = rPoint.mStress;
        rPoint.mStrain[j] -= 2.0 * h;
        rPoint.Evaluate();
        for (unsigned int i = 0; i < 6; ++i) {
            numerical(i, j) = (plus[i] - rPoint.mStress[i]) / (2.0 * h);
            KRATOS_EXPECT_NEAR(analytical(i, j), analytical(j, i), 1.0e-10);
        }
    }
    KRATOS_EXPECT_MATRIX_NEAR(analytical, numerical, 3.0e-7 * rPoint.mMaterial[YOUNG_MODULUS]);
    KRATOS_EXPECT_EQ(rPoint.P(), committed_p);
    KRATOS_EXPECT_VECTOR_NEAR(rPoint.PlasticStrain(), committed_plastic_strain, 0.0);
    rPoint.mStrain = strain;
    rPoint.mValues.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
}

void check_continuation(SwiftMaterialPoint& rPoint, ConstitutiveLaw& rOther)
{
    double other_p = 0.0;
    Vector other_plastic;
    KRATOS_EXPECT_NEAR(rOther.GetValue(ACCUMULATED_PLASTIC_STRAIN, other_p), rPoint.P(), 1.0e-15);
    KRATOS_EXPECT_VECTOR_NEAR(rOther.GetValue(PLASTIC_STRAIN_VECTOR, other_plastic), rPoint.PlasticStrain(), 1.0e-15);
    rPoint.mStrain[0] += 0.005;
    rPoint.mStrain[4] += 0.007;
    rPoint.Evaluate();
    const Vector expected_stress = rPoint.mStress;
    const Matrix expected_tangent = rPoint.mTangent;
    rOther.CalculateMaterialResponseCauchy(rPoint.mValues);
    KRATOS_EXPECT_VECTOR_NEAR(rPoint.mStress, expected_stress, 1.0e-10);
    KRATOS_EXPECT_MATRIX_NEAR(rPoint.mTangent, expected_tangent, 1.0e-10);
    rPoint.Commit();
    rOther.FinalizeMaterialResponseCauchy(rPoint.mValues);
    KRATOS_EXPECT_NEAR(rOther.GetValue(ACCUMULATED_PLASTIC_STRAIN, other_p), rPoint.P(), 1.0e-15);
    KRATOS_EXPECT_VECTOR_NEAR(rOther.GetValue(PLASTIC_STRAIN_VECTOR, other_plastic), rPoint.PlasticStrain(), 1.0e-15);
}

} // namespace

KRATOS_TEST_CASE_IN_SUITE(SwiftJ2FlowStressEquation, KratosConstitutiveLawsFastSuite)
{
    SwiftMaterialPoint point;
    for (const double n : {0.2, 1.0, 2.0}) {
        point.mMaterial[SWIFT_HARDENING_EXPONENT] = n;
        for (const double p : {0.0, 0.001, 0.02, 0.1, 0.5}) {
            point.mLaw.SetValue(ACCUMULATED_PLASTIC_STRAIN, p, point.mrModelPart.GetProcessInfo());
            double flow_stress = 0.0;
            point.mLaw.CalculateValue(point.mValues, YIELD_STRESS, flow_stress);
            KRATOS_EXPECT_NEAR(flow_stress, 1000.0 * std::pow(0.01 + p, n), 1.0e-11);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(SwiftJ2InitialYieldStress, KratosConstitutiveLawsFastSuite)
{
    SwiftMaterialPoint point;
    double flow_stress = 0.0;
    point.mLaw.CalculateValue(point.mValues, YIELD_STRESS, flow_stress);
    // Independent numerical value of 1000*(0.01)^0.2 in MPa.
    KRATOS_EXPECT_NEAR(flow_stress, 398.1071705534972, 1.0e-10);
    KRATOS_EXPECT_EQ(point.P(), 0.0);
}

KRATOS_TEST_CASE_IN_SUITE(SwiftJ2ElasticResponse, KratosConstitutiveLawsFastSuite)
{
    SwiftMaterialPoint point;
    for (unsigned int i = 0; i < 6; ++i) {
        point.mStrain[i] = (static_cast<double>(i) - 2.0) * 1.0e-6;
    }
    point.Evaluate();
    const Matrix elastic = elastic_reference(point.mMaterial);
    const Vector expected = prod(elastic, point.mStrain);
    KRATOS_EXPECT_VECTOR_NEAR(point.mStress, expected, 1.0e-12);
    KRATOS_EXPECT_MATRIX_NEAR(point.mTangent, elastic, 1.0e-10);
    point.Commit();
    KRATOS_EXPECT_EQ(point.P(), 0.0);
}

KRATOS_TEST_CASE_IN_SUITE(SwiftJ2ElasticPlasticTransitionAndPureShear, KratosConstitutiveLawsFastSuite)
{
    SwiftMaterialPoint point;
    const double mu = 210000.0 / 2.6;
    const double yield = 1000.0 * std::pow(0.01, 0.2);
    for (const double factor : {0.5, 0.999, 1.001, 1.1, 2.0, 5.0}) {
        point.mStrain[3] = factor * yield / (std::sqrt(3.0) * mu);
        const double previous_p = point.P();
        point.Evaluate();
        KRATOS_EXPECT_EQ(point.P(), previous_p);
        point.Commit();
        const double q = von_mises(point.mStress);
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
    const double young = point.mMaterial[YOUNG_MODULUS];
    const double nu = point.mMaterial[POISSON_RATIO];
    const double k = point.mMaterial[SWIFT_COEFFICIENT];
    const double eps0 = point.mMaterial[SWIFT_INITIAL_STRAIN];
    const double n = point.mMaterial[SWIFT_HARDENING_EXPONENT];
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
    point.mStrain[3] = gamma; // Only engineering xy shear is prescribed.
    point.Evaluate();
    KRATOS_EXPECT_NEAR(point.mStress[3], tau_ref, stress_tolerance);
    KRATOS_EXPECT_EQ(point.P(), initial_p);
    KRATOS_EXPECT_VECTOR_NEAR(point.PlasticStrain(), initial_plastic_strain, 0.0);

    point.Commit();
    const double p = point.P();
    const Vector plastic_strain = point.PlasticStrain();
    KRATOS_EXPECT_NEAR(point.mStress[3], tau_ref, stress_tolerance);
    KRATOS_EXPECT_NEAR(p, p_ref, strain_tolerance);
    KRATOS_EXPECT_NEAR(plastic_strain[3], gamma_p_ref, strain_tolerance);
    KRATOS_EXPECT_NEAR(plastic_strain[3], sqrt_three * p, strain_tolerance);
    KRATOS_EXPECT_NEAR(p, std::abs(plastic_strain[3]) / sqrt_three, strain_tolerance);
    KRATOS_EXPECT_NEAR(sqrt_three * std::abs(point.mStress[3]), sigma_eq_ref, stress_tolerance);
    KRATOS_EXPECT_NEAR(sqrt_three * std::abs(point.mStress[3]),
        k * std::pow(eps0 + p, n), stress_tolerance);
    for (unsigned int i = 0; i < 6; ++i) {
        if (i != 3) {
            KRATOS_EXPECT_NEAR(point.mStress[i], 0.0, stress_tolerance);
            KRATOS_EXPECT_NEAR(plastic_strain[i], 0.0, strain_tolerance);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(SwiftJ2MonotonicUniaxialStress, KratosConstitutiveLawsFastSuite)
{
    for (const double n : {0.2, 1.0, 2.0}) {
        SwiftMaterialPoint point;
        point.mMaterial[SWIFT_HARDENING_EXPONENT] = n;
        for (const double axial : {0.0001, 0.001, 0.003, 0.006, 0.01, 0.02}) {
            point.SolveUniaxialStress(axial);
            point.Commit();
            const double p = uniaxial_reference_p(point.mMaterial, axial);
            KRATOS_EXPECT_NEAR(point.mStress[0], 210000.0 * (axial - p), 1.0e-7);
            KRATOS_EXPECT_NEAR(point.mStress[1], 0.0, 1.0e-8);
            KRATOS_EXPECT_NEAR(point.mStress[2], 0.0, 1.0e-8);
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
    const double stress = point.mStress[0];
    const double p = point.P();
    const Vector plastic = point.PlasticStrain();
    const double decrement = 0.5 * stress / 210000.0;
    point.SolveUniaxialStress(0.02 - decrement);
    point.Commit();
    KRATOS_EXPECT_NEAR(stress - point.mStress[0], 210000.0 * decrement, 1.0e-8);
    KRATOS_EXPECT_EQ(point.P(), p);
    KRATOS_EXPECT_VECTOR_NEAR(point.PlasticStrain(), plastic, 0.0);
    point.SolveUniaxialStress(p);
    point.Commit();
    KRATOS_EXPECT_NEAR(point.mStress[0], 0.0, 1.0e-8);
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
        const double p = uniaxial_reference_p(compression.mMaterial, axial);
        KRATOS_EXPECT_NEAR(compression.mStress[0], -210000.0 * (axial - p), 1.0e-7);
        KRATOS_EXPECT_NEAR(compression.mStress[1], 0.0, 1.0e-8);
        KRATOS_EXPECT_NEAR(compression.mStress[2], 0.0, 1.0e-8);
        KRATOS_EXPECT_NEAR(compression.P(), p, 1.0e-12);
        const Vector opposite = -tension.PlasticStrain();
        KRATOS_EXPECT_VECTOR_NEAR(compression.PlasticStrain(), opposite, 1.0e-14);
        KRATOS_EXPECT_MATRIX_NEAR(compression.mTangent, tension.mTangent, 1.0e-8);
    }
}

KRATOS_TEST_CASE_IN_SUITE(SwiftJ2StronglyNonlinearHardening, KratosConstitutiveLawsFastSuite)
{
    SwiftMaterialPoint point;
    point.mMaterial[SWIFT_INITIAL_STRAIN] = 1.0e-12;
    point.mMaterial[SWIFT_HARDENING_EXPONENT] = 0.1;
    point.SolveUniaxialStress(0.02);
    point.Commit();
    const double p = uniaxial_reference_p(point.mMaterial, 0.02);
    KRATOS_EXPECT_NEAR(point.P(), p, 1.0e-12);
    KRATOS_EXPECT_NEAR(point.mStress[0], 210000.0 * (0.02 - p), 1.0e-7);
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
            point.mStrain[i] = r_step[i];
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
    point.mStrain[0] = 1.0e-5;
    point.mStrain[4] = 2.0e-5;
    check_tangent(point);
}

KRATOS_TEST_CASE_IN_SUITE(SwiftJ2TangentJustPlastic, KratosConstitutiveLawsFastSuite)
{
    SwiftMaterialPoint point;
    point.mStrain[3] = 1.001 * 1000.0 * std::pow(0.01, 0.2) / (std::sqrt(3.0) * 210000.0 / 2.6);
    check_tangent(point);
}

KRATOS_TEST_CASE_IN_SUITE(SwiftJ2TangentPlasticWithHistory, KratosConstitutiveLawsFastSuite)
{
    for (const double n : {0.2, 1.0, 2.0}) {
        SwiftMaterialPoint point;
        point.mMaterial[SWIFT_HARDENING_EXPONENT] = n;
        point.mStrain[0] = 0.004;
        point.Commit();
        const double strain[] = {0.008, 0.0005, -0.002, 0.012, 0.004, -0.005};
        for (unsigned int i = 0; i < 6; ++i) {
            point.mStrain[i] = strain[i];
        }
        check_tangent(point);
    }
}

KRATOS_TEST_CASE_IN_SUITE(SwiftJ2ClonePreservesHistory, KratosConstitutiveLawsFastSuite)
{
    SwiftMaterialPoint point;
    point.SolveUniaxialStress(0.02);
    point.Commit();
    auto p_clone = point.mLaw.Clone();
    check_continuation(point, *p_clone);
    point.mStrain[0] += 0.01;
    point.Commit();
    double clone_p = 0.0;
    KRATOS_EXPECT_GT(point.P(), p_clone->GetValue(ACCUMULATED_PLASTIC_STRAIN, clone_p));
}

KRATOS_TEST_CASE_IN_SUITE(SwiftJ2SerializationPreservesHistory, KratosConstitutiveLawsFastSuite)
{
    SwiftMaterialPoint point;
    point.SolveUniaxialStress(0.02);
    point.Commit();
    ConstitutiveLaw::Pointer p_original = point.mLaw.Clone();
    StreamSerializer serializer;
    serializer.save("Law", p_original);
    ConstitutiveLaw::Pointer p_loaded;
    serializer.load("Law", p_loaded);
    KRATOS_EXPECT_EQ(p_loaded->Info(), "SmallStrainSwiftJ2Plasticity3DLaw");
    check_continuation(point, *p_loaded);
}

KRATOS_TEST_CASE_IN_SUITE(SwiftJ2MaterialRegistration, KratosConstitutiveLawsFastSuite)
{
    SwiftMaterialPoint point;
    Parameters materials(R"({"properties": [{"model_part_name": "Main", "properties_id": 2,
        "Material": {"constitutive_law": {"name": "SmallStrainSwiftJ2Plasticity3DLaw"},
        "Variables": {"YOUNG_MODULUS": 210000.0, "POISSON_RATIO": 0.3,
        "SWIFT_COEFFICIENT": 1000.0, "SWIFT_INITIAL_STRAIN": 0.01,
        "SWIFT_HARDENING_EXPONENT": 0.2}, "Tables": {}}}]})");
    ReadMaterialsUtility reader(point.mModel);
    reader.ReadMaterials(materials);
    auto& r_material = point.mrModelPart.GetProperties(2);
    auto p_law = r_material[CONSTITUTIVE_LAW];
    KRATOS_EXPECT_EQ(p_law->Info(), "SmallStrainSwiftJ2Plasticity3DLaw");
    KRATOS_EXPECT_EQ(p_law->WorkingSpaceDimension(), 3);
    KRATOS_EXPECT_EQ(p_law->GetStrainSize(), 6);
    KRATOS_EXPECT_EQ(p_law->Check(r_material, *point.mpGeometry, point.mrModelPart.GetProcessInfo()), 0);
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
        const double original = point.mMaterial[*p_variable];
        for (const double invalid : {0.0, -1.0, std::numeric_limits<double>::infinity(),
                std::numeric_limits<double>::quiet_NaN()}) {
            point.mMaterial[*p_variable] = invalid;
            KRATOS_EXPECT_EXCEPTION_IS_THROWN(point.mLaw.Check(point.mMaterial, *point.mpGeometry,
                point.mrModelPart.GetProcessInfo()), "must be finite and positive");
        }
        point.mMaterial[*p_variable] = original;
    }
    for (const double invalid : {-1.0, 0.5, std::numeric_limits<double>::quiet_NaN()}) {
        point.mMaterial[POISSON_RATIO] = invalid;
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(point.mLaw.Check(point.mMaterial, *point.mpGeometry,
            point.mrModelPart.GetProcessInfo()), "POISSON_RATIO must satisfy");
    }
    point.mMaterial[POISSON_RATIO] = 0.3;
    point.mMaterial[SWIFT_HARDENING_EXPONENT] = 1000.0;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(point.Evaluate(), "not representable");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(point.mLaw.SetValue(ACCUMULATED_PLASTIC_STRAIN,
        -0.1, point.mrModelPart.GetProcessInfo()), "p must be finite and non-negative");
}

KRATOS_TEST_CASE_IN_SUITE(SwiftJ2OutputFlagsAndInfinitesimalStrain, KratosConstitutiveLawsFastSuite)
{
    SwiftMaterialPoint point;
    point.mValues.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
    point.mF(0, 0) += 1.0e-5;
    point.mF(0, 1) = 2.0e-5;
    point.mF(1, 0) = 3.0e-5;
    point.Evaluate();
    KRATOS_EXPECT_NEAR(point.mStrain[0], 1.0e-5, 1.0e-15);
    KRATOS_EXPECT_NEAR(point.mStrain[3], 5.0e-5, 1.0e-15);
    const Vector expected = prod(elastic_reference(point.mMaterial), point.mStrain);
    KRATOS_EXPECT_VECTOR_NEAR(point.mStress, expected, 1.0e-10);
    point.mValues.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    point.mStrain = ZeroVector(6);
    point.mStrain[0] = point.mStrain[1] = point.mStrain[2] = 0.001;
    point.Commit();
    KRATOS_EXPECT_EQ(point.P(), 0.0); // Hydrostatic: zero deviator, no unsafe normalization.
    point.mStrain[3] = 0.02;
    point.mStress = ScalarVector(6, -123.0);
    point.mValues.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS, false);
    point.Evaluate();
    KRATOS_EXPECT_VECTOR_NEAR(point.mStress, ScalarVector(6, -123.0), 0.0);
    KRATOS_EXPECT_EQ(point.P(), 0.0);
    const Matrix tangent = point.mTangent;
    point.mValues.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
    point.Commit();
    KRATOS_EXPECT_GT(point.P(), 0.0);
    KRATOS_EXPECT_MATRIX_NEAR(point.mTangent, tangent, 0.0);
}

} // namespace Kratos::Testing
