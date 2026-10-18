//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt

// System includes
#include <array>
#include <limits>
#include <string>
#include <type_traits>

// Project includes
#include "containers/model.h"
#include "geometries/triangle_2d_3.h"
#include "includes/stream_serializer.h"
#include "tests/cpp_tests/constitutive_laws_fast_suite.h"

namespace Kratos::Testing
{
namespace
{

using PlaneStressIntegrator = GenericConstitutiveLawIntegratorPlasticity<
    PlaneStressVonMisesYieldSurface<PlaneStressVonMisesPlasticPotential<3>>>;
using SmallStrainPlaneStressLaw = GenericSmallStrainIsotropicPlasticityPlaneStress<PlaneStressIntegrator>;
using FiniteStrainPlaneStressLaw = GenericFiniteStrainIsotropicPlasticityPlaneStress<PlaneStressIntegrator>;

std::string serialize_plasticity_state(const ConstitutiveLaw& rLaw)
{
    StreamSerializer serializer(Serializer::SERIALIZER_TRACE_ERROR);
    serializer.pGetBuffer()->precision(std::numeric_limits<double>::max_digits10);
    if (rLaw.GetStrainSize() == 3) {
        // Observe the stored base state even when the derived serializer skips it.
        const auto& r_plasticity_base = dynamic_cast<const GenericSmallStrainIsotropicPlasticity<PlaneStressIntegrator>&>(rLaw);
        serializer.save_base("Law", r_plasticity_base);
    } else {
        serializer.save("Law", rLaw);
    }
    return serializer.GetStringRepresentation();
}

void check_plasticity_state(
    ConstitutiveLaw& rOriginal,
    ConstitutiveLaw& rRestored)
{
    Vector original_plastic_strain, restored_plastic_strain;
    rOriginal.GetValue(PLASTIC_STRAIN_VECTOR, original_plastic_strain);
    rRestored.GetValue(PLASTIC_STRAIN_VECTOR, restored_plastic_strain);
    KRATOS_EXPECT_VECTOR_NEAR(original_plastic_strain, restored_plastic_strain, 1.0e-14);

    double original_dissipation = 0.0, restored_dissipation = 0.0;
    rOriginal.GetValue(PLASTIC_DISSIPATION, original_dissipation);
    rRestored.GetValue(PLASTIC_DISSIPATION, restored_dissipation);
    KRATOS_EXPECT_NEAR(original_dissipation, restored_dissipation, 1.0e-14);

    // Threshold has no public GetValue interface. Compare the complete trace state,
    // including its named threshold field, at round-trip floating-point precision.
    const std::string original_state = serialize_plasticity_state(rOriginal);
    KRATOS_EXPECT_NE(original_state.find("\"Threshold\""), std::string::npos);
    KRATOS_EXPECT_EQ(original_state, serialize_plasticity_state(rRestored));
}

ConstitutiveLaw::Pointer round_trip_plasticity_law(
    const ConstitutiveLaw::Pointer& pOriginal,
    const Serializer::TraceType Trace)
{
    StreamSerializer serializer(Trace);
    serializer.pGetBuffer()->precision(std::numeric_limits<double>::max_digits10);
    auto p_other_law = KratosComponents<ConstitutiveLaw>::Get(
        "SmallStrainIsotropicPlasticity3DVonMisesVonMises").Clone();
    const double sentinel = 42.125;
    serializer.save("Law", pOriginal);
    serializer.save("OtherLaw", p_other_law);
    serializer.save("Sentinel", sentinel);

    ConstitutiveLaw::Pointer p_restored, p_other_restored;
    double restored_sentinel = 0.0;
    serializer.load("Law", p_restored);
    serializer.load("OtherLaw", p_other_restored);
    serializer.load("Sentinel", restored_sentinel);
    KRATOS_EXPECT_NE(pOriginal.get(), p_restored.get());
    KRATOS_EXPECT_EQ(p_restored->GetStrainSize(), pOriginal->GetStrainSize());
    check_plasticity_state(*p_other_law, *p_other_restored);
    KRATOS_EXPECT_EQ(restored_sentinel, sentinel);
    return p_restored;
}

template<class TLaw>
void check_plane_stress_plasticity_serialization(
    const std::string& rLawName,
    const Serializer::TraceType Trace,
    const bool PlasticHistory)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("Main");
    auto p_node_1 = r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    auto p_node_2 = r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    auto p_node_3 = r_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    Triangle2D3<Node> geometry(p_node_1, p_node_2, p_node_3);
    auto& r_process_info = r_model_part.GetProcessInfo();
    r_process_info[STEP] = 2;
    r_process_info[NL_ITERATION_NUMBER] = 1;

    Properties properties;
    properties[DENSITY] = 7850.0;
    properties[YOUNG_MODULUS] = 210.0e9;
    properties[POISSON_RATIO] = 0.22;
    properties[YIELD_STRESS] = 250.0e6;
    properties[FRACTURE_ENERGY] = 1.0e8;
    properties[HARDENING_CURVE] = 0;

    Vector strain = ZeroVector(3), stress = ZeroVector(3);
    Matrix tangent = ZeroMatrix(3, 3), deformation_gradient = IdentityMatrix(2);
    ConstitutiveLaw::Parameters parameters(geometry, properties, r_process_info);
    parameters.SetStrainVector(strain);
    parameters.SetStressVector(stress);
    parameters.SetConstitutiveMatrix(tangent);
    parameters.SetDeformationGradientF(deformation_gradient);
    auto& r_options = parameters.GetOptions();
    r_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, std::is_same_v<TLaw, SmallStrainPlaneStressLaw>);
    r_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    const auto set_loading = [&](const double Axial, const double Shear) {
        strain[0] = Axial;
        strain[1] = -0.2 * Axial;
        strain[2] = Shear;
        deformation_gradient = IdentityMatrix(2);
        deformation_gradient(0, 0) += Axial;
        deformation_gradient(1, 1) -= 0.2 * Axial;
        deformation_gradient(0, 1) = Shear;
        parameters.SetDeterminantF(MathUtils<double>::Det(deformation_gradient));
    };

    auto p_original = KratosComponents<ConstitutiveLaw>::Get(rLawName).Clone();
    KRATOS_EXPECT_EQ(p_original->Check(properties, geometry, r_process_info), 0);
    p_original->InitializeMaterial(properties, geometry, Vector(3, 1.0 / 3.0));
    // A committed axial/shear history; the elastic case uses the same path scaled down.
    const double scale = PlasticHistory ? 1.0 : 0.01;
    set_loading(scale * 0.002, scale * 0.0004);
    p_original->CalculateMaterialResponseCauchy(parameters);
    p_original->FinalizeMaterialResponseCauchy(parameters);
    set_loading(scale * 0.003, scale * 0.0008);
    p_original->CalculateMaterialResponseCauchy(parameters);
    p_original->FinalizeMaterialResponseCauchy(parameters);

    Vector plastic_strain;
    double dissipation = 0.0;
    p_original->GetValue(PLASTIC_STRAIN_VECTOR, plastic_strain);
    p_original->GetValue(PLASTIC_DISSIPATION, dissipation);
    if (PlasticHistory) {
        KRATOS_EXPECT_GT(norm_2(plastic_strain), 1.0e-4);
        KRATOS_EXPECT_GT(dissipation, 1.0e-6);
    } else {
        KRATOS_EXPECT_NEAR(norm_2(plastic_strain), 0.0, 1.0e-14);
        KRATOS_EXPECT_NEAR(dissipation, 0.0, 1.0e-14);
    }

    auto p_restored = round_trip_plasticity_law(p_original, Trace);
    auto p_cloned = p_original->Clone();
    TLaw copied_law(dynamic_cast<const TLaw&>(*p_original));
    check_plasticity_state(*p_original, *p_restored);
    check_plasticity_state(*p_original, *p_cloned);
    check_plasticity_state(*p_original, copied_law);

    // Unload, then reload with additional shear. Compare trial responses before
    // finalization, which may replace the computed tangent with the elastic matrix.
    const std::array<std::array<double, 2>, 3> continuation{{
        {{0.0025, 0.0006}}, {{0.0035, 0.0012}}, {{0.0038, -0.0004}}
    }};
    for (const auto& r_loading : continuation) {
        set_loading(scale * r_loading[0], scale * r_loading[1]);
        p_original->CalculateMaterialResponseCauchy(parameters);
        const Vector reference_stress = stress;
        const Matrix reference_tangent = tangent;
        KRATOS_EXPECT_GT(norm_frobenius(reference_tangent), 1.0e9);
        p_original->FinalizeMaterialResponseCauchy(parameters);
        for (auto* p_law : {p_restored.get(), p_cloned.get(), static_cast<ConstitutiveLaw*>(&copied_law)}) {
            set_loading(scale * r_loading[0], scale * r_loading[1]);
            p_law->CalculateMaterialResponseCauchy(parameters);
            // Absolute tolerances: 1e-4 Pa in stress and 1 Pa in the numerical tangent.
            KRATOS_EXPECT_MATRIX_NEAR(
                MathUtils<double>::StressVectorToTensor(stress),
                MathUtils<double>::StressVectorToTensor(reference_stress), 1.0e-4);
            KRATOS_EXPECT_MATRIX_NEAR(tangent, reference_tangent, 1.0);
            p_law->FinalizeMaterialResponseCauchy(parameters);
            check_plasticity_state(*p_original, *p_law);
        }
        p_restored = round_trip_plasticity_law(p_restored, Trace);
        check_plasticity_state(*p_original, *p_restored);
    }
}

} // namespace

KRATOS_TEST_CASE_IN_SUITE(SmallStrainPlaneStressPlasticitySerialization, KratosConstitutiveLawsFastSuite)
{
    for (const auto trace : {Serializer::SERIALIZER_NO_TRACE, Serializer::SERIALIZER_TRACE_ERROR}) {
        for (const bool plastic_history : {true, false}) {
            check_plane_stress_plasticity_serialization<SmallStrainPlaneStressLaw>(
                "SmallStrainIsotropicPlasticityPlaneStressVonMisesVonMises", trace, plastic_history);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(FiniteStrainPlaneStressPlasticitySerialization, KratosConstitutiveLawsFastSuite)
{
    for (const auto trace : {Serializer::SERIALIZER_NO_TRACE, Serializer::SERIALIZER_TRACE_ERROR}) {
        for (const bool plastic_history : {true, false}) {
            check_plane_stress_plasticity_serialization<FiniteStrainPlaneStressLaw>(
                "FiniteStrainIsotropicPlasticityPlaneStressVonMisesVonMises", trace, plastic_history);
        }
    }
}

} // namespace Kratos::Testing
