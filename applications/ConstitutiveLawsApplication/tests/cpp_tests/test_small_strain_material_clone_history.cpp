//    |  /           |
//    ' /   __|  _ \ __|
//    . \  (     |\__ `
//   _|\_\___|_| _|___/
//                   Multi-Physics
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

#include "constitutive_laws_fast_suite.h"
#include "constitutive_laws_application_variables.h"
#include "custom_constitutive/small_strains/damage/small_strain_isotropic_damage_3d.h"
#include "custom_constitutive/small_strains/plasticity/small_strain_j2_plasticity_3d.h"
#include "geometries/tetrahedra_3d_4.h"
#include "includes/kratos_components.h"
#include "includes/mat_variables.h"

namespace Kratos::Testing
{
namespace
{
class CloneHistoryMaterialPoint
{
public:
    CloneHistoryMaterialPoint()
    {
        auto& r_options = mParameters.GetOptions();
        r_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
        r_options.Set(ConstitutiveLaw::COMPUTE_STRESS);
        r_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
        mParameters.SetShapeFunctionsValues(mShapeFunctions);
        mParameters.SetStrainVector(mStrain);
        mParameters.SetStressVector(mStress);
        mParameters.SetConstitutiveMatrix(mTangent);
        mProperties.SetValue(POISSON_RATIO, 0.3);
    }

    void Initialize(ConstitutiveLaw& rLaw)
    {
        rLaw.Check(mProperties, mGeometry, mProcessInfo);
        rLaw.InitializeMaterial(mProperties, mGeometry, mShapeFunctions);
    }

    void Advance(ConstitutiveLaw& rLaw, const double AxialStrain)
    {
        SetStrain(AxialStrain);
        rLaw.InitializeMaterialResponseCauchy(mParameters);
        rLaw.FinalizeMaterialResponseCauchy(mParameters);
    }

    void ExpectSameResponse(ConstitutiveLaw& rOriginal, ConstitutiveLaw& rClone, const double AxialStrain)
    {
        SetStrain(AxialStrain);
        rOriginal.CalculateMaterialResponsePK2(mParameters);
        const Vector stress = mStress;
        const Matrix tangent = mTangent;

        SetStrain(AxialStrain);
        rClone.CalculateMaterialResponsePK2(mParameters);
        KRATOS_EXPECT_VECTOR_NEAR(mStress, stress, 1e-10);
        KRATOS_EXPECT_MATRIX_NEAR(mTangent, tangent, 1e-10);
    }

    Properties mProperties{1};

private:
    void SetStrain(const double AxialStrain)
    {
        noalias(mStrain) = ZeroVector(6);
        mStrain[0] = AxialStrain;
    }

    Tetrahedra3D4<Node> mGeometry{
        Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0),
        Kratos::make_intrusive<Node>(2, 1.0, 0.0, 0.0),
        Kratos::make_intrusive<Node>(3, 0.0, 1.0, 0.0),
        Kratos::make_intrusive<Node>(4, 0.0, 0.0, 1.0)};
    ProcessInfo mProcessInfo;
    const Vector mShapeFunctions = ScalarVector(4, 0.25);
    Vector mStrain = ZeroVector(6);
    Vector mStress = ZeroVector(6);
    Matrix mTangent = ZeroMatrix(6, 6);
    ConstitutiveLaw::Parameters mParameters{mGeometry, mProperties, mProcessInfo};
};

// Read the actual registered law's protected vector before exercising its response.
// A broken clone has an empty vector and an indeterminate accumulated strain.
struct J2CloneHistoryAccess : SmallStrainJ2Plasticity3D
{
    static const Vector& PlasticStrain(const SmallStrainJ2Plasticity3D& rLaw)
    {
        return rLaw.*(&J2CloneHistoryAccess::mPlasticStrain);
    }
};
} // namespace

KRATOS_TEST_CASE_IN_SUITE(SmallStrainJ2Plasticity3DClonePreservesHistory, KratosConstitutiveLawsFastSuite)
{
    CloneHistoryMaterialPoint point;
    point.mProperties.SetValue(YOUNG_MODULUS, 200000.0);
    point.mProperties.SetValue(YIELD_STRESS, 250.0);
    point.mProperties.SetValue(ISOTROPIC_HARDENING_MODULUS, 1000.0);
    point.mProperties.SetValue(EXPONENTIAL_SATURATION_YIELD_STRESS, 250.0);
    point.mProperties.SetValue(HARDENING_EXPONENT, 0.0);

    const auto p_original = KratosComponents<ConstitutiveLaw>::Get("SmallStrainJ2Plasticity3DLaw").Clone();
    point.Initialize(*p_original);
    point.Advance(*p_original, 0.01);
    const auto& r_original = dynamic_cast<const SmallStrainJ2Plasticity3D&>(*p_original);
    const Vector plastic_strain = J2CloneHistoryAccess::PlasticStrain(r_original);
    double accumulated_strain = 0.0;
    p_original->GetValue(ACCUMULATED_PLASTIC_STRAIN, accumulated_strain);
    KRATOS_EXPECT_GT(accumulated_strain, 0.0);
    KRATOS_EXPECT_GT(norm_2(plastic_strain), 0.0);

    const auto p_clone = p_original->Clone();
    const auto& r_clone = dynamic_cast<const SmallStrainJ2Plasticity3D&>(*p_clone);
    const auto& r_clone_plastic_strain = J2CloneHistoryAccess::PlasticStrain(r_clone);
    ASSERT_EQ(r_clone_plastic_strain.size(), plastic_strain.size());
    KRATOS_EXPECT_VECTOR_NEAR(r_clone_plastic_strain, plastic_strain, 1e-12);
    double clone_accumulated_strain = 0.0;
    p_clone->GetValue(ACCUMULATED_PLASTIC_STRAIN, clone_accumulated_strain);
    KRATOS_EXPECT_NEAR(clone_accumulated_strain, accumulated_strain, 1e-12);
    point.ExpectSameResponse(*p_original, *p_clone, 0.009);

    point.Advance(*p_clone, 0.02);
    p_clone->GetValue(ACCUMULATED_PLASTIC_STRAIN, clone_accumulated_strain);
    KRATOS_EXPECT_GT(clone_accumulated_strain, accumulated_strain);
    KRATOS_EXPECT_GT(norm_2(r_clone_plastic_strain - plastic_strain), 0.0);
    double original_accumulated_strain = 0.0;
    p_original->GetValue(ACCUMULATED_PLASTIC_STRAIN, original_accumulated_strain);
    KRATOS_EXPECT_NEAR(original_accumulated_strain, accumulated_strain, 1e-12);
    const auto& r_original_plastic_strain = J2CloneHistoryAccess::PlasticStrain(r_original);
    KRATOS_EXPECT_VECTOR_NEAR(r_original_plastic_strain, plastic_strain, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(SmallStrainIsotropicDamage3DClonePreservesHistory, KratosConstitutiveLawsFastSuite)
{
    CloneHistoryMaterialPoint point;
    point.mProperties.SetValue(YOUNG_MODULUS, 30000.0);
    point.mProperties.SetValue(HARDENING_CURVE, 0);
    Vector stress_limits(2);
    stress_limits[0] = 3.0;
    stress_limits[1] = 4.0;
    point.mProperties.SetValue(STRESS_LIMITS, stress_limits);
    const Vector hardening_parameters = ScalarVector(1, 0.5);
    point.mProperties.SetValue(HARDENING_PARAMETERS, hardening_parameters);

    const auto p_original = KratosComponents<ConstitutiveLaw>::Get("SmallStrainIsotropicDamage3DLaw").Clone();
    point.Initialize(*p_original);
    Vector initial_history;
    p_original->GetValue(INTERNAL_VARIABLES, initial_history);
    ASSERT_EQ(initial_history.size(), 1);
    point.Advance(*p_original, 0.001);
    Vector history;
    p_original->GetValue(INTERNAL_VARIABLES, history);
    ASSERT_EQ(history.size(), 1);
    KRATOS_EXPECT_GT(history[0], initial_history[0]);

    const auto p_clone = p_original->Clone();
    KRATOS_EXPECT_TRUE(dynamic_cast<const SmallStrainIsotropicDamage3D*>(p_clone.get()) != nullptr);
    Vector clone_history;
    p_clone->GetValue(INTERNAL_VARIABLES, clone_history);
    ASSERT_EQ(clone_history.size(), 1);
    ASSERT_NEAR(clone_history[0], history[0], 1e-12);
    point.ExpectSameResponse(*p_original, *p_clone, 0.0008);

    point.Advance(*p_clone, 0.002);
    p_clone->GetValue(INTERNAL_VARIABLES, clone_history);
    KRATOS_EXPECT_GT(clone_history[0], history[0]);
    Vector original_history;
    p_original->GetValue(INTERNAL_VARIABLES, original_history);
    KRATOS_EXPECT_VECTOR_NEAR(original_history, history, 1e-12);
}
} // namespace Kratos::Testing
