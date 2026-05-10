// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Gennady Markelov
//

#include "custom_constitutive/piecewise_linear_moment_capacity_constitutive_law.h"
#include "custom_utilities/check_utilities.hpp"
#include "custom_utilities/ublas_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

using namespace Kratos;

namespace
{

Properties CreateValidProperties()
{
    auto properties = Properties{};
    properties.SetValue(KAPPA_PIECEWISE_LINEAR_LAW, UblasUtilities::CreateVector({0.01, 0.03, 0.05}));
    properties.SetValue(MOMENTUM_PIECEWISE_LINEAR_LAW, UblasUtilities::CreateVector({80.0, 80.0, 128.0}));
    properties.SetValue(YOUNG_MODULUS, 80.0);
    properties.SetValue(POISSON_RATIO, 0.2);
    properties.SetValue(THICKNESS, 1.0);
    properties.SetValue(THICKNESS_EFFECTIVE_Y, 1.0);
    return properties;
}

double CalculateMomentForCurvature(PiecewiseLinearMomentCapacityConstitutiveLaw& rLaw,
                                   const Properties&                             rProperties,
                                   double                                        Curvature)
{
    auto parameters = ConstitutiveLaw::Parameters{};
    // Curvature is expected at index 1 in the constitutive law's strain vector
    auto strain_vector = UblasUtilities::CreateVector({0.0, Curvature, 0.0});
    parameters.SetStrainVector(strain_vector);

    const auto cl_strain_size = rLaw.GetStrainSize();
    auto       stress_vector  = Vector(cl_strain_size, 0.0);
    parameters.SetStressVector(stress_vector);

    auto constitutive_matrix = Matrix(cl_strain_size, cl_strain_size, 0.0);
    parameters.SetConstitutiveMatrix(constitutive_matrix);

    parameters.SetMaterialProperties(rProperties);
    parameters.Set(ConstitutiveLaw::COMPUTE_STRESS);

    rLaw.CalculateMaterialResponsePK2(parameters);
    // Moment stored at index 1 (Mz)
    return parameters.GetStressVector()[1];
}

double CalculateTangentForCurvature(PiecewiseLinearMomentCapacityConstitutiveLaw& rLaw,
                                    const Properties&                             rProperties,
                                    double                                        Curvature)
{
    auto parameters = ConstitutiveLaw::Parameters{};
    // Curvature is expected at index 1 in the constitutive law's strain vector
    auto strain_vector = UblasUtilities::CreateVector({0.0, Curvature, 0.0});
    parameters.SetStrainVector(strain_vector);
    parameters.SetMaterialProperties(rProperties);

    double tangent{};
    rLaw.CalculateValue(parameters, TANGENT_MODULUS, tangent);
    return tangent;
}

void FinalizeForCurvature(PiecewiseLinearMomentCapacityConstitutiveLaw& rLaw, const Properties& rProperties, double Curvature)
{
    auto parameters = ConstitutiveLaw::Parameters{};
    // Curvature is expected at index 1 in the constitutive law's strain vector
    auto strain_vector = UblasUtilities::CreateVector({0.0, Curvature, 0.0});
    parameters.SetStrainVector(strain_vector);

    const auto cl_strain_size = rLaw.GetStrainSize();
    auto       stress_vector  = Vector(cl_strain_size, 0.0);
    parameters.SetStressVector(stress_vector);

    auto constitutive_matrix = Matrix(cl_strain_size, cl_strain_size, 0.0);
    parameters.SetConstitutiveMatrix(constitutive_matrix);

    parameters.SetMaterialProperties(rProperties);
    parameters.Set(ConstitutiveLaw::COMPUTE_STRESS);

    rLaw.FinalizeMaterialResponsePK2(parameters);
}

} // namespace

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(PiecewiseLinearMomentCapacityConstitutiveLaw_CheckOfMomentCapacityLawThrowsWhenCurvatureVectorSizeIsIncorrect,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto law          = PiecewiseLinearMomentCapacityConstitutiveLaw{};
    auto       properties   = CreateValidProperties();
    const auto geometry     = Geometry<Node>{};
    const auto process_info = ProcessInfo{};

    properties.SetValue(KAPPA_PIECEWISE_LINEAR_LAW, UblasUtilities::CreateVector({0.01, 0.03}));

    KRATOS_EXPECT_EXCEPTION_IS_THROWN((void)law.Check(properties, geometry, process_info),
                                      "The number of entries in KAPPA_PIECEWISE_LINEAR_LAW (2) "
                                      "does not match MOMENTUM_PIECEWISE_LINEAR_LAW (3)")
}

KRATOS_TEST_CASE_IN_SUITE(PiecewiseLinearMomentCapacityConstitutiveLaw_MomentCapacityLawReturnsExpectedMomentForAllRegimes,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto       law   = PiecewiseLinearMomentCapacityConstitutiveLaw{};
    const auto props = CreateValidProperties();

    const auto geometry = Geometry<Node>{};
    Vector     dummy_vector;
    law.InitializeMaterial(props, geometry, dummy_vector);

    KRATOS_EXPECT_NEAR(CalculateMomentForCurvature(law, props, 0.005), 40.0, Defaults::absolute_tolerance); // elastic
    KRATOS_EXPECT_NEAR(CalculateMomentForCurvature(law, props, 0.02), 80.0, Defaults::absolute_tolerance); // plastic 1
    KRATOS_EXPECT_NEAR(CalculateMomentForCurvature(law, props, 0.04), 104.0, Defaults::absolute_tolerance); // elasto-plastic
    KRATOS_EXPECT_NEAR(CalculateMomentForCurvature(law, props, 0.07), 128.0, Defaults::absolute_tolerance); // second plateau
    KRATOS_EXPECT_NEAR(CalculateMomentForCurvature(law, props, -0.04), -104.0, Defaults::absolute_tolerance); // sign symmetry
}

KRATOS_TEST_CASE_IN_SUITE(PiecewiseLinearMomentCapacityConstitutiveLaw_MomentCapacityLawReturnsExpectedTangentModulus,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto       law   = PiecewiseLinearMomentCapacityConstitutiveLaw{};
    const auto props = CreateValidProperties();

    const auto geometry = Geometry<Node>{};
    Vector     dummy_vector;
    law.InitializeMaterial(props, geometry, dummy_vector);

    KRATOS_EXPECT_NEAR(CalculateTangentForCurvature(law, props, 0.005), 8000.0, Defaults::absolute_tolerance);
    KRATOS_EXPECT_NEAR(CalculateTangentForCurvature(law, props, 0.02), 0.0, Defaults::absolute_tolerance);
    KRATOS_EXPECT_NEAR(CalculateTangentForCurvature(law, props, 0.04), 2400.0, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(PiecewiseLinearMomentCapacityConstitutiveLaw_UnloadReloadWindow_ElasticResponseInsideWindow,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law        = PiecewiseLinearMomentCapacityConstitutiveLaw{};
    auto properties = CreateValidProperties();
    properties.SetValue(UNRELOAD_MODULUS, 1000.0);

    const auto geometry = Geometry<Node>{};
    Vector     dummy_vector;
    law.InitializeMaterial(properties, geometry, dummy_vector);

    // First, load to a curvature beyond initial zero accumulated curvature
    FinalizeForCurvature(law, properties, 0.04);

    // After finalize the law should have a non-zero accumulated curvature and an elastic window
    // Choose a curvature close to the center so it falls inside the elastic window
    const auto test_kappa = -0.06; // selected to be within the small window created

    // Tangent should equal the UNRELOAD_MODULUS inside the window
    KRATOS_EXPECT_NEAR(CalculateTangentForCurvature(law, properties, test_kappa), 1000.0,
                       Defaults::absolute_tolerance);

    // Verify local derivative (finite difference) equals the unload modulus -> elastic linearity
    const auto eps = 1e-6;
    const auto s1  = CalculateMomentForCurvature(law, properties, test_kappa);
    const auto s2  = CalculateMomentForCurvature(law, properties, test_kappa + eps);
    const auto fd  = (s2 - s1) / eps;
    KRATOS_EXPECT_NEAR(fd, 1000.0, Defaults::relative_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(PiecewiseLinearMomentCapacityConstitutiveLaw_SequenceLoadUnloadReload,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law        = PiecewiseLinearMomentCapacityConstitutiveLaw{};
    auto properties = CreateValidProperties();
    properties.SetValue(UNRELOAD_MODULUS, 8000.0);
    properties.SetValue(KAPPA_PIECEWISE_LINEAR_LAW, UblasUtilities::CreateVector({0.01, 0.03, 0.05}));
    properties.SetValue(MOMENTUM_PIECEWISE_LINEAR_LAW, UblasUtilities::CreateVector({80.0, 100.0, 128.0}));

    const auto geometry = Geometry<Node>{};
    Vector     dummy_vector;
    law.InitializeMaterial(properties, geometry, dummy_vector);

    const auto curvatures = UblasUtilities::CreateVector(
        {0.0, 0.01, 0.02, 0.03, 0.025, 0.015, 0.035, 0.01, 0.00825, 0.0, -0.01, 0.02, 0.022, 0.032});

    const auto expected_moment = UblasUtilities::CreateVector(
        {0.0, 80.0, 90.0, 100.0, 60.0, -20.0, 107.0, -93.0, -107.0, -118.55, -128.0, 112.0, 128.0, 128.0});
    const auto expected_tangent_modulus =
        UblasUtilities::CreateVector({8000.0, 8000.0, 1000.0, 1000.0, 8000.0, 8000.0, 1400.0,
                                      8000.0, 1400.0, 1400.0, 0.0, 8000.0, 0.0, 0.0});
    for (auto i = 0; i < curvatures.size(); i++) {
        const auto moment = CalculateMomentForCurvature(law, properties, curvatures(i));
        KRATOS_EXPECT_NEAR(moment, expected_moment(i), Defaults::relative_tolerance);
        const auto tangent_modulus = CalculateTangentForCurvature(law, properties, curvatures(i));
        KRATOS_EXPECT_NEAR(tangent_modulus, expected_tangent_modulus(i), Defaults::relative_tolerance);
        FinalizeForCurvature(law, properties, curvatures(i));
    }
}

} // namespace Kratos::Testing
