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
    properties.SetValue(STRAINS_OF_PIECEWISE_LINEAR_LAW, UblasUtilities::CreateVector({0.01, 0.03, 0.05}));
    properties.SetValue(STRESSES_OF_PIECEWISE_LINEAR_LAW, UblasUtilities::CreateVector({100.0, 160.0}));
    properties.SetValue(MAX_AXIAL_LOAD_OF_CONSTRUCTION_ELEMENT, 10.0);
    properties.SetValue(MOMENT_CAPACITY_REDUCTION_FACTOR, 0.02);
    properties.SetValue(MINIMUM_MOMENT_CAPACITY_FACTOR, 0.30);
    return properties;
}

double CalculateMomentForCurvature(PiecewiseLinearMomentCapacityConstitutiveLaw& rLaw,
                                   const Properties&                             rProperties,
                                   double                                        Curvature)
{
    auto parameters    = ConstitutiveLaw::Parameters{};
    auto strain_vector = UblasUtilities::CreateVector({Curvature});
    parameters.SetStrainVector(strain_vector);

    auto stress_vector = Vector(1, 0.0);
    parameters.SetStressVector(stress_vector);

    auto constitutive_matrix = Matrix(1, 1, 0.0);
    parameters.SetConstitutiveMatrix(constitutive_matrix);

    parameters.SetMaterialProperties(rProperties);
    parameters.Set(ConstitutiveLaw::COMPUTE_STRESS);

    rLaw.CalculateMaterialResponsePK2(parameters);
    return parameters.GetStressVector()[0];
}

double CalculateTangentForCurvature(PiecewiseLinearMomentCapacityConstitutiveLaw& rLaw,
                                    const Properties&                             rProperties,
                                    double                                        Curvature)
{
    auto parameters    = ConstitutiveLaw::Parameters{};
    auto strain_vector = UblasUtilities::CreateVector({Curvature});
    parameters.SetStrainVector(strain_vector);
    parameters.SetMaterialProperties(rProperties);

    double tangent{};
    rLaw.CalculateValue(parameters, TANGENT_MODULUS, tangent);
    return tangent;
}

} // namespace

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(CheckOfMomentCapacityLawThrowsWhenPropertiesDoesNotHaveMaximumAxialLoad,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto law          = PiecewiseLinearMomentCapacityConstitutiveLaw{};
    auto       properties   = CreateValidProperties();
    const auto geometry     = Geometry<Node>{};
    const auto process_info = ProcessInfo{};

    properties.Erase(MAX_AXIAL_LOAD_OF_CONSTRUCTION_ELEMENT);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto rv = law.Check(properties, geometry, process_info),
        "Error: MAX_AXIAL_LOAD_OF_CONSTRUCTION_ELEMENT does not exist")
}

KRATOS_TEST_CASE_IN_SUITE(CheckOfMomentCapacityLawThrowsWhenCurvatureVectorSizeIsIncorrect,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto law          = PiecewiseLinearMomentCapacityConstitutiveLaw{};
    auto       properties   = CreateValidProperties();
    const auto geometry     = Geometry<Node>{};
    const auto process_info = ProcessInfo{};

    properties.SetValue(STRAINS_OF_PIECEWISE_LINEAR_LAW, UblasUtilities::CreateVector({0.01, 0.03}));

    KRATOS_EXPECT_EXCEPTION_IS_THROWN([[maybe_unused]] const auto rv = law.Check(properties, geometry, process_info), "Error: STRAINS_OF_PIECEWISE_LINEAR_LAW requires exactly 3 values: [kappa_elastic_end, kappa_plastic1_end, kappa_elastoplastic_end]")
}

KRATOS_TEST_CASE_IN_SUITE(MomentCapacityLawReturnsExpectedMomentForAllRegimes, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto       law   = PiecewiseLinearMomentCapacityConstitutiveLaw{};
    const auto props = CreateValidProperties();

    KRATOS_EXPECT_NEAR(CalculateMomentForCurvature(law, props, 0.005), 40.0, Defaults::absolute_tolerance); // elastic
    KRATOS_EXPECT_NEAR(CalculateMomentForCurvature(law, props, 0.02), 80.0, Defaults::absolute_tolerance); // plastic 1
    KRATOS_EXPECT_NEAR(CalculateMomentForCurvature(law, props, 0.04), 104.0, Defaults::absolute_tolerance); // elasto-plastic
    KRATOS_EXPECT_NEAR(CalculateMomentForCurvature(law, props, 0.07), 128.0, Defaults::absolute_tolerance); // second plateau
    KRATOS_EXPECT_NEAR(CalculateMomentForCurvature(law, props, -0.04), -104.0, Defaults::absolute_tolerance); // sign symmetry
}

KRATOS_TEST_CASE_IN_SUITE(MomentCapacityLawReturnsExpectedTangentModulus, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto       law   = PiecewiseLinearMomentCapacityConstitutiveLaw{};
    const auto props = CreateValidProperties();

    KRATOS_EXPECT_NEAR(CalculateTangentForCurvature(law, props, 0.005), 8000.0, Defaults::absolute_tolerance);
    KRATOS_EXPECT_NEAR(CalculateTangentForCurvature(law, props, 0.02), 0.0, Defaults::absolute_tolerance);
    KRATOS_EXPECT_NEAR(CalculateTangentForCurvature(law, props, 0.04), 2400.0, Defaults::absolute_tolerance);
}

} // namespace Kratos::Testing
