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

#include "custom_constitutive/piecewise_linear_moment_capacity_plane_strain_constitutive_law.h"
#include "custom_utilities/check_utilities.hpp"
#include "custom_utilities/registration_utilities.hpp"
#include "custom_utilities/ublas_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "includes/stream_serializer.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

using namespace Kratos;

namespace
{

Properties CreateValidProperties()
{
    auto properties = Properties{};
    properties.SetValue(GEO_KAPPA_PIECEWISE_LINEAR_LAW, UblasUtilities::CreateVector({0.01, 0.03, 0.05}));
    properties.SetValue(GEO_MOMENT_PIECEWISE_LINEAR_LAW, UblasUtilities::CreateVector({80.0, 80.0, 128.0}));
    properties.SetValue(YOUNG_MODULUS, 80.0);
    properties.SetValue(POISSON_RATIO, 0.2);
    properties.SetValue(THICKNESS, 1.0);
    properties.SetValue(THICKNESS_EFFECTIVE_Y, 1.0);
    return properties;
}

double CalculateMomentForCurvature(PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw& rLaw,
                                   const Properties& rProperties,
                                   double            Curvature)
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

double CalculateTangentForCurvature(PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw& rLaw,
                                    const Properties& rProperties,
                                    double            Curvature)
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

void FinalizeForCurvature(PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw& rLaw,
                          const Properties&                                        rProperties,
                          double                                                   Curvature)
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
    // Arrange
    const auto law          = PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw{};
    auto       properties   = CreateValidProperties();
    const auto geometry     = Geometry<Node>{};
    const auto process_info = ProcessInfo{};

    properties.SetValue(GEO_KAPPA_PIECEWISE_LINEAR_LAW, UblasUtilities::CreateVector({0.01, 0.03}));

    // Act & Assert
    KRATOS_EXPECT_EXCEPTION_IS_THROWN((void)law.Check(properties, geometry, process_info),
                                      "The number of entries in GEO_KAPPA_PIECEWISE_LINEAR_LAW (2) "
                                      "does not match GEO_MOMENT_PIECEWISE_LINEAR_LAW (3)")
}

KRATOS_TEST_CASE_IN_SUITE(PiecewiseLinearMomentCapacityConstitutiveLaw_CheckThrowsWhenYoungModulusIsMissing,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto law          = PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw{};
    auto       properties   = CreateValidProperties();
    const auto geometry     = Geometry<Node>{};
    const auto process_info = ProcessInfo{};
    properties.Erase(YOUNG_MODULUS);

    // Act & Assert
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        (void)law.Check(properties, geometry, process_info),
        "YOUNG_MODULUS does not exist in the material properties with Id 0.")
}

KRATOS_TEST_CASE_IN_SUITE(PiecewiseLinearMomentCapacityConstitutiveLaw_CheckThrowsWhenThicknessIsMissing,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto law          = PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw{};
    auto       properties   = CreateValidProperties();
    const auto geometry     = Geometry<Node>{};
    const auto process_info = ProcessInfo{};
    properties.Erase(THICKNESS);

    // Act & Assert
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        (void)law.Check(properties, geometry, process_info),
        "THICKNESS does not exist in the material properties with Id 0.")
}

KRATOS_TEST_CASE_IN_SUITE(PiecewiseLinearMomentCapacityConstitutiveLaw_CheckThrowsWhenThicknessEffectiveYIsMissing,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto law          = PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw{};
    auto       properties   = CreateValidProperties();
    const auto geometry     = Geometry<Node>{};
    const auto process_info = ProcessInfo{};
    properties.Erase(THICKNESS_EFFECTIVE_Y);

    // Act & Assert
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        (void)law.Check(properties, geometry, process_info),
        "THICKNESS_EFFECTIVE_Y does not exist in the material properties with Id 0.")
}

KRATOS_TEST_CASE_IN_SUITE(PiecewiseLinearMomentCapacityConstitutiveLaw_CheckThrowsWhenPoissonRatioIsOutOfRange,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto law          = PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw{};
    auto       properties   = CreateValidProperties();
    const auto geometry     = Geometry<Node>{};
    const auto process_info = ProcessInfo{};
    properties.SetValue(POISSON_RATIO, 0.5); // upper bound exclusive

    // Act & Assert
    KRATOS_EXPECT_EXCEPTION_IS_THROWN((void)law.Check(properties, geometry, process_info),
                                      "POISSON_RATIO in the material properties with Id 0 has an "
                                      "invalid value: 0.5 is out of the range (-1, 0.5).")
}

KRATOS_TEST_CASE_IN_SUITE(PiecewiseLinearMomentCapacityConstitutiveLaw_CheckThrowsWhenUnreloadModulusIsZero,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto law          = PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw{};
    auto       properties   = CreateValidProperties();
    const auto geometry     = Geometry<Node>{};
    const auto process_info = ProcessInfo{};
    properties.SetValue(GEO_UNRELOAD_MODULUS, 0.0); // must be strictly positive

    // Act & Assert
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        (void)law.Check(properties, geometry, process_info),
        "GEO_UNRELOAD_MODULUS in the material properties with Id 0 has "
        "an invalid value: 0 is out of the range (0, -).")
}

KRATOS_TEST_CASE_IN_SUITE(PiecewiseLinearMomentCapacityConstitutiveLaw_CheckPassesWithValidProperties,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto law          = PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw{};
    const auto properties   = CreateValidProperties();
    const auto geometry     = Geometry<Node>{};
    const auto process_info = ProcessInfo{};

    // Act & Assert
    EXPECT_NO_THROW((void)law.Check(properties, geometry, process_info));
}

KRATOS_TEST_CASE_IN_SUITE(PiecewiseLinearMomentCapacityConstitutiveLaw_RequiresFinalizeMaterialResponseReturnsTrue,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto law = PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw{};

    // Act
    const auto requires_finalize = law.RequiresFinalizeMaterialResponse();
    // Assert
    KRATOS_EXPECT_TRUE(requires_finalize)
}

KRATOS_TEST_CASE_IN_SUITE(PiecewiseLinearMomentCapacityConstitutiveLaw_CloneCreatesEquivalentLaw,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto law = PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw{};

    // Act
    const auto p_clone = law.Clone();

    // Assert
    KRATOS_EXPECT_NE(p_clone.get(), nullptr);
    KRATOS_EXPECT_TRUE(
        dynamic_cast<PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw*>(p_clone.get()) != nullptr)
    KRATOS_EXPECT_EQ(p_clone->Info(), law.Info());
}

KRATOS_TEST_CASE_IN_SUITE(PiecewiseLinearMomentCapacityConstitutiveLaw_InfoReturnsClassName,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto law = PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw{};

    // Act
    const auto info = law.Info();

    // Assert
    KRATOS_EXPECT_EQ(info, "PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw");
}

KRATOS_TEST_CASE_IN_SUITE(PiecewiseLinearMomentCapacityConstitutiveLaw_GetLawFeaturesReturnsExpectedDimensions,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto law      = PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw{};
    auto features = ConstitutiveLaw::Features{};

    // Act
    law.GetLawFeatures(features);

    // Assert
    KRATOS_EXPECT_EQ(features.mStrainSize, law.GetStrainSize());
    KRATOS_EXPECT_EQ(features.mSpaceDimension,
                     PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw::space_dimenstion);
}

KRATOS_TEST_CASE_IN_SUITE(PiecewiseLinearMomentCapacityConstitutiveLaw_ComputeConstitutiveTensorReturnsExpectedDiagonal,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto       law        = PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw{};
    const auto properties = CreateValidProperties();
    const auto geometry   = Geometry<Node>{};
    Vector     dummy_vector;
    law.InitializeMaterial(properties, geometry, dummy_vector);

    auto parameters     = ConstitutiveLaw::Parameters{};
    auto strain_vector  = UblasUtilities::CreateVector({0.01, 0.02, 0.03});
    auto stress_vector  = Vector(3, 0.0);
    auto constit_matrix = Matrix(1, 1, 0.0);
    parameters.SetStrainVector(strain_vector);
    parameters.SetStressVector(stress_vector);
    parameters.SetConstitutiveMatrix(constit_matrix);
    parameters.SetMaterialProperties(properties);
    parameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    parameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

    // Act
    law.CalculateMaterialResponseCauchy(parameters);

    // Assert
    const auto& r_constitutive_matrix = parameters.GetConstitutiveMatrix();
    KRATOS_EXPECT_EQ(r_constitutive_matrix.size1(), 3);
    KRATOS_EXPECT_EQ(r_constitutive_matrix.size2(), 3);

    const auto E                        = properties[YOUNG_MODULUS];
    const auto A                        = properties[THICKNESS];
    const auto nu                       = properties[POISSON_RATIO];
    const auto one_minus_nu_squared     = 1.0 - nu * nu;
    const auto expected_axial_stiffness = E * A / one_minus_nu_squared;
    const auto expected_shear_stiffness = (E / (2.0 * (1.0 + nu))) * properties[THICKNESS_EFFECTIVE_Y];
    constexpr auto expected_bending_tangent = 0.0;

    auto expected_constitutive_matrix  = Matrix(3, 3, 0.0);
    expected_constitutive_matrix(0, 0) = expected_axial_stiffness;
    expected_constitutive_matrix(1, 1) = expected_bending_tangent;
    expected_constitutive_matrix(2, 2) = expected_shear_stiffness;
    KRATOS_EXPECT_MATRIX_NEAR(r_constitutive_matrix, expected_constitutive_matrix, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(PiecewiseLinearMomentCapacityConstitutiveLaw_SaveLoadPreservesInternalState,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto law        = PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw{};
    auto properties = CreateValidProperties();
    properties.SetValue(GEO_UNRELOAD_MODULUS, 1000.0);
    const auto geometry = Geometry<Node>{};
    Vector     dummy_vector;
    law.InitializeMaterial(properties, geometry, dummy_vector);

    FinalizeForCurvature(law, properties, 0.04);
    const auto check_curvature = -0.06;
    const auto moment_before_serialize = CalculateMomentForCurvature(law, properties, check_curvature);
    const auto tangent_before_serialize = CalculateTangentForCurvature(law, properties, check_curvature);

    const auto scoped_registration = ScopedSerializerRegistration{
        std::make_pair("PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw",
                       PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw{})};
    auto serializer = StreamSerializer{};

    auto p_law = std::unique_ptr<ConstitutiveLaw>{
        std::make_unique<PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw>(law)};

    // Act
    serializer.save("test_tag", p_law);
    auto p_loaded_law = std::unique_ptr<ConstitutiveLaw>{};
    serializer.load("test_tag", p_loaded_law);

    // Assert
    KRATOS_EXPECT_NE(p_loaded_law.get(), nullptr);
    auto p_loaded_piecewise =
        dynamic_cast<PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw*>(p_loaded_law.get());
    KRATOS_EXPECT_NE(p_loaded_piecewise, nullptr);

    const auto moment_after_serialize =
        CalculateMomentForCurvature(*p_loaded_piecewise, properties, check_curvature);
    const auto tangent_after_serialize =
        CalculateTangentForCurvature(*p_loaded_piecewise, properties, check_curvature);

    KRATOS_EXPECT_NEAR(moment_before_serialize, moment_after_serialize, Defaults::absolute_tolerance);
    KRATOS_EXPECT_NEAR(tangent_before_serialize, tangent_after_serialize, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(PiecewiseLinearMomentCapacityConstitutiveLaw_MomentCapacityLawReturnsExpectedMomentForAllRegimes,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto       law        = PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw{};
    const auto properties = CreateValidProperties();

    const auto geometry = Geometry<Node>{};
    Vector     dummy_vector;
    law.InitializeMaterial(properties, geometry, dummy_vector);
    const auto curvatures = UblasUtilities::CreateVector({0.005, 0.02, 0.04, 0.07, -0.04});

    // Act
    auto calculated_moment = Vector(curvatures.size(), 0.0);
    for (auto i = std::size_t{0}; i < curvatures.size(); i++)
        calculated_moment[i] = CalculateMomentForCurvature(law, properties, curvatures[i]);

    // Assert
    const auto expected_moments = UblasUtilities::CreateVector({40.0, 80.0, 104.0, 128.0, -104.0});
    KRATOS_EXPECT_VECTOR_NEAR(calculated_moment, expected_moments, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(PiecewiseLinearMomentCapacityConstitutiveLaw_MomentCapacityLawReturnsExpectedTangentModulus,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto       law        = PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw{};
    const auto properties = CreateValidProperties();

    const auto geometry = Geometry<Node>{};
    Vector     dummy_vector;
    law.InitializeMaterial(properties, geometry, dummy_vector);
    const auto curvatures = UblasUtilities::CreateVector({0.005, 0.02, 0.04});

    // Act
    auto calculated_tangent_modulus = Vector(curvatures.size(), 0.0);
    for (auto i = std::size_t{0}; i < curvatures.size(); i++)
        calculated_tangent_modulus[i] = CalculateTangentForCurvature(law, properties, curvatures[i]);

    // Assert
    const auto expected_tangent_modulus = UblasUtilities::CreateVector({8000.0, 0.0, 2400.0});
    KRATOS_EXPECT_VECTOR_NEAR(calculated_tangent_modulus, expected_tangent_modulus, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(PiecewiseLinearMomentCapacityConstitutiveLaw_UnloadReloadWindow_ElasticResponseInsideWindow,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto law        = PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw{};
    auto properties = CreateValidProperties();
    properties.SetValue(GEO_UNRELOAD_MODULUS, 1000.0);

    const auto geometry = Geometry<Node>{};
    Vector     dummy_vector;
    law.InitializeMaterial(properties, geometry, dummy_vector);

    // First, load to a curvature beyond initial zero accumulated curvature
    FinalizeForCurvature(law, properties, 0.04);

    // After finalize the law should have a non-zero accumulated curvature and an elastic window
    // Choose a curvature close to the center so it falls inside the elastic window
    const auto test_kappa = -0.06; // selected to be within the small window created

    // Act & Assert
    // Tangent should equal the GEO_UNRELOAD_MODULUS inside the window
    const auto calculated_tangent_modulus = CalculateTangentForCurvature(law, properties, test_kappa);
    constexpr auto expected_tangent_modulus = 1000.0;
    KRATOS_EXPECT_NEAR(calculated_tangent_modulus, expected_tangent_modulus, Defaults::absolute_tolerance);

    // Verify local derivative (finite difference) equals the unload modulus -> elastic linearity
    const auto eps      = 1e-6;
    const auto moment_1 = CalculateMomentForCurvature(law, properties, test_kappa);
    const auto moment_2 = CalculateMomentForCurvature(law, properties, test_kappa + eps);
    const auto fd       = (moment_2 - moment_1) / eps;
    KRATOS_EXPECT_NEAR(fd, expected_tangent_modulus, Defaults::relative_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(PiecewiseLinearMomentCapacityConstitutiveLaw_SequenceLoadUnloadReload,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto law        = PiecewiseLinearMomentCapacityPlaneStrainConstitutiveLaw{};
    auto properties = CreateValidProperties();
    properties.SetValue(GEO_UNRELOAD_MODULUS, 8000.0);
    properties.SetValue(GEO_KAPPA_PIECEWISE_LINEAR_LAW, UblasUtilities::CreateVector({0.01, 0.03, 0.05}));
    properties.SetValue(GEO_MOMENT_PIECEWISE_LINEAR_LAW, UblasUtilities::CreateVector({80.0, 100.0, 128.0}));

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

    // Act & Asset
    for (auto i = std::size_t{0}; i < curvatures.size(); i++) {
        const auto moment = CalculateMomentForCurvature(law, properties, curvatures(i));
        KRATOS_EXPECT_NEAR(moment, expected_moment(i), Defaults::relative_tolerance);
        const auto tangent_modulus = CalculateTangentForCurvature(law, properties, curvatures(i));
        KRATOS_EXPECT_NEAR(tangent_modulus, expected_tangent_modulus(i), Defaults::relative_tolerance);
        FinalizeForCurvature(law, properties, curvatures(i));
    }
}

} // namespace Kratos::Testing
