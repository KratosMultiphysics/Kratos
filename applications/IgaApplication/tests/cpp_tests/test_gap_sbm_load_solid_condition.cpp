//  License: BSD License
//           Kratos default license: kratos/license.txt

#include "containers/model.h"
#include "testing/testing.h"
#include "custom_conditions/gap_sbm_load_solid_condition.h"
#include "test_creation_utility.h"

namespace Kratos::Testing
{
namespace
{
class TestGapSbmLoadSolidCondition : public GapSbmLoadSolidCondition
{
public:
    using GapSbmLoadSolidCondition::GapSbmLoadSolidCondition;

    void InitializeGeometry()
    {
        InitializeMemberVariables();
        InitializeSbmMemberVariables();
    }
};
}

KRATOS_TEST_CASE_IN_SUITE(GapSbmLoadSolidConditionPrescribedForce, KratosIgaFastSuite)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("Test");
    auto p_geometry = TestCreationUtility::GetQuadraturePointGeometryOnCurve(
        r_model_part, 3, IntegrationPoint<3>(0.5, 0.0, 0.0, 2.5));
    auto p_properties = r_model_part.CreateNewProperties(0);
    p_properties->SetValue(THICKNESS, 2.0);
    p_properties->SetValue(YOUNG_MODULUS, 1000.0);
    p_properties->SetValue(POISSON_RATIO, 0.3);
    TestGapSbmLoadSolidCondition condition(1, p_geometry, p_properties);
    condition.SetValue(NEIGHBOUR_GEOMETRIES, std::vector<Geometry<Node>::Pointer>{p_geometry});
    condition.InitializeGeometry();

    array_1d<double, 3> force = ZeroVector(3);
    force[0] = 3.0;
    force[1] = -4.0;
    condition.SetValue(FORCE, force);
    Vector rhs;
    condition.CalculateRightHandSide(rhs, r_model_part.GetProcessInfo());
    double total_x = 0.0;
    double total_y = 0.0;
    for (std::size_t i = 0; i < rhs.size() / 2; ++i) {
        total_x += rhs[2 * i];
        total_y += rhs[2 * i + 1];
    }
    KRATOS_EXPECT_NEAR(total_x, 15.0, 1e-12);
    KRATOS_EXPECT_NEAR(total_y, -20.0, 1e-12);

    condition.SetValue(FORCE, ZeroVector(3));
    condition.CalculateRightHandSide(rhs, r_model_part.GetProcessInfo());
    KRATOS_EXPECT_NEAR(norm_2(rhs), 0.0, 1e-12);
}
} // namespace Kratos::Testing
