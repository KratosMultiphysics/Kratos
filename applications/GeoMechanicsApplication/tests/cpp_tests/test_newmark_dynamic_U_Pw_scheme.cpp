// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//

#include "custom_strategies/schemes/newmark_dynamic_U_Pw_scheme.hpp"
#include "spaces/ublas_space.h"
#include "testing/testing.h"

namespace Kratos::Testing
{

using namespace Kratos;
using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
using LocalSpaceType = UblasSpace<double, Matrix, Vector>;

class NewmarkDynamicUPwSchemeTester
{
public:
    Model mModel;
    NewmarkDynamicUPwScheme<SparseSpaceType, LocalSpaceType> mScheme =
        CreateValidScheme();

    explicit NewmarkDynamicUPwSchemeTester(const bool ThreeDCase = true)
    {
        CreateValidModelPart(ThreeDCase);
    }

    NewmarkDynamicUPwScheme<SparseSpaceType, LocalSpaceType> CreateValidScheme() const
    {
        return NewmarkDynamicUPwScheme<SparseSpaceType, LocalSpaceType>(0.25, 0.5, 0.75);
    }

    void CreateValidModelPart(const bool Add3DDofs)
    {
        auto& result = mModel.CreateModelPart("dummy", 2);
        result.AddNodalSolutionStepVariable(DISPLACEMENT);
        result.AddNodalSolutionStepVariable(VELOCITY);
        result.AddNodalSolutionStepVariable(ACCELERATION);
        result.AddNodalSolutionStepVariable(WATER_PRESSURE);
        result.AddNodalSolutionStepVariable(DT_WATER_PRESSURE);
        result.AddNodalSolutionStepVariable(ROTATION);
        result.AddNodalSolutionStepVariable(ANGULAR_VELOCITY);
        result.AddNodalSolutionStepVariable(ANGULAR_ACCELERATION);

        auto p_node = result.CreateNewNode(0, 0.0, 0.0, 0.0);
        p_node->AddDof(DISPLACEMENT_X);
        p_node->AddDof(DISPLACEMENT_Y);
        p_node->AddDof(WATER_PRESSURE);
        p_node->AddDof(ROTATION_X);
        p_node->AddDof(ROTATION_Y);
        result.GetProcessInfo()[DELTA_TIME] = 4.0;

        p_node->FastGetSolutionStepValue(VELOCITY, 0) =
            Kratos::array_1d<double, 3>{1.0, 2.0, 3.0};
        p_node->FastGetSolutionStepValue(ACCELERATION, 0) =
            Kratos::array_1d<double, 3>{4.0, 5.0, 6.0};
        p_node->FastGetSolutionStepValue(DISPLACEMENT, 0) =
            Kratos::array_1d<double, 3>{0.0, 0.0, 0.0};

        p_node->FastGetSolutionStepValue(VELOCITY, 1) =
            Kratos::array_1d<double, 3>{10.0, 11.0, 12.0};
        p_node->FastGetSolutionStepValue(ACCELERATION, 1) =
            Kratos::array_1d<double, 3>{13.0, 14.0, 15.0};
        p_node->FastGetSolutionStepValue(DISPLACEMENT, 1) =
            Kratos::array_1d<double, 3>{16.0, 17.0, 18.0};

        p_node->FastGetSolutionStepValue(WATER_PRESSURE, 1) = 1.0;
        p_node->FastGetSolutionStepValue(WATER_PRESSURE, 0) = 2.0;

        if (Add3DDofs)
        {
            p_node->AddDof(DISPLACEMENT_Z);
            p_node->AddDof(ROTATION_Z);
        }
    }

    ModelPart& GetModelPart() { return mModel.GetModelPart("dummy"); }
};

KRATOS_TEST_CASE_IN_SUITE(NewmarkDynamicUPwSchemePredictWithFixedAccelerations_UpdatesVariablesDerivatives,
                          KratosGeoMechanicsFastSuite)
{
    NewmarkDynamicUPwSchemeTester tester;

    tester.mScheme.Initialize(tester.GetModelPart()); // This is needed to set the time factors
    ModelPart::DofsArrayType dof_set;
    CompressedMatrix A;
    Vector Dx;
    Vector b;

    tester.GetModelPart().GetNode(0).Fix(ACCELERATION_X);
    tester.GetModelPart().GetNode(0).Fix(ACCELERATION_Y);
    tester.GetModelPart().GetNode(0).Fix(ACCELERATION_Z);

    tester.mScheme.Predict(tester.GetModelPart(), dof_set, A, Dx, b);

    // These expected numbers result from the calculations in UpdateVariablesDerivatives
    const auto expected_displacement = Kratos::array_1d<double, 3>{124, 137, 150};
    const auto actual_displacement =
        tester.GetModelPart().Nodes()[0].FastGetSolutionStepValue(DISPLACEMENT, 0);
    KRATOS_EXPECT_VECTOR_NEAR(actual_displacement, expected_displacement, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(NewmarkDynamicUPwSchemePredictWithFixedVelocities_UpdatesVariablesDerivatives,
                          KratosGeoMechanicsFastSuite)
{
    NewmarkDynamicUPwSchemeTester tester;

    tester.mScheme.Initialize(tester.GetModelPart()); // This is needed to set the time factors
    ModelPart::DofsArrayType dof_set;
    CompressedMatrix A;
    Vector Dx;
    Vector b;

    tester.GetModelPart().GetNode(0).Fix(VELOCITY_X);
    tester.GetModelPart().GetNode(0).Fix(VELOCITY_Y);
    tester.GetModelPart().GetNode(0).Fix(VELOCITY_Z);

    tester.mScheme.Predict(tester.GetModelPart(), dof_set, A, Dx, b);

    // These expected numbers result from the calculations in UpdateVariablesDerivatives
    const auto expected_displacement = Kratos::array_1d<double, 3>{38, 43, 48};
    const auto actual_displacement =
        tester.GetModelPart().Nodes()[0].FastGetSolutionStepValue(DISPLACEMENT, 0);
    KRATOS_EXPECT_VECTOR_NEAR(actual_displacement, expected_displacement, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(NewmarkDynamicUPwSchemePredictWithNoFixedVariables_UpdatesVariablesDerivatives,
                          KratosGeoMechanicsFastSuite)
{
    NewmarkDynamicUPwSchemeTester tester;

    tester.mScheme.Initialize(tester.GetModelPart()); // This is needed to set the time factors
    ModelPart::DofsArrayType dof_set;
    CompressedMatrix A;
    Vector Dx;
    Vector b;

    tester.mScheme.Predict(tester.GetModelPart(), dof_set, A, Dx, b);

    // These expected numbers result from the calculations in UpdateVariablesDerivatives
    const auto expected_displacement = Kratos::array_1d<double, 3>{160, 173, 186};
    const auto actual_displacement =
        tester.GetModelPart().Nodes()[0].FastGetSolutionStepValue(DISPLACEMENT, 0);
    KRATOS_EXPECT_VECTOR_NEAR(actual_displacement, expected_displacement, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(NewmarkDynamicUPwSchemePredictFixedDisplacements_DoesNotUpdateDisplacements,
                          KratosGeoMechanicsFastSuite)
{
    NewmarkDynamicUPwSchemeTester tester;

    tester.mScheme.Initialize(tester.GetModelPart()); // This is needed to set the time factors
    ModelPart::DofsArrayType dof_set;
    CompressedMatrix A;
    Vector Dx;
    Vector b;

    tester.GetModelPart().GetNode(0).Fix(DISPLACEMENT_X);
    tester.GetModelPart().GetNode(0).Fix(DISPLACEMENT_Y);
    tester.GetModelPart().GetNode(0).Fix(DISPLACEMENT_Z);

    tester.mScheme.Predict(tester.GetModelPart(), dof_set, A, Dx, b);

    // Nothing should happen to the displacement, so it stays (0, 0, 0)
    const auto expected_displacement = Kratos::array_1d<double, 3>{0.0, 0.0, 0.0};
    const auto actual_displacement =
        tester.GetModelPart().Nodes()[0].FastGetSolutionStepValue(DISPLACEMENT, 0);
    KRATOS_EXPECT_VECTOR_NEAR(actual_displacement, expected_displacement, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(NewmarkDynamicUPwSchemePredictWithout3DDofs_DoesNotUpdateZDisplacement,
                          KratosGeoMechanicsFastSuite)
{
    const bool three_d_case = false;
    NewmarkDynamicUPwSchemeTester tester(three_d_case);

    tester.mScheme.Initialize(tester.GetModelPart()); // This is needed to set the time factors
    ModelPart::DofsArrayType dof_set;
    CompressedMatrix A;
    Vector Dx;
    Vector b;

    tester.mScheme.Predict(tester.GetModelPart(), dof_set, A, Dx, b);

    // The Z Displacement should stay 0.0, since it is not a dof
    const auto expected_displacement = Kratos::array_1d<double, 3>{160, 173, 0.0};
    const auto actual_displacement =
        tester.GetModelPart().Nodes()[0].FastGetSolutionStepValue(DISPLACEMENT, 0);
    KRATOS_EXPECT_VECTOR_NEAR(actual_displacement, expected_displacement, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(NewmarkDynamicUPwSchemePredict_UpdatesVariablesDerivatives,
                          KratosGeoMechanicsFastSuite)
{
    NewmarkDynamicUPwSchemeTester tester;

    tester.mScheme.Initialize(tester.GetModelPart()); // This is needed to set the time factors
    ModelPart::DofsArrayType dof_set;
    CompressedMatrix A;
    Vector Dx;
    Vector b;

    tester.mScheme.Predict(tester.GetModelPart(), dof_set, A, Dx, b);

    // These expected numbers result from the calculations in UpdateVariablesDerivatives
    const auto expected_acceleration = Kratos::array_1d<double, 3>{13, 14, 15};
    const auto expected_velocity = Kratos::array_1d<double, 3>{62, 67, 72};
    constexpr auto expected_dt_water_pressure = 1.0 / 3.0;

    const auto actual_acceleration =
        tester.GetModelPart().Nodes()[0].FastGetSolutionStepValue(ACCELERATION, 0);
    KRATOS_EXPECT_VECTOR_NEAR(actual_acceleration, expected_acceleration, 1e-6)

    const auto actual_velocity =
        tester.GetModelPart().Nodes()[0].FastGetSolutionStepValue(VELOCITY, 0);
    KRATOS_EXPECT_VECTOR_NEAR(actual_velocity, expected_velocity, 1e-6)

    KRATOS_EXPECT_DOUBLE_EQ(
        tester.GetModelPart().Nodes()[0].FastGetSolutionStepValue(DT_WATER_PRESSURE, 0),
        expected_dt_water_pressure);
}

} // namespace Kratos::Testing