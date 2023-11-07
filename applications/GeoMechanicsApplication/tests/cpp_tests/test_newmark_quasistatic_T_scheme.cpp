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

#include "custom_strategies/schemes/newmark_quasistatic_T_scheme.hpp"
#include "spaces/ublas_space.h"
#include "testing/testing.h"

using namespace Kratos;
using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
using LocalSpaceType = UblasSpace<double, Matrix, Vector>;

namespace Kratos::Testing {

NewmarkQuasistaticTScheme<SparseSpaceType, LocalSpaceType> CreateValidScheme()
{
    NewmarkQuasistaticTScheme<SparseSpaceType, LocalSpaceType> result(0.75);
    return result;
}

ModelPart& CreateValidModelPart(Model& rModel)
{
    auto& result = rModel.CreateModelPart("dummy", 2);
    result.AddNodalSolutionStepVariable(TEMPERATURE);
    result.AddNodalSolutionStepVariable(DT_TEMPERATURE);
    auto p_node = result.CreateNewNode(0, 0.0, 0.0, 0.0);
    p_node->AddDof(TEMPERATURE);

    return result;
}

class SpyElement : public Element
{
public:
    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override
    {
        mFinalized = true;
    }

    bool IsFinalized () const
    {
        return mFinalized;
    }

private:
    bool mFinalized = false;
};

class SpyCondition : public Condition
{
public:
    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override
    {
        mFinalized = true;
    }

    bool IsFinalized () const
    {
        return mFinalized;
    }

private:
    bool mFinalized = false;
};


KRATOS_TEST_CASE_IN_SUITE(CheckBackwardEulerQuasistaticTScheme_WithAllNecessaryParts_Returns0,
                          KratosGeoMechanicsFastSuite)
{
    Model model;
    auto scheme = CreateValidScheme();
    const auto& model_part = CreateValidModelPart(model);

    KRATOS_EXPECT_EQ(scheme.Check(model_part), 0);
}

KRATOS_TEST_CASE_IN_SUITE(ForInvalidTheta_CheckBackwardEulerQuasistaticTScheme_Throws,
                          KratosGeoMechanicsFastSuite)
{
    const int invalid_theta = -2;
    NewmarkQuasistaticTScheme<SparseSpaceType, LocalSpaceType> scheme(invalid_theta);

    Model model;
    auto& model_part = model.CreateModelPart("dummy", 2);
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    model_part.AddNodalSolutionStepVariable(DT_TEMPERATURE);
    auto p_node = model_part.CreateNewNode(0, 0.0, 0.0, 0.0);
    p_node->AddDof(TEMPERATURE);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(scheme.Check(model_part),
                                      "Some of the scheme variables: beta, "
                                      "gamma or theta has an invalid value ")
}

KRATOS_TEST_CASE_IN_SUITE(ForInvalidBufferSize_CheckBackwardEulerQuasistaticTScheme_Throws,
                          KratosGeoMechanicsFastSuite)
{
    NewmarkQuasistaticTScheme<SparseSpaceType, LocalSpaceType> scheme(0.75);

    Model model;
    const int invalid_buffer_size = 1;
    auto& model_part = model.CreateModelPart("dummy", invalid_buffer_size);
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    model_part.AddNodalSolutionStepVariable(DT_TEMPERATURE);
    auto p_node = model_part.CreateNewNode(0, 0.0, 0.0, 0.0);
    p_node->AddDof(TEMPERATURE);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        scheme.Check(model_part),
        "insufficient buffer size. Buffer size should be greater or equal to "
        "2. Current size is ")
}

KRATOS_TEST_CASE_IN_SUITE(ForMissingNodalDof_CheckBackwardEulerQuasistaticTScheme_Throws,
                          KratosGeoMechanicsFastSuite)
{
    NewmarkQuasistaticTScheme<SparseSpaceType, LocalSpaceType> scheme(0.75);

    Model model;
    auto& model_part = model.CreateModelPart("dummy", 2);
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    model_part.AddNodalSolutionStepVariable(DT_TEMPERATURE);
    auto p_node = model_part.CreateNewNode(0, 0.0, 0.0, 0.0);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(scheme.Check(model_part),
                                      "missing TEMPERATURE dof on node ")
}

KRATOS_TEST_CASE_IN_SUITE(ForMissingDtTemperatureSolutionStepVariable_CheckBackwardEulerQuasistaticTScheme_Throws,
                          KratosGeoMechanicsFastSuite)
{
    NewmarkQuasistaticTScheme<SparseSpaceType, LocalSpaceType> scheme(0.75);

    Model model;
    auto& model_part = model.CreateModelPart("dummy", 2);
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    auto p_node = model_part.CreateNewNode(0, 0.0, 0.0, 0.0);
    p_node->AddDof(TEMPERATURE);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        scheme.Check(model_part),
        "DT_TEMPERATURE variable is not allocated for node 0")
}

KRATOS_TEST_CASE_IN_SUITE(ForMissingTemperatureSolutionStepVariable_CheckBackwardEulerQuasistaticTScheme_Throws,
                          KratosGeoMechanicsFastSuite)
{
    NewmarkQuasistaticTScheme<SparseSpaceType, LocalSpaceType> scheme(0.75);

    Model model;
    auto& model_part = model.CreateModelPart("dummy", 2);
    model_part.AddNodalSolutionStepVariable(DT_TEMPERATURE);
    auto p_node = model_part.CreateNewNode(0, 0.0, 0.0, 0.0);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        scheme.Check(model_part),
        "TEMPERATURE variable is not allocated for node 0")
}

KRATOS_TEST_CASE_IN_SUITE(ThermalSchemeUpdate_SetsDtTemperature, KratosGeoMechanicsFastSuite)
{
    auto scheme = CreateValidScheme();
    Model model;
    ModelPart& model_part = CreateValidModelPart(model);

    const double current_temperature = 10.0;
    const double previous_temperature = 5.0;
    const double previous_dt_temperature = 3.0;
    const double theta = 0.75;
    const double delta_time = 2.0;

    model_part.GetProcessInfo()[DELTA_TIME] = delta_time;
    Node& node = model_part.Nodes()[0];
    node.FastGetSolutionStepValue(TEMPERATURE) = current_temperature;
    node.FastGetSolutionStepValue(TEMPERATURE, 1) = previous_temperature;
    node.FastGetSolutionStepValue(DT_TEMPERATURE, 1) = previous_dt_temperature;

    double expected_value = 1.0 / (theta * delta_time) *
                            (current_temperature - previous_temperature -
                             (1.0 - theta) * delta_time * previous_dt_temperature);

    ModelPart::DofsArrayType dof_set;
    CompressedMatrix A;
    Vector Dx;
    Vector b;

    scheme.Initialize(model_part);
    scheme.Predict(model_part, dof_set, A, Dx, b);

    KRATOS_EXPECT_DOUBLE_EQ(node.FastGetSolutionStepValue(DT_TEMPERATURE), expected_value);
}

KRATOS_TEST_CASE_IN_SUITE(InitializeScheme_SetsTimeFactors, KratosGeoMechanicsFastSuite)
{
    auto scheme = CreateValidScheme();
    Model model;
    ModelPart& model_part = CreateValidModelPart(model);
    const double delta_time = 3.0;
    model_part.GetProcessInfo()[DELTA_TIME] = delta_time;

    scheme.Initialize(model_part);

    KRATOS_EXPECT_TRUE(scheme.SchemeIsInitialized())
    const double theta = 0.75;
    KRATOS_EXPECT_DOUBLE_EQ(model_part.GetProcessInfo()[DT_TEMPERATURE_COEFFICIENT],
                            1.0 / (theta * delta_time));
}

KRATOS_TEST_CASE_IN_SUITE(FinalizeSolutionStepActiveEntities_FinalizesOnlyActiveElements,
                          KratosGeoMechanicsFastSuite)
{
    auto scheme = CreateValidScheme();
    Model model;
    ModelPart& model_part = CreateValidModelPart(model);

    auto active_element = Kratos::make_intrusive<SpyElement>();
    active_element->Set(ACTIVE, true);
    active_element->SetId(0);

    auto inactive_element = Kratos::make_intrusive<SpyElement>();
    inactive_element->Set(ACTIVE, false);
    inactive_element->SetId(1);

    model_part.AddElement(active_element);
    model_part.AddElement(inactive_element);

    CompressedMatrix A;
    Vector Dx;
    Vector b;
    scheme.FinalizeSolutionStep(model_part, A, Dx, b);

    KRATOS_EXPECT_TRUE(active_element->IsFinalized())
    KRATOS_EXPECT_FALSE(inactive_element->IsFinalized())
}

KRATOS_TEST_CASE_IN_SUITE(FinalizeSolutionStep_FinalizesOnlyActiveConditions, KratosGeoMechanicsFastSuite)
{
    auto scheme = CreateValidScheme();
    Model model;
    ModelPart& model_part = CreateValidModelPart(model);

    auto active_condition = Kratos::make_intrusive<SpyCondition>();
    active_condition->SetId(0);
    active_condition->Set(ACTIVE, true);

    auto inactive_condition = Kratos::make_intrusive<SpyCondition>();
    inactive_condition->SetId(1);
    inactive_condition->Set(ACTIVE, false);

    model_part.AddCondition(active_condition);
    model_part.AddCondition(inactive_condition);

    CompressedMatrix A;
    Vector Dx;
    Vector b;
    scheme.FinalizeSolutionStep(model_part, A, Dx, b);

    KRATOS_EXPECT_TRUE(active_condition->IsFinalized())
    KRATOS_EXPECT_FALSE(inactive_condition->IsFinalized())
}

} // namespace Kratos::Test