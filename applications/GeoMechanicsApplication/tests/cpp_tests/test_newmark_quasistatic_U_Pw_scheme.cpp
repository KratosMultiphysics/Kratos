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

#include "testing/testing.h"

#include "custom_strategies/schemes/newmark_quasistatic_U_Pw_scheme.hpp"
#include "spaces/ublas_space.h"

namespace {

using namespace Kratos;
using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
using LocalSpaceType = UblasSpace<double, Matrix, Vector>;

class SpyElement : public Element {
public:
    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override
    {
        mSolutionStepFinalized = true;
    }
    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override
    {
        mSolutionStepInitialized = true;
    }
    void InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo) override
    {
        mNonLinIterationInitialized = true;
    }
    void FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo) override
    {
        mNonLinIterationFinalized = true;
    }

    bool IsSolutionStepFinalized() const
    {
        return mSolutionStepFinalized;
    }

    bool IsSolutionStepInitialized() const
    {
        return mSolutionStepInitialized;
    }

    bool IsNonLinIterationInitialized() const
    {
        return mNonLinIterationInitialized;
    }

    bool IsNonLinIterationFinalized() const
    {
        return mNonLinIterationFinalized;
    }

private:
    bool mSolutionStepInitialized = false;
    bool mSolutionStepFinalized = false;
    bool mNonLinIterationInitialized = false;
    bool mNonLinIterationFinalized = false;
};

class SpyCondition : public Condition {
public:
    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override
    {
        mSolutionStepFinalized = true;
    }
    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override
    {
        mSolutionStepInitialized = true;
    }
    void InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo) override
    {
        mNonLinIterationInitialized = true;
    }
    void FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo) override
    {
        mNonLinIterationFinalized = true;
    }

    bool IsSolutionStepFinalized() const
    {
        return mSolutionStepFinalized;
    }

    bool IsSolutionStepInitialized() const
    {
        return mSolutionStepInitialized;
    }

    bool IsNonLinIterationInitialized() const
    {
        return mNonLinIterationInitialized;
    }

    bool IsNonLinIterationFinalized() const
    {
        return mNonLinIterationFinalized;
    }

private:
    bool mSolutionStepInitialized = false;
    bool mSolutionStepFinalized = false;
    bool mNonLinIterationInitialized = false;
    bool mNonLinIterationFinalized = false;
};

} // namespace

namespace Kratos::Testing {

class NewmarkQuasistaticUPwSchemeTester {
public:
    Model mModel;
    ModelPart* mrModelPart = &CreateValidModelPart();
    NewmarkQuasistaticUPwScheme<SparseSpaceType, LocalSpaceType> mScheme =
        CreateValidScheme();

    void Setup()
    {
        mModel.Reset();
        mrModelPart = &CreateValidModelPart();
        mScheme = CreateValidScheme();
    }

    NewmarkQuasistaticUPwScheme<SparseSpaceType, LocalSpaceType> CreateValidScheme() const
    {
        NewmarkQuasistaticUPwScheme<SparseSpaceType, LocalSpaceType> result(0.25, 0.5, 0.75);
        return result;
    }

    ModelPart& CreateValidModelPart()
    {
        auto& result = mModel.CreateModelPart("dummy", 2);
        result.AddNodalSolutionStepVariable(TEMPERATURE);
        result.AddNodalSolutionStepVariable(VELOCITY);
        result.AddNodalSolutionStepVariable(ACCELERATION);
        result.AddNodalSolutionStepVariable(DT_WATER_PRESSURE);
        result.AddNodalSolutionStepVariable(DISPLACEMENT);
        result.AddNodalSolutionStepVariable(WATER_PRESSURE);

        auto p_node = result.CreateNewNode(0, 0.0, 0.0, 0.0);
        p_node->AddDof(DISPLACEMENT_X);
        p_node->AddDof(DISPLACEMENT_Y);
        p_node->AddDof(DISPLACEMENT_Z);
        p_node->AddDof(WATER_PRESSURE);
        result.GetProcessInfo()[DELTA_TIME] = 4.0;

        p_node->FastGetSolutionStepValue(VELOCITY, 1) =
            Kratos::array_1d<double, 3>{1.0, 2.0, 3.0};
        p_node->FastGetSolutionStepValue(ACCELERATION, 1) =
            Kratos::array_1d<double, 3>{4.0, 5.0, 6.0};
        p_node->FastGetSolutionStepValue(DISPLACEMENT, 1) =
            Kratos::array_1d<double, 3>{7.0, 8.0, 9.0};

        p_node->FastGetSolutionStepValue(WATER_PRESSURE, 1) = 1.0;
        p_node->FastGetSolutionStepValue(WATER_PRESSURE) = 2.0;

        return result;
    }

    template <class T>
    std::vector<std::pair<std::function<void()>, std::function<bool(const Kratos::intrusive_ptr<T> rElement)>>> CreateConditionFunctionsAndChecks()
    {
        std::vector<std::pair<std::function<void()>, std::function<bool(const Kratos::intrusive_ptr<T> rElement)>>> functions_and_checks;

        CompressedMatrix A;
        Vector Dx;
        Vector b;

        auto finalize_solution_step = [this, &A, &Dx, &b]() {
            mScheme.FinalizeSolutionStep(*mrModelPart, A, Dx, b);
        };

        auto initialize_function = [this, &A, &Dx, &b]() {
            mScheme.InitializeSolutionStep(*mrModelPart, A, Dx, b);
        };
        auto initialize_non_linear_iteration = [this, &A, &Dx, &b]() {
            mScheme.InitializeNonLinIteration(*mrModelPart, A, Dx, b);
        };

        auto finalize_non_linear_iteration = [this, &A, &Dx, &b]() {
            mScheme.FinalizeNonLinIteration(*mrModelPart, A, Dx, b);
        };

        auto finalize_function_check = [](const Kratos::intrusive_ptr<T> rElement) {
            return rElement->IsSolutionStepFinalized();
        };

        auto initialize_function_check = [](const Kratos::intrusive_ptr<T> rElement) {
            return rElement->IsSolutionStepInitialized();
        };

        auto initialize_non_linear_iteration_check =
            [](const Kratos::intrusive_ptr<T> rCondition) {
                return rCondition->IsNonLinIterationInitialized();
            };

        auto finalize_non_linear_iteration_check = [](const Kratos::intrusive_ptr<T> rCondition) {
            return rCondition->IsNonLinIterationFinalized();
        };

        functions_and_checks.push_back({finalize_solution_step, finalize_function_check});
        functions_and_checks.push_back({initialize_function, initialize_function_check});
        functions_and_checks.push_back({initialize_non_linear_iteration,
                                        initialize_non_linear_iteration_check});
        functions_and_checks.push_back({finalize_non_linear_iteration,
                                        finalize_non_linear_iteration_check});

        return functions_and_checks;
    }
};

KRATOS_TEST_CASE_IN_SUITE(InitializeUPWScheme_SetsTimeFactors, KratosGeoMechanicsFastSuite)
{
    NewmarkQuasistaticUPwSchemeTester tester;
    tester.mrModelPart->GetProcessInfo()[DELTA_TIME] = 4.0;

    tester.mScheme.Initialize(*tester.mrModelPart);

    // These are the expected numbers according to the SetTimeFactors function
    const double expected_dt_pressure_coefficient = 1.0 / 3.0;
    const double expected_velocity_coefficient = 0.5;
    KRATOS_EXPECT_TRUE(tester.mScheme.SchemeIsInitialized())
    KRATOS_EXPECT_DOUBLE_EQ(tester.mrModelPart->GetProcessInfo()[DT_PRESSURE_COEFFICIENT],
                            expected_dt_pressure_coefficient);
    KRATOS_EXPECT_DOUBLE_EQ(tester.mrModelPart->GetProcessInfo()[VELOCITY_COEFFICIENT],
                            expected_velocity_coefficient);
}

KRATOS_TEST_CASE_IN_SUITE(UPWSchemePredict_UpdatesVariablesDerivatives, KratosGeoMechanicsFastSuite)
{
    NewmarkQuasistaticUPwSchemeTester tester;

    tester.mScheme.Initialize(*tester.mrModelPart); // This is needed to set the time factors
    ModelPart::DofsArrayType dof_set;
    CompressedMatrix A;
    Vector Dx;
    Vector b;

    tester.mScheme.Check(*tester.mrModelPart);
    tester.mScheme.Predict(*tester.mrModelPart, dof_set, A, Dx, b);

    auto expected_acceleration = Kratos::array_1d<double, 3>{-6.75, -9.0, -11.25};
    auto expected_velocity = Kratos::array_1d<double, 3>{-4.5, -6.0, -7.5};
    auto expected_dt_water_pressure = 1.0/3.0;

    KRATOS_EXPECT_EQ(tester.mrModelPart->Nodes()[0].FastGetSolutionStepValue(ACCELERATION), expected_acceleration);
    KRATOS_EXPECT_EQ(tester.mrModelPart->Nodes()[0].FastGetSolutionStepValue(VELOCITY), expected_velocity);
    KRATOS_EXPECT_DOUBLE_EQ(tester.mrModelPart->Nodes()[0].FastGetSolutionStepValue(DT_WATER_PRESSURE), expected_dt_water_pressure);
}

KRATOS_TEST_CASE_IN_SUITE(FinalizeSolutionStepActiveEntities_FinalizesOnlyActiveElements,
                          KratosGeoMechanicsFastSuite)
{
    NewmarkQuasistaticUPwSchemeTester tester;
    auto functions_and_checks = tester.CreateConditionFunctionsAndChecks<SpyElement>();

    for (const auto& [function_on_all_elements, function_has_been_called_on_element] :
         functions_and_checks) {
        tester.Setup();
        auto active_element = Kratos::make_intrusive<SpyElement>();
        active_element->Set(ACTIVE, true);
        active_element->SetId(0);

        auto inactive_element = Kratos::make_intrusive<SpyElement>();
        inactive_element->Set(ACTIVE, false);
        inactive_element->SetId(1);

        tester.mrModelPart->AddElement(active_element);
        tester.mrModelPart->AddElement(inactive_element);

        function_on_all_elements();

        KRATOS_EXPECT_TRUE(function_has_been_called_on_element(active_element))
        KRATOS_EXPECT_FALSE(function_has_been_called_on_element(inactive_element))
    }
}

KRATOS_TEST_CASE_IN_SUITE(FinalizeSolutionStep_FinalizesOnlyActiveConditions, KratosGeoMechanicsFastSuite)
{
    NewmarkQuasistaticUPwSchemeTester tester;
    auto functions_and_checks = tester.CreateConditionFunctionsAndChecks<SpyCondition>();

    for (const auto& [function, function_has_been_called_on_condition] : functions_and_checks) {
        tester.Setup();
        auto active_condition = Kratos::make_intrusive<SpyCondition>();
        active_condition->SetId(0);
        active_condition->Set(ACTIVE, true);

        auto inactive_condition = Kratos::make_intrusive<SpyCondition>();
        inactive_condition->SetId(1);
        inactive_condition->Set(ACTIVE, false);

        tester.mrModelPart->AddCondition(active_condition);
        tester.mrModelPart->AddCondition(inactive_condition);

        function();

        KRATOS_EXPECT_TRUE(function_has_been_called_on_condition(active_condition))
        KRATOS_EXPECT_FALSE(function_has_been_called_on_condition(inactive_condition))
    }
}

} // namespace Kratos::Testing
