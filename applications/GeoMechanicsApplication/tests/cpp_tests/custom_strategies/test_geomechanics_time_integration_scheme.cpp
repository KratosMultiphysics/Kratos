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

#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

#include "containers/model.h"
#include "custom_strategies/schemes/geomechanics_time_integration_scheme.hpp"
#include "spaces/ublas_space.h"
#include "tests/cpp_tests/test_utilities/spy_condition.h"
#include "tests/cpp_tests/test_utilities/spy_element.h"

using namespace Kratos;
using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
using LocalSpaceType  = UblasSpace<double, Matrix, Vector>;

namespace Kratos::Testing
{

// We need this class to test all non-abstract functions of the
// GeoMechanicsTimeIntegrationScheme class. We cannot use the
// GeoMechanicsTimeIntegrationScheme class directly, because it is abstract.
class ConcreteGeoMechanicsTimeIntegrationScheme
    : public GeoMechanicsTimeIntegrationScheme<SparseSpaceType, LocalSpaceType>
{
public:
    ConcreteGeoMechanicsTimeIntegrationScheme()
        : GeoMechanicsTimeIntegrationScheme<SparseSpaceType, LocalSpaceType>({}, {})
    {
    }

protected:
    MOCK_METHOD(void, UpdateVariablesDerivatives, (ModelPart&), (override));
};

class GeoMechanicsSchemeTester
{
public:
    Model                                     mModel;
    ConcreteGeoMechanicsTimeIntegrationScheme mScheme;

    void Setup()
    {
        mModel.Reset();
        mModel.CreateModelPart("dummy", 2);
    }

    template <class T, typename AddComponentToModelPartCallable, typename InitializeComponentsInModelPartCallable>
    void TestFunctionCalledOnComponent_IsCalledOnActiveAndInactiveComponents(
        const AddComponentToModelPartCallable&         rAddComponentTo,
        const InitializeComponentsInModelPartCallable& rInitializeComponentsInModelPart)
    {
        typename T::EquationIdVectorType r_equation_id_vector;
        ProcessInfo                      r_process_info;
        typename T::DofsVectorType       r_dofs_vector;

        Setup();
        auto active_component = Kratos::make_intrusive<T>();
        active_component->SetId(0);
        active_component->Set(ACTIVE, true);

        auto inactive_component = Kratos::make_intrusive<T>();
        inactive_component->SetId(1);
        inactive_component->Set(ACTIVE, false);

        auto& r_model_part = mModel.GetModelPart("dummy");
        rAddComponentTo(r_model_part, active_component);
        rAddComponentTo(r_model_part, inactive_component);

        EXPECT_CALL(*active_component, EquationIdVector(testing::_, testing::_)).Times(1);
        mScheme.EquationId(*active_component.get(), r_equation_id_vector, r_process_info);

        EXPECT_CALL(*inactive_component, EquationIdVector(testing::_, testing::_)).Times(1);
        mScheme.EquationId(*inactive_component.get(), r_equation_id_vector, r_process_info);

        EXPECT_CALL(*active_component, GetDofList(testing::_, testing::_)).Times(1);
        mScheme.GetDofList(*active_component.get(), r_dofs_vector, r_process_info);

        EXPECT_CALL(*inactive_component, GetDofList(testing::_, testing::_)).Times(1);
        mScheme.GetDofList(*inactive_component.get(), r_dofs_vector, r_process_info);

        EXPECT_CALL(*inactive_component, Initialize(testing::_)).Times(1);
        EXPECT_CALL(*active_component, Initialize(testing::_)).Times(1);
        rInitializeComponentsInModelPart(mScheme, r_model_part);
    }

    template <class T>
    void TestFunctionCallOnAllComponents_AreOnlyCalledForActiveComponents()
    {
        CompressedMatrix A;
        Vector           Dx;
        Vector           b;

        auto functions_and_checks = CreateFunctionsAndChecksCalledOnAllComponents<T>(A, Dx, b);

        for (const auto& [function_on_all_elements, function_has_been_called_on_element] : functions_and_checks) {
            Setup();
            auto active_component = Kratos::make_intrusive<T>();
            active_component->Set(ACTIVE, true);
            active_component->SetId(0);

            auto inactive_component = Kratos::make_intrusive<T>();
            inactive_component->Set(ACTIVE, false);
            inactive_component->SetId(1);

            AddComponent(active_component);
            AddComponent(inactive_component);

            function_on_all_elements();

            KRATOS_EXPECT_TRUE(function_has_been_called_on_element(active_component))
            KRATOS_EXPECT_FALSE(function_has_been_called_on_element(inactive_component))
        }
    }

    template <class T>
    std::vector<std::pair<std::function<void()>, std::function<bool(const Kratos::intrusive_ptr<T> rElement)>>> CreateFunctionsAndChecksCalledOnAllComponents(
        CompressedMatrix& A, Vector& Dx, Vector& b)
    {
        std::vector<std::pair<std::function<void()>, std::function<bool(const Kratos::intrusive_ptr<T> rElement)>>> functions_and_checks;

        // Create functions that need to be called in the test
        auto finalize_solution_step = [this, &A, &Dx, &b]() {
            mScheme.FinalizeSolutionStep(GetModelPart(), A, Dx, b);
        };

        auto initialize_function = [this, &A, &Dx, &b]() {
            mScheme.InitializeSolutionStep(GetModelPart(), A, Dx, b);
        };
        auto initialize_non_linear_iteration = [this, &A, &Dx, &b]() {
            mScheme.InitializeNonLinIteration(GetModelPart(), A, Dx, b);
        };

        auto finalize_non_linear_iteration = [this, &A, &Dx, &b]() {
            mScheme.FinalizeNonLinIteration(GetModelPart(), A, Dx, b);
        };

        // Create functions that check if the previously mentioned functions have been called
        auto finalize_function_check = [](const Kratos::intrusive_ptr<T> rElement) {
            return rElement->IsSolutionStepFinalized();
        };

        auto initialize_function_check = [](const Kratos::intrusive_ptr<T> rElement) {
            return rElement->IsSolutionStepInitialized();
        };

        auto initialize_non_linear_iteration_check = [](const Kratos::intrusive_ptr<T> rCondition) {
            return rCondition->IsNonLinIterationInitialized();
        };

        auto finalize_non_linear_iteration_check = [](const Kratos::intrusive_ptr<T> rCondition) {
            return rCondition->IsNonLinIterationFinalized();
        };

        functions_and_checks.push_back({finalize_solution_step, finalize_function_check});
        functions_and_checks.push_back({initialize_function, initialize_function_check});
        functions_and_checks.push_back({initialize_non_linear_iteration, initialize_non_linear_iteration_check});
        functions_and_checks.push_back({finalize_non_linear_iteration, finalize_non_linear_iteration_check});

        return functions_and_checks;
    }

    void AddComponent(ModelPart::ElementType::Pointer element)
    {
        GetModelPart().AddElement(element);
    }

    void AddComponent(ModelPart::ConditionType::Pointer condition)
    {
        GetModelPart().AddCondition(condition);
    }

    ModelPart& GetModelPart() { return mModel.GetModelPart("dummy"); }
};

KRATOS_TEST_CASE_IN_SUITE(FunctionCallsOnAllElements_AreOnlyCalledForActiveElements, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    GeoMechanicsSchemeTester tester;
    tester.TestFunctionCallOnAllComponents_AreOnlyCalledForActiveComponents<SpyElement>();
}

KRATOS_TEST_CASE_IN_SUITE(FunctionCallsOnAllConditions_AreOnlyCalledForActiveConditions,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    GeoMechanicsSchemeTester tester;
    tester.TestFunctionCallOnAllComponents_AreOnlyCalledForActiveComponents<SpyCondition>();
}

KRATOS_TEST_CASE_IN_SUITE(FunctionCalledOnElement_IsCalledOnActiveAndInactiveElements, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    GeoMechanicsSchemeTester tester;
    tester.TestFunctionCalledOnComponent_IsCalledOnActiveAndInactiveComponents<SpyElement>(
        [](auto& rModelPart, auto& rElement) { rModelPart.AddElement(rElement); },
        [](auto& rScheme, auto& rModelPart) { rScheme.InitializeElements(rModelPart); });
}

KRATOS_TEST_CASE_IN_SUITE(FunctionCalledOnCondition_IsCalledOnActiveAndInactiveConditions,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    GeoMechanicsSchemeTester tester;
    tester.TestFunctionCalledOnComponent_IsCalledOnActiveAndInactiveComponents<SpyCondition>(
        [](auto& rModelPart, auto& rCondition) { rModelPart.AddCondition(rCondition); },
        [](auto& rScheme, auto& rModelPart) {
        rScheme.SetElementsAreInitialized(); // Precondition for initializing the conditions
        rScheme.InitializeConditions(rModelPart);
    });
}

KRATOS_TEST_CASE_IN_SUITE(InitializeConditions_Throws_IfElementsNotInitialized, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    GeoMechanicsSchemeTester tester;
    tester.Setup();
    auto& model_part = tester.GetModelPart();
    EXPECT_THROW(tester.mScheme.InitializeConditions(model_part), Kratos::Exception);
}

KRATOS_TEST_CASE_IN_SUITE(ForInvalidBufferSize_CheckGeoMechanicsTimeIntegrationScheme_Throws,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    ConcreteGeoMechanicsTimeIntegrationScheme scheme;

    Model         model;
    constexpr int invalid_buffer_size = 1;
    const auto&   model_part          = model.CreateModelPart("dummy", invalid_buffer_size);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(scheme.Check(model_part),
                                      "insufficient buffer size. Buffer size "
                                      "should be greater than or equal to "
                                      "2. Current size is 1")
}

void TestUpdateForNumberOfThreads(int NumberOfThreads)
{
    GeoMechanicsSchemeTester tester;
    tester.Setup();
    CompressedMatrix         A;
    Vector                   Dx = ZeroVector(3);
    Vector                   b;
    ModelPart::DofsArrayType dofs_array;

    ParallelUtilities::SetNumThreads(NumberOfThreads);

    tester.GetModelPart().AddNodalSolutionStepVariable(WATER_PRESSURE);
    tester.GetModelPart().AddNodalSolutionStepVariable(DISPLACEMENT);

    auto p_node = tester.GetModelPart().CreateNewNode(1, 0.0, 0.0, 0.0);
    p_node->AddDof(WATER_PRESSURE);
    p_node->AddDof(DISPLACEMENT_X);
    p_node->AddDof(DISPLACEMENT_Y);

    auto dof_water_pressure                                     = p_node->pGetDof(WATER_PRESSURE);
    dof_water_pressure->GetSolutionStepValue(WATER_PRESSURE, 0) = 42.0;
    dof_water_pressure->SetEquationId(0);
    dofs_array.push_back(dof_water_pressure);
    Dx[0] = 1.0; // Meaning the updated value = 42.0 + 1.0 = 43.0

    auto dof_displacement                                     = p_node->pGetDof(DISPLACEMENT_X);
    dof_displacement->GetSolutionStepValue(DISPLACEMENT_X, 0) = 3.14;
    dof_displacement->SetEquationId(1);
    dofs_array.push_back(dof_displacement);
    Dx[1] = 6.0; // Meaning the updated value = 3.14 + 6.0 = 9.14

    auto dof_inactive_displacement = p_node->pGetDof(DISPLACEMENT_Y);
    dof_inactive_displacement->GetSolutionStepValue(DISPLACEMENT_Y, 0) = 1.0;
    dof_inactive_displacement->SetEquationId(1);
    dof_inactive_displacement->FixDof();
    dofs_array.push_back(dof_inactive_displacement);
    Dx[2] = 3.0; // Meaning the updated value stays 1.0, since the dof_water_pressure is fixed

    tester.mScheme.Update(tester.GetModelPart(), dofs_array, A, Dx, b);

    KRATOS_EXPECT_DOUBLE_EQ(dofs_array.begin()->GetSolutionStepValue(WATER_PRESSURE, 0), 43.0);
    KRATOS_EXPECT_DOUBLE_EQ((dofs_array.begin() + 1)->GetSolutionStepValue(DISPLACEMENT_X, 0), 9.14);
    KRATOS_EXPECT_DOUBLE_EQ((dofs_array.begin() + 2)->GetSolutionStepValue(DISPLACEMENT_Y, 0), 1.0);
}

KRATOS_TEST_CASE_IN_SUITE(GeoMechanicsTimeIntegrationScheme_GivesCorrectDofs_WhenUpdateIsCalled,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    TestUpdateForNumberOfThreads(1);
    TestUpdateForNumberOfThreads(2);
}

} // namespace Kratos::Testing
