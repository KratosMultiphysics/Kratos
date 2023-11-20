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

#include "containers/model.h"
#include "custom_strategies/schemes/geomechanics_time_integration_scheme.hpp"
#include "geo_mechanics_application_variables.h"
#include "spaces/ublas_space.h"
#include "test_utilities/spy_condition.h"
#include "test_utilities/spy_element.h"

using namespace Kratos;
using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
using LocalSpaceType = UblasSpace<double, Matrix, Vector>;

namespace Kratos::Testing {

class GeoMechanicsSchemeTester {
public:
    Model mModel;
    GeoMechanicsTimeIntegrationScheme<SparseSpaceType, LocalSpaceType> mScheme;

    void Setup()
    {
        mModel.Reset();
        mModel.CreateModelPart("dummy", 2);
        mScheme = GeoMechanicsTimeIntegrationScheme<SparseSpaceType, LocalSpaceType>();
    }

    template <class T>
    void TestFunctionCalledOnComponent_IsOnlyCalledWhenComponentIsActive()
    {
        typename T::EquationIdVectorType r_equation_id_vector;
        ProcessInfo r_process_info;
        typename T::DofsVectorType r_dofs_vector;

        auto functions_and_checks = CreateFunctionsAndChecksCalledOnASingleComponent<T>(
            r_equation_id_vector, r_process_info, r_dofs_vector);

        for (const auto& [function_on_component, function_has_been_called_on_component] :
             functions_and_checks) {
            Setup();
            auto active_component = Kratos::make_intrusive<T>();
            active_component->SetId(0);
            active_component->Set(ACTIVE, true);

            auto inactive_component = Kratos::make_intrusive<T>();
            inactive_component->SetId(1);
            inactive_component->Set(ACTIVE, false);

            function_on_component(active_component);
            function_on_component(inactive_component);

            KRATOS_EXPECT_TRUE(function_has_been_called_on_component(active_component))
            KRATOS_EXPECT_FALSE(function_has_been_called_on_component(inactive_component))
        }
    }

    template <class T>
    void TestFunctionCallOnAllComponents_AreOnlyCalledForActiveComponents()
    {
        CompressedMatrix A;
        Vector Dx;
        Vector b;

        auto functions_and_checks =
            CreateFunctionsAndChecksCalledOnAllComponents<T>(A, Dx, b);

        for (const auto& [function_on_all_elements, function_has_been_called_on_element] :
             functions_and_checks) {
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

    template <typename T>
    std::vector<std::pair<std::function<void(const Kratos::intrusive_ptr<T> Component)>, std::function<bool(const Kratos::intrusive_ptr<T> rElement)>>> CreateFunctionsAndChecksCalledOnASingleComponent(
        typename T::EquationIdVectorType& rEquationIdVector,
        ProcessInfo& rProcessInfo,
        typename T::DofsVectorType& rDofsVector)
    {
        std::vector<std::pair<std::function<void(const Kratos::intrusive_ptr<T> Component)>,
                              std::function<bool(const Kratos::intrusive_ptr<T> Component)>>>
            functions_and_checks;

        auto equation_id = [this, &rEquationIdVector,
                            &rProcessInfo](const Kratos::intrusive_ptr<T> Component) {
            mScheme.EquationId(*Component.get(), rEquationIdVector, rProcessInfo);
        };

        auto equation_id_check = [](const Kratos::intrusive_ptr<T> Component) {
            return Component->IsEquationIdRetrieved();
        };

        auto get_dofs_list = [this, &rDofsVector, &rProcessInfo](
                                 const Kratos::intrusive_ptr<T> Component) {
            mScheme.GetDofList(*Component.get(), rDofsVector, rProcessInfo);
        };

        auto get_dofs_list_check = [](const Kratos::intrusive_ptr<T> Component) {
            return Component->IsGetDofListCalled();
        };

        functions_and_checks.push_back({equation_id, equation_id_check});
        functions_and_checks.push_back({get_dofs_list, get_dofs_list_check});
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

    ModelPart& GetModelPart()
    {
        return mModel.GetModelPart("dummy");
    }
};

KRATOS_TEST_CASE_IN_SUITE(FunctionCallsOnAllElements_AreOnlyCalledForActiveElements,
                          KratosGeoMechanicsFastSuite)
{
    GeoMechanicsSchemeTester tester;
    tester.TestFunctionCallOnAllComponents_AreOnlyCalledForActiveComponents<SpyElement>();
}

KRATOS_TEST_CASE_IN_SUITE(FunctionCallsOnAllConditions_AreOnlyCalledForActiveConditions,
                          KratosGeoMechanicsFastSuite)
{
    GeoMechanicsSchemeTester tester;
    tester.TestFunctionCallOnAllComponents_AreOnlyCalledForActiveComponents<SpyCondition>();
}

KRATOS_TEST_CASE_IN_SUITE(FunctionCalledOnCondition_IsOnlyCalledWhenConditionIsActive,
                          KratosGeoMechanicsFastSuite)
{
    GeoMechanicsSchemeTester tester;
    tester.TestFunctionCalledOnComponent_IsOnlyCalledWhenComponentIsActive<SpyCondition>();
}

KRATOS_TEST_CASE_IN_SUITE(FunctionCalledOnElement_IsOnlyCalledWhenElementIsActive,
                          KratosGeoMechanicsFastSuite)
{
    GeoMechanicsSchemeTester tester;
    tester.TestFunctionCalledOnComponent_IsOnlyCalledWhenComponentIsActive<SpyElement>();
}

KRATOS_TEST_CASE_IN_SUITE(ForInvalidBufferSize_CheckGeoMechanicsTimeIntegrationScheme_Throws,
                          KratosGeoMechanicsFastSuite)
{
    GeoMechanicsTimeIntegrationScheme<SparseSpaceType, LocalSpaceType> scheme;

    Model model;
    constexpr int invalid_buffer_size = 1;
    const auto& model_part = model.CreateModelPart("dummy", invalid_buffer_size);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(scheme.Check(model_part),
                                      "insufficient buffer size. Buffer size "
                                      "should be greater than or equal to "
                                      "2. Current size is 1")
}

} // namespace Kratos::Testing