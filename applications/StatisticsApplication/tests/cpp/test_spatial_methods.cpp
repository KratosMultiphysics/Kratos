//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//
//

// System includes
#include <functional>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/checks.h"
#include "includes/variables.h"
#include "testing/testing.h"

// Application includes
#include "custom_methods/spatial_methods.h"
#include "custom_utilities/test_utilities.h"

namespace Kratos
{
namespace Testing
{
namespace
{
template <typename TContainerType, typename TContainerItemType, template <typename T> typename TDataRetrievalFunctor>
void InitializeSpatialMethodVariables(TContainerType& rContainer)
{
    using random_initializer =
        StatisticsApplicationTestUtilities::RandomInitializer<TContainerType, TContainerItemType, TDataRetrievalFunctor>;
    random_initializer::InitializeVariableWithRandomValues(rContainer, DENSITY);
    random_initializer::InitializeVariableWithRandomValues(rContainer, VELOCITY);
}

template <typename TContainerType, typename TContainerItemType, template <typename T> typename TDataRetrievalFunctor>
void RunSpatialSumMethodTest()
{
    Model test_model;
    ModelPart& test_model_part = test_model.CreateModelPart("test_model_part");
    StatisticsApplicationTestUtilities::AddNodalSolutionStepVariables(test_model_part);
    StatisticsApplicationTestUtilities::CreateModelPart(test_model_part);

    TContainerType& r_container =
        MethodsUtilities::GetDataContainer<TContainerType>(test_model_part);

    InitializeSpatialMethodVariables<TContainerType, TContainerItemType, TDataRetrievalFunctor>(
        r_container);
    using statistics_methods =
        SpatialMethods::ContainerSpatialMethods<TContainerType, TContainerItemType, TDataRetrievalFunctor>;

    const double sum_method_density =
        statistics_methods::CalculateSum(test_model_part, DENSITY);
    const array_1d<double, 3>& sum_method_velocity =
        statistics_methods::CalculateSum(test_model_part, VELOCITY);

    double sum_density = 0.0;
    array_1d<double, 3> sum_velocity = ZeroVector(3);
    for (const TContainerItemType& r_item : r_container)
    {
        sum_density += TDataRetrievalFunctor<TContainerItemType>()(r_item, DENSITY);
        sum_velocity += TDataRetrievalFunctor<TContainerItemType>()(r_item, VELOCITY);
    }

    KRATOS_CHECK_NEAR(sum_density, sum_method_density, 1e-12);
    KRATOS_CHECK_VECTOR_NEAR(sum_velocity, sum_method_velocity, 1e-12);
}

template <typename TContainerType, typename TContainerItemType, typename TInputDataType, typename TOutputDataType, template <typename T> typename TDataRetrievalFunctor>
void RunSpatialNormMethodTest(
    const std::string& rNormType,
    const Variable<TInputDataType>& rVariable,
    const std::function<TOutputDataType(const std::string&, const ModelPart&, const Variable<TInputDataType>&)>& rTestMethod,
    const std::function<void(const TOutputDataType&, const TContainerType&, const std::string&, const Variable<TInputDataType>&)>& rCheckMethod)
{
    Model test_model;
    ModelPart& test_model_part = test_model.CreateModelPart("test_model_part");
    StatisticsApplicationTestUtilities::AddNodalSolutionStepVariables(test_model_part);
    StatisticsApplicationTestUtilities::CreateModelPart(test_model_part);

    TContainerType& r_container =
        MethodsUtilities::GetDataContainer<TContainerType>(test_model_part);

    InitializeSpatialMethodVariables<TContainerType, TContainerItemType, TDataRetrievalFunctor>(
        r_container);

    const auto& r_norm_method =
        MethodsUtilities::GetNormMethod<TInputDataType>(rVariable, rNormType);
    const auto& method_output = rTestMethod(rNormType, test_model_part, rVariable);

    rCheckMethod(method_output, r_container, rNormType, rVariable);
}

template <typename TContainerType, typename TContainerItemType, typename TDataType, template <typename T> typename TDataRetrievalFunctor>
void CheckNormSumMethod(const double& rTestValue,
                        const TContainerType& rContainer,
                        const std::string& rNormType,
                        const Variable<TDataType>& rVariable)
{
    KRATOS_TRY

    double sum = 0.0;
    const auto& norm_method =
        MethodsUtilities::GetNormMethod<TDataType>(rVariable, rNormType);
    for (const TContainerItemType& item : rContainer)
    {
        sum += norm_method(TDataRetrievalFunctor<TContainerItemType>()(item, rVariable));
    }

    KRATOS_CHECK_NEAR(rTestValue, sum, 1e-12);

    KRATOS_CATCH("");
}

template <typename TContainerType, typename TContainerItemType, typename TDataType, template <typename T> typename TDataRetrievalFunctor>
void CheckNormMeanMethod(const double& rTestValue,
                         const TContainerType& rContainer,
                         const std::string& rNormType,
                         const Variable<TDataType>& rVariable)
{
    KRATOS_TRY

    double mean = 0.0;
    const auto& norm_method =
        MethodsUtilities::GetNormMethod<TDataType>(rVariable, rNormType);
    for (const TContainerItemType& item : rContainer)
    {
        mean += norm_method(TDataRetrievalFunctor<TContainerItemType>()(item, rVariable));
    }
    mean /= static_cast<double>(rContainer.size());

    KRATOS_CHECK_NEAR(rTestValue, mean, 1e-12);

    KRATOS_CATCH("");
}

template <typename TContainerType, typename TContainerItemType, typename TDataType, template <typename T> typename TDataRetrievalFunctor>
void CheckNormVarianceMethod(const std::tuple<double, double>& rTestValue,
                             const TContainerType& rContainer,
                             const std::string& rNormType,
                             const Variable<TDataType>& rVariable)
{
    KRATOS_TRY

    double mean = 0.0;
    double variance = 0.0;
    const auto& norm_method =
        MethodsUtilities::GetNormMethod<TDataType>(rVariable, rNormType);
    for (const TContainerItemType& item : rContainer)
    {
        mean += norm_method(TDataRetrievalFunctor<TContainerItemType>()(item, rVariable));
        variance += std::pow(
            norm_method(TDataRetrievalFunctor<TContainerItemType>()(item, rVariable)), 2);
    }
    mean /= static_cast<double>(rContainer.size());
    variance = variance / static_cast<double>(rContainer.size()) - std::pow(mean, 2);

    KRATOS_CHECK_NEAR(std::get<0>(rTestValue), mean, 1e-12);
    KRATOS_CHECK_NEAR(std::get<1>(rTestValue), variance, 1e-12);

    KRATOS_CATCH("");
}

template <typename TContainerType, typename TContainerItemType, template <typename T> typename TDataRetrievalFunctor>
void RunSpatialMeanMethodTest()
{
    Model test_model;
    ModelPart& test_model_part = test_model.CreateModelPart("test_model_part");
    StatisticsApplicationTestUtilities::AddNodalSolutionStepVariables(test_model_part);
    StatisticsApplicationTestUtilities::CreateModelPart(test_model_part);

    TContainerType& r_container =
        MethodsUtilities::GetDataContainer<TContainerType>(test_model_part);

    InitializeSpatialMethodVariables<TContainerType, TContainerItemType, TDataRetrievalFunctor>(
        r_container);
    using statistics_methods =
        SpatialMethods::ContainerSpatialMethods<TContainerType, TContainerItemType, TDataRetrievalFunctor>;

    const double mean_method_density =
        statistics_methods::CalculateMean(test_model_part, DENSITY);
    const array_1d<double, 3>& mean_method_velocity =
        statistics_methods::CalculateMean(test_model_part, VELOCITY);

    double mean_density = 0.0;
    array_1d<double, 3> mean_velocity = ZeroVector(3);
    for (const TContainerItemType& r_item : r_container)
    {
        mean_density += TDataRetrievalFunctor<TContainerItemType>()(r_item, DENSITY);
        mean_velocity += TDataRetrievalFunctor<TContainerItemType>()(r_item, VELOCITY);
    }

    mean_density /= r_container.size();
    mean_velocity /= r_container.size();

    KRATOS_CHECK_NEAR(mean_density, mean_method_density, 1e-12);
    KRATOS_CHECK_VECTOR_NEAR(mean_velocity, mean_method_velocity, 1e-12);
}

template <typename TContainerType, typename TContainerItemType, template <typename T> typename TDataRetrievalFunctor>
void RunSpatialVarianceMethodTest()
{
    Model test_model;
    ModelPart& test_model_part = test_model.CreateModelPart("test_model_part");
    StatisticsApplicationTestUtilities::AddNodalSolutionStepVariables(test_model_part);
    StatisticsApplicationTestUtilities::CreateModelPart(test_model_part);

    TContainerType& r_container =
        MethodsUtilities::GetDataContainer<TContainerType>(test_model_part);

    InitializeSpatialMethodVariables<TContainerType, TContainerItemType, TDataRetrievalFunctor>(
        r_container);
    using statistics_methods =
        SpatialMethods::ContainerSpatialMethods<TContainerType, TContainerItemType, TDataRetrievalFunctor>;

    const std::tuple<double, double>& method_density =
        statistics_methods::CalculateVariance(test_model_part, DENSITY);
    const std::tuple<array_1d<double, 3>, array_1d<double, 3>>& method_velocity =
        statistics_methods::CalculateVariance(test_model_part, VELOCITY);

    const double mean_method_density = std::get<0>(method_density);
    const double variance_method_density = std::get<1>(method_density);

    const array_1d<double, 3>& mean_method_velocity = std::get<0>(method_velocity);
    const array_1d<double, 3>& variance_method_velocity = std::get<1>(method_velocity);

    double mean_density = 0.0;
    array_1d<double, 3> mean_velocity = ZeroVector(3);
    double variance_density = 0.0;
    array_1d<double, 3> variance_velocity = ZeroVector(3);
    for (const TContainerItemType& r_item : r_container)
    {
        mean_density += TDataRetrievalFunctor<TContainerItemType>()(r_item, DENSITY);
        mean_velocity += TDataRetrievalFunctor<TContainerItemType>()(r_item, VELOCITY);
        variance_density += MethodsUtilities::RaiseToPower(
            TDataRetrievalFunctor<TContainerItemType>()(r_item, DENSITY), 2);
        variance_velocity += MethodsUtilities::RaiseToPower(
            TDataRetrievalFunctor<TContainerItemType>()(r_item, VELOCITY), 2);
    }

    const int n = r_container.size();
    mean_density /= n;
    mean_velocity /= n;
    variance_density =
        variance_density / n - MethodsUtilities::RaiseToPower(mean_density, 2);
    variance_velocity =
        variance_velocity / n - MethodsUtilities::RaiseToPower(mean_velocity, 2);

    KRATOS_CHECK_NEAR(mean_density, mean_method_density, 1e-12);
    KRATOS_CHECK_VECTOR_NEAR(mean_velocity, mean_method_velocity, 1e-12);
    KRATOS_CHECK_NEAR(variance_density, variance_method_density, 1e-12);
    KRATOS_CHECK_VECTOR_NEAR(variance_velocity, variance_method_velocity, 1e-12);
}

template <typename TContainerType, typename TContainerItemType, template <typename T> typename TDataRetrievalFunctor>
void RunSpatialMinMethodTest(const std::string& rNormType)
{
    Model test_model;
    ModelPart& test_model_part = test_model.CreateModelPart("test_model_part");
    StatisticsApplicationTestUtilities::AddNodalSolutionStepVariables(test_model_part);
    StatisticsApplicationTestUtilities::CreateModelPart(test_model_part);

    TContainerType& r_container =
        MethodsUtilities::GetDataContainer<TContainerType>(test_model_part);

    InitializeSpatialMethodVariables<TContainerType, TContainerItemType, TDataRetrievalFunctor>(
        r_container);
    using statistics_methods =
        SpatialMethods::ContainerSpatialMethods<TContainerType, TContainerItemType, TDataRetrievalFunctor>;

    const double max_value = std::numeric_limits<double>::max();

    if (rNormType == "magnitude" || rNormType == "value")
    {
        double min_density{max_value};
        std::size_t min_density_id{0};

        const std::tuple<double, std::size_t>& method_density =
            statistics_methods::GetMin(rNormType, test_model_part, DENSITY);
        const double min_method_density = std::get<0>(method_density);
        const double min_method_density_id = std::get<1>(method_density);
        for (const TContainerItemType& r_item : r_container)
        {
            const double current_density =
                TDataRetrievalFunctor<TContainerItemType>()(r_item, DENSITY);
            double current_density_norm = 0.0;
            if (rNormType == "magnitude")
            {
                current_density_norm = std::abs(current_density);
            }
            else if (rNormType == "value")
            {
                current_density_norm = current_density;
            }

            if (current_density_norm < min_density)
            {
                min_density = current_density_norm;
                min_density_id = r_item.Id();
            }
        }

        KRATOS_CHECK_NEAR(min_density, min_method_density, 1e-12);
        KRATOS_CHECK_EQUAL(min_density_id, min_method_density_id);
    }

    if (rNormType == "magnitude" || rNormType == "component_x" ||
        rNormType == "component_y" || rNormType == "component_z")
    {
        const std::tuple<double, std::size_t>& method_velocity =
            statistics_methods::GetMin(rNormType, test_model_part, VELOCITY);
        const double min_method_velocity = std::get<0>(method_velocity);
        const double min_method_velocity_id = std::get<1>(method_velocity);

        double min_velocity{max_value};
        std::size_t min_velocity_id{0};

        for (const TContainerItemType& r_item : r_container)
        {
            const array_1d<double, 3>& r_current_velocity =
                TDataRetrievalFunctor<TContainerItemType>()(r_item, VELOCITY);
            double current_velocity_norm = 0.0;
            if (rNormType == "magnitude")
            {
                current_velocity_norm = norm_2(r_current_velocity);
            }
            else if (rNormType == "component_x")
            {
                current_velocity_norm = r_current_velocity[0];
            }
            else if (rNormType == "component_y")
            {
                current_velocity_norm = r_current_velocity[1];
            }
            else if (rNormType == "component_z")
            {
                current_velocity_norm = r_current_velocity[2];
            }

            if (current_velocity_norm < min_velocity)
            {
                min_velocity = current_velocity_norm;
                min_velocity_id = r_item.Id();
            }
        }

        KRATOS_CHECK_NEAR(min_velocity, min_method_velocity, 1e-12);
        KRATOS_CHECK_EQUAL(min_velocity_id, min_method_velocity_id);
    }
}

template <typename TContainerType, typename TContainerItemType, template <typename T> typename TDataRetrievalFunctor>
void RunSpatialMaxMethodTest(const std::string& rNormType)
{
    Model test_model;
    ModelPart& test_model_part = test_model.CreateModelPart("test_model_part");
    StatisticsApplicationTestUtilities::AddNodalSolutionStepVariables(test_model_part);
    StatisticsApplicationTestUtilities::CreateModelPart(test_model_part);

    TContainerType& r_container =
        MethodsUtilities::GetDataContainer<TContainerType>(test_model_part);

    InitializeSpatialMethodVariables<TContainerType, TContainerItemType, TDataRetrievalFunctor>(
        r_container);
    using statistics_methods =
        SpatialMethods::ContainerSpatialMethods<TContainerType, TContainerItemType, TDataRetrievalFunctor>;

    const double min_value = std::numeric_limits<double>::lowest();

    if (rNormType == "magnitude" || rNormType == "value")
    {
        double max_density{min_value};
        std::size_t max_density_id{0};

        const std::tuple<double, std::size_t>& method_density =
            statistics_methods::GetMax(rNormType, test_model_part, DENSITY);
        const double max_method_density = std::get<0>(method_density);
        const std::size_t max_method_density_id = std::get<1>(method_density);

        for (const TContainerItemType& r_item : r_container)
        {
            const double current_density =
                TDataRetrievalFunctor<TContainerItemType>()(r_item, DENSITY);
            double current_density_norm = 0.0;
            if (rNormType == "magnitude")
            {
                current_density_norm = std::abs(current_density);
            }
            else if (rNormType == "value")
            {
                current_density_norm = current_density;
            }

            if (current_density_norm > max_density)
            {
                max_density = current_density_norm;
                max_density_id = r_item.Id();
            }
        }
        KRATOS_CHECK_NEAR(max_density, max_method_density, 1e-12);
        KRATOS_CHECK_EQUAL(max_density_id, max_method_density_id);
    }

    if (rNormType == "magnitude" || rNormType == "component_x" ||
        rNormType == "component_y" || rNormType == "component_z")
    {
        const std::tuple<double, std::size_t>& method_velocity =
            statistics_methods::GetMax(rNormType, test_model_part, VELOCITY);
        const double max_method_velocity = std::get<0>(method_velocity);
        const std::size_t max_method_velocity_id = std::get<1>(method_velocity);

        double max_velocity{min_value};
        std::size_t max_velocity_id{0};

        for (const TContainerItemType& r_item : r_container)
        {
            const array_1d<double, 3>& r_current_velocity =
                TDataRetrievalFunctor<TContainerItemType>()(r_item, VELOCITY);
            double current_velocity_norm = 0.0;
            if (rNormType == "magnitude")
            {
                current_velocity_norm = norm_2(r_current_velocity);
            }
            else if (rNormType == "component_x")
            {
                current_velocity_norm = r_current_velocity[0];
            }
            else if (rNormType == "component_y")
            {
                current_velocity_norm = r_current_velocity[1];
            }
            else if (rNormType == "component_z")
            {
                current_velocity_norm = r_current_velocity[2];
            }

            if (current_velocity_norm > max_velocity)
            {
                max_velocity = current_velocity_norm;
                max_velocity_id = r_item.Id();
            }
        }

        KRATOS_CHECK_NEAR(max_velocity, max_method_velocity, 1e-12);
        KRATOS_CHECK_EQUAL(max_velocity_id, max_method_velocity_id);
    }
}

} // namespace

KRATOS_TEST_CASE_IN_SUITE(SpatialSumMethod, KratosStatisticsFastSuite)
{
    RunSpatialSumMethodTest<ModelPart::NodesContainerType, ModelPart::NodeType,
                            MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>();
    RunSpatialSumMethodTest<ModelPart::ConditionsContainerType, ModelPart::ConditionType,
                            MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>();
    RunSpatialSumMethodTest<ModelPart::ElementsContainerType, ModelPart::ElementType,
                            MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>();
    RunSpatialSumMethodTest<ModelPart::NodesContainerType, ModelPart::NodeType,
                            MethodsUtilities::HistoricalDataValueRetrievalFunctor>();
}

KRATOS_TEST_CASE_IN_SUITE(SpatialNormSumMethod, KratosStatisticsFastSuite)
{
    using node_type = ModelPart::NodeType;
    using nodes_container_type = ModelPart::NodesContainerType;
    using condition_type = ModelPart::ConditionType;
    using conditions_container_type = ModelPart::ConditionsContainerType;
    using element_type = ModelPart::ElementType;
    using elements_container_type = ModelPart::ElementsContainerType;

    // checking for doubles
    using variable_type_scalar = double;
    const auto& test_nodal_historical_method_scalar =
        SpatialMethods::ContainerSpatialMethods<nodes_container_type, node_type,
                                                MethodsUtilities::HistoricalDataValueRetrievalFunctor>::CalculateNormSum<variable_type_scalar>;
    const auto& check_nodal_historical_method_scalar =
        CheckNormSumMethod<nodes_container_type, node_type, variable_type_scalar, MethodsUtilities::HistoricalDataValueRetrievalFunctor>;
    RunSpatialNormMethodTest<nodes_container_type, node_type, variable_type_scalar, double, MethodsUtilities::HistoricalDataValueRetrievalFunctor>(
        "value", DENSITY, test_nodal_historical_method_scalar,
        check_nodal_historical_method_scalar);
    RunSpatialNormMethodTest<nodes_container_type, node_type, variable_type_scalar, double, MethodsUtilities::HistoricalDataValueRetrievalFunctor>(
        "magnitude", DENSITY, test_nodal_historical_method_scalar,
        check_nodal_historical_method_scalar);

    const auto& test_nodal_non_historical_method_scalar =
        SpatialMethods::ContainerSpatialMethods<nodes_container_type, node_type,
                                                MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>::CalculateNormSum<variable_type_scalar>;
    const auto& check_nodal_non_historical_method_scalar =
        CheckNormSumMethod<nodes_container_type, node_type, variable_type_scalar, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>;
    RunSpatialNormMethodTest<nodes_container_type, node_type, variable_type_scalar, double, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "value", DENSITY, test_nodal_non_historical_method_scalar,
        check_nodal_non_historical_method_scalar);
    RunSpatialNormMethodTest<nodes_container_type, node_type, variable_type_scalar, double, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "magnitude", DENSITY, test_nodal_non_historical_method_scalar,
        check_nodal_non_historical_method_scalar);

    const auto& test_condition_non_historical_method_scalar =
        SpatialMethods::ContainerSpatialMethods<conditions_container_type, condition_type,
                                                MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>::CalculateNormSum<variable_type_scalar>;
    const auto& check_condition_non_historical_method_scalar =
        CheckNormSumMethod<conditions_container_type, condition_type, variable_type_scalar, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>;
    RunSpatialNormMethodTest<conditions_container_type, condition_type, variable_type_scalar,
                             double, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "value", DENSITY, test_condition_non_historical_method_scalar,
        check_condition_non_historical_method_scalar);
    RunSpatialNormMethodTest<conditions_container_type, condition_type, variable_type_scalar,
                             double, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "magnitude", DENSITY, test_condition_non_historical_method_scalar,
        check_condition_non_historical_method_scalar);

    const auto& test_element_non_historical_method_scalar =
        SpatialMethods::ContainerSpatialMethods<elements_container_type, element_type,
                                                MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>::CalculateNormSum<variable_type_scalar>;
    const auto& check_element_non_historical_method_scalar =
        CheckNormSumMethod<elements_container_type, element_type, variable_type_scalar, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>;
    RunSpatialNormMethodTest<elements_container_type, element_type, variable_type_scalar,
                             double, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "value", DENSITY, test_element_non_historical_method_scalar,
        check_element_non_historical_method_scalar);
    RunSpatialNormMethodTest<elements_container_type, element_type, variable_type_scalar,
                             double, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "magnitude", DENSITY, test_element_non_historical_method_scalar,
        check_element_non_historical_method_scalar);

    using variable_type_array = array_1d<double, 3>;
    const auto& test_nodal_historical_method_array =
        SpatialMethods::ContainerSpatialMethods<nodes_container_type, node_type,
                                                MethodsUtilities::HistoricalDataValueRetrievalFunctor>::CalculateNormSum<variable_type_array>;
    const auto& check_nodal_historical_method_array =
        CheckNormSumMethod<nodes_container_type, node_type, variable_type_array, MethodsUtilities::HistoricalDataValueRetrievalFunctor>;
    RunSpatialNormMethodTest<nodes_container_type, node_type, variable_type_array, double, MethodsUtilities::HistoricalDataValueRetrievalFunctor>(
        "component_x", VELOCITY, test_nodal_historical_method_array,
        check_nodal_historical_method_array);
    RunSpatialNormMethodTest<nodes_container_type, node_type, variable_type_array, double, MethodsUtilities::HistoricalDataValueRetrievalFunctor>(
        "component_y", VELOCITY, test_nodal_historical_method_array,
        check_nodal_historical_method_array);
    RunSpatialNormMethodTest<nodes_container_type, node_type, variable_type_array, double, MethodsUtilities::HistoricalDataValueRetrievalFunctor>(
        "component_z", VELOCITY, test_nodal_historical_method_array,
        check_nodal_historical_method_array);
    RunSpatialNormMethodTest<nodes_container_type, node_type, variable_type_array, double, MethodsUtilities::HistoricalDataValueRetrievalFunctor>(
        "magnitude", VELOCITY, test_nodal_historical_method_array,
        check_nodal_historical_method_array);

    const auto& test_nodal_non_historical_method_array =
        SpatialMethods::ContainerSpatialMethods<nodes_container_type, node_type,
                                                MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>::CalculateNormSum<variable_type_array>;
    const auto& check_nodal_non_historical_method_array =
        CheckNormSumMethod<nodes_container_type, node_type, variable_type_array, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>;
    RunSpatialNormMethodTest<nodes_container_type, node_type, variable_type_array, double, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "component_x", VELOCITY, test_nodal_non_historical_method_array,
        check_nodal_non_historical_method_array);
    RunSpatialNormMethodTest<nodes_container_type, node_type, variable_type_array, double, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "component_y", VELOCITY, test_nodal_non_historical_method_array,
        check_nodal_non_historical_method_array);
    RunSpatialNormMethodTest<nodes_container_type, node_type, variable_type_array, double, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "component_z", VELOCITY, test_nodal_non_historical_method_array,
        check_nodal_non_historical_method_array);
    RunSpatialNormMethodTest<nodes_container_type, node_type, variable_type_array, double, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "magnitude", VELOCITY, test_nodal_non_historical_method_array,
        check_nodal_non_historical_method_array);

    const auto& test_condition_non_historical_method_array =
        SpatialMethods::ContainerSpatialMethods<conditions_container_type, condition_type,
                                                MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>::CalculateNormSum<variable_type_array>;
    const auto& check_condition_non_historical_method_array =
        CheckNormSumMethod<conditions_container_type, condition_type, variable_type_array, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>;
    RunSpatialNormMethodTest<conditions_container_type, condition_type, variable_type_array,
                             double, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "component_x", VELOCITY, test_condition_non_historical_method_array,
        check_condition_non_historical_method_array);
    RunSpatialNormMethodTest<conditions_container_type, condition_type, variable_type_array,
                             double, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "component_y", VELOCITY, test_condition_non_historical_method_array,
        check_condition_non_historical_method_array);
    RunSpatialNormMethodTest<conditions_container_type, condition_type, variable_type_array,
                             double, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "component_z", VELOCITY, test_condition_non_historical_method_array,
        check_condition_non_historical_method_array);
    RunSpatialNormMethodTest<conditions_container_type, condition_type, variable_type_array,
                             double, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "magnitude", VELOCITY, test_condition_non_historical_method_array,
        check_condition_non_historical_method_array);

    const auto& test_element_non_historical_method_array =
        SpatialMethods::ContainerSpatialMethods<elements_container_type, element_type,
                                                MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>::CalculateNormSum<variable_type_array>;
    const auto& check_element_non_historical_method_array =
        CheckNormSumMethod<elements_container_type, element_type, variable_type_array, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>;
    RunSpatialNormMethodTest<elements_container_type, element_type, variable_type_array,
                             double, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "component_x", VELOCITY, test_element_non_historical_method_array,
        check_element_non_historical_method_array);
    RunSpatialNormMethodTest<elements_container_type, element_type, variable_type_array,
                             double, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "component_y", VELOCITY, test_element_non_historical_method_array,
        check_element_non_historical_method_array);
    RunSpatialNormMethodTest<elements_container_type, element_type, variable_type_array,
                             double, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "component_z", VELOCITY, test_element_non_historical_method_array,
        check_element_non_historical_method_array);
    RunSpatialNormMethodTest<elements_container_type, element_type, variable_type_array,
                             double, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "magnitude", VELOCITY, test_element_non_historical_method_array,
        check_element_non_historical_method_array);
}

KRATOS_TEST_CASE_IN_SUITE(SpatialMeanMethod, KratosStatisticsFastSuite)
{
    RunSpatialMeanMethodTest<ModelPart::NodesContainerType, ModelPart::NodeType,
                             MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>();
    RunSpatialMeanMethodTest<ModelPart::ConditionsContainerType, ModelPart::ConditionType,
                             MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>();
    RunSpatialMeanMethodTest<ModelPart::ElementsContainerType, ModelPart::ElementType,
                             MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>();
    RunSpatialMeanMethodTest<ModelPart::NodesContainerType, ModelPart::NodeType,
                             MethodsUtilities::HistoricalDataValueRetrievalFunctor>();
}

KRATOS_TEST_CASE_IN_SUITE(SpatialNormMeanMethod, KratosStatisticsFastSuite)
{
    using node_type = ModelPart::NodeType;
    using nodes_container_type = ModelPart::NodesContainerType;
    using condition_type = ModelPart::ConditionType;
    using conditions_container_type = ModelPart::ConditionsContainerType;
    using element_type = ModelPart::ElementType;
    using elements_container_type = ModelPart::ElementsContainerType;

    // checking for doubles
    using variable_type_scalar = double;
    const auto& test_nodal_historical_method_scalar =
        SpatialMethods::ContainerSpatialMethods<nodes_container_type, node_type,
                                                MethodsUtilities::HistoricalDataValueRetrievalFunctor>::CalculateNormMean<variable_type_scalar>;
    const auto& check_nodal_historical_method_scalar =
        CheckNormMeanMethod<nodes_container_type, node_type, variable_type_scalar, MethodsUtilities::HistoricalDataValueRetrievalFunctor>;
    RunSpatialNormMethodTest<nodes_container_type, node_type, variable_type_scalar, double, MethodsUtilities::HistoricalDataValueRetrievalFunctor>(
        "value", DENSITY, test_nodal_historical_method_scalar,
        check_nodal_historical_method_scalar);
    RunSpatialNormMethodTest<nodes_container_type, node_type, variable_type_scalar, double, MethodsUtilities::HistoricalDataValueRetrievalFunctor>(
        "magnitude", DENSITY, test_nodal_historical_method_scalar,
        check_nodal_historical_method_scalar);

    const auto& test_nodal_non_historical_method_scalar =
        SpatialMethods::ContainerSpatialMethods<nodes_container_type, node_type,
                                                MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>::CalculateNormMean<variable_type_scalar>;
    const auto& check_nodal_non_historical_method_scalar =
        CheckNormMeanMethod<nodes_container_type, node_type, variable_type_scalar, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>;
    RunSpatialNormMethodTest<nodes_container_type, node_type, variable_type_scalar, double, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "value", DENSITY, test_nodal_non_historical_method_scalar,
        check_nodal_non_historical_method_scalar);
    RunSpatialNormMethodTest<nodes_container_type, node_type, variable_type_scalar, double, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "magnitude", DENSITY, test_nodal_non_historical_method_scalar,
        check_nodal_non_historical_method_scalar);

    const auto& test_condition_non_historical_method_scalar =
        SpatialMethods::ContainerSpatialMethods<conditions_container_type, condition_type,
                                                MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>::CalculateNormMean<variable_type_scalar>;
    const auto& check_condition_non_historical_method_scalar =
        CheckNormMeanMethod<conditions_container_type, condition_type, variable_type_scalar, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>;
    RunSpatialNormMethodTest<conditions_container_type, condition_type, variable_type_scalar,
                             double, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "value", DENSITY, test_condition_non_historical_method_scalar,
        check_condition_non_historical_method_scalar);
    RunSpatialNormMethodTest<conditions_container_type, condition_type, variable_type_scalar,
                             double, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "magnitude", DENSITY, test_condition_non_historical_method_scalar,
        check_condition_non_historical_method_scalar);

    const auto& test_element_non_historical_method_scalar =
        SpatialMethods::ContainerSpatialMethods<elements_container_type, element_type,
                                                MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>::CalculateNormMean<variable_type_scalar>;
    const auto& check_element_non_historical_method_scalar =
        CheckNormMeanMethod<elements_container_type, element_type, variable_type_scalar, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>;
    RunSpatialNormMethodTest<elements_container_type, element_type, variable_type_scalar,
                             double, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "value", DENSITY, test_element_non_historical_method_scalar,
        check_element_non_historical_method_scalar);
    RunSpatialNormMethodTest<elements_container_type, element_type, variable_type_scalar,
                             double, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "magnitude", DENSITY, test_element_non_historical_method_scalar,
        check_element_non_historical_method_scalar);

    using variable_type_array = array_1d<double, 3>;
    const auto& test_nodal_historical_method_array =
        SpatialMethods::ContainerSpatialMethods<nodes_container_type, node_type,
                                                MethodsUtilities::HistoricalDataValueRetrievalFunctor>::CalculateNormMean<variable_type_array>;
    const auto& check_nodal_historical_method_array =
        CheckNormMeanMethod<nodes_container_type, node_type, variable_type_array, MethodsUtilities::HistoricalDataValueRetrievalFunctor>;
    RunSpatialNormMethodTest<nodes_container_type, node_type, variable_type_array, double, MethodsUtilities::HistoricalDataValueRetrievalFunctor>(
        "component_x", VELOCITY, test_nodal_historical_method_array,
        check_nodal_historical_method_array);
    RunSpatialNormMethodTest<nodes_container_type, node_type, variable_type_array, double, MethodsUtilities::HistoricalDataValueRetrievalFunctor>(
        "component_y", VELOCITY, test_nodal_historical_method_array,
        check_nodal_historical_method_array);
    RunSpatialNormMethodTest<nodes_container_type, node_type, variable_type_array, double, MethodsUtilities::HistoricalDataValueRetrievalFunctor>(
        "component_z", VELOCITY, test_nodal_historical_method_array,
        check_nodal_historical_method_array);
    RunSpatialNormMethodTest<nodes_container_type, node_type, variable_type_array, double, MethodsUtilities::HistoricalDataValueRetrievalFunctor>(
        "magnitude", VELOCITY, test_nodal_historical_method_array,
        check_nodal_historical_method_array);

    const auto& test_nodal_non_historical_method_array =
        SpatialMethods::ContainerSpatialMethods<nodes_container_type, node_type,
                                                MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>::CalculateNormMean<variable_type_array>;
    const auto& check_nodal_non_historical_method_array =
        CheckNormMeanMethod<nodes_container_type, node_type, variable_type_array, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>;
    RunSpatialNormMethodTest<nodes_container_type, node_type, variable_type_array, double, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "component_x", VELOCITY, test_nodal_non_historical_method_array,
        check_nodal_non_historical_method_array);
    RunSpatialNormMethodTest<nodes_container_type, node_type, variable_type_array, double, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "component_y", VELOCITY, test_nodal_non_historical_method_array,
        check_nodal_non_historical_method_array);
    RunSpatialNormMethodTest<nodes_container_type, node_type, variable_type_array, double, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "component_z", VELOCITY, test_nodal_non_historical_method_array,
        check_nodal_non_historical_method_array);
    RunSpatialNormMethodTest<nodes_container_type, node_type, variable_type_array, double, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "magnitude", VELOCITY, test_nodal_non_historical_method_array,
        check_nodal_non_historical_method_array);

    const auto& test_condition_non_historical_method_array =
        SpatialMethods::ContainerSpatialMethods<conditions_container_type, condition_type,
                                                MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>::CalculateNormMean<variable_type_array>;
    const auto& check_condition_non_historical_method_array =
        CheckNormMeanMethod<conditions_container_type, condition_type, variable_type_array, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>;
    RunSpatialNormMethodTest<conditions_container_type, condition_type, variable_type_array,
                             double, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "component_x", VELOCITY, test_condition_non_historical_method_array,
        check_condition_non_historical_method_array);
    RunSpatialNormMethodTest<conditions_container_type, condition_type, variable_type_array,
                             double, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "component_y", VELOCITY, test_condition_non_historical_method_array,
        check_condition_non_historical_method_array);
    RunSpatialNormMethodTest<conditions_container_type, condition_type, variable_type_array,
                             double, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "component_z", VELOCITY, test_condition_non_historical_method_array,
        check_condition_non_historical_method_array);
    RunSpatialNormMethodTest<conditions_container_type, condition_type, variable_type_array,
                             double, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "magnitude", VELOCITY, test_condition_non_historical_method_array,
        check_condition_non_historical_method_array);

    const auto& test_element_non_historical_method_array =
        SpatialMethods::ContainerSpatialMethods<elements_container_type, element_type,
                                                MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>::CalculateNormMean<variable_type_array>;
    const auto& check_element_non_historical_method_array =
        CheckNormMeanMethod<elements_container_type, element_type, variable_type_array, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>;
    RunSpatialNormMethodTest<elements_container_type, element_type, variable_type_array,
                             double, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "component_x", VELOCITY, test_element_non_historical_method_array,
        check_element_non_historical_method_array);
    RunSpatialNormMethodTest<elements_container_type, element_type, variable_type_array,
                             double, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "component_y", VELOCITY, test_element_non_historical_method_array,
        check_element_non_historical_method_array);
    RunSpatialNormMethodTest<elements_container_type, element_type, variable_type_array,
                             double, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "component_z", VELOCITY, test_element_non_historical_method_array,
        check_element_non_historical_method_array);
    RunSpatialNormMethodTest<elements_container_type, element_type, variable_type_array,
                             double, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "magnitude", VELOCITY, test_element_non_historical_method_array,
        check_element_non_historical_method_array);
}

KRATOS_TEST_CASE_IN_SUITE(SpatialVarianceMethod, KratosStatisticsFastSuite)
{
    RunSpatialVarianceMethodTest<ModelPart::NodesContainerType, ModelPart::NodeType,
                                 MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>();
    RunSpatialVarianceMethodTest<ModelPart::ConditionsContainerType, ModelPart::ConditionType,
                                 MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>();
    RunSpatialVarianceMethodTest<ModelPart::ElementsContainerType, ModelPart::ElementType,
                                 MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>();
    RunSpatialVarianceMethodTest<ModelPart::NodesContainerType, ModelPart::NodeType,
                                 MethodsUtilities::HistoricalDataValueRetrievalFunctor>();
}

KRATOS_TEST_CASE_IN_SUITE(SpatialNormVarianceMethod, KratosStatisticsFastSuite)
{
    using node_type = ModelPart::NodeType;
    using nodes_container_type = ModelPart::NodesContainerType;
    using condition_type = ModelPart::ConditionType;
    using conditions_container_type = ModelPart::ConditionsContainerType;
    using element_type = ModelPart::ElementType;
    using elements_container_type = ModelPart::ElementsContainerType;

    // checking for doubles
    using variable_type_scalar = double;
    const auto& test_nodal_historical_method_scalar =
        SpatialMethods::ContainerSpatialMethods<nodes_container_type, node_type,
                                                MethodsUtilities::HistoricalDataValueRetrievalFunctor>::CalculateNormVariance<variable_type_scalar>;
    const auto& check_nodal_historical_method_scalar =
        CheckNormVarianceMethod<nodes_container_type, node_type, variable_type_scalar, MethodsUtilities::HistoricalDataValueRetrievalFunctor>;
    RunSpatialNormMethodTest<nodes_container_type, node_type, variable_type_scalar,
                             std::tuple<double, double>, MethodsUtilities::HistoricalDataValueRetrievalFunctor>(
        "value", DENSITY, test_nodal_historical_method_scalar,
        check_nodal_historical_method_scalar);
    RunSpatialNormMethodTest<nodes_container_type, node_type, variable_type_scalar,
                             std::tuple<double, double>, MethodsUtilities::HistoricalDataValueRetrievalFunctor>(
        "magnitude", DENSITY, test_nodal_historical_method_scalar,
        check_nodal_historical_method_scalar);

    const auto& test_nodal_non_historical_method_scalar =
        SpatialMethods::ContainerSpatialMethods<nodes_container_type, node_type,
                                                MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>::CalculateNormVariance<variable_type_scalar>;
    const auto& check_nodal_non_historical_method_scalar =
        CheckNormVarianceMethod<nodes_container_type, node_type, variable_type_scalar, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>;
    RunSpatialNormMethodTest<nodes_container_type, node_type, variable_type_scalar,
                             std::tuple<double, double>, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "value", DENSITY, test_nodal_non_historical_method_scalar,
        check_nodal_non_historical_method_scalar);
    RunSpatialNormMethodTest<nodes_container_type, node_type, variable_type_scalar,
                             std::tuple<double, double>, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "magnitude", DENSITY, test_nodal_non_historical_method_scalar,
        check_nodal_non_historical_method_scalar);

    const auto& test_condition_non_historical_method_scalar =
        SpatialMethods::ContainerSpatialMethods<conditions_container_type, condition_type,
                                                MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>::CalculateNormVariance<variable_type_scalar>;
    const auto& check_condition_non_historical_method_scalar =
        CheckNormVarianceMethod<conditions_container_type, condition_type, variable_type_scalar, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>;
    RunSpatialNormMethodTest<conditions_container_type, condition_type, variable_type_scalar,
                             std::tuple<double, double>, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "value", DENSITY, test_condition_non_historical_method_scalar,
        check_condition_non_historical_method_scalar);
    RunSpatialNormMethodTest<conditions_container_type, condition_type, variable_type_scalar,
                             std::tuple<double, double>, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "magnitude", DENSITY, test_condition_non_historical_method_scalar,
        check_condition_non_historical_method_scalar);

    const auto& test_element_non_historical_method_scalar =
        SpatialMethods::ContainerSpatialMethods<elements_container_type, element_type,
                                                MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>::CalculateNormVariance<variable_type_scalar>;
    const auto& check_element_non_historical_method_scalar =
        CheckNormVarianceMethod<elements_container_type, element_type, variable_type_scalar, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>;
    RunSpatialNormMethodTest<elements_container_type, element_type, variable_type_scalar,
                             std::tuple<double, double>, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "value", DENSITY, test_element_non_historical_method_scalar,
        check_element_non_historical_method_scalar);
    RunSpatialNormMethodTest<elements_container_type, element_type, variable_type_scalar,
                             std::tuple<double, double>, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "magnitude", DENSITY, test_element_non_historical_method_scalar,
        check_element_non_historical_method_scalar);

    using variable_type_array = array_1d<double, 3>;
    const auto& test_nodal_historical_method_array =
        SpatialMethods::ContainerSpatialMethods<nodes_container_type, node_type,
                                                MethodsUtilities::HistoricalDataValueRetrievalFunctor>::CalculateNormVariance<variable_type_array>;
    const auto& check_nodal_historical_method_array =
        CheckNormVarianceMethod<nodes_container_type, node_type, variable_type_array, MethodsUtilities::HistoricalDataValueRetrievalFunctor>;
    RunSpatialNormMethodTest<nodes_container_type, node_type, variable_type_array,
                             std::tuple<double, double>, MethodsUtilities::HistoricalDataValueRetrievalFunctor>(
        "component_x", VELOCITY, test_nodal_historical_method_array,
        check_nodal_historical_method_array);
    RunSpatialNormMethodTest<nodes_container_type, node_type, variable_type_array,
                             std::tuple<double, double>, MethodsUtilities::HistoricalDataValueRetrievalFunctor>(
        "component_y", VELOCITY, test_nodal_historical_method_array,
        check_nodal_historical_method_array);
    RunSpatialNormMethodTest<nodes_container_type, node_type, variable_type_array,
                             std::tuple<double, double>, MethodsUtilities::HistoricalDataValueRetrievalFunctor>(
        "component_z", VELOCITY, test_nodal_historical_method_array,
        check_nodal_historical_method_array);
    RunSpatialNormMethodTest<nodes_container_type, node_type, variable_type_array,
                             std::tuple<double, double>, MethodsUtilities::HistoricalDataValueRetrievalFunctor>(
        "magnitude", VELOCITY, test_nodal_historical_method_array,
        check_nodal_historical_method_array);

    const auto& test_nodal_non_historical_method_array =
        SpatialMethods::ContainerSpatialMethods<nodes_container_type, node_type,
                                                MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>::CalculateNormVariance<variable_type_array>;
    const auto& check_nodal_non_historical_method_array =
        CheckNormVarianceMethod<nodes_container_type, node_type, variable_type_array, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>;
    RunSpatialNormMethodTest<nodes_container_type, node_type, variable_type_array,
                             std::tuple<double, double>, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "component_x", VELOCITY, test_nodal_non_historical_method_array,
        check_nodal_non_historical_method_array);
    RunSpatialNormMethodTest<nodes_container_type, node_type, variable_type_array,
                             std::tuple<double, double>, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "component_y", VELOCITY, test_nodal_non_historical_method_array,
        check_nodal_non_historical_method_array);
    RunSpatialNormMethodTest<nodes_container_type, node_type, variable_type_array,
                             std::tuple<double, double>, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "component_z", VELOCITY, test_nodal_non_historical_method_array,
        check_nodal_non_historical_method_array);
    RunSpatialNormMethodTest<nodes_container_type, node_type, variable_type_array,
                             std::tuple<double, double>, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "magnitude", VELOCITY, test_nodal_non_historical_method_array,
        check_nodal_non_historical_method_array);

    const auto& test_condition_non_historical_method_array =
        SpatialMethods::ContainerSpatialMethods<conditions_container_type, condition_type,
                                                MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>::CalculateNormVariance<variable_type_array>;
    const auto& check_condition_non_historical_method_array =
        CheckNormVarianceMethod<conditions_container_type, condition_type, variable_type_array, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>;
    RunSpatialNormMethodTest<conditions_container_type, condition_type, variable_type_array,
                             std::tuple<double, double>, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "component_x", VELOCITY, test_condition_non_historical_method_array,
        check_condition_non_historical_method_array);
    RunSpatialNormMethodTest<conditions_container_type, condition_type, variable_type_array,
                             std::tuple<double, double>, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "component_y", VELOCITY, test_condition_non_historical_method_array,
        check_condition_non_historical_method_array);
    RunSpatialNormMethodTest<conditions_container_type, condition_type, variable_type_array,
                             std::tuple<double, double>, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "component_z", VELOCITY, test_condition_non_historical_method_array,
        check_condition_non_historical_method_array);
    RunSpatialNormMethodTest<conditions_container_type, condition_type, variable_type_array,
                             std::tuple<double, double>, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "magnitude", VELOCITY, test_condition_non_historical_method_array,
        check_condition_non_historical_method_array);

    const auto& test_element_non_historical_method_array =
        SpatialMethods::ContainerSpatialMethods<elements_container_type, element_type,
                                                MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>::CalculateNormVariance<variable_type_array>;
    const auto& check_element_non_historical_method_array =
        CheckNormVarianceMethod<elements_container_type, element_type, variable_type_array, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>;
    RunSpatialNormMethodTest<elements_container_type, element_type, variable_type_array,
                             std::tuple<double, double>, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "component_x", VELOCITY, test_element_non_historical_method_array,
        check_element_non_historical_method_array);
    RunSpatialNormMethodTest<elements_container_type, element_type, variable_type_array,
                             std::tuple<double, double>, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "component_y", VELOCITY, test_element_non_historical_method_array,
        check_element_non_historical_method_array);
    RunSpatialNormMethodTest<elements_container_type, element_type, variable_type_array,
                             std::tuple<double, double>, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "component_z", VELOCITY, test_element_non_historical_method_array,
        check_element_non_historical_method_array);
    RunSpatialNormMethodTest<elements_container_type, element_type, variable_type_array,
                             std::tuple<double, double>, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "magnitude", VELOCITY, test_element_non_historical_method_array,
        check_element_non_historical_method_array);
}

KRATOS_TEST_CASE_IN_SUITE(SpatialMinMethod, KratosStatisticsFastSuite)
{
    RunSpatialMinMethodTest<ModelPart::NodesContainerType, ModelPart::NodeType, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "magnitude");
    RunSpatialMinMethodTest<ModelPart::ConditionsContainerType, ModelPart::ConditionType, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "magnitude");
    RunSpatialMinMethodTest<ModelPart::ElementsContainerType, ModelPart::ElementType, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "magnitude");
    RunSpatialMinMethodTest<ModelPart::NodesContainerType, ModelPart::NodeType, MethodsUtilities::HistoricalDataValueRetrievalFunctor>(
        "magnitude");

    RunSpatialMinMethodTest<ModelPart::NodesContainerType, ModelPart::NodeType, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "component_x");
    RunSpatialMinMethodTest<ModelPart::ConditionsContainerType, ModelPart::ConditionType, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "component_x");
    RunSpatialMinMethodTest<ModelPart::ElementsContainerType, ModelPart::ElementType, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "component_x");
    RunSpatialMinMethodTest<ModelPart::NodesContainerType, ModelPart::NodeType, MethodsUtilities::HistoricalDataValueRetrievalFunctor>(
        "component_x");

    RunSpatialMinMethodTest<ModelPart::NodesContainerType, ModelPart::NodeType, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "component_y");
    RunSpatialMinMethodTest<ModelPart::ConditionsContainerType, ModelPart::ConditionType, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "component_y");
    RunSpatialMinMethodTest<ModelPart::ElementsContainerType, ModelPart::ElementType, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "component_y");
    RunSpatialMinMethodTest<ModelPart::NodesContainerType, ModelPart::NodeType, MethodsUtilities::HistoricalDataValueRetrievalFunctor>(
        "component_y");

    RunSpatialMinMethodTest<ModelPart::NodesContainerType, ModelPart::NodeType, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "component_z");
    RunSpatialMinMethodTest<ModelPart::ConditionsContainerType, ModelPart::ConditionType, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "component_z");
    RunSpatialMinMethodTest<ModelPart::ElementsContainerType, ModelPart::ElementType, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "component_z");
    RunSpatialMinMethodTest<ModelPart::NodesContainerType, ModelPart::NodeType, MethodsUtilities::HistoricalDataValueRetrievalFunctor>(
        "component_z");
}

KRATOS_TEST_CASE_IN_SUITE(SpatialMaxMethod, KratosStatisticsFastSuite)
{
    RunSpatialMaxMethodTest<ModelPart::NodesContainerType, ModelPart::NodeType, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "value");
    RunSpatialMaxMethodTest<ModelPart::ConditionsContainerType, ModelPart::ConditionType, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "value");
    RunSpatialMaxMethodTest<ModelPart::ElementsContainerType, ModelPart::ElementType, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "value");
    RunSpatialMaxMethodTest<ModelPart::NodesContainerType, ModelPart::NodeType, MethodsUtilities::HistoricalDataValueRetrievalFunctor>(
        "value");

    RunSpatialMaxMethodTest<ModelPart::NodesContainerType, ModelPart::NodeType, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "magnitude");
    RunSpatialMaxMethodTest<ModelPart::ConditionsContainerType, ModelPart::ConditionType, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "magnitude");
    RunSpatialMaxMethodTest<ModelPart::ElementsContainerType, ModelPart::ElementType, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "magnitude");
    RunSpatialMaxMethodTest<ModelPart::NodesContainerType, ModelPart::NodeType, MethodsUtilities::HistoricalDataValueRetrievalFunctor>(
        "magnitude");

    RunSpatialMaxMethodTest<ModelPart::NodesContainerType, ModelPart::NodeType, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "component_x");
    RunSpatialMaxMethodTest<ModelPart::ConditionsContainerType, ModelPart::ConditionType, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "component_x");
    RunSpatialMaxMethodTest<ModelPart::ElementsContainerType, ModelPart::ElementType, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "component_x");
    RunSpatialMaxMethodTest<ModelPart::NodesContainerType, ModelPart::NodeType, MethodsUtilities::HistoricalDataValueRetrievalFunctor>(
        "component_x");

    RunSpatialMaxMethodTest<ModelPart::NodesContainerType, ModelPart::NodeType, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "component_y");
    RunSpatialMaxMethodTest<ModelPart::ConditionsContainerType, ModelPart::ConditionType, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "component_y");
    RunSpatialMaxMethodTest<ModelPart::ElementsContainerType, ModelPart::ElementType, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "component_y");
    RunSpatialMaxMethodTest<ModelPart::NodesContainerType, ModelPart::NodeType, MethodsUtilities::HistoricalDataValueRetrievalFunctor>(
        "component_y");

    RunSpatialMaxMethodTest<ModelPart::NodesContainerType, ModelPart::NodeType, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "component_z");
    RunSpatialMaxMethodTest<ModelPart::ConditionsContainerType, ModelPart::ConditionType, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "component_z");
    RunSpatialMaxMethodTest<ModelPart::ElementsContainerType, ModelPart::ElementType, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>(
        "component_z");
    RunSpatialMaxMethodTest<ModelPart::NodesContainerType, ModelPart::NodeType, MethodsUtilities::HistoricalDataValueRetrievalFunctor>(
        "component_z");
}

} // namespace Testing
} // namespace Kratos