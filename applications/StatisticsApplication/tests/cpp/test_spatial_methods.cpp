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

    double sum_method_density;
    array_1d<double, 3> sum_method_velocity;
    statistics_methods::CalculateSum(sum_method_density, test_model_part, DENSITY);
    statistics_methods::CalculateSum(sum_method_velocity, test_model_part, VELOCITY);

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

    double mean_method_density;
    array_1d<double, 3> mean_method_velocity;
    statistics_methods::CalculateMean(mean_method_density, test_model_part, DENSITY);
    statistics_methods::CalculateMean(mean_method_velocity, test_model_part, VELOCITY);

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

    double mean_method_density, variance_method_density;
    array_1d<double, 3> mean_method_velocity, variance_method_velocity;
    statistics_methods::CalculateVariance(
        mean_method_density, variance_method_density, test_model_part, DENSITY);
    statistics_methods::CalculateVariance(
        mean_method_velocity, variance_method_velocity, test_model_part, VELOCITY);

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

    double min_method_density, min_method_velocity;
    std::size_t min_method_density_id, min_method_velocity_id;
    statistics_methods::GetMin(min_method_density, min_method_density_id,
                               rNormType, test_model_part, DENSITY);
    statistics_methods::GetMin(min_method_velocity, min_method_velocity_id,
                               rNormType, test_model_part, VELOCITY);

    const double max_value = std::numeric_limits<double>::max();

    double min_density{max_value}, min_velocity{max_value};
    std::size_t min_density_id{0}, min_velocity_id{0};

    for (const TContainerItemType& r_item : r_container)
    {
        const double current_density =
            TDataRetrievalFunctor<TContainerItemType>()(r_item, DENSITY);
        if (current_density < min_density)
        {
            min_density = current_density;
            min_density_id = r_item.Id();
        }

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

    KRATOS_CHECK_NEAR(min_density, min_method_density, 1e-12);
    KRATOS_CHECK_EQUAL(min_density_id, min_method_density_id);
    KRATOS_CHECK_NEAR(min_velocity, min_method_velocity, 1e-12);
    KRATOS_CHECK_EQUAL(min_velocity_id, min_method_velocity_id);
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

    double max_method_density, max_method_velocity;
    std::size_t max_method_density_id, max_method_velocity_id;
    statistics_methods::GetMax(max_method_density, max_method_density_id,
                               rNormType, test_model_part, DENSITY);
    statistics_methods::GetMax(max_method_velocity, max_method_velocity_id,
                               rNormType, test_model_part, VELOCITY);

    const double min_value = std::numeric_limits<double>::lowest();

    double max_density{min_value}, max_velocity{min_value};
    std::size_t max_density_id{0}, max_velocity_id{0};

    for (const TContainerItemType& r_item : r_container)
    {
        const double current_density =
            TDataRetrievalFunctor<TContainerItemType>()(r_item, DENSITY);
        if (current_density > max_density)
        {
            max_density = current_density;
            max_density_id = r_item.Id();
        }

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

    KRATOS_CHECK_NEAR(max_density, max_method_density, 1e-12);
    KRATOS_CHECK_EQUAL(max_density_id, max_method_density_id);
    KRATOS_CHECK_NEAR(max_velocity, max_method_velocity, 1e-12);
    KRATOS_CHECK_EQUAL(max_velocity_id, max_method_velocity_id);
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