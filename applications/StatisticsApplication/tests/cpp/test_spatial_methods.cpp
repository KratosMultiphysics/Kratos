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
template <typename TContainerType, typename TContainerItemType>
void InitializeSpatialMethodNonHistoricalVariables(TContainerType& rContainer)
{
    using random_initializer =
        StatisticsApplicationTestUtilities::RandomInitializer<TContainerType, TContainerItemType, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor,
                                                              MethodsUtilities::NonHistoricalDataValueInitializationFunctor>;
    random_initializer::InitializeVariableWithRandomValues(rContainer, DENSITY);
    random_initializer::InitializeVariableWithRandomValues(rContainer, VELOCITY);
    random_initializer::InitializeVariableWithRandomValues(rContainer, EXTERNAL_FORCES_VECTOR);
    random_initializer::InitializeVariableWithRandomValues(
        rContainer, GREEN_LAGRANGE_STRAIN_TENSOR);
}

template <typename TContainerType, typename TContainerItemType>
void InitializeSpatialMethodHistoricalVariables(TContainerType& rContainer)
{
    using random_initializer =
        StatisticsApplicationTestUtilities::RandomInitializer<TContainerType, TContainerItemType, MethodsUtilities::HistoricalDataValueRetrievalFunctor,
                                                              MethodsUtilities::HistoricalDataValueInitializationFunctor>;
    random_initializer::InitializeVariableWithRandomValues(rContainer, DENSITY);
    random_initializer::InitializeVariableWithRandomValues(rContainer, VELOCITY);
    random_initializer::InitializeVariableWithRandomValues(rContainer, EXTERNAL_FORCES_VECTOR);
    random_initializer::InitializeVariableWithRandomValues(
        rContainer, GREEN_LAGRANGE_STRAIN_TENSOR);
}
} // namespace

KRATOS_TEST_CASE_IN_SUITE(SpatialSumNodalNonHistorical, KratosStatisticsFastSuite)
{
    Model test_model;
    ModelPart& test_model_part = test_model.CreateModelPart("test_model_part");
    StatisticsApplicationTestUtilities::CreateModelPart(test_model_part);

    using container_item_type = ModelPart::NodeType;
    using container_type = ModelPart::NodesContainerType;

    container_type& r_container = test_model_part.Nodes();

    InitializeSpatialMethodNonHistoricalVariables<container_type, container_item_type>(r_container);
    using statistics_methods =
        SpatialMethods::ContainerSpatialMethods<container_type, container_item_type, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>;

    double sum_method_density;
    array_1d<double, 3> sum_method_velocity;
    Vector sum_method_external_forces_vector;
    Matrix sum_method_green_lagrange_strain_tensor;
    statistics_methods::CalculateSum(sum_method_density, r_container, DENSITY);
    statistics_methods::CalculateSum(sum_method_velocity, r_container, VELOCITY);
    statistics_methods::CalculateSum(sum_method_external_forces_vector,
                                     r_container, EXTERNAL_FORCES_VECTOR);
    statistics_methods::CalculateSum(sum_method_green_lagrange_strain_tensor,
                                     r_container, GREEN_LAGRANGE_STRAIN_TENSOR);

    double sum_density = 0.0;
    array_1d<double, 3> sum_velocity = ZeroVector(3);
    Vector sum_external_forces_vector(5, 0.0);
    Matrix sum_green_lagrange_strain_tensor(5, 5, 0.0);
    for (const container_item_type& r_item : r_container)
    {
        sum_density += r_item.GetValue(DENSITY);
        sum_velocity += r_item.GetValue(VELOCITY);
        sum_external_forces_vector += r_item.GetValue(EXTERNAL_FORCES_VECTOR);
        sum_green_lagrange_strain_tensor += r_item.GetValue(GREEN_LAGRANGE_STRAIN_TENSOR);
    }

    KRATOS_CHECK_NEAR(sum_density, sum_method_density, 1e-14);
    KRATOS_CHECK_VECTOR_NEAR(sum_velocity, sum_method_velocity, 1e-14);
    KRATOS_CHECK_VECTOR_NEAR(sum_external_forces_vector,
                             sum_method_external_forces_vector, 1e-14);
    KRATOS_CHECK_MATRIX_NEAR(sum_green_lagrange_strain_tensor,
                             sum_method_green_lagrange_strain_tensor, 1e-14);
}

KRATOS_TEST_CASE_IN_SUITE(SpatialSumConditionNonHistorical, KratosStatisticsFastSuite)
{
    Model test_model;
    ModelPart& test_model_part = test_model.CreateModelPart("test_model_part");
    StatisticsApplicationTestUtilities::CreateModelPart(test_model_part);

    using container_item_type = ModelPart::ConditionType;
    using container_type = ModelPart::ConditionsContainerType;

    container_type& r_container = test_model_part.Conditions();

    InitializeSpatialMethodNonHistoricalVariables<container_type, container_item_type>(r_container);
    using statistics_methods =
        SpatialMethods::ContainerSpatialMethods<container_type, container_item_type, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>;

    double sum_method_density;
    array_1d<double, 3> sum_method_velocity;
    Vector sum_method_external_forces_vector;
    Matrix sum_method_green_lagrange_strain_tensor;
    statistics_methods::CalculateSum(sum_method_density, r_container, DENSITY);
    statistics_methods::CalculateSum(sum_method_velocity, r_container, VELOCITY);
    statistics_methods::CalculateSum(sum_method_external_forces_vector,
                                     r_container, EXTERNAL_FORCES_VECTOR);
    statistics_methods::CalculateSum(sum_method_green_lagrange_strain_tensor,
                                     r_container, GREEN_LAGRANGE_STRAIN_TENSOR);

    double sum_density = 0.0;
    array_1d<double, 3> sum_velocity = ZeroVector(3);
    Vector sum_external_forces_vector(5, 0.0);
    Matrix sum_green_lagrange_strain_tensor(5, 5, 0.0);
    for (const container_item_type& r_item : r_container)
    {
        sum_density += r_item.GetValue(DENSITY);
        sum_velocity += r_item.GetValue(VELOCITY);
        sum_external_forces_vector += r_item.GetValue(EXTERNAL_FORCES_VECTOR);
        sum_green_lagrange_strain_tensor += r_item.GetValue(GREEN_LAGRANGE_STRAIN_TENSOR);
    }

    KRATOS_CHECK_NEAR(sum_density, sum_method_density, 1e-14);
    KRATOS_CHECK_VECTOR_NEAR(sum_velocity, sum_method_velocity, 1e-14);
    KRATOS_CHECK_VECTOR_NEAR(sum_external_forces_vector,
                             sum_method_external_forces_vector, 1e-14);
    KRATOS_CHECK_MATRIX_NEAR(sum_green_lagrange_strain_tensor,
                             sum_method_green_lagrange_strain_tensor, 1e-14);
}

KRATOS_TEST_CASE_IN_SUITE(SpatialSumElementNonHistorical, KratosStatisticsFastSuite)
{
    Model test_model;
    ModelPart& test_model_part = test_model.CreateModelPart("test_model_part");
    StatisticsApplicationTestUtilities::CreateModelPart(test_model_part);

    using container_item_type = ModelPart::ElementType;
    using container_type = ModelPart::ElementsContainerType;

    container_type& r_container = test_model_part.Elements();

    InitializeSpatialMethodNonHistoricalVariables<container_type, container_item_type>(r_container);
    using statistics_methods =
        SpatialMethods::ContainerSpatialMethods<container_type, container_item_type, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>;

    double sum_method_density;
    array_1d<double, 3> sum_method_velocity;
    Vector sum_method_external_forces_vector;
    Matrix sum_method_green_lagrange_strain_tensor;
    statistics_methods::CalculateSum(sum_method_density, r_container, DENSITY);
    statistics_methods::CalculateSum(sum_method_velocity, r_container, VELOCITY);
    statistics_methods::CalculateSum(sum_method_external_forces_vector,
                                     r_container, EXTERNAL_FORCES_VECTOR);
    statistics_methods::CalculateSum(sum_method_green_lagrange_strain_tensor,
                                     r_container, GREEN_LAGRANGE_STRAIN_TENSOR);

    double sum_density = 0.0;
    array_1d<double, 3> sum_velocity = ZeroVector(3);
    Vector sum_external_forces_vector(5, 0.0);
    Matrix sum_green_lagrange_strain_tensor(5, 5, 0.0);
    for (const container_item_type& r_item : r_container)
    {
        sum_density += r_item.GetValue(DENSITY);
        sum_velocity += r_item.GetValue(VELOCITY);
        sum_external_forces_vector += r_item.GetValue(EXTERNAL_FORCES_VECTOR);
        sum_green_lagrange_strain_tensor += r_item.GetValue(GREEN_LAGRANGE_STRAIN_TENSOR);
    }

    KRATOS_CHECK_NEAR(sum_density, sum_method_density, 1e-14);
    KRATOS_CHECK_VECTOR_NEAR(sum_velocity, sum_method_velocity, 1e-14);
    KRATOS_CHECK_VECTOR_NEAR(sum_external_forces_vector,
                             sum_method_external_forces_vector, 1e-14);
    KRATOS_CHECK_MATRIX_NEAR(sum_green_lagrange_strain_tensor,
                             sum_method_green_lagrange_strain_tensor, 1e-14);
}

KRATOS_TEST_CASE_IN_SUITE(SpatialMeanNodalNonHistorical, KratosStatisticsFastSuite)
{
    Model test_model;
    ModelPart& test_model_part = test_model.CreateModelPart("test_model_part");
    StatisticsApplicationTestUtilities::CreateModelPart(test_model_part);

    using container_item_type = ModelPart::NodeType;
    using container_type = ModelPart::NodesContainerType;

    container_type& r_container = test_model_part.Nodes();

    InitializeSpatialMethodNonHistoricalVariables<container_type, container_item_type>(r_container);
    using statistics_methods =
        SpatialMethods::ContainerSpatialMethods<container_type, container_item_type, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>;

    double mean_method_density;
    array_1d<double, 3> mean_method_velocity;
    Vector mean_method_external_forces_vector;
    Matrix mean_method_green_lagrange_strain_tensor;
    statistics_methods::CalculateMean(mean_method_density, r_container, DENSITY);
    statistics_methods::CalculateMean(mean_method_velocity, r_container, VELOCITY);
    statistics_methods::CalculateMean(mean_method_external_forces_vector,
                                      r_container, EXTERNAL_FORCES_VECTOR);
    statistics_methods::CalculateMean(mean_method_green_lagrange_strain_tensor,
                                      r_container, GREEN_LAGRANGE_STRAIN_TENSOR);

    double mean_density = 0.0;
    array_1d<double, 3> mean_velocity = ZeroVector(3);
    Vector mean_external_forces_vector(5, 0.0);
    Matrix mean_green_lagrange_strain_tensor(5, 5, 0.0);
    const double n = r_container.size();
    for (const container_item_type& r_item : r_container)
    {
        mean_density += r_item.GetValue(DENSITY) / n;
        mean_velocity += r_item.GetValue(VELOCITY) / n;
        mean_external_forces_vector += r_item.GetValue(EXTERNAL_FORCES_VECTOR) / n;
        mean_green_lagrange_strain_tensor +=
            r_item.GetValue(GREEN_LAGRANGE_STRAIN_TENSOR) / n;
    }

    KRATOS_CHECK_NEAR(mean_density, mean_method_density, 1e-14);
    KRATOS_CHECK_VECTOR_NEAR(mean_velocity, mean_method_velocity, 1e-14);
    KRATOS_CHECK_VECTOR_NEAR(mean_external_forces_vector,
                             mean_method_external_forces_vector, 1e-14);
    KRATOS_CHECK_MATRIX_NEAR(mean_green_lagrange_strain_tensor,
                             mean_method_green_lagrange_strain_tensor, 1e-14);
}

KRATOS_TEST_CASE_IN_SUITE(SpatialMeanConditionNonHistorical, KratosStatisticsFastSuite)
{
    Model test_model;
    ModelPart& test_model_part = test_model.CreateModelPart("test_model_part");
    StatisticsApplicationTestUtilities::CreateModelPart(test_model_part);

    using container_item_type = ModelPart::ConditionType;
    using container_type = ModelPart::ConditionsContainerType;

    container_type& r_container = test_model_part.Conditions();

    InitializeSpatialMethodNonHistoricalVariables<container_type, container_item_type>(r_container);
    using statistics_methods =
        SpatialMethods::ContainerSpatialMethods<container_type, container_item_type, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>;

    double mean_method_density;
    array_1d<double, 3> mean_method_velocity;
    Vector mean_method_external_forces_vector;
    Matrix mean_method_green_lagrange_strain_tensor;
    statistics_methods::CalculateMean(mean_method_density, r_container, DENSITY);
    statistics_methods::CalculateMean(mean_method_velocity, r_container, VELOCITY);
    statistics_methods::CalculateMean(mean_method_external_forces_vector,
                                      r_container, EXTERNAL_FORCES_VECTOR);
    statistics_methods::CalculateMean(mean_method_green_lagrange_strain_tensor,
                                      r_container, GREEN_LAGRANGE_STRAIN_TENSOR);

    double mean_density = 0.0;
    array_1d<double, 3> mean_velocity = ZeroVector(3);
    Vector mean_external_forces_vector(5, 0.0);
    Matrix mean_green_lagrange_strain_tensor(5, 5, 0.0);
    const double n = r_container.size();
    for (const container_item_type& r_item : r_container)
    {
        mean_density += r_item.GetValue(DENSITY) / n;
        mean_velocity += r_item.GetValue(VELOCITY) / n;
        mean_external_forces_vector += r_item.GetValue(EXTERNAL_FORCES_VECTOR) / n;
        mean_green_lagrange_strain_tensor +=
            r_item.GetValue(GREEN_LAGRANGE_STRAIN_TENSOR) / n;
    }

    KRATOS_CHECK_NEAR(mean_density, mean_method_density, 1e-14);
    KRATOS_CHECK_VECTOR_NEAR(mean_velocity, mean_method_velocity, 1e-14);
    KRATOS_CHECK_VECTOR_NEAR(mean_external_forces_vector,
                             mean_method_external_forces_vector, 1e-14);
    KRATOS_CHECK_MATRIX_NEAR(mean_green_lagrange_strain_tensor,
                             mean_method_green_lagrange_strain_tensor, 1e-14);
}

KRATOS_TEST_CASE_IN_SUITE(SpatialMeanElementNonHistorical, KratosStatisticsFastSuite)
{
    Model test_model;
    ModelPart& test_model_part = test_model.CreateModelPart("test_model_part");
    StatisticsApplicationTestUtilities::CreateModelPart(test_model_part);

    using container_item_type = ModelPart::ElementType;
    using container_type = ModelPart::ElementsContainerType;

    container_type& r_container = test_model_part.Elements();

    InitializeSpatialMethodNonHistoricalVariables<container_type, container_item_type>(r_container);
    using statistics_methods =
        SpatialMethods::ContainerSpatialMethods<container_type, container_item_type, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>;

    double mean_method_density;
    array_1d<double, 3> mean_method_velocity;
    Vector mean_method_external_forces_vector;
    Matrix mean_method_green_lagrange_strain_tensor;
    statistics_methods::CalculateMean(mean_method_density, r_container, DENSITY);
    statistics_methods::CalculateMean(mean_method_velocity, r_container, VELOCITY);
    statistics_methods::CalculateMean(mean_method_external_forces_vector,
                                      r_container, EXTERNAL_FORCES_VECTOR);
    statistics_methods::CalculateMean(mean_method_green_lagrange_strain_tensor,
                                      r_container, GREEN_LAGRANGE_STRAIN_TENSOR);

    double mean_density = 0.0;
    array_1d<double, 3> mean_velocity = ZeroVector(3);
    Vector mean_external_forces_vector(5, 0.0);
    Matrix mean_green_lagrange_strain_tensor(5, 5, 0.0);
    const double n = r_container.size();
    for (const container_item_type& r_item : r_container)
    {
        mean_density += r_item.GetValue(DENSITY) / n;
        mean_velocity += r_item.GetValue(VELOCITY) / n;
        mean_external_forces_vector += r_item.GetValue(EXTERNAL_FORCES_VECTOR) / n;
        mean_green_lagrange_strain_tensor +=
            r_item.GetValue(GREEN_LAGRANGE_STRAIN_TENSOR) / n;
    }

    KRATOS_CHECK_NEAR(mean_density, mean_method_density, 1e-14);
    KRATOS_CHECK_VECTOR_NEAR(mean_velocity, mean_method_velocity, 1e-14);
    KRATOS_CHECK_VECTOR_NEAR(mean_external_forces_vector,
                             mean_method_external_forces_vector, 1e-14);
    KRATOS_CHECK_MATRIX_NEAR(mean_green_lagrange_strain_tensor,
                             mean_method_green_lagrange_strain_tensor, 1e-14);
}
} // namespace Testing
} // namespace Kratos