// Project includes
#include "testing/testing.h"
#include "containers/multi_level_data_value_container.h"
#include "utilities/multi_level_data_accesors.h"


namespace Kratos::Testing 
{

KRATOS_TEST_CASE_IN_SUITE(MultiLevelDataValueContainerCopyConstructor, KratosCoreFastSuite)
{
    MultiLevelDataValueContainer container;
    auto p_accesor = Kratos::make_unique<LayeredGaussPointDataAccesor>(3, 4);
    container.AddVariable(TEMPERATURE, std::move(p_accesor));
    container.SetValue(TEMPERATURE, 3.0, 1, 2);

    MultiLevelDataValueContainer container_copy(container);
    KRATOS_EXPECT_TRUE(container_copy.Has(TEMPERATURE));
    KRATOS_EXPECT_NEAR(container_copy.GetValue(TEMPERATURE, 1, 2), 3.0, 1e-12);
    KRATOS_EXPECT_TRUE(container.Has(TEMPERATURE));
    KRATOS_EXPECT_NEAR(container.GetValue(TEMPERATURE, 1, 2), 3.0, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(MultiLevelDataValueContainerMoveConstructor, KratosCoreFastSuite)
{
    MultiLevelDataValueContainer container;
    auto p_accesor = Kratos::make_unique<LayeredGaussPointDataAccesor>(3, 4);
    container.AddVariable(TEMPERATURE, std::move(p_accesor));
    container.SetValue(TEMPERATURE, 3.0, 1, 2);

    MultiLevelDataValueContainer container_copy(std::move(container));
    KRATOS_EXPECT_TRUE(container_copy.Has(TEMPERATURE));
    KRATOS_EXPECT_NEAR(container_copy.GetValue(TEMPERATURE, 1, 2), 3.0, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(MultiLevelDataValueContainerCopyAssignment, KratosCoreFastSuite)
{
    MultiLevelDataValueContainer container;
    auto p_accesor = Kratos::make_unique<LayeredGaussPointDataAccesor>(3, 4);
    container.AddVariable(TEMPERATURE, std::move(p_accesor));
    container.SetValue(TEMPERATURE, 3.0, 1, 2);

    MultiLevelDataValueContainer container_copy;
    container_copy = container;
    KRATOS_EXPECT_TRUE(container_copy.Has(TEMPERATURE));
    KRATOS_EXPECT_NEAR(container_copy.GetValue(TEMPERATURE, 1, 2), 3.0, 1e-12);
    KRATOS_EXPECT_TRUE(container.Has(TEMPERATURE));
    KRATOS_EXPECT_NEAR(container.GetValue(TEMPERATURE, 1, 2), 3.0, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(MultiLevelDataValueContainerMoveAssignment, KratosCoreFastSuite)
{
    MultiLevelDataValueContainer container;
    auto p_accesor = Kratos::make_unique<LayeredGaussPointDataAccesor>(3, 4);
    container.AddVariable(TEMPERATURE, std::move(p_accesor));
    container.SetValue(TEMPERATURE, 3.0, 1, 2);

    MultiLevelDataValueContainer container_copy;
    container_copy = std::move(container);
    KRATOS_EXPECT_TRUE(container_copy.Has(TEMPERATURE));
    KRATOS_EXPECT_NEAR(container_copy.GetValue(TEMPERATURE, 1, 2), 3.0, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(MultiLevelDataValueContainerHas, KratosCoreFastSuite)
{
    MultiLevelDataValueContainer container;
    auto p_accesor = Kratos::make_unique<LayeredGaussPointDataAccesor>(3, 4);
    KRATOS_EXPECT_FALSE(container.Has(TEMPERATURE));
    container.AddVariable(TEMPERATURE, std::move(p_accesor));
    container.SetValue(TEMPERATURE, 2.0, 1, 2);
    KRATOS_EXPECT_TRUE(container.Has(TEMPERATURE));
}

KRATOS_TEST_CASE_IN_SUITE(MultiLevelDataValueContainerRemove, KratosCoreFastSuite)
{
    MultiLevelDataValueContainer container;
    auto p_accesor = Kratos::make_unique<LayeredGaussPointDataAccesor>(3, 4);
    auto p_accesor_2 = Kratos::make_unique<LayeredGaussPointDataAccesor>(3, 4);
    container.AddVariable(TEMPERATURE, std::move(p_accesor));
    container.AddVariable(PRESSURE, std::move(p_accesor_2));
    KRATOS_EXPECT_TRUE(container.Has(TEMPERATURE));
    KRATOS_EXPECT_TRUE(container.Has(PRESSURE));
    container.RemoveVariable(TEMPERATURE);
    KRATOS_EXPECT_FALSE(container.Has(TEMPERATURE));
    KRATOS_EXPECT_TRUE(container.Has(PRESSURE));
}

KRATOS_TEST_CASE_IN_SUITE(MultiLevelDataValueContainerSetValueDouble, KratosCoreFastSuite)
{
    MultiLevelDataValueContainer container;
    auto p_accesor = Kratos::make_unique<LayeredGaussPointDataAccesor>(3, 4);
    container.AddVariable(TEMPERATURE, std::move(p_accesor));
    container.SetValue(TEMPERATURE, 2.0, 1, 2);
    double value = container.GetValue(TEMPERATURE, 1, 2);
    KRATOS_EXPECT_NEAR(value, 2.0, 1e-12);

    auto p_accesor_2 = Kratos::make_unique<LayeredGaussPointDataAccesor>(3, 4);
    container.AddVariable(PRESSURE, std::move(p_accesor_2));
    container.SetValue(PRESSURE, 3.0, 2, 1);
    value = container.GetValue(PRESSURE, 2, 1);
    KRATOS_EXPECT_NEAR(value, 3.0, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(MultiLevelDataValueContainerSetValueArray, KratosCoreFastSuite)
{
    MultiLevelDataValueContainer container;
    auto p_accesor = Kratos::make_unique<LayeredGaussPointDataAccesor>(3, 4);
    array_1d<double, 3> array = ZeroVector(3);
    array[0] = 1.0; array[1] = 2.0; array[2] = 3.0;
    container.AddVariable(VELOCITY, std::move(p_accesor));
    container.SetValue(VELOCITY, array, 1, 2);
    KRATOS_EXPECT_VECTOR_NEAR(container.GetValue(VELOCITY, 1, 2), array, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(MultiLevelDataValueContainerSetValueVector, KratosCoreFastSuite)
{
    MultiLevelDataValueContainer container;
    auto p_accesor = Kratos::make_unique<LayeredGaussPointDataAccesor>(3, 4);
    Vector vector = ZeroVector(3);
    vector[0] = 1.5; vector[1] = 2.5; vector[2] = 3.5;
    container.AddVariable(STRESSES, std::move(p_accesor));
    container.SetValue(STRESSES, vector, 0, 1);
    KRATOS_EXPECT_VECTOR_NEAR(container.GetValue(STRESSES, 0, 1), vector, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(MultiLevelDataValueContainerSetValueMatrix, KratosCoreFastSuite)
{
    MultiLevelDataValueContainer container;
    auto p_accesor = Kratos::make_unique<LayeredGaussPointDataAccesor>(3, 4);
    Matrix matrix = ZeroMatrix(3, 3);
    matrix(0,0) = 5.7;
    container.AddVariable(CAUCHY_STRESS_TENSOR, std::move(p_accesor));
    container.SetValue(CAUCHY_STRESS_TENSOR, matrix, 0, 3);
    KRATOS_EXPECT_MATRIX_NEAR(container.GetValue(CAUCHY_STRESS_TENSOR, 0, 3), matrix, 1e-12);
}


KRATOS_TEST_CASE_IN_SUITE(MultiLevelDataValueContainerSerializer, KratosCoreFastSuite)
{
    MultiLevelDataValueContainer container;
    auto p_accesor = Kratos::make_unique<LayeredGaussPointDataAccesor>(3, 4);
    auto p_accesor_2 = Kratos::make_unique<LayeredGaussPointDataAccesor>(3, 4);
    container.AddVariable(TEMPERATURE, std::move(p_accesor));
    container.AddVariable(VELOCITY, std::move(p_accesor_2));
    container.SetValue(TEMPERATURE, 2.0, 0, 0);
    container.SetValue(TEMPERATURE, 3.0, 1, 1);
    array_1d<double,3> velocity({1.0, 2.0, 3.0});
    container.SetValue(VELOCITY, velocity, 2, 2);

    StreamSerializer serializer;
    serializer.save("Container", container);

    MultiLevelDataValueContainer container_loaded;
    serializer.load("Container", container_loaded);

    KRATOS_EXPECT_TRUE(container_loaded.Has(TEMPERATURE));
    KRATOS_EXPECT_NEAR(container_loaded.GetValue(TEMPERATURE, 0, 0), 2.0, 1e-12);
    KRATOS_EXPECT_NEAR(container_loaded.GetValue(TEMPERATURE, 1, 1), 3.0, 1e-12);
    KRATOS_EXPECT_VECTOR_NEAR(container_loaded.GetValue(VELOCITY, 2, 2), velocity, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(MultiLevelDataValueContainerComponent, KratosCoreFastSuite)
{
    MultiLevelDataValueContainer container;
    auto p_accesor = Kratos::make_unique<LayeredGaussPointDataAccesor>(3, 4);
    container.AddVariable(VELOCITY, std::move(p_accesor));
    array_1d<double, 3> array = ZeroVector(3);
    array[0] = 2.5; array[1] = 3.5; array[2] = 4.5;
    container.SetValue(VELOCITY_X, array[0], 1, 2);
    container.SetValue(VELOCITY_Y, array[1], 1, 2);
    container.SetValue(VELOCITY_Z, array[2], 1, 2);
    KRATOS_EXPECT_NEAR(container.GetValue(VELOCITY_X, 1, 2), array[0], 1e-12);
    KRATOS_EXPECT_NEAR(container.GetValue(VELOCITY_Y, 1, 2), array[1], 1e-12);
    KRATOS_EXPECT_NEAR(container.GetValue(VELOCITY_Z, 1, 2), array[2], 1e-12);
    KRATOS_EXPECT_VECTOR_NEAR(container.GetValue(VELOCITY, 1, 2), array, 1e-12);
}

} // namespace Kratos::Testing
