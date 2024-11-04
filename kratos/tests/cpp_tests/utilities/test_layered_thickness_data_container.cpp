#include "utilities/layered_thickness_data_container.h"
#include "testing/testing.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(LayeredThicknessDataContainerConstructorWithVector, KratosCoreFastSuite)
{
    std::vector<DataValueContainer> data_vector(2);
    data_vector[0].SetValue(TEMPERATURE, 100.0);
    data_vector[1].SetValue(TEMPERATURE, 200.0);

    LayeredThicknessDataContainer container(data_vector);

    KRATOS_EXPECT_EQ(container.Size(), 2);
    KRATOS_EXPECT_NEAR(container.Get(0).GetValue(TEMPERATURE), 100.0, 1e-8);
    KRATOS_EXPECT_NEAR(container.Get(1).GetValue(TEMPERATURE), 200.0, 1e-8);
}

KRATOS_TEST_CASE_IN_SUITE(LayeredThicknessDataContainerAddGetSet, KratosCoreFastSuite)
{
    LayeredThicknessDataContainer container;
    DataValueContainer data1, data2;
    data1.SetValue(TEMPERATURE, 100.0);
    data2.SetValue(TEMPERATURE, 200.0);

    container.Add(data1);
    container.Add(data2);

    KRATOS_EXPECT_EQ(container.Size(), 2);
    const DataValueContainer& data0_ref = container.Get(0);
    DataValueContainer& data1_ref = container.Get(1);
    KRATOS_EXPECT_NEAR(data0_ref.GetValue(TEMPERATURE), 100.0, 1e-8);
    KRATOS_EXPECT_NEAR(data1_ref.GetValue(TEMPERATURE), 200.0, 1e-8);

    DataValueContainer data3;
    data3.SetValue(TEMPERATURE, 300.0);
    container.Set(1, data3);

    KRATOS_EXPECT_NEAR(container.Get(1).GetValue(TEMPERATURE), 300.0, 1e-8);
}

KRATOS_TEST_CASE_IN_SUITE(LayeredThicknessDataContainerSetOutOfRange, KratosCoreFastSuite)
{
    LayeredThicknessDataContainer container;
    DataValueContainer data;
    data.SetValue(TEMPERATURE, 150.0);

    container.Add(data);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(container.Set(1, data), "Index out of range");
}

}
