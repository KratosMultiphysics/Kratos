// Project includes
#include "testing/testing.h"
#include "containers/multi_level_data_value_container.h"


namespace Kratos::Testing 
{

class TestAccesor
{
public:
    TestAccesor(std::size_t num_layers, std::size_t num_gauss_points)
        : mNumLayers(num_layers), mNumGaussPoints(num_gauss_points)
    {
    }

    std::size_t GetIndex(std::vector<std::size_t> const& index) const
    {
        return (index[0] * mNumGaussPoints) + index[1];
    }

    std::size_t Size() const
    {
        return mNumLayers * mNumGaussPoints;
    }

    std::size_t mNumLayers;
    std::size_t mNumGaussPoints;
};

KRATOS_TEST_CASE_IN_SUITE(MultiLevelDataValueContainerSetValue, KratosCoreFastSuite)
{
    int dummy;
    MultiLevelDataValueContainer container;
    TestAccesor accesor(3, 4);
    KRATOS_WATCH("Im starting")
    container.SetValue(TEMPERATURE, accesor, {1,2}, 2.0);
    double value = container.GetValue(TEMPERATURE, accesor, {1,2});
    KRATOS_EXPECT_NEAR(value, 2.0, 1e-12);

    container.SetValue(PRESSURE, accesor, {2,1}, 3.0);
    value = container.GetValue(PRESSURE, accesor, {2,1});
    KRATOS_EXPECT_NEAR(value, 3.0, 1e-12);

    array_1d<double, 3> array = ZeroVector(3);
    array[0] = 1.0; array[1] = 2.0; array[2] = 3.0;
    container.SetValue(VELOCITY, accesor, {1,2}, array);
    KRATOS_WATCH(container.GetValue(VELOCITY, accesor, {1,2}));
    KRATOS_WATCH("Im fine")
}

}

