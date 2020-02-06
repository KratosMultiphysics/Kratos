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

#if !defined(KRATOS_STATISTICS_TEST_UTILITIES_H_INCLUDED)
#define KRATOS_STATISTICS_TEST_UTILITIES_H_INCLUDED

// System includes
#include <random>
#include <string>

// External includes

// Project includes

// Application includes

namespace Kratos
{
namespace StatisticsApplicationTestUtilities
{
template <typename TContainerType, typename TContainerItemType, template <typename T> typename TDataRetrievalFunctor>
class RandomInitializer
{
public:
    void static InitializeVariableWithRandomValues(TContainerType& rContainer,
                                                   const Variable<double>& rVariable)
    {
        std::string seed_str =
            "KratosStatisticsApplicationTestSeed_" + rVariable.Name();
        std::seed_seq seed(seed_str.begin(), seed_str.end());
        std::default_random_engine generator(seed);
        std::uniform_real_distribution<double> distribution(-10.0, 10.0);

        for (TContainerItemType& container_item : rContainer)
        {
            TDataRetrievalFunctor<TContainerItemType>()(container_item, rVariable) =
                distribution(generator);
        }
    }

    void static InitializeVariableWithRandomValues(TContainerType& rContainer,
                                                   const Variable<array_1d<double, 3>>& rVariable)
    {
        std::string seed_str =
            "KratosStatisticsApplicationTestSeed_" + rVariable.Name();
        std::seed_seq seed(seed_str.begin(), seed_str.end());
        std::default_random_engine generator(seed);

        for (int i = 0; i < 3; ++i)
        {
            std::uniform_real_distribution<double> distribution(-10.0, 10.0);

            for (TContainerItemType& container_item : rContainer)
            {
                TDataRetrievalFunctor<TContainerItemType>()(
                    container_item, rVariable)[i] = distribution(generator);
            }
        }
    }
};

void AddNodalSolutionStepVariables(ModelPart& rModelPart);

void CreateModelPart(ModelPart& rModelPart);

} // namespace StatisticsApplicationTestUtilities
} // namespace Kratos

#endif // KRATOS_STATISTICS_TEST_UTILITIES_H_INCLUDED