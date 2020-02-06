//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

#if !defined(KRATOS_SPATIAL_VARIANCE_METHOD_H_INCLUDED)
#define KRATOS_SPATIAL_VARIANCE_METHOD_H_INCLUDED

// System includes
#include <cmath>

// External includes

// Project includes
#include "includes/communicator.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"

// Application includes
#include "custom_utilities/method_utilities.h"

namespace Kratos
{
///@addtogroup RANSApplication
///@{

///@name Kratos Globals
///@{

namespace SpatialMethods
{
template <typename TContainerType, typename TContainerItemType, template <typename T> typename TDataRetrievalFunctor>
class ContainerSpatialMethods
{
public:
    template <typename TDataType>
    void static CalculateSum(TDataType& rSum,
                             const ModelPart& rModelPart,
                             const Variable<TDataType>& rVariable)
    {
        KRATOS_TRY

        const TContainerType& r_container =
            MethodsUtilities::GetDataContainer<TContainerType>(rModelPart);

        const TDataType& r_initial_value =
            TDataRetrievalFunctor<TContainerItemType>()(*r_container.begin(), rVariable);

        rSum = rVariable.Zero();
        MethodsUtilities::DataTypeSizeInitializer(rSum, r_initial_value);

#pragma omp parallel
        {
            TDataType sum = rVariable.Zero();
            MethodsUtilities::DataTypeSizeInitializer(sum, r_initial_value);

#pragma omp for
            for (int i = 0; i < static_cast<int>(r_container.size()); ++i)
            {
                const TContainerItemType& r_item = *(r_container.begin() + i);
                const TDataType& current_value =
                    TDataRetrievalFunctor<TContainerItemType>()(r_item, rVariable);
                MethodsUtilities::DataTypeSizeChecker(current_value, sum);
                sum += current_value;
            }
#pragma omp critical
            {
                rSum += sum;
            }
        }

        rSum = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(rSum);

        KRATOS_CATCH("");
    }

    // special overloaded method for flags
    void static CalculateSum(int& rSum, const ModelPart& rModelPart, const Flags& rVariable)
    {
        const TContainerType& r_container =
            MethodsUtilities::GetDataContainer<TContainerType>(rModelPart);
        int sum = 0;

#pragma omp parallel for reduction(+ : sum)
        for (int i = 0; i < static_cast<int>(r_container.size()); ++i)
        {
            const TContainerItemType& r_item = *(r_container.begin() + i);
            sum += r_item.Is(rVariable);
        }

        rSum = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(sum);
    }

    template <typename TDataType>
    void static CalculateMean(TDataType& rMean,
                              const ModelPart& rModelPart,
                              const Variable<TDataType>& rVariable)
    {
        CalculateSum<TDataType>(rMean, rModelPart, rVariable);
        const TContainerType& r_container =
            MethodsUtilities::GetDataContainer<TContainerType>(rModelPart);

        const unsigned int number_of_items =
            rModelPart.GetCommunicator().GetDataCommunicator().SumAll(
                r_container.size());

        if (number_of_items > 0)
        {
            rMean *= (1.0 / static_cast<double>(number_of_items));
        }
    }

    template <typename TDataType>
    void static CalculateVariance(TDataType& rMean,
                                  TDataType& rVariance,
                                  const ModelPart& rModelPart,
                                  const Variable<TDataType>& rVariable)
    {
        CalculateMean<TDataType>(rMean, rModelPart, rVariable);

        const TContainerType& r_container =
            MethodsUtilities::GetDataContainer<TContainerType>(rModelPart);

        const TDataType& r_initial_value =
            TDataRetrievalFunctor<TContainerItemType>()(*r_container.begin(), rVariable);

        rVariance = rVariable.Zero();
        MethodsUtilities::DataTypeSizeInitializer(rVariance, r_initial_value);

#pragma omp parallel
        {
            TDataType variance = rVariable.Zero();
            MethodsUtilities::DataTypeSizeInitializer(variance, r_initial_value);

#pragma omp for
            for (int i = 0; i < static_cast<int>(r_container.size()); ++i)
            {
                const TContainerItemType& r_item = *(r_container.begin() + i);
                const TDataType& current_value =
                    TDataRetrievalFunctor<TContainerItemType>()(r_item, rVariable);
                MethodsUtilities::DataTypeSizeChecker(current_value, variance);

                variance += MethodsUtilities::RaiseToPower(current_value, 2);
            }
#pragma omp critical
            {
                rVariance += variance;
            }
        }
        rVariance = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(rVariance);
        const unsigned int number_of_items =
            rModelPart.GetCommunicator().GetDataCommunicator().SumAll(
                r_container.size());

        if (number_of_items > 0)
        {
            rVariance *= (1.0 / static_cast<double>(number_of_items));
            rVariance -= MethodsUtilities::RaiseToPower(rMean, 2);
        }
    }

    template <typename TDataType>
    void static GetMax(double& rMax,
                       std::size_t& rId,
                       const std::string& rNormType,
                       const ModelPart& rModelPart,
                       const Variable<TDataType>& rVariable)
    {
        KRATOS_TRY

        const TContainerType& r_container =
            MethodsUtilities::GetDataContainer<TContainerType>(rModelPart);

        rMax = std::numeric_limits<double>::lowest();
        const int multiplication_table = GetMultiplicationTable(rVariable, rNormType);

#pragma omp parallel
        {
            double current_max = std::numeric_limits<double>::lowest();
            std::size_t current_id = 0;
#pragma omp for
            for (int i = 0; i < static_cast<int>(r_container.size()); ++i)
            {
                const TContainerItemType& r_item = *(r_container.begin() + i);
                const TDataType& current_value =
                    TDataRetrievalFunctor<TContainerItemType>()(r_item, rVariable);
                const double value_norm = GetNorm(current_value, multiplication_table);
                if (value_norm > current_max)
                {
                    current_max = value_norm;
                    current_id = r_item.Id();
                }
            }
#pragma omp critical
            {
                if (current_max > rMax)
                {
                    rMax = current_max;
                    rId = current_id;
                }
            }
        }

        rMax = rModelPart.GetCommunicator().GetDataCommunicator().MaxAll(rMax);

        KRATOS_CATCH("");
    }

    template <typename TDataType>
    void static GetMin(double& rMin,
                       std::size_t& rId,
                       const std::string& rNormType,
                       const ModelPart& rModelPart,
                       const Variable<TDataType>& rVariable)
    {
        KRATOS_TRY

        const TContainerType& r_container =
            MethodsUtilities::GetDataContainer<TContainerType>(rModelPart);

        rMin = std::numeric_limits<double>::max();
        const int multiplication_table = GetMultiplicationTable(rVariable, rNormType);

#pragma omp parallel
        {
            double current_min = std::numeric_limits<double>::max();
            std::size_t current_id = 0;
#pragma omp for
            for (int i = 0; i < static_cast<int>(r_container.size()); ++i)
            {
                const TContainerItemType& r_item = *(r_container.begin() + i);
                const TDataType& current_value =
                    TDataRetrievalFunctor<TContainerItemType>()(r_item, rVariable);
                const double value_norm = GetNorm(current_value, multiplication_table);
                if (value_norm < current_min)
                {
                    current_min = value_norm;
                    current_id = r_item.Id();
                }
            }
#pragma omp critical
            {
                if (current_min < rMin)
                {
                    rMin = current_min;
                    rId = current_id;
                }
            }
        }

        rMin = rModelPart.GetCommunicator().GetDataCommunicator().MinAll(rMin);

        KRATOS_CATCH("");
    }

private:
    const int static GetMultiplicationTable(const Variable<double>&, const std::string&)
    {
        return 1;
    }

    const int static GetMultiplicationTable(const Variable<array_1d<double, 3>>& rVariable,
                                            const std::string& rNormType)
    {
        KRATOS_TRY

        int index = -1;
        if (rNormType == "magnitude")
        {
            index = 3;
        }
        else if (rNormType == "component_x")
        {
            index = 0;
        }
        else if (rNormType == "component_y")
        {
            index = 1;
        }
        else if (rNormType == "component_z")
        {
            index = 2;
        }
        else
        {
            KRATOS_ERROR << "Unknown norm type for 3d variable "
                         << rVariable.Name() << ". [ NormType = " << rNormType << " ]\n"
                         << "   Allowed norm types are:\n"
                         << "        magnitude\n"
                         << "        component_x\n"
                         << "        component_y\n"
                         << "        component_z\n";
        }

        return index;

        KRATOS_CATCH("");
    }

    double static GetNorm(const double Value, const int)
    {
        return Value;
    }

    double static GetNorm(const array_1d<double, 3>& rValue, const int Index)
    {
        if (Index < 3)
        {
            return rValue[Index];
        }
        else
        {
            double value_norm = 0.0;
            for (int i = 0; i < 3; ++i)
            {
                value_norm += std::pow(rValue[i], 2);
            }
            return std::sqrt(value_norm);
        }
    }
};
} // namespace SpatialMethods

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_SPATIAL_VARIANCE_METHOD_H_INCLUDED defined
