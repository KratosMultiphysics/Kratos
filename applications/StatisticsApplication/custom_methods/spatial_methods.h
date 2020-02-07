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
#include <functional>
#include <tuple>

// External includes

// Project includes
#include "includes/communicator.h"
#include "includes/define.h"
#include "includes/model_part.h"

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
    TDataType static CalculateSum(const ModelPart& rModelPart,
                                  const Variable<TDataType>& rVariable)
    {
        KRATOS_TRY

        const TContainerType& r_container =
            MethodsUtilities::GetDataContainer<TContainerType>(rModelPart);

        const TDataType& r_initial_value =
            TDataRetrievalFunctor<TContainerItemType>()(*r_container.begin(), rVariable);

        TDataType global_sum = rVariable.Zero();
        MethodsUtilities::DataTypeSizeInitializer(global_sum, r_initial_value);

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
                global_sum += sum;
            }
        }

        global_sum = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(global_sum);
        return global_sum;

        KRATOS_CATCH("");
    }

    // special overloaded method for flags
    int static CalculateSum(const ModelPart& rModelPart, const Flags& rVariable)
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

        sum = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(sum);
        return sum;
    }

    template <typename TDataType>
    TDataType static CalculateMean(const ModelPart& rModelPart,
                                   const Variable<TDataType>& rVariable)
    {
        const TDataType& sum = CalculateSum<TDataType>(rModelPart, rVariable);
        const TContainerType& r_container =
            MethodsUtilities::GetDataContainer<TContainerType>(rModelPart);

        const unsigned int number_of_items =
            rModelPart.GetCommunicator().GetDataCommunicator().SumAll(
                r_container.size());

        if (number_of_items > 0)
        {
            return sum * (1.0 / static_cast<double>(number_of_items));
        }

        return rVariable.Zero();
    }

    template <typename TDataType>
    std::tuple<TDataType, TDataType> static CalculateVariance(const ModelPart& rModelPart,
                                                              const Variable<TDataType>& rVariable)
    {
        TDataType mean = CalculateMean<TDataType>(rModelPart, rVariable);

        const TContainerType& r_container =
            MethodsUtilities::GetDataContainer<TContainerType>(rModelPart);

        const TDataType& r_initial_value =
            TDataRetrievalFunctor<TContainerItemType>()(*r_container.begin(), rVariable);

        TDataType global_variance = rVariable.Zero();
        MethodsUtilities::DataTypeSizeInitializer(global_variance, r_initial_value);

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
                global_variance += variance;
            }
        }
        global_variance =
            rModelPart.GetCommunicator().GetDataCommunicator().SumAll(global_variance);
        const unsigned int number_of_items =
            rModelPart.GetCommunicator().GetDataCommunicator().SumAll(
                r_container.size());

        if (number_of_items > 0)
        {
            global_variance *= (1.0 / static_cast<double>(number_of_items));
            global_variance -= MethodsUtilities::RaiseToPower(mean, 2);
        }

        return std::make_tuple<TDataType, TDataType>(
            std::forward<TDataType>(mean), std::forward<TDataType>(global_variance));
    }

    template <typename TDataType>
    std::tuple<double, std::size_t> static GetMax(const std::string& rNormType,
                                                  const ModelPart& rModelPart,
                                                  const Variable<TDataType>& rVariable)
    {
        KRATOS_TRY

        const TContainerType& r_container =
            MethodsUtilities::GetDataContainer<TContainerType>(rModelPart);

        double global_max = std::numeric_limits<double>::lowest();
        std::size_t global_id = 0;
        const auto& norm_method =
            MethodsUtilities::GetNormMethod<TDataType>(rVariable, rNormType);

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
                const double value_norm = norm_method(current_value);
                if (value_norm > current_max)
                {
                    current_max = value_norm;
                    current_id = r_item.Id();
                }
            }
#pragma omp critical
            {
                if (current_max > global_max)
                {
                    global_max = current_max;
                    global_id = current_id;
                }
            }
        }

        const DataCommunicator& r_communicator =
            rModelPart.GetCommunicator().GetDataCommunicator();
        const std::vector<double> global_max_value_array =
            r_communicator.AllGather(std::vector<double>{{global_max}});
        const std::vector<std::size_t> global_max_id_array =
            r_communicator.AllGather(std::vector<std::size_t>{{global_id}});

        for (std::size_t i = 0; i < global_max_value_array.size(); ++i)
        {
            if (global_max_value_array[i] > global_max)
            {
                global_max = global_max_value_array[i];
                global_id = global_max_id_array[i];
            }
        }

        return std::make_tuple<double, std::size_t>(
            std::forward<double>(global_max), std::forward<std::size_t>(global_id));

        KRATOS_CATCH("");
    }

    template <typename TDataType>
    std::tuple<double, std::size_t> static GetMin(const std::string& rNormType,
                                                  const ModelPart& rModelPart,
                                                  const Variable<TDataType>& rVariable)
    {
        KRATOS_TRY

        const TContainerType& r_container =
            MethodsUtilities::GetDataContainer<TContainerType>(rModelPart);

        double global_min = std::numeric_limits<double>::max();
        std::size_t global_id = 0;
        const auto& norm_method =
            MethodsUtilities::GetNormMethod<TDataType>(rVariable, rNormType);

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
                const double value_norm = norm_method(current_value);
                if (value_norm < current_min)
                {
                    current_min = value_norm;
                    current_id = r_item.Id();
                }
            }
#pragma omp critical
            {
                if (current_min < global_min)
                {
                    global_min = current_min;
                    global_id = current_id;
                }
            }
        }

        const DataCommunicator& r_communicator =
            rModelPart.GetCommunicator().GetDataCommunicator();
        const std::vector<double> global_min_value_array =
            r_communicator.AllGather(std::vector<double>{{global_min}});
        const std::vector<std::size_t> global_min_id_array =
            r_communicator.AllGather(std::vector<std::size_t>{{global_id}});

        for (std::size_t i = 0; i < global_min_value_array.size(); ++i)
        {
            if (global_min_value_array[i] < global_min)
            {
                global_min = global_min_value_array[i];
                global_id = global_min_id_array[i];
            }
        }

        return std::make_tuple<double, std::size_t>(
            std::forward<double>(global_min), std::forward<std::size_t>(global_id));

        KRATOS_CATCH("");
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
