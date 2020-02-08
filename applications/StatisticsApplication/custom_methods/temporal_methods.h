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

#if !defined(KRATOS_TEMPORAL_METHODS_H_INCLUDED)
#define KRATOS_TEMPORAL_METHODS_H_INCLUDED

// System includes
#include <cmath>
#include <vector>

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

namespace TemporalMethods
{
class TemporalMethod
{
public:
    /// Pointer definition of RansApplyFlagProcess
    KRATOS_CLASS_POINTER_DEFINITION(TemporalMethod);

    TemporalMethod(ModelPart& rModelPart) : mrModelPart(rModelPart)
    {
    }

    virtual void CalculateStatistics(const double DeltaTime)
    {
        KRATOS_TRY

        KRATOS_ERROR << "Calling base class Execute method. Please implement "
                        "it in derrived class.";

        KRATOS_CATCH("");
    }

    virtual void InitializeStatisticsMethod()
    {
        KRATOS_TRY

        KRATOS_ERROR << "Calling base class initialize. Please implement "
                        "it in derrived class.";

        KRATOS_CATCH("");
    }

    ModelPart& GetModelPart() const
    {
        return mrModelPart;
    }

private:
    ModelPart& mrModelPart;
};

template <typename TContainerType, typename TContainerItemType, template <typename T> typename TDataRetrievalFunctor, template <typename T> typename TDataStorageFunctor>
class ContainerTemporalMethods
{
public:
    template <typename TDataType>
    class VarianceMethod : public TemporalMethod
    {
    public:
        /// Pointer definition of RansApplyFlagProcess
        KRATOS_CLASS_POINTER_DEFINITION(VarianceMethod);

        VarianceMethod(ModelPart& rModelPart,
                       const Variable<TDataType>& rInputVariable,
                       const Variable<TDataType>& rOutputMeanVariable,
                       const Variable<TDataType>& rOutputVarianceVariable)
            : TemporalMethod(rModelPart),
              mrInputVariable(rInputVariable),
              mrOutputMeanVariable(rOutputMeanVariable),
              mrOutputVarianceVariable(rOutputVarianceVariable),
              mTotalTime(0.0),
              mInitialized(false)
        {
        }

        void CalculateStatistics(const double DeltaTime) override
        {
            if (!mInitialized)
                InitializeStatisticsMethod();

            TContainerType& r_container =
                MethodsUtilities::GetDataContainer<TContainerType>(this->GetModelPart());

            const int number_of_items = r_container.size();
#pragma omp parallel for
            for (int i = 0; i < number_of_items; ++i)
            {
                TContainerItemType& r_item = *(r_container.begin() + i);
                const TDataType& r_input_value =
                    TDataRetrievalFunctor<TContainerItemType>()(r_item, mrInputVariable);
                TDataType& r_output_mean_value =
                    TDataStorageFunctor<TContainerItemType>()(r_item, mrOutputMeanVariable);
                TDataType& r_output_variance_value =
                    TDataStorageFunctor<TContainerItemType>()(r_item, mrOutputVarianceVariable);

                MethodsUtilities::DataTypeSizeChecker(r_input_value, r_output_mean_value);

                CalculateMeanAndVariance(r_output_mean_value, r_output_variance_value,
                                         r_input_value, DeltaTime);
            }
        }

        void InitializeStatisticsMethod() override
        {
            TContainerType& r_container =
                MethodsUtilities::GetDataContainer<TContainerType>(this->GetModelPart());
            if (r_container.size() > 0)
            {
                const int number_of_items = r_container.size();
#pragma omp parallel for
                for (int i = 0; i < number_of_items; ++i)
                {
                    TContainerItemType& r_item = *(r_container.begin() + i);
                    const TDataType& r_input_value =
                        TDataRetrievalFunctor<TContainerItemType>()(r_item, mrInputVariable);
                    TDataType& r_output_mean_value =
                        TDataStorageFunctor<TContainerItemType>()(r_item, mrOutputMeanVariable);
                    TDataType& r_output_variance_value =
                        TDataStorageFunctor<TContainerItemType>()(
                            r_item, mrOutputVarianceVariable);
                    r_output_mean_value = mrOutputMeanVariable.Zero();
                    r_output_variance_value = mrOutputVarianceVariable.Zero();
                    MethodsUtilities::DataTypeSizeInitializer<TDataType>(
                        r_output_mean_value, r_input_value);
                    MethodsUtilities::DataTypeSizeInitializer<TDataType>(
                        r_output_variance_value, r_input_value);
                }
            }
            mTotalTime = 0.0;
            mInitialized = true;
        }

    private:
        const Variable<TDataType>& mrInputVariable;
        const Variable<TDataType>& mrOutputMeanVariable;
        const Variable<TDataType>& mrOutputVarianceVariable;
        double mTotalTime;
        bool mInitialized;

        void CalculateMeanAndVariance(TDataType& rMean,
                                      TDataType& rVariance,
                                      const TDataType& rNewDataPoint,
                                      const double DeltaTime)
        {
            const double new_total_time = mTotalTime + DeltaTime;
            const TDataType new_mean =
                (rMean * mTotalTime + rNewDataPoint) * (1.0 / new_total_time);
            rVariance =
                ((rVariance + MethodsUtilities::RaiseToPower<TDataType>(new_mean - rMean, 2)) * mTotalTime +
                 MethodsUtilities::RaiseToPower<TDataType>(rNewDataPoint - new_mean, 2)) *
                (1.0 / new_total_time);
            rMean = new_mean;
            mTotalTime = new_total_time;
        }
    };
};

using NodeType = ModelPart::NodeType;
using ElementType = ModelPart::ElementType;
using ConditionType = ModelPart::ConditionType;

using NodesContainerType = ModelPart::NodesContainerType;
using ElementsContainerType = ModelPart::ElementsContainerType;
using ConditionsContainerType = ModelPart::ConditionsContainerType;

template <template <typename T> typename TDataStorageFunctor>
class HistoricalTemporalMethods
    : public ContainerTemporalMethods<NodesContainerType, NodeType, MethodsUtilities::HistoricalDataValueRetrievalFunctor, TDataStorageFunctor>
{
};

class HistoricalInputHistoricalOutputTemporalMethods
    : public HistoricalTemporalMethods<MethodsUtilities::HistoricalDataValueRetrievalFunctor>
{
};

class HistoricalInputNonHistoricalOutputTemporalMethods
    : public HistoricalTemporalMethods<MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>
{
};

class NodalNonHistoricalTemporalMethods
    : public ContainerTemporalMethods<NodesContainerType, NodeType, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>
{
};

class ConditionNonHistoricalTemporalMethods
    : public ContainerTemporalMethods<ConditionsContainerType, ConditionType, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>
{
};

class ElementNonHistoricalTemporalMethods
    : public ContainerTemporalMethods<ElementsContainerType, ElementType, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>
{
};

} // namespace TemporalMethods
} // namespace Kratos
#endif // KRATOS_TEMPORAL_METHODS_H_INCLUDED