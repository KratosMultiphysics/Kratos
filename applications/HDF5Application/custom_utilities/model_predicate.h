//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: HDF5Application/license.txt
//
//  Main author:     Máté Kelemen
//

#ifndef KRATOS_HDF5_MODEL_PREDICATE_H
#define KRATOS_HDF5_MODEL_PREDICATE_H

// Internal includes
#include "basic_pipes.h"

// Core includes
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "utilities/interval_utility.h"

// STL includes


namespace Kratos::HDF5 {


/// @brief Base class for functors that take a @ref Model and return a bool.
struct ModelPredicate : Pipes::Traits<const Model&, bool>
{
    virtual ~ModelPredicate() {}

    virtual bool operator()(const Model& rModel) const = 0;
}; // struct ModelPredicate


template <class TPipe, class TPipeFactory = Pipes::Factory<TPipe>>
class PipedModelPredicate : public ModelPredicate
{
public:
    PipedModelPredicate()
        : mPipe(TPipeFactory::Make(Parameters()))
    {}

    PipedModelPredicate(const Parameters& rParameters)
        : mPipe(TPipeFactory::Make(rParameters))
    {}

    PipedModelPredicate(PipedModelPredicate&& rOther) noexcept = default;

    PipedModelPredicate(const PipedModelPredicate& rOther) = default;

    bool operator()(const Model& rModel) const override
    {return mPipe(rModel);}

private:
    TPipe mPipe;
}; // class PipedModelPredicate


/**
 *  @brief Check whether @ref TIME in a @ref ModelPart is within an interval.
 *  @details Model => ModelPart => ProcessInfo => TIME => IntervalUtility::IsIninterval.
 *           Required parameters (other settings ignored):
 *           @code
 *           {
 *              "model_part_name" : "",
 *              "interval" : ["Begin", "End"]
 *           }
 *           @endcode
 *  @note See @ref IntervalUtility for details.
 */
using TimeIntervalPredicate = PipedModelPredicate<Pipes::Pipeline<
    Pipes::ModelPartFromModel,
    Pipes::ProcessInfoFromModelPart,
    Pipes::TimeFromProcessInfo,
    Pipes::IntervalPredicate<double>
>>;


/**
 *  @brief Check whether @ref STEP in a @ref ModelPart is within an interval.
 *
 *  @details Model => ModelPart => ProcessInfo => STEP => DiscreteIntervalUtility::IsIninterval.
 *           Required parameters (other settings ignored):
 *           @code
 *           {
 *              "model_part_name" : "",
 *              "interval" : ["Begin", "End"]
 *           }
 *           @endcode
 *
 *  @note See @ref DiscreteIntervalUtility for details.
 */
using StepIntervalPredicate = PipedModelPredicate<Pipes::Pipeline<
    Pipes::ModelPartFromModel,
    Pipes::ProcessInfoFromModelPart,
    Pipes::StepFromProcessInfo,
    Pipes::IntervalPredicate<double>
>>;


/**
 *  @brief Check whether @ref TIME in a @ref ModelPart is within a cyclic interval.
 *
 *  @details Model => ModelPart => ProcessInfo => TIME => Modulo => IntervalUtility::IsInInterval.
 *           Required parameters (other settings ignored):
 *           @code
 *           {
 *              "model_part_name" : "",
 *              "mod" : 0,
 *              "interval" : ["Begin", "End"]
 *           }
 *           @endcode
 *
 *           Example with @code {"mod" : 12.0, "interval" : [3.0, 6.0]} @endcode
 *           @code
 *           TIME:   12          24          36          48          60
 *           ... ----|-----------|-----------|-----------|-----------|---- ...
 *                      ++++        ++++        ++++        ++++        ++
 *           @endcode
 *
 *  @note See @ref IntervalUtility for details.
 */
 using PeriodicTimeIntervalPredicate = PipedModelPredicate<Pipes::Pipeline<
    Pipes::ModelPartFromModel,
    Pipes::ProcessInfoFromModelPart,
    Pipes::TimeFromProcessInfo,
    Pipes::Modulo<double>,
    Pipes::IntervalPredicate<double>
>>;


/**
 *  @brief Check whether @ref STEP in a @ref ModelPart is within a cyclic interval.
 *
 *  @details Model => ModelPart => ProcessInfo => STEP => Modulo => DiscreteIntervalUtility::IsInInterval.
 *           Required parameters (other settings ignored):
 *           @code
 *           {
 *              "model_part_name" : "",
 *              "mod" : 0,
 *              "interval" : ["Begin", "End"]
 *           }
 *           @endcode
 *
 *           Example with @code {"mod" : 12, "interval" : [3, 6]} @endcode
 *           @code
 *           TIME:   12          24          36          48          60
 *           ... ----|-----------|-----------|-----------|-----------|---- ...
 *                      ++++        ++++        ++++        ++++        ++
 *           @endcode
 *
 *  @note See @ref DiscreteIntervalUtility for details.
 */
 using PeriodicStepIntervalPredicate = PipedModelPredicate<Pipes::Pipeline<
    Pipes::ModelPartFromModel,
    Pipes::ProcessInfoFromModelPart,
    Pipes::StepFromProcessInfo,
    Pipes::Modulo<int>,
    Pipes::IntervalPredicate<int>
>>;


} // namespace Kratos::HDF5

#endif
