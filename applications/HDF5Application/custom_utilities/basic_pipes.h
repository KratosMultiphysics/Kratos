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

#pragma once

// Internal includes
#include "containers/variable.h"
#include "includes/process_info.h"
#include "pipe.h"

// Core includes
#include "containers/model.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "utilities/interval_utility.h"

// STL includes
#include <optional>
#include <functional> // std::reference_wrapper
#include <cmath> // std::fmod


namespace Kratos::Pipes {


/// @brief Get a @ref ModelPart from a @ref Model by name.
/// @note Constructible from @ref Parameters with a "model_part_name" string entry.
class ModelPartFromModel : public Traits<const Model&, const ModelPart&>
{
public:
    ModelPartFromModel() = default;

    ModelPartFromModel(const std::string& rModelPartName);

    ModelPartFromModel(std::string&& rModelPartName) noexcept;

    ModelPartFromModel(const Parameters& rParameters);

    ModelPartFromModel(ModelPartFromModel&& rOther) noexcept = default;

    ModelPartFromModel(const ModelPartFromModel& rOther) = default;

    const ModelPart& operator()(const Model& rModel) const
    {return rModel.GetModelPart(mModelPartName);}

private:
    std::string mModelPartName;
}; // class ModelPartFromModel


/// @brief Get the @ref ProcessInfo of a @ref ModelPart.
/// @note No-op constructible from @ref Parameters.
struct ProcessInfoFromModelPart : public Traits<const ModelPart&, const ProcessInfo&>
{
    ProcessInfoFromModelPart() noexcept = default;

    ProcessInfoFromModelPart(const Parameters& rParameters) noexcept {}

    ProcessInfoFromModelPart(ProcessInfoFromModelPart&& rOther) noexcept = default;

    ProcessInfoFromModelPart(const ProcessInfoFromModelPart& rOther) = default;

    const ProcessInfo& operator()(const ModelPart& rModelPart) const
    {return rModelPart.GetProcessInfo();}
}; // struct ProcessInfoFromModelPart


/// @brief Get a variable from @ref ProcessInfo.
/// @note Constructible from @ref Parameters with a "process_info_variable" entry.
/// @warning This type is invalid if default-constructed.
template <class TVariable>
class VariableFromProcessInfo : public Traits<const ProcessInfo&, typename TVariable::Type>
{
public:
    VariableFromProcessInfo() = default;

    VariableFromProcessInfo(const TVariable& rVariable)
        : mVariable({rVariable})
    {}

    VariableFromProcessInfo(const Parameters& rParameters);

    VariableFromProcessInfo(VariableFromProcessInfo&& rOther) noexcept = default;

    VariableFromProcessInfo(const VariableFromProcessInfo& rOther) = default;

    typename TVariable::Type operator()(const ProcessInfo& rProcessInfo) const
    {
        KRATOS_ERROR_IF_NOT(bool(mVariable)) << "uninitialized variable in VariableFromProcessInfo";
        return rProcessInfo[mVariable.value().get()];
    }

private:
    std::optional<std::reference_wrapper<const TVariable>> mVariable;
}; // class VariableFromProcessInfo


/// @brief Get @ref TIME from a @ref ProcessInfo.
/// @note Constructible from @ref Parameters without any requirements.
struct TimeFromProcessInfo : public VariableFromProcessInfo<decltype(TIME)>
{
    TimeFromProcessInfo();

    TimeFromProcessInfo(const Parameters& rParameters);

    TimeFromProcessInfo(TimeFromProcessInfo&& rOther) = default;

    TimeFromProcessInfo(const TimeFromProcessInfo& rOther) = default;
}; // struct TimeFromProcessInfo


/// @brief Get @ref STEP from a @ref ProcessInfo.
/// @note Constructible from @ref Parameters without any requirements.
struct StepFromProcessInfo : public VariableFromProcessInfo<decltype(STEP)>
{
    StepFromProcessInfo();

    StepFromProcessInfo(const Parameters& rParameters);

    StepFromProcessInfo(StepFromProcessInfo&& rOther) = default;

    StepFromProcessInfo(const StepFromProcessInfo& rOther) = default;
}; // struct TimeFromProcessInfo


/// @brief Pipe wrapper for @ref Detail::IntervalUtility.
/// @note Constructible from @ref Parameters (passed on to @ref Detail::IntervalUtility).
template <class TValue>
class IntervalPredicate : public Traits<TValue,bool>
{
public:
    IntervalPredicate() = default;

    IntervalPredicate(TValue begin, TValue end);

    IntervalPredicate(const Parameters& rParameters);

    IntervalPredicate(IntervalPredicate&& rOther) noexcept = default;

    IntervalPredicate(const IntervalPredicate& rOther) = default;

    bool operator()(TValue Value) const
    {return mInterval.IsInInterval(Value);}

private:
    ::Kratos::Detail::IntervalUtility<TValue> mInterval;
}; // class IntervalPredicate


/// @brief Compue the mod of the input.
/// @note Constructible from @ref Parameters with a "mod" entry (@a int or @a double).
template <class TValue>
class Modulo : public Traits<TValue,TValue>
{
public:
    Modulo() noexcept;

    Modulo(TValue modulo) noexcept;

    Modulo(const Parameters& rParameters);

    Modulo(Modulo&& rOther) noexcept = default;

    Modulo(const Modulo& rOther) = default;

    TValue operator()(TValue Value) const
    {
        if constexpr (std::is_integral_v<TValue>) {
            return Value % mModulo;
        } else {
            return std::fmod(Value, mModulo);
        }
    }

private:
    TValue mModulo;
}; // class Modulo


} // namespace Kratos::Pipes


// Template definitions
#include <custom_utilities/basic_pipes_impl.h>
