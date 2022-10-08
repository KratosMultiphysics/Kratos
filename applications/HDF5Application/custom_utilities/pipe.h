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

// Core includes
#include "includes/kratos_parameters.h"

// STL includes
#include <type_traits> // std::true_type, std::false_type
#include <utility> // std::move


namespace Kratos
{

/**
 *  @brief Unix-inspired pipes in C++.
 *
 *  @details A pipe is a modular sequence of operations that takes an immutable input and produces an output.
 *           Modularity is provided via static polymorphism (templates) to avoid the overhead incurred by virtual
 *           functions. Therefore, there is no base class for a @a pipe; rather, a class must satisfy a set of
 *           requirements to be considered a @a Pipe. These requirements are:
 *              - has member type @a InputType
 *              - has member type @a OutputType
 *              - has member function operator() with the following signature:
 *                @code OutputType operator()(@a InputType) const @endcode
 *              - must not have side-effects (not checked)
 *              - must not mutate its input (not checked)
 *              - must be default, move, and copy constructible (not checked)
 *              - [optional]: should be constructible from @ref Parameters (required for @ref Factory).
 *
 *  @details Pipes are meant to be stitched together to form a pipeline capable of performing more complex
 *           tasks. Two operators can be used for this purpose:
 *              - @ref Kratos::Pipes::operator| takes two (possibly different type of) pipe instances and produces a new pipe
 *                that combines the two. For example:
 *                @code auto compound_pipe = pipe_1 | pipe_2; @endcode
 *                The resulting @a compound_pipe takes an input, feeds it to @a pipe_1, then feed the result
 *                into @a pipe_2, then returns its result. The benefit of this approach over @a operator>>
 *                is that the resulting @a compound_pipe is a pipe as well, and can be passed around to other
 *                functions / classes, or combined with other pipes.
 *              - @ref Kratos::Pipes::operator>> feeds an input into a pipe. It is equivalent to calling @a operator(), but is
 *                more readable when pushing through multiple pipes. Example:
 *                @code auto result = input >> pipe_1 >> pipe_2; @endcode
 *                The benefit of this approach over @a operator| is that no compound pipe instance is constructed.
 *
 *  @details When defining a new pipe class, make sure to derive from @ref Traits, which defines the necessary type aliases
 *           for your new class to be considered a pipe. The integral constant @ref IsPipe checks whether a type is a valid pipe.
 *           Usage: @code IsPipe<MyPipe>::value @endcode. A type must satisfy the pipe requirements in order to be usable in
 *           @ref operator| and @ref operator>>.
 *
 *  @details Pipes that are constructible from @ref Parameters can be used with @ref Factory. @ref Factory recursively
 *           constructs each pipe segment by passing the same @ref Parameters instance. This can be very handy for constructing
 *           longer pipes or using pipes in generic templates.
 */
namespace Pipes {


/**
 *  @brief Metaclass containing type information every pipe must inherit from.
 *  @tparam TInput: cv-qualified argument type of the pipe's operator().
 *  @tparam TOutput: qualified output type of the pipe's operator().
 *  @note Every pipe must be default, move, and copy constructible + have @a OutputType operator()(@a InputType) const.
 */
template <class TInput, class TOutput>
struct Traits
{
    /// @brief cv-qualified argument type of operator()
    using InputType = TInput;

    /// @brief qualified return type of operator()
    using OutputType = TOutput;
}; // struct Traits


/**
 *  @brief Bool constant checking whether @a TPipe satisfies the requirements of a pipe.
 *  @details Pipe requirements:
 *              - has non-void member aliases @a InputType and @a OutputType
 *              - has member function with the following signature
 *                @a OutputType operator()(InputType) const
 *              - is default, move, and copy constructible @todo @matekelemen
 */
template <class TPipe>
using IsPipe = std::integral_constant<bool,
    // Has operator() with the following signature: OutputType operator()(InputType) const
    std::is_same_v<decltype(std::declval<const TPipe>().operator()(std::declval<typename TPipe::InputType>())), typename TPipe::OutputType>

    // Has non-void InputType
    && !std::is_same_v<typename TPipe::InputType,void>

    // Has non-void OutputType
    && !std::is_same_v<typename TPipe::OutputType,void>
>;


/// @brief Operator for calling operator() of the pipe.
template <class TInput, class TPipe, std::enable_if_t<IsPipe<TPipe>::value && std::is_convertible_v<TInput,typename TPipe::InputType>, bool> = true>
typename TPipe::OutputType operator>>(TInput&& rInput, const TPipe& rPipe)
{
    return rPipe(rInput);
}


/**
 *  @brief A pipe that takes the output of one pipe and feeds its result into another.
 *  @tparam TInputPipe: input pipe type (left hand side).
 *  @tparam TOutputPipe: output pipe type (right hand side).
 *  @details @ref CompoundPipe can be constructed by copying or moving its input- and output
 *           pipes, which are stored by @a value in the resulting @ref CompoundPipe.
 */
template <class TInputPipe, class TOutputPipe>
class CompoundPipe : public Traits<typename TInputPipe::InputType, typename TOutputPipe::OutputType>
{
public:
    using InputPipe = TInputPipe;

    using OutputPipe = TOutputPipe;

    /// @brief Default construct the input- and output pipes.
    CompoundPipe()
        : mInputPipe(),
          mOutputPipe()
    {}

    /// @brief Construct from a @ref Parameters object.
    /// @details Enabled only if both the input- and the output pipes
    ///          are constructible from a @ref Parameters instance.
    template <class TInPipe = TInputPipe, class TOutPipe = TOutputPipe>
    CompoundPipe(const Parameters& rParameters,
                 std::enable_if_t<std::is_constructible_v<TInPipe,Parameters>
                 && std::is_constructible_v<TOutPipe,Parameters>>* = 0)
        : mInputPipe(rParameters),
          mOutputPipe(rParameters)
    {}

    /// @brief Move construct the input- and output pipes.
    CompoundPipe(TInputPipe&& rInputPipe, TOutputPipe&& rOutputPipe) noexcept
        : mInputPipe(std::move(rInputPipe)),
          mOutputPipe(std::move(rOutputPipe))
    {}

    /// @brief Move construct the input pipe, copy construct the output pipe.
    CompoundPipe(TInputPipe&& rInputPipe, const TOutputPipe& rOutputPipe) noexcept
        : mInputPipe(std::move(rInputPipe)),
          mOutputPipe(rOutputPipe)
    {}

    /// @brief Copy construct the input pipe, move construct the output pipe.
    CompoundPipe(const TInputPipe& rInputPipe, TOutputPipe&& rOutputPipe) noexcept
        : mInputPipe(rInputPipe),
          mOutputPipe(std::move(rOutputPipe))
    {}

    /// @brief Copy construct the input- and output pipes.
    CompoundPipe(const TInputPipe& rInputPipe, const TOutputPipe& rOutputPipe)
        : mInputPipe(rInputPipe),
          mOutputPipe(rOutputPipe)
    {}

    CompoundPipe(CompoundPipe&& rOther) noexcept = default;

    CompoundPipe(const CompoundPipe& rOther) = default;

    /// @brief Feed the result of the input pipe into the output pipe.
    typename CompoundPipe::OutputType operator()(typename CompoundPipe::InputType Input) const
    {return mOutputPipe(mInputPipe(Input));}

private:
    TInputPipe mInputPipe;

    TOutputPipe mOutputPipe;
}; // class CompoundPipe


/// @brief Construct a pipe that takes the output of an input pipe and feeds it into an output pipe.
template <class TInputPipe,
          class TOutputPipe,
          std::enable_if_t<
            IsPipe<TInputPipe>::value && IsPipe<TOutputPipe>::value,
            bool> = true
          >
CompoundPipe<TInputPipe,TOutputPipe> operator|(TInputPipe&& rInputPipe, TOutputPipe&& rOutputPipe)
{
    return CompoundPipe<TInputPipe,TOutputPipe>(std::move(rInputPipe), std::move(rOutputPipe));
}


/// @brief Construct a pipe that takes the output of an input pipe and feeds it into an output pipe.
template <class TInputPipe,
          class TOutputPipe,
          std::enable_if_t<
            IsPipe<TInputPipe>::value && IsPipe<TOutputPipe>::value,
            bool> = true
          >
CompoundPipe<TInputPipe,TOutputPipe> operator|(const TInputPipe& rInputPipe, TOutputPipe&& rOutputPipe)
{
    return CompoundPipe<TInputPipe,TOutputPipe>(rInputPipe, std::move(rOutputPipe));
}


/// @brief Construct a pipe that takes the output of an input pipe and feeds it into an output pipe.
template <class TInputPipe,
          class TOutputPipe,
          std::enable_if_t<
            IsPipe<TInputPipe>::value && IsPipe<TOutputPipe>::value,
            bool> = true
          >
CompoundPipe<TInputPipe,TOutputPipe> operator|(TInputPipe&& rInputPipe, const TOutputPipe& rOutputPipe)
{
    return CompoundPipe<TInputPipe,TOutputPipe>(std::move(rInputPipe), rOutputPipe);
}


/// @brief Construct a pipe that takes the output of an input pipe and feeds it into an output pipe.
template <class TInputPipe,
          class TOutputPipe,
          std::enable_if_t<
            IsPipe<TInputPipe>::value && IsPipe<TOutputPipe>::value,
            bool> = true
          >
CompoundPipe<TInputPipe,TOutputPipe> operator|(const TInputPipe& rInputPipe, const TOutputPipe& rOutputPipe)
{
    return CompoundPipe<TInputPipe,TOutputPipe>(rInputPipe, rOutputPipe);
}


/// @brief Convenience type alias for complex pipes.
template <class ...TPipes>
using Pipeline = decltype((... | std::declval<TPipes>()));


/// @brief Bool constant checking whether a type is a @ref CompoundPipe.
template <class TPipe>
using IsCompoundPipe = std::integral_constant<bool,
    IsPipe<TPipe>::value
    && IsPipe<typename TPipe::InputPipe>::value
    && IsPipe<typename TPipe::OutputPipe>::value
>;


/**
 *  @brief A factory for creating pipes from a single @ref Parameters object.
 *  @note Every pipe segment within @a TPipe must be constructible from @ref Parameters too.
 */
template <class TPipe, std::enable_if_t<std::is_constructible_v<TPipe,Parameters>, bool> = true>
class Factory
{
private:
    // Leaf pipe factory.
    template <class TSimplePipe, class Enable = std::enable_if_t<IsCompoundPipe<TSimplePipe>::value,int>>
    struct Impl
    {
        static TSimplePipe Make(const Parameters& rParameters)
        {return TSimplePipe(rParameters);}
    }; // struct Impl

    // Compound pipe factory
    template <class TCompoundPipe>
    struct Impl<TCompoundPipe,void>
    {
        static TCompoundPipe Make(const Parameters& rParameters)
        {
            auto input_pipe = Factory<typename TCompoundPipe::InputPipe>::Make(rParameters);
            auto output_pipe = Factory<typename TCompoundPipe::OutputPipe>::Make(rParameters);
            return TCompoundPipe(std::move(input_pipe), std::move(output_pipe));
        }
    }; // struct Impl

public:
    static TPipe Make(const Parameters& rParameters)
    {return Impl<TPipe>::Make(rParameters);}
}; // class Factory


} // namespace Pipes
} // namespace Kratos
