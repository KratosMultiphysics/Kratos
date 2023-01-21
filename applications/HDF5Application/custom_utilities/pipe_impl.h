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

// --- HDF5 Includes ---
// Included from "custom_utilities/pipe.h"
#include "custom_utilities/pipe.h" // unnecessary include to get language servers working


namespace Kratos::Pipes
{


/// @brief Operator for calling operator() of the pipe.
template <class TInput, class TPipe, std::enable_if_t<IsPipe<TPipe>::value && std::is_convertible_v<TInput,typename TPipe::InputType>, bool>>
typename TPipe::OutputType operator>>(TInput&& rInput, const TPipe& rPipe)
{
    return rPipe(rInput);
}


template <class TInputPipe, class TOutputPipe>
CompoundPipe<TInputPipe,TOutputPipe>::CompoundPipe()
    : mInputPipe(),
      mOutputPipe()
{
}


namespace Detail {


/**
 *  @brief Class for constructing compound and simple pipes.
 *  @details A lot of heavy machinery is required to make the compile-time choice
 *           between constructing a compound- or simple pipe, but at the end of the
 *           day, it all happens at compile-time and does not affect runtime performance.
 */
template <class TPipe>
class Factory
{
private:
    /// @details Factory::Make needs to be recursively called, so all pipe factories must be friends.
    template <class>
    friend class Factory;

    /// @brief Dummy struct for overload resolution.
    struct Simple {};

    /// @brief Derived dummy struct for overload resolution.
    struct Compound : Simple {};

    /// @brief Dummy struct for SFINAE.
    template <class>
    struct Dummy {using Type = int;};

    /**
     *  @brief Pipe factory for @ref CompoundPipe s.
     *  @details Ignoring all the machinery to make SFINAE possible, this function
     *           just constructs an input- and output pipe, then assembles them into
     *           a @ref CompoundPipe.
     *  @param[out] pPipe: pointer to a preallocated compound pipe the output will be assigned to.
     *  @param[in] rParameters: input parameters to construct the pipeline with. Must consist of
     *                          an array of subparameters to construct individual segments with.
     *  @param[inout] rSubParamOffset: index of the subparameters within @a rParameters to construct
     *                                 the input pipe with. This index is updated and will point to
     *                                 one past the last subparameter used for constructing the
     *                                 output pipe.
     */
    template <class TCompoundPipe, typename Dummy<decltype(TCompoundPipe::mInputPipe)>::Type = 0>
    void Make(TCompoundPipe* pPipe, const Parameters& rParameters, std::size_t& rSubParamOffset, Compound)
    {
        static_assert(std::is_same_v<TPipe, TCompoundPipe>);
        // Construct with placement new at the specified location.
        // This avoids move/copy assigning the pipe.
        new(pPipe) TCompoundPipe(
            Factory<typename TPipe::InputPipe>::Make(rParameters, rSubParamOffset),
            Factory<typename TPipe::OutputPipe>::Make(rParameters, rSubParamOffset)
        );
    }

    /**
     *  @brief Pipe factory for simple pipes.
     *  @details Ignoring all the machinery to make SFINAE possible, this function
     *           just constructs a simple (non-compound) pipe.
     *  @param[out] pPipe: pointer to a preallocated compound pipe the output will be assigned to.
     *  @param[in] rParameters: input parameters to construct the pipeline with. Must consist of
     *                          an array of subparameters to construct individual segments with.
     *  @param[inout] rSubParamOffset: index of the subparameters within @a rParameters to construct
     *                                 the segment with. This index is incremented once in this function.
     */
    template <class TSimplePipe>
    void Make(TSimplePipe* pPipe, const Parameters& rParameters, std::size_t& rSubParamOffset, Simple)
    {
        static_assert(std::is_same_v<TPipe, TSimplePipe>);
        KRATOS_ERROR_IF_NOT(rSubParamOffset < rParameters.size())
            << "No subparameter found for constructing "
            << typeid(TSimplePipe).name()
            << " pipe at index "
            << rSubParamOffset
            << " from parameters: "
            << rParameters;

        // Construct with placement new at the specified location.
        // This avoids move/copy assigning the pipe.
        KRATOS_TRY
        new(pPipe) TSimplePipe(rParameters[rSubParamOffset]);
        KRATOS_CATCH(
            "Failed to construct a "
            << typeid(TSimplePipe).name()
            << " pipe segment at index "
            << rSubParamOffset << " with parameters: "
            << rParameters);
        ++rSubParamOffset;
    }

    /**
     *  @brief Pipe factory for simple or compound pipes.
     *  @param[in] rParameters: @ref Parameters consisting of an array of subparameters to
     *                          construct the pipe or pipeline with.
     *  @param[inout] rSubParamOffset: index of the subparameters within @a rParameters to begin
     *                                 constructing the pipe with. This index is updated and will
     *                                 point to one past the last subparameter used for constructing
     *                                 the last pipe in the pipeline.
     */
    static TPipe Make(const Parameters& rParameters, std::size_t& rSubParamOffset)
    {
        // Allocate memory on the stack to later construct the pipe into.
        // This avoids default constructing, then later copy/move assigning
        // it from a different place.
        char buffer[sizeof(TPipe)];
        TPipe* p_pipe = reinterpret_cast<TPipe*>(buffer);
        Factory().Make(
            p_pipe,
            rParameters,
            rSubParamOffset,
            Compound() // assume a compound pipe, the compiler will cast it up to Simple if that's not the case
        );
        return *reinterpret_cast<TPipe*>(p_pipe);
    }

public:
    /// @brief Construct a simple or compound pipe from a @ref Parameters consisting of an array of subparameters.
    /// @param[in] rParameters: @ref Parameters consisting of an array of subparameters to construct the pipe or
    ///                         pipeline with.
    static TPipe Make(const Parameters& rParameters)
    {
        std::size_t sub_param_offset = 0;
        return Factory::Make(rParameters, sub_param_offset);
    }
}; // struct Factory


} // namespace Detail


template <class TInputPipe, class TOutputPipe>
template <class TInPipe, class TOutPipe>
CompoundPipe<TInputPipe,TOutputPipe>::CompoundPipe(const Parameters& rParameters,
                                                   std::enable_if_t<std::is_constructible_v<TInPipe,Parameters>
                                                       && std::is_constructible_v<TOutPipe,Parameters>>*)
    : CompoundPipe(Detail::Factory<CompoundPipe>::Make(rParameters))
{
}


template <class TInputPipe, class TOutputPipe>
CompoundPipe<TInputPipe,TOutputPipe>::CompoundPipe(TInputPipe&& rInputPipe, TOutputPipe&& rOutputPipe) noexcept
    : mInputPipe(std::move(rInputPipe)),
      mOutputPipe(std::move(rOutputPipe))
{
}


template <class TInputPipe, class TOutputPipe>
CompoundPipe<TInputPipe,TOutputPipe>::CompoundPipe(const TInputPipe& rInputPipe, const TOutputPipe& rOutputPipe)
    : mInputPipe(rInputPipe),
      mOutputPipe(rOutputPipe)
{
}


template <class TInputPipe, class TOutputPipe>
inline typename CompoundPipe<TInputPipe,TOutputPipe>::OutputType
CompoundPipe<TInputPipe,TOutputPipe>::operator()(typename CompoundPipe::InputType Input) const
{
    return mOutputPipe(mInputPipe(Input));
}


template <class TInputPipe,
          class TOutputPipe,
          std::enable_if_t<
            IsPipe<TInputPipe>::value && IsPipe<TOutputPipe>::value,
            bool>
          >
CompoundPipe<TInputPipe,TOutputPipe> operator|(TInputPipe&& rInputPipe, TOutputPipe&& rOutputPipe)
{
    return CompoundPipe<TInputPipe,TOutputPipe>(std::move(rInputPipe), std::move(rOutputPipe));
}


template <class TInputPipe,
          class TOutputPipe,
          std::enable_if_t<
            IsPipe<TInputPipe>::value && IsPipe<TOutputPipe>::value,
            bool>
          >
CompoundPipe<TInputPipe,TOutputPipe> operator|(const TInputPipe& rInputPipe, const TOutputPipe& rOutputPipe)
{
    return CompoundPipe<TInputPipe,TOutputPipe>(rInputPipe, rOutputPipe);
}


} // namespace Kratos::Pipes
