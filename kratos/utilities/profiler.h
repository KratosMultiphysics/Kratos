//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Máté Kelemen
//

#pragma once

// Project includes
#include "includes/code_location.h"

// System includes
#include <unordered_map>
#include <filesystem>
#include <chrono>
#include <thread>
#include <optional>
#include <list>
#include <mutex>


namespace Kratos::Internals {


template <class TTimeUnit>
class KRATOS_API(KRATOS_CORE) Profiler
{
private:
    using TimeUnit = TTimeUnit;

    using Duration = TimeUnit;

    using Clock = std::chrono::high_resolution_clock;

    /// @brief Class for identifying a profiled scope and aggregating its stats.
    class Item
    {
    public:
        explicit Item(CodeLocation&& rLocation);

    private:
        Item(std::size_t CallCount,
             Duration CumulativeDuration,
             Duration MinDuration,
             Duration MaxDuration,
             CodeLocation&& rLocation);

        Item& operator+=(const Item& rOther);

    private:
        friend class Profiler;

        unsigned mRecursionLevel;

        std::size_t mCallCount;

        Duration mCumulative;

        Duration mMin;

        Duration mMax;

        CodeLocation mLocation;
    }; // class Item

    struct SourceLocationHash
    {
        std::size_t operator()(const CodeLocation& rArgument) const
        {
            return std::hash<std::string>()(rArgument.GetFileName() + rArgument.GetFunctionName());
        }
    };

    struct SourceLocationEquality
    {
        bool operator()(const CodeLocation& rLhs,
                        const CodeLocation& rRhs) const
        {
            return (rLhs.GetFileName() == rRhs.GetFileName())
                   && (rLhs.GetFunctionName() == rRhs.GetFunctionName());
        }
    };

public:
    /// @brief RAII wrapper for updating an @ref Item.
    class Scope
    {
    public:
        ~Scope();

    private:
        Scope(Item& rItem);

        Scope(Item& rItem, std::chrono::high_resolution_clock::time_point Begin);

        Scope(Scope&&) = delete;

        Scope(const Scope&) = delete;

        Scope& operator=(Scope&&) = delete;

        Scope& operator=(const Scope&) = delete;

    private:
        friend class Profiler;

        /// @brief A raw reference to the associated @ref Item this @ref Scope updates.
        /// @details @see Profiler::mItemContainerMap for guarantees why this reference
        ///          doesn't become dangling.
        Item& mrItem;

        const std::chrono::high_resolution_clock::time_point mBegin;
    }; // class Scope

    using ItemMap = std::unordered_map<
        CodeLocation,
        Item,
        SourceLocationHash,
        SourceLocationEquality
    >;

public:
    Profiler();

    Profiler(Profiler&& rOther) = default;

    Profiler(std::filesystem::path&& r_outputPath);

    ~Profiler();

    Profiler& operator=(Profiler&& rOther) = default;

    [[nodiscard]] Item& Create(CodeLocation&& rItem);

    [[nodiscard]] Scope Profile(Item& rItem);

    /// @brief Collect results from all threads into a single map.
    ItemMap Aggregate() const;

    void Write(std::ostream& rStream) const;

private:
    Profiler(const Profiler&) = delete;

    Profiler& operator=(const Profiler&) = delete;

private:
    /** @brief Container for storing profiled items persistently.
     *  @details The container storing profiled @ref Item s needs to satisfy two requirements:
     *           @code
     *           - it must be thread-safe
     *           - items must not be moved once created.
     *           @endcode
     *           With these requirements in mind, one solution is a hash table
     *           associating thread ids (initialized with the maximum number of threads
     *           supported by the system upon construction) to linked lists of profiled
     *           @ref Item s. As a result, each thread can construct and insert a thread-local
     *           static item the first time they encounter a specific profile macro, and safely
     *           reference it later on, eliminating the need to go through the expensive process
     *           of finding the related @ref Item in the container again.
     */
    std::unordered_map<std::thread::id,std::list<Item>> mItemContainerMap;

    /// @brief @ref Item for measuring the total lifetime of the @ref Profiler.
    Item mItem;

    /// @brief @ref Scope measuring the total lifetime of the @ref Profiler.
    std::unique_ptr<Scope> mpScope;

    /// @brief Path to the output file to write the results to upon destruction.
    std::filesystem::path mOutputPath;
}; // class Profiler


template <class T>
KRATOS_API(KRATOS_CORE) std::ostream& operator<<(std::ostream& rStream, const Profiler<T>& rProfiler);


template <class TTimeUnit>
class KRATOS_API(KRATOS_CORE) ProfilerSingleton
{
public:
    static Profiler<TTimeUnit>& Get() noexcept;

private:
    static std::optional<Profiler<TTimeUnit>> mProfiler;

    static std::mutex mMutex;
}; // class ProfilerSingleton


} // namespace Kratos::Internals


#if defined(KRATOS_ENABLE_PROFILING)
    #define KRATOS_DEFINE_SCOPE_PROFILER(KRATOS_TIME_UNIT, CODE_LOCATION)                                                     \
        thread_local static auto& KRATOS_STATIC_PROFILER_REF = Kratos::Internals::ProfilerSingleton<KRATOS_TIME_UNIT>::Get(); \
        thread_local static auto& KRATOS_SCOPE_PROFILED_ITEM = KRATOS_STATIC_PROFILER_REF.Create(CODE_LOCATION);              \
        const auto KRATOS_SCOPE_PROFILER = KRATOS_STATIC_PROFILER_REF.Profile(KRATOS_SCOPE_PROFILED_ITEM)

    #define KRATOS_PROFILE_SCOPE_MILLI(CODE_LOCATION) KRATOS_DEFINE_SCOPE_PROFILER(std::chrono::milliseconds, CODE_LOCATION)

    #define KRATOS_PROFILE_SCOPE_MICRO(CODE_LOCATION) KRATOS_DEFINE_SCOPE_PROFILER(std::chrono::microseconds, CODE_LOCATION)

    #define KRATOS_PROFILE_SCOPE_NANO(CODE_LOCATION) KRATOS_DEFINE_SCOPE_PROFILER(std::chrono::nanoseconds, CODE_LOCATION)

    #define KRATOS_PROFILE_SCOPE(CODE_LOCATION) KRATOS_PROFILE_SCOPE_MICRO(CODE_LOCATION)

#else
    #define KRATOS_PROFILE_SCOPE_MILLI(CODE_LOCATION)

    #define KRATOS_PROFILE_SCOPE_MICRO(CODE_LOCATION)

    #define KRATOS_PROFILE_SCOPE_NANO(CODE_LOCATION)

    #define KRATOS_PROFILE_SCOPE(CODE_LOCATION)

#endif


// Definitions of "inlined" functions.
#include "utilities/profiler_impl.h"
