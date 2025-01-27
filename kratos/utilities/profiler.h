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
#include "includes/kratos_parameters.h"

// System includes
#include <unordered_map>
#include <filesystem>
#include <chrono>
#include <thread>
#include <optional>
#include <list>
#include <mutex>
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <vector>
#include <fstream>
#include <sstream>
#include <atomic>
#include <limits>

namespace Kratos::Internals {

namespace {
template <class TTimeUnit>
std::string GetTimeUnit()
{
    KRATOS_ERROR << "Unsupported time unit";
}

template <>
std::string GetTimeUnit<std::chrono::milliseconds>()
{
    return "ms";
}

template <>
std::string GetTimeUnit<std::chrono::microseconds>()
{
    return "us";
}

template <>
std::string GetTimeUnit<std::chrono::nanoseconds>()
{
    return "ns";
}
} // unnamed namespace

template <class TTimeUnit>
class Profiler
{
private:
    /// @brief Absolute time type.
    using TimeUnit = TTimeUnit;

    /// @brief Relative time type.
    using Duration = TimeUnit;

    /// @brief  Clock type used for measuring durations.
    using Clock = std::chrono::high_resolution_clock;

    /// @brief Class for identifying a profiled scope and aggregating its stats.
    class Item
    {
    public:
        explicit Item(CodeLocation&& rLocation)
            : Item(0,                                                               // <== .mCallCount
                  Duration(0),                                                     // <== .mCumulative
                  Duration(std::numeric_limits<typename Duration::rep>::max()),    // <== .mMin
                  Duration(0),                                                     // <== .mMax
                  std::move(rLocation))                                            // <== .mLocation
        {
        }

    private:
        Item(std::size_t CallCount,
             Duration CumulativeDuration,
             Duration MinDuration,
             Duration MaxDuration,
             CodeLocation&& rLocation)
            : mRecursionLevel(0)
            , mCallCount(CallCount)
            , mCumulative(CumulativeDuration)
            , mMin(MinDuration)
            , mMax(MaxDuration)
            , mLocation(std::move(rLocation))
        {
        }

        /// @brief Aggregate profiled data from another @ref Item in the same scope.
        Item& operator+=(const Item& rOther)
        {
            mCallCount += rOther.mCallCount;
            mCumulative += rOther.mCumulative;
            mMin = std::min(mMin, rOther.mMin);
            mMax = std::max(mMax, rOther.mMax);
            return *this;
        }

    private:
        friend class Profiler;

        /// @brief Counter for keeping track of recursive calls.
        /// @details Recursive function calls are aggregated onto the top
        ///          level call. To do that, the @ref Item must keep track
        ///          of its recursion depth.
        unsigned mRecursionLevel;

        /// @brief Counter tracking total number of calls to a function during the program's entire execution time.
        std::size_t mCallCount;

        /// @brief Counter summing the duration of each call to the profiled scope.
        Duration mCumulative;

        /// @brief Minimum time spent in the profiled scope.
        Duration mMin;

        /// @brief Maximum time spent in the profiled scope.
        Duration mMax;

        /// @brief Source information about the profiled scope.
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
        ~Scope()
        {
            const auto end = Clock::now();
            const auto duration = std::chrono::duration_cast<Duration>(end - mBegin);
            if (0 == mrItem.mRecursionLevel) {
                ++mrItem.mCallCount;
                mrItem.mCumulative += duration;
                mrItem.mMin = std::min(mrItem.mMin, duration);
                mrItem.mMax = std::max(mrItem.mMax, duration);
            }
            --mrItem.mRecursionLevel;
        }

    private:
        Scope(Item& rItem)
            : mrItem(rItem)
            , mBegin(Clock::now())
        {
            ++mrItem.mRecursionLevel;
        }

        Scope(Item& rItem, std::chrono::high_resolution_clock::time_point Begin)
            : mrItem(rItem)
            , mBegin(Begin)
        {
            ++mrItem.mRecursionLevel;
        }

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
    Profiler()
        : Profiler("kratos_profiler_output_" + GetTimeUnit<T>() + ".json")
    {
    }

    Profiler(std::filesystem::path&& r_outputPath)
        : mItemContainerMap()
        , mItem(KRATOS_CODE_LOCATION)
        , mpScope()
        , mOutputPath(std::move(r_outputPath))
    {
        // "Reserve" thread map to avoid bucket moving later on.
        const auto number_of_threads = std::thread::hardware_concurrency();
        std::vector<std::thread> threads;
        threads.reserve(number_of_threads);
        std::atomic<std::size_t> thread_counter = 0;
        for (std::size_t i_thread=0; i_thread<number_of_threads; ++i_thread) {
            threads.emplace_back([i_thread, &thread_counter, this](){
                while (thread_counter < i_thread) {} // <== wait until the previous thread finishes
                mItemContainerMap.emplace(std::this_thread::get_id(), std::list<Item> {});
                ++thread_counter;
            });
        }
        for (auto& r_thread : threads)
            r_thread.join();
        threads.clear();

        // Measure the total lifetime of the profiler
        mpScope.reset(new Scope(mItem));
    }

    Profiler(Profiler&& rOther) = default;

    ~Profiler()
    {
        std::ofstream file(mOutputPath);
        file << *this;
    }

    Profiler& operator=(Profiler&& rOther) = default;

    [[nodiscard]] Item& Create(CodeLocation&& rItem)
    {
        auto& r_list = mItemContainerMap[std::this_thread::get_id()];
        r_list.emplace_back(std::move(rItem));
        return r_list.back();
    }

    [[nodiscard]] Scope Profile(Item& rItem)
    {
        return Scope(rItem);
    }

    /// @brief Collect results from all threads into a single map.
    ItemMap Aggregate() const
    {
        KRATOS_TRY

        // Aggregate maps from all threads
        ItemMap output;
        for (const auto& r_threadMapPair : mItemContainerMap) {
            for (const auto& r_item : r_threadMapPair.second) {
                auto it = output.find(r_item.mLocation);
                if (it == output.end()) {
                    it = output.emplace(r_item.mLocation, r_item).first;
                } else {
                    it->second += r_item;
                }
            } // for item in vector
        } // for (location, vector) in map

        return output;

        KRATOS_CATCH("")
    }

    void Write(std::ostream& rStream) const
    {
        KRATOS_TRY

        // Time the profiler's scope without changing its state
        Item profiler_item(this->mItem);
        --profiler_item.mCallCount;
        --profiler_item.mRecursionLevel;
        {Scope(profiler_item, this->mpScope->mBegin);} // Force update the copied item

        auto aggregate_map = this->Aggregate();

        // Sort items based on their total duration
        std::vector<const Item*> items;
        items.reserve(aggregate_map.size());
        for (const auto& r_pair : aggregate_map) {
            items.push_back(&r_pair.second);
        }

        std::sort(items.begin(),
                  items.end(),
                  [](const auto& rpLeft, const auto& rpRight)
                    {return rpLeft->mCumulative < rpRight->mCumulative;});

        // Start with metadata
        Parameters root;
        {
            Parameters object;
            object.AddString("name", "total");
            object.AddString("timeUnit", GetTimeUnit<T>());
            object.AddInt("total", profiler_item.mCumulative.count());
            root.AddValue("meta", std::move(object));
        }

        // Add all items
        root.AddEmptyArray("results");
        Parameters results = root["results"];
        for (const typename Profiler<T>::Item* p_item : items) {
            Parameters result;
            const auto& r_location = p_item->mLocation;
            result.AddString("file", std::string(r_location.GetFileName()));
            result.AddInt("line", int(r_location.GetLineNumber()));
            result.AddString("signature", std::string(r_location.GetFunctionName()));
            result.AddString("function", std::string(r_location.CleanFunctionName()));
            result.AddInt("callCount", p_item->mCallCount);

            std::stringstream stream;
            stream << std::chrono::duration_cast<TimeUnit>(p_item->mCumulative).count();
            result.AddString("total", stream.str());

            stream.str("");
            stream << std::chrono::duration_cast<TimeUnit>(p_item->mMin).count();
            result.AddString("min", stream.str());

            stream.str("");
            stream << std::chrono::duration_cast<TimeUnit>(p_item->mMax).count();
            result.AddString("max", stream.str());

            results.Append(result);
        }

        rStream << root.PrettyPrintJsonString() << std::endl;

        KRATOS_CATCH("")
    }

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
std::ostream& operator<<(std::ostream& rStream, const Profiler<T>& rProfiler)
{
    rProfiler.Write(rStream);
    return rStream;
}

template <class TTimeUnit>
class ProfilerSingleton
{
public:
    static Profiler<TTimeUnit>& Get() noexcept
    {
        std::scoped_lock<std::mutex> lock(mMutex);
        if (!mProfiler.has_value())
            mProfiler.emplace();

        return mProfiler.value();
    }

private:
    static std::optional<Profiler<TTimeUnit>> mProfiler;
    static std::mutex mMutex;
}; // class ProfilerSingleton

template <class T>
std::optional<Profiler<T>> ProfilerSingleton<T>::mProfiler;

template <class T>
std::mutex ProfilerSingleton<T>::mMutex;

// Template instantiations
template class Profiler<std::chrono::milliseconds>;
template class ProfilerSingleton<std::chrono::milliseconds>;
template std::ostream& operator<<(std::ostream&, const Profiler<std::chrono::milliseconds>&);

template class Profiler<std::chrono::microseconds>;
template class ProfilerSingleton<std::chrono::microseconds>;
template std::ostream& operator<<(std::ostream&, const Profiler<std::chrono::microseconds>&);

template class Profiler<std::chrono::nanoseconds>;
template class ProfilerSingleton<std::chrono::nanoseconds>;
template std::ostream& operator<<(std::ostream&, const Profiler<std::chrono::nanoseconds>&);

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
