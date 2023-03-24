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

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "utilities/profiler.h"

// System includes
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <thread>
#include <vector>
#include <tuple>
#include <fstream>
#include <sstream>
#include <atomic>


namespace Kratos::Internals {


template <class T>
Profiler<T>::Item::Item(CodeLocation&& rLocation)
    : Item(0,
           Duration(0),
           Duration(0),
           Duration(0),
           std::move(rLocation))
{
}


template <class T>
Profiler<T>::Item::Item(std::size_t CallCount,
                        Duration CumulativeDuration,
                        Duration MinDuration,
                        Duration MaxDuration,
                        CodeLocation&& rLocation)
    : mRecursionLevel(0),
      mCallCount(CallCount),
      mCumulative(CumulativeDuration),
      mMin(MinDuration),
      mMax(MaxDuration),
      mLocation(std::move(rLocation))
{
}


template <class T>
typename Profiler<T>::Item& Profiler<T>::Item::operator+=(const Item& rOther)
{
    mCallCount += rOther.mCallCount;
    mCumulative += rOther.mCumulative;
    mMin = std::min(mMin, rOther.mMin);
    mMax = std::max(mMax, rOther.mMax);
    return *this;
}


template <class T>
Profiler<T>::Profiler()
    : Profiler("kratos_profiler_output.json")
{
}


template <class T>
Profiler<T>::Profiler(std::filesystem::path&& r_outputPath)
    : mItemContainerMap(),
      mItem(KRATOS_CODE_LOCATION),
      mpScope(),
      mOutputPath(std::move(r_outputPath))
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


template <class T>
typename Profiler<T>::Item& Profiler<T>::Create(CodeLocation&& r_item)
{
    auto& r_list = mItemContainerMap[std::this_thread::get_id()];
    r_list.emplace_back(std::move(r_item));
    return r_list.back();
}


namespace {
template <class TTimeUnit>
std::string GetTimeUnit()
{
    KRATOS_ERROR << "Unknown time unit";
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


template <class T>
typename Profiler<T>::ItemMap Profiler<T>::Aggregate() const
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


template <class T>
void Profiler<T>::Write(std::ostream& rStream) const
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
        result.AddString("function", std::string(r_location.GetFunctionName()));
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


template <class T>
std::ostream& operator<<(std::ostream& rStream, const Profiler<T>& rProfiler)
{
    rProfiler.Write(rStream);
    return rStream;
}


template <class T>
Profiler<T>::~Profiler()
{
    std::ofstream file(mOutputPath);
    file << *this;
}


template <class T>
Profiler<T>& ProfilerSingleton<T>::Get() noexcept
{
    std::scoped_lock<std::mutex> lock(mMutex);
    if (!mProfiler.has_value())
        mProfiler.emplace();

    return mProfiler.value();
}


template <class T>
std::optional<Profiler<T>> ProfilerSingleton<T>::mProfiler;


template <class T>
std::mutex ProfilerSingleton<T>::mMutex;


template class Profiler<std::chrono::milliseconds>;
template class ProfilerSingleton<std::chrono::milliseconds>;
template std::ostream& operator<<(std::ostream&, const Profiler<std::chrono::milliseconds>&);


template class Profiler<std::chrono::microseconds>;
template class ProfilerSingleton<std::chrono::microseconds>;
template std::ostream& operator<<(std::ostream&, const Profiler<std::chrono::microseconds>&);


template class Profiler<std::chrono::nanoseconds>;
template class ProfilerSingleton<std::chrono::nanoseconds>;
template std::ostream& operator<<(std::ostream&, const Profiler<std::chrono::nanoseconds>&);


} // namespace cie::utils
