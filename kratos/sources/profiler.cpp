// --- Core Includes ---
#include "includes/kratos_parameters.h"
#include "utilities/profiler.h"

// --- STL Includes ---
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
Profiler<T>::Item::Item(CodeLocation&& r_location)
    : mRecursionLevel(0),
      mCallCount(0),
      mTime(0),
      mLocation(std::move(r_location))
{
}


template <class T>
typename Profiler<T>::Item& Profiler<T>::Item::operator+=(const Item& r_rhs)
{
    mCallCount += r_rhs.mCallCount;
    mTime += r_rhs.mTime;
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
      mpScope(),
      mOutputPath(std::move(r_outputPath))
{
    // "Reserve" thread map to avoid bucket moving later on.
    const auto numberOfThreads = std::thread::hardware_concurrency();
    std::vector<std::thread> threads;
    threads.reserve(numberOfThreads);
    std::atomic<int> thread_counter = 0;
    for (int i_thread=0; i_thread<numberOfThreads; ++i_thread) {
        threads.emplace_back([i_thread, &thread_counter, this](){
            while (thread_counter < i_thread) {} // <== wait until the previous thread finishes
            mItemContainerMap.emplace(std::this_thread::get_id(), std::list<Item> {});
            ++thread_counter;
        });
    }
    for (auto& r_thread : threads)
        r_thread.join();
    threads.clear();

    // Insert first item measuring the total runtime
    auto& r_item = this->Create(CodeLocation());
    mpScope.reset(new Scope(r_item)); // <== no, this doesn't immediately become a dangling reference
}


template <class T>
typename Profiler<T>::Item& Profiler<T>::Create(CodeLocation&& r_item)
{
    auto& r_list = mItemContainerMap[std::this_thread::get_id()];
    r_list.emplace_back(std::move(r_item));
    return r_list.back();
}


struct SourceLocationHash
{
    std::size_t operator()(const CodeLocation& r_argument) const
    {
        std::string string(r_argument.GetFileName());
        string.append(std::to_string(r_argument.GetLineNumber()));
        return std::hash<std::string>()(string);
    }
};


struct SourceLocationEquals
{
    bool operator()(const CodeLocation& r_lhs,
                    const CodeLocation& r_rhs) const
    {
        return (std::string(r_lhs.GetFileName()) == std::string(r_rhs.GetFileName())) && (r_lhs.GetLineNumber() == r_rhs.GetLineNumber());
    }
};


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
Profiler<T>::~Profiler()
{
    // Stop the total timer
    mpScope.reset();
    const auto totalTime = mItemContainerMap[std::this_thread::get_id()].front().mTime;
    mItemContainerMap[std::this_thread::get_id()].pop_front();

    std::unordered_map<
        CodeLocation,
        Item,
        SourceLocationHash,
        SourceLocationEquals
    > aggregateMap;

    // Combine maps from all threads
    for (const auto& r_threadMapPair : mItemContainerMap)
    {
        for (const auto& r_item : r_threadMapPair.second)
        {
            auto it = aggregateMap.find(r_item.mLocation);
            if (it == aggregateMap.end())
                it = aggregateMap.emplace(r_item.mLocation, r_item).first;
            else
                it->second += r_item;
        } // for item in vector
    } // for (location, vector) in map

    // Sort items based on their total duration
    using Numeric = typename Duration::rep;
    std::vector<std::tuple<CodeLocation,std::size_t,Numeric>> items;
    items.reserve(aggregateMap.size());
    for (const auto& r_pair : aggregateMap)
        items.emplace_back(r_pair.first, r_pair.second.mCallCount, r_pair.second.mTime.count());
    std::sort(items.begin(),
              items.end(),
              [](const auto& r_lhs, const auto& r_rhs)
                {return std::get<Numeric>(r_lhs) < std::get<Numeric>(r_rhs);});

    Parameters root;

    // Add metadata
    {
        Parameters object;
        object.AddString("name", "total");
        object.AddString("timeUnit", GetTimeUnit<T>());
        object.AddInt("time", totalTime.count());
        root.AddValue("meta", std::move(object));
    }

    // Add all items
    root.AddEmptyArray("results");
    Parameters results = root["results"];
    for (const auto& r_item : items)
    {
        Parameters result;
        const auto& r_location = std::get<0>(r_item);
        result.AddString("file", std::string(r_location.GetFileName()));
        result.AddInt("line", int(r_location.GetLineNumber()));
        result.AddString("function", std::string(std::get<0>(r_item).GetFunctionName()));
        result.AddInt("callCount", std::get<1>(r_item));

        std::stringstream stream;
        stream << std::get<2>(r_item);
        result.AddString("time", stream.str());
        results.Append(result);
    }

    std::ofstream file(mOutputPath);
    file << root.WriteJsonString() << std::endl;
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


template class Profiler<std::chrono::microseconds>;
template class ProfilerSingleton<std::chrono::microseconds>;


template class Profiler<std::chrono::nanoseconds>;
template class ProfilerSingleton<std::chrono::nanoseconds>;


} // namespace cie::utils
