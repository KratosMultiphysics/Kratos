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
#include <fstream>
#include <sstream>
#include <atomic>
#include <vector>

namespace Kratos::Internals {

template <class TTimeUnit>
class Profiler
{
private:
    using TimeUnit = TTimeUnit;
    using Duration = TimeUnit;
    using Clock = std::chrono::high_resolution_clock;

    /// @brief Class for identifying a profiled scope and aggregating its stats.
    class Item
    {
    public:
        Item(CodeLocation&& rLocation);
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
        std::size_t operator()(const CodeLocation& r_argument) const
        {
            std::string string(r_argument.GetFileName());
            string.append(std::to_string(r_argument.GetLineNumber()));
            return std::hash<std::string>()(string);
        }
    };

    struct SourceLocationEquality
    {
        bool operator()(const CodeLocation& r_lhs,
                        const CodeLocation& r_rhs) const
        {
            return (std::string(r_lhs.GetFileName()) == std::string(r_rhs.GetFileName())) && 
                   (r_lhs.GetLineNumber() == r_rhs.GetLineNumber());
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

        Item& mrItem;
        const std::chrono::high_resolution_clock::time_point mBegin;
    }; // class Scope

    using ItemMap = std::unordered_map<CodeLocation, Item, SourceLocationHash, SourceLocationEquality>;

public:
    Profiler();
    Profiler(Profiler&& rOther) = default;
    Profiler(std::filesystem::path&& r_outputPath);
    ~Profiler();
    Profiler& operator=(Profiler&& rOther) = default;

    [[nodiscard]] Item& Create(CodeLocation&& rItem);
    [[nodiscard]] Scope Profile(Item& rItem);

    ItemMap Aggregate() const;
    void Write(std::ostream& rStream) const;

private:
    Profiler(const Profiler&) = delete;
    Profiler& operator=(const Profiler&) = delete;

private:
    std::unordered_map<std::thread::id, std::list<Item>> mItemContainerMap;
    Item mItem;
    std::unique_ptr<Scope> mpScope;
    std::filesystem::path mOutputPath;
}; // class Profiler

template <class T>
std::ostream& operator<<(std::ostream& rStream, const Profiler<T>& rProfiler);

template <class TTimeUnit>
class ProfilerSingleton
{
public:
    static Profiler<TTimeUnit>& Get() noexcept;

private:
    static std::optional<Profiler<TTimeUnit>> mProfiler;
    static std::mutex mMutex;
}; // class ProfilerSingleton

template <class TTimeUnit>
std::string GetTimeUnit();

template <>
std::string GetTimeUnit<std::chrono::milliseconds>();

template <>
std::string GetTimeUnit<std::chrono::microseconds>();

template <>
std::string GetTimeUnit<std::chrono::nanoseconds>();

// Definitions of inline functions.
template <class T>
Profiler<T>::Item::Item(CodeLocation&& rLocation)
    : Item(0, Duration(0), Duration(0), Duration(0), std::move(rLocation))
{
}

template <class T>
Profiler<T>::Item::Item(std::size_t CallCount,
                        Duration CumulativeDuration,
                        Duration MinDuration,
                        Duration MaxDuration,
                        CodeLocation&& rLocation)
    : mRecursionLevel(0), mCallCount(CallCount), mCumulative(CumulativeDuration), mMin(MinDuration), mMax(MaxDuration), mLocation(std::move(rLocation))
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
    : mItemContainerMap(), mItem(KRATOS_CODE_LOCATION), mpScope(), mOutputPath(std::move(r_outputPath))
{
    const auto number_of_threads = std::thread::hardware_concurrency();
    std::vector<std::thread> threads;
    threads.reserve(number_of_threads);
    std::atomic<std::size_t> thread_counter = 0;
    for (std::size_t i_thread = 0; i_thread < number_of_threads; ++i_thread) {
        threads.emplace_back([i_thread, &thread_counter, this]() {
            while (thread_counter < i_thread) {}
            mItemContainerMap.emplace(std::this_thread::get_id(), std::list<Item>{});
            ++thread_counter;
        });
    }
    for (auto& r_thread : threads)
        r_thread.join();
    threads.clear();

    mpScope.reset(new Scope(mItem));
}

template <class T>
typename Profiler<T>::Item& Profiler<T>::Create(CodeLocation&& r_item)
{
    auto& r_list = mItemContainerMap[std::this_thread::get_id()];
    r_list.emplace_back(std::move(r_item));
    return r_list.back();
}

template <class T>
typename Profiler<T>::ItemMap Profiler<T>::Aggregate() const
{
    ItemMap output;
    for (const auto& r_threadMapPair : mItemContainerMap) {
        for (const auto& r_item : r_threadMapPair.second) {
            auto it = output.find(r_item.mLocation);
            if (it == output.end()) {
                it = output.emplace(r_item.mLocation, r_item).first;
            } else {
                it->second += r_item;
            }
        }
    }
    return output;
}

template <class T>
void Profiler<T>::Write(std::ostream& rStream) const
{
    Item profiler_item(this->mItem);
    --profiler_item.mCallCount;
    --profiler_item.mRecursionLevel;
    { Scope(profiler_item, this->mpScope->mBegin); }

    auto aggregate_map = this->Aggregate();

    std::vector<const Item*> items;
    items.reserve(aggregate_map.size());
    for (const auto& r_pair : aggregate_map) {
        items.push_back(&r_pair.second);
    }

    std::sort(items.begin(), items.end(), [](const auto& rpLeft, const auto& rpRight) {
        return rpLeft->mCumulative < rpRight->mCumulative;
    });

    Parameters root;
    {
        Parameters object;
        object.AddString("name", "total");
        object.AddString("timeUnit", GetTimeUnit<T>());
        object.AddInt("total", profiler_item.mCumulative.count());
        root.AddValue("meta", std::move(object));
    }

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

} // namespace Kratos::Internals
