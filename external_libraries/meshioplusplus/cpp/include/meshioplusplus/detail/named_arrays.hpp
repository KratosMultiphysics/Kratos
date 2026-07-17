//  ██████   ██████ ██████████  █████████  █████   █████ █████    ███████
// ░░██████ ██████ ░░███░░░░░█ ███░░░░░███░░███   ░░███ ░░███   ███░░░░░███      ███         ███
//  ░███░█████░███  ░███  █ ░ ░███    ░░░  ░███    ░███  ░███  ███     ░░███    ░███        ░███
//  ░███░░███ ░███  ░██████   ░░█████████  ░███████████  ░███ ░███      ░███ ███████████ ███████████
//  ░███ ░░░  ░███  ░███░░█    ░░░░░░░░███ ░███░░░░░███  ░███ ░███      ░███░░░░░███░░░ ░░░░░███░░░
//  ░███      ░███  ░███ ░   █ ███    ░███ ░███    ░███  ░███ ░░███     ███     ░███        ░███
//  █████     █████ ██████████░░█████████  █████   █████ █████ ░░░███████░      ░░░         ░░░
// ░░░░░     ░░░░░ ░░░░░░░░░░  ░░░░░░░░░  ░░░░░   ░░░░░ ░░░░░    ░░░░░░░
//
//
//  License:         MIT License
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//
#pragma once

/**
 * @file named_arrays.hpp
 * @brief Insertion-ordered `name -> NDArray` (and `name -> vector<NDArray>`)
 * containers with O(1) name lookup.
 *
 * Storage backbone for the NATIVE mesh backend's point/cell/field data (and
 * reused by the KRATOS backend's field data): a contiguous
 * `std::vector<std::pair<name, value>>` preserving insertion order plus an
 * `std::unordered_map<name, index>` kept in sync for O(1) lookup — the same
 * vector-plus-access-map pattern Kratos's CoSimIO uses for its entity
 * containers. The uniform mesh API's observable name order is *sorted*
 * (see `mesh_api.hpp`), which `SortedNames()` provides regardless of
 * insertion order — matching what `detail::sorted_keys` produces for the
 * MESHIO backend's `unordered_map` storage.
 */

// System includes
#include <algorithm>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/ndarray.hpp"

namespace meshioplusplus {
namespace detail {

/** @brief Insertion-ordered `name -> T` store with O(1) lookup by name. */
template <class T>
class NamedItems {
public:
    /** @brief Insert-or-assign `rValue` under `rName`. */
    void Set(std::string name, T value) {
        auto it = mIndex.find(name);
        if (it != mIndex.end()) {
            mItems[it->second].second = std::move(value);
            return;
        }
        mIndex.emplace(name, mItems.size());
        mItems.emplace_back(std::move(name), std::move(value));
    }
    /** @brief Whether an entry named @p rName exists. */
    bool Has(const std::string& rName) const { return mIndex.count(rName) > 0; }
    /** @brief The entry named @p rName; throws `std::out_of_range` if absent. */
    const T& Get(const std::string& rName) const {
        auto it = mIndex.find(rName);
        if (it == mIndex.end())
            throw std::out_of_range("meshio++: no data array named '" + rName + "'");
        return mItems[it->second].second;
    }
    /** @brief Mutable access to the entry named @p rName, created if absent. */
    T& GetOrCreate(const std::string& rName) {
        auto it = mIndex.find(rName);
        if (it != mIndex.end())
            return mItems[it->second].second;
        mIndex.emplace(rName, mItems.size());
        mItems.emplace_back(rName, T{});
        return mItems.back().second;
    }
    /** @brief Number of entries. */
    std::size_t Size() const { return mItems.size(); }
    /** @brief All names, sorted ascending (the API's observable order). */
    std::vector<std::string> SortedNames() const {
        std::vector<std::string> names;
        names.reserve(mItems.size());
        for (const auto& kv : mItems)
            names.push_back(kv.first);
        std::sort(names.begin(), names.end());
        return names;
    }
    /** @brief The underlying insertion-ordered items (for direct iteration). */
    const std::vector<std::pair<std::string, T>>& Items() const { return mItems; }

private:
    std::vector<std::pair<std::string, T>> mItems;
    std::unordered_map<std::string, std::size_t> mIndex;
};

using NamedArrays = NamedItems<NDArray>;
using NamedArrayLists = NamedItems<std::vector<NDArray>>;

}  // namespace detail
}  // namespace meshioplusplus
