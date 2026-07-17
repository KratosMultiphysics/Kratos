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
 * @file map_order.hpp
 * @brief Deterministic iteration order for the `Mesh` data maps.
 *
 * `Mesh::point_data`/`cell_data`/`field_data` are `std::unordered_map` (O(1)
 * lookup, no ordering guarantee), but their key order is observable — it drives
 * Python dict key order and the on-disk field/variable order of several writers
 * (VTU, XDMF, Exodus, Tecplot, HMF) as well as medit's "first int field"
 * selection. `sorted_keys` recovers that order explicitly at each such
 * consumption site, decoupling "how we store" from "how we emit" so output stays
 * byte-identical regardless of the storage container.
 */

// System includes
#include <algorithm>
#include <vector>

namespace meshioplusplus {
namespace detail {

/**
 * @brief Collect a map's keys in sorted order.
 * @tparam Map An associative container (ordered or unordered).
 * @param m The map whose keys to enumerate.
 * @return The keys of @p m sorted ascending.
 */
template <class Map>
std::vector<typename Map::key_type> sorted_keys(const Map& rM) {
    std::vector<typename Map::key_type> keys;
    keys.reserve(rM.size());
    for (const auto& kv : rM)
        keys.push_back(kv.first);
    std::sort(keys.begin(), keys.end());
    return keys;
}

}  // namespace detail
}  // namespace meshioplusplus
