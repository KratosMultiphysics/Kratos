//     __ __      __  __________  ___                ___            __  _           
//    / //_/___ _/ / / /  _/ __ \/   |  ____  ____  / (_)________ _/ /_(_)___  ____ 
//   / ,< / __ `/ /_/ // // /_/ / /| | / __ \/ __ \/ / / ___/ __ `/ __/ / __ \/ __ \
//  / /| / /_/ / __  // // ____/ ___ |/ /_/ / /_/ / / / /__/ /_/ / /_/ / / /_/ / / /
// /_/ |_\__,_/_/ /_/___/_/   /_/  |_/ .___/ .___/_/_/\___/\__,_/\__/_/\____/_/ /_/ 
//                                  /_/   /_/                                       
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes
#include <vector>
#include <string>
#include <limits>

// External includes
#include "kaHIP_interface.h"

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"

namespace Kratos
{

///@addtogroup KaHIPApplication
///@{

///@name Kratos Classes
///@{

/**
 * @class KaHIPPartitioner
 * @ingroup KaHIPApplication
 * @brief Low-level wrapper around the KaHIP @c kaffpa() serial partitioning API.
 * @details Accepts a graph in CSR format (the output of KaHIPCSRConverter) and calls
 *          @c kaffpa() to partition it into k blocks. Supports:
 *
 *          - All six KaHIP preconfiguration modes (FAST, ECO, STRONG, and social variants)
 *          - Configurable imbalance tolerance (fraction, e.g. 0.03 = 3 %)
 *          - Multi-trial mode: runs KaHIP multiple times with different seeds and keeps
 *            the partition with the minimum edge cut
 *          - Optional suppression of KaHIP's stdout output
 *
 *          Thread safety: @c kaffpa() is not thread-safe. If partitioning multiple
 *          graphs concurrently, protect each call with an external mutex.
 *
 * @note Index type: KaHIP's @c kahip_idx is int32_t by default (int64_t if compiled
 *       with -D64BITMODE=ON). Call @c kahip_sizeof_idx() at runtime to check.
 *
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(KAHIP_APPLICATION) KaHIPPartitioner
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of KaHIPPartitioner
    KRATOS_CLASS_POINTER_DEFINITION(KaHIPPartitioner);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Construct from Parameters.
     * @details Accepted keys and defaults:
     * @code{.json}
     * {
     *     "preconfiguration": "eco",       // "fast"|"eco"|"strong"|"fastsocial"|"ecosocial"|"strongsocial"
     *     "imbalance":         0.03,        // allowed imbalance fraction
     *     "seed":              0,           // base random seed
     *     "suppress_output":   true,        // suppress KaHIP stdout
     *     "num_trials":        1            // repeat with seeds 0..num_trials-1; keep best
     * }
     * @endcode
     * @param rSettings  Configuration parameters (validated and assigned defaults)
     */
    explicit KaHIPPartitioner(Parameters rSettings);

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Partition a pre-built CSR graph into @p NumPartitions blocks.
     * @details Calls @c kaffpa() once per trial (see "num_trials"), seeds are
     *          base_seed, base_seed+1, …, base_seed+num_trials-1. Returns the
     *          partition with the minimum edge cut across all trials.
     *
     *          All input arrays must be consistent:
     *          - @p rXAdj has size @p n + 1
     *          - @p rAdjncy has size @p rXAdj[n]
     *          - @p rVWgt has size @p n (may be empty → uniform weights)
     *          - @p rAdjcWgt has size @p rXAdj[n] (may be empty → uniform weights)
     *
     * @param n              Number of graph vertices
     * @param rXAdj          CSR row offsets
     * @param rAdjncy        CSR column indices (0-indexed)
     * @param rVWgt          Vertex weights (empty = uniform)
     * @param rAdjcWgt       Edge weights (empty = uniform)
     * @param NumPartitions  Number of blocks k
     * @return Partition assignment vector of size @p n; part[i] ∈ [0, k)
     */
    std::vector<int> PartitionGraph(
        int n,
        std::vector<kahip_idx>& rXAdj,
        std::vector<kahip_idx>& rAdjncy,
        std::vector<int>& rVWgt,
        std::vector<kahip_idx>& rAdjcWgt,
        int NumPartitions);

    /**
     * @brief Return the KaHIP mode integer corresponding to the current preconfiguration.
     */
    int GetMode() const { return mMode; }

    /**
     * @brief Return the configured imbalance tolerance.
     */
    double GetImbalance() const { return mImbalance; }

    /**
     * @brief Return the default Parameters schema for this class.
     * @details Use with ValidateAndAssignDefaults() to fill partial user settings.
     */
    static Parameters GetDefaultParameters();

    ///@}

private:
    ///@name Member Variables
    ///@{

    /// KaHIP mode integer (FAST=0, ECO=1, STRONG=2, FASTSOCIAL=3, ECOSOCIAL=4, STRONGSOCIAL=5)
    int mMode;

    /// Allowed imbalance fraction (e.g. 0.03 = 3 %)
    double mImbalance;

    /// Base random seed; successive trials use seed, seed+1, seed+2, …
    int mSeed;

    /// When true, kaffpa() does not print to stdout
    bool mSuppressOutput;

    /// Number of independent trials; the best (lowest edge cut) result is returned
    int mNumTrials;

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Convert a preconfiguration string to the KaHIP mode integer.
     * @throws KRATOS_ERROR for unknown preconfiguration strings.
     */
    static int PreconfigurationToMode(const std::string& rPreconfiguration);

    ///@}

}; // class KaHIPPartitioner

///@}

} // namespace Kratos
