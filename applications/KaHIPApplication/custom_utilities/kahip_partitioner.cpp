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

// System includes
#include <algorithm>
#include <limits>

// External includes
#include "kaHIP_interface.h"

// Project includes
#include "custom_utilities/kahip_partitioner.h"
#include "includes/exception.h"

namespace Kratos 
{

KaHIPPartitioner::KaHIPPartitioner(Parameters rSettings)
{
    rSettings.ValidateAndAssignDefaults(GetDefaultParameters());

    mMode           = PreconfigurationToMode(rSettings["preconfiguration"].GetString());
    mImbalance      = rSettings["imbalance"].GetDouble();
    mSeed           = rSettings["seed"].GetInt();
    mSuppressOutput = rSettings["suppress_output"].GetBool();
    mNumTrials      = rSettings["num_trials"].GetInt();

    KRATOS_ERROR_IF(mImbalance <= 0.0 || mImbalance >= 1.0)
        << "KaHIPPartitioner: imbalance must be in (0, 1), got " << mImbalance << std::endl;

    KRATOS_ERROR_IF(mNumTrials < 1)
        << "KaHIPPartitioner: num_trials must be >= 1, got " << mNumTrials << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

std::vector<int> KaHIPPartitioner::PartitionGraph(
    const int NGraphVertices,
    std::vector<kahip_idx>& rXAdj,
    std::vector<kahip_idx>& rAdjncy,
    std::vector<int>& rVWgt,
    std::vector<kahip_idx>& rAdjcWgt,
    const int NumPartitions
    )
{
    KRATOS_ERROR_IF(NGraphVertices <= 0)
        << "KaHIPPartitioner: number of vertices n must be > 0, got " << NGraphVertices << std::endl;

    KRATOS_ERROR_IF(NumPartitions < 1)
        << "KaHIPPartitioner: NumPartitions must be >= 1, got " << NumPartitions << std::endl;

    KRATOS_ERROR_IF(static_cast<int>(rXAdj.size()) != NGraphVertices + 1)
        << "KaHIPPartitioner: xadj size " << rXAdj.size()
        << " != n+1 = " << (NGraphVertices + 1) << std::endl;

    // kaffpa takes non-const pointers; make local copies of the const parameters
    int n_vertices   = NGraphVertices;
    int n_parts      = NumPartitions;

    // Pointers for optional weight arrays (nullptr triggers KaHIP uniform weighting)
    int*        p_vwgt    = rVWgt.empty()    ? nullptr : rVWgt.data();
    kahip_idx*  p_adjcwgt = rAdjcWgt.empty() ? nullptr : rAdjcWgt.data();

    // Storage for the best result across all trials
    std::vector<int> best_part(NGraphVertices);
    kahip_idx        best_cut = std::numeric_limits<kahip_idx>::max();

    std::vector<int> trial_part(NGraphVertices);
    double imbalance = mImbalance; // kaffpa takes a pointer, must be non-const

    for (int trial = 0; trial < mNumTrials; ++trial) {
        const int seed = mSeed + trial;
        kahip_idx edgecut = 0;

        kaffpa(
            &n_vertices,
            p_vwgt,
            rXAdj.data(),
            p_adjcwgt,
            rAdjncy.data(),
            &n_parts,
            &imbalance,
            mSuppressOutput,
            seed,
            mMode,
            &edgecut,
            trial_part.data());

        KRATOS_INFO_IF("KaHIPPartitioner", !mSuppressOutput)
            << "  Trial " << trial << " (seed=" << seed
            << "): edge cut = " << edgecut << std::endl;

        if (edgecut < best_cut) {
            best_cut = edgecut;
            best_part = trial_part;
        }
    }

    KRATOS_INFO("KaHIPPartitioner")
        << "Partitioned graph (" << NGraphVertices << " nodes, " << NumPartitions
        << " parts): best edge cut = " << best_cut
        << " (over " << mNumTrials << " trial(s))" << std::endl;

    return best_part;
}

/***********************************************************************************/
/***********************************************************************************/

Parameters KaHIPPartitioner::GetDefaultParameters()
{
    Parameters defaults(R"({
        "preconfiguration": "eco",
        "imbalance":         0.03,
        "seed":              0,
        "suppress_output":   true,
        "num_trials":        1
    })");
    return defaults;
}

/***********************************************************************************/
/***********************************************************************************/

int KaHIPPartitioner::PreconfigurationToMode(const std::string& rPreconfiguration)
{
    if (rPreconfiguration == "fast")         return FAST;
    if (rPreconfiguration == "eco")          return ECO;
    if (rPreconfiguration == "strong")       return STRONG;
    if (rPreconfiguration == "fastsocial")   return FASTSOCIAL;
    if (rPreconfiguration == "ecosocial")    return ECOSOCIAL;
    if (rPreconfiguration == "strongsocial") return STRONGSOCIAL;

    KRATOS_ERROR << "KaHIPPartitioner: unknown preconfiguration '" << rPreconfiguration
                 << "'. Valid values: fast, eco, strong, fastsocial, ecosocial, strongsocial." << std::endl;
}

} // namespace Kratos
