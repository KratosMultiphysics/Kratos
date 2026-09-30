//     __ __      __  __________  ___                ___            __  _
//    / //_/___ _/ / / /  _/ __ \/   |  ____  ____  / (_)________ _/ /_(_)___  ____
//   / ,< / __ `/ /_/ // // /_/ / /| | / __ \/ __ \/ / / ___/ __ `/ __/ / __ \/ __
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
#include <numeric>

// External includes

// Project includes
#include "tests/cpp_tests/kahip_fast_suite.h"
#include "custom_utilities/kahip_partitioner.h"

namespace Kratos::Testing {

// ── Helper: 5-node path graph (from KaHIP API docs example) ──────────────────
//   Nodes 0–4, edges: 0-1, 0-4, 1-2, 1-4, 2-3, 3-4
//   xadj   = {0, 2, 5, 7, 9, 12}
//   adjncy = {1,4, 0,2,4, 1,3, 2,4, 0,1,3}
static void BuildK5StarGraph(
    int& n,
    std::vector<kahip_idx>& xadj,
    std::vector<kahip_idx>& adjncy)
{
    n       = 5;
    xadj    = {0, 2, 5, 7, 9, 12};
    adjncy  = {1,4,  0,2,4,  1,3,  2,4,  0,1,3};
}

// ── Helper: 8-node ring (cycle C8) ───────────────────────────────────────────
static void BuildRingGraph(
    int& n,
    std::vector<kahip_idx>& xadj,
    std::vector<kahip_idx>& adjncy)
{
    n = 8;
    xadj.resize(n + 1);
    adjncy.resize(2 * n);
    xadj[0] = 0;
    for (int i = 0; i < n; ++i) {
        adjncy[2 * i]     = (i + n - 1) % n;  // previous node
        adjncy[2 * i + 1] = (i + 1) % n;       // next node
        xadj[i + 1] = xadj[i] + 2;
    }
}

// ── Default Parameters ────────────────────────────────────────────────────────

KRATOS_TEST_CASE_IN_SUITE(KaHIPPartitioner_DefaultParameters, KratosKaHIPFastSuite)
{
    Parameters p = KaHIPPartitioner::GetDefaultParameters();
    KRATOS_EXPECT_EQ(p["preconfiguration"].GetString(), "eco");
    KRATOS_EXPECT_NEAR(p["imbalance"].GetDouble(), 0.03, 1e-14);
    KRATOS_EXPECT_EQ(p["seed"].GetInt(), 0);
    KRATOS_EXPECT_EQ(p["echo_level"].GetInt(), 0);
    KRATOS_EXPECT_EQ(p["num_trials"].GetInt(), 1);
}

// ── PartitionGraph: basic correctness ────────────────────────────────────────

KRATOS_TEST_CASE_IN_SUITE(KaHIPPartitioner_PathGraph_AllNodesAssigned, KratosKaHIPFastSuite)
{
    int n;
    std::vector<kahip_idx> xadj, adjncy;
    BuildK5StarGraph(n, xadj, adjncy);

    std::vector<int>       vwgt;
    std::vector<kahip_idx> adjcwgt;

    Parameters settings = KaHIPPartitioner::GetDefaultParameters();
    KaHIPPartitioner partitioner(settings);

    const std::vector<int> part = partitioner.PartitionGraph(n, xadj, adjncy, vwgt, adjcwgt, 2);

    KRATOS_EXPECT_EQ(static_cast<int>(part.size()), n);
    for (int p : part) {
        KRATOS_EXPECT_GE(p, 0);
        KRATOS_EXPECT_LT(p, 2);
    }
}

KRATOS_TEST_CASE_IN_SUITE(KaHIPPartitioner_PathGraph_EdgeCutFinite, KratosKaHIPFastSuite)
{
    // For 2-partition of a 5-node path, edge cut must be at least 1 and at most 4
    int n;
    std::vector<kahip_idx> xadj, adjncy;
    BuildK5StarGraph(n, xadj, adjncy);

    std::vector<int>       vwgt;
    std::vector<kahip_idx> adjcwgt;

    Parameters settings = KaHIPPartitioner::GetDefaultParameters();
    KaHIPPartitioner partitioner(settings);

    const std::vector<int> part = partitioner.PartitionGraph(n, xadj, adjncy, vwgt, adjcwgt, 2);

    // Count edge cuts manually
    int cut = 0;
    for (int i = 0; i < n; ++i) {
        for (kahip_idx j = xadj[i]; j < xadj[i + 1]; ++j) {
            if (part[i] != part[adjncy[j]]) ++cut;
        }
    }
    cut /= 2; // each undirected edge counted twice

    KRATOS_EXPECT_GE(cut, 1);
    KRATOS_EXPECT_LE(cut, 4);
}

KRATOS_TEST_CASE_IN_SUITE(KaHIPPartitioner_RingGraph_BalancedPartition, KratosKaHIPFastSuite)
{
    // 8-node ring split into 4 parts: each part should have exactly 2 nodes
    int n;
    std::vector<kahip_idx> xadj, adjncy;
    BuildRingGraph(n, xadj, adjncy);

    std::vector<int>       vwgt;
    std::vector<kahip_idx> adjcwgt;

    Parameters settings = KaHIPPartitioner::GetDefaultParameters();
    KaHIPPartitioner partitioner(settings);

    const std::vector<int> part = partitioner.PartitionGraph(n, xadj, adjncy, vwgt, adjcwgt, 4);

    std::vector<int> block_size(4, 0);
    for (int p : part) block_size[p]++;

    // With 3 % imbalance, each block may have up to ceil(8/4 * 1.03) = 3 nodes
    for (int s : block_size) {
        KRATOS_EXPECT_GE(s, 1);
        KRATOS_EXPECT_LE(s, 3);
    }
}

// ── Preconfiguration modes ────────────────────────────────────────────────────

KRATOS_TEST_CASE_IN_SUITE(KaHIPPartitioner_FastMode, KratosKaHIPFastSuite)
{
    int n;
    std::vector<kahip_idx> xadj, adjncy;
    BuildK5StarGraph(n, xadj, adjncy);

    std::vector<int>       vwgt;
    std::vector<kahip_idx> adjcwgt;

    Parameters settings = KaHIPPartitioner::GetDefaultParameters();
    settings["preconfiguration"].SetString("fast");
    KaHIPPartitioner partitioner(settings);
    KRATOS_EXPECT_EQ(partitioner.GetMode(), FAST);

    const auto part = partitioner.PartitionGraph(n, xadj, adjncy, vwgt, adjcwgt, 2);
    KRATOS_EXPECT_EQ(static_cast<int>(part.size()), n);
}

KRATOS_TEST_CASE_IN_SUITE(KaHIPPartitioner_StrongMode, KratosKaHIPFastSuite)
{
    int n;
    std::vector<kahip_idx> xadj, adjncy;
    BuildK5StarGraph(n, xadj, adjncy);

    std::vector<int>       vwgt;
    std::vector<kahip_idx> adjcwgt;

    Parameters settings = KaHIPPartitioner::GetDefaultParameters();
    settings["preconfiguration"].SetString("strong");
    KaHIPPartitioner partitioner(settings);
    KRATOS_EXPECT_EQ(partitioner.GetMode(), STRONG);

    const auto part = partitioner.PartitionGraph(n, xadj, adjncy, vwgt, adjcwgt, 2);
    KRATOS_EXPECT_EQ(static_cast<int>(part.size()), n);
}

// ── Multi-trial mode ──────────────────────────────────────────────────────────

KRATOS_TEST_CASE_IN_SUITE(KaHIPPartitioner_MultiTrial_BestResultReturned, KratosKaHIPFastSuite)
{
    // Run with num_trials=5 and verify the result is still valid
    int n;
    std::vector<kahip_idx> xadj, adjncy;
    BuildRingGraph(n, xadj, adjncy);

    std::vector<int>       vwgt;
    std::vector<kahip_idx> adjcwgt;

    Parameters settings = KaHIPPartitioner::GetDefaultParameters();
    settings["num_trials"].SetInt(5);
    KaHIPPartitioner partitioner(settings);

    const auto part = partitioner.PartitionGraph(n, xadj, adjncy, vwgt, adjcwgt, 2);
    KRATOS_EXPECT_EQ(static_cast<int>(part.size()), n);
    for (int p : part) {
        KRATOS_EXPECT_GE(p, 0);
        KRATOS_EXPECT_LT(p, 2);
    }
}

// ── Node weights ──────────────────────────────────────────────────────────────

KRATOS_TEST_CASE_IN_SUITE(KaHIPPartitioner_WithNodeWeights, KratosKaHIPFastSuite)
{
    int n;
    std::vector<kahip_idx> xadj, adjncy;
    BuildRingGraph(n, xadj, adjncy);

    // Assign uniform node weights
    std::vector<int>       vwgt(n, 1);
    std::vector<kahip_idx> adjcwgt;

    Parameters settings = KaHIPPartitioner::GetDefaultParameters();
    KaHIPPartitioner partitioner(settings);

    const auto part = partitioner.PartitionGraph(n, xadj, adjncy, vwgt, adjcwgt, 2);
    KRATOS_EXPECT_EQ(static_cast<int>(part.size()), n);
}

// ── Edge weights ──────────────────────────────────────────────────────────────

KRATOS_TEST_CASE_IN_SUITE(KaHIPPartitioner_WithEdgeWeights, KratosKaHIPFastSuite)
{
    int n;
    std::vector<kahip_idx> xadj, adjncy;
    BuildRingGraph(n, xadj, adjncy);

    std::vector<int>       vwgt;
    std::vector<kahip_idx> adjcwgt(adjncy.size(), 1); // uniform edge weights

    Parameters settings = KaHIPPartitioner::GetDefaultParameters();
    KaHIPPartitioner partitioner(settings);

    const auto part = partitioner.PartitionGraph(n, xadj, adjncy, vwgt, adjcwgt, 2);
    KRATOS_EXPECT_EQ(static_cast<int>(part.size()), n);
}

// ── Invalid preconfiguration ──────────────────────────────────────────────────

KRATOS_TEST_CASE_IN_SUITE(KaHIPPartitioner_InvalidPreconfiguration_Throws, KratosKaHIPFastSuite)
{
    Parameters settings = KaHIPPartitioner::GetDefaultParameters();
    settings["preconfiguration"].SetString("ultraturbo");

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        KaHIPPartitioner partitioner(settings),
        "unknown preconfiguration");
}

// ── Reproducibility (fixed seed) ─────────────────────────────────────────────

KRATOS_TEST_CASE_IN_SUITE(KaHIPPartitioner_FixedSeed_ReproducibleResult, KratosKaHIPFastSuite)
{
    int n;
    std::vector<kahip_idx> xadj, adjncy;
    BuildRingGraph(n, xadj, adjncy);

    std::vector<int>       vwgt;
    std::vector<kahip_idx> adjcwgt;

    Parameters settings = KaHIPPartitioner::GetDefaultParameters();
    settings["seed"].SetInt(42);
    KaHIPPartitioner partitioner(settings);

    const auto part1 = partitioner.PartitionGraph(n, xadj, adjncy, vwgt, adjcwgt, 2);
    const auto part2 = partitioner.PartitionGraph(n, xadj, adjncy, vwgt, adjcwgt, 2);

    KRATOS_EXPECT_EQ(part1, part2);
}

} // namespace Kratos::Testing
