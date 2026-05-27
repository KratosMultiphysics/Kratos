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

// External includes

// Project includes
#include "tests/cpp_tests/kahip_fast_suite.h"
#include "custom_utilities/kahip_csr_converter.h"

namespace Kratos::Testing {

/**
 * @brief Helper: build a ConnectivitiesContainerType for a path graph 0–1–2–3–4.
 * @details Each node lists its 1-indexed neighbours.
 *          Path:  1-2-3-4-5  (IDs are 1-indexed in Kratos)
 */
static IO::ConnectivitiesContainerType MakePathGraph()
{
    // 5-node path: 1—2—3—4—5
    IO::ConnectivitiesContainerType conn(5);
    conn[0] = {2};           // node 1 → node 2
    conn[1] = {1, 3};        // node 2 → nodes 1, 3
    conn[2] = {2, 4};        // node 3 → nodes 2, 4
    conn[3] = {3, 5};        // node 4 → nodes 3, 5
    conn[4] = {4};           // node 5 → node 4
    return conn;
}

/**
 * @brief Helper: build a complete graph K4 on nodes 1–4.
 */
static IO::ConnectivitiesContainerType MakeK4Graph()
{
    IO::ConnectivitiesContainerType conn(4);
    conn[0] = {2, 3, 4};
    conn[1] = {1, 3, 4};
    conn[2] = {1, 2, 4};
    conn[3] = {1, 2, 3};
    return conn;
}

// ── ConvertToCSRFormat ────────────────────────────────────────────────────────

KRATOS_TEST_CASE_IN_SUITE(KaHIPCSRConverter_PathGraph_XAdjSize, KratosKaHIPFastSuite)
{
    const auto conn = MakePathGraph();
    std::vector<kahip_idx> xadj, adjncy;
    KaHIPCSRConverter::ConvertToCSRFormat(conn, xadj, adjncy);

    // xadj must have n+1 entries
    KRATOS_EXPECT_EQ(static_cast<int>(xadj.size()), 6);
    KRATOS_EXPECT_EQ(xadj[0], 0);
}

KRATOS_TEST_CASE_IN_SUITE(KaHIPCSRConverter_PathGraph_AdjncySize, KratosKaHIPFastSuite)
{
    // Path 1-2-3-4-5: 4 undirected edges → 8 directed entries
    const auto conn = MakePathGraph();
    std::vector<kahip_idx> xadj, adjncy;
    KaHIPCSRConverter::ConvertToCSRFormat(conn, xadj, adjncy);

    KRATOS_EXPECT_EQ(static_cast<int>(adjncy.size()), 8);
    KRATOS_EXPECT_EQ(xadj[5], 8);
}

KRATOS_TEST_CASE_IN_SUITE(KaHIPCSRConverter_PathGraph_ZeroIndexed, KratosKaHIPFastSuite)
{
    // Kratos uses 1-indexed IDs; CSR output must be 0-indexed
    const auto conn = MakePathGraph();
    std::vector<kahip_idx> xadj, adjncy;
    KaHIPCSRConverter::ConvertToCSRFormat(conn, xadj, adjncy);

    // Node 0 (Kratos ID 1) has a single neighbour: node 1 (Kratos ID 2)
    KRATOS_EXPECT_EQ(xadj[0], 0);
    KRATOS_EXPECT_EQ(xadj[1], 1);
    KRATOS_EXPECT_EQ(adjncy[0], 1); // 0-indexed ID of Kratos node 2
}

KRATOS_TEST_CASE_IN_SUITE(KaHIPCSRConverter_PathGraph_MonotoneXAdj, KratosKaHIPFastSuite)
{
    const auto conn = MakePathGraph();
    std::vector<kahip_idx> xadj, adjncy;
    KaHIPCSRConverter::ConvertToCSRFormat(conn, xadj, adjncy);

    for (std::size_t i = 1; i < xadj.size(); ++i) {
        KRATOS_EXPECT_GE(xadj[i], xadj[i - 1]);
    }
}

// ── ConvertToCSRFormatWithWeights ─────────────────────────────────────────────

KRATOS_TEST_CASE_IN_SUITE(KaHIPCSRConverter_UniformWeights, KratosKaHIPFastSuite)
{
    const auto conn = MakeK4Graph();
    std::vector<kahip_idx> xadj, adjncy, adjcwgt;
    std::vector<int> vwgt;
    KaHIPCSRConverter::ConvertToCSRFormatWithWeights(conn, xadj, adjncy, vwgt, adjcwgt);

    // 4 nodes → vwgt of size 4, all 1
    KRATOS_EXPECT_EQ(static_cast<int>(vwgt.size()), 4);
    for (int w : vwgt) KRATOS_EXPECT_EQ(w, 1);

    // K4 has 6 undirected edges → 12 directed entries → adjcwgt of size 12, all 1
    KRATOS_EXPECT_EQ(static_cast<int>(adjcwgt.size()), 12);
    for (auto w : adjcwgt) KRATOS_EXPECT_EQ(w, 1);
}

// ── ValidateCSRGraph ──────────────────────────────────────────────────────────

KRATOS_TEST_CASE_IN_SUITE(KaHIPCSRConverter_ValidatePathGraph, KratosKaHIPFastSuite)
{
    const auto conn = MakePathGraph();
    std::vector<kahip_idx> xadj, adjncy;
    KaHIPCSRConverter::ConvertToCSRFormat(conn, xadj, adjncy);

    KRATOS_EXPECT_TRUE(KaHIPCSRConverter::ValidateCSRGraph(5, xadj, adjncy));
}

KRATOS_TEST_CASE_IN_SUITE(KaHIPCSRConverter_ValidateK4Graph, KratosKaHIPFastSuite)
{
    const auto conn = MakeK4Graph();
    std::vector<kahip_idx> xadj, adjncy;
    KaHIPCSRConverter::ConvertToCSRFormat(conn, xadj, adjncy);

    KRATOS_EXPECT_TRUE(KaHIPCSRConverter::ValidateCSRGraph(4, xadj, adjncy));
}

KRATOS_TEST_CASE_IN_SUITE(KaHIPCSRConverter_ValidateRejectsAsymmetric, KratosKaHIPFastSuite)
{
    // Manually construct an asymmetric graph (edge 0→1 but not 1→0)
    // n=2, only one directed edge
    std::vector<kahip_idx> xadj   = {0, 1, 1};
    std::vector<kahip_idx> adjncy = {1};

    // ValidateCSRGraph should detect the missing reverse edge and return false
    KRATOS_EXPECT_FALSE(KaHIPCSRConverter::ValidateCSRGraph(2, xadj, adjncy));
}

KRATOS_TEST_CASE_IN_SUITE(KaHIPCSRConverter_ValidateRejectsSelfLoop, KratosKaHIPFastSuite)
{
    // n=3, node 1 has a self-loop: adjncy entry = 1 for node 1
    std::vector<kahip_idx> xadj   = {0, 1, 3, 4};
    // node 0 → node 2, node 1 → self (1) and node 2, node 2 → node 1
    std::vector<kahip_idx> adjncy = {2, 1, 2, 1};

    KRATOS_EXPECT_FALSE(KaHIPCSRConverter::ValidateCSRGraph(3, xadj, adjncy));
}

} // namespace Kratos::Testing
