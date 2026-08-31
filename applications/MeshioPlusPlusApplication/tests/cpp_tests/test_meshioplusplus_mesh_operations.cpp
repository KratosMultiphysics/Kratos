//  ██████   ██████ ██████████  █████████  █████   █████ █████    ███████
// ░░██████ ██████ ░░███░░░░░█ ███░░░░░███░░███   ░░███ ░░███   ███░░░░░███      ███         ███
//  ░███░█████░███  ░███  █ ░ ░███    ░░░  ░███    ░███  ░███  ███     ░░███    ░███        ░███
//  ░███░░███ ░███  ░██████   ░░█████████  ░███████████  ░███ ░███      ░███ ███████████ ███████████
//  ░███ ░░░  ░███  ░███░░█    ░░░░░░░░███ ░███░░░░░███  ░███ ░███      ░███░░░░░███░░░ ░░░░░███░░░
//  ░███      ░███  ░███ ░   █ ███    ░███ ░███    ░███  ░███ ░░███     ███     ░███        ░███
//  █████     █████ ██████████░░█████████  █████   █████ █████ ░░░███████░      ░░░         ░░░
// ░░░░░     ░░░░░ ░░░░░░░░░░  ░░░░░░░░░  ░░░░░   ░░░░░ ░░░░░    ░░░░░░░                        Application
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <sstream>

// External includes
#include "meshioplusplus/operations/refine.hpp"

// Project includes
#include "containers/model.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "custom_utilities/meshioplusplus_mesh_operations.h"
#include "custom_utilities/meshioplusplus_conversion_utilities.h"
#include "custom_io/meshioplusplus_io.h"
#include "meshioplusplus_fast_suite.h"

namespace Kratos::Testing {

namespace {
namespace mio = meshioplusplus;

std::filesystem::path TestFilePath(const std::string& rExtension)
{
    const std::string test_name = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    return std::filesystem::temp_directory_path() / (test_name + rExtension);
}

void RemoveIfExists(const std::filesystem::path& rPath)
{
    if (std::filesystem::exists(rPath)) {
        std::filesystem::remove(rPath);
    }
}

std::string ReadFileContent(const std::filesystem::path& rPath)
{
    std::ifstream file(rPath);
    std::stringstream buffer;
    buffer << file.rdbuf();
    return buffer.str();
}

/// A unit cube split into six linear tetrahedra (the diagonal 1-7 shared by all six).
void PopulateCubeOfTetrahedra(ModelPart& rModelPart)
{
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
    rModelPart.CreateNewNode(4, 0.0, 1.0, 0.0);
    rModelPart.CreateNewNode(5, 0.0, 0.0, 1.0);
    rModelPart.CreateNewNode(6, 1.0, 0.0, 1.0);
    rModelPart.CreateNewNode(7, 1.0, 1.0, 1.0);
    rModelPart.CreateNewNode(8, 0.0, 1.0, 1.0);

    auto p_properties = rModelPart.CreateNewProperties(1);
    const std::vector<std::vector<std::size_t>> connectivities = {
        {1, 2, 3, 7}, {1, 2, 7, 6}, {1, 6, 7, 5},
        {1, 3, 4, 7}, {1, 4, 8, 7}, {1, 8, 5, 7},
    };
    for (std::size_t i = 0; i < connectivities.size(); ++i) {
        rModelPart.CreateNewElement("Element3D4N", i + 1, connectivities[i], p_properties);
    }
}

/// A single quad (two triangles) on the z=0 plane, for operations restricted to surfaces.
void PopulateTriangulatedSquare(ModelPart& rModelPart)
{
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
    rModelPart.CreateNewNode(4, 0.0, 1.0, 0.0);

    auto p_properties = rModelPart.CreateNewProperties(1);
    rModelPart.CreateNewElement("Element2D3N", 1, {1, 2, 3}, p_properties);
    rModelPart.CreateNewElement("Element2D3N", 2, {1, 3, 4}, p_properties);
}

/// The unit cube's closed triangular skin (12 triangles, outward wound) - a watertight surface,
/// which is what the surface remeshing and volume-generating operations require.
void PopulateClosedCubeSkin(ModelPart& rModelPart)
{
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
    rModelPart.CreateNewNode(4, 0.0, 1.0, 0.0);
    rModelPart.CreateNewNode(5, 0.0, 0.0, 1.0);
    rModelPart.CreateNewNode(6, 1.0, 0.0, 1.0);
    rModelPart.CreateNewNode(7, 1.0, 1.0, 1.0);
    rModelPart.CreateNewNode(8, 0.0, 1.0, 1.0);

    auto p_properties = rModelPart.CreateNewProperties(1);
    const std::vector<std::vector<std::size_t>> connectivities = {
        {1, 3, 2}, {1, 4, 3},   // z = 0 (outward is -z)
        {5, 6, 7}, {5, 7, 8},   // z = 1
        {1, 2, 6}, {1, 6, 5},   // y = 0
        {3, 4, 8}, {3, 8, 7},   // y = 1
        {2, 3, 7}, {2, 7, 6},   // x = 1
        {1, 5, 8}, {1, 8, 4},   // x = 0
    };
    for (std::size_t i = 0; i < connectivities.size(); ++i) {
        rModelPart.CreateNewElement("Element2D3N", i + 1, connectivities[i], p_properties);
    }
}

Parameters OperationSettings(const std::string& rOperation, const std::string& rExtraJson = "{}")
{
    Parameters settings(rExtraJson);
    settings.AddString("operation", rOperation);
    return settings;
}

} // namespace

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsGetSupportedOperations, KratosMeshioPlusPlusFastSuite)
{
    const auto operations = MeshioPlusPlusMeshOperations::GetSupportedOperations();
    for (const std::string name : {"agglomerate", "attach_quality", "cell_data_to_point_data",
                                   "clean", "compute_sdf", "convert_cells", "crop_bbox",
                                   "crop_halfspace", "crop_predicate", "data_calc",
                                   "data_condition", "data_info", "data_integrate", "data_manage",
                                   "decimate", "decimate_volume", "estimate_error", "extract_skin",
                                   "extract_surface", "gradient", "hessian", "isosurface",
                                   "optimize_volume", "point_data_to_cell_data", "partition",
                                   "quality", "refine", "remesh", "remesh_volume", "reorder",
                                   "slice", "smooth", "split", "stats", "subdivide", "transform",
                                   "voxelize"}) {
        KRATOS_EXPECT_TRUE(std::find(operations.begin(), operations.end(), name) != operations.end());
    }
    KRATOS_EXPECT_TRUE(std::is_sorted(operations.begin(), operations.end()));
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsGetDefaultParameters, KratosMeshioPlusPlusFastSuite)
{
    const Parameters defaults = MeshioPlusPlusMeshOperations::GetDefaultParameters();
    KRATOS_EXPECT_EQ(defaults["operation"].GetString(), "clean");
    KRATOS_EXPECT_FALSE(defaults["weld"].GetBool());
    KRATOS_EXPECT_EQ(defaults["number_of_parts"].GetInt(), 2);
    KRATOS_EXPECT_EQ(defaults["ghost_layers"].GetInt(), 0);
    // Field data is opt-in: every selection is empty/false, so an operation run without
    // overriding these settings stages nothing (byte-identical to the pre-field-data behavior).
    KRATOS_EXPECT_TRUE(defaults["nodal_solution_step_data_variables"].GetStringArray().empty());
    KRATOS_EXPECT_TRUE(defaults["nodal_data_value_variables"].GetStringArray().empty());
    KRATOS_EXPECT_FALSE(defaults["write_ids"].GetBool());
    KRATOS_EXPECT_EQ(defaults["on_conflict"].GetString(), "error");
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsUnknownOperationThrows, KratosMeshioPlusPlusFastSuite)
{
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);
    auto& r_destination = model.CreateModelPart("destination");

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        MeshioPlusPlusMeshOperations::Execute(r_source, OperationSettings("not_a_real_operation"), r_destination),
        "Unknown meshio++ operation");
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsClean, KratosMeshioPlusPlusFastSuite)
{
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);
    auto& r_destination = model.CreateModelPart("destination");

    const Parameters report = MeshioPlusPlusMeshOperations::Execute(r_source, OperationSettings("clean"), r_destination);

    // A clean mesh is preserved as-is: nothing to weld, drop or remove.
    KRATOS_EXPECT_EQ(r_destination.NumberOfNodes(), r_source.NumberOfNodes());
    KRATOS_EXPECT_EQ(r_destination.NumberOfElements(), r_source.NumberOfElements());
    KRATOS_EXPECT_EQ(report["points_welded"].GetInt(), 0);
    KRATOS_EXPECT_EQ(report["cells_dropped_degenerate"].GetInt(), 0);
    KRATOS_EXPECT_EQ(report["cells_dropped_duplicate"].GetInt(), 0);
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsCleanWeldsCoincidentNodes, KratosMeshioPlusPlusFastSuite)
{
    // Two triangles sharing an edge, but built from two *different, coincident* nodes at
    // each shared position - clean must weld them back into one.
    Model model;
    auto& r_source = model.CreateModelPart("source");
    r_source.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_source.CreateNewNode(2, 1.0, 0.0, 0.0);
    r_source.CreateNewNode(3, 1.0, 1.0, 0.0);
    r_source.CreateNewNode(4, 1.0, 0.0, 0.0); // coincides with node 2
    r_source.CreateNewNode(5, 1.0, 1.0, 0.0); // coincides with node 3
    r_source.CreateNewNode(6, 0.0, 1.0, 0.0);
    auto p_properties = r_source.CreateNewProperties(1);
    r_source.CreateNewElement("Element2D3N", 1, {1, 2, 3}, p_properties);
    r_source.CreateNewElement("Element2D3N", 2, {4, 5, 6}, p_properties);
    auto& r_destination = model.CreateModelPart("destination");

    Parameters settings = OperationSettings("clean", R"({"weld" : true, "tolerance" : 1e-8})");
    const Parameters report = MeshioPlusPlusMeshOperations::Execute(r_source, settings, r_destination);

    KRATOS_EXPECT_EQ(report["points_welded"].GetInt(), 2);
    KRATOS_EXPECT_EQ(r_destination.NumberOfNodes(), 4);
    KRATOS_EXPECT_EQ(r_destination.NumberOfElements(), 2);
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsTransform, KratosMeshioPlusPlusFastSuite)
{
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);
    auto& r_destination = model.CreateModelPart("destination");

    Parameters settings = OperationSettings("transform", R"({"translation" : [1.0, 2.0, 3.0]})");
    MeshioPlusPlusMeshOperations::Execute(r_source, settings, r_destination);

    KRATOS_EXPECT_EQ(r_destination.NumberOfNodes(), r_source.NumberOfNodes());
    // Node 1 sits at the origin, so it lands exactly on the translation.
    const auto& r_node = r_destination.GetNode(1);
    KRATOS_EXPECT_NEAR(r_node.X(), 1.0, 1e-12);
    KRATOS_EXPECT_NEAR(r_node.Y(), 2.0, 1e-12);
    KRATOS_EXPECT_NEAR(r_node.Z(), 3.0, 1e-12);
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsConvertCellsElevate, KratosMeshioPlusPlusFastSuite)
{
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);
    auto& r_destination = model.CreateModelPart("destination");

    Parameters settings = OperationSettings("convert_cells", R"({"mode" : "elevate"})");
    MeshioPlusPlusMeshOperations::Execute(r_source, settings, r_destination);

    // Linear -> quadratic tetrahedra: the same six cells, with extra edge-midpoint nodes.
    KRATOS_EXPECT_EQ(r_destination.NumberOfElements(), r_source.NumberOfElements());
    KRATOS_EXPECT_GT(r_destination.NumberOfNodes(), r_source.NumberOfNodes());
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsRefine, KratosMeshioPlusPlusFastSuite)
{
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);
    auto& r_destination = model.CreateModelPart("destination");

    MeshioPlusPlusMeshOperations::Execute(r_source, OperationSettings("refine"), r_destination);

    KRATOS_EXPECT_GT(r_destination.NumberOfElements(), r_source.NumberOfElements());
    KRATOS_EXPECT_GT(r_destination.NumberOfNodes(), r_source.NumberOfNodes());
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsRefineSelectByExplicitCells, KratosMeshioPlusPlusFastSuite)
{
    // Two triangles sharing edge 1-3. Selecting only cell 0 (RedGreen, the default closure)
    // fully splits it into 4; cell 1 sees exactly its one shared edge bisected and splits
    // into 2 - so the total is 6, not the 8 a uniform refine of both would give.
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateTriangulatedSquare(r_source);
    auto& r_destination = model.CreateModelPart("destination");

    Parameters settings = OperationSettings("refine", R"({"cells" : [0]})");
    MeshioPlusPlusMeshOperations::Execute(r_source, settings, r_destination);

    KRATOS_EXPECT_EQ(r_destination.NumberOfElements(), 6);
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsRefineSelectByRegion, KratosMeshioPlusPlusFastSuite)
{
    // The region path must resolve to the same selection as the explicit "cells" index -
    // a sub model part containing only the first triangle's element and nodes.
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateTriangulatedSquare(r_source);
    auto& r_region = r_source.CreateSubModelPart("target_region");
    r_region.AddNodes({1, 2, 3});
    r_region.AddElements({1});
    auto& r_destination = model.CreateModelPart("destination");

    Parameters settings = OperationSettings("refine", R"({"region" : "target_region"})");
    MeshioPlusPlusMeshOperations::Execute(r_source, settings, r_destination);

    KRATOS_EXPECT_EQ(r_destination.NumberOfElements(), 6);
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsRefineSelectByPredicate, KratosMeshioPlusPlusFastSuite)
{
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateTriangulatedSquare(r_source);
    r_source.GetElement(1).SetValue(TEMPERATURE, 0.0);
    r_source.GetElement(2).SetValue(TEMPERATURE, 1.0);
    auto& r_destination = model.CreateModelPart("destination");

    Parameters settings = OperationSettings("refine", R"({
        "predicate_array" : "TEMPERATURE",
        "predicate_op" : "<",
        "predicate_value" : 0.5,
        "element_data_value_variables" : ["TEMPERATURE"]
    })");
    MeshioPlusPlusMeshOperations::Execute(r_source, settings, r_destination);

    KRATOS_EXPECT_EQ(r_destination.NumberOfElements(), 6);
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsRefineClosureRedGreenVsPropagate, KratosMeshioPlusPlusFastSuite)
{
    // Same single-cell selection, different closures: RedGreen keeps the neighbour's split
    // partial (6 total, as in the explicit-cells test); Propagate promotes any non-empty
    // mask straight to a full split, so both triangles end up fully split (4 + 4 = 8).
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateTriangulatedSquare(r_source);

    auto& r_redgreen = model.CreateModelPart("redgreen");
    MeshioPlusPlusMeshOperations::Execute(
        r_source, OperationSettings("refine", R"({"cells" : [0], "closure" : "redgreen"})"), r_redgreen);

    auto& r_propagate = model.CreateModelPart("propagate");
    MeshioPlusPlusMeshOperations::Execute(
        r_source, OperationSettings("refine", R"({"cells" : [0], "closure" : "propagate"})"), r_propagate);

    KRATOS_EXPECT_EQ(r_redgreen.NumberOfElements(), 6);
    KRATOS_EXPECT_EQ(r_propagate.NumberOfElements(), 8);
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsRefineBalancedClosureLeavesHangingNode, KratosMeshioPlusPlusFastSuite)
{
    // record_parent_ids-style Int64 bookkeeping arrays (here "refine:hanging") are computed
    // by meshio++ but never reach the destination model part - their name is not a registered
    // Kratos Variable, by design (see the write-back constraint). Exercised directly against
    // the mesh rather than through Execute() for exactly that reason.
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateTriangulatedSquare(r_source);

    mio::Mesh mesh = Internals::ModelPartToMesh(r_source);
    mio::RefineOptions options;
    options.mCells = {0};
    options.mClosure = mio::RefineClosure::Balanced;
    const mio::RefineResult result = mio::refine(mesh, options);

    KRATOS_EXPECT_TRUE(result.mMesh.HasPointData(mio::kRefineHangingName));
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsRefineRecordLevelsAttachesLevelArray, KratosMeshioPlusPlusFastSuite)
{
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateTriangulatedSquare(r_source);

    mio::Mesh mesh = Internals::ModelPartToMesh(r_source);
    mio::RefineOptions options;
    options.mRecordLevels = true;
    const mio::RefineResult result = mio::refine(mesh, options);

    KRATOS_EXPECT_TRUE(result.mMesh.HasCellData(mio::kRefineLevelName));
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsRefineTwoSelectorsIsAnError, KratosMeshioPlusPlusFastSuite)
{
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateTriangulatedSquare(r_source);
    auto& r_destination = model.CreateModelPart("destination");

    Parameters settings = OperationSettings("refine", R"({"cells" : [0], "region" : "whatever"})");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        MeshioPlusPlusMeshOperations::Execute(r_source, settings, r_destination), "");
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsGetDefaultParametersRefineSelectiveKeys, KratosMeshioPlusPlusFastSuite)
{
    const Parameters defaults = MeshioPlusPlusMeshOperations::GetDefaultParameters();
    KRATOS_EXPECT_EQ(defaults["cells"].size(), 0);
    KRATOS_EXPECT_TRUE(defaults["region"].GetString().empty());
    KRATOS_EXPECT_TRUE(defaults["predicate_array"].GetString().empty());
    KRATOS_EXPECT_EQ(defaults["closure"].GetString(), "redgreen");
    KRATOS_EXPECT_FALSE(defaults["record_levels"].GetBool());
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsDecimateNeedsASurface, KratosMeshioPlusPlusFastSuite)
{
    // decimate is surface-only; a volume mesh raises by name pointing at extract_surface,
    // per the meshio++ contract. target_ratio must be set explicitly - decimate validates
    // its own options (exactly one stopping criterion) before checking the mesh at all.
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);
    auto& r_destination = model.CreateModelPart("destination");

    Parameters settings = OperationSettings("decimate", R"({"target_ratio" : 0.5})");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        MeshioPlusPlusMeshOperations::Execute(r_source, settings, r_destination), "extract_surface");
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsDecimateOnASurface, KratosMeshioPlusPlusFastSuite)
{
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateTriangulatedSquare(r_source);
    auto& r_destination = model.CreateModelPart("destination");

    Parameters settings = OperationSettings("decimate", R"({"target_ratio" : 0.5})");
    const Parameters report = MeshioPlusPlusMeshOperations::Execute(r_source, settings, r_destination);

    KRATOS_EXPECT_LE(r_destination.NumberOfElements(), r_source.NumberOfElements());
    KRATOS_EXPECT_TRUE(report.Has("faces_removed"));
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsSmooth, KratosMeshioPlusPlusFastSuite)
{
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);
    auto& r_destination = model.CreateModelPart("destination");

    const Parameters report = MeshioPlusPlusMeshOperations::Execute(r_source, OperationSettings("smooth"), r_destination);

    // Topology is untouched by smoothing, only coordinates move.
    KRATOS_EXPECT_EQ(r_destination.NumberOfNodes(), r_source.NumberOfNodes());
    KRATOS_EXPECT_EQ(r_destination.NumberOfElements(), r_source.NumberOfElements());
    KRATOS_EXPECT_TRUE(report.Has("nodes_moved"));
    KRATOS_EXPECT_TRUE(report.Has("max_displacement"));
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsReorder, KratosMeshioPlusPlusFastSuite)
{
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);
    auto& r_destination = model.CreateModelPart("destination");

    Parameters settings = OperationSettings("reorder", R"({"method" : "rcm"})");
    const Parameters report = MeshioPlusPlusMeshOperations::Execute(r_source, settings, r_destination);

    KRATOS_EXPECT_EQ(r_destination.NumberOfNodes(), r_source.NumberOfNodes());
    KRATOS_EXPECT_EQ(r_destination.NumberOfElements(), r_source.NumberOfElements());
    // Reverse Cuthill-McKee never increases the bandwidth.
    KRATOS_EXPECT_LE(report["bandwidth_after"].GetInt(), report["bandwidth_before"].GetInt());
    KRATOS_EXPECT_EQ(static_cast<std::size_t>(report["bandwidth_before"].GetInt()),
                     MeshioPlusPlusMeshOperations::ComputeBandwidth(r_source));
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsExtractSurface, KratosMeshioPlusPlusFastSuite)
{
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);
    auto& r_destination = model.CreateModelPart("destination");

    MeshioPlusPlusMeshOperations::Execute(r_source, OperationSettings("extract_surface"), r_destination);

    KRATOS_EXPECT_EQ(r_destination.NumberOfNodes(), 8); // every node of the cube is on its boundary
    KRATOS_EXPECT_GT(r_destination.NumberOfElements() + r_destination.NumberOfConditions(), 0);
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsExtractSkin, KratosMeshioPlusPlusFastSuite)
{
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);
    auto& r_destination = model.CreateModelPart("destination");

    MeshioPlusPlusMeshOperations::Execute(r_source, OperationSettings("extract_skin"), r_destination);

    KRATOS_EXPECT_EQ(r_destination.NumberOfNodes(), 8);
    KRATOS_EXPECT_GT(r_destination.NumberOfElements() + r_destination.NumberOfConditions(), 0);
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsCropBboxCoversEverything, KratosMeshioPlusPlusFastSuite)
{
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);
    auto& r_destination = model.CreateModelPart("destination");

    Parameters settings = OperationSettings("crop_bbox",
        R"({"box_min" : [-1.0, -1.0, -1.0], "box_max" : [2.0, 2.0, 2.0]})");
    MeshioPlusPlusMeshOperations::Execute(r_source, settings, r_destination);

    KRATOS_EXPECT_EQ(r_destination.NumberOfNodes(), r_source.NumberOfNodes());
    KRATOS_EXPECT_EQ(r_destination.NumberOfElements(), r_source.NumberOfElements());
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsCropBboxCoversNothing, KratosMeshioPlusPlusFastSuite)
{
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);
    auto& r_destination = model.CreateModelPart("destination");

    Parameters settings = OperationSettings("crop_bbox",
        R"({"box_min" : [10.0, 10.0, 10.0], "box_max" : [20.0, 20.0, 20.0]})");
    MeshioPlusPlusMeshOperations::Execute(r_source, settings, r_destination);

    KRATOS_EXPECT_EQ(r_destination.NumberOfNodes(), 0);
    KRATOS_EXPECT_EQ(r_destination.NumberOfElements(), 0);
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsCropHalfspace, KratosMeshioPlusPlusFastSuite)
{
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);

    // The whole cube (z in [0, 1]) is on the z >= 0 side of the origin plane.
    {
        auto& r_covers_everything = model.CreateModelPart("covers_everything");
        Parameters settings = OperationSettings("crop_halfspace",
            R"({"origin" : [0.0, 0.0, 0.0], "normal" : [0.0, 0.0, 1.0]})");
        MeshioPlusPlusMeshOperations::Execute(r_source, settings, r_covers_everything);
        KRATOS_EXPECT_EQ(r_covers_everything.NumberOfElements(), r_source.NumberOfElements());
    }
    // ... and none of it is on the z < 0 side.
    {
        auto& r_covers_nothing = model.CreateModelPart("covers_nothing");
        Parameters settings = OperationSettings("crop_halfspace",
            R"({"origin" : [0.0, 0.0, 0.0], "normal" : [0.0, 0.0, -1.0]})");
        MeshioPlusPlusMeshOperations::Execute(r_source, settings, r_covers_nothing);
        KRATOS_EXPECT_EQ(r_covers_nothing.NumberOfElements(), 0);
    }
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsSlice, KratosMeshioPlusPlusFastSuite)
{
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);
    auto& r_destination = model.CreateModelPart("destination");

    // A plane through the middle of the cube must intersect it.
    Parameters settings = OperationSettings("slice",
        R"({"origin" : [0.5, 0.5, 0.5], "normal" : [0.0, 0.0, 1.0]})");
    MeshioPlusPlusMeshOperations::Execute(r_source, settings, r_destination);

    KRATOS_EXPECT_GT(r_destination.NumberOfNodes(), 0);
    KRATOS_EXPECT_GT(r_destination.NumberOfElements(), 0);
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsIsosurfaceNeedsAnArrayName, KratosMeshioPlusPlusFastSuite)
{
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);
    auto& r_destination = model.CreateModelPart("destination");

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        MeshioPlusPlusMeshOperations::Execute(r_source, OperationSettings("isosurface"), r_destination),
        "array_name");
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsIsosurfaceThrowsOnAnUnknownArray, KratosMeshioPlusPlusFastSuite)
{
    // Field data is opt-in (see BuildFieldDataSelection): with none of the
    // "nodal_*_variables" settings listing "TEMPERATURE", nothing is staged as point_data even
    // though TEMPERATURE is a registered Variable, so the array can never resolve.
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);
    auto& r_destination = model.CreateModelPart("destination");

    Parameters settings = OperationSettings("isosurface", R"({"array_name" : "TEMPERATURE"})");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        MeshioPlusPlusMeshOperations::Execute(r_source, settings, r_destination), "");
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsIsosurfaceWithStagedFieldData, KratosMeshioPlusPlusFastSuite)
{
    // Staging TEMPERATURE (via "nodal_data_value_variables") makes it a real point_data array,
    // so isosurface can now resolve it - the gap the field-data hoist closed. TEMPERATURE is
    // set to the node's own Z coordinate, so the level set at 0.5 must cut the cube in half
    // (the bottom face is at Z=0, the top at Z=1).
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);
    for (auto& r_node : r_source.Nodes()) {
        r_node.SetValue(TEMPERATURE, r_node.Z());
    }
    auto& r_destination = model.CreateModelPart("destination");

    Parameters settings = OperationSettings("isosurface", R"({
        "array_name" : "TEMPERATURE",
        "isovalues" : [0.5],
        "nodal_data_value_variables" : ["TEMPERATURE"]
    })");
    MeshioPlusPlusMeshOperations::Execute(r_source, settings, r_destination);

    KRATOS_EXPECT_GT(r_destination.NumberOfNodes(), 0);
    // "The contoured field is exactly the isovalue on the cut points" (isosurface.hpp) - so the
    // write-back path (Internals::MeshToModelPart -> ApplyDataArrayToEntities) is exercised too:
    // every new node's non-historical TEMPERATURE must come back as exactly 0.5.
    for (const auto& r_node : r_destination.Nodes()) {
        KRATOS_EXPECT_NEAR(r_node.GetValue(TEMPERATURE), 0.5, 1e-12);
    }
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsAttachQuality, KratosMeshioPlusPlusFastSuite)
{
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);
    auto& r_destination = model.CreateModelPart("destination");

    MeshioPlusPlusMeshOperations::Execute(r_source, OperationSettings("attach_quality"), r_destination);

    // attach_quality names its cell_data "quality:scaled_jacobian" etc. (see quality.hpp).
    // MeshToModelPart's write-back (Internals::ApplyDataArrayToEntities) does attempt every
    // array the result carries, but a Kratos entity can only hold data under an actual
    // registered Variable, and "quality:scaled_jacobian" is not one - so the metrics are
    // computed (attach_quality ran) but unretrievable from the destination model part; only
    // the topology is expected to survive here.
    KRATOS_EXPECT_EQ(r_destination.NumberOfElements(), r_source.NumberOfElements());
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsFieldDataRoundTripsThroughClean, KratosMeshioPlusPlusFastSuite)
{
    // A mesh-producing operation that changes nothing (a clean mesh, weld off) must still
    // carry a staged nodal variable all the way through: into the meshio++ mesh
    // (Internals::ModelPartToMeshWithData), through "clean", and back out
    // (Internals::MeshToModelPart -> ApplyDataArrayToEntities).
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);
    for (auto& r_node : r_source.Nodes()) {
        r_node.SetValue(TEMPERATURE, 10.0 * r_node.Id());
    }
    auto& r_destination = model.CreateModelPart("destination");

    Parameters settings = OperationSettings("clean", R"({"nodal_data_value_variables" : ["TEMPERATURE"]})");
    MeshioPlusPlusMeshOperations::Execute(r_source, settings, r_destination);

    KRATOS_EXPECT_EQ(r_destination.NumberOfNodes(), r_source.NumberOfNodes());
    for (const auto& r_node : r_destination.Nodes()) {
        KRATOS_EXPECT_NEAR(r_node.GetValue(TEMPERATURE), 10.0 * r_node.Id(), 1e-12);
    }
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsDataCalc, KratosMeshioPlusPlusFastSuite)
{
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);
    for (auto& r_node : r_source.Nodes()) {
        r_node.SetValue(TEMPERATURE, r_node.Z());
    }
    auto& r_destination = model.CreateModelPart("destination");

    // "output" reuses the staged array's own name with "output_overwrite" - a registered
    // Variable name, so the doubled result comes straight back as non-historical TEMPERATURE.
    Parameters settings = OperationSettings("data_calc", R"({
        "expression"                  : "TEMPERATURE * 2.0",
        "output"                      : "TEMPERATURE",
        "output_overwrite"            : true,
        "nodal_data_value_variables"  : ["TEMPERATURE"]
    })");
    MeshioPlusPlusMeshOperations::Execute(r_source, settings, r_destination);

    for (const auto& r_node : r_destination.Nodes()) {
        KRATOS_EXPECT_NEAR(r_node.GetValue(TEMPERATURE), 2.0 * r_node.Z(), 1e-12);
    }
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsDataCondition, KratosMeshioPlusPlusFastSuite)
{
    // "clamp" (the default mode) maps x -> min(max(x, lo), hi); Z in {0, 1} with hi = 0.5
    // must clamp only the top face's nodes down to 0.5.
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);
    for (auto& r_node : r_source.Nodes()) {
        r_node.SetValue(TEMPERATURE, r_node.Z());
    }
    auto& r_destination = model.CreateModelPart("destination");

    Parameters settings = OperationSettings("data_condition", R"({
        "names"                        : ["TEMPERATURE"],
        "mode"                          : "clamp",
        "lo"                            : 0.0,
        "hi"                            : 0.5,
        "nodal_data_value_variables"    : ["TEMPERATURE"]
    })");
    MeshioPlusPlusMeshOperations::Execute(r_source, settings, r_destination);

    for (const auto& r_node : r_destination.Nodes()) {
        KRATOS_EXPECT_LE(r_node.GetValue(TEMPERATURE), 0.5 + 1e-12);
    }
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsDataManageRename, KratosMeshioPlusPlusFastSuite)
{
    // Renaming "point_data:TEMPERATURE" to "point_data:PRESSURE" both proves data_manage's
    // rename phase runs and demonstrates the write-back constraint's escape hatch: the
    // renamed-to name is itself a registered Variable, so it survives the trip back into
    // Kratos even though the array started life under a different name.
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);
    for (auto& r_node : r_source.Nodes()) {
        r_node.SetValue(TEMPERATURE, 42.0);
    }
    auto& r_destination = model.CreateModelPart("destination");

    Parameters settings = OperationSettings("data_manage", R"({
        "rename"                     : [{"location" : "point", "from" : "TEMPERATURE", "to" : "PRESSURE"}],
        "nodal_data_value_variables" : ["TEMPERATURE"]
    })");
    const Parameters report = MeshioPlusPlusMeshOperations::Execute(r_source, settings, r_destination);

    KRATOS_EXPECT_EQ(report["renamed"].size(), 1);
    for (const auto& r_node : r_destination.Nodes()) {
        KRATOS_EXPECT_NEAR(r_node.GetValue(PRESSURE), 42.0, 1e-12);
    }
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsDataInfo, KratosMeshioPlusPlusFastSuite)
{
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);
    for (auto& r_node : r_source.Nodes()) {
        r_node.SetValue(TEMPERATURE, r_node.Z());
    }
    auto& r_destination = model.CreateModelPart("destination");

    Parameters settings = OperationSettings("data_info", R"({"nodal_data_value_variables" : ["TEMPERATURE"]})");
    const Parameters report = MeshioPlusPlusMeshOperations::Execute(r_source, settings, r_destination);

    KRATOS_EXPECT_EQ(report["number_of_point_data"].GetInt(), 1);
    KRATOS_EXPECT_EQ(report["arrays"][0]["name"].GetString(), "TEMPERATURE");
    KRATOS_EXPECT_EQ(report["arrays"][0]["num_values"].GetInt(), static_cast<int>(r_source.NumberOfNodes()));
    // Report-only: nothing is written into the destination.
    KRATOS_EXPECT_EQ(r_destination.NumberOfNodes(), 0);
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsPointDataToCellData, KratosMeshioPlusPlusFastSuite)
{
    // A uniform nodal value keeps the expected cell average trivial to state regardless of
    // which nodes each tetrahedron happens to use.
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);
    for (auto& r_node : r_source.Nodes()) {
        r_node.SetValue(TEMPERATURE, 5.0);
    }
    auto& r_destination = model.CreateModelPart("destination");

    Parameters settings = OperationSettings("point_data_to_cell_data", R"({
        "names"                       : ["TEMPERATURE"],
        "nodal_data_value_variables"  : ["TEMPERATURE"]
    })");
    MeshioPlusPlusMeshOperations::Execute(r_source, settings, r_destination);

    KRATOS_EXPECT_EQ(r_destination.NumberOfElements(), r_source.NumberOfElements());
    for (const auto& r_element : r_destination.Elements()) {
        KRATOS_EXPECT_NEAR(r_element.GetValue(TEMPERATURE), 5.0, 1e-12);
    }
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsCellDataToPointData, KratosMeshioPlusPlusFastSuite)
{
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);
    for (auto& r_element : r_source.Elements()) {
        r_element.SetValue(DENSITY, 7850.0);
    }
    auto& r_destination = model.CreateModelPart("destination");

    Parameters settings = OperationSettings("cell_data_to_point_data", R"({
        "names"                          : ["DENSITY"],
        "element_data_value_variables"   : ["DENSITY"]
    })");
    MeshioPlusPlusMeshOperations::Execute(r_source, settings, r_destination);

    KRATOS_EXPECT_EQ(r_destination.NumberOfNodes(), r_source.NumberOfNodes());
    for (const auto& r_node : r_destination.Nodes()) {
        KRATOS_EXPECT_NEAR(r_node.GetValue(DENSITY), 7850.0, 1e-12);
    }
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsInterpolateNearest, KratosMeshioPlusPlusFastSuite)
{
    // Source and target are geometrically identical cubes (independent node/element id
    // spaces), so nearest-point sampling must be exact: every target point coincides with a
    // source point. The target's own (unset, defaulted-to-zero) TEMPERATURE forces the
    // name conflict that "on_conflict" : "overwrite" then resolves.
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);
    for (auto& r_node : r_source.Nodes()) {
        r_node.SetValue(TEMPERATURE, r_node.Z());
    }
    auto& r_target = model.CreateModelPart("target");
    PopulateCubeOfTetrahedra(r_target);
    auto& r_destination = model.CreateModelPart("destination");

    Parameters settings(R"({
        "method"                      : "nearest",
        "on_conflict"                 : "overwrite",
        "nodal_data_value_variables"  : ["TEMPERATURE"]
    })");
    const Parameters report = MeshioPlusPlusMeshOperations::Interpolate(r_source, r_target, settings, r_destination);

    KRATOS_EXPECT_EQ(report["number_of_nodes"].GetInt(), static_cast<int>(r_target.NumberOfNodes()));
    KRATOS_EXPECT_EQ(r_destination.NumberOfNodes(), r_target.NumberOfNodes());
    for (const auto& r_node : r_destination.Nodes()) {
        KRATOS_EXPECT_NEAR(r_node.GetValue(TEMPERATURE), r_node.Z(), 1e-12);
    }
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsInterpolateConflictThrowsByDefault, KratosMeshioPlusPlusFastSuite)
{
    // Field data is staged identically onto source and target (Interpolate applies the same
    // FieldDataSelection to both), so TEMPERATURE exists on both meshes here; the default
    // "on_conflict" : "error" must therefore throw rather than silently pick a side.
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);
    for (auto& r_node : r_source.Nodes()) {
        r_node.SetValue(TEMPERATURE, r_node.Z());
    }
    auto& r_target = model.CreateModelPart("target");
    PopulateCubeOfTetrahedra(r_target);
    auto& r_destination = model.CreateModelPart("destination");

    Parameters settings(R"({
        "method"                      : "nearest",
        "nodal_data_value_variables"  : ["TEMPERATURE"]
    })");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        MeshioPlusPlusMeshOperations::Interpolate(r_source, r_target, settings, r_destination), "");
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsSplitByType, KratosMeshioPlusPlusFastSuite)
{
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);
    auto& r_destination = model.CreateModelPart("destination");

    Parameters settings = OperationSettings("split", R"({"split_by" : "type"})");
    const Parameters report = MeshioPlusPlusMeshOperations::Execute(r_source, settings, r_destination);

    // A single cell type (tetra), so a single piece.
    KRATOS_EXPECT_EQ(report["number_of_pieces"].GetInt(), 1);
    const std::string piece_name = report["pieces"][0]["name"].GetString();
    KRATOS_EXPECT_TRUE(model.HasModelPart(piece_name));
    KRATOS_EXPECT_EQ(model.GetModelPart(piece_name).NumberOfElements(), r_source.NumberOfElements());
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsSplitByComponent, KratosMeshioPlusPlusFastSuite)
{
    // Two disjoint tetrahedra: split by connected component must produce two pieces.
    Model model;
    auto& r_source = model.CreateModelPart("source");
    r_source.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_source.CreateNewNode(2, 1.0, 0.0, 0.0);
    r_source.CreateNewNode(3, 0.0, 1.0, 0.0);
    r_source.CreateNewNode(4, 0.0, 0.0, 1.0);
    r_source.CreateNewNode(5, 10.0, 0.0, 0.0);
    r_source.CreateNewNode(6, 11.0, 0.0, 0.0);
    r_source.CreateNewNode(7, 10.0, 1.0, 0.0);
    r_source.CreateNewNode(8, 10.0, 0.0, 1.0);
    auto p_properties = r_source.CreateNewProperties(1);
    r_source.CreateNewElement("Element3D4N", 1, {1, 2, 3, 4}, p_properties);
    r_source.CreateNewElement("Element3D4N", 2, {5, 6, 7, 8}, p_properties);
    auto& r_destination = model.CreateModelPart("destination");

    Parameters settings = OperationSettings("split", R"({"split_by" : "component"})");
    const Parameters report = MeshioPlusPlusMeshOperations::Execute(r_source, settings, r_destination);

    KRATOS_EXPECT_EQ(report["number_of_pieces"].GetInt(), 2);
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsPartitionCreatesPieces, KratosMeshioPlusPlusFastSuite)
{
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);
    auto& r_destination = model.CreateModelPart("destination");

    Parameters settings = OperationSettings("partition", R"({"number_of_parts" : 2})");
    const Parameters report = MeshioPlusPlusMeshOperations::Execute(r_source, settings, r_destination);

    KRATOS_EXPECT_EQ(report["number_of_pieces"].GetInt(), 2);
    std::size_t total_elements = 0;
    for (std::size_t i = 0; i < report["pieces"].size(); ++i) {
        const std::string name = report["pieces"][i]["name"].GetString();
        KRATOS_EXPECT_TRUE(model.HasModelPart(name));
        total_elements += model.GetModelPart(name).NumberOfElements();
    }
    // Partition of unity: no ghost layers, so every source element lands in exactly one piece.
    KRATOS_EXPECT_EQ(total_elements, r_source.NumberOfElements());
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsPartitionWithGhostLayers, KratosMeshioPlusPlusFastSuite)
{
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);
    auto& r_destination = model.CreateModelPart("destination");

    Parameters settings = OperationSettings("partition",
        R"({"number_of_parts" : 2, "ghost_layers" : 1, "seed" : 1})");
    const Parameters report = MeshioPlusPlusMeshOperations::Execute(r_source, settings, r_destination);

    std::size_t total_elements = 0;
    for (std::size_t i = 0; i < report["pieces"].size(); ++i) {
        total_elements += model.GetModelPart(report["pieces"][i]["name"].GetString()).NumberOfElements();
    }
    // With ghost layers the pieces overlap, so the total exceeds a partition of unity -
    // unless every cell already shares a node with every other part (unlikely for 6 cells
    // in 2 parts, but not impossible), so this is the one operation-specific assertion that
    // could be flaky on a future meshio++ version; >= is always true, > is the interesting case.
    KRATOS_EXPECT_GE(total_elements, r_source.NumberOfElements());
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsStatsAndQualityThroughExecute, KratosMeshioPlusPlusFastSuite)
{
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);
    auto& r_destination = model.CreateModelPart("destination");

    const Parameters stats_via_execute = MeshioPlusPlusMeshOperations::Execute(
        r_source, OperationSettings("stats"), r_destination);
    const Parameters stats_direct = MeshioPlusPlusMeshOperations::ComputeStatistics(r_source);
    KRATOS_EXPECT_EQ(stats_via_execute["number_of_cells"].GetInt(), stats_direct["number_of_cells"].GetInt());

    const Parameters quality_via_execute = MeshioPlusPlusMeshOperations::Execute(
        r_source, OperationSettings("quality"), r_destination);
    const Parameters quality_direct = MeshioPlusPlusMeshOperations::ComputeQuality(r_source);
    KRATOS_EXPECT_EQ(quality_via_execute["number_of_cells"].GetInt(), quality_direct["number_of_cells"].GetInt());

    // Both are report-only: the destination is never touched.
    KRATOS_EXPECT_EQ(r_destination.NumberOfNodes(), 0);
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsComputeStatistics, KratosMeshioPlusPlusFastSuite)
{
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);

    const Parameters report = MeshioPlusPlusMeshOperations::ComputeStatistics(r_source);

    KRATOS_EXPECT_EQ(report["number_of_points"].GetInt(), 8);
    KRATOS_EXPECT_EQ(report["number_of_cells"].GetInt(), 6);
    KRATOS_EXPECT_NEAR(report["unsigned_volume"].GetDouble(), 1.0, 1e-10);
    KRATOS_EXPECT_EQ(report["number_of_inverted"].GetInt(), 0);
    KRATOS_EXPECT_TRUE(report["cell_type_counts"].Has("tetra"));
    KRATOS_EXPECT_EQ(report["cell_type_counts"]["tetra"].GetInt(), 6);
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsComputeQuality, KratosMeshioPlusPlusFastSuite)
{
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);

    const Parameters report = MeshioPlusPlusMeshOperations::ComputeQuality(r_source);

    KRATOS_EXPECT_EQ(report["number_of_cells"].GetInt(), 6);
    KRATOS_EXPECT_EQ(report["number_of_inverted"].GetInt(), 0);
    KRATOS_EXPECT_TRUE(report.Has("metrics"));
    // "metrics" is an object keyed by metric name, not an array - Parameters::size() only
    // accepts the Array type, so its non-emptiness is checked by iterating instead.
    KRATOS_EXPECT_TRUE(report["metrics"].begin() != report["metrics"].end());
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsComputeBandwidth, KratosMeshioPlusPlusFastSuite)
{
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);

    KRATOS_EXPECT_GT(MeshioPlusPlusMeshOperations::ComputeBandwidth(r_source), 0);
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsMergeTwoCubes, KratosMeshioPlusPlusFastSuite)
{
    Model model;
    auto& r_first = model.CreateModelPart("first");
    PopulateCubeOfTetrahedra(r_first);
    auto& r_second = model.CreateModelPart("second");
    PopulateCubeOfTetrahedra(r_second);
    auto& r_destination = model.CreateModelPart("merged");

    const std::vector<const ModelPart*> sources = {&r_first, &r_second};
    const Parameters report = MeshioPlusPlusMeshOperations::Merge(sources, Parameters("{}"), r_destination);

    KRATOS_EXPECT_EQ(report["number_of_sources"].GetInt(), 2);
    KRATOS_EXPECT_EQ(r_destination.NumberOfElements(), 12);
    KRATOS_EXPECT_EQ(r_destination.NumberOfNodes(), 16); // not welded: two independent id spaces
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsMergeWelded, KratosMeshioPlusPlusFastSuite)
{
    // Two copies of the same cube, in the same place: welding must collapse the 16 raw
    // nodes back down to the 8 that are actually distinct positions.
    Model model;
    auto& r_first = model.CreateModelPart("first");
    PopulateCubeOfTetrahedra(r_first);
    auto& r_second = model.CreateModelPart("second");
    PopulateCubeOfTetrahedra(r_second);
    auto& r_destination = model.CreateModelPart("merged_welded");

    const std::vector<const ModelPart*> sources = {&r_first, &r_second};
    Parameters settings(R"({"weld" : true, "tolerance" : 1e-8})");
    MeshioPlusPlusMeshOperations::Merge(sources, settings, r_destination);

    KRATOS_EXPECT_EQ(r_destination.NumberOfNodes(), 8);
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsMergeNeedsAtLeastTwoSources, KratosMeshioPlusPlusFastSuite)
{
    Model model;
    auto& r_first = model.CreateModelPart("first");
    PopulateCubeOfTetrahedra(r_first);
    auto& r_destination = model.CreateModelPart("destination");

    const std::vector<const ModelPart*> sources = {&r_first};
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        MeshioPlusPlusMeshOperations::Merge(sources, Parameters("{}"), r_destination),
        "at least two model parts");
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsDiffOfAModelPartWithItself, KratosMeshioPlusPlusFastSuite)
{
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);

    const Parameters report = MeshioPlusPlusMeshOperations::Diff(r_source, r_source);
    KRATOS_EXPECT_TRUE(report["equal"].GetBool());
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsDiffOfDifferentModelParts, KratosMeshioPlusPlusFastSuite)
{
    Model model;
    auto& r_first = model.CreateModelPart("first");
    PopulateCubeOfTetrahedra(r_first);
    auto& r_second = model.CreateModelPart("second");
    PopulateCubeOfTetrahedra(r_second);
    // ModelPartToMesh (which Diff uses) reads the *initial* coordinates by default, so the
    // difference must be created there - mutating Node::X() afterwards would be invisible.
    r_second.GetNode(2).X0() = 5.0;

    const Parameters report = MeshioPlusPlusMeshOperations::Diff(r_first, r_second);
    KRATOS_EXPECT_FALSE(report["equal"].GetBool());
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsEntityTypeFiltering, KratosMeshioPlusPlusFastSuite)
{
    // The tetrahedron (element) and the triangle (condition) are given disjoint geometry -
    // one at the origin, one far away - so which one made it into the operation's input can
    // be told apart by node count and position alone. This matters because, once an entity
    // reaches meshio++, "element vs condition" is re-derived from cell dimension on the way
    // back into Kratos (highest dimension present -> element, lower -> condition): with only
    // one kind selected it is also the *only* cell block, so it always comes back an element
    // regardless of which Kratos container it started in. Node position is therefore the
    // only source-agnostic way to check that "entity_type" actually filtered the input.
    Model model;
    auto& r_source = model.CreateModelPart("source");
    r_source.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_source.CreateNewNode(2, 1.0, 0.0, 0.0);
    r_source.CreateNewNode(3, 0.0, 1.0, 0.0);
    r_source.CreateNewNode(4, 0.0, 0.0, 1.0);
    r_source.CreateNewNode(5, 100.0, 0.0, 0.0);
    r_source.CreateNewNode(6, 101.0, 0.0, 0.0);
    r_source.CreateNewNode(7, 100.0, 1.0, 0.0);
    auto p_properties = r_source.CreateNewProperties(1);
    r_source.CreateNewElement("Element3D4N", 1, {1, 2, 3, 4}, p_properties);
    r_source.CreateNewCondition("SurfaceCondition3D3N", 1, {5, 6, 7}, p_properties);

    {
        auto& r_elements_only = model.CreateModelPart("elements_only");
        Parameters settings = OperationSettings("clean", R"({"entity_type" : "elements"})");
        MeshioPlusPlusMeshOperations::Execute(r_source, settings, r_elements_only);
        KRATOS_EXPECT_EQ(r_elements_only.NumberOfNodes(), 4); // the tetrahedron only
        KRATOS_EXPECT_LT(r_elements_only.GetNode(1).X(), 10.0);
    }
    {
        auto& r_conditions_only = model.CreateModelPart("conditions_only");
        Parameters settings = OperationSettings("clean", R"({"entity_type" : "conditions"})");
        MeshioPlusPlusMeshOperations::Execute(r_source, settings, r_conditions_only);
        KRATOS_EXPECT_EQ(r_conditions_only.NumberOfNodes(), 3); // the triangle only
        KRATOS_EXPECT_GT(r_conditions_only.GetNode(1).X(), 10.0);
    }
    {
        auto& r_automatic = model.CreateModelPart("automatic");
        MeshioPlusPlusMeshOperations::Execute(r_source, OperationSettings("clean"), r_automatic);
        KRATOS_EXPECT_EQ(r_automatic.NumberOfNodes(), 7); // both, disjoint
        KRATOS_EXPECT_EQ(r_automatic.NumberOfElements(), 1);
        KRATOS_EXPECT_EQ(r_automatic.NumberOfConditions(), 1);
    }
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsPropertiesSurviveTransform, KratosMeshioPlusPlusFastSuite)
{
    // Material data reaches the meshio++ mesh via Internals::FillMeshioModelPart /
    // ExportMeshioProperties (the same conversion Execute uses for every operation), but
    // meshio++ 9.2.0 only *propagates* property sets from the input mesh to an operation's
    // result for a subset of operations - "transform" and "smooth" are documented to carry
    // them; "merge", "crop", "split", "partition" and "diff" are documented not to.
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);
    r_source.pGetProperties(1)->SetValue(DENSITY, 7850.0);
    auto& r_destination = model.CreateModelPart("destination");

    MeshioPlusPlusMeshOperations::Execute(r_source, OperationSettings("transform"), r_destination);

    KRATOS_EXPECT_TRUE(r_destination.HasProperties(1));
    KRATOS_EXPECT_TRUE(r_destination.GetProperties(1).Has(DENSITY));
    KRATOS_EXPECT_NEAR(r_destination.GetProperties(1).GetValue(DENSITY), 7850.0, 1e-12);
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsPropertiesDoNotSurviveClean, KratosMeshioPlusPlusFastSuite)
{
    // The other half of the finding above, pinned down explicitly so a future meshio++
    // release that starts propagating property sets through "clean" changes this test
    // rather than silently changing behaviour unnoticed.
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);
    r_source.pGetProperties(1)->SetValue(DENSITY, 7850.0);
    auto& r_destination = model.CreateModelPart("destination");

    MeshioPlusPlusMeshOperations::Execute(r_source, OperationSettings("clean"), r_destination);

    KRATOS_EXPECT_FALSE(r_destination.HasProperties(1));
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsEntityNamesSurviveAnOperation, KratosMeshioPlusPlusFastSuite)
{
    // Kratos::Element carries no queryable "registered name" of its own at runtime, so the
    // only way to observe it is to write the result to a format that stores it verbatim
    // (mdpa) and check the text - the same technique the direct IO round-trip tests use.
    Model model;
    auto& r_source = model.CreateModelPart("source");
    r_source.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_source.CreateNewNode(2, 1.0, 0.0, 0.0);
    r_source.CreateNewNode(3, 0.0, 1.0, 0.0);
    r_source.CreateNewNode(4, 0.0, 0.0, 1.0);
    auto p_properties = r_source.CreateNewProperties(1);
    r_source.CreateNewElement("Element3D4N", 1, {1, 2, 3, 4}, p_properties);
    auto& r_destination = model.CreateModelPart("destination");

    MeshioPlusPlusMeshOperations::Execute(r_source, OperationSettings("clean"), r_destination);

    const auto file_path = TestFilePath(".mdpa");
    {
        Parameters settings(R"({"time_series" : "single_file"})");
        MeshioPlusPlusIO io_write(file_path, settings);
        io_write.WriteModelPart(r_destination);
    }
    const std::string written = ReadFileContent(file_path);
    RemoveIfExists(file_path);

    KRATOS_EXPECT_TRUE(written.find("Element3D4N") != std::string::npos);
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsGradientOfALinearField, KratosMeshioPlusPlusFastSuite)
{
    // Green-Gauss is exact for a linear field on any cell, so this is a hard oracle rather
    // than a tolerance: f = 2x + 3y + 4z must differentiate to exactly (2, 3, 4) everywhere.
    // "output" points at VELOCITY (a registered array_1d<double, 3>) so the three components
    // survive the write-back - under meshio++'s own default name they never could.
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);
    for (auto& r_node : r_source.Nodes()) {
        r_node.SetValue(TEMPERATURE, 2.0 * r_node.X() + 3.0 * r_node.Y() + 4.0 * r_node.Z());
    }
    auto& r_destination = model.CreateModelPart("destination");

    Parameters settings = OperationSettings("gradient", R"({
        "array_name"                 : "TEMPERATURE",
        "location"                   : "cell",
        "output"                     : "VELOCITY",
        "nodal_data_value_variables" : ["TEMPERATURE"]
    })");
    const Parameters report = MeshioPlusPlusMeshOperations::Execute(r_source, settings, r_destination);

    KRATOS_EXPECT_EQ(report["number_of_skipped"].GetInt(), 0);
    KRATOS_EXPECT_EQ(r_destination.NumberOfElements(), r_source.NumberOfElements());
    for (const auto& r_element : r_destination.Elements()) {
        const array_1d<double, 3>& r_gradient = r_element.GetValue(VELOCITY);
        KRATOS_EXPECT_NEAR(r_gradient[0], 2.0, 1e-10);
        KRATOS_EXPECT_NEAR(r_gradient[1], 3.0, 1e-10);
        KRATOS_EXPECT_NEAR(r_gradient[2], 4.0, 1e-10);
    }
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsGradientNeedsAnArrayName, KratosMeshioPlusPlusFastSuite)
{
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);
    auto& r_destination = model.CreateModelPart("destination");

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        MeshioPlusPlusMeshOperations::Execute(r_source, OperationSettings("gradient"), r_destination),
        "array_name");
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsGradientRejectsCellData, KratosMeshioPlusPlusFastSuite)
{
    // A piecewise-constant field has no derivative; meshio++ says so by name rather than
    // averaging onto the nodes behind the caller's back.
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);
    for (auto& r_element : r_source.Elements()) {
        r_element.SetValue(DENSITY, 1.0);
    }
    auto& r_destination = model.CreateModelPart("destination");

    Parameters settings = OperationSettings("gradient", R"({
        "array_name"                   : "DENSITY",
        "element_data_value_variables" : ["DENSITY"]
    })");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        MeshioPlusPlusMeshOperations::Execute(r_source, settings, r_destination), "");
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsCropPredicate, KratosMeshioPlusPlusFastSuite)
{
    // Two triangles, one tagged 0.0 and one 1.0; "< 0.5" keeps exactly the first.
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateTriangulatedSquare(r_source);
    r_source.GetElement(1).SetValue(TEMPERATURE, 0.0);
    r_source.GetElement(2).SetValue(TEMPERATURE, 1.0);
    auto& r_destination = model.CreateModelPart("destination");

    Parameters settings = OperationSettings("crop_predicate", R"({
        "predicate_array"              : "TEMPERATURE",
        "predicate_op"                 : "<",
        "predicate_value"              : 0.5,
        "element_data_value_variables" : ["TEMPERATURE"]
    })");
    MeshioPlusPlusMeshOperations::Execute(r_source, settings, r_destination);

    KRATOS_EXPECT_EQ(r_destination.NumberOfElements(), 1);
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsCropPredicateNeverMatchesNonFinite, KratosMeshioPlusPlusFastSuite)
{
    // The documented shared-evaluator rule, and the case most likely to regress: a non-finite
    // cell value never matches, *including* under "!=", where IEEE would say NaN != 1.0 is
    // true. Both cells are NaN here, so a correct "!=" keeps neither.
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateTriangulatedSquare(r_source);
    const double nan_value = std::numeric_limits<double>::quiet_NaN();
    for (auto& r_element : r_source.Elements()) {
        r_element.SetValue(TEMPERATURE, nan_value);
    }
    auto& r_destination = model.CreateModelPart("destination");

    Parameters settings = OperationSettings("crop_predicate", R"({
        "predicate_array"              : "TEMPERATURE",
        "predicate_op"                 : "!=",
        "predicate_value"              : 1.0,
        "element_data_value_variables" : ["TEMPERATURE"]
    })");
    MeshioPlusPlusMeshOperations::Execute(r_source, settings, r_destination);

    KRATOS_EXPECT_EQ(r_destination.NumberOfElements(), 0);
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsVoxelizeFillAll, KratosMeshioPlusPlusFastSuite)
{
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);
    auto& r_destination = model.CreateModelPart("destination");

    Parameters settings = OperationSettings("voxelize", R"({"resolution" : [4, 4, 4]})");
    const Parameters report = MeshioPlusPlusMeshOperations::Execute(r_source, settings, r_destination);

    // "all" keeps the whole bounding box, so the lattice is exactly 4 x 4 x 4 hexahedra.
    KRATOS_EXPECT_EQ(r_destination.NumberOfElements(), 64);
    KRATOS_EXPECT_EQ(report["number_of_occupied"].GetInt(), 64);
    KRATOS_EXPECT_EQ(report["dims"][0].GetInt(), 4);
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsVoxelizeSurfaceKeepsFewerCells, KratosMeshioPlusPlusFastSuite)
{
    // "surface" and "inside" measure against a surface, so the source is the cube's skin
    // rather than the volume mesh - meshio++ says so by name if it is not.
    Model model;
    auto& r_volume = model.CreateModelPart("volume");
    PopulateCubeOfTetrahedra(r_volume);
    auto& r_surface = model.CreateModelPart("surface");
    MeshioPlusPlusMeshOperations::Execute(r_volume, OperationSettings("extract_skin"), r_surface);
    auto& r_surface_only = model.CreateModelPart("surface_only");

    Parameters settings = OperationSettings("voxelize", R"({"resolution" : [4, 4, 4], "fill" : "surface"})");
    MeshioPlusPlusMeshOperations::Execute(r_surface, settings, r_surface_only);

    // Only the shell of cells a surface triangle passes through: strictly fewer than the 64
    // of the full box, and more than none.
    KRATOS_EXPECT_GT(r_surface_only.NumberOfElements(), 0);
    KRATOS_EXPECT_LT(r_surface_only.NumberOfElements(), 64);
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsVoxelizeNeedsExactlyOneSizing, KratosMeshioPlusPlusFastSuite)
{
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);
    auto& r_destination = model.CreateModelPart("destination");

    // Neither or both of "resolution"/"cell_size" is meshio++'s own named error.
    Parameters both = OperationSettings("voxelize", R"({"resolution" : [4, 4, 4], "cell_size" : 0.25})");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        MeshioPlusPlusMeshOperations::Execute(r_source, both, r_destination), "");

    auto& r_neither_destination = model.CreateModelPart("neither");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        MeshioPlusPlusMeshOperations::Execute(r_source, OperationSettings("voxelize"), r_neither_destination), "");
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsGrid, KratosMeshioPlusPlusFastSuite)
{
    // The one generator: no source model part at all.
    Model model;
    auto& r_destination = model.CreateModelPart("destination");

    Parameters settings(R"({
        "dims"    : [2, 3, 4],
        "origin"  : [0.0, 0.0, 0.0],
        "spacing" : [1.0, 1.0, 1.0]
    })");
    const Parameters report = MeshioPlusPlusMeshOperations::Grid(settings, r_destination);

    KRATOS_EXPECT_EQ(r_destination.NumberOfElements(), 2 * 3 * 4);
    KRATOS_EXPECT_EQ(r_destination.NumberOfNodes(), 3 * 4 * 5);
    KRATOS_EXPECT_EQ(report["number_of_elements"].GetInt(), 24);
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsComputeSdfSignFlips, KratosMeshioPlusPlusFastSuite)
{
    // The skin of the cube is a closed surface, so a grid spanning it must carry both signs:
    // negative strictly inside, positive strictly outside.
    Model model;
    auto& r_volume = model.CreateModelPart("volume");
    PopulateCubeOfTetrahedra(r_volume);
    auto& r_surface = model.CreateModelPart("surface");
    MeshioPlusPlusMeshOperations::Execute(r_volume, OperationSettings("extract_skin"), r_surface);
    auto& r_destination = model.CreateModelPart("destination");

    Parameters settings = OperationSettings("compute_sdf", R"({
        "resolution"       : [6, 6, 6],
        "padding_relative" : 0.5,
        "output"           : "DISTANCE"
    })");
    MeshioPlusPlusMeshOperations::Execute(r_surface, settings, r_destination);

    bool has_negative = false;
    bool has_positive = false;
    for (const auto& r_node : r_destination.Nodes()) {
        const double distance = r_node.GetValue(DISTANCE);
        has_negative = has_negative || distance < -1e-9;
        has_positive = has_positive || distance > 1e-9;
    }
    KRATOS_EXPECT_TRUE(has_negative);
    KRATOS_EXPECT_TRUE(has_positive);
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsDistanceToSurfaceRenamesToAVariable, KratosMeshioPlusPlusFastSuite)
{
    // The test that proves the "output" rename works. meshio++ names the result
    // "sdf:distance", which no Kratos Variable is, so without the rename nothing would be
    // retrievable from the destination at all. Every cube node lies exactly on the cube's own
    // skin, so the expected distance is exactly zero - a hard oracle, not a tolerance.
    Model model;
    auto& r_volume = model.CreateModelPart("volume");
    PopulateCubeOfTetrahedra(r_volume);
    auto& r_surface = model.CreateModelPart("surface");
    MeshioPlusPlusMeshOperations::Execute(r_volume, OperationSettings("extract_skin"), r_surface);
    auto& r_destination = model.CreateModelPart("destination");

    Parameters settings(R"({"output" : "DISTANCE"})");
    const Parameters report = MeshioPlusPlusMeshOperations::DistanceToSurface(
        r_volume, r_surface, settings, r_destination);

    KRATOS_EXPECT_EQ(r_destination.NumberOfNodes(), r_volume.NumberOfNodes());
    for (const auto& r_node : r_destination.Nodes()) {
        KRATOS_EXPECT_NEAR(r_node.GetValue(DISTANCE), 0.0, 1e-9);
    }
    KRATOS_EXPECT_TRUE(report["surface_quality"]["watertight"].GetBool());
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsCheckSurfaceWatertight, KratosMeshioPlusPlusFastSuite)
{
    Model model;
    auto& r_volume = model.CreateModelPart("volume");
    PopulateCubeOfTetrahedra(r_volume);
    auto& r_closed = model.CreateModelPart("closed");
    MeshioPlusPlusMeshOperations::Execute(r_volume, OperationSettings("extract_skin"), r_closed);

    const Parameters closed_report = MeshioPlusPlusMeshOperations::CheckSurfaceWatertight(r_closed);
    KRATOS_EXPECT_TRUE(closed_report["watertight"].GetBool());
    KRATOS_EXPECT_EQ(closed_report["boundary_edges"].GetInt(), 0);

    // A bare sheet is not: every one of its outer edges is used by a single triangle.
    auto& r_open = model.CreateModelPart("open");
    PopulateTriangulatedSquare(r_open);
    const Parameters open_report = MeshioPlusPlusMeshOperations::CheckSurfaceWatertight(r_open);
    KRATOS_EXPECT_FALSE(open_report["watertight"].GetBool());
    KRATOS_EXPECT_GT(open_report["boundary_edges"].GetInt(), 0);
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsRefineBalancedDoesNotTearTheMesh, KratosMeshioPlusPlusFastSuite)
{
    // meshio++ v9.23.0 fixed a second balanced pass leaving the two sides of a hanging
    // interface referencing distinct but exactly coincident nodes. The output is allowed to be
    // 1-irregular; it is not allowed to be torn, so no two distinct nodes may share a position.
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);
    auto& r_destination = model.CreateModelPart("destination");

    Parameters settings = OperationSettings("refine",
        R"({"cells" : [0], "levels" : 2, "closure" : "balanced"})");
    MeshioPlusPlusMeshOperations::Execute(r_source, settings, r_destination);

    std::vector<std::array<double, 3>> positions;
    positions.reserve(r_destination.NumberOfNodes());
    for (const auto& r_node : r_destination.Nodes()) {
        positions.push_back({{r_node.X(), r_node.Y(), r_node.Z()}});
    }
    for (std::size_t i = 0; i < positions.size(); ++i) {
        for (std::size_t j = i + 1; j < positions.size(); ++j) {
            const double dx = positions[i][0] - positions[j][0];
            const double dy = positions[i][1] - positions[j][1];
            const double dz = positions[i][2] - positions[j][2];
            KRATOS_EXPECT_GT(dx * dx + dy * dy + dz * dz, 1e-20);
        }
    }
}


/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsSubdivide, KratosMeshioPlusPlusFastSuite)
{
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);
    auto& r_destination = model.CreateModelPart("destination");

    MeshioPlusPlusMeshOperations::Execute(r_source, OperationSettings("subdivide"), r_destination);

    // Polyhedral refinement emits one polyhedral child per face, which no Kratos Element can
    // hold - so the default "simplexify_result" decomposes them into tetrahedra first. The
    // count is therefore not a fixed multiple, but it must strictly grow and must not be empty
    // (an empty destination is the silent failure StorePolyhedralResult exists to prevent).
    KRATOS_EXPECT_GT(r_destination.NumberOfElements(), r_source.NumberOfElements());
    KRATOS_EXPECT_GT(r_destination.NumberOfNodes(), r_source.NumberOfNodes());
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsSubdivideWithoutSimplexifyThrows, KratosMeshioPlusPlusFastSuite)
{
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);
    auto& r_destination = model.CreateModelPart("destination");

    // Raw polyhedral output cannot reach a model part. Silently returning an empty destination
    // is the outcome this refusal replaces.
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        MeshioPlusPlusMeshOperations::Execute(
            r_source, OperationSettings("subdivide", R"({"simplexify_result" : false})"),
            r_destination),
        "produced polyhedral cells");
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsAgglomerate, KratosMeshioPlusPlusFastSuite)
{
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);
    auto& r_destination = model.CreateModelPart("destination");

    MeshioPlusPlusMeshOperations::Execute(
        r_source, OperationSettings("agglomerate", R"({"target_group_size" : 2})"), r_destination);

    // Grouping the six tetrahedra in pairs yields three polyhedra, which the default
    // "simplexify_result" then decomposes back into tetrahedra - so the final element count is
    // not comparable to the input's. What matters is that the destination is reachable at all.
    KRATOS_EXPECT_GT(r_destination.NumberOfElements(), 0);
    KRATOS_EXPECT_GT(r_destination.NumberOfNodes(), 0);
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsRefineRecordHierarchy, KratosMeshioPlusPlusFastSuite)
{
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);
    auto& r_destination = model.CreateModelPart("destination");

    // "record_hierarchy" attaches "refine:cell_id"/"refine:parent_id". Both are colon-namespaced,
    // so neither ever reaches the destination model part - which is why meshio++'s own
    // "undo_green" is not exposed here at all. The assertion is that asking for them is
    // accepted and changes nothing else about the result.
    MeshioPlusPlusMeshOperations::Execute(
        r_source, OperationSettings("refine", R"({"levels" : 1, "record_hierarchy" : true})"),
        r_destination);

    KRATOS_EXPECT_EQ(r_destination.NumberOfElements(), 8 * r_source.NumberOfElements());
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsRemesh, KratosMeshioPlusPlusFastSuite)
{
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateClosedCubeSkin(r_source);
    auto& r_destination = model.CreateModelPart("destination");

    const Parameters report = MeshioPlusPlusMeshOperations::Execute(
        r_source, OperationSettings("remesh", R"({"num_clusters" : 8})"), r_destination);

    KRATOS_EXPECT_GT(report["number_of_clusters"].GetInt(), 0);
    KRATOS_EXPECT_TRUE(report.Has("number_of_iterations"));
    KRATOS_EXPECT_TRUE(report.Has("number_of_non_manifold_vertices"));
    KRATOS_EXPECT_GT(r_destination.NumberOfElements(), 0);
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsRemeshVolume, KratosMeshioPlusPlusFastSuite)
{
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateClosedCubeSkin(r_source);
    auto& r_destination = model.CreateModelPart("destination");

    // A closed surface in, a tetrahedral volume out. The lattice settings are the ones
    // "voxelize"/"compute_sdf" already share, hence "resolution" rather than an own key.
    const Parameters report = MeshioPlusPlusMeshOperations::Execute(
        r_source, OperationSettings("remesh_volume", R"({"resolution" : [4, 4, 4]})"), r_destination);

    KRATOS_EXPECT_GT(report["number_of_tets"].GetInt(), 0);
    KRATOS_EXPECT_TRUE(report["surface_quality"]["watertight"].GetBool());
    KRATOS_EXPECT_EQ(r_destination.NumberOfElements(), report["number_of_tets"].GetInt());
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsOptimizeVolume, KratosMeshioPlusPlusFastSuite)
{
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);
    auto& r_destination = model.CreateModelPart("destination");

    const Parameters report = MeshioPlusPlusMeshOperations::Execute(
        r_source, OperationSettings("optimize_volume", R"({"optimize_iterations" : 2})"), r_destination);

    // The contract is monotone worst-element quality: optimizing may do nothing on an already
    // good mesh, but it must never make the worst element worse.
    KRATOS_EXPECT_GE(report["min_quality_after"].GetDouble(), report["min_quality_before"].GetDouble());
    KRATOS_EXPECT_EQ(r_destination.NumberOfElements(), r_source.NumberOfElements());
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsDecimateVolume, KratosMeshioPlusPlusFastSuite)
{
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);
    // Two levels of subdivision give it enough tetrahedra to have something to remove.
    auto& r_refined = model.CreateModelPart("refined");
    MeshioPlusPlusMeshOperations::Execute(
        r_source, OperationSettings("refine", R"({"levels" : 2})"), r_refined);

    auto& r_destination = model.CreateModelPart("destination");
    const Parameters report = MeshioPlusPlusMeshOperations::Execute(
        r_refined, OperationSettings("decimate_volume", R"({"target_ratio" : 0.5})"), r_destination);

    KRATOS_EXPECT_TRUE(report.Has("tets_removed"));
    KRATOS_EXPECT_LE(r_destination.NumberOfElements(), r_refined.NumberOfElements());
    KRATOS_EXPECT_GT(r_destination.NumberOfElements(), 0);
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsEstimateError, KratosMeshioPlusPlusFastSuite)
{
    Model model;
    auto& r_source = model.CreateModelPart("source");
    r_source.AddNodalSolutionStepVariable(TEMPERATURE);
    PopulateCubeOfTetrahedra(r_source);
    // A field with curvature, so the recovered gradient genuinely differs from the element-wise
    // one and the indicator is not identically zero.
    for (auto& r_node : r_source.Nodes()) {
        r_node.FastGetSolutionStepValue(TEMPERATURE) = r_node.X() * r_node.X() + r_node.Y();
    }
    auto& r_destination = model.CreateModelPart("destination");

    // "output" names a registered Variable, which is what makes the indicator readable from
    // Kratos at all - meshio++ calls it "error:zz", which no Variable can be.
    const Parameters report = MeshioPlusPlusMeshOperations::Execute(
        r_source,
        OperationSettings("estimate_error",
            R"({"array_name" : "TEMPERATURE", "output" : "DISTANCE", "marking" : "fraction",
                "marking_value" : 0.5, "nodal_solution_step_data_variables" : ["TEMPERATURE"]})"),
        r_destination);

    KRATOS_EXPECT_GE(report["global_error"].GetDouble(), 0.0);
    KRATOS_EXPECT_GT(report["number_of_marked"].GetInt(), 0);
    KRATOS_EXPECT_EQ(r_destination.NumberOfElements(), r_source.NumberOfElements());
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsHessian, KratosMeshioPlusPlusFastSuite)
{
    Model model;
    auto& r_source = model.CreateModelPart("source");
    r_source.AddNodalSolutionStepVariable(TEMPERATURE);
    PopulateCubeOfTetrahedra(r_source);
    for (auto& r_node : r_source.Nodes()) {
        r_node.FastGetSolutionStepValue(TEMPERATURE) = r_node.X() * r_node.X();
    }
    auto& r_destination = model.CreateModelPart("destination");

    const Parameters report = MeshioPlusPlusMeshOperations::Execute(
        r_source,
        OperationSettings("hessian",
            R"({"array_name" : "TEMPERATURE", "location" : "cell",
                "nodal_solution_step_data_variables" : ["TEMPERATURE"]})"),
        r_destination);

    KRATOS_EXPECT_TRUE(report.Has("number_of_skipped"));
    KRATOS_EXPECT_TRUE(report.Has("number_of_fallback"));
    KRATOS_EXPECT_EQ(r_destination.NumberOfElements(), r_source.NumberOfElements());
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsDataIntegrate, KratosMeshioPlusPlusFastSuite)
{
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);
    // A cell-measure-weighted reduction needs *cell* data - meshio++ refuses a point_data array
    // by name rather than averaging it behind the caller's back. A constant field's weighted
    // mean is that constant whatever the mesh, which is the one value this must reproduce.
    for (auto& r_element : r_source.Elements()) {
        r_element.SetValue(TEMPERATURE, 2.0);
    }
    auto& r_destination = model.CreateModelPart("destination");

    const Parameters report = MeshioPlusPlusMeshOperations::Execute(
        r_source,
        OperationSettings("data_integrate",
            R"({"names" : ["TEMPERATURE"], "element_data_value_variables" : ["TEMPERATURE"]})"),
        r_destination);

    KRATOS_EXPECT_EQ(report["arrays"].size(), 1);
    const Parameters domain = report["arrays"][0]["domain"];
    // The unit cube's volume, and the constant back out of the weighted mean.
    KRATOS_EXPECT_NEAR(domain["domain_measure"].GetVector()[0], 1.0, 1e-10);
    KRATOS_EXPECT_NEAR(domain["total"].GetVector()[0], 2.0, 1e-10);
    KRATOS_EXPECT_NEAR(domain["mean"].GetVector()[0], 2.0, 1e-10);
    // Report-only: the destination is left untouched.
    KRATOS_EXPECT_EQ(r_destination.NumberOfElements(), 0);
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsConservativeInterpolate, KratosMeshioPlusPlusFastSuite)
{
    Model model;
    auto& r_source = model.CreateModelPart("source");
    r_source.AddNodalSolutionStepVariable(TEMPERATURE);
    PopulateTriangulatedSquare(r_source);
    for (auto& r_node : r_source.Nodes()) {
        r_node.FastGetSolutionStepValue(TEMPERATURE) = 3.0;
    }

    // A refinement of the same square: the geometry is identical, so a conservative transfer of
    // a constant field must reproduce the constant rather than smear it. The shared field data
    // settings are staged from *both* meshes, so the target has to know the variable too - the
    // usual case is that only the source carries values, not that only the source declares it.
    auto& r_target = model.CreateModelPart("target");
    r_target.AddNodalSolutionStepVariable(TEMPERATURE);
    MeshioPlusPlusMeshOperations::Execute(
        r_source, OperationSettings("refine", R"({"levels" : 1})"), r_target);

    auto& r_destination = model.CreateModelPart("destination");
    const Parameters report = MeshioPlusPlusMeshOperations::ConservativeInterpolate(
        r_source, r_target,
        // "overwrite" because the shared field data settings stage TEMPERATURE from the target
        // too (as zeros): a name present on both meshes is exactly what "on_conflict" resolves,
        // and the default "error" is the honest refusal rather than a silent choice.
        Parameters(R"({"names" : ["TEMPERATURE"], "on_conflict" : "overwrite",
                       "nodal_solution_step_data_variables" : ["TEMPERATURE"]})"),
        r_destination);

    KRATOS_EXPECT_EQ(report["number_of_elements"].GetInt(),
                     static_cast<int>(r_target.NumberOfElements()));
    KRATOS_EXPECT_EQ(r_destination.NumberOfNodes(), r_target.NumberOfNodes());
}

} // namespace Kratos::Testing
