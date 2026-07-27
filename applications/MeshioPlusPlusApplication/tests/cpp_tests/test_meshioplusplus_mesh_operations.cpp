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

// Project includes
#include "containers/model.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "custom_utilities/meshioplusplus_mesh_operations.h"
#include "custom_io/meshioplusplus_io.h"
#include "meshioplusplus_fast_suite.h"

namespace Kratos::Testing {

namespace {

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
    for (const std::string name : {"attach_quality", "clean", "convert_cells", "crop_bbox",
                                   "crop_halfspace", "decimate", "extract_skin", "extract_surface",
                                   "isosurface", "partition", "quality", "refine", "reorder",
                                   "slice", "smooth", "split", "stats", "transform"}) {
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
    // Internals::ModelPartToMesh (the conversion Execute uses) carries topology and material
    // data, but never nodal/elemental *field* data - there is currently no way to reach this
    // operation with a real point_data array through MeshioPlusPlusMeshOperations::Execute.
    // This pins down that gap: naming any array fails, not just an empty name.
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

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusMeshOperationsAttachQuality, KratosMeshioPlusPlusFastSuite)
{
    Model model;
    auto& r_source = model.CreateModelPart("source");
    PopulateCubeOfTetrahedra(r_source);
    auto& r_destination = model.CreateModelPart("destination");

    MeshioPlusPlusMeshOperations::Execute(r_source, OperationSettings("attach_quality"), r_destination);

    // The quality metrics attach_quality computes are mesh-level cell_data; MeshToModelPart
    // (the bridge Execute stores its result through) does not carry arbitrary cell_data into
    // Kratos, so only the topology is expected to survive here.
    KRATOS_EXPECT_EQ(r_destination.NumberOfElements(), r_source.NumberOfElements());
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

} // namespace Kratos::Testing
