//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes
#include <filesystem>
#include <fstream>
#include <sstream>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "input_output/meshioplusplus_io.h"
#include "testing/testing.h"

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

/// Creates 5 nodes, 2 tetrahedra elements and 1 triangle condition
void PopulateTetrahedraModelPart(ModelPart& rModelPart)
{
    auto p_properties = rModelPart.CreateNewProperties(1);
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
    rModelPart.CreateNewNode(4, 0.0, 0.0, 1.0);
    rModelPart.CreateNewNode(5, 1.0, 1.0, 1.0);
    rModelPart.CreateNewElement("Element3D4N", 1, {1, 2, 3, 4}, p_properties);
    rModelPart.CreateNewElement("Element3D4N", 2, {2, 3, 4, 5}, p_properties);
    rModelPart.CreateNewCondition("SurfaceCondition3D3N", 1, {1, 2, 3}, p_properties);
}

/// Creates 4 nodes and 2 triangle elements (a planar surface mesh)
void PopulateTriangleModelPart(ModelPart& rModelPart)
{
    auto p_properties = rModelPart.CreateNewProperties(1);
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
    rModelPart.CreateNewNode(4, 0.0, 1.0, 0.0);
    rModelPart.CreateNewElement("Element2D3N", 1, {1, 2, 3}, p_properties);
    rModelPart.CreateNewElement("Element2D3N", 2, {1, 3, 4}, p_properties);
}

void ExpectSameMesh(const ModelPart& rExpected, const ModelPart& rActual)
{
    KRATOS_EXPECT_EQ(rExpected.NumberOfNodes(), rActual.NumberOfNodes());
    KRATOS_EXPECT_EQ(rExpected.NumberOfElements(), rActual.NumberOfElements());
    KRATOS_EXPECT_EQ(rExpected.NumberOfConditions(), rActual.NumberOfConditions());

    auto it_expected = rExpected.NodesBegin();
    auto it_actual = rActual.NodesBegin();
    for (; it_expected != rExpected.NodesEnd(); ++it_expected, ++it_actual) {
        KRATOS_EXPECT_NEAR(it_expected->X(), it_actual->X(), 1e-12);
        KRATOS_EXPECT_NEAR(it_expected->Y(), it_actual->Y(), 1e-12);
        KRATOS_EXPECT_NEAR(it_expected->Z(), it_actual->Z(), 1e-12);
    }

    auto it_expected_elem = rExpected.ElementsBegin();
    auto it_actual_elem = rActual.ElementsBegin();
    for (; it_expected_elem != rExpected.ElementsEnd(); ++it_expected_elem, ++it_actual_elem) {
        const auto& r_expected_geom = it_expected_elem->GetGeometry();
        const auto& r_actual_geom = it_actual_elem->GetGeometry();
        KRATOS_EXPECT_EQ(r_expected_geom.size(), r_actual_geom.size());
        for (std::size_t i = 0; i < r_expected_geom.size(); ++i) {
            KRATOS_EXPECT_EQ(r_expected_geom[i].Id(), r_actual_geom[i].Id());
        }
    }
}

void WriteReadRoundTrip(const std::string& rExtension)
{
    Model model;
    auto& r_write_model_part = model.CreateModelPart("write");
    auto& r_read_model_part = model.CreateModelPart("read");
    PopulateTetrahedraModelPart(r_write_model_part);

    const auto file_path = TestFilePath(rExtension);
    {
        Parameters settings(R"({"time_series" : "single_file"})");
        MeshioPlusPlusIO io_write(file_path, settings);
        io_write.WriteModelPart(r_write_model_part);
    }
    {
        MeshioPlusPlusIO io_read(file_path);
        io_read.ReadModelPart(r_read_model_part);
    }
    RemoveIfExists(file_path);

    ExpectSameMesh(r_write_model_part, r_read_model_part);
}

std::string ReadFileContent(const std::filesystem::path& rPath)
{
    std::ifstream file(rPath);
    std::stringstream buffer;
    buffer << file.rdbuf();
    return buffer.str();
}

std::size_t CountOccurrences(const std::string& rContent, const std::string& rToken)
{
    std::size_t count = 0;
    std::size_t position = rContent.find(rToken);
    while (position != std::string::npos) {
        ++count;
        position = rContent.find(rToken, position + rToken.size());
    }
    return count;
}

} // namespace

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusIOFormatIntrospection, KratosCoreFastSuite)
{
    // The self-contained formats must always be present
    const auto supported_formats = MeshioPlusPlusIO::GetSupportedFormats();
    for (const std::string name : {"vtu", "vtk", "gmsh", "stl", "obj", "xdmf", "abaqus", "medit",
                                   "ensight", "vtp", "triangle", "svg", "tikz"}) {
        KRATOS_EXPECT_TRUE(std::find(supported_formats.begin(), supported_formats.end(), name) != supported_formats.end());
    }

    // openfoam is read-only
    const auto read_formats = MeshioPlusPlusIO::GetSupportedReadFormats();
    const auto write_formats = MeshioPlusPlusIO::GetSupportedWriteFormats();
    KRATOS_EXPECT_TRUE(std::find(read_formats.begin(), read_formats.end(), "openfoam") != read_formats.end());
    KRATOS_EXPECT_TRUE(std::find(write_formats.begin(), write_formats.end(), "openfoam") == write_formats.end());

    // svg and tikz are write-only
    for (const std::string name : {"svg", "tikz"}) {
        KRATOS_EXPECT_TRUE(std::find(write_formats.begin(), write_formats.end(), name) != write_formats.end());
        KRATOS_EXPECT_TRUE(std::find(read_formats.begin(), read_formats.end(), name) == read_formats.end());
    }

    // Name <-> enum round trips
    KRATOS_EXPECT_EQ(static_cast<int>(MeshioPlusPlusIO::FormatFromString("gmsh")),
                     static_cast<int>(MeshioPlusPlusIO::Format::GMSH));
    KRATOS_EXPECT_EQ(MeshioPlusPlusIO::FormatName(MeshioPlusPlusIO::Format::GMSH), "gmsh");
    KRATOS_EXPECT_EQ(static_cast<int>(MeshioPlusPlusIO::FormatFromString("auto")),
                     static_cast<int>(MeshioPlusPlusIO::Format::AUTOMATIC));

    // Name <-> enum round trips of the newer formats
    KRATOS_EXPECT_EQ(static_cast<int>(MeshioPlusPlusIO::FormatFromString("ensight")),
                     static_cast<int>(MeshioPlusPlusIO::Format::ENSIGHT));
    KRATOS_EXPECT_EQ(MeshioPlusPlusIO::FormatName(MeshioPlusPlusIO::Format::ENSIGHT), "ensight");
    KRATOS_EXPECT_EQ(static_cast<int>(MeshioPlusPlusIO::FormatFromString("vtp")),
                     static_cast<int>(MeshioPlusPlusIO::Format::VTP));
    KRATOS_EXPECT_EQ(MeshioPlusPlusIO::FormatName(MeshioPlusPlusIO::Format::TIKZ), "tikz");

    // Extension resolution
    KRATOS_EXPECT_EQ(static_cast<int>(MeshioPlusPlusIO::ResolveFormat("mesh.vtu")),
                     static_cast<int>(MeshioPlusPlusIO::Format::VTU));
    KRATOS_EXPECT_EQ(static_cast<int>(MeshioPlusPlusIO::ResolveFormat("mesh.msh")),
                     static_cast<int>(MeshioPlusPlusIO::Format::GMSH));
    KRATOS_EXPECT_EQ(static_cast<int>(MeshioPlusPlusIO::ResolveFormat("mesh.case")),
                     static_cast<int>(MeshioPlusPlusIO::Format::ENSIGHT));
    KRATOS_EXPECT_EQ(static_cast<int>(MeshioPlusPlusIO::ResolveFormat("mesh.vtp")),
                     static_cast<int>(MeshioPlusPlusIO::Format::VTP));
    // .node/.ele keep resolving to tetgen; triangle needs the "format" setting
    KRATOS_EXPECT_EQ(static_cast<int>(MeshioPlusPlusIO::ResolveFormat("mesh.node")),
                     static_cast<int>(MeshioPlusPlusIO::Format::TETGEN));
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(MeshioPlusPlusIO::ResolveFormat("mesh.unknown_extension"),
                                      "Cannot resolve a format");

    // Unknown names throw with the supported list
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(MeshioPlusPlusIO::FormatFromString("not_a_format"),
                                      "Unknown format");

    // Availability is consistent: self-contained always, optional ones report their state
    KRATOS_EXPECT_TRUE(MeshioPlusPlusIO::IsFormatAvailable(MeshioPlusPlusIO::Format::VTU));
    KRATOS_EXPECT_TRUE(MeshioPlusPlusIO::IsFormatAvailable(MeshioPlusPlusIO::Format::AUTOMATIC));
}

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusIOWriteReadVtu, KratosCoreFastSuite)
{
    WriteReadRoundTrip(".vtu");
}

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusIOWriteReadGmsh, KratosCoreFastSuite)
{
    WriteReadRoundTrip(".msh");
}

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusIOWriteReadMed, KratosCoreFastSuite)
{
    if (!MeshioPlusPlusIO::IsFormatAvailable(MeshioPlusPlusIO::Format::MED)) {
        GTEST_SKIP() << "med requires an HDF5-enabled build";
    }
    WriteReadRoundTrip(".med");
}

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusIOWriteReadEnsight, KratosCoreFastSuite)
{
    Model model;
    auto& r_write_model_part = model.CreateModelPart("write");
    auto& r_read_model_part = model.CreateModelPart("read");
    PopulateTetrahedraModelPart(r_write_model_part);

    // EnSight writes a .case file plus a sibling .geo file
    const auto file_path = TestFilePath(".case");
    auto geo_path = file_path;
    geo_path.replace_extension(".geo");
    {
        Parameters settings(R"({"time_series" : "single_file"})");
        MeshioPlusPlusIO io_write(file_path, settings);
        io_write.WriteModelPart(r_write_model_part);
    }
    KRATOS_EXPECT_TRUE(std::filesystem::exists(geo_path));
    {
        MeshioPlusPlusIO io_read(file_path);
        io_read.ReadModelPart(r_read_model_part);
    }
    RemoveIfExists(file_path);
    RemoveIfExists(geo_path);

    ExpectSameMesh(r_write_model_part, r_read_model_part);
}

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusIOWriteReadVtp, KratosCoreFastSuite)
{
    // vtp is a surface (PolyData) format: round trip a triangle mesh
    Model model;
    auto& r_write_model_part = model.CreateModelPart("write");
    auto& r_read_model_part = model.CreateModelPart("read");
    PopulateTriangleModelPart(r_write_model_part);

    const auto file_path = TestFilePath(".vtp");
    {
        Parameters settings(R"({"time_series" : "single_file"})");
        MeshioPlusPlusIO io_write(file_path, settings);
        io_write.WriteModelPart(r_write_model_part);
    }
    {
        MeshioPlusPlusIO io_read(file_path);
        io_read.ReadModelPart(r_read_model_part);
    }
    RemoveIfExists(file_path);

    ExpectSameMesh(r_write_model_part, r_read_model_part);
}

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusIOTriangleFormat, KratosCoreFastSuite)
{
    // Shewchuk Triangle .node/.ele pair: 4 nodes, 2 triangles (1-based indices).
    // The .node/.ele extensions resolve to tetgen, so the format is explicit.
    const auto node_path = TestFilePath(".node");
    auto ele_path = node_path;
    ele_path.replace_extension(".ele");
    {
        std::ofstream node_file(node_path);
        node_file << R"(4 2 0 0
1 0.0 0.0
2 1.0 0.0
3 1.0 1.0
4 0.0 1.0
)";
        std::ofstream ele_file(ele_path);
        ele_file << R"(2 3 0
1 1 2 3
2 1 3 4
)";
    }

    Model model;
    auto& r_model_part = model.CreateModelPart("read");
    {
        Parameters settings(R"({"format" : "triangle"})");
        MeshioPlusPlusIO io_read(node_path, settings);
        io_read.ReadModelPart(r_model_part);
    }
    RemoveIfExists(node_path);
    RemoveIfExists(ele_path);

    KRATOS_EXPECT_EQ(r_model_part.NumberOfNodes(), 4);
    KRATOS_EXPECT_EQ(r_model_part.NumberOfElements(), 2);
    for (const auto& r_node : r_model_part.Nodes()) {
        KRATOS_EXPECT_NEAR(r_node.Z(), 0.0, 1e-12); // 2D points are zero-padded
    }

    // Writing is not supported from Kratos (nodes are staged as 3D points)
    auto& r_write_model_part = model.CreateModelPart("write");
    PopulateTriangleModelPart(r_write_model_part);
    const auto write_path = TestFilePath("_write.node");
    Parameters write_settings(R"({"format" : "triangle", "time_series" : "single_file"})");
    MeshioPlusPlusIO io_write(write_path, write_settings);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(io_write.WriteModelPart(r_write_model_part),
                                      "can only write 2D points");
}

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusIOWriteSvgAndTikz, KratosCoreFastSuite)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("write");
    PopulateTriangleModelPart(r_model_part);

    for (const std::string extension : {".svg", ".tikz"}) {
        const auto file_path = TestFilePath(extension);
        {
            Parameters settings(R"({"time_series" : "single_file"})");
            MeshioPlusPlusIO io_write(file_path, settings);
            io_write.WriteModelPart(r_model_part);
        }
        KRATOS_EXPECT_TRUE(std::filesystem::exists(file_path));
        KRATOS_EXPECT_TRUE(std::filesystem::file_size(file_path) > 0);

        // svg/tikz are write-only
        auto& r_read = model.CreateModelPart("read" + extension.substr(1));
        MeshioPlusPlusIO io_read(file_path);
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(io_read.ReadModelPart(r_read), "write-only");
        RemoveIfExists(file_path);
    }
}

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusIOStlSkinOfVolumeMesh, KratosCoreFastSuite)
{
    // 2 tetrahedra sharing one face: the boundary skin holds 6 triangles
    Model model;
    auto& r_model_part = model.CreateModelPart("write");
    PopulateTetrahedraModelPart(r_model_part);

    const auto file_path = TestFilePath(".stl");
    auto count_read_triangles = [&model, &file_path](const std::string& rName) {
        auto& r_read = model.CreateModelPart(rName);
        MeshioPlusPlusIO(file_path).ReadModelPart(r_read);
        return r_read.NumberOfElements();
    };

    { // default: the extracted boundary skin of the volume mesh is written
        Parameters settings(R"({"time_series" : "single_file"})");
        MeshioPlusPlusIO(file_path, settings).WriteModelPart(r_model_part);
        KRATOS_EXPECT_EQ(count_read_triangles("read_skin"), 6);
    }
    { // the skin path is also honored by the ascii/binary override
        Parameters settings(R"({"time_series" : "single_file", "file_format" : "binary"})");
        MeshioPlusPlusIO(file_path, settings).WriteModelPart(r_model_part);
        KRATOS_EXPECT_EQ(count_read_triangles("read_skin_binary"), 6);
    }
    { // "skin" off: legacy behavior, only the existing triangle condition is written
        Parameters settings(R"({"time_series" : "single_file", "skin" : false})");
        MeshioPlusPlusIO(file_path, settings).WriteModelPart(r_model_part);
        KRATOS_EXPECT_EQ(count_read_triangles("read_no_skin"), 1);
    }
    RemoveIfExists(file_path);
}

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusIOReadGmshPhysicalTagsAsSubModelParts, KratosCoreFastSuite)
{
    // Gmsh 2.2 ASCII: one tetrahedron with physical tag 1, one triangle with
    // physical tag 2. The physical tags must become sub model parts.
    const auto file_path = TestFilePath(".msh");
    {
        std::ofstream file(file_path);
        file << R"($MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
4
1 0 0 0
2 1 0 0
3 0 1 0
4 0 0 1
$EndNodes
$Elements
2
1 4 2 1 1 1 2 3 4
2 2 2 2 2 1 2 3
$EndElements
)";
    }

    Model model;
    auto& r_model_part = model.CreateModelPart("read");
    {
        MeshioPlusPlusIO io_read(file_path);
        io_read.ReadModelPart(r_model_part);
    }
    RemoveIfExists(file_path);

    KRATOS_EXPECT_EQ(r_model_part.NumberOfNodes(), 4);
    KRATOS_EXPECT_EQ(r_model_part.NumberOfElements(), 1);   // the tetrahedron (highest dimension)
    KRATOS_EXPECT_EQ(r_model_part.NumberOfConditions(), 1); // the triangle

    KRATOS_EXPECT_TRUE(r_model_part.HasSubModelPart("gmsh_physical_1"));
    KRATOS_EXPECT_TRUE(r_model_part.HasSubModelPart("gmsh_physical_2"));
    const auto& r_volume = r_model_part.GetSubModelPart("gmsh_physical_1");
    const auto& r_surface = r_model_part.GetSubModelPart("gmsh_physical_2");
    KRATOS_EXPECT_EQ(r_volume.NumberOfElements(), 1);
    KRATOS_EXPECT_EQ(r_volume.NumberOfNodes(), 4);
    KRATOS_EXPECT_EQ(r_surface.NumberOfConditions(), 1);
    KRATOS_EXPECT_EQ(r_surface.NumberOfNodes(), 3);
}

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusIOFileSeries, KratosCoreFastSuite)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("write");
    PopulateTetrahedraModelPart(r_model_part);
    auto& r_process_info = r_model_part.GetProcessInfo();

    const auto file_path = TestFilePath(".vtu");
    auto labeled_path = [&file_path](const std::string& rLabel) {
        auto path = file_path;
        path.replace_extension();
        path += "_" + rLabel + ".vtu";
        return path;
    };

    {
        MeshioPlusPlusIO io_write(file_path); // "time_series" default: automatic -> file series for vtu
        r_process_info[STEP] = 1;
        io_write.WriteModelPart(r_model_part);
        r_process_info[STEP] = 2;
        io_write.WriteModelPart(r_model_part);
    }

    KRATOS_EXPECT_TRUE(std::filesystem::exists(labeled_path("1")));
    KRATOS_EXPECT_TRUE(std::filesystem::exists(labeled_path("2")));
    KRATOS_EXPECT_FALSE(std::filesystem::exists(file_path)); // only labeled files are written

    RemoveIfExists(labeled_path("1"));
    RemoveIfExists(labeled_path("2"));
}

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusIOSingleFileOverwrites, KratosCoreFastSuite)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("write");
    PopulateTetrahedraModelPart(r_model_part);

    const auto file_path = TestFilePath(".vtu");
    {
        Parameters settings(R"({"time_series" : "single_file"})");
        MeshioPlusPlusIO io_write(file_path, settings);
        io_write.WriteModelPart(r_model_part);
        io_write.WriteModelPart(r_model_part); // overwrites the same file
    }
    KRATOS_EXPECT_TRUE(std::filesystem::exists(file_path));
    RemoveIfExists(file_path);
}

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusIOXdmfTimeSeriesAppend, KratosCoreFastSuite)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("write");
    r_model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    PopulateTetrahedraModelPart(r_model_part);
    auto& r_process_info = r_model_part.GetProcessInfo();

    const auto file_path = TestFilePath(".xdmf");
    Parameters settings(R"({
        "output_control_type"                : "time",
        "xdmf_data_format"                   : "XML",
        "nodal_solution_step_data_variables" : ["TEMPERATURE"]
    })");

    {
        MeshioPlusPlusIO io_write(file_path, settings.Clone());
        for (int step = 1; step <= 3; ++step) {
            r_process_info[TIME] = 0.1 * step;
            for (auto& r_node : r_model_part.Nodes()) {
                r_node.FastGetSolutionStepValue(TEMPERATURE) = step * r_node.Id();
            }
            io_write.WriteModelPart(r_model_part); // first call writes the mesh, later ones append
        }
    }

    std::string content = ReadFileContent(file_path);
    KRATOS_EXPECT_EQ(CountOccurrences(content, "<Time Value="), 3);
    KRATOS_EXPECT_EQ(CountOccurrences(content, "Name=\"mesh\""), 1);
    KRATOS_EXPECT_EQ(CountOccurrences(content, "Name=\"TEMPERATURE\""), 3);
    KRATOS_EXPECT_EQ(CountOccurrences(content, "CollectionType=\"Temporal\""), 1);

    { // A NEW IO on the same file must detect the existing "buffer" and extend it
        MeshioPlusPlusIO io_append(file_path, settings.Clone());
        r_process_info[TIME] = 0.4;
        io_append.WriteModelPart(r_model_part);
    }

    content = ReadFileContent(file_path);
    KRATOS_EXPECT_EQ(CountOccurrences(content, "<Time Value="), 4);
    KRATOS_EXPECT_EQ(CountOccurrences(content, "Name=\"mesh\""), 1); // the mesh is not rewritten

    RemoveIfExists(file_path);
}

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusIOWritePointData, KratosCoreFastSuite)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("write");
    PopulateTetrahedraModelPart(r_model_part);
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.SetValue(PRESSURE, 2.0 * r_node.Id());
        r_node.SetValue(VELOCITY, array_1d<double, 3>{1.0, 2.0, 3.0});
    }

    const auto file_path = TestFilePath(".vtu");
    {
        Parameters settings(R"({
            "time_series"                : "single_file",
            "nodal_data_value_variables" : ["PRESSURE", "VELOCITY"]
        })");
        MeshioPlusPlusIO io_write(file_path, settings);
        io_write.WriteModelPart(r_model_part);
    }

    // The variable names must appear as point data arrays in the vtu file
    const std::string content = ReadFileContent(file_path);
    KRATOS_EXPECT_TRUE(content.find("PRESSURE") != std::string::npos);
    KRATOS_EXPECT_TRUE(content.find("VELOCITY") != std::string::npos);
    RemoveIfExists(file_path);
}

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusIOInvalidSettings, KratosCoreFastSuite)
{
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        MeshioPlusPlusIO("mesh.vtu", Parameters(R"({"format" : "not_a_format"})")),
        "Unknown format");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        MeshioPlusPlusIO("mesh.vtu", Parameters(R"({"time_series" : "wrong"})")),
        "Unknown \"time_series\" setting");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        MeshioPlusPlusIO("mesh.vtu", Parameters(R"({"entity_type" : "wrong"})")),
        "Unknown \"entity_type\" setting");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        MeshioPlusPlusIO("mesh.vtu", Parameters(R"({"file_format" : "wrong"})")),
        "Unknown \"file_format\" setting");

    // Reading a non-existent file must fail with a clear message
    Model model;
    auto& r_model_part = model.CreateModelPart("read");
    MeshioPlusPlusIO io_read("definitely_missing_file.vtu");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(io_read.ReadModelPart(r_model_part), "does not exist");
}

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusIOEntityType, KratosCoreFastSuite)
{
    Model model;
    auto& r_write_model_part = model.CreateModelPart("write");
    PopulateTetrahedraModelPart(r_write_model_part); // 2 tetrahedra + 1 triangle condition
    const auto file_path = TestFilePath(".vtu");

    { // only elements: the read-back file holds the 2 tetrahedra
        Parameters settings(R"({"time_series" : "single_file", "entity_type" : "element"})");
        MeshioPlusPlusIO(file_path, settings).WriteModelPart(r_write_model_part);
        auto& r_read = model.CreateModelPart("read_element");
        MeshioPlusPlusIO(file_path).ReadModelPart(r_read);
        KRATOS_EXPECT_EQ(r_read.NumberOfElements(), 2);
        KRATOS_EXPECT_EQ(r_read.NumberOfConditions(), 0);
    }
    { // only conditions: the file holds the triangle (its highest dimension)
        Parameters settings(R"({"time_series" : "single_file", "entity_type" : "condition"})");
        MeshioPlusPlusIO(file_path, settings).WriteModelPart(r_write_model_part);
        auto& r_read = model.CreateModelPart("read_condition");
        MeshioPlusPlusIO(file_path).ReadModelPart(r_read);
        KRATOS_EXPECT_EQ(r_read.NumberOfElements(), 1); // the triangle is now the max dimension
        KRATOS_EXPECT_EQ(r_read.NumberOfConditions(), 0);
    }
    { // automatic (default): both
        Parameters settings(R"({"time_series" : "single_file"})");
        MeshioPlusPlusIO(file_path, settings).WriteModelPart(r_write_model_part);
        auto& r_read = model.CreateModelPart("read_automatic");
        MeshioPlusPlusIO(file_path).ReadModelPart(r_read);
        KRATOS_EXPECT_EQ(r_read.NumberOfElements(), 2);
        KRATOS_EXPECT_EQ(r_read.NumberOfConditions(), 1);
    }
    RemoveIfExists(file_path);
}

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusIODeformedConfiguration, KratosCoreFastSuite)
{
    Model model;
    auto& r_write_model_part = model.CreateModelPart("write");
    PopulateTetrahedraModelPart(r_write_model_part);
    for (auto& r_node : r_write_model_part.Nodes()) {
        r_node.X() = r_node.X0() + 0.5; // move the mesh
    }
    const auto file_path = TestFilePath(".vtu");

    { // default: the initial configuration (X0) is written
        Parameters settings(R"({"time_series" : "single_file"})");
        MeshioPlusPlusIO(file_path, settings).WriteModelPart(r_write_model_part);
        auto& r_read = model.CreateModelPart("read_initial");
        MeshioPlusPlusIO(file_path).ReadModelPart(r_read);
        auto it_write = r_write_model_part.NodesBegin();
        for (const auto& r_node : r_read.Nodes()) {
            KRATOS_EXPECT_NEAR(r_node.X(), it_write->X0(), 1e-12);
            ++it_write;
        }
    }
    { // deformed: the current coordinates are written
        Parameters settings(R"({"time_series" : "single_file", "write_deformed_configuration" : true})");
        MeshioPlusPlusIO(file_path, settings).WriteModelPart(r_write_model_part);
        auto& r_read = model.CreateModelPart("read_deformed");
        MeshioPlusPlusIO(file_path).ReadModelPart(r_read);
        auto it_write = r_write_model_part.NodesBegin();
        for (const auto& r_node : r_read.Nodes()) {
            KRATOS_EXPECT_NEAR(r_node.X(), it_write->X0() + 0.5, 1e-12);
            ++it_write;
        }
    }
    RemoveIfExists(file_path);
}

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusIOWriteIds, KratosCoreFastSuite)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("write");
    PopulateTetrahedraModelPart(r_model_part);

    const auto file_path = TestFilePath(".vtu");
    Parameters settings(R"({"time_series" : "single_file", "write_ids" : true})");
    MeshioPlusPlusIO(file_path, settings).WriteModelPart(r_model_part);

    const std::string content = ReadFileContent(file_path);
    KRATOS_EXPECT_TRUE(content.find("KRATOS_NODE_ID") != std::string::npos);
    KRATOS_EXPECT_TRUE(content.find("KRATOS_ELEMENT_ID") != std::string::npos);
    KRATOS_EXPECT_TRUE(content.find("KRATOS_CONDITION_ID") != std::string::npos);
    KRATOS_EXPECT_TRUE(content.find("PROPERTIES_ID") != std::string::npos);
    RemoveIfExists(file_path);
}

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusIOFlags, KratosCoreFastSuite)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("write");
    PopulateTetrahedraModelPart(r_model_part);
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.Set(TO_ERASE, true); // defined on every node -> no -1 values
    }
    r_model_part.ElementsBegin()->Set(ACTIVE, false);

    const auto file_path = TestFilePath(".xdmf");
    Parameters settings(R"({
        "xdmf_data_format" : "XML",
        "nodal_flags"      : ["TO_ERASE"],
        "element_flags"    : ["ACTIVE"],
        "condition_flags"  : ["TO_ERASE"]
    })");
    MeshioPlusPlusIO(file_path, settings).WriteModelPart(r_model_part);

    const std::string content = ReadFileContent(file_path);
    KRATOS_EXPECT_TRUE(content.find("Name=\"TO_ERASE\"") != std::string::npos);
    KRATOS_EXPECT_TRUE(content.find("Name=\"ACTIVE\"") != std::string::npos);
    KRATOS_EXPECT_EQ(CountOccurrences(content, "Center=\"Cell\""), 2); // ACTIVE + condition TO_ERASE
    RemoveIfExists(file_path);
}

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusIOCellData, KratosCoreFastSuite)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("write");
    PopulateTetrahedraModelPart(r_model_part);
    for (auto& r_element : r_model_part.Elements()) {
        r_element.SetValue(PRESSURE, 2.0 * r_element.Id());
    }
    for (auto& r_condition : r_model_part.Conditions()) {
        r_condition.SetValue(TEMPERATURE, 3.0 * r_condition.Id());
    }

    const auto file_path = TestFilePath(".vtu");
    Parameters settings(R"({
        "time_series"                    : "single_file",
        "element_data_value_variables"   : ["PRESSURE"],
        "condition_data_value_variables" : ["TEMPERATURE"]
    })");
    MeshioPlusPlusIO(file_path, settings).WriteModelPart(r_model_part);

    const std::string content = ReadFileContent(file_path);
    KRATOS_EXPECT_TRUE(content.find("PRESSURE") != std::string::npos);
    KRATOS_EXPECT_TRUE(content.find("TEMPERATURE") != std::string::npos);
    RemoveIfExists(file_path);
}

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusIOVariableTypes, KratosCoreFastSuite)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("write");
    PopulateTetrahedraModelPart(r_model_part);
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.SetValue(PARTITION_INDEX, static_cast<int>(r_node.Id()));      // int
        r_node.SetValue(BDF_COEFFICIENTS, Vector(4, 1.0));                    // Vector (4 components)
        r_node.SetValue(VELOCITY, array_1d<double, 3>{1.0, 2.0, 3.0});        // array_1d<double, 3>
    }

    const auto file_path = TestFilePath(".vtu");
    Parameters settings(R"({
        "time_series"                : "single_file",
        "nodal_data_value_variables" : ["PARTITION_INDEX", "BDF_COEFFICIENTS", "VELOCITY",
                                        "GREEN_LAGRANGE_STRAIN_TENSOR"]
    })");
    // The Matrix variable is unsupported: it must only warn, not throw
    MeshioPlusPlusIO(file_path, settings).WriteModelPart(r_model_part);

    const std::string content = ReadFileContent(file_path);
    KRATOS_EXPECT_TRUE(content.find("PARTITION_INDEX") != std::string::npos);
    KRATOS_EXPECT_TRUE(content.find("BDF_COEFFICIENTS") != std::string::npos);
    KRATOS_EXPECT_TRUE(content.find("VELOCITY") != std::string::npos);
    KRATOS_EXPECT_TRUE(content.find("GREEN_LAGRANGE_STRAIN_TENSOR") == std::string::npos);
    RemoveIfExists(file_path);
}

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusIOFileFormatOverride, KratosCoreFastSuite)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("write");
    PopulateTetrahedraModelPart(r_model_part);
    const auto file_path = TestFilePath(".vtu");

    { // ascii vtu
        Parameters settings(R"({"time_series" : "single_file", "file_format" : "ascii"})");
        MeshioPlusPlusIO(file_path, settings).WriteModelPart(r_model_part);
        const std::string content = ReadFileContent(file_path);
        KRATOS_EXPECT_TRUE(content.find("format=\"ascii\"") != std::string::npos);
    }
    { // binary vtu
        Parameters settings(R"({"time_series" : "single_file", "file_format" : "binary"})");
        MeshioPlusPlusIO(file_path, settings).WriteModelPart(r_model_part);
        const std::string content = ReadFileContent(file_path);
        KRATOS_EXPECT_TRUE(content.find("format=\"ascii\"") == std::string::npos);
    }
    RemoveIfExists(file_path);

    { // formats without an ascii/binary variant warn and use their default
        const auto su2_path = TestFilePath(".su2");
        Parameters settings(R"({"time_series" : "single_file", "file_format" : "ascii"})");
        MeshioPlusPlusIO(su2_path, settings).WriteModelPart(r_model_part);
        KRATOS_EXPECT_TRUE(std::filesystem::exists(su2_path));
        RemoveIfExists(su2_path);
    }
}

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusIONamePrefixPostfix, KratosCoreFastSuite)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("write");
    PopulateTetrahedraModelPart(r_model_part);
    r_model_part.GetProcessInfo()[STEP] = 7;

    const auto file_path = TestFilePath(".vtu");
    Parameters settings(R"({"custom_name_prefix" : "pre_", "custom_name_postfix" : "_post"})");
    MeshioPlusPlusIO(file_path, settings).WriteModelPart(r_model_part); // file series (automatic)

    auto expected_path = file_path.parent_path() /
        ("pre_" + file_path.stem().string() + "_post_7" + file_path.extension().string());
    KRATOS_EXPECT_TRUE(std::filesystem::exists(expected_path));
    RemoveIfExists(expected_path);
}

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusIOOutputSubModelParts, KratosCoreFastSuite)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("write");
    PopulateTetrahedraModelPart(r_model_part);
    auto& r_volume = r_model_part.CreateSubModelPart("Volume");
    r_volume.AddNodes({1, 2, 3, 4, 5});
    r_volume.AddElements({1, 2});
    auto& r_surface = r_model_part.CreateSubModelPart("Surface");
    r_surface.AddNodes({1, 2, 3});
    r_surface.AddConditions({1});
    r_model_part.GetProcessInfo()[STEP] = 1;

    { // file series: one file per target and step
        const auto file_path = TestFilePath(".vtu");
        Parameters settings(R"({"output_sub_model_parts" : true})");
        MeshioPlusPlusIO(file_path, settings).WriteModelPart(r_model_part);

        auto labeled = [&file_path](const std::string& rSuffix) {
            return file_path.parent_path() / (file_path.stem().string() + rSuffix + "_1.vtu");
        };
        KRATOS_EXPECT_TRUE(std::filesystem::exists(labeled("")));
        KRATOS_EXPECT_TRUE(std::filesystem::exists(labeled("_write_Volume")));
        KRATOS_EXPECT_TRUE(std::filesystem::exists(labeled("_write_Surface")));
        RemoveIfExists(labeled(""));
        RemoveIfExists(labeled("_write_Volume"));
        RemoveIfExists(labeled("_write_Surface"));
    }
    { // transient XDMF: one time-series file per target
        const auto file_path = TestFilePath(".xdmf");
        Parameters settings(R"({"output_sub_model_parts" : true, "xdmf_data_format" : "XML"})");
        MeshioPlusPlusIO io_write(file_path, settings);
        io_write.WriteModelPart(r_model_part);
        io_write.WriteModelPart(r_model_part);

        auto target = [&file_path](const std::string& rSuffix) {
            return file_path.parent_path() / (file_path.stem().string() + rSuffix + ".xdmf");
        };
        for (const std::string suffix : {"", "_write_Volume", "_write_Surface"}) {
            KRATOS_EXPECT_TRUE(std::filesystem::exists(target(suffix)));
            KRATOS_EXPECT_EQ(CountOccurrences(ReadFileContent(target(suffix)), "<Time Value="), 2);
            RemoveIfExists(target(suffix));
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusIOGaussPointVariables, KratosCoreFastSuite)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("write");
    PopulateTetrahedraModelPart(r_model_part);

    const auto file_path = TestFilePath(".vtu");
    // The generic core entities return no integration point results: the
    // in-elements variable must be skipped with a warning (not throw), and the
    // extrapolation process must run producing (zero-valued) nodal data.
    Parameters settings(R"({
        "time_series"                                 : "single_file",
        "gauss_point_variables_in_elements"           : ["CAUCHY_STRESS_VECTOR"],
        "gauss_point_variables_extrapolated_to_nodes" : ["TEMPERATURE"]
    })");
    MeshioPlusPlusIO(file_path, settings).WriteModelPart(r_model_part);

    const std::string content = ReadFileContent(file_path);
    KRATOS_EXPECT_TRUE(content.find("CAUCHY_STRESS_VECTOR") == std::string::npos); // skipped
    KRATOS_EXPECT_TRUE(content.find("TEMPERATURE") != std::string::npos);          // extrapolated nodal data
    RemoveIfExists(file_path);
}

} // namespace Kratos::Testing
