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
//

// System includes
#include <fstream>
#include <sstream>
#include <filesystem>

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "includes/variables.h"
#include "input_output/unv_output.h"

namespace Kratos::Testing {

namespace {

/// A single element record parsed from a UNV 2412 dataset.
struct UnvElementRecord {
    int fe_id = 0;
    std::vector<int> connectivity;
};

/// A parsed UNV dataset block: its id and the whitespace tokens between the delimiters.
struct UnvBlock {
    int id = 0;
    std::vector<std::string> tokens;
};

/// Light UNV parser: splits the file into "-1 <id> ... -1" delimited blocks.
std::vector<UnvBlock> ParseUnvBlocks(const std::string& rFileName)
{
    std::ifstream input_file(rFileName);
    std::vector<std::string> all_tokens;
    std::string token;
    while (input_file >> token) {
        all_tokens.push_back(token);
    }

    std::vector<UnvBlock> blocks;
    std::size_t i = 0;
    while (i < all_tokens.size()) {
        if (all_tokens[i] == "-1") {
            ++i;
            if (i >= all_tokens.size() || all_tokens[i] == "-1") {
                continue;
            }
            UnvBlock block;
            block.id = std::stoi(all_tokens[i]);
            ++i;
            while (i < all_tokens.size() && all_tokens[i] != "-1") {
                block.tokens.push_back(all_tokens[i]);
                ++i;
            }
            ++i; // skip the closing "-1"
            blocks.push_back(block);
        } else {
            ++i;
        }
    }
    return blocks;
}

/// Returns the first block with the given dataset id (fails the test if not found is caller's job).
const UnvBlock* GetBlock(const std::vector<UnvBlock>& rBlocks, const int Id)
{
    for (const auto& r_block : rBlocks) {
        if (r_block.id == Id) {
            return &r_block;
        }
    }
    return nullptr;
}

/// Parses the element records of a UNV 2412 dataset block.
std::vector<UnvElementRecord> ParseElements(const UnvBlock& rBlock)
{
    std::vector<UnvElementRecord> records;
    const auto& t = rBlock.tokens;
    std::size_t i = 0;
    while (i + 6 <= t.size()) {
        const int fe_id = std::stoi(t[i + 1]);
        const int number_of_nodes = std::stoi(t[i + 5]);
        i += 6;
        const bool is_beam = (fe_id == 11 || fe_id == 21 || fe_id == 22 || fe_id == 23 || fe_id == 24);
        if (is_beam) {
            i += 3; // skip the beam orientation record
        }
        UnvElementRecord record;
        record.fe_id = fe_id;
        for (int k = 0; k < number_of_nodes; ++k) {
            record.connectivity.push_back(std::stoi(t[i]));
            ++i;
        }
        records.push_back(record);
    }
    return records;
}

/// Writes the mesh of the "Main" model part to a UNV file and returns the parsed blocks.
std::vector<UnvBlock> WriteAndParse(Model& rModel, const bool Decompose)
{
    Parameters parameters(R"({
        "model_part_name"                 : "Main",
        "save_output_files_in_folder"     : false,
        "decompose_quadratic_into_linear" : false
    })");
    parameters["decompose_quadratic_into_linear"].SetBool(Decompose);

    const std::string file_name = "Main.unv";
    if (std::filesystem::exists(file_name)) {
        std::filesystem::remove(file_name);
    }
    UnvOutput unv_output(rModel, parameters);
    unv_output.InitializeOutputFile();
    unv_output.WriteMesh();
    return ParseUnvBlocks(file_name);
}

} // anonymous namespace

KRATOS_TEST_CASE_IN_SUITE(UnvOutputLinearTriangle, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");
    auto p_prop = r_model_part.CreateNewProperties(1);
    r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    r_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    r_model_part.CreateNewElement("Element2D3N", 1, std::vector<ModelPart::IndexType>{1, 2, 3}, p_prop);

    const auto blocks = WriteAndParse(current_model, false);

    const UnvBlock* p_nodes = GetBlock(blocks, 2411);
    KRATOS_EXPECT_NE(p_nodes, nullptr);
    KRATOS_EXPECT_EQ(p_nodes->tokens.size() / 7, 3u); // 4 header + 3 coord tokens per node

    const UnvBlock* p_elements = GetBlock(blocks, 2412);
    KRATOS_EXPECT_NE(p_elements, nullptr);
    const auto elements = ParseElements(*p_elements);
    KRATOS_EXPECT_EQ(elements.size(), 1u);
    KRATOS_EXPECT_EQ(elements[0].fe_id, 41); // Plane Stress Linear Triangle
    const std::vector<int> expected_connectivity{1, 2, 3};
    KRATOS_EXPECT_EQ(elements[0].connectivity, expected_connectivity);
}

KRATOS_TEST_CASE_IN_SUITE(UnvOutputQuadraticTetrahedra, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");
    auto p_prop = r_model_part.CreateNewProperties(1);
    for (std::size_t i = 1; i <= 10; ++i) {
        r_model_part.CreateNewNode(i, static_cast<double>(i), 0.0, 0.0);
    }
    r_model_part.CreateNewElement("Element3D10N", 1, std::vector<ModelPart::IndexType>{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, p_prop);

    const auto blocks = WriteAndParse(current_model, false);
    const UnvBlock* p_elements = GetBlock(blocks, 2412);
    KRATOS_EXPECT_NE(p_elements, nullptr);
    const auto elements = ParseElements(*p_elements);
    KRATOS_EXPECT_EQ(elements.size(), 1u);
    KRATOS_EXPECT_EQ(elements[0].fe_id, 118); // Solid Parabolic Tetrahedron
    // UNV order from Kratos indices [0,4,1,5,2,6,7,8,9,3] -> node ids (1-based)
    const std::vector<int> expected_connectivity{1, 5, 2, 6, 3, 7, 8, 9, 10, 4};
    KRATOS_EXPECT_EQ(elements[0].connectivity, expected_connectivity);
}

KRATOS_TEST_CASE_IN_SUITE(UnvOutputQuadraticHexahedra, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");
    auto p_prop = r_model_part.CreateNewProperties(1);
    std::vector<ModelPart::IndexType> connectivity;
    for (std::size_t i = 1; i <= 20; ++i) {
        r_model_part.CreateNewNode(i, static_cast<double>(i), 0.0, 0.0);
        connectivity.push_back(i);
    }
    r_model_part.CreateNewElement("Element3D20N", 1, connectivity, p_prop);

    const auto blocks = WriteAndParse(current_model, false);
    const UnvBlock* p_elements = GetBlock(blocks, 2412);
    KRATOS_EXPECT_NE(p_elements, nullptr);
    const auto elements = ParseElements(*p_elements);
    KRATOS_EXPECT_EQ(elements.size(), 1u);
    KRATOS_EXPECT_EQ(elements[0].fe_id, 116); // Solid Parabolic Brick
    // UNV order from Kratos indices [0,8,1,9,2,10,3,11,12,13,14,15,4,16,5,17,6,18,7,19]
    const std::vector<int> expected_connectivity{1, 9, 2, 10, 3, 11, 4, 12, 13, 14, 15, 16, 5, 17, 6, 18, 7, 19, 8, 20};
    KRATOS_EXPECT_EQ(elements[0].connectivity, expected_connectivity);
}

KRATOS_TEST_CASE_IN_SUITE(UnvOutputBeamRecord, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");
    auto p_prop = r_model_part.CreateNewProperties(1);
    r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    r_model_part.CreateNewElement("Element2D2N", 1, std::vector<ModelPart::IndexType>{1, 2}, p_prop); // Line2D2

    const auto blocks = WriteAndParse(current_model, false);
    const UnvBlock* p_elements = GetBlock(blocks, 2412);
    KRATOS_EXPECT_NE(p_elements, nullptr);
    // The parser skips the orientation record for beams; connectivity must be the two nodes.
    const auto elements = ParseElements(*p_elements);
    KRATOS_EXPECT_EQ(elements.size(), 1u);
    KRATOS_EXPECT_EQ(elements[0].fe_id, 21); // Linear beam
    const std::vector<int> expected_connectivity{1, 2};
    KRATOS_EXPECT_EQ(elements[0].connectivity, expected_connectivity);
    // Header(6) + orientation(3) + connectivity(2) = 11 tokens
    KRATOS_EXPECT_EQ(p_elements->tokens.size(), 11u);
}

KRATOS_TEST_CASE_IN_SUITE(UnvOutputDegradeQuadratic9Quad, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");
    auto p_prop = r_model_part.CreateNewProperties(1);
    for (std::size_t i = 1; i <= 9; ++i) {
        r_model_part.CreateNewNode(i, static_cast<double>(i), 0.0, 0.0);
    }
    r_model_part.CreateNewElement("Element2D9N", 1, std::vector<ModelPart::IndexType>{1, 2, 3, 4, 5, 6, 7, 8, 9}, p_prop);

    const auto blocks = WriteAndParse(current_model, false);
    const auto elements = ParseElements(*GetBlock(blocks, 2412));
    KRATOS_EXPECT_EQ(elements.size(), 1u);
    KRATOS_EXPECT_EQ(elements[0].fe_id, 45); // Plane Stress Parabolic Quadrilateral (degraded)
    // Center node (Kratos index 8, id 9) is dropped.
    const std::vector<int> expected_connectivity{1, 5, 2, 6, 3, 7, 4, 8};
    KRATOS_EXPECT_EQ(elements[0].connectivity, expected_connectivity);
}

KRATOS_TEST_CASE_IN_SUITE(UnvOutputDegradeHexahedra27, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");
    auto p_prop = r_model_part.CreateNewProperties(1);
    std::vector<ModelPart::IndexType> connectivity;
    for (std::size_t i = 1; i <= 27; ++i) {
        r_model_part.CreateNewNode(i, static_cast<double>(i), 0.0, 0.0);
        connectivity.push_back(i);
    }
    r_model_part.CreateNewElement("Element3D27N", 1, connectivity, p_prop);

    const auto blocks = WriteAndParse(current_model, false);
    const auto elements = ParseElements(*GetBlock(blocks, 2412));
    KRATOS_EXPECT_EQ(elements.size(), 1u);
    KRATOS_EXPECT_EQ(elements[0].fe_id, 116); // degraded to 20-node parabolic brick
    KRATOS_EXPECT_EQ(elements[0].connectivity.size(), 20u);
}

KRATOS_TEST_CASE_IN_SUITE(UnvOutputUnsupportedPyramidThrows, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");
    auto p_prop = r_model_part.CreateNewProperties(1);
    for (std::size_t i = 1; i <= 5; ++i) {
        r_model_part.CreateNewNode(i, static_cast<double>(i), 0.0, 0.0);
    }
    r_model_part.CreateNewElement("Element3D5N", 1, std::vector<ModelPart::IndexType>{1, 2, 3, 4, 5}, p_prop); // Pyramid3D5

    UnvOutput unv_output(r_model_part, "test_unv_pyramid");
    unv_output.InitializeOutputFile();
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        unv_output.WriteMesh(),
        "cannot be represented in UNV");
}

KRATOS_TEST_CASE_IN_SUITE(UnvOutputNodalAndElementResults, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");
    r_model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    r_model_part.SetBufferSize(1);
    auto p_prop = r_model_part.CreateNewProperties(1);
    r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    r_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    auto p_element = r_model_part.CreateNewElement("Element2D3N", 1, std::vector<ModelPart::IndexType>{1, 2, 3}, p_prop);
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.FastGetSolutionStepValue(TEMPERATURE) = 20.0;
    }
    p_element->SetValue(PRESSURE, 5.0);
    r_model_part.GetProcessInfo()[TIME] = 1.0;

    // The (Model, Parameters) constructor derives the file name from the model part name.
    const std::string name = "Main";
    if (std::filesystem::exists(name + ".unv")) {
        std::filesystem::remove(name + ".unv");
    }
    Parameters parameters(R"({
        "model_part_name"                    : "Main",
        "save_output_files_in_folder"        : false,
        "nodal_solution_step_data_variables" : ["TEMPERATURE"],
        "element_data_value_variables"       : ["PRESSURE"]
    })");
    UnvOutput unv_output(current_model, parameters);
    unv_output.InitializeOutputFile();
    unv_output.WriteMesh();
    unv_output.PrintOutput();

    const auto blocks = ParseUnvBlocks(name + ".unv");
    // Two result datasets (2414): one nodal (location 1) and one element (location 2)
    int nodal_results = 0, element_results = 0;
    for (const auto& r_block : blocks) {
        if (r_block.id == 2414 && r_block.tokens.size() >= 3) {
            if (r_block.tokens[2] == "1") ++nodal_results;
            if (r_block.tokens[2] == "2") ++element_results;
        }
    }
    KRATOS_EXPECT_EQ(nodal_results, 1);
    KRATOS_EXPECT_EQ(element_results, 1);
}

KRATOS_TEST_CASE_IN_SUITE(UnvOutputSubModelPartGroups, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");
    auto p_prop = r_model_part.CreateNewProperties(1);
    r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    r_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    r_model_part.CreateNewElement("Element2D3N", 1, std::vector<ModelPart::IndexType>{1, 2, 3}, p_prop);
    auto& r_sub = r_model_part.CreateSubModelPart("Sub");
    r_sub.AddNodes({1, 2, 3});
    r_sub.AddElement(r_model_part.pGetElement(1));
    r_model_part.GetProcessInfo()[TIME] = 1.0;

    // The (Model, Parameters) constructor derives the file name from the model part name.
    const std::string name = "Main";
    if (std::filesystem::exists(name + ".unv")) {
        std::filesystem::remove(name + ".unv");
    }
    Parameters parameters(R"({
        "model_part_name"                 : "Main",
        "save_output_files_in_folder"     : false,
        "output_sub_model_parts"          : true,
        "decompose_quadratic_into_linear" : false
    })");
    UnvOutput unv_output(current_model, parameters);
    unv_output.InitializeOutputFile();
    unv_output.WriteMesh();

    const auto blocks = ParseUnvBlocks(name + ".unv");
    const UnvBlock* p_groups = GetBlock(blocks, 2467);
    KRATOS_EXPECT_NE(p_groups, nullptr);
    // Record 1 header: group number, 6 zeros, number of entities (1 element)
    KRATOS_EXPECT_EQ(std::stoi(p_groups->tokens[7]), 1);
    // Entity tuple: type code 8 (finite element), tag 1
    KRATOS_EXPECT_EQ(std::stoi(p_groups->tokens[9]), 8);
    KRATOS_EXPECT_EQ(std::stoi(p_groups->tokens[10]), 1);
}

KRATOS_TEST_CASE_IN_SUITE(UnvOutputFlagsIdsAndGaussPoints, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");
    auto p_prop = r_model_part.CreateNewProperties(1);
    r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    r_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    r_model_part.CreateNewElement("Element2D3N", 1, std::vector<ModelPart::IndexType>{1, 2, 3}, p_prop);
    r_model_part.GetNode(1).Set(ACTIVE, true);
    r_model_part.GetProcessInfo()[TIME] = 1.0;

    const std::string name = "Main";
    if (std::filesystem::exists(name + ".unv")) {
        std::filesystem::remove(name + ".unv");
    }
    Parameters parameters(R"({
        "model_part_name"                   : "Main",
        "save_output_files_in_folder"       : false,
        "write_ids"                         : true,
        "nodal_flags"                       : ["ACTIVE"],
        "gauss_point_variables_in_elements" : ["TEMPERATURE"]
    })");
    UnvOutput unv_output(current_model, parameters);
    unv_output.InitializeOutputFile();
    unv_output.WriteMesh();
    // Must not throw at runtime, even for a generic element with no integration point data.
    unv_output.PrintOutput();

    const auto blocks = ParseUnvBlocks(name + ".unv");
    // Datasets expected: nodes, elements and several results (flag, gauss, and 3 id datasets)
    KRATOS_EXPECT_NE(GetBlock(blocks, 2411), nullptr);
    KRATOS_EXPECT_NE(GetBlock(blocks, 2412), nullptr);
    int result_datasets = 0;
    for (const auto& r_block : blocks) {
        if (r_block.id == 2414) ++result_datasets;
    }
    // 1 nodal flag + 1 gauss-in-elements + 3 id datasets (node/element/condition)
    KRATOS_EXPECT_EQ(result_datasets, 5);
}

KRATOS_TEST_CASE_IN_SUITE(UnvOutputDecomposeQuadraticTetrahedra, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");
    auto p_prop = r_model_part.CreateNewProperties(1);
    for (std::size_t i = 1; i <= 10; ++i) {
        r_model_part.CreateNewNode(i, static_cast<double>(i), 0.0, 0.0);
    }
    r_model_part.CreateNewElement("Element3D10N", 1, std::vector<ModelPart::IndexType>{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, p_prop);

    const auto blocks = WriteAndParse(current_model, true);
    const auto elements = ParseElements(*GetBlock(blocks, 2412));
    // Tet10 decomposes into 8 linear tetrahedra (FE 111)
    KRATOS_EXPECT_EQ(elements.size(), 8u);
    for (const auto& r_element : elements) {
        KRATOS_EXPECT_EQ(r_element.fe_id, 111);
        KRATOS_EXPECT_EQ(r_element.connectivity.size(), 4u);
    }
    // First sub-tet uses Kratos indices [0,4,6,7] -> node ids 1,5,7,8
    const std::vector<int> expected_first{1, 5, 7, 8};
    KRATOS_EXPECT_EQ(elements[0].connectivity, expected_first);
}

KRATOS_TEST_CASE_IN_SUITE(UnvOutputDecomposeQuadraticTriangle, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");
    auto p_prop = r_model_part.CreateNewProperties(1);
    for (std::size_t i = 1; i <= 6; ++i) {
        r_model_part.CreateNewNode(i, static_cast<double>(i), 0.0, 0.0);
    }
    r_model_part.CreateNewElement("Element2D6N", 1, std::vector<ModelPart::IndexType>{1, 2, 3, 4, 5, 6}, p_prop);

    const auto blocks = WriteAndParse(current_model, true);
    const auto elements = ParseElements(*GetBlock(blocks, 2412));
    // Tri6 decomposes into 4 linear triangles (FE 41)
    KRATOS_EXPECT_EQ(elements.size(), 4u);
    for (const auto& r_element : elements) {
        KRATOS_EXPECT_EQ(r_element.fe_id, 41);
        KRATOS_EXPECT_EQ(r_element.connectivity.size(), 3u);
    }
    // First sub-triangle uses Kratos indices [0,3,5] -> node ids 1,4,6
    const std::vector<int> expected_first{1, 4, 6};
    KRATOS_EXPECT_EQ(elements[0].connectivity, expected_first);
}

KRATOS_TEST_CASE_IN_SUITE(UnvOutputDecomposeIsDefaultOff, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");
    auto p_prop = r_model_part.CreateNewProperties(1);
    for (std::size_t i = 1; i <= 10; ++i) {
        r_model_part.CreateNewNode(i, static_cast<double>(i), 0.0, 0.0);
    }
    r_model_part.CreateNewElement("Element3D10N", 1, std::vector<ModelPart::IndexType>{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, p_prop);

    // The lightweight (ModelPart, name) constructor uses the default, which keeps the quadratic element.
    const std::string file_name = "test_unv_default.unv";
    if (std::filesystem::exists(file_name)) {
        std::filesystem::remove(file_name);
    }
    UnvOutput unv_output(r_model_part, "test_unv_default");
    unv_output.InitializeOutputFile();
    unv_output.WriteMesh();

    const auto blocks = ParseUnvBlocks(file_name);
    const auto elements = ParseElements(*GetBlock(blocks, 2412));
    KRATOS_EXPECT_EQ(elements.size(), 1u); // not decomposed by default
    KRATOS_EXPECT_EQ(elements[0].fe_id, 118); // Solid Parabolic Tetrahedron
}

KRATOS_TEST_CASE_IN_SUITE(UnvOutputDeformationFactorScalesDisplacement, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");
    r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_model_part.SetBufferSize(1);
    auto p_prop = r_model_part.CreateNewProperties(1);
    r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    r_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    r_model_part.CreateNewElement("Element2D3N", 1, std::vector<ModelPart::IndexType>{1, 2, 3}, p_prop);
    array_1d<double, 3> displacement;
    displacement[0] = 1000.0;
    displacement[1] = 1000.0;
    displacement[2] = 1000.0;
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.FastGetSolutionStepValue(DISPLACEMENT) = displacement;
    }
    r_model_part.GetProcessInfo()[TIME] = 1.0;

    const std::string name = "Main";
    if (std::filesystem::exists(name + ".unv")) {
        std::filesystem::remove(name + ".unv");
    }
    Parameters parameters(R"({
        "model_part_name"                    : "Main",
        "save_output_files_in_folder"        : false,
        "deformation_factor"                 : 0.001,
        "nodal_solution_step_data_variables" : ["DISPLACEMENT"]
    })");
    UnvOutput unv_output(current_model, parameters);
    unv_output.InitializeOutputFile();
    unv_output.WriteMesh();
    unv_output.PrintOutput();

    const auto blocks = ParseUnvBlocks(name + ".unv");
    // Find the DISPLACEMENT results dataset (2414, location "1") and read the first node's values.
    const UnvBlock* p_displacement = nullptr;
    for (const auto& r_block : blocks) {
        if (r_block.id == 2414 && !r_block.tokens.empty() && r_block.tokens[0] == "DISPLACEMENT") {
            p_displacement = &r_block;
            break;
        }
    }
    KRATOS_EXPECT_NE(p_displacement, nullptr);
    // The dataset header occupies tokens [0..30]; token 31 is node 1's id and tokens 32-34 its values.
    // 1000.0 * 0.001 (deformation_factor) = 1.0 for each component.
    const std::size_t first_value_index = 32;
    KRATOS_EXPECT_NEAR(std::stod(p_displacement->tokens[first_value_index]), 1.0, 1e-9);
    KRATOS_EXPECT_NEAR(std::stod(p_displacement->tokens[first_value_index + 1]), 1.0, 1e-9);
    KRATOS_EXPECT_NEAR(std::stod(p_displacement->tokens[first_value_index + 2]), 1.0, 1e-9);
}

} // namespace Kratos::Testing
