// KRATOS  __  __          _    _                _ _           _   _
//        |  \/  | ___  __| |  / \   _ __  _ __ | (_) ___ __ _| |_(_) ___  _ ___
//        | |\/| |/ _ \/ _` | / _ \ | '_ \| '_ \| | |/ __/ _` | __| |/ _ \| '_  |
//        | |  | |  __/ (_| |/ ___ \| |_) | |_) | | | (_| (_| | |_| | (_) | | | |
//        |_|  |_|\___|\__,_/_/   \_\ .__/| .__/|_|_|\___\__,_|\__|_|\___/|_| |_|
//                                  |_|   |_|
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/hexahedra_3d_8.h"
#include "processes/structured_mesh_generator_process.h"
#include "testing/testing.h"
#include "custom_io/med_model_part_io.h"
#include "custom_utilities/med_testing_utilities.h"

namespace Kratos::Testing {

namespace { // helpers namespace

void MedWriteReadModelPart(
    const std::filesystem::path& rFileName,
    const std::function<void(ModelPart&)>& rPopulateFunction)
{
    Model model;
    auto& test_model_part_write = model.CreateModelPart("test_write");
    auto& test_model_part_read = model.CreateModelPart("test_read");

    auto full_name = rFileName;
    full_name.replace_extension(".med");

    rPopulateFunction(test_model_part_write);

    { // encapsulating to ensure memory (aka file handle) is freed
        MedModelPartIO io_write(full_name, IO::WRITE);
        io_write.WriteModelPart(test_model_part_write);
    }
    { // encapsulating to ensure memory (aka file handle) is freed
        MedModelPartIO io_read(full_name);
        io_read.ReadModelPart(test_model_part_read);
    }

    // remove before checking ModelParts, as would be left over if comparison fails
    if (std::filesystem::exists(full_name)) {
        std::filesystem::remove(full_name);
    }

    MedTestingUtilities::CheckModelPartsAreEqual(test_model_part_write, test_model_part_read);
}

} // helpers namespace

KRATOS_TEST_CASE_IN_SUITE(WriteReadMedEmpty, KratosMedFastSuite)
{
    MedWriteReadModelPart(this->Name(), [](ModelPart& rModelPart){
        // deliberately do not create any entities
    });
}

KRATOS_TEST_CASE_IN_SUITE(WriteReadMedNodes, KratosMedFastSuite)
{
    MedWriteReadModelPart(this->Name(), [](ModelPart& rModelPart){
        int node_id = 0;
        for (int x=0; x<20; ++x) {
            for (int y=0; y<10; ++y) {
                for (int z=0; z<15; ++z) {
                    rModelPart.CreateNewNode(++node_id, x, y, z);
                }
            }
        }
    });
}

KRATOS_TEST_CASE_IN_SUITE(WriteReadMedNodesNonConsecutiveNodeIds, KratosMedFastSuite)
{
    MedWriteReadModelPart(this->Name(), [](ModelPart& rModelPart){
        int node_id = 1;
        for (int x=0; x<20; ++x) {
            for (int y=0; y<10; ++y) {
                for (int z=0; z<15; ++z) {
                    rModelPart.CreateNewNode(node_id, x, y, z);
                    node_id += 5;
                }
            }
        }
    });
}

KRATOS_TEST_CASE_IN_SUITE(WriteRead2DLineMesh, KratosMedFastSuite)
{
    MedWriteReadModelPart(this->Name(), [](ModelPart& rModelPart){
        constexpr std::size_t num_nodes = 10;
        for (std::size_t i=0; i<num_nodes; ++i) {
            rModelPart.CreateNewNode(i+1, i,0,0);
        }
        for (std::size_t i=0; i<num_nodes-1; ++i) {
            rModelPart.CreateNewGeometry("Line3D2", std::vector<ModelPart::IndexType>{i+1, i+2});
        }
    });
}

KRATOS_TEST_CASE_IN_SUITE(WriteRead2DLineMeshNonConsecutiveNodeIds, KratosMedFastSuite)
{
    MedWriteReadModelPart(this->Name(), [](ModelPart& rModelPart){
        constexpr std::size_t num_nodes = 10;
        for (std::size_t i=0; i<num_nodes; ++i) {
            rModelPart.CreateNewNode(i*10+1, i,0,0)->Id();
        }
        for (std::size_t i=0; i<num_nodes-1; ++i) {
            rModelPart.CreateNewGeometry("Line3D2", std::vector<ModelPart::IndexType>{i*10+1, (i+1)*10+1});
        }
    });
}

KRATOS_TEST_CASE_IN_SUITE(WriteRead2DTriangularMesh, KratosMedFastSuite)
{
    MedWriteReadModelPart(this->Name(), [](ModelPart& rModelPart){
        auto p_point_1 = Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0);
        auto p_point_2 = Kratos::make_intrusive<Node>(2, 0.0, 1.0, 0.0);
        auto p_point_3 = Kratos::make_intrusive<Node>(3, 1.0, 1.0, 0.0);
        auto p_point_4 = Kratos::make_intrusive<Node>(4, 1.0, 0.0, 0.0);

        Quadrilateral2D4<Node> geometry(p_point_1, p_point_2, p_point_3, p_point_4);

        Parameters mesher_parameters(R"(
        {
            "number_of_divisions": 7,
            "element_name": "Element3D3N",
            "create_skin_sub_model_part" : false,
            "create_body_sub_model_part" : false
        })");

        StructuredMeshGeneratorProcess(geometry, rModelPart, mesher_parameters).Execute();
        // create geometries!
        MedTestingUtilities::AddGeometriesFromElements(rModelPart);
    });
}

KRATOS_TEST_CASE_IN_SUITE(WriteRead2DQuadrilateralMesh, KratosMedFastSuite)
{
    MedWriteReadModelPart(this->Name(), [](ModelPart& rModelPart){
        constexpr std::size_t num_quads = 10;
        int node_id = 0;
        for (std::size_t i=0; i<num_quads+1; ++i) {
            rModelPart.CreateNewNode(++node_id, i*1,0,0);
            rModelPart.CreateNewNode(++node_id, i*1,0.5,0);
        }

        for (std::size_t i=0; i<num_quads; ++i) {
            rModelPart.CreateNewGeometry("Quadrilateral3D4", std::vector<ModelPart::IndexType>{
                (i*2)+1,
                (i*2)+3,
                (i*2)+4,
                (i*2)+2
            });
        }
    });
}

KRATOS_TEST_CASE_IN_SUITE(WriteRead2DTriAndQuadMesh, KratosMedFastSuite)
{
    MedWriteReadModelPart(this->Name(), [](ModelPart& rModelPart){
        constexpr std::size_t num_quads = 10;
        int node_id = 0;
        for (std::size_t i=0; i<num_quads+1; ++i) {
            rModelPart.CreateNewNode(++node_id, i*1,0,0);
            rModelPart.CreateNewNode(++node_id, i*1,0.5,0);
        }

        // first writing all quads
        for (std::size_t i=0; i<num_quads; ++i) {
            rModelPart.CreateNewGeometry("Quadrilateral3D4", std::vector<ModelPart::IndexType>{
                (i*2)+1,
                (i*2)+3,
                (i*2)+4,
                (i*2)+2
            });
        }

        // then writing all triangles
        for (std::size_t i=0; i<num_quads; ++i) {
            rModelPart.CreateNewNode(++node_id, i*1+0.5,1,0);
            rModelPart.CreateNewGeometry("Triangle3D3", std::vector<ModelPart::IndexType>{
                (i*2)+4,
                static_cast<ModelPart::IndexType>(node_id),
                (i*2)+2
            });
        }
    });
}

KRATOS_TEST_CASE_IN_SUITE(WriteRead3DTetraMesh, KratosMedFastSuite)
{
    MedWriteReadModelPart(this->Name(), [](ModelPart& rModelPart){
        const double max_x = 1.0;
        const double min_x = 0.0;
        const double max_y = 1.0;
        const double min_y = 0.0;
        const double max_z = 1.0;
        const double min_z = 0.0;
        auto p_point_1 = Kratos::make_intrusive<Node>(1, min_x, min_y, min_z);
        auto p_point_2 = Kratos::make_intrusive<Node>(2, max_x, min_y, min_z);
        auto p_point_3 = Kratos::make_intrusive<Node>(3, max_x, max_y, min_z);
        auto p_point_4 = Kratos::make_intrusive<Node>(4, min_x, max_y, min_z);
        auto p_point_5 = Kratos::make_intrusive<Node>(5, min_x, min_y, max_z);
        auto p_point_6 = Kratos::make_intrusive<Node>(6, max_x, min_y, max_z);
        auto p_point_7 = Kratos::make_intrusive<Node>(7, max_x, max_y, max_z);
        auto p_point_8 = Kratos::make_intrusive<Node>(8, min_x, max_y, max_z);
        Hexahedra3D8<Node> geometry(p_point_1, p_point_2, p_point_3, p_point_4, p_point_5, p_point_6, p_point_7, p_point_8);

        Parameters mesher_parameters(R"(
        {
            "number_of_divisions": 5,
            "element_name": "Element3D4N",
            "create_skin_sub_model_part" : false,
            "create_body_sub_model_part" : false
        })");

        StructuredMeshGeneratorProcess(geometry, rModelPart, mesher_parameters).Execute();
        // create geometries!
        MedTestingUtilities::AddGeometriesFromElements(rModelPart);
    });
}

KRATOS_TEST_CASE_IN_SUITE(WriteRead3DHexa, KratosMedFastSuite)
{
    MedWriteReadModelPart(this->Name(), [](ModelPart& rModelPart){
        const double max_x = 1.0;
        const double min_x = 0.0;
        const double max_y = 1.0;
        const double min_y = 0.0;
        const double max_z = 1.0;
        const double min_z = 0.0;
        auto p_point_1 = Kratos::make_intrusive<Node>(1, min_x, min_y, min_z);
        auto p_point_2 = Kratos::make_intrusive<Node>(2, max_x, min_y, min_z);
        auto p_point_3 = Kratos::make_intrusive<Node>(3, max_x, max_y, min_z);
        auto p_point_4 = Kratos::make_intrusive<Node>(4, min_x, max_y, min_z);
        auto p_point_5 = Kratos::make_intrusive<Node>(5, min_x, min_y, max_z);
        auto p_point_6 = Kratos::make_intrusive<Node>(6, max_x, min_y, max_z);
        auto p_point_7 = Kratos::make_intrusive<Node>(7, max_x, max_y, max_z);
        auto p_point_8 = Kratos::make_intrusive<Node>(8, min_x, max_y, max_z);
        auto p_geom = Kratos::make_shared<Hexahedra3D8<Node>>(p_point_1, p_point_2, p_point_3, p_point_4, p_point_5, p_point_6, p_point_7, p_point_8);

        rModelPart.AddGeometry(p_geom);
    });
}

} // Kratos::Testing
