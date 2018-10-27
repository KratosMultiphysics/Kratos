//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//
//

// System includes


// External includes


// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "processes/mesh_coarsening_process.h"
#include "processes/structured_mesh_generator_process.h"
#include "processes/find_nodal_neighbours_process.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/hexahedra_3d_8.h"
#include "includes/gid_io.h"


namespace Kratos {
    namespace Testing {

        KRATOS_TEST_CASE_IN_SUITE(Structured2DMeshCoarseningProcess, KratosCoreFastSuite)
        {
            Model current_model;

            Node<3>::Pointer p_point1(new Node<3>(1, 0.00, 0.00, 0.00));
            Node<3>::Pointer p_point2(new Node<3>(2, 0.00, 10.00, 0.00));
            Node<3>::Pointer p_point3(new Node<3>(3, 10.00, 10.00, 0.00));
            Node<3>::Pointer p_point4(new Node<3>(4, 10.00, 0.00, 0.00));

            Quadrilateral2D4<Node<3> > geometry(p_point1, p_point2, p_point3, p_point4);

            ModelPart& model_part = current_model.CreateModelPart("Test");

            Parameters mesher_parameters(R"(
            {
                "number_of_divisions":1,
                "element_name": "Element2D3N"
            }  )");

            std::size_t number_of_divisions = mesher_parameters["number_of_divisions"].GetInt();

            StructuredMeshGeneratorProcess(geometry, model_part, mesher_parameters).Execute();

            // GidIO<> gid_io("c:/temp/coarsening/coarsening_test", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
            // gid_io.InitializeMesh(0.00);
            // gid_io.WriteMesh(model_part.GetMesh());
            // gid_io.FinalizeMesh();

            ModelPart& r_skin_model_part = model_part.GetSubModelPart("Skin");

            for (auto i_node = r_skin_model_part.NodesBegin(); i_node != r_skin_model_part.NodesEnd(); i_node++)
                i_node->Set(BOUNDARY);

            FindNodalNeighboursProcess(model_part).Execute();

            MeshCoarseningProcess(model_part).Execute();

            // gid_io.InitializeMesh(1.00);
            // gid_io.WriteMesh(model_part.GetMesh());
            // gid_io.FinalizeMesh();


            //KRATOS_WATCH(model_part.NumberOfNodes());
            //KRATOS_WATCH(model_part.NumberOfElements());
            //std::size_t number_of_nodes = (number_of_divisions + 1) * (number_of_divisions + 1);
            //KRATOS_CHECK_EQUAL(model_part.NumberOfNodes(), number_of_nodes);
            //KRATOS_CHECK_EQUAL(model_part.NumberOfElements(), number_of_divisions * number_of_divisions * 2);

            double total_area = 0.00;
            for (auto i_element = model_part.ElementsBegin(); i_element != model_part.ElementsEnd(); i_element++) {
                double element_area = i_element->GetGeometry().Area();
                KRATOS_CHECK_GREATER(element_area, 1.00/number_of_divisions) << " for element #" << i_element->Id() << " with nodes ["
                    << i_element->GetGeometry()[0].Id()
                    << "," << i_element->GetGeometry()[1].Id()
                    << "," << i_element->GetGeometry()[2].Id() << "] with area : " << element_area << std::endl << *i_element;
                total_area += element_area;
            }
            KRATOS_CHECK_NEAR(total_area, 100., 1e-6) << " : " << total_area << " != 100" << std::endl;
        }

        KRATOS_TEST_CASE_IN_SUITE(PerturbedStructured2DMeshCoarseningProcess, KratosCoreFastSuite)
        {
            Model current_model;

            Node<3>::Pointer p_point1(new Node<3>(1, 0.00, 0.00, 0.00));
            Node<3>::Pointer p_point2(new Node<3>(2, 0.00, 10.00, 0.00));
            Node<3>::Pointer p_point3(new Node<3>(3, 10.00, 10.00, 0.00));
            Node<3>::Pointer p_point4(new Node<3>(4, 10.00, 0.00, 0.00));

            Quadrilateral2D4<Node<3> > geometry(p_point1, p_point2, p_point3, p_point4);

            ModelPart& model_part = current_model.CreateModelPart("Test");

            Parameters mesher_parameters(R"(
            {
                "number_of_divisions":10,
                "element_name": "Element2D3N"
            }  )");

            std::size_t number_of_divisions = mesher_parameters["number_of_divisions"].GetInt();

            StructuredMeshGeneratorProcess(geometry, model_part, mesher_parameters).Execute();
            for (std::size_t i = 0; i < model_part.NumberOfNodes(); i++)
                model_part.GetNode(i + 1).Coordinates()[i % 3] += .1 / number_of_divisions;

            double original_mesh_area = 0.00;
            for (auto i_element = model_part.ElementsBegin(); i_element != model_part.ElementsEnd(); i_element++)
                original_mesh_area += i_element->GetGeometry().Area();

            // GidIO<> gid_io("c:/temp/coarsening/perturbed_coarsening_test", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
            // gid_io.InitializeMesh(0.00);
            // gid_io.WriteMesh(model_part.GetMesh());
            // gid_io.FinalizeMesh();

            ModelPart& r_skin_model_part = model_part.GetSubModelPart("Skin");

            for (auto i_node = r_skin_model_part.NodesBegin(); i_node != r_skin_model_part.NodesEnd(); i_node++)
                i_node->Set(BOUNDARY);

            FindNodalNeighboursProcess(model_part).Execute();

            MeshCoarseningProcess(model_part).Execute();

            // gid_io.InitializeMesh(1.00);
            // gid_io.WriteMesh(model_part.GetMesh());
            // gid_io.FinalizeMesh();


            //std::size_t number_of_nodes = (number_of_divisions + 1) * (number_of_divisions + 1);
            //KRATOS_CHECK_EQUAL(model_part.NumberOfNodes(), number_of_nodes);
            //KRATOS_CHECK_EQUAL(model_part.NumberOfElements(), number_of_divisions * number_of_divisions * 2);

            double total_area = 0.00;
            for (auto i_element = model_part.ElementsBegin(); i_element != model_part.ElementsEnd(); i_element++) {
                double element_area = i_element->GetGeometry().Area();
                KRATOS_CHECK_GREATER(element_area, 1.00 / number_of_divisions) << " for element #" << i_element->Id() << " with nodes ["
                    << i_element->GetGeometry()[0].Id()
                    << "," << i_element->GetGeometry()[1].Id()
                    << "," << i_element->GetGeometry()[2].Id() << "] with area : " << element_area << std::endl << *i_element;
                total_area += element_area;
            }
            KRATOS_CHECK_NEAR(total_area, original_mesh_area, 1e-6);
        }

        KRATOS_TEST_CASE_IN_SUITE(RedistributedStructured2DMeshCoarseningProcess, KratosCoreFastSuite)
        {
            Model current_model;

            Node<3>::Pointer p_point1(new Node<3>(1, 0.00, 0.00, 0.00));
            Node<3>::Pointer p_point2(new Node<3>(2, 0.00, 10.00, 0.00));
            Node<3>::Pointer p_point3(new Node<3>(3, 10.00, 10.00, 0.00));
            Node<3>::Pointer p_point4(new Node<3>(4, 10.00, 0.00, 0.00));

            Quadrilateral2D4<Node<3> > geometry(p_point1, p_point2, p_point3, p_point4);

            ModelPart& model_part = current_model.CreateModelPart("Test");

            Parameters mesher_parameters(R"(
            {
                "number_of_divisions":10,
                "element_name": "Element2D3N"
            }  )");

            std::size_t number_of_divisions = mesher_parameters["number_of_divisions"].GetInt();

            StructuredMeshGeneratorProcess(geometry, model_part, mesher_parameters).Execute();
            for (std::size_t i = 0; i < model_part.NumberOfNodes(); i++)
                model_part.GetNode(i + 1).Coordinates()[i % 2] += 4. / number_of_divisions;

            double original_mesh_area = 0.00;
            for (auto i_element = model_part.ElementsBegin(); i_element != model_part.ElementsEnd(); i_element++)
                original_mesh_area += i_element->GetGeometry().Area();

            // GidIO<> gid_io("c:/temp/coarsening/redistributed_coarsening_test", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
            // gid_io.InitializeMesh(0.00);
            // gid_io.WriteMesh(model_part.GetMesh());
            // gid_io.FinalizeMesh();

            ModelPart& r_skin_model_part = model_part.GetSubModelPart("Skin");

            for (auto i_node = r_skin_model_part.NodesBegin(); i_node != r_skin_model_part.NodesEnd(); i_node++)
                i_node->Set(BOUNDARY);

            FindNodalNeighboursProcess(model_part).Execute();

            MeshCoarseningProcess(model_part).Execute();

            FindNodalNeighboursProcess(model_part).Execute();

            MeshCoarseningProcess(model_part).Execute();


            // gid_io.InitializeMesh(1.00);
            // gid_io.WriteMesh(model_part.GetMesh());
            // gid_io.FinalizeMesh();


            //std::size_t number_of_nodes = (number_of_divisions + 1) * (number_of_divisions + 1);
            //KRATOS_CHECK_EQUAL(model_part.NumberOfNodes(), number_of_nodes);
            //KRATOS_CHECK_EQUAL(model_part.NumberOfElements(), number_of_divisions * number_of_divisions * 2);

            double total_area = 0.00;
            for (auto i_element = model_part.ElementsBegin(); i_element != model_part.ElementsEnd(); i_element++) {
                double element_area = i_element->GetGeometry().Area();
                KRATOS_CHECK_GREATER(element_area, 0.01 / number_of_divisions) << " for element #" << i_element->Id() << " with nodes ["
                    << i_element->GetGeometry()[0].Id()
                    << "," << i_element->GetGeometry()[1].Id()
                    << "," << i_element->GetGeometry()[2].Id() << "] with area : " << element_area << std::endl << *i_element;
                total_area += element_area;
            }
            KRATOS_CHECK_NEAR(total_area, original_mesh_area, 1e-6);
        }

        KRATOS_TEST_CASE_IN_SUITE(Structured3DMeshCoarseningProcess, KratosCoreFastSuite)
        {
            Model current_model;

            Node<3>::Pointer p_point1(new Node<3>(1, 0.00, 0.00, 0.00));
            Node<3>::Pointer p_point2(new Node<3>(2, 10.00, 0.00, 0.00));
            Node<3>::Pointer p_point3(new Node<3>(3, 10.00, 10.00, 0.00));
            Node<3>::Pointer p_point4(new Node<3>(4, 0.00, 10.00, 0.00));
            Node<3>::Pointer p_point5(new Node<3>(5, 0.00, 0.00, 10.00));
            Node<3>::Pointer p_point6(new Node<3>(6, 10.00, 0.00, 10.00));
            Node<3>::Pointer p_point7(new Node<3>(7, 10.00, 10.00, 10.00));
            Node<3>::Pointer p_point8(new Node<3>(8, 0.00, 10.00, 10.00));

            Hexahedra3D8<Node<3> > geometry(p_point1, p_point2, p_point3, p_point4, p_point5, p_point6, p_point7, p_point8);

            ModelPart& model_part = current_model.CreateModelPart("Test");

            Parameters mesher_parameters(R"(
            {
                "number_of_divisions":10,
                "element_name": "Element3D4N"
            }  )");

            //std::size_t number_of_divisions = mesher_parameters["number_of_divisions"].GetInt();

            StructuredMeshGeneratorProcess(geometry, model_part, mesher_parameters).Execute();

            // GidIO<> gid_io("c:/temp/coarsening/coarsening_3d_test", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
            // gid_io.InitializeMesh(0.00);
            // gid_io.WriteMesh(model_part.GetMesh());
            // gid_io.FinalizeMesh();

            ModelPart& r_skin_model_part = model_part.GetSubModelPart("Skin");

            for (auto i_node = r_skin_model_part.NodesBegin(); i_node != r_skin_model_part.NodesEnd(); i_node++)
                i_node->Set(BOUNDARY);

            FindNodalNeighboursProcess(model_part).Execute();

            MeshCoarseningProcess(model_part).Execute();

            // gid_io.InitializeMesh(1.00);
            // gid_io.WriteMesh(model_part.GetMesh());
            // gid_io.FinalizeMesh();


            //KRATOS_WATCH(model_part.NumberOfNodes());
            //KRATOS_WATCH(model_part.NumberOfElements());
            KRATOS_CHECK_NOT_EQUAL(model_part.NumberOfNodes(), 0);
            KRATOS_CHECK_NOT_EQUAL(model_part.NumberOfElements(), 0);


            const double cube_volume = 1000.;
            const double avarage_element_volume = cube_volume / model_part.NumberOfElements();
            double total_volume = 0.00;
            for (auto i_element = model_part.ElementsBegin(); i_element != model_part.ElementsEnd(); i_element++) {
                double element_volume = i_element->GetGeometry().Volume();
                KRATOS_CHECK_GREATER(element_volume, avarage_element_volume / 1000) << " for element #" << i_element->Id() << " with nodes ["
                    << i_element->GetGeometry()[0].Id()
                    << "," << i_element->GetGeometry()[1].Id()
                    << "," << i_element->GetGeometry()[2].Id() << "] with volume : " << element_volume << std::endl << *i_element;
                total_volume += element_volume;
            }
            KRATOS_CHECK_NEAR(total_volume, cube_volume, 1.E-6) << "with total_volume = " << total_volume;
        }

    }
}  // namespace Kratos.
