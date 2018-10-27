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
#include "includes/kernel.h"
#include "processes/structured_mesh_generator_process.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/hexahedra_3d_8.h"


namespace Kratos {
    namespace Testing {

        KRATOS_TEST_CASE_IN_SUITE(StructuredMeshGeneratorProcessHexahedra, KratosCoreFastSuite)
        {
            Kernel kernel;

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

            ModelPart& model_part = current_model.CreateModelPart("Generated");

            Parameters mesher_parameters(R"(
            {
                "number_of_divisions":10,
                "element_name": "Element3D4N"
            }  )");

            std::size_t number_of_divisions = mesher_parameters["number_of_divisions"].GetInt();

            StructuredMeshGeneratorProcess(geometry, model_part, mesher_parameters).Execute();
            std::size_t number_of_nodes = (number_of_divisions + 1) * (number_of_divisions + 1) * (number_of_divisions + 1);
            std::size_t number_of_elements = number_of_divisions * number_of_divisions * number_of_divisions * 6;
            KRATOS_CHECK_EQUAL(model_part.NumberOfNodes(), number_of_nodes);
            KRATOS_CHECK_EQUAL(model_part.NumberOfElements(), number_of_elements) << " Number of elements = " << model_part.NumberOfElements() ;

            double total_volume = 0.00;
            for (auto i_element = model_part.ElementsBegin(); i_element != model_part.ElementsEnd(); i_element++) {
                double element_volume = i_element->GetGeometry().Volume();
                KRATOS_CHECK_GREATER(element_volume, 0.00) << " for element #" << i_element->Id() << " with nodes ["
                    << i_element->GetGeometry()[0].Id()
                    << "," << i_element->GetGeometry()[1].Id()
                    << "," << i_element->GetGeometry()[2].Id()
                    << "," << i_element->GetGeometry()[3].Id() << "] with volume : " << element_volume << std::endl << *i_element;
                total_volume += element_volume;
            }
            KRATOS_CHECK_NEAR(total_volume, 1000., 1.E-6) << "with total_volume = " << total_volume;

            KRATOS_CHECK(model_part.HasSubModelPart("Skin"));

            KRATOS_CHECK_EQUAL(model_part.GetSubModelPart("Skin").NumberOfNodes(), 602);
            KRATOS_CHECK_EQUAL(model_part.GetSubModelPart("Skin").NumberOfElements(), 0);
        }

        KRATOS_TEST_CASE_IN_SUITE(StructuredMeshGeneratorProcessQuadrilateral, KratosCoreFastSuite)
        {
            Kernel kernel;

            Model current_model;

            Node<3>::Pointer p_point1(new Node<3>(1, 0.00, 0.00, 0.00));
            Node<3>::Pointer p_point2(new Node<3>(2, 0.00, 10.00, 0.00));
            Node<3>::Pointer p_point3(new Node<3>(3, 10.00, 10.00, 0.00));
            Node<3>::Pointer p_point4(new Node<3>(4, 10.00, 0.00, 0.00));

            Quadrilateral2D4<Node<3> > geometry(p_point1, p_point2, p_point3, p_point4);

            ModelPart& model_part = current_model.CreateModelPart("Generated");

            Parameters mesher_parameters(R"(
            {
                "number_of_divisions":10,
                "element_name": "Element2D3N",
                "create_skin_sub_model_part": false
            }  )");

            std::size_t number_of_divisions = mesher_parameters["number_of_divisions"].GetInt();

            StructuredMeshGeneratorProcess(geometry, model_part, mesher_parameters).Execute();
            std::size_t number_of_nodes = (number_of_divisions + 1) * (number_of_divisions + 1);
            KRATOS_CHECK_EQUAL(model_part.NumberOfNodes(), number_of_nodes);
            KRATOS_CHECK_EQUAL(model_part.NumberOfElements(), number_of_divisions * number_of_divisions * 2);

            double total_area = 0.00;
            for (auto i_element = model_part.ElementsBegin(); i_element != model_part.ElementsEnd(); i_element++) {
                double element_area = i_element->GetGeometry().Area();
                KRATOS_CHECK_GREATER(element_area, 0.00) << " for element #" << i_element->Id() << " with nodes ["
                    << i_element->GetGeometry()[0].Id()
                    << "," << i_element->GetGeometry()[1].Id()
                    << "," << i_element->GetGeometry()[2].Id() << "] with area : " << element_area << std::endl << *i_element;
                total_area += element_area;
            }
            KRATOS_CHECK_NEAR(total_area, 100., 1.E-6) << "with total_area = " << total_area;

            KRATOS_CHECK_IS_FALSE(model_part.HasSubModelPart("Skin"));
        }
    }
}  // namespace Kratos.


