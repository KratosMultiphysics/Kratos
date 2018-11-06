//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:   BSD License
//      Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

// System includes

// External includes

// Project includes
#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "testing/testing.h"
#include "includes/kratos_flags.h"
#include "containers/model.h"
#include "meshing_application.h"

// Utilities
#include "custom_utilities/uniform_refinement_utility.h"

namespace Kratos
{
    namespace Testing
    {
        typedef Node<3> NodeType;

        double Tolerance = 1e-6;

        double DistanceFunction(ModelPart::NodesContainerType::iterator pNode)
        {
            return 1.0 * pNode->X() + 2.0 * pNode->Y();
        }

        /**
         * Checks the correct refining utility with triangles
         */
        KRATOS_TEST_CASE_IN_SUITE(UniformRefineTrianglesUtility, KratosMeshingApplicationFastSuite)
        {
            Model this_model;
            ModelPart& this_model_part = this_model.CreateModelPart("Main", 2);

            this_model_part.AddNodalSolutionStepVariable(DISTANCE);

            Properties::Pointer p_properties = this_model_part.pGetProperties(0);

            // Creating the sub model parts
            ModelPart& r_sub_model_part_1 = this_model_part.CreateSubModelPart("BodySubModelPart");
            ModelPart& r_sub_model_part_2 = this_model_part.CreateSubModelPart("SkinSubModelPart");

            // Creating the nodes
            NodeType::Pointer p_node_1 = this_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
            NodeType::Pointer p_node_2 = this_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
            NodeType::Pointer p_node_3 = this_model_part.CreateNewNode(3, 2.0, 0.0, 0.0);
            NodeType::Pointer p_node_4 = this_model_part.CreateNewNode(4, 0.0, 1.0, 0.0);
            NodeType::Pointer p_node_5 = this_model_part.CreateNewNode(5, 1.0, 1.0, 0.0);
            NodeType::Pointer p_node_6 = this_model_part.CreateNewNode(6, 2.0, 1.0, 0.0);

            // Creating the elements
            std::vector<NodeType::Pointer> triangle_nodes(3);
            triangle_nodes[0] = p_node_1;
            triangle_nodes[1] = p_node_2;
            triangle_nodes[2] = p_node_4;
            Element::Pointer p_elem_1 = this_model_part.CreateNewElement("Element2D3N", 1, PointerVector<NodeType>{triangle_nodes}, p_properties);

            triangle_nodes[0] = p_node_2;
            triangle_nodes[1] = p_node_5;
            triangle_nodes[2] = p_node_4;
            Element::Pointer p_elem_2 = this_model_part.CreateNewElement("Element2D3N", 2, PointerVector<NodeType>{triangle_nodes}, p_properties);

            triangle_nodes[0] = p_node_2;
            triangle_nodes[1] = p_node_6;
            triangle_nodes[2] = p_node_5;
            Element::Pointer p_elem_3 = this_model_part.CreateNewElement("Element2D3N", 3, PointerVector<NodeType>{triangle_nodes}, p_properties);

            triangle_nodes[0] = p_node_2;
            triangle_nodes[1] = p_node_3;
            triangle_nodes[2] = p_node_6;
            Element::Pointer p_elem_4 = this_model_part.CreateNewElement("Element2D3N", 4, PointerVector<NodeType>{triangle_nodes}, p_properties);

            // Adding elements and its nodes to the sub model part
            r_sub_model_part_1.AddNode(p_node_1);
            r_sub_model_part_1.AddNode(p_node_2);
            r_sub_model_part_1.AddNode(p_node_3);
            r_sub_model_part_1.AddNode(p_node_4);
            r_sub_model_part_1.AddNode(p_node_5);
            r_sub_model_part_1.AddNode(p_node_6);
            r_sub_model_part_1.AddElement(p_elem_1);
            r_sub_model_part_1.AddElement(p_elem_2);
            r_sub_model_part_1.AddElement(p_elem_3);
            r_sub_model_part_1.AddElement(p_elem_4);

            // Creating the conditions
            std::vector<NodeType::Pointer> line_nodes(2);
            line_nodes[0] = p_node_1;
            line_nodes[1] = p_node_4;
            Condition::Pointer p_cond_1 = this_model_part.CreateNewCondition("Condition2D2N", 1, PointerVector<NodeType>{line_nodes}, p_properties);

            // Adding conditions and its nodes to the sub model part
            r_sub_model_part_2.AddNode(p_node_1);
            r_sub_model_part_2.AddNode(p_node_2);
            r_sub_model_part_2.AddCondition(p_cond_1);

            // Set the variables
            for (std::size_t i_node = 0; i_node < this_model_part.Nodes().size(); ++i_node) {
                auto it_node = this_model_part.Nodes().begin() + i_node;
                double& distance = it_node->FastGetSolutionStepValue(DISTANCE);
                distance = DistanceFunction(it_node);
            }

            // Store the initial values
            std::vector<unsigned int> initial_nodes(3);
            std::vector<unsigned int> initial_elems(3);
            std::vector<unsigned int> initial_conds(3);
            initial_nodes[0] = this_model_part.NumberOfNodes();
            initial_elems[0] = this_model_part.NumberOfElements();
            initial_conds[0] = this_model_part.NumberOfConditions();
            initial_nodes[1] = r_sub_model_part_1.NumberOfNodes();
            initial_elems[1] = r_sub_model_part_1.NumberOfElements();
            initial_conds[1] = r_sub_model_part_1.NumberOfConditions();
            initial_nodes[2] = r_sub_model_part_2.NumberOfNodes();
            initial_nodes[2] = r_sub_model_part_2.NumberOfElements();
            initial_conds[2] = r_sub_model_part_2.NumberOfConditions();

            // Execute the utility
            int refinement_level = 2;
            UniformRefinementUtility uniform_refinement(this_model_part);
            uniform_refinement.Refine(refinement_level);

            // Check the number of entities in the main model part
            unsigned int final_nodes = (std::pow(2,refinement_level)+1) * ((0.5*initial_nodes[0]-1)*std::pow(2,refinement_level)+1);
            KRATOS_CHECK_EQUAL(final_nodes, this_model_part.NumberOfNodes());

            unsigned int final_elems = initial_elems[0] * std::pow(4,refinement_level);
            KRATOS_CHECK_EQUAL(final_elems, this_model_part.NumberOfElements());

            unsigned int final_conds = initial_conds[0] * std::pow(2,refinement_level);
            KRATOS_CHECK_EQUAL(final_conds, this_model_part.NumberOfConditions());

            // Check the number of entities in the first sub model part
            // final_nodes = (std::pow(2,refinement_level)+1) * ((0.5*initial_nodes[1]-1)*std::pow(2,refinement_level)+1);
            // KRATOS_CHECK_EQUAL(final_nodes, r_sub_model_part_1.NumberOfNodes());

            final_elems = initial_elems[1] * std::pow(4,refinement_level);
            KRATOS_CHECK_EQUAL(final_elems, r_sub_model_part_1.NumberOfElements());

            final_conds = initial_conds[1] * std::pow(2,refinement_level);
            KRATOS_CHECK_EQUAL(final_conds, r_sub_model_part_1.NumberOfConditions());

            // Check the number of entities in the second sub model part
            // final_nodes = (std::pow(2,refinement_level)+1);
            // KRATOS_CHECK_EQUAL(final_nodes, r_sub_model_part_2.NumberOfNodes());

            final_elems = initial_elems[2] * std::pow(4,refinement_level);
            KRATOS_CHECK_EQUAL(final_elems, r_sub_model_part_2.NumberOfElements());

            final_conds = initial_conds[2] * std::pow(2,refinement_level);
            KRATOS_CHECK_EQUAL(final_conds, r_sub_model_part_2.NumberOfConditions());

            // Check the variables interpolation
            for (std::size_t i_node = 0; i_node < this_model_part.Nodes().size(); ++i_node) {
                auto it_node = this_model_part.Nodes().begin() + i_node;
                double& distance = it_node->FastGetSolutionStepValue(DISTANCE);
                double value = DistanceFunction(it_node);
                KRATOS_CHECK_NEAR(value, distance, Tolerance);
            }
        } // UniformRefineTrianglesUtility


        /**
         * Checks the correct refining utility with quadrilaterals
         */
        KRATOS_TEST_CASE_IN_SUITE(UniformRefineQuadrilateralsUtility, KratosMeshingApplicationFastSuite)
        {
            Model this_model;
            ModelPart& this_model_part = this_model.CreateModelPart("Main", 2);

            this_model_part.AddNodalSolutionStepVariable(VELOCITY);

            Properties::Pointer p_properties = this_model_part.pGetProperties(0);

            // Creating the sub model parts
            ModelPart& r_sub_model_part_1 = this_model_part.CreateSubModelPart("BodySubModelPart");

            // Creating the nodes
            NodeType::Pointer p_node_1 = this_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
            NodeType::Pointer p_node_2 = this_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
            NodeType::Pointer p_node_3 = this_model_part.CreateNewNode(3, 2.0, 0.0, 0.0);
            NodeType::Pointer p_node_4 = this_model_part.CreateNewNode(4, 0.0, 1.0, 0.0);
            NodeType::Pointer p_node_5 = this_model_part.CreateNewNode(5, 1.0, 1.0, 0.0);
            NodeType::Pointer p_node_6 = this_model_part.CreateNewNode(6, 2.0, 1.0, 0.0);

            // Creating the elements
            std::vector<NodeType::Pointer> quadrilateral_nodes(4);
            quadrilateral_nodes[0] = p_node_1;
            quadrilateral_nodes[1] = p_node_2;
            quadrilateral_nodes[2] = p_node_5;
            quadrilateral_nodes[3] = p_node_4;
            Element::Pointer p_elem_1 = this_model_part.CreateNewElement("Element2D4N", 1, PointerVector<NodeType>{quadrilateral_nodes}, p_properties);

            quadrilateral_nodes[0] = p_node_2;
            quadrilateral_nodes[1] = p_node_3;
            quadrilateral_nodes[2] = p_node_6;
            quadrilateral_nodes[3] = p_node_5;
            Element::Pointer p_elem_2 = this_model_part.CreateNewElement("Element2D4N", 2, PointerVector<NodeType>{quadrilateral_nodes}, p_properties);

            // Adding elements and its nodes to the sub model part
            r_sub_model_part_1.AddNode(p_node_1);
            r_sub_model_part_1.AddNode(p_node_2);
            r_sub_model_part_1.AddNode(p_node_3);
            r_sub_model_part_1.AddNode(p_node_4);
            r_sub_model_part_1.AddNode(p_node_5);
            r_sub_model_part_1.AddNode(p_node_6);
            r_sub_model_part_1.AddElement(p_elem_1);
            r_sub_model_part_1.AddElement(p_elem_2);

            // Setting the variables
            for (std::size_t i_node = 0; i_node < this_model_part.Nodes().size(); ++i_node) {
                auto it_node = this_model_part.Nodes().begin() + i_node;
                double& velocity_x = it_node->FastGetSolutionStepValue(VELOCITY_X);
                double& velocity_y = it_node->FastGetSolutionStepValue(VELOCITY_Y);
                velocity_x = DistanceFunction(it_node);
                velocity_y = DistanceFunction(it_node);
            }

            // Store the initial values
            std::vector<unsigned int> initial_nodes(3);
            std::vector<unsigned int> initial_elems(3);
            initial_nodes[0] = this_model_part.NumberOfNodes();
            initial_elems[0] = this_model_part.NumberOfElements();
            initial_nodes[1] = r_sub_model_part_1.NumberOfNodes();
            initial_elems[1] = r_sub_model_part_1.NumberOfElements();

            // Execute the utility
            int refinement_level = 3;
            UniformRefinementUtility uniform_refinement(this_model_part);
            uniform_refinement.Refine(refinement_level);

            // Check the number of entities in the main model part
            unsigned int final_nodes = (std::pow(2,refinement_level)+1) * ((0.5*initial_nodes[0]-1)*std::pow(2,refinement_level)+1);
            KRATOS_CHECK_EQUAL(final_nodes, this_model_part.NumberOfNodes());

            unsigned int final_elems = initial_elems[0] * std::pow(4,refinement_level);
            KRATOS_CHECK_EQUAL(final_elems, this_model_part.NumberOfElements());

            // Check the number of entities in the first sub model part
            // final_nodes = (std::pow(2,refinement_level)+1) * ((0.5*initial_nodes[1]-1)*std::pow(2,refinement_level)+1);
            // KRATOS_CHECK_EQUAL(final_nodes, r_sub_model_part_1.NumberOfNodes());

            final_elems = initial_elems[1] * std::pow(4,refinement_level);
            KRATOS_CHECK_EQUAL(final_elems, r_sub_model_part_1.NumberOfElements());

            // Check the variables interpolation
            for (std::size_t i_node = 0; i_node < this_model_part.Nodes().size(); ++i_node) {
                auto it_node = this_model_part.Nodes().begin() + i_node;
                const double& velocity_x = it_node->FastGetSolutionStepValue(VELOCITY_X);
                const double& velocity_y = it_node->FastGetSolutionStepValue(VELOCITY_Y);
                double value = DistanceFunction(it_node);
                KRATOS_CHECK_NEAR(value, velocity_x, Tolerance);
                KRATOS_CHECK_NEAR(value, velocity_y, Tolerance);
            }
        } // UniformRefineQuadrilateralsUtility


        /**
         * Checks the correct refining utility with hexahedrons
         */
        KRATOS_TEST_CASE_IN_SUITE(UniformRefineHexahedronsUtility, KratosMeshingApplicationFastSuite)
        {
            Model this_model;
            ModelPart& this_model_part = this_model.CreateModelPart("Main", 2);

            this_model_part.AddNodalSolutionStepVariable(VELOCITY);

            Properties::Pointer p_properties = this_model_part.pGetProperties(0);

            // Creating the sub model parts
            ModelPart& r_sub_model_part_1 = this_model_part.CreateSubModelPart("BodySubModelPart");
            ModelPart& r_sub_model_part_2 = this_model_part.CreateSubModelPart("SkinSubModelPart");

            // Creating the nodes
            NodeType::Pointer p_node_1 = this_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
            NodeType::Pointer p_node_2 = this_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
            NodeType::Pointer p_node_3 = this_model_part.CreateNewNode(3, 1.0, 1.0, 0.0);
            NodeType::Pointer p_node_4 = this_model_part.CreateNewNode(4, 0.0, 1.0, 0.0);
            NodeType::Pointer p_node_5 = this_model_part.CreateNewNode(5, 0.0, 0.0, 1.0);
            NodeType::Pointer p_node_6 = this_model_part.CreateNewNode(6, 1.0, 0.0, 1.0);
            NodeType::Pointer p_node_7 = this_model_part.CreateNewNode(7, 1.0, 1.0, 1.0);
            NodeType::Pointer p_node_8 = this_model_part.CreateNewNode(8, 0.0, 1.0, 1.0);

            // Creating the elements
            PointerVector<NodeType> hexahedra_nodes(8);
            hexahedra_nodes(0) = p_node_1;
            hexahedra_nodes(1) = p_node_2;
            hexahedra_nodes(2) = p_node_3;
            hexahedra_nodes(3) = p_node_4;
            hexahedra_nodes(4) = p_node_5;
            hexahedra_nodes(5) = p_node_6;
            hexahedra_nodes(6) = p_node_7;
            hexahedra_nodes(7) = p_node_8;
            Element::Pointer p_elem_1 = this_model_part.CreateNewElement("Element3D8N", 1, hexahedra_nodes, p_properties);

            // Creating the conditions
            PointerVector<NodeType> quadrilateral_nodes(4);
            quadrilateral_nodes(0) = p_node_1;
            quadrilateral_nodes(1) = p_node_2;
            quadrilateral_nodes(2) = p_node_3;
            quadrilateral_nodes(3) = p_node_4;
            Condition::Pointer p_cond_1 = this_model_part.CreateNewCondition("SurfaceCondition3D4N", 1, quadrilateral_nodes, p_properties);

            r_sub_model_part_1.AddElement(p_elem_1);
            r_sub_model_part_2.AddCondition(p_cond_1);

            // Get the initial values
            unsigned int initial_elements = r_sub_model_part_1.NumberOfElements();
            unsigned int initial_conditions = r_sub_model_part_2.NumberOfConditions();

            // Execute the utility
            int refinement_level = 2;
            UniformRefinementUtility uniform_refinement(this_model_part);
            uniform_refinement.Refine(refinement_level);

            // Check the number of elements (tetrahedrons)
            unsigned int final_elements = initial_elements * std::pow(8, refinement_level);
            KRATOS_CHECK_EQUAL(final_elements, r_sub_model_part_1.NumberOfElements());

            // Check the number of conditions (quadrilaterals)
            unsigned int final_conditions = initial_conditions * std::pow(4, refinement_level);
            KRATOS_CHECK_EQUAL(final_conditions, r_sub_model_part_2.NumberOfConditions());
        } // UniformRefineHexahedronsUtility



        /**
         * Checks the correct refining utility with tetrahedrons
         */
        KRATOS_TEST_CASE_IN_SUITE(UniformRefineTetrahedronsUtility, KratosMeshingApplicationFastSuite)
        {
            Model this_model;
            ModelPart& this_model_part = this_model.CreateModelPart("Main", 2);

            this_model_part.AddNodalSolutionStepVariable(VELOCITY);

            Properties::Pointer p_properties = this_model_part.pGetProperties(0);

            // Creating the sub model parts
            ModelPart& r_sub_model_part_1 = this_model_part.CreateSubModelPart("BodySubModelPart");
            ModelPart& r_sub_model_part_2 = this_model_part.CreateSubModelPart("SkinSubModelPart");

            // Creating the nodes
            NodeType::Pointer p_node_1 = this_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
            NodeType::Pointer p_node_2 = this_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
            NodeType::Pointer p_node_3 = this_model_part.CreateNewNode(3, 1.0, 1.0, 0.0);
            NodeType::Pointer p_node_4 = this_model_part.CreateNewNode(4, 0.0, 1.0, 0.0);
            NodeType::Pointer p_node_5 = this_model_part.CreateNewNode(5, 0.0, 0.0, 1.0);
            NodeType::Pointer p_node_6 = this_model_part.CreateNewNode(6, 1.0, 0.0, 1.0);
            NodeType::Pointer p_node_7 = this_model_part.CreateNewNode(7, 1.0, 1.0, 1.0);
            NodeType::Pointer p_node_8 = this_model_part.CreateNewNode(8, 0.0, 1.0, 1.0);

            // Creating the elements
            PointerVector<NodeType> tetrahedra_nodes_1(4);
            tetrahedra_nodes_1(0) = p_node_1;
            tetrahedra_nodes_1(1) = p_node_2;
            tetrahedra_nodes_1(2) = p_node_3;
            tetrahedra_nodes_1(3) = p_node_5;

            PointerVector<NodeType> tetrahedra_nodes_2(4);
            tetrahedra_nodes_2(0) = p_node_2;
            tetrahedra_nodes_2(1) = p_node_6;
            tetrahedra_nodes_2(2) = p_node_3;
            tetrahedra_nodes_2(3) = p_node_5;

            PointerVector<NodeType> tetrahedra_nodes_3(4);
            tetrahedra_nodes_3(0) = p_node_3;
            tetrahedra_nodes_3(1) = p_node_6;
            tetrahedra_nodes_3(2) = p_node_7;
            tetrahedra_nodes_3(3) = p_node_5;

            PointerVector<NodeType> tetrahedra_nodes_4(4);
            tetrahedra_nodes_4(0) = p_node_3;
            tetrahedra_nodes_4(1) = p_node_7;
            tetrahedra_nodes_4(2) = p_node_8;
            tetrahedra_nodes_4(3) = p_node_5;

            PointerVector<NodeType> tetrahedra_nodes_5(4);
            tetrahedra_nodes_5(0) = p_node_3;
            tetrahedra_nodes_5(1) = p_node_8;
            tetrahedra_nodes_5(2) = p_node_4;
            tetrahedra_nodes_5(3) = p_node_5;

            PointerVector<NodeType> tetrahedra_nodes_6(4);
            tetrahedra_nodes_6(0) = p_node_1;
            tetrahedra_nodes_6(1) = p_node_3;
            tetrahedra_nodes_6(2) = p_node_4;
            tetrahedra_nodes_6(3) = p_node_5;

            Element::Pointer p_elem_1 = this_model_part.CreateNewElement("Element3D4N", 1, tetrahedra_nodes_1, p_properties);
            Element::Pointer p_elem_2 = this_model_part.CreateNewElement("Element3D4N", 2, tetrahedra_nodes_2, p_properties);
            Element::Pointer p_elem_3 = this_model_part.CreateNewElement("Element3D4N", 3, tetrahedra_nodes_3, p_properties);
            Element::Pointer p_elem_4 = this_model_part.CreateNewElement("Element3D4N", 4, tetrahedra_nodes_4, p_properties);
            Element::Pointer p_elem_5 = this_model_part.CreateNewElement("Element3D4N", 5, tetrahedra_nodes_5, p_properties);
            Element::Pointer p_elem_6 = this_model_part.CreateNewElement("Element3D4N", 6, tetrahedra_nodes_6, p_properties);

            // Creating the conditions
            PointerVector<NodeType> triangle_nodes_1(3);
            triangle_nodes_1(0) = p_node_1;
            triangle_nodes_1(1) = p_node_5;
            triangle_nodes_1(2) = p_node_4;

            PointerVector<NodeType> triangle_nodes_2(3);
            triangle_nodes_2(0) = p_node_4;
            triangle_nodes_2(1) = p_node_5;
            triangle_nodes_2(2) = p_node_8;

            Condition::Pointer p_cond_1 = this_model_part.CreateNewCondition("SurfaceCondition3D3N", 1, triangle_nodes_1, p_properties);
            Condition::Pointer p_cond_2 = this_model_part.CreateNewCondition("SurfaceCondition3D3N", 2, triangle_nodes_2, p_properties);

            r_sub_model_part_1.AddElement(p_elem_1);
            r_sub_model_part_1.AddElement(p_elem_2);
            r_sub_model_part_1.AddElement(p_elem_3);
            r_sub_model_part_2.AddCondition(p_cond_1);
            r_sub_model_part_2.AddCondition(p_cond_2);

            // Get the initial values
            unsigned int initial_elements = r_sub_model_part_1.NumberOfElements();
            unsigned int initial_conditions = r_sub_model_part_2.NumberOfConditions();

            // Execute the utility
            int refinement_level = 2;
            UniformRefinementUtility uniform_refinement(this_model_part);
            uniform_refinement.Refine(refinement_level);

            // Check the number of elements (tetrahedrons)
            unsigned int final_elements = initial_elements * std::pow(8, refinement_level);
            KRATOS_CHECK_EQUAL(final_elements, r_sub_model_part_1.NumberOfElements());

            // Check the number of conditions (quadrilaterals)
            unsigned int final_conditions = initial_conditions * std::pow(4, refinement_level);
            KRATOS_CHECK_EQUAL(final_conditions, r_sub_model_part_2.NumberOfConditions());
        } // UniformRefineTetrahedronsUtility

    } // namespace Testing
} // namespace Kratos
