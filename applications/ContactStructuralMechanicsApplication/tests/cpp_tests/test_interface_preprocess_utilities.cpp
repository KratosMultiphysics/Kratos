// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "containers/model.h"
#include "custom_utilities/interface_preprocess.h"

namespace Kratos {
    namespace Testing {

        typedef Node<3> NodeType;
        typedef Geometry<NodeType> GeometryType;

        void CreateSimple2DModelPart(ModelPart& rModelPart)
        {
            Properties::Pointer p_elem_prop = rModelPart.CreateNewProperties(0);

            // First we create the nodes
            rModelPart.CreateNewNode(1, 0.0 , 0.0 , 0.0);
            rModelPart.CreateNewNode(2, 1.0 , 0.0 , 0.0);
            rModelPart.CreateNewNode(3, 1.0 , 1.0 , 0.0);
            rModelPart.CreateNewNode(4, 0.0 , 1.0 , 0.0);
            rModelPart.CreateNewNode(5, 2.0 , 0.0 , 0.0);
            rModelPart.CreateNewNode(6, 2.0 , 1.0 , 0.0);

            // Set flag
            for (auto& r_node : rModelPart.Nodes()) {
                if (r_node.Y() == 0.0) {
                    r_node.Set(INTERFACE);
                }
            }

            // Now we create the elements
            rModelPart.CreateNewElement("Element2D3N", 1, {{1,2,3}}, p_elem_prop);
            rModelPart.CreateNewElement("Element2D3N", 2, {{1,3,4}}, p_elem_prop);
            rModelPart.CreateNewElement("Element2D3N", 3, {{2,5,3}}, p_elem_prop);
            rModelPart.CreateNewElement("Element2D3N", 4, {{5,6,3}}, p_elem_prop);
        }

        void CreateSimple3DModelPart(ModelPart& rModelPart)
        {
            Properties::Pointer p_elem_prop = rModelPart.CreateNewProperties(0);

            // First we create the nodes
            rModelPart.CreateNewNode(1 , 0.0 , 1.0 , 1.0);
            rModelPart.CreateNewNode(2 , 0.0 , 1.0 , 0.0);
            rModelPart.CreateNewNode(3 , 0.0 , 0.0 , 1.0);
            rModelPart.CreateNewNode(4 , 1.0 , 1.0 , 1.0);
            rModelPart.CreateNewNode(5 , 0.0 , 0.0 , 0.0);
            rModelPart.CreateNewNode(6 , 1.0 , 1.0 , 0.0);

            rModelPart.CreateNewNode(7 , 1.0 , 0.0 , 1.0);
            rModelPart.CreateNewNode(8 , 1.0 , 0.0 , 0.0);
            rModelPart.CreateNewNode(9 , 2.0 , 1.0 , 1.0);
            rModelPart.CreateNewNode(10 , 2.0 , 1.0 , 0.0);
            rModelPart.CreateNewNode(11 , 2.0 , 0.0 , 1.0);
            rModelPart.CreateNewNode(12 , 2.0 , 0.0 , 0.0);

            // Set flag
            for (auto& r_node : rModelPart.Nodes()) {
                if (r_node.Z() == 0.0) {
                    r_node.Set(INTERFACE);
                }
            }

            // Now we create the elements
            rModelPart.CreateNewElement("Element3D4N", 1, {{12,10,8,9}}, p_elem_prop);
            rModelPart.CreateNewElement("Element3D4N", 2, {{4,6,9,7}}, p_elem_prop);
            rModelPart.CreateNewElement("Element3D4N", 3, {{11,7,9,8}}, p_elem_prop);
            rModelPart.CreateNewElement("Element3D4N", 4, {{5,3,8,6}}, p_elem_prop);
            rModelPart.CreateNewElement("Element3D4N", 5, {{4,6,7,3}}, p_elem_prop);
            rModelPart.CreateNewElement("Element3D4N", 6, {{2,3,5,6}}, p_elem_prop);
            rModelPart.CreateNewElement("Element3D4N", 7, {{10,9,6,8}}, p_elem_prop);
            rModelPart.CreateNewElement("Element3D4N", 8, {{7,8,3,6}}, p_elem_prop);
            rModelPart.CreateNewElement("Element3D4N", 9, {{7,8,6,9}}, p_elem_prop);
            rModelPart.CreateNewElement("Element3D4N", 10, {{4,1,6,3}}, p_elem_prop);
            rModelPart.CreateNewElement("Element3D4N", 11, {{9,12,11,8}}, p_elem_prop);
            rModelPart.CreateNewElement("Element3D4N", 12, {{3,2,1,6}}, p_elem_prop);
        }

        /**
        * Checks the correct work of the InterfacePreprocessCondition
        * Test triangle
        */
        KRATOS_TEST_CASE_IN_SUITE(InterfacePreprocessCondition2D, KratosContactStructuralMechanicsFastSuite)
        {
            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);
            ProcessInfo& r_current_process_info = r_model_part.GetProcessInfo();
            r_current_process_info[DOMAIN_SIZE] = 2;

            CreateSimple2DModelPart(r_model_part);

            auto utility = InterfacePreprocessCondition(r_model_part);
            utility.GenerateInterfacePart(r_model_part);

            KRATOS_CHECK(r_model_part.NumberOfConditions() == 2);
        }

        /**
        * Checks the correct work of the InterfacePreprocessCondition
        * Test tetrahedra
        */
        KRATOS_TEST_CASE_IN_SUITE(InterfacePreprocessCondition3D, KratosContactStructuralMechanicsFastSuite)
        {
            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);
            ProcessInfo& r_current_process_info = r_model_part.GetProcessInfo();
            r_current_process_info[DOMAIN_SIZE] = 3;

            CreateSimple3DModelPart(r_model_part);

            auto utility = InterfacePreprocessCondition(r_model_part);
            utility.GenerateInterfacePart(r_model_part);

            KRATOS_CHECK(r_model_part.NumberOfConditions() == 4);
        }

    } // namespace Testing
} // namespace Kratos
