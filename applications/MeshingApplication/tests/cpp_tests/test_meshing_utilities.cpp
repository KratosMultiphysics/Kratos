// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                       license: MeshingApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "utilities/cpp_tests_utilities.h"
#include "custom_utilities/meshing_utilities.h"

namespace Kratos {
    namespace Testing {

        typedef Node<3> NodeType;
        typedef Geometry<NodeType> GeometryType;

        void CreateDummy2DNoModelPartPropertiesModelPart(ModelPart& rModelPart)
        {
            Properties::Pointer p_elem_prop = Kratos::make_shared<Properties>(0);

            // First we create the nodes
            rModelPart.CreateNewNode(1, 0.0 , 0.0 , 0.0);
            rModelPart.CreateNewNode(2, 1.0 , 0.0 , 0.0);
            rModelPart.CreateNewNode(3, 1.0 , 1.0 , 0.0);
            rModelPart.CreateNewNode(4, 0.0 , 1.0 , 0.0);
            rModelPart.CreateNewNode(5, 2.0 , 0.0 , 0.0);
            rModelPart.CreateNewNode(6, 2.0 , 1.0 , 0.0);

            // Now we create the elements
            rModelPart.CreateNewElement("Element2D3N", 1, {{1,2,3}}, p_elem_prop);
            rModelPart.CreateNewElement("Element2D3N", 2, {{1,3,4}}, p_elem_prop);
            rModelPart.CreateNewElement("Element2D3N", 3, {{2,5,3}}, p_elem_prop);
            rModelPart.CreateNewElement("Element2D3N", 4, {{5,6,3}}, p_elem_prop);
        }

        /**
        * Checks the correct work of the EnsureModelPartOwnsProperties
        * Test triangle
        */
        KRATOS_TEST_CASE_IN_SUITE(EnsureModelPartOwnsProperties, KratosMeshingApplicationFastSuite)
        {
            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);
            ProcessInfo& r_current_process_info = r_model_part.GetProcessInfo();
            r_current_process_info[DOMAIN_SIZE] = 2;

            CreateDummy2DNoModelPartPropertiesModelPart(r_model_part);

            KRATOS_CHECK(r_model_part.NumberOfProperties() == 0);

            MeshingUtilities::EnsureModelPartOwnsProperties(r_model_part);

            KRATOS_CHECK(r_model_part.NumberOfProperties() == 1);
        }

        /**
        * Checks the correct work of the BlockThresholdSizeElements
        * Test triangle
        */
        KRATOS_TEST_CASE_IN_SUITE(BlockThresholdSizeElements2D, KratosMeshingApplicationFastSuite)
        {
            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);
            ProcessInfo& r_current_process_info = r_model_part.GetProcessInfo();
            r_current_process_info[DOMAIN_SIZE] = 2;

            CppTestsUtilities::Create2DGeometry(r_model_part, "Element2D3N");

            Parameters parameters = Parameters(R"(
            {
                "minimal_size" : 2.0,
                "maximal_size" : 10.0
            })" );

            MeshingUtilities::BlockThresholdSizeElements(r_model_part, parameters);

            for (auto& r_element: r_model_part.Elements()) {
                KRATOS_CHECK(r_element.Is(BLOCKED));
//                 KRATOS_WATCH(r_element.GetValue(ELEMENT_H))
            }
        }

        /**
        * Checks the correct work of the BlockThresholdSizeElements
        * Test tetrahedra
        */
        KRATOS_TEST_CASE_IN_SUITE(BlockThresholdSizeElements3D, KratosMeshingApplicationFastSuite)
        {
            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);
            ProcessInfo& r_current_process_info = r_model_part.GetProcessInfo();
            r_current_process_info[DOMAIN_SIZE] = 3;

            CppTestsUtilities::Create3DGeometry(r_model_part, "Element3D4N");

            Parameters parameters = Parameters(R"(
            {
                "minimal_size" : 2.0,
                "maximal_size" : 10.0
            })" );

            MeshingUtilities::BlockThresholdSizeElements(r_model_part, parameters);

            for (auto& r_element: r_model_part.Elements()) {
                KRATOS_CHECK(r_element.Is(BLOCKED));
//                 KRATOS_WATCH(r_element.GetValue(ELEMENT_H))
            }
        }

    } // namespace Testing
} // namespace Kratos
