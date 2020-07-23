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
#include "containers/model.h"
#include "custom_utilities/interface_preprocess.h"
#include "utilities/cpp_tests_utilities.h"

namespace Kratos {
    namespace Testing {

        typedef Node<3> NodeType;
        typedef Geometry<NodeType> GeometryType;

        void CreateSimple2DModelPart(ModelPart& rModelPart)
        {
            CppTestsUtilities::Create2DGeometry(rModelPart, "Element2D3N");

            // Set flag
            for (auto& r_node : rModelPart.Nodes()) {
                if (r_node.Y() == 0.0) {
                    r_node.Set(INTERFACE);
                }
            }
        }

        void CreateSimple3DModelPart(ModelPart& rModelPart)
        {
            CppTestsUtilities::Create3DGeometry(rModelPart, "Element3D4N");

            // Set flag
            for (auto& r_node : rModelPart.Nodes()) {
                if (r_node.Z() == 0.0) {
                    r_node.Set(INTERFACE);
                }
            }
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
