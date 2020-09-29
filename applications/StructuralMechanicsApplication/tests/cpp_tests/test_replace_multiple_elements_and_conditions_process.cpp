// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Armin Geiser
//


// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/triangle_3d_3.h"
#include "containers/model.h"

/* Processes */
#include "custom_response_functions/adjoint_processes/replace_multiple_elements_and_conditions_process.h"
#include "utilities/compare_elements_and_conditions_utility.h"

namespace Kratos
{
    namespace Testing
    {
        typedef Node<3> NodeType;

        void ReplaceMultipleProcessCreateModelPart(ModelPart& ThisModelPart)
        {
            Properties::Pointer p_prop = ThisModelPart.CreateNewProperties(0);

            NodeType::Pointer p_node_1 = ThisModelPart.CreateNewNode(1, 0.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_2 = ThisModelPart.CreateNewNode(2, 1.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_3 = ThisModelPart.CreateNewNode(3, 1.0 , 1.0 , 0.0);

            std::vector<NodeType::Pointer> element_nodes_0 (3);
            element_nodes_0[0] = p_node_1;
            element_nodes_0[1] = p_node_2;
            element_nodes_0[2] = p_node_3;
            Triangle3D3 <NodeType> triangle_0( PointerVector<NodeType>{element_nodes_0} );
            ThisModelPart.CreateNewElement("Element3D3N", 1, triangle_0, p_prop);

            std::vector<NodeType::Pointer> element_nodes_1 (3);
            element_nodes_1[0] = p_node_2;
            element_nodes_1[1] = p_node_3;
            element_nodes_1[2] = p_node_1;
            Triangle3D3 <NodeType> triangle_1( PointerVector<NodeType>{element_nodes_1} );
            ThisModelPart.CreateNewElement("ShellThinElement3D3N", 2, triangle_1, p_prop);

            ThisModelPart.CreateNewCondition("SurfaceCondition3D3N", 1, triangle_0, p_prop);
            ThisModelPart.CreateNewCondition("SurfaceLoadCondition3D3N", 2, triangle_1, p_prop);
        }

        KRATOS_TEST_CASE_IN_SUITE(ReplaceMultipleElementsAndConditionsProcess1, KratosStructuralMechanicsFastSuite)
        {
            // all types are explicitly defined in the replacement table
            Model current_model;
            ModelPart& this_model_part =  current_model.CreateModelPart("Main");

            ReplaceMultipleProcessCreateModelPart(this_model_part);

            Parameters parameters = Parameters(R"(
            {
                "element_name_table"    : {
                    "Element3D3N" : "ShellThinElement3D3N",
                    "ShellThinElement3D3N": "ShellThinElement3D3N"
                },
                "condition_name_table"    : {
                    "SurfaceCondition3D3N" : "SurfaceLoadCondition3D3N",
                    "SurfaceLoadCondition3D3N": "SurfaceLoadCondition3D3N"
                }
            })" );

            ReplaceMultipleElementsAndConditionsProcess replacement_process(this_model_part, parameters);
            replacement_process.Execute();

            std::string component_name;

            CompareElementsAndConditionsUtility::GetRegisteredName(this_model_part.GetElement(1), component_name);
            KRATOS_CHECK_EQUAL(component_name, "ShellThinElement3D3N");

            CompareElementsAndConditionsUtility::GetRegisteredName(this_model_part.GetElement(2), component_name);
            KRATOS_CHECK_EQUAL(component_name, "ShellThinElement3D3N");

            CompareElementsAndConditionsUtility::GetRegisteredName(this_model_part.GetCondition(1), component_name);
            KRATOS_CHECK_EQUAL(component_name, "SurfaceLoadCondition3D3N");

            CompareElementsAndConditionsUtility::GetRegisteredName(this_model_part.GetCondition(2), component_name);
            KRATOS_CHECK_EQUAL(component_name, "SurfaceLoadCondition3D3N");
        }

        // Does not work because KRATOS_CHECK_EXCEPTION_IS_THROWN works only for errors thrown in serial regions
        // KRATOS_TEST_CASE_IN_SUITE(ReplaceMultipleElementsAndConditionsProcess2, KratosStructuralMechanicsFastSuite)
        // {
        //     Model current_model;
        //     ModelPart& this_model_part =  current_model.CreateModelPart("Main");

        //     ReplaceMultipleProcessCreateModelPart(this_model_part);

        //     Parameters parameters = Parameters(R"(
        //     {
        //         "element_name_table"    : {
        //             "ShellThinElement3D3N" : "ShellThinElement3D3N"
        //         }
        //     })" );

        //     ReplaceMultipleElementsAndConditionsProcess replacement_process(this_model_part, parameters);

        //     KRATOS_CHECK_EXCEPTION_IS_THROWN(
        //         replacement_process.Execute(),
        //         "Element3D3N was not defined in the replacement table!"
        //     );

        // }

        KRATOS_TEST_CASE_IN_SUITE(ReplaceMultipleElementsAndConditionsProcess3, KratosStructuralMechanicsFastSuite)
        {
            // undefiend types are ignored by the replacement process
            Model current_model;
            ModelPart& this_model_part =  current_model.CreateModelPart("Main");

            ReplaceMultipleProcessCreateModelPart(this_model_part);

            Parameters parameters = Parameters(R"(
            {
                "element_name_table"    : {
                    "Element3D3N" : "ShellThinElement3D3N"
                },
                "condition_name_table"    : {
                    "SurfaceCondition3D3N" : "SurfaceLoadCondition3D3N"
                },
                "ignore_undefined_types" : true
            })" );

            ReplaceMultipleElementsAndConditionsProcess replacement_process(this_model_part, parameters);
            replacement_process.Execute();

            std::string component_name;

            CompareElementsAndConditionsUtility::GetRegisteredName(this_model_part.GetElement(1), component_name);
            KRATOS_CHECK_EQUAL(component_name, "ShellThinElement3D3N");

            CompareElementsAndConditionsUtility::GetRegisteredName(this_model_part.GetElement(2), component_name);
            KRATOS_CHECK_EQUAL(component_name, "ShellThinElement3D3N");

            CompareElementsAndConditionsUtility::GetRegisteredName(this_model_part.GetCondition(1), component_name);
            KRATOS_CHECK_EQUAL(component_name, "SurfaceLoadCondition3D3N");

            CompareElementsAndConditionsUtility::GetRegisteredName(this_model_part.GetCondition(2), component_name);
            KRATOS_CHECK_EQUAL(component_name, "SurfaceLoadCondition3D3N");
        }


        KRATOS_TEST_CASE_IN_SUITE(ReplaceMultipleElementsAndConditionsProcess4, KratosStructuralMechanicsFastSuite)
        {
            // ignored types are explicitly defined in a list
            Model current_model;
            ModelPart& this_model_part =  current_model.CreateModelPart("Main");

            ReplaceMultipleProcessCreateModelPart(this_model_part);

            Parameters parameters = Parameters(R"(
            {
                "element_name_table"    : {
                    "Element3D3N" : "ShellThinElement3D3N"
                },
                "condition_name_table"    : {
                    "SurfaceCondition3D3N" : "SurfaceLoadCondition3D3N"
                },
                "ignore_elements" : [
                    "ShellThinElement3D3N"
                ],
                "ignore_conditions" : [
                    "SurfaceLoadCondition3D3N"
                ]
            })" );

            ReplaceMultipleElementsAndConditionsProcess replacement_process(this_model_part, parameters);
            replacement_process.Execute();

            std::string component_name;

            CompareElementsAndConditionsUtility::GetRegisteredName(this_model_part.GetElement(1), component_name);
            KRATOS_CHECK_EQUAL(component_name, "ShellThinElement3D3N");

            CompareElementsAndConditionsUtility::GetRegisteredName(this_model_part.GetElement(2), component_name);
            KRATOS_CHECK_EQUAL(component_name, "ShellThinElement3D3N");

            CompareElementsAndConditionsUtility::GetRegisteredName(this_model_part.GetCondition(1), component_name);
            KRATOS_CHECK_EQUAL(component_name, "SurfaceLoadCondition3D3N");

            CompareElementsAndConditionsUtility::GetRegisteredName(this_model_part.GetCondition(2), component_name);
            KRATOS_CHECK_EQUAL(component_name, "SurfaceLoadCondition3D3N");
        }

    } // namespace Testing
}  // namespace Kratos.
