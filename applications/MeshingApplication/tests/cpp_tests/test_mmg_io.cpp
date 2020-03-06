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

#ifdef INCLUDE_MMG
// System includes
#include<iostream>

// External includes

// Project includes
#include "testing/testing.h"
#include "includes/checks.h"
#include "includes/kratos_filesystem.h"
#include "includes/kratos_flags.h"
#include "containers/model.h"
#include "utilities/cpp_tests_utilities.h"
#include "meshing_application_variables.h"
#include "custom_io/mmg_io.h"

namespace Kratos
{
    namespace Testing
    {
        /**
        * Checks the correct work of the level set MMG IO
        * Test triangle
        */
        KRATOS_TEST_CASE_IN_SUITE(MMGIO1, KratosMeshingApplicationFastSuite)
        {
            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);

            Properties::Pointer p_elem_prop = r_model_part.CreateNewProperties(0);

            auto& process_info = r_model_part.GetProcessInfo();
            process_info[STEP] = 1;
            process_info[NL_ITERATION_NUMBER] = 1;

            CppTestsUtilities::Create2DGeometry(r_model_part, "Element2D3N");

            // Set DISTANCE and other variables
            Vector ref_metric(3);
            ref_metric[0] = 1.0;
            ref_metric[1] = 1.0;
            ref_metric[2] = 0.0;
            for (std::size_t i_node = 0; i_node < r_model_part.Nodes().size(); ++i_node) {
                auto it_node = r_model_part.Nodes().begin() + i_node;
                it_node->SetValue(METRIC_TENSOR_2D, ref_metric);
            }

            // Auxiliar submodelpart to check colors
            auto& r_sub = r_model_part.CreateSubModelPart("AuxiliarSubModelPart");

            // Now we create the conditions
            r_sub.AddNode(r_model_part.pGetNode(1));
            r_sub.AddNode(r_model_part.pGetNode(2));
            auto p_cond = r_sub.CreateNewCondition("LineCondition2D2N", 1, {{1,2}}, p_elem_prop);

            // Compute read/write
            Parameters params = Parameters(R"({ "echo_level" : 0 })" );
            MmgIO<MMGLibrary::MMG2D> mmg_io(FilesystemExtensions::JoinPaths({FilesystemExtensions::CurrentWorkingDirectory(), "mmg_output_2d"}));
            mmg_io.WriteModelPart(r_model_part);

            Model this_aux_model;
            ModelPart& r_aux_model_part = this_aux_model.CreateModelPart("Main", 2);

            mmg_io.ReadModelPart(r_aux_model_part);

            KRATOS_CHECK_EQUAL(r_model_part.NumberOfNodes(), r_aux_model_part.NumberOfNodes());
            KRATOS_CHECK_EQUAL(r_model_part.NumberOfElements(), r_aux_model_part.NumberOfElements());

            for (auto& r_sub_model_part_name : r_model_part.GetSubModelPartNames()) {
                KRATOS_CHECK(r_aux_model_part.HasSubModelPart(r_sub_model_part_name));
            }

            Kratos::filesystem::remove((FilesystemExtensions::JoinPaths({FilesystemExtensions::CurrentWorkingDirectory(), "mmg_output_2d.mesh"})).c_str());
            Kratos::filesystem::remove((FilesystemExtensions::JoinPaths({FilesystemExtensions::CurrentWorkingDirectory(), "mmg_output_2d.sol"})).c_str());
            Kratos::filesystem::remove((FilesystemExtensions::JoinPaths({FilesystemExtensions::CurrentWorkingDirectory(), "mmg_output_2d.json"})).c_str());
            Kratos::filesystem::remove((FilesystemExtensions::JoinPaths({FilesystemExtensions::CurrentWorkingDirectory(), "mmg_output_2d.cond.ref.json"})).c_str());
            Kratos::filesystem::remove((FilesystemExtensions::JoinPaths({FilesystemExtensions::CurrentWorkingDirectory(), "mmg_output_2d.elem.ref.json"})).c_str());
        }

        /**
        * Checks the correct work of the level set MMG IO
        * Test tetrahedra
        */
        KRATOS_TEST_CASE_IN_SUITE(MMGIO2, KratosMeshingApplicationFastSuite)
        {
            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);

            r_model_part.AddNodalSolutionStepVariable(DISTANCE);
            r_model_part.AddNodalSolutionStepVariable(DISTANCE_GRADIENT);

            Properties::Pointer p_elem_prop = r_model_part.CreateNewProperties(0);

            auto& process_info = r_model_part.GetProcessInfo();
            process_info[STEP] = 1;
            process_info[NL_ITERATION_NUMBER] = 1;

            CppTestsUtilities::Create3DGeometry(r_model_part, "Element3D4N");

            // Set DISTANCE and other variables
            array_1d<double, 6> ref_metric = ZeroVector(6);
            ref_metric[0] = 1.0;
            ref_metric[1] = 1.0;
            ref_metric[2] = 1.0;
            for (std::size_t i_node = 0; i_node < r_model_part.Nodes().size(); ++i_node) {
                auto it_node = r_model_part.Nodes().begin() + i_node;
                it_node->SetValue(METRIC_TENSOR_3D, ref_metric);
            }

            // Auxiliar submodelpart to check colors
            auto& r_sub = r_model_part.CreateSubModelPart("AuxiliarSubModelPart");

            // Now we create the conditions
            r_sub.AddNode(r_model_part.pGetNode(1));
            r_sub.AddNode(r_model_part.pGetNode(2));
            r_sub.AddNode(r_model_part.pGetNode(3));
            auto p_cond = r_sub.CreateNewCondition("SurfaceCondition3D3N", 1, {{1,2,3}}, p_elem_prop);

            // Compute read/write
            Parameters params = Parameters(R"({ "echo_level" : 0 })" );
            MmgIO<MMGLibrary::MMG3D> mmg_io(FilesystemExtensions::JoinPaths({FilesystemExtensions::CurrentWorkingDirectory(), "mmg_output_3d"}));
            mmg_io.WriteModelPart(r_model_part);

            Model this_aux_model;
            ModelPart& r_aux_model_part = this_aux_model.CreateModelPart("Main", 2);

            mmg_io.ReadModelPart(r_aux_model_part);

            KRATOS_CHECK_EQUAL(r_model_part.NumberOfNodes(), r_aux_model_part.NumberOfNodes());
            KRATOS_CHECK_EQUAL(r_model_part.NumberOfElements(), r_aux_model_part.NumberOfElements());

            for (auto& r_sub_model_part_name : r_model_part.GetSubModelPartNames()) {
                KRATOS_CHECK(r_aux_model_part.HasSubModelPart(r_sub_model_part_name));
            }

            Kratos::filesystem::remove((FilesystemExtensions::JoinPaths({FilesystemExtensions::CurrentWorkingDirectory(), "mmg_output_3d.mesh"})).c_str());
            Kratos::filesystem::remove((FilesystemExtensions::JoinPaths({FilesystemExtensions::CurrentWorkingDirectory(), "mmg_output_3d.sol"})).c_str());
            Kratos::filesystem::remove((FilesystemExtensions::JoinPaths({FilesystemExtensions::CurrentWorkingDirectory(), "mmg_output_3d.json"})).c_str());
            Kratos::filesystem::remove((FilesystemExtensions::JoinPaths({FilesystemExtensions::CurrentWorkingDirectory(), "mmg_output_3d.cond.ref.json"})).c_str());
            Kratos::filesystem::remove((FilesystemExtensions::JoinPaths({FilesystemExtensions::CurrentWorkingDirectory(), "mmg_output_3d.elem.ref.json"})).c_str());
        }
    } // namespace Testing
}  // namespace Kratos.
#endif
