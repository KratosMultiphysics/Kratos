//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:   BSD License
//      Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
#ifdef INCLUDE_PMMG
// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "includes/kratos_flags.h"
// #include "includes/gid_io.h"
#include "containers/model.h"
#include "meshing_application_variables.h"

/* Processes */
#include "custom_processes/pmmg_process.h"
#include "processes/find_nodal_h_process.h"

namespace Kratos
{
    namespace Testing
    {
        typedef Node<3> NodeType;

        /**
        * Checks the correct work of the level set PMMG process
        * Test tetrahedra
        */
        KRATOS_TEST_CASE_IN_SUITE(PMMGProcess1, KratosMeshingApplicationFastSuite)
        {
            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);

            r_model_part.AddNodalSolutionStepVariable(DISTANCE);
            r_model_part.AddNodalSolutionStepVariable(DISTANCE_GRADIENT);

            Properties::Pointer p_elem_prop = r_model_part.CreateNewProperties(0);

            auto& process_info = r_model_part.GetProcessInfo();
            process_info[STEP] = 1;
            process_info[NL_ITERATION_NUMBER] = 1;

            // First we create the nodes
            r_model_part.CreateNewNode(1 , 0.0 , 1.0 , 1.0);
            r_model_part.CreateNewNode(2 , 0.0 , 1.0 , 0.0);
            r_model_part.CreateNewNode(3 , 0.0 , 0.0 , 1.0);
            r_model_part.CreateNewNode(4 , 1.0 , 1.0 , 1.0);
            r_model_part.CreateNewNode(5 , 0.0 , 0.0 , 0.0);
            r_model_part.CreateNewNode(6 , 1.0 , 1.0 , 0.0);

            r_model_part.CreateNewNode(7 , 1.0 , 0.0 , 1.0);
            r_model_part.CreateNewNode(8 , 1.0 , 0.0 , 0.0);
            r_model_part.CreateNewNode(9 , 2.0 , 1.0 , 1.0);
            r_model_part.CreateNewNode(10 , 2.0 , 1.0 , 0.0);
            r_model_part.CreateNewNode(11 , 2.0 , 0.0 , 1.0);
            r_model_part.CreateNewNode(12 , 2.0 , 0.0 , 0.0);

            // Now we create the elements
            r_model_part.CreateNewElement("Element3D4N", 1, {{12,10,8,9}}, p_elem_prop);
            r_model_part.CreateNewElement("Element3D4N", 2, {{4,6,9,7}}, p_elem_prop);
            r_model_part.CreateNewElement("Element3D4N", 3, {{11,7,9,8}}, p_elem_prop);
            r_model_part.CreateNewElement("Element3D4N", 4, {{5,3,8,6}}, p_elem_prop);
            r_model_part.CreateNewElement("Element3D4N", 5, {{4,6,7,3}}, p_elem_prop);
            r_model_part.CreateNewElement("Element3D4N", 6, {{2,3,5,6}}, p_elem_prop);
            r_model_part.CreateNewElement("Element3D4N", 7, {{10,9,6,8}}, p_elem_prop);
            r_model_part.CreateNewElement("Element3D4N", 8, {{7,8,3,6}}, p_elem_prop);
            r_model_part.CreateNewElement("Element3D4N", 9, {{7,8,6,9}}, p_elem_prop);
            r_model_part.CreateNewElement("Element3D4N", 10, {{4,1,6,3}}, p_elem_prop);
            r_model_part.CreateNewElement("Element3D4N", 11, {{9,12,11,8}}, p_elem_prop);
            r_model_part.CreateNewElement("Element3D4N", 12, {{3,2,1,6}}, p_elem_prop);

            // We set the flag to check that is transfered
            for (auto& i_elem : r_model_part.Elements())
                i_elem.Set(ACTIVE, true);

            // Set DISTANCE and other variables
            array_1d<double, 6> ref_metric = ZeroVector(6);
            ref_metric[0] = 1.0;
            ref_metric[1] = 1.0;
            ref_metric[2] = 1.0;
            for (std::size_t i_node = 0; i_node < r_model_part.Nodes().size(); ++i_node) {
                auto it_node = r_model_part.Nodes().begin() + i_node;
                it_node->SetValue(METRIC_TENSOR_3D, ref_metric);
            }

            // Compute remesh
            Parameters params = Parameters(R"({ "echo_level" : 0 })" );
            ParMmgProcess<PMMGLibrary::PMMG3D> parmmg_process(r_model_part, params);
            parmmg_process.Execute();

            // Compute NodalH
            auto process = FindNodalHProcess<FindNodalHSettings::SaveAsNonHistoricalVariable>(r_model_part);
            process.Execute();

//             // DEBUG
//             GiDIODebugMMG(r_model_part);

            double max = 0.0;
            for (auto& i_node : r_model_part.Nodes())
                if (i_node.GetValue(NODAL_H) > max)
                    max = i_node.GetValue(NODAL_H);

            const double tolerance = 1.0e-2;
            KRATOS_CHECK_LESS_EQUAL(std::abs(max - 1.0/std::sqrt(2.0))/max, tolerance);

            for (auto& i_elem : r_model_part.Elements())
                KRATOS_CHECK(i_elem.Is(ACTIVE));
        }
    } // namespace Testing
}  // namespace Kratos.
#endif
