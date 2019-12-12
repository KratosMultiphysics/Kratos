// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "includes/kratos_flags.h"
// #include "includes/gid_io.h"
#include "contact_structural_mechanics_application.h"
#include "custom_processes/contact_search_wrapper_process.h"

namespace Kratos
{
    namespace Testing
    {
        typedef Node<3> NodeType;

//         void GiDIOSearchDebug(ModelPart& rModelPart)
//         {
//             // Activating
//             for (auto& r_cond : rModelPart.Conditions()) {
//                 r_cond.Set(ACTIVE);
//             }
//
//             GidIO<> gid_io("TEST_SEARCH_PROCESS", GiD_PostBinary, SingleFile, WriteUndeformed,  WriteConditions);
//             const int nl_iter = rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER];
//             const double label = static_cast<double>(nl_iter);
//
//             gid_io.InitializeMesh(label);
//             gid_io.WriteMesh(rModelPart.GetMesh());
//             gid_io.WriteNodeMesh(rModelPart.GetMesh());
//             gid_io.FinalizeMesh();
//             gid_io.InitializeResults(label, rModelPart.GetMesh());
//             gid_io.WriteNodalFlags(ACTIVE, "ACTIVE", rModelPart.Nodes(), label);
//             gid_io.WriteNodalFlags(SLAVE, "SLAVE", rModelPart.Nodes(), label);
//             gid_io.WriteNodalResults(DISPLACEMENT, rModelPart.Nodes(), label, 0);
//             gid_io.WriteNodalResultsNonHistorical(NORMAL, rModelPart.Nodes(), label);
//             gid_io.WriteNodalResultsNonHistorical(TANGENT_XI, rModelPart.Nodes(), label);
//             gid_io.WriteNodalResultsNonHistorical(TANGENT_ETA, rModelPart.Nodes(), label);
//             gid_io.WriteNodalResultsNonHistorical(NORMAL_GAP, rModelPart.Nodes(), label);
//             gid_io.WriteNodalResultsNonHistorical(NODAL_AREA, rModelPart.Nodes(), label);
//             gid_io.WriteNodalResultsNonHistorical(AUXILIAR_COORDINATES, rModelPart.Nodes(), label);
//             gid_io.WriteNodalResults(NORMAL, rModelPart.Nodes(), label, 0);
//         }

        /**
         * This method can be used to create a plane/cylinder condition set
         */
        void CreatePlaneCilynderProblemForSearch(
            ModelPart& rModelPart,
            const std::size_t NumberOfDivisions,
            const double Lenght,
            const double Radius,
            const double Angle,
            const double Slope = 0.0,
            const bool MoveMesh = false
            )
        {
            Properties::Pointer p_cond_prop = rModelPart.CreateNewProperties(0);
            ModelPart& r_contact_model_part = rModelPart.CreateSubModelPart("Contact");

            double x, y;

            // Creating the base geometry
            std::size_t id_node = 0;
            const double dx = Lenght/static_cast<double>(NumberOfDivisions);
            for (std::size_t i = 0; i < NumberOfDivisions + 1; ++i) {
                x = dx * i;
                y = Slope * dx * i;
                id_node++;
                NodeType::Pointer p_node_1 = r_contact_model_part.CreateNewNode(id_node, x , y , 0.0);
                p_node_1->Set(SLAVE, true);
                p_node_1->Set(MASTER, false);
                id_node++;
                NodeType::Pointer p_node_2 = r_contact_model_part.CreateNewNode(id_node, x , y , 1.0);
                p_node_2->Set(SLAVE, true);
                p_node_2->Set(MASTER, false);
            }

            std::size_t id_cond = 0;
            std::vector<Condition::Pointer> slave_conds;
            for (std::size_t i = 0; i < NumberOfDivisions; i++) {
                id_cond++;
                const std::size_t ref_id = (2 * i)+1;
                Condition::Pointer pcond = r_contact_model_part.CreateNewCondition("SurfaceCondition3D4N", id_cond, {{ref_id, ref_id + 1, ref_id + 3, ref_id + 2}}, p_cond_prop);
                pcond->Set(SLAVE, true);
                pcond->Set(MASTER, false);
                slave_conds.push_back(pcond);
            }

            // Creating the base circle
            x = 0.0;
            std::size_t count = 0;
            const double dtheta = Angle/static_cast<double>(NumberOfDivisions);
            while (x < Lenght && count * dtheta < Globals::Pi/2.0) {
                x = Radius * std::sin(count * dtheta);
                y = Radius * (1.0 - std::cos(count * dtheta));
                id_node++;
                NodeType::Pointer p_node_1 = r_contact_model_part.CreateNewNode(id_node, x, y , 0.0);
                p_node_1->Set(SLAVE, false);
                p_node_1->Set(MASTER, true);
                id_node++;
                NodeType::Pointer p_node_2 = r_contact_model_part.CreateNewNode(id_node, x, y , 1.0);
                p_node_2->Set(SLAVE, false);
                p_node_2->Set(MASTER, true);
                count++;
            }

            // Adding master conditions
            IndexSet this_set;
            std::vector<Condition::Pointer> master_conds;
            for (std::size_t i = 0; i < count - 1; i++) {
                id_cond++;
                this_set.AddId(id_cond);
                const std::size_t ref_id = (2 * (i + NumberOfDivisions + 1)+1);
                Condition::Pointer pcond = r_contact_model_part.CreateNewCondition("SurfaceCondition3D4N", id_cond, {{ref_id +2, ref_id + 3, ref_id + 1, ref_id}}, p_cond_prop);
                pcond->Set(SLAVE, false);
                pcond->Set(MASTER, true);
                master_conds.push_back(pcond);
            }

            // We compute the normals
            MortarUtilities::ComputeNodesMeanNormalModelPart(r_contact_model_part);

            // We move mesh in order to test the dynamic search
            if (MoveMesh) {
                for (auto& inode : r_contact_model_part.Nodes()) {
                    if (inode.Is(MASTER)) {
                        inode.FastGetSolutionStepValue(DISPLACEMENT_X) = 0.1;
                        inode.Coordinates() += inode.FastGetSolutionStepValue(DISPLACEMENT);
                    }
                }
            }
        }

        /**
         * Checks the correct work of the search process
         * Test KDTree
         */
        KRATOS_TEST_CASE_IN_SUITE(SearchProcessKDTree, KratosContactStructuralMechanicsFastSuite)
        {
            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);
            r_model_part.CreateSubModelPart("ComputingContact");

            r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            r_model_part.AddNodalSolutionStepVariable(WEIGHTED_GAP);
            r_model_part.AddNodalSolutionStepVariable(VECTOR_LAGRANGE_MULTIPLIER);
            r_model_part.AddNodalSolutionStepVariable(WEIGHTED_VECTOR_RESIDUAL);
            r_model_part.AddNodalSolutionStepVariable(NODAL_H);
            r_model_part.AddNodalSolutionStepVariable(NORMAL);

            auto& process_info = r_model_part.GetProcessInfo();
            process_info[STEP] = 1;
            process_info[NL_ITERATION_NUMBER] = 1;
            process_info[DELTA_TIME] = 1.0;
            process_info[DISTANCE_THRESHOLD] = 1.0;

            // First we create the nodes
            const std::size_t number_of_divisions = 8;
            const double lenght = 4.0;
            const double radius = 6.0;
            const double angle = Globals::Pi/6;
            const double slope = 0.0;

            // We create our problem
            CreatePlaneCilynderProblemForSearch(r_model_part, number_of_divisions, lenght, radius, angle, slope);

            // Assign NodalH
            for (auto& r_node : r_model_part.Nodes()) {
                r_node.FastGetSolutionStepValue(NODAL_H) = 1.0;
            }

//             // DEBUG
//             GiDIOSearchDebug(r_model_part);

            // We compute the search
            Parameters this_parameters = Parameters(R"(
            {
                "simple_search"                        : true,
                "search_factor"                        : 3.5,
                "normal_orientation_threshold"         : 0.0,
                "type_search"                          : "InRadius",
                "check_gap"                            : "MappingCheck"
            })" );

            auto search_process = ContactSearchWrapperProcess(r_model_part, this_parameters);

            search_process.ExecuteInitialize();
            search_process.ExecuteInitializeSolutionStep();

//             // DEBUG
//             GiDIOSearchDebug(r_model_part);

//             // DEBUG
//             for (auto& r_cond : r_model_part.Conditions()) {
//                 if (r_cond.Is(SLAVE)) {
//                     auto p_indexes_pairs = r_cond.GetValue(INDEX_MAP);
//                     KRATOS_WATCH(r_cond.Id())
//                     KRATOS_WATCH(p_indexes_pairs->size())
//                     for (auto it_pair = p_indexes_pairs->begin(); it_pair != p_indexes_pairs->end(); ++it_pair ) {
//                         const IndexType master_id = p_indexes_pairs->GetId(it_pair);
//                         KRATOS_WATCH(master_id)
//                     }
//                 }
//             }

            // Check results
            std::size_t aux_index;
            auto& r_cond_1 = r_model_part.GetCondition(1);
            auto p_indexes_pairs_1 = r_cond_1.GetValue(INDEX_MAP);
            KRATOS_CHECK_EQUAL(p_indexes_pairs_1->size(), 2);
            aux_index = p_indexes_pairs_1->GetId(p_indexes_pairs_1->begin());
            KRATOS_CHECK(aux_index == 10 || aux_index == 9);
            auto it_pair_check = p_indexes_pairs_1->begin();
            it_pair_check++;
            aux_index = p_indexes_pairs_1->GetId(it_pair_check);
            KRATOS_CHECK(aux_index == 10 || aux_index == 9);

            auto& r_cond_2 = r_model_part.GetCondition(2);
            auto p_indexes_pairs_2 = r_cond_2.GetValue(INDEX_MAP);
            KRATOS_CHECK_EQUAL(p_indexes_pairs_2->size(), 2);
            aux_index = p_indexes_pairs_2->GetId(p_indexes_pairs_2->begin());
            KRATOS_CHECK(aux_index == 10 || aux_index == 11);
            it_pair_check = p_indexes_pairs_2->begin();
            it_pair_check++;
            aux_index = p_indexes_pairs_2->GetId(it_pair_check);
            KRATOS_CHECK(aux_index == 10 || aux_index == 11);

            auto& r_cond_3 = r_model_part.GetCondition(3);
            auto p_indexes_pairs_3 = r_cond_3.GetValue(INDEX_MAP);
            KRATOS_CHECK_EQUAL(p_indexes_pairs_3->size(), 2);
            aux_index = p_indexes_pairs_3->GetId(p_indexes_pairs_3->begin());
            KRATOS_CHECK(aux_index == 12 || aux_index == 11);
            it_pair_check = p_indexes_pairs_3->begin();
            it_pair_check++;
            aux_index = p_indexes_pairs_3->GetId(it_pair_check);
            KRATOS_CHECK(aux_index == 12 || aux_index == 11);

            auto& r_cond_4 = r_model_part.GetCondition(4);
            auto p_indexes_pairs_4 = r_cond_4.GetValue(INDEX_MAP);
            KRATOS_CHECK_EQUAL(p_indexes_pairs_4->size(), 3);
            aux_index = p_indexes_pairs_4->GetId(p_indexes_pairs_4->begin());
            KRATOS_CHECK(aux_index == 12 || aux_index == 13 || aux_index == 14);
            it_pair_check = p_indexes_pairs_4->begin();
            it_pair_check++;
            aux_index = p_indexes_pairs_4->GetId(it_pair_check);
            KRATOS_CHECK(aux_index == 12 || aux_index == 13 || aux_index == 14);
            it_pair_check++;
            aux_index = p_indexes_pairs_4->GetId(it_pair_check);
            KRATOS_CHECK(aux_index == 12 || aux_index == 13 || aux_index == 14);

            auto& r_cond_5 = r_model_part.GetCondition(5);
            auto p_indexes_pairs_5 = r_cond_5.GetValue(INDEX_MAP);
            KRATOS_CHECK_EQUAL(p_indexes_pairs_5->size(), 2);
            aux_index = p_indexes_pairs_5->GetId(p_indexes_pairs_5->begin());
            KRATOS_CHECK(aux_index == 14 || aux_index == 15);
            it_pair_check = p_indexes_pairs_5->begin();
            it_pair_check++;
            aux_index = p_indexes_pairs_5->GetId(it_pair_check);
            KRATOS_CHECK(aux_index == 14 || aux_index == 15);

            auto& r_cond_6 = r_model_part.GetCondition(6);
            auto p_indexes_pairs_6 = r_cond_6.GetValue(INDEX_MAP);
            KRATOS_CHECK_EQUAL(p_indexes_pairs_6->size(), 2);
            aux_index = p_indexes_pairs_6->GetId(p_indexes_pairs_6->begin());
            KRATOS_CHECK(aux_index == 16 || aux_index == 15);
            it_pair_check = p_indexes_pairs_6->begin();
            it_pair_check++;
            aux_index = p_indexes_pairs_6->GetId(it_pair_check);
            KRATOS_CHECK(aux_index == 16 || aux_index == 15);

            auto& r_cond_7 = r_model_part.GetCondition(7);
            auto p_indexes_pairs_7 = r_cond_7.GetValue(INDEX_MAP);
            KRATOS_CHECK_EQUAL(p_indexes_pairs_7->size(), 0);

            auto& r_cond_8 = r_model_part.GetCondition(8);
            auto p_indexes_pairs_8 = r_cond_8.GetValue(INDEX_MAP);
            KRATOS_CHECK_EQUAL(p_indexes_pairs_8->size(), 0);
        }

        /**
         * Checks the correct work of the search process
         * Test KDTree
         */
        KRATOS_TEST_CASE_IN_SUITE(SearchProcessKDTreeWithOBB, KratosContactStructuralMechanicsFastSuite)
        {
            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);
            r_model_part.CreateSubModelPart("ComputingContact");

            r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            r_model_part.AddNodalSolutionStepVariable(WEIGHTED_GAP);
            r_model_part.AddNodalSolutionStepVariable(VECTOR_LAGRANGE_MULTIPLIER);
            r_model_part.AddNodalSolutionStepVariable(WEIGHTED_VECTOR_RESIDUAL);
            r_model_part.AddNodalSolutionStepVariable(NODAL_H);
            r_model_part.AddNodalSolutionStepVariable(NORMAL);

            auto& process_info = r_model_part.GetProcessInfo();
            process_info[STEP] = 1;
            process_info[NL_ITERATION_NUMBER] = 1;
            process_info[DELTA_TIME] = 1.0;
            process_info[DISTANCE_THRESHOLD] = 1.0;

            // First we create the nodes
            const std::size_t number_of_divisions = 8;
            const double lenght = 4.0;
            const double radius = 6.0;
            const double angle = Globals::Pi/6;
            const double slope = 0.0;

            // We create our problem
            CreatePlaneCilynderProblemForSearch(r_model_part, number_of_divisions, lenght, radius, angle, slope);

            // Assign NodalH
            for (auto& r_node : r_model_part.Nodes()) {
                r_node.FastGetSolutionStepValue(NODAL_H) = 1.0;
            }

//             // DEBUG
//             GiDIOSearchDebug(r_model_part);

            // We compute the search
            Parameters this_parameters = Parameters(R"(
            {
                "simple_search"                        : true,
                "search_factor"                        : 3.5,
                "normal_orientation_threshold"         : 0.0,
                "type_search"                          : "InRadiusWithOBB",
                "check_gap"                            : "MappingCheck",
                "octree_search_parameters" : {
                    "bounding_box_factor"    : 0.1,
                    "debug_obb"              : false
                    }
            })" );

            auto search_process = ContactSearchWrapperProcess(r_model_part, this_parameters);

            search_process.ExecuteInitialize();
            search_process.ExecuteInitializeSolutionStep();

//             // DEBUG
//             GiDIOSearchDebug(r_model_part);

//             // DEBUG
//             for (auto& r_cond : r_model_part.Conditions()) {
//                 if (r_cond.Is(SLAVE)) {
//                     auto p_indexes_pairs = r_cond.GetValue(INDEX_MAP);
//                     KRATOS_WATCH(r_cond.Id())
//                     KRATOS_WATCH(p_indexes_pairs->size())
//                     for (auto it_pair = p_indexes_pairs->begin(); it_pair != p_indexes_pairs->end(); ++it_pair ) {
//                         const IndexType master_id = p_indexes_pairs->GetId(it_pair);
//                         KRATOS_WATCH(master_id)
//                     }
//                 }
//             }

            // Check results
            std::size_t aux_index;
            auto& r_cond_1 = r_model_part.GetCondition(1);
            auto p_indexes_pairs_1 = r_cond_1.GetValue(INDEX_MAP);
            KRATOS_CHECK_EQUAL(p_indexes_pairs_1->size(), 2);
            aux_index = p_indexes_pairs_1->GetId(p_indexes_pairs_1->begin());
            KRATOS_CHECK(aux_index == 10 || aux_index == 9);
            auto it_pair_check = p_indexes_pairs_1->begin();
            it_pair_check++;
            aux_index = p_indexes_pairs_1->GetId(it_pair_check);
            KRATOS_CHECK(aux_index == 10 || aux_index == 9);

            auto& r_cond_2 = r_model_part.GetCondition(2);
            auto p_indexes_pairs_2 = r_cond_2.GetValue(INDEX_MAP);
            KRATOS_CHECK_EQUAL(p_indexes_pairs_2->size(), 2);
            aux_index = p_indexes_pairs_2->GetId(p_indexes_pairs_2->begin());
            KRATOS_CHECK(aux_index == 10 || aux_index == 11);
            it_pair_check = p_indexes_pairs_2->begin();
            it_pair_check++;
            aux_index = p_indexes_pairs_2->GetId(it_pair_check);
            KRATOS_CHECK(aux_index == 10 || aux_index == 11);

            auto& r_cond_3 = r_model_part.GetCondition(3);
            auto p_indexes_pairs_3 = r_cond_3.GetValue(INDEX_MAP);
            KRATOS_CHECK_EQUAL(p_indexes_pairs_3->size(), 2);
            aux_index = p_indexes_pairs_3->GetId(p_indexes_pairs_3->begin());
            KRATOS_CHECK(aux_index == 12 || aux_index == 11);
            it_pair_check = p_indexes_pairs_3->begin();
            it_pair_check++;
            aux_index = p_indexes_pairs_3->GetId(it_pair_check);
            KRATOS_CHECK(aux_index == 12 || aux_index == 11);

            auto& r_cond_4 = r_model_part.GetCondition(4);
            auto p_indexes_pairs_4 = r_cond_4.GetValue(INDEX_MAP);
            KRATOS_CHECK_EQUAL(p_indexes_pairs_4->size(), 3);
            aux_index = p_indexes_pairs_4->GetId(p_indexes_pairs_4->begin());
            KRATOS_CHECK(aux_index == 12 || aux_index == 13 || aux_index == 14);
            it_pair_check = p_indexes_pairs_4->begin();
            it_pair_check++;
            aux_index = p_indexes_pairs_4->GetId(it_pair_check);
            KRATOS_CHECK(aux_index == 12 || aux_index == 13 || aux_index == 14);
            it_pair_check++;
            aux_index = p_indexes_pairs_4->GetId(it_pair_check);
            KRATOS_CHECK(aux_index == 12 || aux_index == 13 || aux_index == 14);

            auto& r_cond_5 = r_model_part.GetCondition(5);
            auto p_indexes_pairs_5 = r_cond_5.GetValue(INDEX_MAP);
            KRATOS_CHECK_EQUAL(p_indexes_pairs_5->size(), 1);
            KRATOS_CHECK_EQUAL(p_indexes_pairs_5->GetId(p_indexes_pairs_5->begin()), 14);

            auto& r_cond_6 = r_model_part.GetCondition(6);
            auto p_indexes_pairs_6 = r_cond_6.GetValue(INDEX_MAP);
            KRATOS_CHECK_EQUAL(p_indexes_pairs_6->size(), 0);

            auto& r_cond_7 = r_model_part.GetCondition(7);
            auto p_indexes_pairs_7 = r_cond_7.GetValue(INDEX_MAP);
            KRATOS_CHECK_EQUAL(p_indexes_pairs_7->size(), 0);

            auto& r_cond_8 = r_model_part.GetCondition(8);
            auto p_indexes_pairs_8 = r_cond_8.GetValue(INDEX_MAP);
            KRATOS_CHECK_EQUAL(p_indexes_pairs_8->size(), 0);
        }

        /**
         * Checks the correct work of the search process
         * Test Octree
         */
        KRATOS_TEST_CASE_IN_SUITE(SearchProcessOctree, KratosContactStructuralMechanicsFastSuite)
        {
            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);
            r_model_part.CreateSubModelPart("ComputingContact");

            r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            r_model_part.AddNodalSolutionStepVariable(WEIGHTED_GAP);
            r_model_part.AddNodalSolutionStepVariable(VECTOR_LAGRANGE_MULTIPLIER);
            r_model_part.AddNodalSolutionStepVariable(WEIGHTED_VECTOR_RESIDUAL);
            r_model_part.AddNodalSolutionStepVariable(NODAL_H);
            r_model_part.AddNodalSolutionStepVariable(NORMAL);

            auto& process_info = r_model_part.GetProcessInfo();
            process_info[STEP] = 1;
            process_info[NL_ITERATION_NUMBER] = 1;
            process_info[DELTA_TIME] = 1.0;
            process_info[DISTANCE_THRESHOLD] = 1.0;

            // First we create the nodes
            const std::size_t number_of_divisions = 8;
            const double lenght = 4.0;
            const double radius = 6.0;
            const double angle = Globals::Pi/6;
            const double slope = 0.0;

            // We create our problem
            CreatePlaneCilynderProblemForSearch(r_model_part, number_of_divisions, lenght, radius, angle, slope);

            // Assign NodalH
            for (auto& r_node : r_model_part.Nodes()) {
                r_node.FastGetSolutionStepValue(NODAL_H) = 1.0;
            }

//             // DEBUG
//             GiDIOSearchDebug(r_model_part);

            // We compute the search
            Parameters this_parameters = Parameters(R"(
            {
                "simple_search"                        : true,
                "search_factor"                        : 3.5,
                "normal_orientation_threshold"         : 0.0,
                "type_search"                          : "OctreeWithOBB",
                "check_gap"                            : "MappingCheck",
                "octree_search_parameters" : {
                    "bounding_box_factor"    : 0.1,
                    "debug_obb"              : false,
                    "OBB_intersection_type"  : "SeparatingAxisTheorem"
                    }
            })" );

            auto search_process = ContactSearchWrapperProcess(r_model_part, this_parameters);

            search_process.ExecuteInitialize();
            search_process.ExecuteInitializeSolutionStep();

//             // DEBUG
//             GiDIOSearchDebug(r_model_part);

//             // DEBUG
//             for (auto& r_cond : r_model_part.Conditions()) {
//                 if (r_cond.Is(SLAVE)) {
//                     auto p_indexes_pairs = r_cond.GetValue(INDEX_MAP);
//                     KRATOS_WATCH(r_cond.Id())
//                     KRATOS_WATCH(p_indexes_pairs->size())
//                     for (auto it_pair = p_indexes_pairs->begin(); it_pair != p_indexes_pairs->end(); ++it_pair ) {
//                         const IndexType master_id = p_indexes_pairs->GetId(it_pair);
//                         KRATOS_WATCH(master_id)
//                     }
//                 }
//             }

            // Check results
            std::size_t aux_index;
            auto& r_cond_1 = r_model_part.GetCondition(1);
            auto p_indexes_pairs_1 = r_cond_1.GetValue(INDEX_MAP);
            KRATOS_CHECK_EQUAL(p_indexes_pairs_1->size(), 2);
            aux_index = p_indexes_pairs_1->GetId(p_indexes_pairs_1->begin());
            KRATOS_CHECK(aux_index == 10 || aux_index == 9);
            auto it_pair_check = p_indexes_pairs_1->begin();
            it_pair_check++;
            aux_index = p_indexes_pairs_1->GetId(it_pair_check);
            KRATOS_CHECK(aux_index == 10 || aux_index == 9);

            auto& r_cond_2 = r_model_part.GetCondition(2);
            auto p_indexes_pairs_2 = r_cond_2.GetValue(INDEX_MAP);
            KRATOS_CHECK_EQUAL(p_indexes_pairs_2->size(), 2);
            aux_index = p_indexes_pairs_2->GetId(p_indexes_pairs_2->begin());
            KRATOS_CHECK(aux_index == 10 || aux_index == 11);
            it_pair_check = p_indexes_pairs_2->begin();
            it_pair_check++;
            aux_index = p_indexes_pairs_2->GetId(it_pair_check);
            KRATOS_CHECK(aux_index == 10 || aux_index == 11);

            auto& r_cond_3 = r_model_part.GetCondition(3);
            auto p_indexes_pairs_3 = r_cond_3.GetValue(INDEX_MAP);
            KRATOS_CHECK_EQUAL(p_indexes_pairs_3->size(), 2);
            aux_index = p_indexes_pairs_3->GetId(p_indexes_pairs_3->begin());
            KRATOS_CHECK(aux_index == 12 || aux_index == 11);
            it_pair_check = p_indexes_pairs_3->begin();
            it_pair_check++;
            aux_index = p_indexes_pairs_3->GetId(it_pair_check);
            KRATOS_CHECK(aux_index == 12 || aux_index == 11);

            auto& r_cond_4 = r_model_part.GetCondition(4);
            auto p_indexes_pairs_4 = r_cond_4.GetValue(INDEX_MAP);
            KRATOS_CHECK_EQUAL(p_indexes_pairs_4->size(), 3);
            aux_index = p_indexes_pairs_4->GetId(p_indexes_pairs_4->begin());
            KRATOS_CHECK(aux_index == 12 || aux_index == 13 || aux_index == 14);
            it_pair_check = p_indexes_pairs_4->begin();
            it_pair_check++;
            aux_index = p_indexes_pairs_4->GetId(it_pair_check);
            KRATOS_CHECK(aux_index == 12 || aux_index == 13 || aux_index == 14);
            it_pair_check++;
            aux_index = p_indexes_pairs_4->GetId(it_pair_check);
            KRATOS_CHECK(aux_index == 12 || aux_index == 13 || aux_index == 14);

            auto& r_cond_5 = r_model_part.GetCondition(5);
            auto p_indexes_pairs_5 = r_cond_5.GetValue(INDEX_MAP);
            KRATOS_CHECK_EQUAL(p_indexes_pairs_5->size(), 1);
            KRATOS_CHECK_EQUAL(p_indexes_pairs_5->GetId(p_indexes_pairs_5->begin()), 14);

            auto& r_cond_6 = r_model_part.GetCondition(6);
            auto p_indexes_pairs_6 = r_cond_6.GetValue(INDEX_MAP);
            KRATOS_CHECK_EQUAL(p_indexes_pairs_6->size(), 0);

            auto& r_cond_7 = r_model_part.GetCondition(7);
            auto p_indexes_pairs_7 = r_cond_7.GetValue(INDEX_MAP);
            KRATOS_CHECK_EQUAL(p_indexes_pairs_7->size(), 0);

            auto& r_cond_8 = r_model_part.GetCondition(8);
            auto p_indexes_pairs_8 = r_cond_8.GetValue(INDEX_MAP);
            KRATOS_CHECK_EQUAL(p_indexes_pairs_8->size(), 0);
        }

    } // namespace Testing
}  // namespace Kratos.
