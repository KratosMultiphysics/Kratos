// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:	   BSD License
//				   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:  Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "contact_structural_mechanics_application_variables.h"
// #include "includes/gid_io.h"
#include "utilities/variable_utils.h"
#include "processes/simple_mortar_mapper_process.h"

namespace Kratos
{
    namespace Testing
    {
        typedef Node<3> NodeType;

//         void GiDIOGapDebug(ModelPart& rModelPart)
//         {
//             GidIO<> gid_io("TEST_WEIGHTED_GAP", GiD_PostBinary, SingleFile, WriteUndeformed,  WriteConditionsOnly);
//             const int nl_iter = rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER];
//             const double label = static_cast<double>(nl_iter);
//
//             gid_io.InitializeMesh(label);
//             gid_io.WriteMesh(rModelPart.GetMesh());
//             gid_io.FinalizeMesh();
//             gid_io.InitializeResults(label, rModelPart.GetMesh());
//             gid_io.WriteNodalFlags(ACTIVE, "ACTIVE", rModelPart.Nodes(), label);
//             gid_io.WriteNodalFlags(SLAVE, "SLAVE", rModelPart.Nodes(), label);
//             gid_io.WriteNodalFlags(ISOLATED, "ISOLATED", rModelPart.Nodes(), label);
//             gid_io.WriteNodalResults(DISPLACEMENT, rModelPart.Nodes(), label, 0);
//             gid_io.WriteNodalResults(WEIGHTED_GAP, rModelPart.Nodes(), label, 0);
//             gid_io.WriteNodalResults(WEIGHTED_SLIP, rModelPart.Nodes(), label, 0);
//             gid_io.WriteNodalResultsNonHistorical(NORMAL_GAP, rModelPart.Nodes(), label);
//             gid_io.WriteNodalResultsNonHistorical(NODAL_AREA, rModelPart.Nodes(), label);
//             gid_io.WriteNodalResultsNonHistorical(AUXILIAR_COORDINATES, rModelPart.Nodes(), label);
//             gid_io.WriteNodalResults(NORMAL, rModelPart.Nodes(), label, 0);
//         }

        /**
         * This method can be used to create a plane/cylinder condition set
         */
        void CreateNewProblem(
            ModelPart& rModelPart,
            const std::size_t NumberOfDivisions,
            const double Lenght,
            const double Radius = 0.0,
            const double Angle = 0.0,
            const double Slope = 0.0,
            const double MoveMesh = 0.0
            )
        {
            rModelPart.CreateSubModelPart("SlaveModelPart");
            ModelPart& r_slave_model_part = rModelPart.GetSubModelPart("SlaveModelPart");
            rModelPart.CreateSubModelPart("MasterModelPart");
            ModelPart& r_master_model_part = rModelPart.GetSubModelPart("MasterModelPart");

            Properties::Pointer p_cond_prop = rModelPart.CreateNewProperties(0);

            double x, y;

            // Creating the base geometry
            std::size_t id_node = 0;
            double dx = Lenght/static_cast<double>(NumberOfDivisions);
            for (std::size_t i = 0; i < NumberOfDivisions + 1; ++i) {
                x = dx * i;
                y = Slope * dx * i;
                id_node++;
                NodeType::Pointer p_node_1 = rModelPart.CreateNewNode(id_node, x , y , 0.0);
                r_slave_model_part.AddNode(p_node_1);
                p_node_1->Set(SLAVE, true);
                p_node_1->Set(MASTER, false);
                p_node_1->Set(ACTIVE, true);
                id_node++;
                NodeType::Pointer p_node_2 = rModelPart.CreateNewNode(id_node, x , y , 1.0);
                r_slave_model_part.AddNode(p_node_2);
                p_node_2->Set(SLAVE, true);
                p_node_2->Set(MASTER, false);
                p_node_2->Set(ACTIVE, true);
            }

            std::size_t r_cond = 0;
            std::vector<Condition::Pointer> slave_conds;
            for (std::size_t i = 0; i < NumberOfDivisions; i++) {
                r_cond++;
                const std::size_t ref_id = (2 * i)+1;
                Condition::Pointer p_cond = rModelPart.CreateNewCondition("Condition3D4N", r_cond, {{ref_id, ref_id + 1, ref_id + 3, ref_id + 2}}, p_cond_prop);
                r_slave_model_part.AddCondition(p_cond);
                p_cond->Set(SLAVE, true);
                p_cond->Set(MASTER, false);
                slave_conds.push_back(p_cond);
            }

            // Values
            x = 0.0;
            std::size_t count = 0;
            if (Radius > 0.0) {
                // Creating the base circle
                const double dtheta = Angle/static_cast<double>(NumberOfDivisions);
                while (x < Lenght && count * dtheta < Globals::Pi/2.0) {
                    x = Radius * std::sin(count * dtheta);
                    y = Radius * (1.0 - std::cos(count * dtheta));
                    id_node++;
                    NodeType::Pointer p_node_1 = rModelPart.CreateNewNode(id_node, x, y , 0.0);
                    r_master_model_part.AddNode(p_node_1);
                    p_node_1->Set(SLAVE, false);
                    p_node_1->Set(MASTER, true);
                    p_node_1->Set(ACTIVE, false);
                    id_node++;
                    NodeType::Pointer p_node_2 = rModelPart.CreateNewNode(id_node, x, y , 1.0);
                    r_master_model_part.AddNode(p_node_2);
                    p_node_2->Set(SLAVE, false);
                    p_node_2->Set(MASTER, true);
                    p_node_2->Set(ACTIVE, false);
                    count++;
                }
            } else {
                // Creating the second plane
                dx = 0.25 * Lenght/static_cast<double>(NumberOfDivisions);
                for (std::size_t i = 0; i < NumberOfDivisions + 1; ++i) {
                    x = dx * i;
                    y = Slope * dx * i + 1.0e-3;
                    id_node++;
                    NodeType::Pointer p_node_1 = rModelPart.CreateNewNode(id_node, x, y, 0.0);
                    r_master_model_part.AddNode(p_node_1);
                    p_node_1->Set(SLAVE, false);
                    p_node_1->Set(MASTER, true);
                    p_node_1->Set(ACTIVE, false);
                    id_node++;
                    NodeType::Pointer p_node_2 = rModelPart.CreateNewNode(id_node, x, y, 1.0);
                    r_master_model_part.AddNode(p_node_2);
                    p_node_2->Set(SLAVE, false);
                    p_node_2->Set(MASTER, true);
                    p_node_2->Set(ACTIVE, false);
                    count++;
                }
            }

            // Adding map
            IndexSet this_set;
            std::vector<Condition::Pointer> master_conds;
            for (std::size_t i = 0; i < count - 1; i++) {
                r_cond++;
                this_set.AddId(r_cond);
                const std::size_t ref_id = (2 * (i + NumberOfDivisions + 1)+1);
                Condition::Pointer p_cond = rModelPart.CreateNewCondition("Condition3D4N", r_cond, {{ref_id +2, ref_id + 3, ref_id + 1, ref_id}}, p_cond_prop);
                r_master_model_part.AddCondition(p_cond);
                p_cond->Set(SLAVE, false);
                p_cond->Set(MASTER, true);
                master_conds.push_back(p_cond);
            }

            // We compute the normals
            MortarUtilities::ComputeNodesMeanNormalModelPart(rModelPart);

            // We compute the normal gap to compare with the weighted gap
            // We add the index SetScalarVar
            for(auto& r_cond : rModelPart.Conditions()) {
                if (r_cond.Is(SLAVE))
                    r_cond.SetValue(INDEX_SET, Kratos::make_shared<IndexSet>(this_set));
            }

            // We set the auxiliar Coordinates
            for(auto& r_node : rModelPart.Nodes()) {
                if (r_node.Is(MASTER))
                    r_node.SetValue(AUXILIAR_COORDINATES, r_node.Coordinates());
                else
                    r_node.SetValue(AUXILIAR_COORDINATES, ZeroVector(3));
            }

            // We set the mapper parameters
            Parameters mapping_parameters = Parameters(R"({"distance_threshold" : 1.0e24, "origin_variable_historical" : false,
        "destination_variable_historical" : false})" );
            mapping_parameters["distance_threshold"].SetDouble(rModelPart.GetProcessInfo()[DISTANCE_THRESHOLD]);
            typedef SimpleMortarMapperProcess<3, 4, Variable<array_1d<double, 3>>> MapperType;
            MapperType mapper = MapperType(r_master_model_part, r_slave_model_part, AUXILIAR_COORDINATES, mapping_parameters);
            mapper.Execute();

            // We compute now the normal gap and set the nodes under certain threshold as active
            for(auto& r_node : rModelPart.Nodes()) {
                if (r_node.Is(SLAVE)) {
                    // We compute the gap
                    const array_1d<double, 3>& r_normal = r_node.FastGetSolutionStepValue(NORMAL);
                    const array_1d<double, 3>& r_auxiliar_coordinates = r_node.GetValue(AUXILIAR_COORDINATES);
                    const array_1d<double, 3>& r_components_gap = ( r_node.Coordinates() - r_auxiliar_coordinates);
                    const double gap = inner_prod(r_components_gap, - r_normal);
                    r_node.SetValue(NORMAL_GAP, gap);
                } else
                    r_node.SetValue(NORMAL_GAP, 0.0);
            }

            // We set the database
            auto& r_process_info = rModelPart.GetProcessInfo();
            ModelPart& r_computing_contact_model_part = rModelPart.GetSubModelPart("ComputingContact");
            for (auto& r_slave_cond : slave_conds) {
                for (auto& r_master_cond : master_conds) {
                    r_cond++;
                    Condition::Pointer p_auxiliar_condition = r_computing_contact_model_part.CreateNewCondition("ALMFrictionalMortarContactCondition3D4N", r_cond, r_slave_cond->GetGeometry(), p_cond_prop);
                    // We set the geometrical values
                    p_auxiliar_condition->SetValue(PAIRED_GEOMETRY, r_master_cond->pGetGeometry());
                    p_auxiliar_condition->SetValue(NORMAL, r_slave_cond->GetValue(NORMAL));
                    p_auxiliar_condition->SetValue(PAIRED_NORMAL, r_master_cond->GetValue(NORMAL));
                    // We activate the condition and initialize it
                    p_auxiliar_condition->Set(ACTIVE, true);
                    p_auxiliar_condition->Initialize();
                    p_auxiliar_condition->InitializeSolutionStep(r_process_info);
                }
            }

            // We move mesh in order to test the slip
            if (MoveMesh > 0.0) {
                for (auto& r_node : rModelPart.Nodes()) {
                    if (r_node.Is(MASTER)) {
                        r_node.FastGetSolutionStepValue(DISPLACEMENT_X) = MoveMesh;
                        r_node.Coordinates() += r_node.FastGetSolutionStepValue(DISPLACEMENT);
                    }
                }

                // We set the auxiliar Coordinates
                for(auto& r_node : rModelPart.Nodes()) {
                    if (r_node.Is(MASTER))
                        r_node.SetValue(AUXILIAR_COORDINATES, r_node.Coordinates());
                    else
                        r_node.SetValue(AUXILIAR_COORDINATES, ZeroVector(3));
                }

                mapper.Execute();

                // We compute now the normal gap and set the nodes under certain threshold as active
                for(auto& r_node : rModelPart.Nodes()) {
                    if (r_node.Is(SLAVE)) {
                        // We compute the gap
                        const array_1d<double, 3>& r_normal = r_node.FastGetSolutionStepValue(NORMAL);
                        const array_1d<double, 3>& r_auxiliar_coordinates = r_node.GetValue(AUXILIAR_COORDINATES);
                        const array_1d<double, 3>& r_components_gap = ( r_node.Coordinates() - r_auxiliar_coordinates);
                        const double gap = inner_prod(r_components_gap, - r_normal);
                        r_node.SetValue(NORMAL_GAP, gap);
                    } else
                        r_node.SetValue(NORMAL_GAP, 0.0);
                }
            }
        }

        /**
        * Checks the correct work of the weighted gap computation
        * Test 1
        */
        KRATOS_TEST_CASE_IN_SUITE(WeightedGap1, KratosContactStructuralMechanicsFastSuite)
        {
            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Contact", 3);
            r_model_part.CreateSubModelPart("ComputingContact");

            r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            r_model_part.AddNodalSolutionStepVariable(WEIGHTED_GAP);
            r_model_part.AddNodalSolutionStepVariable(WEIGHTED_SLIP);
            r_model_part.AddNodalSolutionStepVariable(VECTOR_LAGRANGE_MULTIPLIER);
            r_model_part.AddNodalSolutionStepVariable(NORMAL);

            auto& r_process_info = r_model_part.GetProcessInfo();
            r_process_info[STEP] = 1;
            r_process_info[NL_ITERATION_NUMBER] = 1;
            r_process_info[DELTA_TIME] = 1.0;
            r_process_info[DISTANCE_THRESHOLD] = 1.0;

            // First we create the nodes
            const std::size_t number_of_divisions = 8;
            const double lenght = 4.0;
            const double radius = 6.0;
            const double angle = Globals::Pi/6;
            const double slope = 0.0;

            // We create our problem
            CreateNewProblem(r_model_part, number_of_divisions, lenght, radius, angle, slope);

            // We compute the explicit contribution
            VariableUtils().SetScalarVar<Variable<double>>(WEIGHTED_GAP, 0.0, r_model_part.Nodes());
            for (auto& r_cond : r_model_part.GetSubModelPart("ComputingContact").Conditions())
                r_cond.AddExplicitContribution(r_process_info);

            // DEBUG
//             GiDIOGapDebug(r_model_part);

            const double tolerance = 1.0e-4;
            for (auto& r_node : r_model_part.Nodes()) {
                if (r_node.Is(SLAVE)) {
                    if (std::abs(r_node.FastGetSolutionStepValue(WEIGHTED_GAP)) > 0.0) {
                        const double normal_gap = r_node.GetValue(NORMAL_GAP);
                        const double weighted_gap_corrected = r_node.FastGetSolutionStepValue(WEIGHTED_GAP)/r_node.GetValue(NODAL_AREA);
                        KRATOS_CHECK_LESS_EQUAL(std::abs(weighted_gap_corrected - normal_gap)/std::abs(normal_gap), tolerance);
                    }
                }
            }
        }

        /**
        * Checks the correct work of the weighted gap computation
        * Test 2
        */
        KRATOS_TEST_CASE_IN_SUITE(WeightedGap2, KratosContactStructuralMechanicsFastSuite)
        {
            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Contact", 3);
            r_model_part.CreateSubModelPart("ComputingContact");

            r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            r_model_part.AddNodalSolutionStepVariable(WEIGHTED_GAP);
            r_model_part.AddNodalSolutionStepVariable(WEIGHTED_SLIP);
            r_model_part.AddNodalSolutionStepVariable(VECTOR_LAGRANGE_MULTIPLIER);
            r_model_part.AddNodalSolutionStepVariable(NORMAL);

            auto& r_process_info = r_model_part.GetProcessInfo();
            r_process_info[STEP] = 1;
            r_process_info[NL_ITERATION_NUMBER] = 1;
            r_process_info[DELTA_TIME] = 1.0;
            r_process_info[DISTANCE_THRESHOLD] = 1.0;

            // First we create the nodes
            const std::size_t number_of_divisions = 8;
            const double lenght = 4.0;
            const double radius = 6.0;
            const double angle = Globals::Pi/6;
            const double slope = 0.0;

            // We create our problem
            CreateNewProblem(r_model_part, number_of_divisions, lenght, radius, angle, slope);

            // We compute the explicit contribution
            const array_1d<double, 3> zero_vector = ZeroVector(3);;
            VariableUtils().SetScalarVar<Variable<double>>(WEIGHTED_GAP, 0.0, r_model_part.Nodes());
            VariableUtils().SetVectorVar(WEIGHTED_SLIP, zero_vector, r_model_part.Nodes());
            for (auto& r_cond : r_model_part.GetSubModelPart("ComputingContact").Conditions())
                r_cond.AddExplicitContribution(r_process_info);

            // DEBUG
//             GiDIOGapDebug(r_model_part);

            const double tolerance = 1.0e-4;
            for (auto& r_node : r_model_part.Nodes()) {
                if (r_node.Is(SLAVE)) {
                    if (std::abs(r_node.FastGetSolutionStepValue(WEIGHTED_GAP)) > 0.0) {
                        const double normal_gap = r_node.GetValue(NORMAL_GAP);
                        const double weighted_gap_corrected = r_node.FastGetSolutionStepValue(WEIGHTED_GAP)/r_node.GetValue(NODAL_AREA);
                        KRATOS_CHECK_LESS_EQUAL(std::abs(weighted_gap_corrected - normal_gap)/std::abs(normal_gap), tolerance);
                    }
                    KRATOS_CHECK_LESS_EQUAL(norm_2(r_node.FastGetSolutionStepValue(WEIGHTED_SLIP)), tolerance);
                }
            }
        }

        /**
        * Checks the correct work of the weighted gap computation (with objective slip)
        * Test 3
        */
        KRATOS_TEST_CASE_IN_SUITE(WeightedGap3, KratosContactStructuralMechanicsFastSuite2)
        {
            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Contact", 3);
            r_model_part.CreateSubModelPart("ComputingContact");

            r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            r_model_part.AddNodalSolutionStepVariable(WEIGHTED_GAP);
            r_model_part.AddNodalSolutionStepVariable(WEIGHTED_SLIP);
            r_model_part.AddNodalSolutionStepVariable(VECTOR_LAGRANGE_MULTIPLIER);
            r_model_part.AddNodalSolutionStepVariable(NORMAL);

            auto& r_process_info = r_model_part.GetProcessInfo();
            r_process_info[STEP] = 1;
            r_process_info[NL_ITERATION_NUMBER] = 1;
            r_process_info[DELTA_TIME] = 1.0;
            r_process_info[DISTANCE_THRESHOLD] = 1.0;

            // First we create the nodes
            const std::size_t number_of_divisions = 8;
            const double lenght = 4.0;
            const double radius = 0.0;
            const double angle = 0.0;
            const double slope = 0.0;

            // We create our problem
            const double delta_x = 0.1;
            CreateNewProblem(r_model_part, number_of_divisions, lenght, radius, angle, slope, delta_x);

            // We compute the explicit contribution
            const array_1d<double, 3> zero_vector = ZeroVector(3);;
            VariableUtils().SetScalarVar<Variable<double>>(WEIGHTED_GAP, 0.0, r_model_part.Nodes());
            VariableUtils().SetVectorVar(WEIGHTED_SLIP, zero_vector, r_model_part.Nodes());
            for (auto& r_cond : r_model_part.GetSubModelPart("ComputingContact").Conditions()) {
                r_cond.AddExplicitContribution(r_process_info);
                r_cond.FinalizeSolutionStep(r_process_info);
            }

//             // DEBUG
//             GiDIOGapDebug(r_model_part);

            const double tolerance = 1.0e-4;
            array_1d<double, 3> slip = ZeroVector(3);
            slip[0] = delta_x;
            std::size_t counter = 0;
            for (auto& r_node : r_model_part.Nodes()) {
                if (r_node.Is(SLAVE)) {
                    const double normal_gap = r_node.GetValue(NORMAL_GAP);
                    const double weighted_gap_corrected = r_node.FastGetSolutionStepValue(WEIGHTED_GAP)/r_node.GetValue(NODAL_AREA);
                    if (std::abs(weighted_gap_corrected - normal_gap)/std::abs(normal_gap) < tolerance) {
                        if (norm_2(r_node.FastGetSolutionStepValue(WEIGHTED_SLIP)) > 0.0) {
                            const array_1d<double, 3> weighted_slip_corrected = r_node.FastGetSolutionStepValue(WEIGHTED_SLIP)/r_node.GetValue(NODAL_AREA);
                            KRATOS_WATCH(slip)
                            KRATOS_WATCH(weighted_slip_corrected)
//                             KRATOS_CHECK_LESS_EQUAL(norm_2(weighted_slip_corrected - slip)/norm_2(slip), tolerance);
                            ++counter;
                        }
                    }
                }
            }
            KRATOS_CHECK(counter != 0);
        }

        /**
        * Checks the correct work of the weighted gap computation (with objective slip)
        * Test 4
        */
        KRATOS_TEST_CASE_IN_SUITE(WeightedGap4, KratosContactStructuralMechanicsFastSuite)
        {
            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Contact", 3);
            r_model_part.CreateSubModelPart("ComputingContact");

            r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            r_model_part.AddNodalSolutionStepVariable(WEIGHTED_GAP);
            r_model_part.AddNodalSolutionStepVariable(WEIGHTED_SLIP);
            r_model_part.AddNodalSolutionStepVariable(VECTOR_LAGRANGE_MULTIPLIER);
            r_model_part.AddNodalSolutionStepVariable(NORMAL);

            auto& r_process_info = r_model_part.GetProcessInfo();
            r_process_info[STEP] = 1;
            r_process_info[NL_ITERATION_NUMBER] = 1;
            r_process_info[DELTA_TIME] = 1.0;
            r_process_info[DISTANCE_THRESHOLD] = 1.0;

            // First we create the nodes
            const std::size_t number_of_divisions = 8;
            const double lenght = 4.0;
            const double radius = 0.0;
            const double angle = 0.0;
            const double slope = 0.0;

            // We create our problem
            const double delta_x = 0.1;
            CreateNewProblem(r_model_part, number_of_divisions, lenght, radius, angle, slope, delta_x);

            // We compute the explicit contribution
            const array_1d<double, 3> zero_vector = ZeroVector(3);;
            VariableUtils().SetScalarVar<Variable<double>>(WEIGHTED_GAP, 0.0, r_model_part.Nodes());
            VariableUtils().SetVectorVar(WEIGHTED_SLIP, zero_vector, r_model_part.Nodes());
            for (auto& r_cond : r_model_part.GetSubModelPart("ComputingContact").Conditions()) {
                r_cond.Set(MODIFIED, true);
                r_cond.AddExplicitContribution(r_process_info);
                r_cond.FinalizeSolutionStep(r_process_info);
            }

//             // DEBUG
//             GiDIOGapDebug(r_model_part);

            const double tolerance = 1.0e-4;
            array_1d<double, 3> slip = ZeroVector(3);
            slip[0] = delta_x;
            std::size_t counter = 0;
            for (auto& r_node : r_model_part.Nodes()) {
                if (r_node.Is(SLAVE)) {
                    const double normal_gap = r_node.GetValue(NORMAL_GAP);
                    const double weighted_gap_corrected = r_node.FastGetSolutionStepValue(WEIGHTED_GAP)/r_node.GetValue(NODAL_AREA);
                    if (std::abs(weighted_gap_corrected - normal_gap)/std::abs(normal_gap) < tolerance) {
                        if (norm_2(r_node.FastGetSolutionStepValue(WEIGHTED_SLIP)) > 0.0) {
                            const array_1d<double, 3> weighted_slip_corrected = r_node.FastGetSolutionStepValue(WEIGHTED_SLIP)/r_node.GetValue(NODAL_AREA);
                            KRATOS_CHECK_LESS_EQUAL(norm_2(weighted_slip_corrected - slip)/norm_2(slip), tolerance);
                            ++counter;
                        }
                    }
                }
            }
            KRATOS_CHECK(counter != 0);
        }

        /**
        * Checks the correct work of the weighted gap computation (with non-objective slip)
        * Test 5
        */
        KRATOS_TEST_CASE_IN_SUITE(WeightedGap5, KratosContactStructuralMechanicsFastSuite)
        {
            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Contact", 3);
            r_model_part.CreateSubModelPart("ComputingContact");

            r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            r_model_part.AddNodalSolutionStepVariable(WEIGHTED_GAP);
            r_model_part.AddNodalSolutionStepVariable(WEIGHTED_SLIP);
            r_model_part.AddNodalSolutionStepVariable(VECTOR_LAGRANGE_MULTIPLIER);
            r_model_part.AddNodalSolutionStepVariable(NORMAL);

            auto& r_process_info = r_model_part.GetProcessInfo();
            r_process_info[STEP] = 1;
            r_process_info[NL_ITERATION_NUMBER] = 1;
            r_process_info[DELTA_TIME] = 1.0;
            r_process_info[DISTANCE_THRESHOLD] = 1.0;

            // First we create the nodes
            const std::size_t number_of_divisions = 8;
            const double lenght = 4.0;
            const double radius = 6.0;
            const double angle = Globals::Pi/6;
            const double slope = 0.0;

            // We create our problem
            const double delta_x = 0.1;
            CreateNewProblem(r_model_part, number_of_divisions, lenght, radius, angle, slope, delta_x);

            // We compute the explicit contribution
            const array_1d<double, 3> zero_vector = ZeroVector(3);;
            VariableUtils().SetScalarVar<Variable<double>>(WEIGHTED_GAP, 0.0, r_model_part.Nodes());
            VariableUtils().SetVectorVar(WEIGHTED_SLIP, zero_vector, r_model_part.Nodes());
            for (auto& r_cond : r_model_part.GetSubModelPart("ComputingContact").Conditions()) {
                r_cond.Set(MODIFIED, true);
                r_cond.AddExplicitContribution(r_process_info);
                r_cond.FinalizeSolutionStep(r_process_info);
            }

//             // DEBUG
//             GiDIOGapDebug(r_model_part);

            const double tolerance = 1.0e-4;
            array_1d<double, 3> slip = ZeroVector(3);
            slip[0] = delta_x;
            std::size_t counter = 0;
            for (auto& r_node : r_model_part.Nodes()) {
                if (r_node.Is(SLAVE)) {
                    const double normal_gap = r_node.GetValue(NORMAL_GAP);
                    const double weighted_gap_corrected = r_node.FastGetSolutionStepValue(WEIGHTED_GAP)/r_node.GetValue(NODAL_AREA);
                    if (std::abs(weighted_gap_corrected - normal_gap)/std::abs(normal_gap) < tolerance) {
                        if (norm_2(r_node.FastGetSolutionStepValue(WEIGHTED_SLIP)) > 0.0) {
                            const array_1d<double, 3> weighted_slip_corrected = r_node.FastGetSolutionStepValue(WEIGHTED_SLIP)/r_node.GetValue(NODAL_AREA);
                            KRATOS_CHECK_LESS_EQUAL(norm_2(weighted_slip_corrected - slip)/norm_2(slip), tolerance);
                            ++counter;
                        }
                    }
                }
            }
            KRATOS_CHECK(counter != 0);
        }

    } // namespace Testing
}  // namespace Kratos.
