// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Franziska Wahl
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "includes/model_part.h"
#include "includes/model_part_io.h"
#include "includes/convection_diffusion_settings.h"
#include "includes/variables.h"
#include "includes/dof.h"
#include "processes/find_nodal_neighbours_process.h"

// Application includes
#include "convection_diffusion_application.h"
#include "convection_diffusion_application_variables.h"
#include "custom_utilities/embedded_mls_constraint_process.h"


namespace Kratos::Testing
{
    void SetEmbeddedConstraintProcessTestModelPart(ModelPart &rModelPart)
    {
        // Set buffer size
        rModelPart.SetBufferSize(1);

        // Set convection diffusion settings
        ConvectionDiffusionSettings::Pointer p_conv_dff_set = Kratos::make_shared<ConvectionDiffusionSettings>();
        p_conv_dff_set->SetDensityVariable(DENSITY);
        p_conv_dff_set->SetDiffusionVariable(CONDUCTIVITY);
        p_conv_dff_set->SetUnknownVariable(TEMPERATURE);
        p_conv_dff_set->SetVolumeSourceVariable(HEAT_FLUX);
        p_conv_dff_set->SetSurfaceSourceVariable(FACE_HEAT_FLUX);
        p_conv_dff_set->SetProjectionVariable(PROJECTED_SCALAR1);
        p_conv_dff_set->SetConvectionVariable(CONVECTION_VELOCITY);
        p_conv_dff_set->SetMeshVelocityVariable(MESH_VELOCITY);
        p_conv_dff_set->SetVelocityVariable(VELOCITY);
        p_conv_dff_set->SetSpecificHeatVariable(SPECIFIC_HEAT);
        p_conv_dff_set->SetReactionVariable(REACTION_FLUX);
        rModelPart.GetProcessInfo().SetValue(CONVECTION_DIFFUSION_SETTINGS, p_conv_dff_set);

        // Variables addition
        rModelPart.AddNodalSolutionStepVariable(DENSITY);
        rModelPart.AddNodalSolutionStepVariable(CONDUCTIVITY);
        rModelPart.AddNodalSolutionStepVariable(TEMPERATURE);
        rModelPart.AddNodalSolutionStepVariable(HEAT_FLUX);
        rModelPart.AddNodalSolutionStepVariable(FACE_HEAT_FLUX);
        rModelPart.AddNodalSolutionStepVariable(PROJECTED_SCALAR1);
        rModelPart.AddNodalSolutionStepVariable(CONVECTION_VELOCITY);
        rModelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
        rModelPart.AddNodalSolutionStepVariable(VELOCITY);
        rModelPart.AddNodalSolutionStepVariable(SPECIFIC_HEAT);
        rModelPart.AddNodalSolutionStepVariable(REACTION_FLUX);
        rModelPart.AddNodalSolutionStepVariable(DISTANCE);

        // Create a fake properties container
        auto p_elem_prop = rModelPart.CreateNewProperties(0);
    }

    KRATOS_TEST_CASE_IN_SUITE(EmbeddedMLSConstraintProcess, KratosConvectionDiffusionFastSuite)
    {
        // Get test model part
        Model model;
        auto &r_test_model_part = model.CreateModelPart("TestModelPart");
        SetEmbeddedConstraintProcessTestModelPart(r_test_model_part);

        // Read mdpa and set domain size
        ModelPartIO(std::filesystem::path(__FILE__).parent_path().string()+"/test_embedded_mls_constraint_process/embedded_mls_constraints_test", IO::READ).ReadModelPart(r_test_model_part);
        auto& r_process_info = r_test_model_part.GetProcessInfo();
        r_process_info[DOMAIN_SIZE] = 2;

        // Add TEMPERATURE Dof and set DISTANCE values
        const double dist_y = 0.8889;
        for (auto& r_node : r_test_model_part.Nodes()) {
            r_node.AddDof(TEMPERATURE);
            const double dist = dist_y - r_node[1];
            r_node.FastGetSolutionStepValue(DISTANCE) = dist;
        }

        // Run find neighbors processes
        FindGlobalNodalNeighboursProcess find_nodal_neighbors_process(r_test_model_part);
        find_nodal_neighbors_process.Execute();

        // Run embedded MLS constraint process
        Parameters mls_constraint_params( R"(
        {
            "model_part_name"                        : "TestModelPart",
            "mls_extension_operator_order"           : 1,
            "deactivate_negative_elements"           : true,
            "deactivate_intersected_elements"        : false
        }  )" );
        EmbeddedMLSConstraintProcess mls_constraint_process(model, mls_constraint_params);
        mls_constraint_process.Execute();

        // Define expected master-slave constraints
        std::unordered_map< std::size_t, std::vector<std::size_t> > mastersIds_of_slaveId = {
            { 10, {16,23,17,33,24,21} },
            { 11, {16,17,23,24,21,25,29,33,34} },
            { 14, {17,21,16,23,24,29,25,33,34,36,32,31} },
            { 20, {21,25,17,24,29,32,31,16,23,33,34,36,39,37,38} },
            { 26, {25,31,21,29,32,37,38,17,24,34,36,39,42,44} },
            { 35, {31,38,25,32,37,42,21,29,36,39,44,47} }
        };
        std::unordered_map< std::size_t, std::vector<bool> > mastersIds_of_slaveId_found = {
            { 10, {false,false,false,false,false,false} },
            { 11, {false,false,false,false,false,false,false,false,false} },
            { 14, {false,false,false,false,false,false,false,false,false,false,false,false} },
            { 20, {false,false,false,false,false,false,false,false,false,false,false,false,false,false,false} },
            { 26, {false,false,false,false,false,false,false,false,false,false,false,false,false,false} },
            { 35, {false,false,false,false,false,false,false,false,false,false,false,false} }
        };

        std::vector< Dof<double>::Pointer > dof_list_slave, dof_list_master;

        // Get the array of constraints from the modeler
        const auto &r_constraints_array = r_test_model_part.MasterSlaveConstraints();
        const size_t number_of_constraints = static_cast<int>(r_constraints_array.size());
        for (std::size_t i = 0; i < number_of_constraints; ++i) {
            const auto it_constraints = r_constraints_array.begin() + i;

            // Gets list of Dof involved on every element
            it_constraints->GetDofList(dof_list_slave, dof_list_master, r_process_info);
            const std::size_t slave_id = (*dof_list_slave.begin())->GetId();
            const std::size_t master_id = (*dof_list_master.begin())->GetId();

            // Check for constraint in unordered map and set to 1 if found, error if not found
            const auto& master_ids = mastersIds_of_slaveId[slave_id];
            const auto it_master = std::find(master_ids.begin(), master_ids.end(), master_id);
            if (it_master != master_ids.end()) {
                mastersIds_of_slaveId_found[slave_id][it_master - master_ids.begin()] = true;
            } else {
                KRATOS_ERROR << "Master node not found in list of expected contraint masters.";
            }
        }

        // Check in unordered map if all constraints where found
        for (const auto& m : mastersIds_of_slaveId_found) {
            const auto& master_ids_found = m.second;
            const auto it_found = std::find(master_ids_found.begin(), master_ids_found.end(), false);
            if (it_found != master_ids_found.end()) {
                KRATOS_ERROR << "Expected master-slave-constraint was not found.";
            }
        }
    }

} // namespace Kratos::Testing
