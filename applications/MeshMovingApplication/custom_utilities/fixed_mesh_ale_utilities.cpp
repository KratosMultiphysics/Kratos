//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

// System includes

// External includes

// Application includes
#include "custom_elements/structural_meshmoving_element.h"
#include "custom_utilities/fixed_mesh_ale_utilities.h"
#include "custom_utilities/mesh_velocity_calculation.h"
#include "custom_utilities/move_mesh_utilities.h"

// Project includes
#include "includes/mesh_moving_variables.h"
#include "utilities/binbased_fast_point_locator.h"
#include "utilities/variable_utils.h"
#include "utilities/parallel_utilities.h"

namespace Kratos
{
    /* Public functions *******************************************************/

    FixedMeshALEUtilities::FixedMeshALEUtilities(
        ModelPart &rVirtualModelPart,
        ModelPart &rStructureModelPart) :
        mrVirtualModelPart(rVirtualModelPart),
        mrStructureModelPart(rStructureModelPart)
    {
        // Get the default settings
        auto default_parameters = this->GetDefaultParameters();

        // Set default embedded nodal variable settings
        mEmbeddedNodalVariableSettings = default_parameters["embedded_nodal_variable_settings"];

        // Set default linear solver pointer
        this->SetLinearSolverPointer(default_parameters["linear_solver_settings"]);

        // Check the structure model part
        if (mrStructureModelPart.GetBufferSize() < 2) {
            (mrStructureModelPart.GetRootModelPart()).SetBufferSize(2);
            KRATOS_WARNING("FixedMeshALEUtilities") << "Structure model part buffer size is 1. Setting buffer size to 2." << std::endl;
        }
    }

    FixedMeshALEUtilities::FixedMeshALEUtilities(
        Model &rModel,
        Parameters &rParameters) :
        mrVirtualModelPart(rModel.GetModelPart(rParameters["virtual_model_part_name"].GetString())),
        mrStructureModelPart(rModel.GetModelPart(rParameters["structure_model_part_name"].GetString()))
    {
        // Validate with default parameters
        rParameters.ValidateAndAssignDefaults(this->GetDefaultParameters());

        // Save the embedded nodal variable settings
        mEmbeddedNodalVariableSettings = rParameters["embedded_nodal_variable_settings"];

        // Set the linear solver pointer
        this->SetLinearSolverPointer(rParameters["linear_solver_settings"]);

        // Check the structure model part
        if (mrStructureModelPart.GetBufferSize() < 2) {
            (mrStructureModelPart.GetRootModelPart()).SetBufferSize(2);
            KRATOS_WARNING("FixedMeshALEUtilities") << "Structure model part buffer size is 1. Setting buffer size to 2." << std::endl;
        }
    }

    void FixedMeshALEUtilities::Initialize(ModelPart &rOriginModelPart)
    {
        // Fill the virtual model part as a copy of the origin mesh
        this->FillVirtualModelPart(rOriginModelPart);

        // Create and initialize the mesh moving strategy
        this->SetMeshMovingStrategy();
    }

    void FixedMeshALEUtilities::FillVirtualModelPart(ModelPart &rOriginModelPart)
    {
        // Check that the origin model part has nodes and elements to be copied
        KRATOS_ERROR_IF(rOriginModelPart.NumberOfNodes() == 0) << "Origin model part has no nodes.";
        KRATOS_ERROR_IF(rOriginModelPart.NumberOfElements() == 0) << "Origin model part has no elements.";

        // Check that the origin model part has the required variables
        KRATOS_ERROR_IF(!rOriginModelPart.HasNodalSolutionStepVariable(MESH_VELOCITY)) << "Missing required MESH_VELOCITY variable in origin model part." << std::endl;
        KRATOS_ERROR_IF(!rOriginModelPart.HasNodalSolutionStepVariable(MESH_DISPLACEMENT)) << "Missing required MESH_DISPLACEMENT variable in origin model part." << std::endl;

        // Save the selected origin model part
        mpOriginModelPart = &rOriginModelPart;

        // Set the origin process info to the virtual model part
        mrVirtualModelPart.SetProcessInfo(rOriginModelPart.GetProcessInfo());

        // Add the required varibles to the virtual model part
        mrVirtualModelPart.AddNodalSolutionStepVariable(VELOCITY);
        mrVirtualModelPart.AddNodalSolutionStepVariable(PRESSURE);
        mrVirtualModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);
        mrVirtualModelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
        mrVirtualModelPart.AddNodalSolutionStepVariable(MESH_DISPLACEMENT);

        // Set the buffer size in the virtual model part
        mrVirtualModelPart.SetBufferSize(rOriginModelPart.GetBufferSize());

        // Copy the origin model part nodes
        auto &r_nodes_array = rOriginModelPart.NodesArray();
        for(auto &it_node : r_nodes_array){
            // Create a copy of the origin model part node and add DOFs
            auto p_node = mrVirtualModelPart.CreateNewNode(it_node->Id(), *it_node, 0);
            p_node->pAddDof(MESH_DISPLACEMENT_X);
            p_node->pAddDof(MESH_DISPLACEMENT_Y);
            p_node->pAddDof(MESH_DISPLACEMENT_Z);
        }

        // Copy the origin model part elements
        this->CreateVirtualModelPartElements(rOriginModelPart);

        // Check that the nodes and elements have been correctly copied
        KRATOS_ERROR_IF(rOriginModelPart.NumberOfNodes() != mrVirtualModelPart.NumberOfNodes())
            << "Origin and virtual model part have different number of nodes.";
        KRATOS_ERROR_IF(rOriginModelPart.NumberOfElements() != mrVirtualModelPart.NumberOfElements())
            << "Origin and virtual model part have different number of elements.";
    }

    void FixedMeshALEUtilities::SetVirtualMeshValuesFromOriginMesh()
    {
        auto virt_nodes_begin = mrVirtualModelPart.NodesBegin();
        auto orig_nodes_begin = mpOriginModelPart->NodesBegin();
        const unsigned int buffer_size = mpOriginModelPart->GetBufferSize();

        IndexPartition<std::size_t>( static_cast<int>(mpOriginModelPart->NumberOfNodes()) ).for_each(
            [&]( std::size_t index )
            {
                auto it_virt_node       = virt_nodes_begin + index;
                const auto it_orig_node = orig_nodes_begin + index;

                for (unsigned int step = 1; step < buffer_size; ++step) {
                    it_virt_node->FastGetSolutionStepValue(PRESSURE, step) = it_orig_node->FastGetSolutionStepValue(PRESSURE, step);
                    noalias(it_virt_node->FastGetSolutionStepValue(VELOCITY, step)) = it_orig_node->FastGetSolutionStepValue(VELOCITY, step);
                }
            } );
    }

    void FixedMeshALEUtilities::ComputeMeshMovement(const double DeltaTime)
    {
        // Initialize the PRESSURE and VELOCITY virtual mesh values
        this->InitializeVirtualMeshValues();

        // Initialize the MESH_DISPLACEMENT fixity
        this->InitializeMeshDisplacementFixity();

        // Get the MESH_DISPLACEMENT fixity from the origin model part
        this->SetMeshDisplacementFixityFromOriginModelPart();

        // Set embedded MESH_DISPLACEMENT from the immersed structure
        this->SetEmbeddedNodalMeshDisplacement();

        // Set the mesh moving strategy
        this->SolveMeshMovementStrategy(DeltaTime);
    }

    void FixedMeshALEUtilities::UndoMeshMovement()
    {
        // Revert the MESH_DISPLACEMENT fixity
        this->RevertMeshDisplacementFixity();

        // Revert the virtual mesh movement to its original configuration
        auto &r_nodes = mrVirtualModelPart.Nodes();
        VariableUtils().UpdateCurrentToInitialConfiguration(r_nodes);
    }

    template <unsigned int TDim>
    void FixedMeshALEUtilities::ProjectVirtualValues(
        ModelPart& rOriginModelPart,
        const unsigned int BufferSize)
    {
        // Check that the virtual model part has elements
        KRATOS_ERROR_IF(mrVirtualModelPart.NumberOfNodes() == 0) << "Virtual model part has no nodes.";
        KRATOS_ERROR_IF(mrVirtualModelPart.NumberOfElements() == 0) << "Virtual model part has no elements.";

        // Set the binbased fast point locator utility
        BinBasedFastPointLocator<TDim> bin_based_point_locator(mrVirtualModelPart);
        bin_based_point_locator.UpdateSearchDatabase();

        // Search the origin model part nodes in the virtual mesh elements and
        // interpolate the values in the virtual element to the origin model part node
        block_for_each( rOriginModelPart.Nodes(),
            [&]( Node<3>& rNode )
            {
                // Find the origin model part node in the virtual mesh
                Vector aux_N;
                Element::Pointer p_elem = nullptr;
                const bool is_found = bin_based_point_locator.FindPointOnMeshSimplified(rNode.Coordinates(), aux_N, p_elem);

                // Check if the node is found
                if (is_found){
                    // Initialize MESH_VELOCITY
                    noalias(rNode.FastGetSolutionStepValue(MESH_VELOCITY)) = ZeroVector(3);

                    // Initialize historical data
                    for (unsigned int i_step = 1; i_step < BufferSize; ++i_step) {
                        rNode.FastGetSolutionStepValue(PRESSURE, i_step) = 0.0;
                        noalias(rNode.FastGetSolutionStepValue(VELOCITY, i_step)) = ZeroVector(3);
                    }

                    // Interpolate the origin model part nodal values
                    const auto &r_geom = p_elem->GetGeometry();
                    for (std::size_t i_virt_node = 0; i_virt_node < r_geom.PointsNumber(); ++i_virt_node){
                        // Project MESH_VELOCITY
                        const auto i_virt_v_mesh = r_geom[i_virt_node].FastGetSolutionStepValue(MESH_VELOCITY);
                        rNode.FastGetSolutionStepValue(MESH_VELOCITY) += aux_N(i_virt_node) * i_virt_v_mesh;

                        // Project historical data
                        for (unsigned int i_step = 1; i_step < BufferSize; ++i_step){
                            const auto &i_virt_v = r_geom[i_virt_node].FastGetSolutionStepValue(VELOCITY, i_step);
                            const double &i_virt_p = r_geom[i_virt_node].FastGetSolutionStepValue(PRESSURE, i_step);
                            rNode.FastGetSolutionStepValue(PRESSURE, i_step) += aux_N(i_virt_node) * i_virt_p;
                            noalias(rNode.FastGetSolutionStepValue(VELOCITY, i_step)) += aux_N(i_virt_node) * i_virt_v;
                        }
                    }
                } else {
                    KRATOS_WARNING("FixedMeshALEUtilities")
                        << "Origin model part node " << rNode.Id() << " has not been found in any virtual model part element. Origin node coordinates: (" << rNode.X() << " , " << rNode.Y() << " , " << rNode.Z() << ")" << std::endl;
                } // else ( is_found )
            } );

    }

    /* Protected functions *******************************************************/

    void FixedMeshALEUtilities::CreateVirtualModelPartElements(const ModelPart &rOriginModelPart)
    {
        auto &r_elems = rOriginModelPart.Elements();
        ModelPart::ElementsContainerType new_elements_container;
        for(auto &elem : r_elems) {
            // Set the array of virtual nodes to create the element from the original ids.
            NodesArrayType new_nodes_array;
            const auto &r_orig_geom = elem.GetGeometry();
            for (unsigned int i = 0; i < r_orig_geom.PointsNumber(); ++i){
                new_nodes_array.push_back(mrVirtualModelPart.pGetNode(r_orig_geom[i].Id()));
            }
            auto p_new_geom = r_orig_geom.Create(new_nodes_array);

            // Create a structural mesh moving element with the same Id() but the virtual mesh nodes
            auto p_elem = Kratos::make_intrusive<StructuralMeshMovingElement>(elem.Id(), p_new_geom, elem.pGetProperties());
            new_elements_container.push_back(p_elem);
        }
        mrVirtualModelPart.AddElements(new_elements_container.begin(), new_elements_container.end());
    }

    /* Private functions *******************************************************/

    Parameters FixedMeshALEUtilities::GetDefaultParameters()
    {
        Parameters default_parameters(R"(
        {
            "virtual_model_part_name": "",
            "structure_model_part_name": "",
            "linear_solver_settings": {
                "solver_type": "cg",
                "tolerance": 1.0e-8,
                "max_iteration": 1000
            },
            "embedded_nodal_variable_settings": {
                "gradient_penalty_coefficient": 0.0,
                "linear_solver_settings": {
                    "preconditioner_type": "amg",
                    "solver_type": "amgcl",
                    "smoother_type": "ilu0",
                    "krylov_type": "cg",
                    "max_iteration": 1000,
                    "verbosity": 0,
                    "tolerance": 1e-8,
                    "scaling": false,
                    "block_size": 1,
                    "use_block_matrices_if_possible": true
                }
            }
        })");
        return default_parameters;
    }

    void FixedMeshALEUtilities::SetLinearSolverPointer(const Parameters &rLinearSolverSettings)
    {
        FixedMeshALEUtilities::LinearSolverFactoryType linear_solver_factory;
        mpLinearSolver = linear_solver_factory.Create(rLinearSolverSettings);
    }

    void FixedMeshALEUtilities::SetMeshMovingStrategy()
    {
        // Set the scheme and buildar and solver pointers
        SchemePointerType p_scheme = Kratos::make_shared<SchemeType>();
        BuilderAndSolverPointerType p_builder_and_solver = Kratos::make_shared<BuilderAndSolverType>(mpLinearSolver);

        // Set a linear strategy to solve the mesh moving problem
        const unsigned int echo_level = 0;
        const bool compute_reactions = false;
        const bool calculate_norm_dx_flag = false;
        const bool reform_dof_set_at_each_step = false;
        mpMeshMovingStrategy = Kratos::make_shared<StrategyType>(
            mrVirtualModelPart,
            p_scheme,
            p_builder_and_solver,
            compute_reactions,
            reform_dof_set_at_each_step,
            calculate_norm_dx_flag);

        mpMeshMovingStrategy->Check();
        mpMeshMovingStrategy->Initialize();
        mpMeshMovingStrategy->SetEchoLevel(echo_level);
    }

    void FixedMeshALEUtilities::InitializeVirtualMeshValues()
    {
        // Initialize the MESH_DISPLACEMENT and MESH_VELOCITY values
        // Note that both positions of the buffer are initialized to zero. This is important
        // in case the CloneTimeStep() is done in the virtual model part, since the method
        // assumes that the mesh is in the origin position when computing the MESH_VELOCITY.
        block_for_each( mrVirtualModelPart.Nodes(),
            []( Node<3>& rNode )
            {
                rNode.FastGetSolutionStepValue(MESH_VELOCITY, 0) = ZeroVector(3);
                rNode.FastGetSolutionStepValue(MESH_VELOCITY, 1) = ZeroVector(3);
                rNode.FastGetSolutionStepValue(MESH_DISPLACEMENT, 0) = ZeroVector(3);
                rNode.FastGetSolutionStepValue(MESH_DISPLACEMENT, 1) = ZeroVector(3);
            } );
    }

    void FixedMeshALEUtilities::InitializeMeshDisplacementFixity()
    {
        block_for_each( mrVirtualModelPart.Nodes(),
            []( Node<3>& rNode )
            {
                rNode.Free(MESH_DISPLACEMENT_X);
                rNode.Free(MESH_DISPLACEMENT_Y);
                rNode.Free(MESH_DISPLACEMENT_Z);
            } );
    }

    void FixedMeshALEUtilities::SetMeshDisplacementFixityFromOriginModelPart()
    {
        IndexPartition<std::size_t>( static_cast<size_t>(mrVirtualModelPart.NumberOfNodes()) ).for_each(
            [this]( std::size_t index )
            {
                auto it_node = mrVirtualModelPart.NodesBegin() + index;
                const auto it_orig_node = mpOriginModelPart->NodesBegin() + index;
                if (it_orig_node->IsFixed(MESH_DISPLACEMENT_X)) {
                    it_node->Fix(MESH_DISPLACEMENT_X);
                }
                if (it_orig_node->IsFixed(MESH_DISPLACEMENT_Y)) {
                    it_node->Fix(MESH_DISPLACEMENT_Y);
                }
                if (it_orig_node->IsFixed(MESH_DISPLACEMENT_Z)) {
                    it_node->Fix(MESH_DISPLACEMENT_Z);
                }
            }
        );
    }

    void FixedMeshALEUtilities::SetEmbeddedNodalMeshDisplacement()
    {
        // Initialize the DISPLACEMENT variable values
        block_for_each( mrVirtualModelPart.Nodes(),
        []( Node<3>& rNode )
        {
            rNode.FastGetSolutionStepValue(DISPLACEMENT, 0) = ZeroVector(3);
            rNode.FastGetSolutionStepValue(DISPLACEMENT, 1) = ZeroVector(3);
        } );

        // Place the structure in its previous position
        block_for_each( mrStructureModelPart.Nodes(),
        []( Node<3>& rNode )
        {
            const auto d_1 = rNode.FastGetSolutionStepValue(DISPLACEMENT, 1);
            rNode.X() = rNode.X0() + d_1[0];
            rNode.Y() = rNode.Y0() + d_1[1];
            rNode.Z() = rNode.Z0() + d_1[2];
        } );

        // Compute the DISPLACEMENT increment from the structure model part and save it in the origin mesh MESH_DISPLACEMENT
        const unsigned int buff_pos_0 = 0;
        FixedMeshALEUtilities::EmbeddedNodalVariableProcessArrayType emb_nod_var_from_skin_proc_array_0(
            mrVirtualModelPart,
            mrStructureModelPart,
            mEmbeddedNodalVariableSettings["linear_solver_settings"],
            DISPLACEMENT,
            DISPLACEMENT,
            mEmbeddedNodalVariableSettings["gradient_penalty_coefficient"].GetDouble(),
            buff_pos_0);
        emb_nod_var_from_skin_proc_array_0.Execute();
        emb_nod_var_from_skin_proc_array_0.Clear();

        const unsigned int buff_pos_1 = 1;
        FixedMeshALEUtilities::EmbeddedNodalVariableProcessArrayType emb_nod_var_from_skin_proc_array_1(
            mrVirtualModelPart,
            mrStructureModelPart,
            mEmbeddedNodalVariableSettings["linear_solver_settings"],
            DISPLACEMENT,
            DISPLACEMENT,
            mEmbeddedNodalVariableSettings["gradient_penalty_coefficient"].GetDouble(),
            buff_pos_1);
        emb_nod_var_from_skin_proc_array_1.Execute();
        emb_nod_var_from_skin_proc_array_1.Clear();

        // In the intersected elements, set the MESH_DISPLACEMENT as the increment of displacement and fix it
        // Note that this assumes that the flag INTERFACE has been set by the embedded nodal variable process
        for (int i_elem = 0; i_elem < static_cast<int>(mrVirtualModelPart.NumberOfElements()); ++i_elem) {
            auto it_elem = mrVirtualModelPart.ElementsBegin() + i_elem;
            if (it_elem->Is(INTERFACE)) {
                auto &r_geom = it_elem->GetGeometry();
                for (unsigned int i_node = 0; i_node < r_geom.PointsNumber(); ++i_node) {
                    auto &r_node = r_geom[i_node];
                    if (r_node.Is(VISITED)) {
                        const auto &r_d_0 = r_node.FastGetSolutionStepValue(DISPLACEMENT, 0);
                        const auto &r_d_1 = r_node.FastGetSolutionStepValue(DISPLACEMENT, 1);
                        auto &r_mesh_disp = r_node.FastGetSolutionStepValue(MESH_DISPLACEMENT, 0);
                        if (!r_node.IsFixed(MESH_DISPLACEMENT_X)) {
                            r_node.Fix(MESH_DISPLACEMENT_X);
                            r_mesh_disp[0] = r_d_0[0] - r_d_1[0];
                        }
                        if (!r_node.IsFixed(MESH_DISPLACEMENT_Y)) {
                            r_node.Fix(MESH_DISPLACEMENT_Y);
                            r_mesh_disp[1] = r_d_0[1] - r_d_1[1];
                        }
                        if (!r_node.IsFixed(MESH_DISPLACEMENT_Z)) {
                            r_node.Fix(MESH_DISPLACEMENT_Z);
                            r_mesh_disp[2] = r_d_0[2] - r_d_1[2];
                        }
                    }
                }
            }
        }

        // Revert the structure movement to its current position
        block_for_each( mrStructureModelPart.Nodes(),
        []( Node<3>& rNode )
        {
            const auto d_0 = rNode.FastGetSolutionStepValue(DISPLACEMENT, 0);
            rNode.X() = rNode.X0() + d_0[0];
            rNode.Y() = rNode.Y0() + d_0[1];
            rNode.Z() = rNode.Z0() + d_0[2];
        } );
    }

    void FixedMeshALEUtilities::SolveMeshMovementStrategy(const double DeltaTime)
    {
        // Set the time increment in the virtual model part so the MESH_VELOCITY can be computed
        mrVirtualModelPart.GetProcessInfo()[DELTA_TIME] = DeltaTime;

        // Solve the mesh problem
        mpMeshMovingStrategy->Solve();
        TimeDiscretization::BDF1 time_disc_BDF1;
        MeshVelocityCalculation::CalculateMeshVelocities(mrVirtualModelPart, time_disc_BDF1);
        MoveMeshUtilities::MoveMesh(mrVirtualModelPart.Nodes());

        // Check that the moved virtual mesh has no negative Jacobian elements
#ifdef KRATOS_DEBUG
        for (const auto& it_elem : mrVirtualModelPart.ElementsArray()) {
            KRATOS_ERROR_IF((it_elem->GetGeometry()).Area() < 0.0) << "Element " << it_elem->Id() << " in virtual model part has negative jacobian." << std::endl;
        }
#endif
    }

    void FixedMeshALEUtilities::RevertMeshDisplacementFixity()
    {
        block_for_each( mrVirtualModelPart.Nodes(),
        []( Node<3>& rNode )
        {
            if (rNode.IsFixed(MESH_DISPLACEMENT_X)) {
                rNode.Free(MESH_DISPLACEMENT_X);
            }
            if (rNode.IsFixed(MESH_DISPLACEMENT_Y)) {
                rNode.Free(MESH_DISPLACEMENT_Y);
            }
            if (rNode.IsFixed(MESH_DISPLACEMENT_Z)) {
                rNode.Free(MESH_DISPLACEMENT_Z);
            }
        } );
    }

    /* External functions *****************************************************/

    /// output stream function
    inline std::ostream& operator << (
        std::ostream& rOStream,
        const FixedMeshALEUtilities& rThis)
    {
        rThis.PrintData(rOStream);
        return rOStream;
    }

    template void KRATOS_API(MESH_MOVING_APPLICATION) FixedMeshALEUtilities::ProjectVirtualValues<2>(ModelPart &rOriginModelPart, const unsigned int BufferSize);
    template void KRATOS_API(MESH_MOVING_APPLICATION) FixedMeshALEUtilities::ProjectVirtualValues<3>(ModelPart &rOriginModelPart, const unsigned int BufferSize);
}
