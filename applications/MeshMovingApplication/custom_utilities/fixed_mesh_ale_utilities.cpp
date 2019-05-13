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

namespace Kratos
{
    /* Public functions *******************************************************/

    FixedMeshALEUtilities::FixedMeshALEUtilities(
        ModelPart &rVirtualModelPart,
        ModelPart &rStructureModelPart,
        const std::string LevelSetType) :
        mrVirtualModelPart(rVirtualModelPart),
        mrStructureModelPart(rStructureModelPart),
        mLevelSetType(LevelSetType)
    {
        // Get the default settings
        auto default_parameters = this->GetDefaultSettings();

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
        mrStructureModelPart(rModel.GetModelPart(rParameters["structure_model_part_name"].GetString())),
        mLevelSetType(rParameters["level_set_type"].GetString())
    {
        // Validate with default parameters
        rParameters.ValidateAndAssignDefaults(this->GetDefaultSettings());

        // Check the input level set type
        if (mLevelSetType != "continuous" && mLevelSetType != "discontinuous") {
            KRATOS_ERROR << "Provided level set type is: " << mLevelSetType << ". Only \"continuous\" and \"discontinuous\" types are supported.";
        }

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
        KRATOS_ERROR_IF(!rOriginModelPart.HasNodalSolutionStepVariable(DISPLACEMENT)) << "Missing required DISPLACEMENT variable in origin model part." << std::endl;
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
        auto &r_nodes_array = rOriginModelPart.NodesArray();
        for(auto &it_node : r_nodes_array){
            // Find the origin model part node in the virtual mesh
            Vector aux_N;
            Element::Pointer p_elem;
            const bool is_found = bin_based_point_locator.FindPointOnMeshSimplified(it_node->Coordinates(), aux_N, p_elem);

            // Check if the node is found
            if (is_found){
                // Initialize historical data
                // The current step values are also set as a prediction
                auto &r_mesh_vel = it_node->GetSolutionStepValue(MESH_VELOCITY);
                r_mesh_vel = ZeroVector(3);
                if (!it_node->IsFixed(PRESSURE)) {
                    for (unsigned int i_step = 0; i_step < BufferSize; ++i_step){
                        it_node->GetSolutionStepValue(PRESSURE, i_step) = 0.0;
                    }
                }
                if (!it_node->IsFixed(VELOCITY_X)) {
                    for (unsigned int i_step = 0; i_step < BufferSize; ++i_step){
                        it_node->GetSolutionStepValue(VELOCITY_X, i_step) = 0.0;
                    }
                }
                if (!it_node->IsFixed(VELOCITY_Y)) {
                    for (unsigned int i_step = 0; i_step < BufferSize; ++i_step){
                        it_node->GetSolutionStepValue(VELOCITY_Y, i_step) = 0.0;
                    }
                }
                if (!it_node->IsFixed(VELOCITY_Z)) {
                    for (unsigned int i_step = 0; i_step < BufferSize; ++i_step){
                        it_node->GetSolutionStepValue(VELOCITY_Z, i_step) = 0.0;
                    }
                }

                // Interpolate the origin model part nodal values
                const auto &r_geom = p_elem->GetGeometry();
                for (std::size_t i_virt_node = 0; i_virt_node < r_geom.PointsNumber(); ++i_virt_node){
                    r_mesh_vel += aux_N(i_virt_node) * r_geom[i_virt_node].GetSolutionStepValue(MESH_VELOCITY);
                    if (!it_node->IsFixed(PRESSURE)) {
                        for (unsigned int i_step = 0; i_step < BufferSize; ++i_step){
                            it_node->GetSolutionStepValue(PRESSURE, i_step) += aux_N(i_virt_node) * r_geom[i_virt_node].GetSolutionStepValue(PRESSURE, i_step);
                        }
                    }
                    if (!it_node->IsFixed(VELOCITY_X)) {
                        for (unsigned int i_step = 0; i_step < BufferSize; ++i_step){
                            it_node->GetSolutionStepValue(VELOCITY_X, i_step) += aux_N(i_virt_node) * r_geom[i_virt_node].GetSolutionStepValue(VELOCITY_X, i_step);
                        }
                    }
                    if (!it_node->IsFixed(VELOCITY_Y)) {
                        for (unsigned int i_step = 0; i_step < BufferSize; ++i_step){
                            it_node->GetSolutionStepValue(VELOCITY_Y, i_step) += aux_N(i_virt_node) * r_geom[i_virt_node].GetSolutionStepValue(VELOCITY_Y, i_step);
                        }
                    }
                    if (!it_node->IsFixed(VELOCITY_Z)) {
                        for (unsigned int i_step = 0; i_step < BufferSize; ++i_step){
                            it_node->GetSolutionStepValue(VELOCITY_Z, i_step) += aux_N(i_virt_node) * r_geom[i_virt_node].GetSolutionStepValue(VELOCITY_Z, i_step);
                        }
                    }
                }
            } else {
                KRATOS_WARNING("FixedMeshALEUtilities")
                    << "Origin model part node " << it_node->Id() << " has not been found in any virtual model part element. Origin node coordinates: (" << it_node->X() << " , " << it_node->Y() << " , " << it_node->Z() << ")" << std::endl;
            }
        }
    }

    /* Protected functions *******************************************************/

    FixedMeshALEUtilities::FixedMeshALEUtilities(
        ModelPart &rVirtualModelPart,
        ModelPart &rStructureModelPart) :
        mrVirtualModelPart(rVirtualModelPart),
        mrStructureModelPart(rStructureModelPart)
    {
        if (mrStructureModelPart.GetBufferSize() < 2)
        {
            (mrStructureModelPart.GetRootModelPart()).SetBufferSize(2);
            KRATOS_WARNING("FixedMeshALEUtilities") << "Structure model part buffer size is 1. Setting buffer size to 2." << std::endl;
        }
    }

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

    Parameters FixedMeshALEUtilities::GetDefaultSettings()
    {
        Parameters default_parameters(R"(
        {
            "virtual_model_part_name": "",
            "structure_model_part_name": "",
            "level_set_type": "",
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
        const unsigned int echo_level = 2;
        const bool compute_reactions = false;
        const bool calculate_norm_dx_flag = false;
        const bool reform_dof_set_at_each_step = false;
        mpMeshMovingStrategy = Kratos::make_shared<StrategyType>(
            mrVirtualModelPart,
            p_scheme,
            mpLinearSolver,
            p_builder_and_solver,
            compute_reactions,
            reform_dof_set_at_each_step,
            calculate_norm_dx_flag);

        mpMeshMovingStrategy->Check();
        mpMeshMovingStrategy->Initialize();
        mpMeshMovingStrategy->SetEchoLevel(echo_level);
    }

    const Vector FixedMeshALEUtilities::SetDistancesVector(ModelPart::ElementIterator ItElem) const
    {
        auto &r_geom = ItElem->GetGeometry();
        Vector nodal_distances(r_geom.PointsNumber());

        if (mLevelSetType == "continuous"){
            // Continuous nodal distance function case
            for (unsigned int i_node = 0; i_node < r_geom.PointsNumber(); ++i_node) {
                nodal_distances[i_node] = r_geom[i_node].FastGetSolutionStepValue(DISTANCE);
            }
        } else if (mLevelSetType == "discontinuous") {
            // Discontinuous elemental distance function case
            nodal_distances = ItElem->GetValue(ELEMENTAL_DISTANCES);
        } else {
            KRATOS_ERROR << "Level set type must be either 'continuous' or 'discontinuous'. Got " << mLevelSetType;
        }

        return nodal_distances;
    }

    inline bool FixedMeshALEUtilities::IsSplit(const Vector &rDistances) const
    {
        unsigned int n_pos = 0, n_neg = 0;
        for (double dist : rDistances) {
            if(dist >= 0) {
                ++n_pos;
            } else {
                ++n_neg;
            }
        }
        if (n_pos > 0 && n_neg > 0) {
            return true;
        }
        return false;
    }

    void FixedMeshALEUtilities::InitializeVirtualMeshValues()
    {
        // Copy the PRESSURE and VELOCITY values from the origin mesh to the virtual one
        VariableUtils::Pointer p_var_utils = Kratos::make_shared<VariableUtils>();
        for (unsigned int i_step = 0; i_step <  mpOriginModelPart->GetBufferSize(); ++i_step) {
            p_var_utils->CopyModelPartNodalVar<Variable<double>>(PRESSURE, *mpOriginModelPart, mrVirtualModelPart, i_step);
            p_var_utils->CopyModelPartNodalVar<Variable<array_1d<double,3>>>(VELOCITY, *mpOriginModelPart, mrVirtualModelPart, i_step);
        }

        // Initialize the DISPLACEMENT, MESH_DISPLACMENT and MESH_VELOCITY values
        p_var_utils->SetHistoricalVariableToZero(MESH_VELOCITY, mrVirtualModelPart.Nodes());
        p_var_utils->SetHistoricalVariableToZero(MESH_DISPLACEMENT, mrVirtualModelPart.Nodes());
    }

    void FixedMeshALEUtilities::InitializeMeshDisplacementFixity()
    {
        #pragma omp parallel for
        for(int i_fl = 0; i_fl < static_cast<int>(mrVirtualModelPart.NumberOfNodes()); ++i_fl) {
            // Free all the MESH_DISPLACEMENT DOFs
            auto it_node = mrVirtualModelPart.NodesBegin() + i_fl;
            it_node->Free(MESH_DISPLACEMENT_X);
            it_node->Free(MESH_DISPLACEMENT_Y);
            it_node->Free(MESH_DISPLACEMENT_Z);
        }
    }

    void FixedMeshALEUtilities::SetMeshDisplacementFixityFromOriginModelPart()
    {
        #pragma omp parallel for
        for(int i_fl = 0; i_fl < static_cast<int>(mrVirtualModelPart.NumberOfNodes()); ++i_fl) {
            // Check origin model part mesh displacement fixity
            auto it_node = mrVirtualModelPart.NodesBegin() + i_fl;
            const auto it_orig_node = mpOriginModelPart->NodesBegin() + i_fl;
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
    }

    void FixedMeshALEUtilities::SetEmbeddedNodalMeshDisplacement()
    {
        // Initialize the DISPLACEMENT variable values
        #pragma omp parallel for
        for (int i_node = 0; i_node < static_cast<int>(mpOriginModelPart->NumberOfNodes()); ++i_node) {
            auto it_node = mpOriginModelPart->NodesBegin() + i_node;
            it_node->FastGetSolutionStepValue(DISPLACEMENT, 0) = ZeroVector(3);
            it_node->FastGetSolutionStepValue(DISPLACEMENT, 1) = ZeroVector(3);
        }

        // Compute the DISPLACEMENT increment from the structure model part and save it in the origin mesh MESH_DISPLACEMENT
        const unsigned int buff_pos_0 = 0;
        FixedMeshALEUtilities::EmbeddedNodalVariableProcessArrayType emb_nod_var_from_skin_proc_array_0(
            *mpOriginModelPart,
            mrStructureModelPart,
            mEmbeddedNodalVariableSettings["linear_solver_settings"],
            DISPLACEMENT,
            DISPLACEMENT,
            mEmbeddedNodalVariableSettings["gradient_penalty_coefficient"].GetDouble(),
            mLevelSetType,
            buff_pos_0);
        emb_nod_var_from_skin_proc_array_0.Execute();
        emb_nod_var_from_skin_proc_array_0.Clear();

        const unsigned int buff_pos_1 = 1;
        FixedMeshALEUtilities::EmbeddedNodalVariableProcessArrayType emb_nod_var_from_skin_proc_array_1(
            *mpOriginModelPart,
            mrStructureModelPart,
            mEmbeddedNodalVariableSettings["linear_solver_settings"],
            DISPLACEMENT,
            DISPLACEMENT,
            mEmbeddedNodalVariableSettings["gradient_penalty_coefficient"].GetDouble(),
            mLevelSetType,
            buff_pos_1);
        emb_nod_var_from_skin_proc_array_1.Execute();
        emb_nod_var_from_skin_proc_array_1.Clear();

        // In the intersected elements, set the MESH_DISPLACEMENT as the increment of displacement and fix it
        for (int i_elem = 0; i_elem < static_cast<int>(mpOriginModelPart->NumberOfElements()); ++i_elem) {
            auto it_elem = mpOriginModelPart->ElementsBegin() + i_elem;
            if (this->IsSplit(this->SetDistancesVector(it_elem))) {
                const auto &r_geom = it_elem->GetGeometry();
                auto &r_virt_geom = (mrVirtualModelPart.ElementsBegin() + i_elem)->GetGeometry();
                for (unsigned int i_node = 0; i_node < r_virt_geom.PointsNumber(); ++i_node) {
                    const auto &r_d_0 = r_geom[i_node].FastGetSolutionStepValue(DISPLACEMENT, 0);
                    const auto &r_d_1 = r_geom[i_node].FastGetSolutionStepValue(DISPLACEMENT, 1);
                    noalias(r_virt_geom[i_node].FastGetSolutionStepValue(MESH_DISPLACEMENT, 0)) = r_d_0 - r_d_1;
                    r_virt_geom[i_node].Fix(MESH_DISPLACEMENT_X);
                    r_virt_geom[i_node].Fix(MESH_DISPLACEMENT_Y);
                    r_virt_geom[i_node].Fix(MESH_DISPLACEMENT_Z);
                }
            }
        }
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
        for (const auto it_elem : mrVirtualModelPart.ElementsArray()) {
            KRATOS_ERROR_IF((it_elem->GetGeometry()).Area() < 0.0) << "Element " << it_elem->Id() << " in virtual model part has negative jacobian." << std::endl;
        }
#endif
    }

    void FixedMeshALEUtilities::RevertMeshDisplacementFixity()
    {
        #pragma omp parallel for
        for (int i_node = 0; i_node < static_cast<int>(mrVirtualModelPart.NumberOfNodes()); ++i_node) {
            auto it_node = mrVirtualModelPart.NodesBegin() + i_node;
            if (it_node->IsFixed(MESH_DISPLACEMENT_X)) {
                it_node->Free(MESH_DISPLACEMENT_X);
            }
            if (it_node->IsFixed(MESH_DISPLACEMENT_Y)) {
                it_node->Free(MESH_DISPLACEMENT_Y);
            }
            if (it_node->IsFixed(MESH_DISPLACEMENT_Z)) {
                it_node->Free(MESH_DISPLACEMENT_Z);
            }
        }
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

    template void FixedMeshALEUtilities::ProjectVirtualValues<2>(ModelPart &rOriginModelPart, const unsigned int BufferSize);
    template void FixedMeshALEUtilities::ProjectVirtualValues<3>(ModelPart &rOriginModelPart, const unsigned int BufferSize);
}
