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
#include "custom_elements/laplacian_meshmoving_element.h"
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
        // Set the linear solver pointer
        Parameters linear_solver_settings(R"({
            "solver_type": "cg",
            "tolerance": 1.0e-8,
            "max_iteration": 200
        })");
        this->SetLinearSolverPointer(linear_solver_settings);

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
        Parameters default_parameters(R"(
        {
            "virtual_model_part_name": "",
            "structure_model_part_name": "",
            "level_set_type": "",
            "linear_solver_settings": {
                "solver_type": "cg",
                "tolerance": 1.0e-8,
                "max_iteration": 200
            }
        }  )");
        rParameters.ValidateAndAssignDefaults(default_parameters);

        // Check the input level set type
        if (mLevelSetType != "continuous" || mLevelSetType != "discontinuous") {
            KRATOS_ERROR << "Provided level set type is: " << mLevelSetType << ". Only \"continuous\" and \"discontinuous\" types are supported.";
        }

        // Set the linear solver pointer
        this->SetLinearSolverPointer(rParameters["linear_solver_settings"]);

        // Check the structure model part
        if (mrStructureModelPart.GetBufferSize() < 2) {
            (mrStructureModelPart.GetRootModelPart()).SetBufferSize(2);
            KRATOS_WARNING("FixedMeshALEUtilities") << "Structure model part buffer size is 1. Setting buffer size to 2." << std::endl;
        }
    }

    void FixedMeshALEUtilities::FillVirtualModelPart(ModelPart& rOriginModelPart)
    {
        // Check that the origin model part has nodes and elements to be copied
        KRATOS_ERROR_IF(rOriginModelPart.NumberOfNodes() == 0) << "Origin model part has no nodes.";
        KRATOS_ERROR_IF(rOriginModelPart.NumberOfElements() == 0) << "Origin model part has no elements.";

        // Save the selected origin model part
        mpOriginModelPart = &rOriginModelPart;

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
            auto p_node = mrVirtualModelPart.CreateNewNode(it_node->Id(),*it_node, 0);
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
        // Get the MESH_DISPLACEMENT fixity from the origin model part
        this->SetMeshDisplacementFixityFromOriginModelPart();

        // Set embedded MESH_DISPLACEMENT from the immersed structure
        this->SetEmbeddedNodalMeshDisplacement();

        // Set the mesh moving strategy
        this->SetAndSolveMeshMovementStrategy(DeltaTime);
    }

    void FixedMeshALEUtilities::UndoMeshMovement()
    {
        auto &r_nodes = mrVirtualModelPart.Nodes();
        VariableUtils().UpdateCurrentToInitialConfiguration(r_nodes);
    }

    template <unsigned int TDim>
    void FixedMeshALEUtilities::ProjectVirtualValues(
        ModelPart& rOriginModelPart,
        const unsigned int BufferSize)
    {
        // Check that the virtual model part has elements
        KRATOS_ERROR_IF(mrVirtualModelPart.NumberOfElements() == 0) << "Virtual model part has no elements.";

        // Set the binbased fast point locator utility
        BinBasedFastPointLocator<TDim> bin_based_point_locator(mrVirtualModelPart);
        bin_based_point_locator.UpdateSearchDatabase();

        // Search the origin model part nodes in the virtual mesh elements and
        // interpolate the values in the virtual element to the origin model part node
        auto &r_nodes_array = rOriginModelPart.NodesArray();
        for(auto it_node : r_nodes_array){
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
                auto &r_geom = p_elem->GetGeometry();
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
                KRATOS_WARNING("ExplicitMeshMovingUtility") << "Origin model part node " << it_node->Id() << " has not been found in any virtual model part element." << std::endl;
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
        for(auto &elem : r_elems) {
            // Set the array of virtual nodes to create the element from the original ids.
            NodesArrayType new_nodes_array;
            auto &r_orig_geom = elem.GetGeometry();
            for (unsigned int i = 0; i < r_orig_geom.PointsNumber(); ++i){
                new_nodes_array.push_back(mrVirtualModelPart.pGetNode(r_orig_geom[i].Id()));
            }
            auto p_new_geom = r_orig_geom.Create(new_nodes_array);

            // Create a Laplacian mesh moving element with the same Id() but the virtual mesh nodes
            auto p_elem = Kratos::make_shared<LaplacianMeshMovingElement>(elem.Id(), p_new_geom, elem.pGetProperties());
            mrVirtualModelPart.AddElement(p_elem);
        }
    }

    /* Private functions *******************************************************/

    void FixedMeshALEUtilities::SetLinearSolverPointer(const Parameters &rLinearSolverSettings)
    {
        FixedMeshALEUtilities::LinearSolverFactoryType linear_solver_factory;
        mpLinearSolver = linear_solver_factory.Create(rLinearSolverSettings);
    }

    const Vector FixedMeshALEUtilities::SetDistancesVector(ModelPart::ElementIterator ItElem)
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

    inline bool FixedMeshALEUtilities::IsSplit(const Vector &rDistances)
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
        // Compute the DISPLACEMENT from the structure model part and save it in the origin mesh MESH_DISPLACEMENT
        FixedMeshALEUtilities::EmbeddedNodalVariableProcessArrayType emb_nod_var_from_skin_proc_array(
            *mpOriginModelPart,
            mrStructureModelPart,
            DISPLACEMENT,
            MESH_DISPLACEMENT,
            mLevelSetType);
        emb_nod_var_from_skin_proc_array.Execute();
        emb_nod_var_from_skin_proc_array.Clear();

        // TODO: ADD OpenMP PARALLELIZATION
        // In the intersected elements, set the MESH_DISPLACEMENT as the increment of displacement and fix it
        // Note that we take advantage of the fact that the virtual mesh nodes and the origin mesh nodes
        // have the same ids.
        for (int i_elem = 0; i_elem < static_cast<int>(mpOriginModelPart->NumberOfElements()); ++i_elem) {
            auto it_elem = mpOriginModelPart->ElementsBegin() + i_elem;
            if (this->IsSplit(this->SetDistancesVector(it_elem))) {
                const auto &r_geom = it_elem->GetGeometry();
                for (auto &r_node : r_geom) {
                    const auto d_0 = r_node.FastGetSolutionStepValue(MESH_DISPLACEMENT, 0);
                    const auto d_1 = r_node.FastGetSolutionStepValue(MESH_DISPLACEMENT, 1);
                    auto p_virt_node = mrVirtualModelPart.pGetNode(r_node.Id());
                    noalias(p_virt_node->FastGetSolutionStepValue(MESH_DISPLACEMENT, 0)) = d_0 - d_1;
                    p_virt_node->Fix(MESH_DISPLACEMENT_X);
                    p_virt_node->Fix(MESH_DISPLACEMENT_Y);
                    p_virt_node->Fix(MESH_DISPLACEMENT_Z);
                }
            }
        }
    }

    void FixedMeshALEUtilities::SetAndSolveMeshMovementStrategy(const double DeltaTime)
    {
        // Set the Laplacian mesh moving strategy
        auto p_mesh_moving_strategy = Kratos::make_shared<LaplacianMeshMovingStrategyType>(mrVirtualModelPart, mpLinearSolver);
        p_mesh_moving_strategy->Check();
        p_mesh_moving_strategy->Initialize();

        // Set the time increment in the virtual model part so the MESH_VELOCITY can be computed
        mrVirtualModelPart.GetProcessInfo()[DELTA_TIME] = DeltaTime;

        // Solve the mesh problem
        p_mesh_moving_strategy->Solve();
        p_mesh_moving_strategy->UpdateReferenceMesh();
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
