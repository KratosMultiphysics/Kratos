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

// System includes
#include <set>

// External includes

// Project includes
#include "custom_processes/mmg_process.h"
#include "containers/model.h"
#include "utilities/assign_unique_model_part_collection_tag_utility.h"
// We indlude the internal variable interpolation process
#include "custom_processes/nodal_values_interpolation_process.h"
#include "custom_processes/internal_variables_interpolation_process.h"
#include "processes/fast_transfer_between_model_parts_process.h"
// Include the point locator
#include "utilities/binbased_fast_point_locator.h"
// Include the spatial containers needed for search
#include "spatial_containers/spatial_containers.h" // kd-tree
#include "includes/gid_io.h"
#include "includes/model_part_io.h"


// NOTE: The following contains the license of the MMG library
/* =============================================================================
**  Copyright (c) Bx INP/Inria/UBordeaux/UPMC, 2004- .
**
**  mmg is free software: you can redistribute it and/or modify it
**  under the terms of the GNU Lesser General Public License as published
**  by the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  mmg is distributed in the hope that it will be useful, but WITHOUT
**  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
**  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
**  License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License and of the GNU General Public License along with mmg (in
**  files COPYING.LESSER and COPYING). If not, see
**  <http://www.gnu.org/licenses/>. Please read their terms carefully and
**  use this copy of the mmg distribution only if you accept them.
** =============================================================================
*/

namespace Kratos
{

/************************************* CONSTRUCTOR *********************************/
/***********************************************************************************/

template<MMGLibrary TMMGLibrary>
MmgProcess<TMMGLibrary>::MmgProcess(
    ModelPart& rThisModelPart,
    Parameters ThisParameters
    ):mrThisModelPart(rThisModelPart),
      mThisParameters(ThisParameters)
{
    Parameters default_parameters = GetDefaultParameters();
    mThisParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

    mStdStringFilename = mThisParameters["filename"].GetString();
    mEchoLevel = mThisParameters["echo_level"].GetInt();

    mFilename = new char [mStdStringFilename.length() + 1];
    std::strcpy (mFilename, mStdStringFilename.c_str());

    // The framework type
    mFramework = ConvertFramework(mThisParameters["framework"].GetString());

    // The discretization type
    mDiscretization = ConvertDiscretization(mThisParameters["discretization_type"].GetString());

    if ( mDiscretization == DiscretizationOption::ISOSURFACE ){
        mRemoveRegions = mThisParameters["isosurface_parameters"]["remove_regions"].GetBool();
    } else{
        mRemoveRegions = false;
    }

    mpRefElement.clear();
    mpRefCondition.clear();
}

/*************************************** EXECUTE ***********************************/
/***********************************************************************************/

template<MMGLibrary TMMGLibrary>
void MmgProcess<TMMGLibrary>::Execute()
{
    KRATOS_TRY;

    // We execute all the necessary steps
    ExecuteInitialize();
    ExecuteInitializeSolutionStep();
    ExecuteFinalize();

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibrary TMMGLibrary>
void MmgProcess<TMMGLibrary>::ExecuteInitialize()
{
    KRATOS_TRY;

    if( mRemoveRegions ){
        // the conditions are re-creted in the process
        mrThisModelPart.Conditions().clear();
        KRATOS_INFO("MmgProcess") << "Conditions were cleared" << std::endl;
    }

    /* We restart the MMG mesh and solution */
    mMmmgUtilities.SetEchoLevel(mEchoLevel);
    mMmmgUtilities.SetDiscretization(mDiscretization);
    mMmmgUtilities.SetRemoveRegions(mRemoveRegions);
    mMmmgUtilities.InitMesh();

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibrary TMMGLibrary>
void MmgProcess<TMMGLibrary>::ExecuteBeforeSolutionLoop()
{
    KRATOS_TRY;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibrary TMMGLibrary>
void MmgProcess<TMMGLibrary>::ExecuteInitializeSolutionStep()
{
    KRATOS_TRY;

    const bool safe_to_file = mThisParameters["save_external_files"].GetBool();

    /* We print the original model part */
    KRATOS_INFO_IF("", mEchoLevel > 0) <<
    "//---------------------------------------------------//" << std::endl <<
    "//---------------------------------------------------//" << std::endl <<
    "//---------------  BEFORE REMESHING   ---------------//" << std::endl <<
    "//---------------------------------------------------//" << std::endl <<
    "//---------------------------------------------------//" << std::endl <<
    std::endl << mrThisModelPart << std::endl;

    // We initialize the mesh and solution data
    InitializeMeshData();

    // We retrieve the data form the Kratos model part to fill sol
    if (mDiscretization == DiscretizationOption::STANDARD) {
        InitializeSolDataMetric();
    } else if (mDiscretization == DiscretizationOption::ISOSURFACE) {
        InitializeSolDataDistance();
    } else {
        KRATOS_ERROR << "Discretization type: " << static_cast<int>(mDiscretization) << " not fully implemented" << std::endl;
    }

    // Check if the number of given entities match with mesh size
    mMmmgUtilities.CheckMeshData();

    // Save to file
    if (safe_to_file) SaveSolutionToFile(false);

    // We execute the remeshing
    ExecuteRemeshing();

    /* We print the resulting model part */
    KRATOS_INFO_IF("", mEchoLevel > 0) <<
    "//---------------------------------------------------//" << std::endl <<
    "//---------------------------------------------------//" << std::endl <<
    "//---------------   AFTER REMESHING   ---------------//" << std::endl <<
    "//---------------------------------------------------//" << std::endl <<
    "//---------------------------------------------------//" << std::endl <<
    std::endl << mrThisModelPart << std::endl;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibrary TMMGLibrary>
void MmgProcess<TMMGLibrary>::ExecuteFinalizeSolutionStep()
{
    KRATOS_TRY;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibrary TMMGLibrary>
void MmgProcess<TMMGLibrary>::ExecuteBeforeOutputStep()
{
    KRATOS_TRY;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibrary TMMGLibrary>
void MmgProcess<TMMGLibrary>::ExecuteAfterOutputStep()
{
    KRATOS_TRY;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibrary TMMGLibrary>
void MmgProcess<TMMGLibrary>::ExecuteFinalize()
{
    KRATOS_TRY;

    // We release the memory
    FreeMemory();

    if( mRemoveRegions ){
        // nodes not belonging to an element are removed
        CleanSuperfluousNodes();
    }

    KRATOS_CATCH("");
}

/************************************* OPERATOR() **********************************/
/***********************************************************************************/

template<MMGLibrary TMMGLibrary>
void MmgProcess<TMMGLibrary>::operator()()
{
    Execute();
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibrary TMMGLibrary>
void MmgProcess<TMMGLibrary>::InitializeMeshData()
{
    // We create a list of submodelparts to later reassign flags after remesh
    CreateAuxiliarSubModelPartForFlags();

    // Before computing colors we do some check and throw a warning to get the user informed
    const std::vector<std::string> sub_model_part_names = AssignUniqueModelPartCollectionTagUtility::GetRecursiveSubModelPartNames(mrThisModelPart);

    for (auto sub_model_part_name : sub_model_part_names) {
        ModelPart& r_sub_model_part = AssignUniqueModelPartCollectionTagUtility::GetRecursiveSubModelPart(mrThisModelPart, sub_model_part_name);

        KRATOS_WARNING_IF("MmgProcess", mEchoLevel > 0 && (r_sub_model_part.NumberOfNodes() > 0 && (r_sub_model_part.NumberOfConditions() == 0 && r_sub_model_part.NumberOfElements() == 0))) <<
        "The submodelpart: " << sub_model_part_name << " contains only nodes and no geometries (conditions/elements)." << std::endl <<
        "It is not guaranteed that the submodelpart will be preserved." << std::endl <<
        "PLEASE: Add some \"dummy\" conditions to the submodelpart to preserve it" << std::endl;
    }

    // First we compute the colors
    mColors.clear();
    ColorsMapType nodes_colors, cond_colors, elem_colors;
    AssignUniqueModelPartCollectionTagUtility model_part_collections(mrThisModelPart);
    model_part_collections.ComputeTags(nodes_colors, cond_colors, elem_colors, mColors);

    /////////* MESH FILE */////////
    // Build mesh in MMG5 format //

    // Iterate over components
    NodesArrayType& r_nodes_array = mrThisModelPart.Nodes();
    ConditionsArrayType& r_conditions_array = mrThisModelPart.Conditions();
    ElementsArrayType& r_elements_array = mrThisModelPart.Elements();

    /* Manually set of the mesh */
    MMGMeshInfo<TMMGLibrary> mmg_mesh_info;
    if (TMMGLibrary == MMGLibrary::MMG2D) { // 2D
        mmg_mesh_info.NumberOfLines = r_conditions_array.size();
        mmg_mesh_info.NumberOfTriangles = r_elements_array.size();
    } else if (TMMGLibrary == MMGLibrary::MMG3D) { // 3D
        /* Conditions */
        std::size_t num_tri = 0, num_quad = 0;
        #pragma omp parallel for reduction(+:num_tri,num_quad)
        for(int i = 0; i < static_cast<int>(r_conditions_array.size()); ++i) {
            auto it_cond = r_conditions_array.begin() + i;

            if ((it_cond->GetGeometry()).GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle3D3) { // Triangles
                num_tri += 1;
            } else if ((it_cond->GetGeometry()).GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4) { // Quadrilaterals
                num_quad += 1;
            } else
                KRATOS_WARNING("MmgProcess") << "WARNING: YOUR GEOMETRY CONTAINS " << it_cond->GetGeometry().PointsNumber() <<" NODES THAT CAN NOT BE REMESHED" << std::endl;
        }

        mmg_mesh_info.NumberOfTriangles = num_tri;
        mmg_mesh_info.NumberOfQuadrilaterals = num_quad;

        KRATOS_INFO_IF("MmgProcess", ((num_tri + num_quad) < r_conditions_array.size()) && mEchoLevel > 0) <<
        "Number of Conditions: " << r_conditions_array.size() << " Number of Triangles: " << num_tri << " Number of Quadrilaterals: " << num_quad << std::endl;

        /* Elements */
        std::size_t num_tetra = 0, num_prisms = 0;
        #pragma omp parallel for reduction(+:num_tetra,num_prisms)
        for(int i = 0; i < static_cast<int>(r_elements_array.size()); ++i) {
            auto it_elem = r_elements_array.begin() + i;

            if ((it_elem->GetGeometry()).GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4) { // Tetrahedron
                num_tetra += 1;
            } else if ((it_elem->GetGeometry()).GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Prism3D6) { // Prisms
                num_prisms += 1;
            } else
                KRATOS_WARNING("MmgProcess") << "WARNING: YOUR GEOMETRY CONTAINS " << it_elem->GetGeometry().PointsNumber() <<" NODES CAN NOT BE REMESHED" << std::endl;
        }

        mmg_mesh_info.NumberOfTetrahedra = num_tetra;
        mmg_mesh_info.NumberOfPrism = num_prisms;

        KRATOS_INFO_IF("MmgProcess", ((num_tetra + num_prisms) < r_elements_array.size()) && mEchoLevel > 0) <<
        "Number of Elements: " << r_elements_array.size() << " Number of Tetrahedron: " << num_tetra << " Number of Prisms: " << num_prisms << std::endl;
    } else { // Surfaces
        mmg_mesh_info.NumberOfLines = r_conditions_array.size();
        mmg_mesh_info.NumberOfTriangles = r_elements_array.size();
    }

    mmg_mesh_info.NumberOfNodes = r_nodes_array.size();
    mMmmgUtilities.SetMeshSize(mmg_mesh_info);

    /* Nodes */
    // We copy the DOF from the first node (after we release, to avoid problem with previous conditions)
    mDofs = r_nodes_array.begin()->GetDofs();
    for (typename NodeType::DofsContainerType::const_iterator it_dof = mDofs.begin(); it_dof != mDofs.end(); ++it_dof)
        it_dof->FreeDof();

    if (mFramework == FrameworkEulerLagrange::LAGRANGIAN){ // NOTE: The code is repeated due to performance reasons
        #pragma omp parallel for firstprivate(nodes_colors)
        for(int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i) {
            auto it_node = r_nodes_array.begin() + i;

            mMmmgUtilities.SetNodes(it_node->X0(), it_node->Y0(), it_node->Z0(), nodes_colors[it_node->Id()], i + 1);

            bool blocked = false;
            if (it_node->IsDefined(BLOCKED))
                blocked = it_node->Is(BLOCKED);
            if (blocked)
                mMmmgUtilities.BlockNode(i + 1);

            // RESETING THE ID OF THE NODES (important for non consecutive meshes)
            it_node->SetId(i + 1);
        }
    }
    else {
        #pragma omp parallel for firstprivate(nodes_colors)
        for(int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i) {
            auto it_node = r_nodes_array.begin() + i;

            mMmmgUtilities.SetNodes(it_node->X(), it_node->Y(), it_node->Z(), nodes_colors[it_node->Id()], i + 1);

            bool blocked = false;
            if (it_node->IsDefined(BLOCKED))
                blocked = it_node->Is(BLOCKED);
            if (blocked)
                mMmmgUtilities.BlockNode(i + 1);

            // RESETING THE ID OF THE NODES (important for non consecutive meshes)
            it_node->SetId(i + 1);
        }
    }

    /* Conditions */
    #pragma omp parallel for firstprivate(cond_colors)
    for(int i = 0; i < static_cast<int>(r_conditions_array.size()); ++i)  {
        auto it_cond = r_conditions_array.begin() + i;

        mMmmgUtilities.SetConditions(it_cond->GetGeometry(), cond_colors[it_cond->Id()], i + 1);

        bool blocked = false;
        if (it_cond->IsDefined(BLOCKED))
            blocked = it_cond->Is(BLOCKED);
        if (blocked)
            mMmmgUtilities.BlockCondition(i + 1);

        // RESETING THE ID OF THE CONDITIONS (important for non consecutive meshes)
        it_cond->SetId(i + 1);
    }

    /* Elements */
    #pragma omp parallel for firstprivate(elem_colors)
    for(int i = 0; i < static_cast<int>(r_elements_array.size()); ++i) {
        auto it_elem = r_elements_array.begin() + i;

        mMmmgUtilities.SetElements(it_elem->GetGeometry(), elem_colors[it_elem->Id()], i + 1);

        bool blocked = false;
        if (it_elem->IsDefined(BLOCKED))
            blocked = it_elem->Is(BLOCKED);
        if (blocked)
            mMmmgUtilities.BlockElement(i + 1);

        // RESETING THE ID OF THE ELEMENTS (important for non consecutive meshes)
        it_elem->SetId(i + 1);
    }

    // Create auxiliar colors maps
    ColorsMapType aux_ref_cond;
    for(int i = 0; i < static_cast<int>(r_conditions_array.size()); ++i)  {
        auto it_cond = r_conditions_array.begin() + i;
        const IndexType cond_id = it_cond->Id();
        const IndexType color = cond_colors[cond_id];
        if (!(aux_ref_cond.find(color) != aux_ref_cond.end()))
            aux_ref_cond.insert (IndexPairType(color,cond_id));
    }
    ColorsMapType aux_ref_elem;
    for(int i = 0; i < static_cast<int>(r_elements_array.size()); ++i) {
        auto it_elem = r_elements_array.begin() + i;
        const IndexType elem_id = it_elem->Id();
        const IndexType color = elem_colors[elem_id];
        if (!(aux_ref_elem.find(color) != aux_ref_elem.end()))
            aux_ref_elem.insert (IndexPairType(color,elem_id));
    }

    /* We clone the first condition and element of each type (we will assume that each sub model part has just one kind of condition, in my opinion it is quite reccomended to create more than one sub model part if you have more than one element or condition) */
    // First we add the main model part
    if (r_conditions_array.size() > 0) {
        const std::string type_name = (Dimension == 2) ? "Condition2D2N" : (TMMGLibrary == MMGLibrary::MMG3D) ? "Condition3D" : "Condition3D2N";
        Condition const& r_clone_condition = KratosComponents<Condition>::Get(type_name);
        mpRefCondition[0] = r_clone_condition.Create(0, r_clone_condition.GetGeometry(), r_conditions_array.begin()->pGetProperties());
    }
    if (r_elements_array.size() > 0) {
        mpRefElement[0] = r_elements_array.begin()->Create(0, r_elements_array.begin()->GetGeometry(), r_elements_array.begin()->pGetProperties());
    }

    // Now we add the reference elements and conditions
    for (auto& ref_cond : aux_ref_cond) {
        Condition::Pointer p_cond = mrThisModelPart.pGetCondition(ref_cond.second);
        mpRefCondition[ref_cond.first] = p_cond->Create(0, p_cond->GetGeometry(), p_cond->pGetProperties());
    }
    for (auto& ref_elem : aux_ref_elem) {
        Element::Pointer p_elem = mrThisModelPart.pGetElement(ref_elem.second);
        mpRefElement[ref_elem.first] = p_elem->Create(0, p_elem->GetGeometry(), p_elem->pGetProperties());
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibrary TMMGLibrary>
void MmgProcess<TMMGLibrary>::InitializeSolDataMetric()
{
    ////////* SOLUTION FILE *////////

    // Iterate in the nodes
    NodesArrayType& r_nodes_array = mrThisModelPart.Nodes();

    // Set size of the solution
    mMmmgUtilities.SetSolSizeTensor(r_nodes_array.size());

    const Variable<TensorArrayType>& r_tensor_variable = KratosComponents<Variable<TensorArrayType>>::Get("METRIC_TENSOR_"+std::to_string(Dimension)+"D");

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i) {
        auto it_node = r_nodes_array.begin() + i;

        KRATOS_DEBUG_ERROR_IF_NOT(it_node->Has(r_tensor_variable)) << "METRIC_TENSOR_" + std::to_string(Dimension) + "D  not defined for node " << it_node->Id() << std::endl;

        // We get the metric
        const TensorArrayType& metric = it_node->GetValue(r_tensor_variable);

        // We set the metric
        mMmmgUtilities.SetMetricTensor(metric, i + 1);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibrary TMMGLibrary>
void MmgProcess<TMMGLibrary>::InitializeSolDataDistance()
{
    ////////* SOLUTION FILE for ISOSURFACE*////////
    // Iterate in the nodes
    NodesArrayType& r_nodes_array = mrThisModelPart.Nodes();

    // Set size of the solution
    mMmmgUtilities.SetSolSizeScalar(r_nodes_array.size() );

    // GEtting variable for scalar filed
    const std::string& r_isosurface_variable_name = mThisParameters["isosurface_parameters"]["isosurface_variable"].GetString();
    const bool nonhistorical_variable = mThisParameters["isosurface_parameters"]["nonhistorical_variable"].GetBool();
    const Variable<double>& r_scalar_variable = KratosComponents<Variable<double>>::Get(r_isosurface_variable_name);

    // Auxiliar value
    double isosurface_value = 0.0;

    // We iterate over the nodes
    #pragma omp parallel for firstprivate(isosurface_value)
    for(int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i) {
        auto it_node = r_nodes_array.begin() + i;

        if (nonhistorical_variable) {
            KRATOS_DEBUG_ERROR_IF_NOT(it_node->Has(r_scalar_variable)) << r_isosurface_variable_name << " field not found as a non-historical variable " << std::endl;

            // We get the isosurface value (non-historical variable)
            isosurface_value = it_node->GetValue( r_scalar_variable );
        } else {
            KRATOS_DEBUG_ERROR_IF_NOT(it_node->SolutionStepsDataHas(r_scalar_variable)) << r_isosurface_variable_name << " field not found as a historical variable " << std::endl;

            // We get the isosurface value (historical variable)
            isosurface_value = it_node->FastGetSolutionStepValue( r_scalar_variable );
        }

        // We set the isosurface variable
        mMmmgUtilities.SetMetricScalar(isosurface_value, i + 1);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibrary TMMGLibrary>
void MmgProcess<TMMGLibrary>::ExecuteRemeshing()
{
    // Getting the parameters
    const bool save_to_file = mThisParameters["save_external_files"].GetBool();

    // We initialize some values
    const SizeType step_data_size = mrThisModelPart.GetNodalSolutionStepDataSize();
    const SizeType buffer_size   = mrThisModelPart.NodesBegin()->GetBufferSize();

    mThisParameters["step_data_size"].SetInt(step_data_size);
    mThisParameters["buffer_size"].SetInt(buffer_size);

    KRATOS_INFO_IF("MmgProcess", mEchoLevel > 0) << "Step data size: " << step_data_size << " Buffer size: " << buffer_size << std::endl;

    ////////* MMG LIBRARY CALL *////////
    KRATOS_INFO_IF("MmgProcess", mEchoLevel > 0) << "////////* MMG LIBRARY CALL *////////" << std::endl;

    // Calling the library functions
    if (mDiscretization == DiscretizationOption::STANDARD) {
        mMmmgUtilities.MMGLibCallMetric(mThisParameters);
    } else if (mDiscretization == DiscretizationOption::ISOSURFACE) {
        mMmmgUtilities.MMGLibCallIsoSurface(mThisParameters);
    } else {
        KRATOS_ERROR << "Discretization type: " << static_cast<int>(mDiscretization) << " not fully implemented" << std::endl;
    }

    /* Save to file */
    if (save_to_file) SaveSolutionToFile(true);

    // Some information
    MMGMeshInfo<TMMGLibrary> mmg_mesh_info;
    mMmmgUtilities.PrintAndGetMmgMeshInfo(mmg_mesh_info);

    ////////* EMPTY AND BACKUP THE MODEL PART *////////
    Model& owner_model = mrThisModelPart.GetModel();
    ModelPart& r_old_model_part = owner_model.CreateModelPart(mrThisModelPart.Name()+"_Old", mrThisModelPart.GetBufferSize());

    // First we empty the model part
    NodesArrayType& r_nodes_array = mrThisModelPart.Nodes();

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i)
        (r_nodes_array.begin() + i)->Set(TO_ERASE, true);
    r_old_model_part.AddNodes( mrThisModelPart.NodesBegin(), mrThisModelPart.NodesEnd() );
    mrThisModelPart.RemoveNodesFromAllLevels(TO_ERASE);

    ConditionsArrayType& r_conditions_array = mrThisModelPart.Conditions();

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(r_conditions_array.size()); ++i)
        (r_conditions_array.begin() + i)->Set(TO_ERASE, true);
    r_old_model_part.AddConditions( mrThisModelPart.ConditionsBegin(), mrThisModelPart.ConditionsEnd() );
    mrThisModelPart.RemoveConditionsFromAllLevels(TO_ERASE);

    ElementsArrayType& r_elements_array = mrThisModelPart.Elements();

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(r_elements_array.size()); ++i)
        (r_elements_array.begin() + i)->Set(TO_ERASE, true);
    r_old_model_part.AddElements( mrThisModelPart.ElementsBegin(), mrThisModelPart.ElementsEnd() );
    mrThisModelPart.RemoveElementsFromAllLevels(TO_ERASE);

    // Create a new model part // TODO: Use a different kind of element for each submodelpart (in order to be able of remeshing more than one kind o element or condition)
    std::unordered_map<IndexType, IndexVectorType> color_nodes, first_color_cond, second_color_cond, first_color_elem, second_color_elem;

    // The tempotal store of
    ConditionsArrayType created_conditions_vector;
    ElementsArrayType created_elements_vector;

    // Auxiliar values
    int ref, is_required;

    /* NODES */ // TODO: ADD OMP
    for (IndexType i_node = 1; i_node <= mmg_mesh_info.NumberOfNodes; ++i_node) {
        NodeType::Pointer p_node = mMmmgUtilities.CreateNode(mrThisModelPart, i_node, ref, is_required);

        // Set the DOFs in the nodes
        for (auto& r_dof : mDofs)
            p_node->pAddDof(r_dof);

        if (ref != 0) color_nodes[static_cast<IndexType>(ref)].push_back(i_node);// NOTE: ref == 0 is the MainModelPart
    }

    /* CONDITIONS */ // TODO: ADD OMP
    if (mpRefCondition.size() > 0) {
        IndexType cond_id = 1;

        IndexType counter_first_cond = 0;
        const IndexVectorType first_condition_to_remove = mMmmgUtilities.CheckFirstTypeConditions();
        for (IndexType i_cond = 1; i_cond <= mmg_mesh_info.NumberFirstTypeConditions(); ++i_cond) {
            bool skip_creation = false;
            if (counter_first_cond < first_condition_to_remove.size()) {
                if (first_condition_to_remove[counter_first_cond] == i_cond) {
                    skip_creation = true;
                    counter_first_cond += 1;
                }
            }

            Condition::Pointer p_condition = mMmmgUtilities.CreateFirstTypeCondition(mrThisModelPart, mpRefCondition, cond_id, ref, is_required, skip_creation);

            if (p_condition.get() != nullptr) {
                created_conditions_vector.push_back(p_condition);
//                 mrThisModelPart.AddCondition(p_condition);
                if (ref != 0) first_color_cond[static_cast<IndexType>(ref)].push_back(cond_id);// NOTE: ref == 0 is the MainModelPart
                cond_id += 1;
            }
        }

        IndexType counter_second_cond = 0;
        const IndexVectorType second_condition_to_remove = mMmmgUtilities.CheckSecondTypeConditions();
        for (IndexType i_cond = 1; i_cond <= mmg_mesh_info.NumberSecondTypeConditions(); ++i_cond) {
            bool skip_creation = false;
            if (counter_second_cond < second_condition_to_remove.size()) {
                if (second_condition_to_remove[counter_second_cond] == i_cond) {
                    skip_creation = true;
                    counter_second_cond += 1;
                }
            }
            Condition::Pointer p_condition = mMmmgUtilities.CreateSecondTypeCondition(mrThisModelPart, mpRefCondition, cond_id, ref, is_required, skip_creation);

            if (p_condition.get() != nullptr) {
                created_conditions_vector.push_back(p_condition);
//                 mrThisModelPart.AddCondition(p_condition);
                if (ref != 0) second_color_cond[static_cast<IndexType>(ref)].push_back(cond_id);// NOTE: ref == 0 is the MainModelPart
                cond_id += 1;
            }
        }
    }

    /* ELEMENTS */ // TODO: ADD OMP
    if (mpRefElement.size() > 0) {
        IndexType elem_id = 1;

        IndexType counter_first_elem = 0;
        const IndexVectorType first_elements_to_remove = mMmmgUtilities.CheckFirstTypeElements();
        for (IndexType i_elem = 1; i_elem <= mmg_mesh_info.NumberFirstTypeElements(); ++i_elem) {
            bool skip_creation = false;
            if (counter_first_elem < first_elements_to_remove.size()) {
                if (first_elements_to_remove[counter_first_elem] == i_elem) {
                    skip_creation = true;
                    counter_first_elem += 1;
                }
            }

            Element::Pointer p_element = mMmmgUtilities.CreateFirstTypeElement(mrThisModelPart, mpRefElement, elem_id, ref, is_required, skip_creation);

            if (p_element.get() != nullptr) {
                created_elements_vector.push_back(p_element);
//                 mrThisModelPart.AddElement(p_element);
                if (ref != 0) first_color_elem[static_cast<IndexType>(ref)].push_back(elem_id);// NOTE: ref == 0 is the MainModelPart
                elem_id += 1;
            }
        }

        IndexType counter_second_elem = 0;
        const IndexVectorType second_elements_to_remove = mMmmgUtilities.CheckSecondTypeElements();
        for (IndexType i_elem = 1; i_elem <= mmg_mesh_info.NumberSecondTypeElements(); ++i_elem) {
            bool skip_creation = false;
            if (counter_second_elem < second_elements_to_remove.size()) {
                if (second_elements_to_remove[counter_second_elem] == i_elem) {
                    skip_creation = true;
                    counter_second_elem += 1;
                }
            }

            Element::Pointer p_element = mMmmgUtilities.CreateSecondTypeElement(mrThisModelPart, mpRefElement, elem_id, ref, is_required,skip_creation);

            if (p_element.get() != nullptr) {
                created_elements_vector.push_back(p_element);
//                 mrThisModelPart.AddElement(p_element);
                if (ref != 0) second_color_elem[static_cast<IndexType>(ref)].push_back(elem_id);// NOTE: ref == 0 is the MainModelPart
                elem_id += 1;
            }
        }
    }

    // Finally we add the conditions and elements to the main model part
    mrThisModelPart.AddConditions(created_conditions_vector.begin(), created_conditions_vector.end());
    mrThisModelPart.AddElements(created_elements_vector.begin(), created_elements_vector.end());

    // We add nodes, conditions and elements to the sub model parts
    for (auto & color_list : mColors) {
        const IndexType key = color_list.first;

        if (key != 0) {// NOTE: key == 0 is the MainModelPart
            for (auto sub_model_part_name : color_list.second) {
                ModelPart& r_sub_model_part = AssignUniqueModelPartCollectionTagUtility::GetRecursiveSubModelPart(mrThisModelPart, sub_model_part_name);

                if (color_nodes.find(key) != color_nodes.end()) r_sub_model_part.AddNodes(color_nodes[key]);
                if (first_color_cond.find(key) != first_color_cond.end()) r_sub_model_part.AddConditions(first_color_cond[key]);
                if (second_color_cond.find(key) != second_color_cond.end()) r_sub_model_part.AddConditions(second_color_cond[key]);
                if (first_color_elem.find(key) != first_color_elem.end()) r_sub_model_part.AddElements(first_color_elem[key]);
                if (second_color_elem.find(key) != second_color_elem.end()) r_sub_model_part.AddElements(second_color_elem[key]);
            }
        }
    }

    // TODO: Add OMP
    // NOTE: We add the nodes from the elements and conditions to the respective submodelparts
    const std::vector<std::string> sub_model_part_names = AssignUniqueModelPartCollectionTagUtility::GetRecursiveSubModelPartNames(mrThisModelPart);

    for (auto sub_model_part_name : sub_model_part_names) {
        ModelPart& r_sub_model_part = AssignUniqueModelPartCollectionTagUtility::GetRecursiveSubModelPart(mrThisModelPart, sub_model_part_name);

        std::unordered_set<IndexType> node_ids;

        ConditionsArrayType& sub_conditions_array = r_sub_model_part.Conditions();
        const SizeType sub_num_conditions = sub_conditions_array.end() - sub_conditions_array.begin();

        for(IndexType i = 0; i < sub_num_conditions; ++i)  {
            auto it_cond = sub_conditions_array.begin() + i;
            auto& cond_geom = it_cond->GetGeometry();

            for (SizeType i_node = 0; i_node < cond_geom.size(); ++i_node)
                node_ids.insert(cond_geom[i_node].Id());
        }

        ElementsArrayType& sub_elements_array = r_sub_model_part.Elements();
        const SizeType sub_num_elements = sub_elements_array.end() - sub_elements_array.begin();

        for(IndexType i = 0; i < sub_num_elements; ++i) {
            auto it_elem = sub_elements_array.begin() + i;
            auto& elem_geom = it_elem->GetGeometry();

            for (SizeType i_node = 0; i_node < elem_geom.size(); ++i_node)
                node_ids.insert(elem_geom[i_node].Id());
        }

        IndexVectorType vector_ids;
        std::copy(node_ids.begin(), node_ids.end(), std::back_inserter(vector_ids));
        r_sub_model_part.AddNodes(vector_ids);
    }

    ///* Free memory */
    //FreeMemory();

    /* After that we reorder nodes, conditions and elements: */
    ReorderAllIds();

    /* We assign flags and clear the auxiliar model parts created to reassing the flags */
    AssignAndClearAuxiliarSubModelPartForFlags();

    /* Unmoving the original mesh to be able to interpolate */
    if (mFramework == FrameworkEulerLagrange::LAGRANGIAN) {
        NodesArrayType& old_nodes_array = r_old_model_part.Nodes();

        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(old_nodes_array.size()); ++i) {
            auto it_node = old_nodes_array.begin() + i;
            noalias(it_node->Coordinates()) = it_node->GetInitialPosition().Coordinates();
        }
    }

    // We create an auxiliar mesh for debugging purposes
    if (mThisParameters["debug_result_mesh"].GetBool()) {
        CreateDebugPrePostRemeshOutput(r_old_model_part);
    }

    /* We interpolate all the values */
    Parameters interpolate_parameters = Parameters(R"({})" );
    interpolate_parameters.AddValue("echo_level", mThisParameters["echo_level"]);
    interpolate_parameters.AddValue("framework", mThisParameters["framework"]);
    interpolate_parameters.AddValue("max_number_of_searchs", mThisParameters["max_number_of_searchs"]);
    interpolate_parameters.AddValue("step_data_size", mThisParameters["step_data_size"]);
    interpolate_parameters.AddValue("buffer_size", mThisParameters["buffer_size"]);
    interpolate_parameters.AddValue("interpolate_non_historical", mThisParameters["interpolate_non_historical"]);
    interpolate_parameters.AddValue("extrapolate_contour_values", mThisParameters["extrapolate_contour_values"]);
    interpolate_parameters.AddValue("surface_elements", mThisParameters["surface_elements"]);
    interpolate_parameters.AddValue("search_parameters", mThisParameters["search_parameters"]);
    if (TMMGLibrary == MMGLibrary::MMGS) interpolate_parameters["surface_elements"].SetBool(true);
    NodalValuesInterpolationProcess<Dimension> InterpolateNodalValues = NodalValuesInterpolationProcess<Dimension>(r_old_model_part, mrThisModelPart, interpolate_parameters);
    InterpolateNodalValues.Execute();

    /* We initialize elements and conditions */
    if (mThisParameters["initialize_entities"].GetBool())
        InitializeElementsAndConditions();

    /* We do some operations related with the Lagrangian framework */
    if (mFramework == FrameworkEulerLagrange::LAGRANGIAN) {
        // If we remesh during non linear iteration we just move to the previous displacement, to the last displacement otherwise
        const IndexType step = mThisParameters["remesh_at_non_linear_iteration"].GetBool() ? 1 : 0;

        /* We move the mesh */
        r_nodes_array = mrThisModelPart.Nodes();
        const auto it_node_begin = r_nodes_array.begin();

        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i) {
            auto it_node = it_node_begin + i;

            noalias(it_node->Coordinates())  = it_node->GetInitialPosition().Coordinates();
            noalias(it_node->Coordinates()) += it_node->FastGetSolutionStepValue(DISPLACEMENT, step);
        }

        /* We interpolate the internal variables */
        InternalVariablesInterpolationProcess internal_variables_interpolation = InternalVariablesInterpolationProcess(r_old_model_part, mrThisModelPart, mThisParameters["internal_variables_parameters"]);
        internal_variables_interpolation.Execute();
    }

    // We set to zero the variables contained on the elements and conditions
    if (r_old_model_part.Conditions().size() > 0)
        SetToZeroEntityData(mrThisModelPart.Conditions(), r_old_model_part.Conditions());
    if (r_old_model_part.Elements().size() > 0)
        SetToZeroEntityData(mrThisModelPart.Elements(), r_old_model_part.Elements());

    // Finally remove old model part
    owner_model.DeleteModelPart(mrThisModelPart.Name()+"_Old");

    /* We clean conditions with duplicated geometries (this is an error on fluid simulations) */
    if (mFramework == FrameworkEulerLagrange::EULERIAN) {
        ClearConditionsDuplicatedGeometries();
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibrary TMMGLibrary>
void MmgProcess<TMMGLibrary>::ReorderAllIds()
{
    // Iterate over nodes
    NodesArrayType& r_nodes_array = mrThisModelPart.Nodes();
    const auto it_node_begin = r_nodes_array.begin();
    for(IndexType i = 0; i < r_nodes_array.size(); ++i)
        (it_node_begin + i)->SetId(i + 1);

    // Iterate over conditions
    ConditionsArrayType& r_conditions_array = mrThisModelPart.Conditions();
    const auto it_cond_begin = r_conditions_array.begin();
    for(IndexType i = 0; i < r_conditions_array.size(); ++i)
        (it_cond_begin + i)->SetId(i + 1);

    // Iterate over elements
    ElementsArrayType& r_elements_array = mrThisModelPart.Elements();
    const auto it_elem_begin = r_elements_array.begin();
    for(IndexType i = 0; i < r_elements_array.size(); ++i)
        (it_elem_begin + i)->SetId(i + 1);
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibrary TMMGLibrary>
void MmgProcess<TMMGLibrary>::InitializeElementsAndConditions()
{
    // Iterate over conditions
    ConditionsArrayType& r_conditions_array = mrThisModelPart.Conditions();
    const auto it_cond_begin = r_conditions_array.begin();
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(r_conditions_array.size()); ++i)
        (it_cond_begin + i)->Initialize();

    // Iterate over elements
    ElementsArrayType& r_elements_array = mrThisModelPart.Elements();
    const auto it_elem_begin = r_elements_array.begin();
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(r_elements_array.size()); ++i)
        (it_elem_begin + i)->Initialize();
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibrary TMMGLibrary>
void MmgProcess<TMMGLibrary>::SaveSolutionToFile(const bool PostOutput)
{
    /* GET RESULTS */

    const IndexType step = mrThisModelPart.GetProcessInfo()[STEP];

    // Automatically save the mesh
    mMmmgUtilities.OutputMesh(mStdStringFilename, PostOutput, step);

    // Automatically save the solution
    mMmmgUtilities.OutputSol(mStdStringFilename, PostOutput, step);

    // Save the mesh in an .mdpa format
    const bool save_mdpa_file = mThisParameters["save_mdpa_file"].GetBool();
    if(save_mdpa_file) OutputMdpa();
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibrary TMMGLibrary>
void MmgProcess<TMMGLibrary>::FreeMemory()
{
    // Free the MMG structures
    mMmmgUtilities.FreeAll();

    // Free filename (NOTE: Problems with more that one iteration)
//     free(mFilename);
//     mFilename = nullptr;

    // Free reference std::unordered_map
    mpRefElement.clear();
    mpRefCondition.clear();
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibrary TMMGLibrary>
void MmgProcess<TMMGLibrary>::OutputMdpa()
{
    std::ofstream output_file;
    ModelPartIO model_part_io("output", IO::WRITE);
    model_part_io.WriteModelPart(mrThisModelPart);
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibrary TMMGLibrary>
void MmgProcess<TMMGLibrary>::CreateAuxiliarSubModelPartForFlags()
{
    ModelPart& r_auxiliar_model_part = mrThisModelPart.CreateSubModelPart("AUXILIAR_MODEL_PART_TO_LATER_REMOVE");

    const auto& flags = KratosComponents<Flags>::GetComponents();

    for (auto& flag : flags) {
        const std::string name_sub_model = "FLAG_"+flag.first;
        if (name_sub_model.find("NOT") == std::string::npos) { // Avoiding inactive flags
            r_auxiliar_model_part.CreateSubModelPart(name_sub_model);
            ModelPart& auxiliar_sub_model_part = r_auxiliar_model_part.GetSubModelPart(name_sub_model);
            FastTransferBetweenModelPartsProcess transfer_process = FastTransferBetweenModelPartsProcess(auxiliar_sub_model_part, mrThisModelPart, FastTransferBetweenModelPartsProcess::EntityTransfered::ALL, *(flag.second));
            transfer_process.Execute();
            // If the number of elements transfered is 0 we remove the model part
            if (auxiliar_sub_model_part.NumberOfNodes() == 0
            && auxiliar_sub_model_part.NumberOfElements() == 0
            && auxiliar_sub_model_part.NumberOfConditions() == 0) {
                r_auxiliar_model_part.RemoveSubModelPart(name_sub_model);
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibrary TMMGLibrary>
void MmgProcess<TMMGLibrary>::AssignAndClearAuxiliarSubModelPartForFlags()
{
    const auto& flags = KratosComponents<Flags>::GetComponents();

    ModelPart& auxiliar_model_part = mrThisModelPart.GetSubModelPart("AUXILIAR_MODEL_PART_TO_LATER_REMOVE");
    for (auto& flag : flags) {
        const std::string name_sub_model = "FLAG_"+flag.first;
        if (auxiliar_model_part.HasSubModelPart(name_sub_model)) {
            ModelPart& auxiliar_sub_model_part = auxiliar_model_part.GetSubModelPart(name_sub_model);
            VariableUtils().SetFlag(*(flag.second), true, auxiliar_sub_model_part.Nodes());
            VariableUtils().SetFlag(*(flag.second), true, auxiliar_sub_model_part.Conditions());
            VariableUtils().SetFlag(*(flag.second), true, auxiliar_sub_model_part.Elements());
        }
    }

    mrThisModelPart.RemoveSubModelPart("AUXILIAR_MODEL_PART_TO_LATER_REMOVE");
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibray TMMGLibray>
void MmgProcess<TMMGLibray>::ClearConditionsDuplicatedGeometries()
{
    // Next check that the conditions are oriented accordingly to do so begin by putting all of the conditions in a set
    typedef std::unordered_map<DenseVector<IndexType>, std::vector<IndexType>, KeyHasherRange<DenseVector<IndexType>>, KeyComparorRange<DenseVector<IndexType>> > HashMapType;
    HashMapType faces_map;

    // Iterate over conditions
    ConditionsArrayType& r_conditions_array = mrThisModelPart.Conditions();

    // Reset flag
    const auto it_cond_begin = r_conditions_array.begin();
    const int number_of_conditions = static_cast<int>(r_conditions_array.size());
    #pragma omp parallel for
    for(int i = 0; i < number_of_conditions; ++i) {
        const auto it_cond = it_cond_begin + i;
        it_cond->Reset(TO_ERASE);
    }

    // Create map
    for(auto& r_cond : r_conditions_array) {

        GeometryType& r_geom = r_cond.GetGeometry();
        DenseVector<IndexType> ids(r_geom.size());

        for(IndexType i = 0; i < ids.size(); ++i) {
            ids[i] = r_geom[i].Id();
        }

        //*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
        std::sort(ids.begin(), ids.end());

        // Insert a pointer to the condition identified by the hash value ids
        HashMapType::iterator it_face = faces_map.find(ids);
        if(it_face != faces_map.end() ) { // Already defined vector
            (it_face->second).push_back(r_cond.Id());
        } else {
            std::vector<IndexType> aux_cond_id(1);
            aux_cond_id[0] = r_cond.Id();
            faces_map.insert( HashMapType::value_type(std::pair<DenseVector<IndexType>, std::vector<IndexType>>({ids, aux_cond_id})) );
        }
    }

    // We set the flag
    SizeType counter = 1;
    for (auto& i_pair : faces_map) {
        const auto& r_pairs = i_pair.second;
        for (auto i_id : r_pairs) {
            auto p_cond = mrThisModelPart.pGetCondition(i_id);
            if (p_cond->Is(MARKER) && counter < r_pairs.size()) { // Only remove dummy conditions repeated
                p_cond->Set(TO_ERASE);
                KRATOS_INFO_IF("MmgProcess", mEchoLevel > 2) << "Condition created ID:\t" << i_id << " will be removed" << std::endl;
                ++counter;
            }
            counter = 1;
        }
    }

    // We remove the conditions marked to be removed
    mrThisModelPart.RemoveConditionsFromAllLevels(TO_ERASE);
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibrary TMMGLibrary>
void MmgProcess<TMMGLibrary>::CreateDebugPrePostRemeshOutput(ModelPart& rOldModelPart)
{
    Model& owner_model = mrThisModelPart.GetModel();
    ModelPart& auxiliar_model_part = owner_model.CreateModelPart(mrThisModelPart.Name()+"_Auxiliar", mrThisModelPart.GetBufferSize());
    ModelPart& copy_old_model_part = owner_model.CreateModelPart(mrThisModelPart.Name()+"_Old_Copy", mrThisModelPart.GetBufferSize());

    Properties::Pointer p_prop_1 = auxiliar_model_part.pGetProperties(1);
    Properties::Pointer p_prop_2 = auxiliar_model_part.pGetProperties(2);

    // We just transfer nodes and elements
    // Current model part
    FastTransferBetweenModelPartsProcess transfer_process_current = FastTransferBetweenModelPartsProcess(auxiliar_model_part, mrThisModelPart, FastTransferBetweenModelPartsProcess::EntityTransfered::NODESANDELEMENTS);
    transfer_process_current.Set(MODIFIED); // We replicate, not transfer
    transfer_process_current.Execute();

    ElementsArrayType& elements_array_1 = auxiliar_model_part.Elements();

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(elements_array_1.size()); ++i) {
        auto it_elem = elements_array_1.begin() + i;
        it_elem->SetProperties(p_prop_1);
    }
    // Old model part
    FastTransferBetweenModelPartsProcess transfer_process_old = FastTransferBetweenModelPartsProcess(copy_old_model_part, rOldModelPart, FastTransferBetweenModelPartsProcess::EntityTransfered::NODESANDELEMENTS);
    transfer_process_current.Set(MODIFIED); // We replicate, not transfer
    transfer_process_old.Execute();

    ElementsArrayType& elements_array_2 = copy_old_model_part.Elements();

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(elements_array_2.size()); ++i) {
        auto it_elem = elements_array_2.begin() + i;
        it_elem->SetProperties(p_prop_2);
    }

    // Reorder ids to ensure be consecuent
    NodesArrayType& auxiliar_nodes_array = auxiliar_model_part.Nodes();
    const SizeType auxiliar_number_of_nodes = (auxiliar_nodes_array.end() - 1)->Id();
    NodesArrayType& copy_old_nodes_array = copy_old_model_part.Nodes();

    for(IndexType i = 0; i < copy_old_nodes_array.size(); ++i) {
        auto it_node = copy_old_nodes_array.begin() + i;
        it_node->SetId(auxiliar_number_of_nodes + i + 1);
    }

    // Last transfer
    FastTransferBetweenModelPartsProcess transfer_process_last = FastTransferBetweenModelPartsProcess(auxiliar_model_part, copy_old_model_part, FastTransferBetweenModelPartsProcess::EntityTransfered::NODESANDELEMENTS);
    transfer_process_last.Set(MODIFIED);
    transfer_process_last.Execute();

    const int step = mrThisModelPart.GetProcessInfo()[STEP];
    const double label = static_cast<double>(step);
    GidIO<> gid_io("BEFORE_AND_AFTER_MMG_MESH_STEP=" + std::to_string(step), GiD_PostBinary, SingleFile, WriteUndeformed,  WriteElementsOnly);

    gid_io.InitializeMesh(label);
    gid_io.WriteMesh(auxiliar_model_part.GetMesh());
    gid_io.FinalizeMesh();
    gid_io.InitializeResults(label, auxiliar_model_part.GetMesh());

    // Remove auxiliar model parts
    owner_model.DeleteModelPart(mrThisModelPart.Name()+"_Auxiliar");
    owner_model.DeleteModelPart(mrThisModelPart.Name()+"_Old_Copy");
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibrary TMMGLibrary>
void MmgProcess<TMMGLibrary>::CleanSuperfluousNodes()
{
    const int initial_num = mrThisModelPart.Nodes().size();

    auto& r_nodes_array = mrThisModelPart.Nodes();
    const auto& r_elem_array = mrThisModelPart.Elements();

    // marking all nodes as "superfluous"
    #pragma omp parallel for
    for( int i_node = 0; i_node < static_cast<int>(r_nodes_array.size()); ++i_node ){

        auto node = r_nodes_array.begin() + i_node;
        node->Set(TO_ERASE, true);
    }

    // saving the nodes that belong to an element
    #pragma omp parallel for
    for( int i_elem = 0; i_elem < static_cast<int>(r_elem_array.size()); ++i_elem ){

        const auto elem = r_elem_array.begin() + i_elem;
        auto& r_geom = elem->GetGeometry();

        for (unsigned int i = 0; i < r_geom.size(); ++i){
            r_geom[i].Set(TO_ERASE, false);
        }
    }

    mrThisModelPart.RemoveNodesFromAllLevels(TO_ERASE);
    const int final_num = mrThisModelPart.Nodes().size();
    KRATOS_INFO("MmgProcess") << "In total " << (initial_num - final_num) <<" superfluous nodes were cleared" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibrary TMMGLibrary>
Parameters MmgProcess<TMMGLibrary>::GetDefaultParameters()
{
    Parameters default_parameters = Parameters(R"(
    {
        "filename"                             : "out",
        "discretization_type"                  : "Standard",
        "isosurface_parameters"                :
        {
            "isosurface_variable"              : "DISTANCE",
            "nonhistorical_variable"           : false,
            "remove_regions"                   : false
        },
        "framework"                            : "Eulerian",
        "internal_variables_parameters"        :
        {
            "allocation_size"                      : 1000,
            "bucket_size"                          : 4,
            "search_factor"                        : 2,
            "interpolation_type"                   : "LST",
            "internal_variable_interpolation_list" :[]
        },
        "force_sizes"                          :
        {
            "force_min"                           : false,
            "minimal_size"                        : 0.1,
            "force_max"                           : false,
            "maximal_size"                        : 10.0
        },
        "advanced_parameters"                  :
        {
            "force_hausdorff_value"               : false,
            "hausdorff_value"                     : 0.0001,
            "no_move_mesh"                        : false,
            "no_surf_mesh"                        : false,
            "no_insert_mesh"                      : false,
            "no_swap_mesh"                        : false,
            "normal_regularization_mesh"          : false,
            "deactivate_detect_angle"             : false,
            "force_gradation_value"               : false,
            "gradation_value"                     : 1.3
        },
        "save_external_files"                  : false,
        "save_mdpa_file"                       : false,
        "max_number_of_searchs"                : 1000,
        "interpolate_non_historical"           : true,
        "extrapolate_contour_values"           : true,
        "surface_elements"                     : false,
        "search_parameters"                    : {
            "allocation_size"                     : 1000,
            "bucket_size"                         : 4,
            "search_factor"                       : 2.0
        },
        "echo_level"                           : 3,
        "debug_result_mesh"                    : false,
        "step_data_size"                       : 0,
        "initialize_entities"                  : true,
        "remesh_at_non_linear_iteration"       : false,
        "buffer_size"                          : 0
    })" );

    return default_parameters;
}

/***********************************************************************************/
/***********************************************************************************/

template class MmgProcess<MMGLibrary::MMG2D>;
template class MmgProcess<MMGLibrary::MMG3D>;
template class MmgProcess<MMGLibrary::MMGS>;

}// namespace Kratos.
