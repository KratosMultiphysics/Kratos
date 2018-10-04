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
#include <unordered_set>

// External includes
// The includes related with the MMG library
#include "mmg/libmmg.h"
#include "mmg/mmg2d/libmmg2d.h"
#include "mmg/mmg3d/libmmg3d.h"
#include "mmg/mmgs/libmmgs.h"

// Project includes
#include "custom_processes/mmg_process.h"
#include "utilities/sub_model_parts_list_utility.h"
#include "utilities/variable_utils.h"
// We indlude the internal variable interpolation process
#include "custom_processes/nodal_values_interpolation_process.h"
#include "custom_processes/internal_variables_interpolation_process.h"
#include "processes/fast_transfer_between_model_parts_process.h"
// Include the point locator
#include "utilities/binbased_fast_point_locator.h"
// Include the spatial containers needed for search
#include "spatial_containers/spatial_containers.h" // kd-tree
#include "includes/io.h"
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
// The member variables related with the MMG library
MMG5_pMesh mmgMesh;
MMG5_pSol  mmgSol;

/************************************* CONSTRUCTOR *********************************/
/***********************************************************************************/

template<MMGLibray TMMGLibray>
MmgProcess<TMMGLibray>::MmgProcess(
    ModelPart& rThisModelPart,
    Parameters ThisParameters
    ):mrThisModelPart(rThisModelPart),
      mThisParameters(ThisParameters)
{
    Parameters DefaultParameters = Parameters(R"(
    {
        "filename"                             : "out",
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
            "deactivate_detect_angle"             : false,
            "force_gradation_value"               : false,
            "gradation_value"                     : 1.3
        },
        "save_external_files"                  : false,
        "save_mdpa_file"                       : false,
        "max_number_of_searchs"                : 1000,
        "interpolate_non_historical"           : true,
        "extrapolate_contour_values"           : true,
        "search_parameters"                    : {
            "allocation_size"                     : 1000,
            "bucket_size"                         : 4,
            "search_factor"                       : 2.0
        },
        "echo_level"                           : 3,
        "debug_result_mesh"                    : false,
        "step_data_size"                       : 0,
        "remesh_at_non_linear_iteration"       : false,
        "buffer_size"                          : 0
    })" );

    mThisParameters.ValidateAndAssignDefaults(DefaultParameters);

    mStdStringFilename = mThisParameters["filename"].GetString();
    mEchoLevel = mThisParameters["echo_level"].GetInt();

    mFilename = new char [mStdStringFilename.length() + 1];
    std::strcpy (mFilename, mStdStringFilename.c_str());

    mFramework = ConvertFramework(mThisParameters["framework"].GetString());

    mpRefElement.clear();
    mpRefCondition.clear();
}

/*************************************** EXECUTE ***********************************/
/***********************************************************************************/

template<MMGLibray TMMGLibray>
void MmgProcess<TMMGLibray>::Execute()
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

template<MMGLibray TMMGLibray>
void MmgProcess<TMMGLibray>::ExecuteInitialize()
{
    KRATOS_TRY;

    /* We restart the MMG mesh and solution */
    InitMesh();

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibray TMMGLibray>
void MmgProcess<TMMGLibray>::ExecuteBeforeSolutionLoop()
{
    KRATOS_TRY;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibray TMMGLibray>
void MmgProcess<TMMGLibray>::ExecuteInitializeSolutionStep()
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
    InitializeSolData();

    // Check if the number of given entities match with mesh size
    CheckMeshData();

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

template<MMGLibray TMMGLibray>
void MmgProcess<TMMGLibray>::ExecuteFinalizeSolutionStep()
{
    KRATOS_TRY;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibray TMMGLibray>
void MmgProcess<TMMGLibray>::ExecuteBeforeOutputStep()
{
    KRATOS_TRY;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibray TMMGLibray>
void MmgProcess<TMMGLibray>::ExecuteAfterOutputStep()
{
    KRATOS_TRY;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibray TMMGLibray>
void MmgProcess<TMMGLibray>::ExecuteFinalize()
{
    KRATOS_TRY;

    // We release the memory
    FreeMemory();

    KRATOS_CATCH("");
}

/************************************* OPERATOR() **********************************/
/***********************************************************************************/

template<MMGLibray TMMGLibray>
void MmgProcess<TMMGLibray>::operator()()
{
    Execute();
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibray TMMGLibray>
void MmgProcess<TMMGLibray>::InitializeMeshData()
{
    // We create a list of submodelparts to later reassign flags after remesh
    CreateAuxiliarSubModelPartForFlags();

    // First we compute the colors
    mColors.clear();
    ColorsMapType nodes_colors, cond_colors, elem_colors;
    SubModelPartsListUtility sub_model_parts_list(mrThisModelPart);
    sub_model_parts_list.ComputeSubModelPartsList(nodes_colors, cond_colors, elem_colors, mColors);

    /////////* MESH FILE */////////
    // Build mesh in MMG5 format //

    // Iterate over components
    NodesArrayType& nodes_array = mrThisModelPart.Nodes();
    ConditionsArrayType& conditions_array = mrThisModelPart.Conditions();
    ElementsArrayType& elements_array = mrThisModelPart.Elements();

    /* Manually set of the mesh */
    array_1d<SizeType, ConditionsArraySize> num_array_conditions;
    array_1d<SizeType, ElementsArraySize> num_array_elements;
    if (TMMGLibray == MMGLibray::MMG2D) { // 2D
        num_array_conditions[0] = conditions_array.size();
        num_array_elements[0]   = elements_array.size();
    } else if (TMMGLibray == MMGLibray::MMG3D) { // 3D
        /* Elements */
        std::size_t num_tetra = 0, num_prisms = 0;
        #pragma omp parallel for reduction(+:num_tetra,num_prisms)
        for(int i = 0; i < static_cast<int>(elements_array.size()); ++i) {
            auto it_elem = elements_array.begin() + i;

            if ((it_elem->GetGeometry()).GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4) { // Tetrahedron
                num_tetra += 1;
            } else if ((it_elem->GetGeometry()).GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Prism3D6) { // Prisms
                num_prisms += 1;
            } else
                KRATOS_WARNING("MmgProcess") << "WARNING: YOUR GEOMETRY CONTAINS " << it_elem->GetGeometry().PointsNumber() <<" NODES CAN NOT BE REMESHED" << std::endl;
        }

        num_array_elements[0] = num_tetra;  // Tetrahedron
        num_array_elements[1] = num_prisms; // Prisms

        KRATOS_INFO_IF("MmgProcess", ((num_tetra + num_tetra) < elements_array.size()) && mEchoLevel > 0) <<
        "Number of Elements: " << elements_array.size() << " Number of Tetrahedron: " << num_tetra << " Number of Prisms: " << num_tetra << std::endl;

        /* Conditions */
        std::size_t num_tri = 0, num_quad = 0;
        #pragma omp parallel for reduction(+:num_tri,num_quad)
        for(int i = 0; i < static_cast<int>(conditions_array.size()); ++i) {
            auto it_cond = conditions_array.begin() + i;

            if ((it_cond->GetGeometry()).GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle3D3) { // Triangles
                num_tri += 1;
            } else if ((it_cond->GetGeometry()).GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4) { // Quadrilaterals
                num_quad += 1;
            } else
                KRATOS_WARNING("MmgProcess") << "WARNING: YOUR GEOMETRY CONTAINS CERTAIN TYPE THAT CAN NOT BE REMESHED" << std::endl;
        }

        num_array_conditions[0] = num_tri;  // Triangles
        num_array_conditions[1] = num_quad; // Quadrilaterals

        KRATOS_INFO_IF("MmgProcess", ((num_tri + num_quad) < conditions_array.size()) && mEchoLevel > 0) <<
        "Number of Conditions: " << conditions_array.size() << " Number of Triangles: " << num_tri << " Number of Quadrilaterals: " << num_quad << std::endl;
    } else { // Surfaces
        num_array_conditions[0] = conditions_array.size();
        num_array_elements[0]   = elements_array.size();
    }

    SetMeshSize(nodes_array.size(), num_array_elements, num_array_conditions);

    /* Nodes */
    // We copy the DOF from the first node (after we release, to avoid problem with previous conditions)
    mDofs = nodes_array.begin()->GetDofs();
    for (typename NodeType::DofsContainerType::const_iterator it_dof = mDofs.begin(); it_dof != mDofs.end(); ++it_dof)
        it_dof->FreeDof();

    if (mFramework == FrameworkEulerLagrange::LAGRANGIAN){ // NOTE: The code is repeated due to performance reasons
        #pragma omp parallel for firstprivate(nodes_colors)
        for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) {
            auto it_node = nodes_array.begin() + i;

            SetNodes(it_node->X0(), it_node->Y0(), it_node->Z0(), nodes_colors[it_node->Id()], i + 1);

            bool blocked = false;
            if (it_node->IsDefined(BLOCKED))
                blocked = it_node->Is(BLOCKED);
            if (Dimension == 3 && blocked)
                BlockNode(i + 1);

            // RESETING THE ID OF THE NODES (important for non consecutive meshes)
            it_node->SetId(i + 1);
        }
    }
    else {
        #pragma omp parallel for firstprivate(nodes_colors)
        for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) {
            auto it_node = nodes_array.begin() + i;

            SetNodes(it_node->X(), it_node->Y(), it_node->Z(), nodes_colors[it_node->Id()], i + 1);

            bool blocked = false;
            if (it_node->IsDefined(BLOCKED))
                blocked = it_node->Is(BLOCKED);
            if (Dimension == 3 && blocked)
                BlockNode(i + 1);

            // RESETING THE ID OF THE NODES (important for non consecutive meshes)
            it_node->SetId(i + 1);
        }
    }

    /* Conditions */
    #pragma omp parallel for firstprivate(cond_colors)
    for(int i = 0; i < static_cast<int>(conditions_array.size()); ++i)  {
        auto it_cond = conditions_array.begin() + i;
        SetConditions(it_cond->GetGeometry(), cond_colors[it_cond->Id()], i + 1);
    }

    /* Elements */
    #pragma omp parallel for firstprivate(elem_colors)
    for(int i = 0; i < static_cast<int>(elements_array.size()); ++i) {
        auto it_elem = elements_array.begin() + i;
        SetElements(it_elem->GetGeometry(), elem_colors[it_elem->Id()], i + 1);
    }

    // Create auxiliar colors maps
    ColorsMapType aux_ref_cond;
    for(int i = 0; i < static_cast<int>(conditions_array.size()); ++i)  {
        auto it_cond = conditions_array.begin() + i;
        const IndexType cond_id = it_cond->Id();
        const IndexType color = cond_colors[cond_id];
        if (!(aux_ref_cond.find(color) != aux_ref_cond.end()))
            aux_ref_cond.insert (IndexPairType(color,cond_id));
    }
    ColorsMapType aux_ref_elem;
    for(int i = 0; i < static_cast<int>(elements_array.size()); ++i) {
        auto it_elem = elements_array.begin() + i;
        const IndexType elem_id = it_elem->Id();
        const IndexType color = elem_colors[elem_id];
        if (!(aux_ref_elem.find(color) != aux_ref_elem.end()))
            aux_ref_elem.insert (IndexPairType(color,elem_id));
    }

    /* We clone the first condition and element of each type (we will assume that each sub model part has just one kind of condition, in my opinion it is quite reccomended to create more than one sub model part if you have more than one element or condition) */
    // First we add the main model part
    if (conditions_array.size() > 0) {
        const std::string type_name = (Dimension == 2) ? "Condition2D2N" : (TMMGLibray == MMGLibray::MMG3D) ? "Condition3D" : "Condition3D2N";
        Condition const& r_clone_condition = KratosComponents<Condition>::Get(type_name);
        mpRefCondition[0] = r_clone_condition.Create(0, r_clone_condition.GetGeometry(), conditions_array.begin()->pGetProperties());
    }
    if (elements_array.size() > 0) {
        mpRefElement[0] = elements_array.begin()->Create(0, elements_array.begin()->GetGeometry(), elements_array.begin()->pGetProperties());
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

template<MMGLibray TMMGLibray>
void MmgProcess<TMMGLibray>::InitializeSolData()
{
    ////////* SOLUTION FILE *////////

    // Iterate in the nodes
    NodesArrayType& nodes_array = mrThisModelPart.Nodes();

    // Set size of the solution
    SetSolSizeTensor(nodes_array.size());

    const Variable<TensorArrayType>& tensor_variable = KratosComponents<Variable<TensorArrayType>>::Get("METRIC_TENSOR_"+std::to_string(Dimension)+"D");

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) {
        auto it_node = nodes_array.begin() + i;

        KRATOS_DEBUG_ERROR_IF_NOT(it_node->Has(tensor_variable)) << "METRIC_TENSOR_" + std::to_string(Dimension) + "D  not defined for node " << it_node->Id() << std::endl;

        // We get the metric
        const TensorArrayType& metric = it_node->GetValue(tensor_variable);

        // We set the metric
        SetMetricTensor(metric, i + 1);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibray TMMGLibray>
void MmgProcess<TMMGLibray>::ExecuteRemeshing()
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

    MMGLibCall();

    const SizeType number_of_nodes = mmgMesh->np;
    array_1d<SizeType, 2> n_conditions;
    if (TMMGLibray == MMGLibray::MMG2D) { // 2D
        n_conditions[0] = mmgMesh->na;
        n_conditions[1] = 0;
    } else if (TMMGLibray == MMGLibray::MMG3D) { // 3D
        n_conditions[0] = mmgMesh->nt;
        n_conditions[1] = mmgMesh->nquad;
    }
    array_1d<SizeType, 2> n_elements;
    if (TMMGLibray == MMGLibray::MMG2D) { // 2D
        n_elements[0] = mmgMesh->nt;
        n_elements[1] = 0;
    } else if (TMMGLibray == MMGLibray::MMG3D) { // 3D
        n_elements[0] = mmgMesh->ne;
        n_elements[1] = mmgMesh->nprism;
    }

    KRATOS_INFO_IF("MmgProcess", mEchoLevel > 0) << "\tNodes created: " << number_of_nodes << std::endl;
    if (TMMGLibray == MMGLibray::MMG2D) { // 2D
        KRATOS_INFO_IF("MmgProcess", mEchoLevel > 0) <<
        "Conditions created: " << n_conditions[0] << "\n" <<
        "Elements created: " << n_elements[0] << std::endl;
    } else if (TMMGLibray == MMGLibray::MMG3D) { // 3D
        KRATOS_INFO_IF("MmgProcess", mEchoLevel > 0) <<
        "Conditions created: " << n_conditions[0] + n_conditions[1] << "\n\tTriangles: " << n_conditions[0] << "\tQuadrilaterals: " << n_conditions[1] << "\n" <<
        "Elements created: " << n_elements[0] + n_elements[1] << "\n\tTetrahedron: " << n_elements[0] << "\tPrisms: " << n_elements[1] << std::endl;
    } else { // Surfaces
    }

    ////////* EMPTY AND BACKUP THE MODEL PART *////////

    ModelPart r_old_model_part;

    // First we empty the model part
    NodesArrayType& nodes_array = mrThisModelPart.Nodes();

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i)
        (nodes_array.begin() + i)->Set(TO_ERASE, true);
    r_old_model_part.AddNodes( mrThisModelPart.NodesBegin(), mrThisModelPart.NodesEnd() );
    mrThisModelPart.RemoveNodesFromAllLevels(TO_ERASE);

    ConditionsArrayType& conditions_array = mrThisModelPart.Conditions();

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(conditions_array.size()); ++i)
        (conditions_array.begin() + i)->Set(TO_ERASE, true);
    r_old_model_part.AddConditions( mrThisModelPart.ConditionsBegin(), mrThisModelPart.ConditionsEnd() );
    mrThisModelPart.RemoveConditionsFromAllLevels(TO_ERASE);

    ElementsArrayType& elements_array = mrThisModelPart.Elements();

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(elements_array.size()); ++i)
        (elements_array.begin() + i)->Set(TO_ERASE, true);
    r_old_model_part.AddElements( mrThisModelPart.ElementsBegin(), mrThisModelPart.ElementsEnd() );
    mrThisModelPart.RemoveElementsFromAllLevels(TO_ERASE);

    // Create a new model part // TODO: Use a different kind of element for each submodelpart (in order to be able of remeshing more than one kind o element or condition)
    std::unordered_map<IndexType, IndexVectorType> color_nodes, color_cond_0, color_cond_1, color_elem_0, color_elem_1;

    // The tempotal store of
    ConditionsArrayType created_conditions_vector;
    ElementsArrayType created_elements_vector;

    // Auxiliar values
    int ref, is_required;

    /* NODES */ // TODO: ADD OMP
    for (IndexType i_node = 1; i_node <= number_of_nodes; ++i_node) {
        NodeType::Pointer p_node = CreateNode(i_node, ref, is_required);

        // Set the DOFs in the nodes
        for (typename NodeType::DofsContainerType::const_iterator it_dof = mDofs.begin(); it_dof != mDofs.end(); ++it_dof)
            p_node->pAddDof(*it_dof);

        if (ref != 0) color_nodes[static_cast<IndexType>(ref)].push_back(i_node);// NOTE: ref == 0 is the MainModelPart
    }

    /* CONDITIONS */ // TODO: ADD OMP
    if (mpRefCondition.size() > 0) {
        IndexType cond_id = 1;

        IndexType counter_cond_0 = 0;
        const IndexVectorType condition_to_remove_0 = CheckConditions0();
        for (IndexType i_cond = 1; i_cond <= n_conditions[0]; ++i_cond) {
            bool skip_creation = false;
            if (counter_cond_0 < condition_to_remove_0.size()) {
                if (condition_to_remove_0[counter_cond_0] == i_cond) {
                    skip_creation = true;
                    counter_cond_0 += 1;
                }
            }
            ConditionType::Pointer p_condition = CreateCondition0(cond_id, ref, is_required, skip_creation);

            if (p_condition != nullptr) {
                created_conditions_vector.push_back(p_condition);
//                 mrThisModelPart.AddCondition(p_condition);
                if (ref != 0) color_cond_0[static_cast<IndexType>(ref)].push_back(cond_id);// NOTE: ref == 0 is the MainModelPart
                cond_id += 1;
            }
        }

        IndexType counter_cond_1 = 0;
        const IndexVectorType condition_to_remove_1 = CheckConditions1();
        for (IndexType i_cond = 1; i_cond <= n_conditions[1]; ++i_cond)
        {
            bool skip_creation = false;
            if (counter_cond_1 < condition_to_remove_1.size()) {
                if (condition_to_remove_1[counter_cond_1] == i_cond) {
                    skip_creation = true;
                    counter_cond_1 += 1;
                }
            }
            ConditionType::Pointer p_condition = CreateCondition1(cond_id, ref, is_required, skip_creation);

            if (p_condition != nullptr) {
                created_conditions_vector.push_back(p_condition);
//                 mrThisModelPart.AddCondition(p_condition);
                if (ref != 0) color_cond_1[static_cast<IndexType>(ref)].push_back(cond_id);// NOTE: ref == 0 is the MainModelPart
                cond_id += 1;
            }
        }
    }

    /* ELEMENTS */ // TODO: ADD OMP
    if (mpRefElement.size() > 0) {
        IndexType elem_id = 1;

        IndexType counter_elem_0 = 0;
        const IndexVectorType elements_to_remove_0 = CheckElements0();
        for (IndexType i_elem = 1; i_elem <= n_elements[0]; ++i_elem) {
            bool skip_creation = false;
            if (counter_elem_0 < elements_to_remove_0.size()) {
                if (elements_to_remove_0[counter_elem_0] == i_elem) {
                    skip_creation = true;
                    counter_elem_0 += 1;
                }
            }

            ElementType::Pointer p_element = CreateElement0(elem_id, ref, is_required, skip_creation);

            if (p_element != nullptr) {
                created_elements_vector.push_back(p_element);
//                 mrThisModelPart.AddElement(p_element);
                if (ref != 0) color_elem_0[static_cast<IndexType>(ref)].push_back(elem_id);// NOTE: ref == 0 is the MainModelPart
                elem_id += 1;
            }
        }

        IndexType counter_elem_1 = 0;
        const IndexVectorType elements_to_remove_1 = CheckElements1();
        for (IndexType i_elem = 1; i_elem <= n_elements[1]; ++i_elem) {
            bool skip_creation = false;
            if (counter_elem_1 < elements_to_remove_1.size()) {
                if (elements_to_remove_1[counter_elem_1] == i_elem) {
                    skip_creation = true;
                    counter_elem_1 += 1;
                }
            }

            ElementType::Pointer p_element = CreateElement1(elem_id, ref, is_required,skip_creation);

            if (p_element != nullptr) {
                created_elements_vector.push_back(p_element);
//                 mrThisModelPart.AddElement(p_element);
                if (ref != 0) color_elem_1[static_cast<IndexType>(ref)].push_back(elem_id);// NOTE: ref == 0 is the MainModelPart
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
                ModelPart& r_sub_model_part = SubModelPartsListUtility::GetRecursiveSubModelPart(mrThisModelPart, sub_model_part_name);

                if (color_nodes.find(key) != color_nodes.end()) r_sub_model_part.AddNodes(color_nodes[key]);
                if (color_cond_0.find(key) != color_cond_0.end()) r_sub_model_part.AddConditions(color_cond_0[key]);
                if (color_cond_1.find(key) != color_cond_1.end()) r_sub_model_part.AddConditions(color_cond_1[key]);
                if (color_elem_0.find(key) != color_elem_0.end()) r_sub_model_part.AddElements(color_elem_0[key]);
                if (color_elem_1.find(key) != color_elem_1.end()) r_sub_model_part.AddElements(color_elem_1[key]);
            }
        }
    }

    // TODO: Add OMP
    // NOTE: We add the nodes from the elements and conditions to the respective submodelparts
    const std::vector<std::string> sub_model_part_names = SubModelPartsListUtility::GetRecursiveSubModelPartNames(mrThisModelPart);

    for (auto sub_model_part_name : sub_model_part_names) {
        ModelPart& r_sub_model_part = SubModelPartsListUtility::GetRecursiveSubModelPart(mrThisModelPart, sub_model_part_name);

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

    /* Save to file */
    if (save_to_file) SaveSolutionToFile(true);

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
    Parameters InterpolateParameters = Parameters(R"({})" );
    InterpolateParameters.AddValue("echo_level", mThisParameters["echo_level"]);
    InterpolateParameters.AddValue("framework", mThisParameters["framework"]);
    InterpolateParameters.AddValue("max_number_of_searchs", mThisParameters["max_number_of_searchs"]);
    InterpolateParameters.AddValue("step_data_size", mThisParameters["step_data_size"]);
    InterpolateParameters.AddValue("buffer_size", mThisParameters["buffer_size"]);
    InterpolateParameters.AddValue("interpolate_non_historical", mThisParameters["interpolate_non_historical"]);
    InterpolateParameters.AddValue("extrapolate_contour_values", mThisParameters["extrapolate_contour_values"]);
    InterpolateParameters.AddValue("search_parameters", mThisParameters["search_parameters"]);
    NodalValuesInterpolationProcess<Dimension> InterpolateNodalValues = NodalValuesInterpolationProcess<Dimension>(r_old_model_part, mrThisModelPart, InterpolateParameters);
    InterpolateNodalValues.Execute();

    /* We initialize elements and conditions */
    InitializeElementsAndConditions();

    /* We do some operations related with the Lagrangian framework */
    if (mFramework == FrameworkEulerLagrange::LAGRANGIAN) {
        // If we remesh during non linear iteration we just move to the previous displacement, to the last displacement otherwise
        const IndexType step = mThisParameters["remesh_at_non_linear_iteration"].GetBool() ? 1 : 0;

        /* We move the mesh */
        nodes_array = mrThisModelPart.Nodes();

        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) {
            auto it_node = nodes_array.begin() + i;

            noalias(it_node->Coordinates())  = it_node->GetInitialPosition().Coordinates();
            noalias(it_node->Coordinates()) += it_node->FastGetSolutionStepValue(DISPLACEMENT, step);
        }

        /* We interpolate the internal variables */
        InternalVariablesInterpolationProcess InternalVariablesInterpolation = InternalVariablesInterpolationProcess(r_old_model_part, mrThisModelPart, mThisParameters["internal_variables_parameters"]);
        InternalVariablesInterpolation.Execute();
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibray TMMGLibray>
void MmgProcess<TMMGLibray>::ReorderAllIds()
{
    NodesArrayType& nodes_array = mrThisModelPart.Nodes();
    for(SizeType i = 0; i < nodes_array.size(); ++i)
        (nodes_array.begin() + i)->SetId(i + 1);

    ConditionsArrayType& condition_array = mrThisModelPart.Conditions();
    for(SizeType i = 0; i < condition_array.size(); ++i)
        (condition_array.begin() + i)->SetId(i + 1);

    ElementsArrayType& element_array = mrThisModelPart.Elements();
    for(SizeType i = 0; i < element_array.size(); ++i)
        (element_array.begin() + i)->SetId(i + 1);
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibray TMMGLibray>
void MmgProcess<TMMGLibray>::InitializeElementsAndConditions()
{
    ConditionsArrayType& condition_array = mrThisModelPart.Conditions();
    for(SizeType i = 0; i < condition_array.size(); ++i)
        (condition_array.begin() + i)->Initialize();

    ElementsArrayType& element_array = mrThisModelPart.Elements();
    for(SizeType i = 0; i < element_array.size(); ++i)
        (element_array.begin() + i)->Initialize();
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibray TMMGLibray>
IndexVectorType MmgProcess<TMMGLibray>::CheckNodes()
{
    DoubleVectorMapType node_map;

    IndexVectorType nodes_to_remove_ids;

    DoubleVectorType coords(Dimension);

    NodesArrayType& nodes_array = mrThisModelPart.Nodes();

    for(SizeType i = 0; i < nodes_array.size(); ++i) {
        auto it_node = nodes_array.begin() + i;

        const array_1d<double, 3>& coordinates = it_node->Coordinates();

        for(IndexType i_coord = 0; i_coord < Dimension; i_coord++)
            coords[i_coord] = coordinates[i_coord];

        node_map[coords] += 1;

        if (node_map[coords] > 1) {
            nodes_to_remove_ids.push_back(it_node->Id());
            KRATOS_WARNING_IF("MmgProcess", mEchoLevel > 0) << "The mode " << it_node->Id() <<  " is repeated"<< std::endl;
        }
    }

    return nodes_to_remove_ids;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
IndexVectorType MmgProcess<MMGLibray::MMG2D>::CheckConditions0()
{
    IndexVectorMapType edge_map;

    IndexVectorType ids(2);

    IndexVectorType conditions_to_remove;

    // Iterate in the conditions
    for(int i = 0; i < mmgMesh->na; ++i) {
        int edge_0, edge_1, prop_id, is_ridge, is_required;

        if (MMG2D_Get_edge(mmgMesh, &edge_0, &edge_1, &prop_id, &is_ridge, &is_required) != 1 )
            exit(EXIT_FAILURE);

        ids[0] = edge_0;
        ids[1] = edge_1;

        //*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
        std::sort(ids.begin(), ids.end());

        edge_map[ids] += 1;

        if (edge_map[ids] > 1)
            conditions_to_remove.push_back(i + 1);
    }

    return conditions_to_remove;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
IndexVectorType MmgProcess<MMGLibray::MMG3D>::CheckConditions0()
{
    IndexVectorMapType triangle_map;

    IndexVectorType ids_triangles(3);

    IndexVectorType conditions_to_remove;

    for(int i = 0; i < mmgMesh->nt; ++i) {
        int vertex_0, vertex_1, vertex_2, prop_id, is_required;

        if (MMG3D_Get_triangle(mmgMesh, &vertex_0, &vertex_1, &vertex_2, &prop_id, &is_required) != 1 )
            exit(EXIT_FAILURE);

        ids_triangles[0] = vertex_0;
        ids_triangles[1] = vertex_1;
        ids_triangles[2] = vertex_2;

        //*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
        std::sort(ids_triangles.begin(), ids_triangles.end());

        triangle_map[ids_triangles] += 1;

        if (triangle_map[ids_triangles] > 1)
            conditions_to_remove.push_back(i + 1);
    }

    return conditions_to_remove;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
IndexVectorType MmgProcess<MMGLibray::MMGS>::CheckConditions0()
{
    IndexVectorMapType edge_map;

    IndexVectorType ids(2);

    IndexVectorType conditions_to_remove;

    // Iterate in the conditions
    for(int i = 0; i < mmgMesh->na; ++i) {
        int edge_0, edge_1, prop_id, is_ridge, is_required;

        if (MMGS_Get_edge(mmgMesh, &edge_0, &edge_1, &prop_id, &is_ridge, &is_required) != 1 )
            exit(EXIT_FAILURE);

        ids[0] = edge_0;
        ids[1] = edge_1;

        //*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
        std::sort(ids.begin(), ids.end());

        edge_map[ids] += 1;

        if (edge_map[ids] > 1)
            conditions_to_remove.push_back(i + 1);
    }

    return conditions_to_remove;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
IndexVectorType MmgProcess<MMGLibray::MMG2D>::CheckConditions1()
{
    IndexVectorType conditions_to_remove(0);

    return conditions_to_remove;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
IndexVectorType MmgProcess<MMGLibray::MMG3D>::CheckConditions1()
{
    IndexVectorMapType quadrilateral_map;

    IndexVectorType ids_quadrialteral(4);

    IndexVectorType conditions_to_remove;

    for(int i = 0; i < mmgMesh->nquad; ++i) {
        int vertex_0, vertex_1, vertex_2, vertex_3, prop_id, is_required;

        if (MMG3D_Get_quadrilateral(mmgMesh, &vertex_0, &vertex_1, &vertex_2, &vertex_3, &prop_id, &is_required) != 1 )
            exit(EXIT_FAILURE);

        ids_quadrialteral[0] = vertex_0;
        ids_quadrialteral[1] = vertex_1;
        ids_quadrialteral[2] = vertex_2;
        ids_quadrialteral[3] = vertex_3;

        //*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
        std::sort(ids_quadrialteral.begin(), ids_quadrialteral.end());

        quadrilateral_map[ids_quadrialteral] += 1;

        if (quadrilateral_map[ids_quadrialteral] > 1)
            conditions_to_remove.push_back(i + 1);
    }

    return conditions_to_remove;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
IndexVectorType MmgProcess<MMGLibray::MMGS>::CheckConditions1()
{
    IndexVectorType conditions_to_remove(0);

    return conditions_to_remove;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
IndexVectorType MmgProcess<MMGLibray::MMG2D>::CheckElements0()
{
    IndexVectorMapType triangle_map;

    IndexVectorType ids_triangles(3);

    IndexVectorType elements_to_remove;

    // Iterate in the elements
    for(int i = 0; i < mmgMesh->nt; ++i) {
        int vertex_0, vertex_1, vertex_2, prop_id, is_required;

        if (MMG2D_Get_triangle(mmgMesh, &vertex_0, &vertex_1, &vertex_2, &prop_id, &is_required) != 1 )
            exit(EXIT_FAILURE);

        ids_triangles[0] = vertex_0;
        ids_triangles[1] = vertex_1;
        ids_triangles[2] = vertex_2;

        //*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
        std::sort(ids_triangles.begin(), ids_triangles.end());

        triangle_map[ids_triangles] += 1;

        if (triangle_map[ids_triangles] > 1)
            elements_to_remove.push_back(i + 1);
    }

    return elements_to_remove;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
IndexVectorType MmgProcess<MMGLibray::MMG3D>::CheckElements0()
{
    IndexVectorMapType triangle_map;

    IndexVectorType ids_tetrahedron(4);

    IndexVectorType elements_to_remove;

    for(int i = 0; i < mmgMesh->ne; ++i) {
        int vertex_0, vertex_1, vertex_2, vertex_3, prop_id, is_required;

        if (MMG3D_Get_tetrahedron(mmgMesh, &vertex_0, &vertex_1, &vertex_2, &vertex_3, &prop_id, &is_required) != 1 )
            exit(EXIT_FAILURE);

        ids_tetrahedron[0] = vertex_0;
        ids_tetrahedron[1] = vertex_1;
        ids_tetrahedron[2] = vertex_2;
        ids_tetrahedron[3] = vertex_3;

        //*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
        std::sort(ids_tetrahedron.begin(), ids_tetrahedron.end());

        triangle_map[ids_tetrahedron] += 1;

        if (triangle_map[ids_tetrahedron] > 1)
            elements_to_remove.push_back(i + 1);
    }

    return elements_to_remove;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
IndexVectorType MmgProcess<MMGLibray::MMGS>::CheckElements0()
{
    IndexVectorMapType triangle_map;

    IndexVectorType ids_triangles(3);

    IndexVectorType elements_to_remove;

    // Iterate in the elements
    for(int i = 0; i < mmgMesh->nt; ++i) {
        int vertex_0, vertex_1, vertex_2, prop_id, is_required;

        if (MMGS_Get_triangle(mmgMesh, &vertex_0, &vertex_1, &vertex_2, &prop_id, &is_required) != 1 )
            exit(EXIT_FAILURE);

        ids_triangles[0] = vertex_0;
        ids_triangles[1] = vertex_1;
        ids_triangles[2] = vertex_2;

        //*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
        std::sort(ids_triangles.begin(), ids_triangles.end());

        triangle_map[ids_triangles] += 1;

        if (triangle_map[ids_triangles] > 1)
            elements_to_remove.push_back(i + 1);
    }

    return elements_to_remove;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
IndexVectorType MmgProcess<MMGLibray::MMG2D>::CheckElements1()
{
    IndexVectorType elements_to_remove(0);
    return elements_to_remove;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
IndexVectorType MmgProcess<MMGLibray::MMG3D>::CheckElements1()
{
    IndexVectorMapType prism_map;

    IndexVectorType ids_prisms(6);

    IndexVectorType elements_to_remove;

    for(int i = 0; i < mmgMesh->nprism; ++i) {
        int vertex_0, vertex_1, vertex_2, vertex_3, vertex_4, vertex_5, prop_id, is_required;

        if (MMG3D_Get_prism(mmgMesh, &vertex_0, &vertex_1, &vertex_2, &vertex_3, &vertex_4, &vertex_5, &prop_id, &is_required) != 1 )
            exit(EXIT_FAILURE);

        ids_prisms[0] = vertex_0;
        ids_prisms[1] = vertex_1;
        ids_prisms[2] = vertex_2;
        ids_prisms[3] = vertex_3;
        ids_prisms[4] = vertex_4;
        ids_prisms[5] = vertex_5;

        //*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
        std::sort(ids_prisms.begin(), ids_prisms.end());

        prism_map[ids_prisms] += 1;

        if (prism_map[ids_prisms] > 1)
            elements_to_remove.push_back(i + 1);
    }

    return elements_to_remove;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
IndexVectorType MmgProcess<MMGLibray::MMGS>::CheckElements1()
{
    IndexVectorType elements_to_remove(0);
    return elements_to_remove;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMG2D>::BlockNode(IndexType iNode)
{
    if (MMG2D_Set_requiredVertex(mmgMesh, iNode) != 1 )
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/


template<>
void MmgProcess<MMGLibray::MMG3D>::BlockNode(IndexType iNode)
{
    if (MMG3D_Set_requiredVertex(mmgMesh, iNode) != 1 )
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMGS>::BlockNode(IndexType iNode)
{
    if (MMGS_Set_requiredVertex(mmgMesh, iNode) != 1 )
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
NodeType::Pointer MmgProcess<MMGLibray::MMG2D>::CreateNode(
    IndexType iNode,
    int& Ref,
    int& IsRequired
    )
{
    double coord_0, coord_1;
    int is_corner;

    if (MMG2D_Get_vertex(mmgMesh, &coord_0, &coord_1, &Ref, &is_corner, &IsRequired) != 1 )
        exit(EXIT_FAILURE);

    NodeType::Pointer p_node = mrThisModelPart.CreateNewNode(iNode, coord_0, coord_1, 0.0);

    return p_node;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
NodeType::Pointer MmgProcess<MMGLibray::MMG3D>::CreateNode(
    IndexType iNode,
    int& Ref,
    int& IsRequired
    )
{
    double coord_0, coord_1, coord_2;
    int is_corner;

    if (MMG3D_Get_vertex(mmgMesh, &coord_0, &coord_1, &coord_2, &Ref, &is_corner, &IsRequired) != 1 )
        exit(EXIT_FAILURE);

    NodeType::Pointer p_node = mrThisModelPart.CreateNewNode(iNode, coord_0, coord_1, coord_2);

    return p_node;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
NodeType::Pointer MmgProcess<MMGLibray::MMGS>::CreateNode(
    IndexType iNode,
    int& Ref,
    int& IsRequired
    )
{
    double coord_0, coord_1, coord_2;
    int is_corner;

    if (MMGS_Get_vertex(mmgMesh, &coord_0, &coord_1, &coord_2, &Ref, &is_corner, &IsRequired) != 1 )
        exit(EXIT_FAILURE);

    NodeType::Pointer p_node = mrThisModelPart.CreateNewNode(iNode, coord_0, coord_1, coord_2);

    return p_node;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
ConditionType::Pointer MmgProcess<MMGLibray::MMG2D>::CreateCondition0(
    const IndexType CondId,
    int& PropId,
    int& IsRequired,
    bool SkipCreation
    )
{
    // Sometimes MMG creates conditions where there are not, then we skip
    if (mpRefCondition[PropId] == nullptr) return nullptr;

    // We create the default one
    ConditionType::Pointer p_condition = nullptr;

    int edge_0, edge_1, is_ridge;

    if (MMG2D_Get_edge(mmgMesh, &edge_0, &edge_1, &PropId, &is_ridge, &IsRequired) != 1 )
        exit(EXIT_FAILURE);

    // FIXME: This is not the correct solution to the problem, I asked in the MMG Forum
    if (edge_0 == 0) SkipCreation = true;
    if (edge_1 == 0) SkipCreation = true;

    if (SkipCreation == false) {
        std::vector<NodeType::Pointer> condition_nodes (2);
        condition_nodes[0] = mrThisModelPart.pGetNode(edge_0);
        condition_nodes[1] = mrThisModelPart.pGetNode(edge_1);

        p_condition = mpRefCondition[PropId]->Create(CondId, PointerVector<NodeType>{condition_nodes}, mpRefCondition[PropId]->pGetProperties());
    } else if (mEchoLevel > 2)
        KRATOS_INFO("MmgProcess") << "Condition creation avoided" << std::endl;

    return p_condition;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
ConditionType::Pointer MmgProcess<MMGLibray::MMG3D>::CreateCondition0(
    const IndexType CondId,
    int& PropId,
    int& IsRequired,
    bool SkipCreation
    )
{
    // Sometimes MMG creates conditions where there are not, then we skip
    if (mpRefCondition[PropId] == nullptr) return nullptr;

    // We create the default one
    ConditionType::Pointer p_condition = nullptr;

    int vertex_0, vertex_1, vertex_2;

    if (MMG3D_Get_triangle(mmgMesh, &vertex_0, &vertex_1, &vertex_2, &PropId, &IsRequired) != 1 )
        exit(EXIT_FAILURE);

    // FIXME: This is not the correct solution to the problem, I asked in the MMG Forum
    if (vertex_0 == 0) SkipCreation = true;
    if (vertex_1 == 0) SkipCreation = true;
    if (vertex_2 == 0) SkipCreation = true;

    if (SkipCreation == false) {
        std::vector<NodeType::Pointer> condition_nodes (3);
        condition_nodes[0] = mrThisModelPart.pGetNode(vertex_0);
        condition_nodes[1] = mrThisModelPart.pGetNode(vertex_1);
        condition_nodes[2] = mrThisModelPart.pGetNode(vertex_2);

        p_condition = mpRefCondition[PropId]->Create(CondId, PointerVector<NodeType>{condition_nodes}, mpRefCondition[PropId]->pGetProperties());
    } else if (mEchoLevel > 2)
        KRATOS_WARNING("MmgProcess") << "Condition creation avoided" << std::endl;

    return p_condition;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
ConditionType::Pointer MmgProcess<MMGLibray::MMGS>::CreateCondition0(
    const IndexType CondId,
    int& PropId,
    int& IsRequired,
    bool SkipCreation
    )
{
    // Sometimes MMG creates conditions where there are not, then we skip
    if (mpRefCondition[PropId] == nullptr) return nullptr;

    // We create the default one
    ConditionType::Pointer p_condition = nullptr;

    int edge_0, edge_1, is_ridge;

    if (MMGS_Get_edge(mmgMesh, &edge_0, &edge_1, &PropId, &is_ridge, &IsRequired) != 1 )
        exit(EXIT_FAILURE);

    // FIXME: This is not the correct solution to the problem, I asked in the MMG Forum
    if (edge_0 == 0) SkipCreation = true;
    if (edge_1 == 0) SkipCreation = true;

    if (SkipCreation == false) {
        std::vector<NodeType::Pointer> condition_nodes (2);
        condition_nodes[0] = mrThisModelPart.pGetNode(edge_0);
        condition_nodes[1] = mrThisModelPart.pGetNode(edge_1);

        p_condition = mpRefCondition[PropId]->Create(CondId, PointerVector<NodeType>{condition_nodes}, mpRefCondition[PropId]->pGetProperties());
    } else if (mEchoLevel > 2)
        KRATOS_INFO("MmgProcess") << "Condition creation avoided" << std::endl;

    return p_condition;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
ConditionType::Pointer MmgProcess<MMGLibray::MMG2D>::CreateCondition1(
    const IndexType CondId,
    int& PropId,
    int& IsRequired,
    bool SkipCreation
    )
{
    return nullptr;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
ConditionType::Pointer MmgProcess<MMGLibray::MMG3D>::CreateCondition1(
    const IndexType CondId,
    int& PropId,
    int& IsRequired,
    bool SkipCreation
    )
{
    ConditionType::Pointer p_condition = nullptr;

    int vertex_0, vertex_1, vertex_2, vertex_3;

    if (MMG3D_Get_quadrilateral(mmgMesh, &vertex_0, &vertex_1, &vertex_2, &vertex_3, &PropId, &IsRequired) != 1 )
        exit(EXIT_FAILURE);

    // FIXME: This is not the correct solution to the problem, I asked in the MMG Forum
    if (vertex_0 == 0) SkipCreation = true;
    if (vertex_1 == 0) SkipCreation = true;
    if (vertex_2 == 0) SkipCreation = true;
    if (vertex_3 == 0) SkipCreation = true;

    if (SkipCreation == false) {
        std::vector<NodeType::Pointer> condition_nodes (4);
        condition_nodes[0] = mrThisModelPart.pGetNode(vertex_0);
        condition_nodes[1] = mrThisModelPart.pGetNode(vertex_1);
        condition_nodes[2] = mrThisModelPart.pGetNode(vertex_2);
        condition_nodes[3] = mrThisModelPart.pGetNode(vertex_3);

        p_condition = mpRefCondition[PropId]->Create(CondId, PointerVector<NodeType>{condition_nodes}, mpRefCondition[PropId]->pGetProperties());
    } else if (mEchoLevel > 2)
        KRATOS_WARNING("MmgProcess") << "Condition creation avoided" << std::endl;

    return p_condition;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
ConditionType::Pointer MmgProcess<MMGLibray::MMGS>::CreateCondition1(
    const IndexType CondId,
    int& PropId,
    int& IsRequired,
    bool SkipCreation
    )
{
    return nullptr;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
ElementType::Pointer MmgProcess<MMGLibray::MMG2D>::CreateElement0(
    const IndexType ElemId,
    int& PropId,
    int& IsRequired,
    bool SkipCreation
    )
{
    ElementType::Pointer p_element = nullptr;

    int vertex_0, vertex_1, vertex_2;

    if (MMG2D_Get_triangle(mmgMesh, &vertex_0, &vertex_1, &vertex_2, &PropId, &IsRequired) != 1 )
        exit(EXIT_FAILURE);

    // FIXME: This is not the correct solution to the problem, I asked in the MMG Forum
    if (vertex_0 == 0) SkipCreation = true;
    if (vertex_1 == 0) SkipCreation = true;
    if (vertex_2 == 0) SkipCreation = true;

    if (SkipCreation == false) {
        std::vector<NodeType::Pointer> element_nodes (3);
        element_nodes[0] = mrThisModelPart.pGetNode(vertex_0);
        element_nodes[1] = mrThisModelPart.pGetNode(vertex_1);
        element_nodes[2] = mrThisModelPart.pGetNode(vertex_2);

        p_element = mpRefElement[PropId]->Create(ElemId, PointerVector<NodeType>{element_nodes}, mpRefElement[PropId]->pGetProperties());
    } else if (mEchoLevel > 2)
        KRATOS_WARNING("MmgProcess") << "Element creation avoided" << std::endl;

    return p_element;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
ElementType::Pointer MmgProcess<MMGLibray::MMG3D>::CreateElement0(
    const IndexType ElemId,
    int& PropId,
    int& IsRequired,
    bool SkipCreation
    )
{
    ElementType::Pointer p_element = nullptr;

    int vertex_0, vertex_1, vertex_2, vertex_3;

    if (MMG3D_Get_tetrahedron(mmgMesh, &vertex_0, &vertex_1, &vertex_2, &vertex_3, &PropId, &IsRequired) != 1 )
        exit(EXIT_FAILURE);

    // FIXME: This is not the correct solution to the problem, I asked in the MMG Forum
    if (vertex_0 == 0) SkipCreation = true;
    if (vertex_1 == 0) SkipCreation = true;
    if (vertex_2 == 0) SkipCreation = true;
    if (vertex_3 == 0) SkipCreation = true;

    if (SkipCreation == false) {
        std::vector<NodeType::Pointer> element_nodes (4);
        element_nodes[0] = mrThisModelPart.pGetNode(vertex_0);
        element_nodes[1] = mrThisModelPart.pGetNode(vertex_1);
        element_nodes[2] = mrThisModelPart.pGetNode(vertex_2);
        element_nodes[3] = mrThisModelPart.pGetNode(vertex_3);

        p_element = mpRefElement[PropId]->Create(ElemId, PointerVector<NodeType>{element_nodes}, mpRefElement[PropId]->pGetProperties());
    } else if (mEchoLevel > 2)
        KRATOS_WARNING("MmgProcess") << "Element creation avoided" << std::endl;

    return p_element;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
ElementType::Pointer MmgProcess<MMGLibray::MMGS>::CreateElement0(
    const IndexType ElemId,
    int& PropId,
    int& IsRequired,
    bool SkipCreation
    )
{
    ElementType::Pointer p_element = nullptr;

    int vertex_0, vertex_1, vertex_2;

    if (MMGS_Get_triangle(mmgMesh, &vertex_0, &vertex_1, &vertex_2, &PropId, &IsRequired) != 1 )
        exit(EXIT_FAILURE);

    // FIXME: This is not the correct solution to the problem, I asked in the MMG Forum
    if (vertex_0 == 0) SkipCreation = true;
    if (vertex_1 == 0) SkipCreation = true;
    if (vertex_2 == 0) SkipCreation = true;

    if (SkipCreation == false) {
        std::vector<NodeType::Pointer> element_nodes (3);
        element_nodes[0] = mrThisModelPart.pGetNode(vertex_0);
        element_nodes[1] = mrThisModelPart.pGetNode(vertex_1);
        element_nodes[2] = mrThisModelPart.pGetNode(vertex_2);

        p_element = mpRefElement[PropId]->Create(ElemId, PointerVector<NodeType>{element_nodes}, mpRefElement[PropId]->pGetProperties());
    } else if (mEchoLevel > 2)
        KRATOS_WARNING("MmgProcess") << "Element creation avoided" << std::endl;

    return p_element;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
ElementType::Pointer MmgProcess<MMGLibray::MMG2D>::CreateElement1(
    const IndexType ElemId,
    int& PropId,
    int& IsRequired,
    bool SkipCreation
    )
{
    return nullptr;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
ElementType::Pointer MmgProcess<MMGLibray::MMG3D>::CreateElement1(
    const IndexType ElemId,
    int& PropId,
    int& IsRequired,
    bool SkipCreation
    )
{
    ElementType::Pointer p_element = nullptr;

    int vertex_0, vertex_1, vertex_2, vertex_3, vertex_4, vertex_5;

    if (MMG3D_Get_prism(mmgMesh, &vertex_0, &vertex_1, &vertex_2, &vertex_3, &vertex_4, &vertex_5, &PropId, &IsRequired) != 1 )
        exit(EXIT_FAILURE);

    // FIXME: This is not the correct solution to the problem, I asked in the MMG Forum
    if (vertex_0 == 0) SkipCreation = true;
    if (vertex_1 == 0) SkipCreation = true;
    if (vertex_2 == 0) SkipCreation = true;
    if (vertex_3 == 0) SkipCreation = true;
    if (vertex_4 == 0) SkipCreation = true;
    if (vertex_5 == 0) SkipCreation = true;

    if (SkipCreation == false) {
        std::vector<NodeType::Pointer> element_nodes (6);
        element_nodes[0] = mrThisModelPart.pGetNode(vertex_0);
        element_nodes[1] = mrThisModelPart.pGetNode(vertex_1);
        element_nodes[2] = mrThisModelPart.pGetNode(vertex_2);
        element_nodes[3] = mrThisModelPart.pGetNode(vertex_3);
        element_nodes[4] = mrThisModelPart.pGetNode(vertex_4);
        element_nodes[5] = mrThisModelPart.pGetNode(vertex_5);

        p_element = mpRefElement[PropId]->Create(ElemId, PointerVector<NodeType>{element_nodes}, mpRefElement[PropId]->pGetProperties());
    } else if (mEchoLevel > 2)
        KRATOS_WARNING("MmgProcess") << "Element creation avoided" << std::endl;

    return p_element;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
ElementType::Pointer MmgProcess<MMGLibray::MMGS>::CreateElement1(
    const IndexType ElemId,
    int& PropId,
    int& IsRequired,
    bool SkipCreation
    )
{
    return nullptr;
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibray TMMGLibray>
void MmgProcess<TMMGLibray>::SaveSolutionToFile(const bool PostOutput)
{
    /* GET RESULTS */

    const IndexType step = mrThisModelPart.GetProcessInfo()[STEP];

    // Automatically save the mesh
    OutputMesh(PostOutput, step);

    // Automatically save the solution
    OutputSol(PostOutput, step);

    // Save the mesh in an .mdpa format
    const bool save_mdpa_file = mThisParameters["save_mdpa_file"].GetBool();
    if(save_mdpa_file) OutputMdpa();
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibray TMMGLibray>
void MmgProcess<TMMGLibray>::FreeMemory()
{
    // Free the MMG structures
    FreeAll();

    // Free filename (NOTE: Problems with more that one iteration)
//     free(mFilename);
//     mFilename = nullptr;

    // Free reference std::unordered_map
    mpRefElement.clear();
    mpRefCondition.clear();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMG2D>::InitMesh()
{
    mmgMesh = nullptr;
    mmgSol = nullptr;

    // We init the MMG mesh and sol
    MMG2D_Init_mesh( MMG5_ARG_start, MMG5_ARG_ppMesh, &mmgMesh, MMG5_ARG_ppMet, &mmgSol, MMG5_ARG_end);

    InitVerbosity();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMG3D>::InitMesh()
{
    mmgMesh = nullptr;
    mmgSol = nullptr;

    // We init the MMG mesh and sol
    MMG3D_Init_mesh( MMG5_ARG_start, MMG5_ARG_ppMesh, &mmgMesh, MMG5_ARG_ppMet, &mmgSol, MMG5_ARG_end);

    InitVerbosity();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMGS>::InitMesh()
{
    mmgMesh = nullptr;
    mmgSol = nullptr;

    // We init the MMG mesh and sol
    MMGS_Init_mesh( MMG5_ARG_start, MMG5_ARG_ppMesh, &mmgMesh, MMG5_ARG_ppMet, &mmgSol, MMG5_ARG_end);

    InitVerbosity();
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibray TMMGLibray>
void MmgProcess<TMMGLibray>::InitVerbosity()
{
    /* We set the MMG verbosity */
    int verbosity_mmg;
    if (mEchoLevel == 0)
        verbosity_mmg = -1;
    else if (mEchoLevel == 1)
        verbosity_mmg = 0; // NOTE: This way just the essential info from MMG will be printed, but the custom message will appear
    else if (mEchoLevel == 2)
        verbosity_mmg = 3;
    else if (mEchoLevel == 3)
        verbosity_mmg = 5;
    else
        verbosity_mmg = 10;

    InitVerbosityParameter(verbosity_mmg);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMG2D>::InitVerbosityParameter(const IndexType VerbosityMMG)
{
    if ( !MMG2D_Set_iparameter(mmgMesh,mmgSol,MMG2D_IPARAM_verbose, VerbosityMMG) )
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMG3D>::InitVerbosityParameter(const IndexType VerbosityMMG)
{
    if ( !MMG3D_Set_iparameter(mmgMesh,mmgSol,MMG3D_IPARAM_verbose, VerbosityMMG) )
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMGS>::InitVerbosityParameter(const IndexType VerbosityMMG)
{
    if ( !MMGS_Set_iparameter(mmgMesh,mmgSol,MMGS_IPARAM_verbose, VerbosityMMG) )
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMG2D>::SetMeshSize(
    const SizeType NumNodes,
    const array_1d<SizeType, ElementsArraySize>& NumArrayElements,
    const array_1d<SizeType, ConditionsArraySize>& NumArrayConditions
    )
{
    //Give the size of the mesh: NumNodes vertices, num_elements triangles, num_conditions edges (2D)
    if ( MMG2D_Set_meshSize(mmgMesh, NumNodes, NumArrayElements[0], NumArrayConditions[0]) != 1 )
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMG3D>::SetMeshSize(
    const SizeType NumNodes,
    const array_1d<SizeType, ElementsArraySize>& NumArrayElements,  // NOTE: We do this tricky thing to take into account the prisms
    const array_1d<SizeType, ConditionsArraySize>& NumArrayConditions // NOTE: We do this tricky thing to take into account the quadrilaterals
    )
{
    //Give the size of the mesh: NumNodes Vertex, num_elements tetra and prism, NumArrayConditions triangles and quadrilaterals, 0 edges (3D)
    if ( MMG3D_Set_meshSize(mmgMesh, NumNodes, NumArrayElements[0], NumArrayElements[1], NumArrayConditions[0], NumArrayConditions[1], 0) != 1 )
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMGS>::SetMeshSize(
    const SizeType NumNodes,
    const array_1d<SizeType, ElementsArraySize>& NumArrayElements,
    const array_1d<SizeType, ConditionsArraySize>& NumArrayConditions
    )
{
    //Give the size of the mesh: NumNodes vertices, num_elements triangles, num_conditions edges (2D)
    if ( MMGS_Set_meshSize(mmgMesh, NumNodes, NumArrayElements[0], NumArrayConditions[0]) != 1 )
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMG2D>::SetSolSizeScalar(const SizeType NumNodes)
{
    if ( MMG2D_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,NumNodes,MMG5_Scalar) != 1 )
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMG3D>::SetSolSizeScalar(const SizeType NumNodes)
{
    if ( MMG3D_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,NumNodes,MMG5_Scalar) != 1 )
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMGS>::SetSolSizeScalar(const SizeType NumNodes)
{
    if ( MMGS_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,NumNodes,MMG5_Scalar) != 1 )
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMG2D>::SetSolSizeVector(const SizeType NumNodes)
{
    if ( MMG2D_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,NumNodes,MMG5_Vector) != 1 )
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMG3D>::SetSolSizeVector(const SizeType NumNodes)
{
    if ( MMG3D_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,NumNodes,MMG5_Vector) != 1 )
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMGS>::SetSolSizeVector(const SizeType NumNodes)
{
    if ( MMGS_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,NumNodes,MMG5_Vector) != 1 )
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMG2D>::SetSolSizeTensor(const SizeType NumNodes)
{
    if ( MMG2D_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,NumNodes,MMG5_Tensor) != 1 )
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMG3D>::SetSolSizeTensor(const SizeType NumNodes)
{
    if ( MMG3D_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,NumNodes,MMG5_Tensor) != 1 )
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMGS>::SetSolSizeTensor(const SizeType NumNodes)
{
    if ( MMGS_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,NumNodes,MMG5_Tensor) != 1 )
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMG2D>::CheckMeshData()
{
    if ( MMG2D_Chk_meshData(mmgMesh, mmgSol) != 1 )
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMG3D>::CheckMeshData()
{
    if ( MMG3D_Chk_meshData(mmgMesh, mmgSol) != 1 )
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMGS>::CheckMeshData()
{
    if ( MMGS_Chk_meshData(mmgMesh, mmgSol) != 1 )
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMG2D>::OutputMesh(
    const bool PostOutput,
    const IndexType Step
    )
{
    std::string mesh_name;
    if (PostOutput)
        mesh_name = mStdStringFilename+"_step="+std::to_string(Step)+".o.mesh";
    else
        mesh_name = mStdStringFilename+"_step="+std::to_string(Step)+".mesh";

    auto  mesh_file = new char [mesh_name.length() + 1];
    std::strcpy (mesh_file, mesh_name.c_str());

    // a)  Give the ouptut mesh name using MMG2D_Set_outputMeshName (by default, the mesh is saved in the "mesh.o.mesh" file
    MMG2D_Set_outputMeshName(mmgMesh,mesh_file);

    // b) function calling
    KRATOS_INFO_IF("MmgProcess", MMG2D_saveMesh(mmgMesh,mesh_file) != 1) << "UNABLE TO SAVE MESH" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMG3D>::OutputMesh(
    const bool PostOutput,
    const IndexType Step
    )
{
    std::string mesh_name;
    if (PostOutput)
        mesh_name = mStdStringFilename+"_step="+std::to_string(Step)+".o.mesh";
    else
        mesh_name = mStdStringFilename+"_step="+std::to_string(Step)+".mesh";

    auto  mesh_file = new char [mesh_name.length() + 1];
    std::strcpy (mesh_file, mesh_name.c_str());

    // a)  Give the ouptut mesh name using MMG3D_Set_outputMeshName (by default, the mesh is saved in the "mesh.o.mesh" file
    MMG3D_Set_outputMeshName(mmgMesh,mesh_file);

    // b) function calling
    KRATOS_INFO_IF("MmgProcess", MMG3D_saveMesh(mmgMesh,mesh_file) != 1) << "UNABLE TO SAVE MESH" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMGS>::OutputMesh(
    const bool PostOutput,
    const IndexType Step
    )
{
    std::string mesh_name;
    if (PostOutput)
        mesh_name = mStdStringFilename+"_step="+std::to_string(Step)+".o.mesh";
    else
        mesh_name = mStdStringFilename+"_step="+std::to_string(Step)+".mesh";

    auto  mesh_file = new char [mesh_name.length() + 1];
    std::strcpy (mesh_file, mesh_name.c_str());

    // a)  Give the ouptut mesh name using MMGS_Set_outputMeshName (by default, the mesh is saved in the "mesh.o.mesh" file
    MMGS_Set_outputMeshName(mmgMesh,mesh_file);

    // b) function calling
    KRATOS_INFO_IF("MmgProcess", MMGS_saveMesh(mmgMesh,mesh_file) != 1) << "UNABLE TO SAVE MESH" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibray TMMGLibray>
void MmgProcess<TMMGLibray>::OutputMdpa()
{
    std::ofstream output_file;
    ModelPartIO model_part_io("output", IO::WRITE);
    model_part_io.WriteModelPart(mrThisModelPart);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMG2D>::OutputSol(
    const bool PostOutput,
    const IndexType Step
    )
{
    std::string sol_name;
    if (PostOutput)
        sol_name = mStdStringFilename+"_step="+std::to_string(Step)+".o.sol";
    else
        sol_name = mStdStringFilename+"_step="+std::to_string(Step)+".sol";

    auto  sol_file = new char [sol_name.length() + 1];
    std::strcpy (sol_file, sol_name.c_str());

    // a)  Give the ouptut sol name using MMG2D_Set_outputSolName (by default, the mesh is saved in the "mesh.o.sol" file
    MMG2D_Set_outputSolName(mmgMesh, mmgSol, sol_file);

    // b) Function calling
    KRATOS_INFO_IF("MmgProcess", MMG2D_saveSol(mmgMesh, mmgSol, sol_file) != 1) << "UNABLE TO SAVE SOL" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMG3D>::OutputSol(
    const bool PostOutput,
    const IndexType Step
    )
{
    std::string sol_name;
    if (PostOutput)
        sol_name = mStdStringFilename+"_step="+std::to_string(Step)+".o.sol";
    else
        sol_name = mStdStringFilename+"_step="+std::to_string(Step)+".sol";

    auto  sol_file = new char [sol_name.length() + 1];
    std::strcpy (sol_file, sol_name.c_str());

    // a)  Give the ouptut sol name using MMG3D_Set_outputSolName (by default, the mesh is saved in the "mesh.o.sol" file
    MMG3D_Set_outputSolName(mmgMesh, mmgSol, sol_file);

    // b) Function calling
    KRATOS_INFO_IF("MmgProcess", MMG3D_saveSol(mmgMesh,mmgSol, sol_file) != 1)<< "UNABLE TO SAVE SOL" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMGS>::OutputSol(
    const bool PostOutput,
    const IndexType Step
    )
{
    std::string sol_name;
    if (PostOutput)
        sol_name = mStdStringFilename+"_step="+std::to_string(Step)+".o.sol";
    else
        sol_name = mStdStringFilename+"_step="+std::to_string(Step)+".sol";

    auto  sol_file = new char [sol_name.length() + 1];
    std::strcpy (sol_file, sol_name.c_str());

    // a)  Give the ouptut sol name using MMGS_Set_outputSolName (by default, the mesh is saved in the "mesh.o.sol" file
    MMGS_Set_outputSolName(mmgMesh, mmgSol, sol_file);

    // b) Function calling
    KRATOS_INFO_IF("MmgProcess", MMGS_saveSol(mmgMesh,mmgSol, sol_file) != 1)<< "UNABLE TO SAVE SOL" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMG2D>::MMGLibCall()
{
    KRATOS_TRY;

    /* Advanced configurations */
    // Global hausdorff value (default value = 0.01) applied on the whole boundary
    if (mThisParameters["advanced_parameters"]["force_hausdorff_value"].GetBool()) {
        if ( MMG2D_Set_dparameter(mmgMesh,mmgSol,MMG2D_DPARAM_hausd, mThisParameters["advanced_parameters"]["hausdorff_value"].GetDouble()) != 1 )
            KRATOS_ERROR << "Unable to set the Hausdorff parameter" << std::endl;
    }

    // Avoid/allow point relocation
    if ( MMG2D_Set_iparameter(mmgMesh,mmgSol,MMG2D_IPARAM_nomove, static_cast<int>(mThisParameters["advanced_parameters"]["no_move_mesh"].GetBool())) != 1 )
        KRATOS_ERROR << "Unable to fix the nodes" << std::endl;

    // Avoid/allow surface modifications
    if ( MMG2D_Set_iparameter(mmgMesh,mmgSol,MMG2D_IPARAM_nosurf, static_cast<int>(mThisParameters["advanced_parameters"]["no_surf_mesh"].GetBool())) != 1 )
        KRATOS_ERROR << "Unable to set no surfacic modifications" << std::endl;

    // Don't insert nodes on mesh
    if ( MMG2D_Set_iparameter(mmgMesh,mmgSol,MMG2D_IPARAM_noinsert, static_cast<int>(mThisParameters["advanced_parameters"]["no_insert_mesh"].GetBool())) != 1 )
        KRATOS_ERROR << "Unable to set no insertion/suppression point" << std::endl;

    // Don't swap mesh
    if ( MMG2D_Set_iparameter(mmgMesh,mmgSol,MMG2D_IPARAM_noswap, static_cast<int>(mThisParameters["advanced_parameters"]["no_swap_mesh"].GetBool())) != 1 )
        KRATOS_ERROR << "Unable to set no edge flipping" << std::endl;

    // Set the angle detection
    const bool deactivate_detect_angle = mThisParameters["advanced_parameters"]["deactivate_detect_angle"].GetBool();
    if ( deactivate_detect_angle) {
        if ( MMG2D_Set_iparameter(mmgMesh,mmgSol,MMG2D_IPARAM_angle, static_cast<int>(!deactivate_detect_angle)) != 1 )
            KRATOS_ERROR << "Unable to set the angle detection on" << std::endl;
    }

    // Set the gradation
    if (mThisParameters["advanced_parameters"]["force_gradation_value"].GetBool()) {
        if ( MMG2D_Set_dparameter(mmgMesh,mmgSol,MMG2D_DPARAM_hgrad, mThisParameters["advanced_parameters"]["gradation_value"].GetDouble()) != 1 )
            KRATOS_ERROR << "Unable to set gradation" << std::endl;
    }

    // Minimal edge size
    if (mThisParameters["force_sizes"]["force_min"].GetBool()) {
        if ( MMG2D_Set_dparameter(mmgMesh,mmgSol,MMG2D_DPARAM_hmin, mThisParameters["force_sizes"]["minimal_size"].GetDouble()) != 1 )
            KRATOS_ERROR << "Unable to set the minimal edge size " << std::endl;
    }

    // Minimal edge size
    if (mThisParameters["force_sizes"]["force_max"].GetBool()) {
        if ( MMG2D_Set_dparameter(mmgMesh,mmgSol,MMG2D_DPARAM_hmax, mThisParameters["force_sizes"]["maximal_size"].GetDouble()) != 1 ) {
            KRATOS_ERROR << "Unable to set the maximal edge size " << std::endl;
        }
    }

    const int ier = MMG2D_mmg2dlib(mmgMesh, mmgSol);

    if ( ier == MMG5_STRONGFAILURE )
        KRATOS_ERROR << "ERROR: BAD ENDING OF MMG2DLIB: UNABLE TO SAVE MESH. ier: " << ier << std::endl;
    else if ( ier == MMG5_LOWFAILURE )
        KRATOS_ERROR << "ERROR: BAD ENDING OF MMG2DLIB. ier: " << ier << std::endl;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMG3D>::MMGLibCall()
{
    KRATOS_TRY;

    /* Advanced configurations */
    // Global hausdorff value (default value = 0.01) applied on the whole boundary
    if (mThisParameters["advanced_parameters"]["force_hausdorff_value"].GetBool()) {
        if ( MMG3D_Set_dparameter(mmgMesh,mmgSol,MMG3D_DPARAM_hausd, mThisParameters["advanced_parameters"]["hausdorff_value"].GetDouble()) != 1 )
            KRATOS_ERROR << "Unable to set the Hausdorff parameter" << std::endl;
    }

    // Avoid/allow point relocation
    if ( MMG3D_Set_iparameter(mmgMesh,mmgSol,MMG3D_IPARAM_nomove, static_cast<int>(mThisParameters["advanced_parameters"]["no_move_mesh"].GetBool())) != 1 )
        KRATOS_ERROR << "Unable to fix the nodes" << std::endl;

    // Avoid/allow surface modifications
    if ( MMG3D_Set_iparameter(mmgMesh,mmgSol,MMG3D_IPARAM_nosurf, static_cast<int>(mThisParameters["advanced_parameters"]["no_surf_mesh"].GetBool())) != 1 )
        KRATOS_ERROR << "Unable to set no surfacic modifications" << std::endl;

    // Don't insert nodes on mesh
    if ( MMG3D_Set_iparameter(mmgMesh,mmgSol,MMG3D_IPARAM_noinsert, static_cast<int>(mThisParameters["advanced_parameters"]["no_insert_mesh"].GetBool())) != 1 )
        KRATOS_ERROR << "Unable to set no insertion/suppression point" << std::endl;

    // Don't swap mesh
    if ( MMG3D_Set_iparameter(mmgMesh,mmgSol,MMG3D_IPARAM_noswap, static_cast<int>(mThisParameters["advanced_parameters"]["no_swap_mesh"].GetBool())) != 1 )
        KRATOS_ERROR << "Unable to set no edge flipping" << std::endl;

    // Set the angle detection
    const bool deactivate_detect_angle = mThisParameters["advanced_parameters"]["deactivate_detect_angle"].GetBool();
    if ( deactivate_detect_angle) {
        if ( MMG3D_Set_iparameter(mmgMesh,mmgSol,MMG3D_IPARAM_angle, static_cast<int>(!deactivate_detect_angle)) != 1 )
            KRATOS_ERROR << "Unable to set the angle detection on" << std::endl;
    }

    // Set the gradation
    if (mThisParameters["advanced_parameters"]["force_gradation_value"].GetBool()) {
        if ( MMG3D_Set_dparameter(mmgMesh,mmgSol,MMG3D_DPARAM_hgrad, mThisParameters["advanced_parameters"]["gradation_value"].GetDouble()) != 1 )
            KRATOS_ERROR << "Unable to set gradation" << std::endl;
    }

    // Minimal edge size
    if (mThisParameters["force_sizes"]["force_min"].GetBool()) {
        if ( MMG3D_Set_dparameter(mmgMesh,mmgSol,MMG3D_DPARAM_hmin, mThisParameters["force_sizes"]["minimal_size"].GetDouble()) != 1 )
            KRATOS_ERROR << "Unable to set the minimal edge size " << std::endl;
    }

    // Minimal edge size
    if (mThisParameters["force_sizes"]["force_max"].GetBool()) {
        if ( MMG3D_Set_dparameter(mmgMesh,mmgSol,MMG3D_DPARAM_hmax, mThisParameters["force_sizes"]["maximal_size"].GetDouble()) != 1 ) {
            KRATOS_ERROR << "Unable to set the maximal edge size " << std::endl;
        }
    }

    const int ier = MMG3D_mmg3dlib(mmgMesh, mmgSol);

    if ( ier == MMG5_STRONGFAILURE )
        KRATOS_ERROR << "ERROR: BAD ENDING OF MMG3DLIB: UNABLE TO SAVE MESH. ier: " << ier << std::endl;
    else if ( ier == MMG5_LOWFAILURE )
        KRATOS_ERROR << "ERROR: BAD ENDING OF MMG3DLIB. ier: " << ier << std::endl;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMGS>::MMGLibCall()
{
    KRATOS_TRY;

    /* Advanced configurations */
    // Global hausdorff value (default value = 0.01) applied on the whole boundary
    if (mThisParameters["advanced_parameters"]["force_hausdorff_value"].GetBool()) {
        if ( MMGS_Set_dparameter(mmgMesh,mmgSol,MMGS_DPARAM_hausd, mThisParameters["advanced_parameters"]["hausdorff_value"].GetDouble()) != 1 )
            KRATOS_ERROR << "Unable to set the Hausdorff parameter" << std::endl;
    }

    // Avoid/allow point relocation
    if ( MMGS_Set_iparameter(mmgMesh,mmgSol,MMGS_IPARAM_nomove, static_cast<int>(mThisParameters["advanced_parameters"]["no_move_mesh"].GetBool())) != 1 )
        KRATOS_ERROR << "Unable to fix the nodes" << std::endl;

    // Don't insert nodes on mesh
    if ( MMGS_Set_iparameter(mmgMesh,mmgSol,MMGS_IPARAM_noinsert, static_cast<int>(mThisParameters["advanced_parameters"]["no_insert_mesh"].GetBool())) != 1 )
        KRATOS_ERROR << "Unable to set no insertion/suppression point" << std::endl;

    // Don't swap mesh
    if ( MMGS_Set_iparameter(mmgMesh,mmgSol,MMGS_IPARAM_noswap, static_cast<int>(mThisParameters["advanced_parameters"]["no_swap_mesh"].GetBool())) != 1 )
        KRATOS_ERROR << "Unable to set no edge flipping" << std::endl;

    // Set the angle detection
    const bool deactivate_detect_angle = mThisParameters["advanced_parameters"]["deactivate_detect_angle"].GetBool();
    if ( deactivate_detect_angle) {
        if ( MMGS_Set_iparameter(mmgMesh,mmgSol,MMGS_IPARAM_angle, static_cast<int>(!deactivate_detect_angle)) != 1 )
            KRATOS_ERROR << "Unable to set the angle detection on" << std::endl;
    }

    // Set the gradation
    if (mThisParameters["advanced_parameters"]["force_gradation_value"].GetBool()) {
        if ( MMGS_Set_dparameter(mmgMesh,mmgSol,MMGS_DPARAM_hgrad, mThisParameters["advanced_parameters"]["gradation_value"].GetDouble()) != 1 )
            KRATOS_ERROR << "Unable to set gradation" << std::endl;
    }

    // Minimal edge size
    if (mThisParameters["force_sizes"]["force_min"].GetBool()) {
        if ( MMGS_Set_dparameter(mmgMesh,mmgSol,MMGS_DPARAM_hmin, mThisParameters["force_sizes"]["minimal_size"].GetDouble()) != 1 )
            KRATOS_ERROR << "Unable to set the minimal edge size " << std::endl;
    }

    // Minimal edge size
    if (mThisParameters["force_sizes"]["force_max"].GetBool()) {
        if ( MMGS_Set_dparameter(mmgMesh,mmgSol,MMGS_DPARAM_hmax, mThisParameters["force_sizes"]["maximal_size"].GetDouble()) != 1 ) {
            KRATOS_ERROR << "Unable to set the maximal edge size " << std::endl;
        }
    }

    const int ier = MMGS_mmgslib(mmgMesh, mmgSol);

    if ( ier == MMG5_STRONGFAILURE )
        KRATOS_ERROR << "ERROR: BAD ENDING OF MMGSLIB: UNABLE TO SAVE MESH. ier: " << ier << std::endl;
    else if ( ier == MMG5_LOWFAILURE )
        KRATOS_ERROR << "ERROR: BAD ENDING OF MMGSLIB. ier: " << ier << std::endl;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMG2D>::FreeAll()
{
    MMG2D_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,MMG5_ARG_end);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMG3D>::FreeAll()
{
    MMG3D_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,MMG5_ARG_end);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMGS>::FreeAll()
{
    MMGS_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,MMG5_ARG_end);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMG2D>::SetNodes(
    const double X,
    const double Y,
    const double Z,
    const IndexType Color,
    const IndexType Index
    )
{
    if ( MMG2D_Set_vertex(mmgMesh, X, Y, Color, Index) != 1 )
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMG3D>::SetNodes(
    const double X,
    const double Y,
    const double Z,
    const IndexType Color,
    const IndexType Index
    )
{
    if ( MMG3D_Set_vertex(mmgMesh, X, Y, Z, Color, Index) != 1 )
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMGS>::SetNodes(
    const double X,
    const double Y,
    const double Z,
    const IndexType Color,
    const IndexType Index
    )
{
    if ( MMGS_Set_vertex(mmgMesh, X, Y, Z, Color, Index) != 1 )
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMG2D>::SetConditions(
    GeometryType & Geom,
    const IndexType Color,
    const IndexType Index
    )
{
    if (Geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Point2D) // Point
        KRATOS_ERROR << "ERROR:: Nodal condition, will be meshed with the node. Condition existence after meshing not guaranteed" << std::endl;
    else if (Geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Line2D2) { // Line
        const IndexType id_1 = Geom[0].Id(); // First node id
        const IndexType id_2 = Geom[1].Id(); // Second node id

        if ( MMG2D_Set_edge(mmgMesh, id_1, id_2, Color, Index) != 1 )
            exit(EXIT_FAILURE);

        // Set fixed boundary
        bool blocked_1 = false;
        if (Geom[0].IsDefined(BLOCKED))
            blocked_1 = Geom[0].Is(BLOCKED);
        bool blocked_2 = false;
        if (Geom[1].IsDefined(BLOCKED))
            blocked_2 = Geom[1].Is(BLOCKED);

        if ((blocked_1 && blocked_2))
            if ( MMG2D_Set_requiredEdge(mmgMesh, Index) != 1 )
                exit(EXIT_FAILURE);
    } else {
        const IndexType size_geometry = Geom.size();
        KRATOS_ERROR << "ERROR: I DO NOT KNOW WHAT IS THIS. Size: " << size_geometry << " Type: " << Geom.GetGeometryType() << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMG3D>::SetConditions(
    GeometryType & Geom,
    const IndexType Color,
    const IndexType Index
    )
{
    if (Geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Point3D) // Point
        KRATOS_ERROR << "ERROR:: Nodal condition, will be meshed with the node. Condition existence after meshing not guaranteed" << std::endl;
    else if (Geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Line3D2) { // Line
        KRATOS_ERROR << "Kratos_Line3D2 remeshing pending to be implemented" << std::endl;
//         const IndexType id1 = Geom[0].Id(); // First node id
//         const IndexType id2 = Geom[1].Id(); // Second node id
//
//         if ( MMG3D_Set_edge(mmgMesh, id1, id2, Color, Index) != 1 )
//             exit(EXIT_FAILURE);
//
//         // Set fixed boundary
//         bool blocked_1 = false;
//         if (Geom[0].IsDefined(BLOCKED))
//             blocked_1 = Geom[0].Is(BLOCKED);
//         bool blocked_2 = false;
//         if (Geom[1].IsDefined(BLOCKED))
//             blocked_2 = Geom[1].Is(BLOCKED);
//
//         if ((blocked_1 && blocked_2))
//             if ( MMG3D_Set_requiredEdge(mmgMesh, Index) != 1 )
//                 exit(EXIT_FAILURE);
    } else if (Geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle3D3) {// Triangle
        const IndexType id_1 = Geom[0].Id(); // First node Id
        const IndexType id_2 = Geom[1].Id(); // Second node Id
        const IndexType id_3 = Geom[2].Id(); // Third node Id

        if ( MMG3D_Set_triangle(mmgMesh, id_1, id_2, id_3, Color, Index) != 1 )
            exit(EXIT_FAILURE);

        // Set fixed boundary
        bool blocked_1 = false;
        if (Geom[0].IsDefined(BLOCKED))
            blocked_1 = Geom[0].Is(BLOCKED);
        bool blocked_2 = false;
        if (Geom[1].IsDefined(BLOCKED))
            blocked_2 = Geom[1].Is(BLOCKED);
        bool blocked_3 = false;
        if (Geom[2].IsDefined(BLOCKED))
            blocked_3 = Geom[2].Is(BLOCKED);

        if ((blocked_1 && blocked_2 && blocked_3))
            if ( MMG3D_Set_requiredTriangle(mmgMesh, Index) != 1 )
                exit(EXIT_FAILURE);
    } else if (Geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4) { // Quadrilaterals
        const IndexType id_1 = Geom[0].Id(); // First node Id
        const IndexType id_2 = Geom[1].Id(); // Second node Id
        const IndexType id_3 = Geom[2].Id(); // Third node Id
        const IndexType id_4 = Geom[3].Id(); // Fourth node Id

        if ( MMG3D_Set_quadrilateral(mmgMesh, id_1, id_2, id_3, id_4, Color, Index) != 1 )
            exit(EXIT_FAILURE);
    } else {
        const SizeType size_geometry = Geom.size();
        KRATOS_ERROR << "ERROR: I DO NOT KNOW WHAT IS THIS. Size: " << size_geometry << " Type: " << Geom.GetGeometryType() << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMGS>::SetConditions(
    GeometryType & Geom,
    const IndexType Color,
    const IndexType Index
    )
{
    if (Geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Point3D) // Point
        KRATOS_ERROR << "ERROR:: Nodal condition, will be meshed with the node. Condition existence after meshing not guaranteed" << std::endl;
    else if (Geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Line3D2) { // Line
        const IndexType id_1 = Geom[0].Id(); // First node id
        const IndexType id_2 = Geom[1].Id(); // Second node id

        if ( MMGS_Set_edge(mmgMesh, id_1, id_2, Color, Index) != 1 )
            exit(EXIT_FAILURE);

        // Set fixed boundary
        bool blocked_1 = false;
        if (Geom[0].IsDefined(BLOCKED))
            blocked_1 = Geom[0].Is(BLOCKED);
        bool blocked_2 = false;
        if (Geom[1].IsDefined(BLOCKED))
            blocked_2 = Geom[1].Is(BLOCKED);

        if ((blocked_1 && blocked_2))
            if ( MMGS_Set_requiredEdge(mmgMesh, Index) != 1 )
                exit(EXIT_FAILURE);
    } else {
        const IndexType size_geometry = Geom.size();
        KRATOS_ERROR << "ERROR: I DO NOT KNOW WHAT IS THIS. Size: " << size_geometry << " Type: " << Geom.GetGeometryType() << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMG2D>::SetElements(
    GeometryType & Geom,
    const IndexType Color,
    const IndexType Index
    )
{
    const IndexType id_1 = Geom[0].Id(); // First node Id
    const IndexType id_2 = Geom[1].Id(); // Second node Id
    const IndexType id_3 = Geom[2].Id(); // Third node Id

    if ( MMG2D_Set_triangle(mmgMesh, id_1, id_2, id_3, Color, Index) != 1 )
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMG3D>::SetElements(
    GeometryType & Geom,
    const IndexType Color,
    const IndexType Index
    )
{
    const IndexType id_1 = Geom[0].Id(); // First node Id
    const IndexType id_2 = Geom[1].Id(); // Second node Id
    const IndexType id_3 = Geom[2].Id(); // Third node Id
    const IndexType id_4 = Geom[3].Id(); // Fourth node Id

    if (Geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4) { // Tetrahedron
        if ( MMG3D_Set_tetrahedron(mmgMesh, id_1, id_2, id_3, id_4, Color, Index) != 1 )
            exit(EXIT_FAILURE);
    } else if (Geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Prism3D6) { // Prisms
        const IndexType id_5 = Geom[4].Id(); // 5th node Id
        const IndexType id_6 = Geom[5].Id(); // 6th node Id

        if ( MMG3D_Set_prism(mmgMesh, id_1, id_2, id_3, id_4, id_5, id_6, Color, Index) != 1 )
            exit(EXIT_FAILURE);
    } else if (Geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Hexahedra3D8) { // Hexaedron
//         const IndexType id_5 = Geom[4].Id(); // 5th node Id
//         const IndexType id_6 = Geom[5].Id(); // 6th node Id
//         const IndexType id_6 = Geom[7].Id(); // 7th node Id
//         const IndexType id_6 = Geom[8].Id(); // 8th node Id

        const SizeType size_geometry = Geom.size();
        KRATOS_ERROR << "ERROR: HEXAEDRON NON IMPLEMENTED IN THE LIBRARY " << size_geometry << std::endl;
    } else {
        const SizeType size_geometry = Geom.size();
        KRATOS_ERROR << "ERROR: I DO NOT KNOW WHAT IS THIS. Size: " << size_geometry << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMGS>::SetElements(
    GeometryType & Geom,
    const IndexType Color,
    const IndexType Index
    )
{
    const IndexType id_1 = Geom[0].Id(); // First node Id
    const IndexType id_2 = Geom[1].Id(); // Second node Id
    const IndexType id_3 = Geom[2].Id(); // Third node Id

    if ( MMGS_Set_triangle(mmgMesh, id_1, id_2, id_3, Color, Index) != 1 )
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMG2D>::SetMetricScalar(
    const double Metric,
    const IndexType NodeId
    )
{
    if ( MMG2D_Set_scalarSol(mmgSol, Metric, NodeId) != 1 )
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMG3D>::SetMetricScalar(
    const double Metric,
    const IndexType NodeId
    )
{
    if ( MMG3D_Set_scalarSol(mmgSol, Metric, NodeId) != 1 )
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMGS>::SetMetricScalar(
    const double Metric,
    const IndexType NodeId
    )
{
    if ( MMGS_Set_scalarSol(mmgSol, Metric, NodeId) != 1 )
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMG2D>::SetMetricVector(
    const array_1d<double, 2>& Metric,
    const IndexType NodeId
    )
{
    if ( MMG2D_Set_vectorSol(mmgSol, Metric[0], Metric[1], NodeId) != 1 )
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMG3D>::SetMetricVector(
    const array_1d<double, 3>& Metric,
    const IndexType NodeId
    )
{
    if ( MMG3D_Set_vectorSol(mmgSol, Metric[0], Metric[1], Metric[2], NodeId) != 1 )
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMGS>::SetMetricVector(
    const array_1d<double, 3>& Metric,
    const IndexType NodeId
    )
{
    if ( MMGS_Set_vectorSol(mmgSol, Metric[0], Metric[1], Metric[2], NodeId) != 1 )
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMG2D>::SetMetricTensor(
    const array_1d<double, 3>& Metric,
    const IndexType NodeId
    )
{
    // The order is XX, XY, YY
    if ( MMG2D_Set_tensorSol(mmgSol, Metric[0], Metric[2], Metric[1], NodeId) != 1 )
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMG3D>::SetMetricTensor(
    const array_1d<double, 6>& Metric,
    const IndexType NodeId
    )
{
    // The order is XX, XY, XZ, YY, YZ, ZZ
    if ( MMG3D_Set_tensorSol(mmgSol, Metric[0], Metric[3], Metric[5], Metric[1], Metric[4], Metric[2], NodeId) != 1 )
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MmgProcess<MMGLibray::MMGS>::SetMetricTensor(
    const array_1d<double, 6>& Metric,
    const IndexType NodeId
    )
{
    // The order is XX, XY, XZ, YY, YZ, ZZ
    if ( MMGS_Set_tensorSol(mmgSol, Metric[0], Metric[3], Metric[5], Metric[1], Metric[4], Metric[2], NodeId) != 1 )
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibray TMMGLibray>
void MmgProcess<TMMGLibray>::CreateAuxiliarSubModelPartForFlags()
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

template<MMGLibray TMMGLibray>
void MmgProcess<TMMGLibray>::AssignAndClearAuxiliarSubModelPartForFlags()
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
void MmgProcess<TMMGLibray>::CreateDebugPrePostRemeshOutput(ModelPart& rOldModelPart)
{
    ModelPart auxiliar_model_part("auxiliar_thing");

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
    ModelPart copy_old_model_part("old_model_part_copy");
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
}

/***********************************************************************************/
/***********************************************************************************/

template class MmgProcess<MMGLibray::MMG2D>;
template class MmgProcess<MMGLibray::MMG3D>;
template class MmgProcess<MMGLibray::MMGS>;

}// namespace Kratos.
