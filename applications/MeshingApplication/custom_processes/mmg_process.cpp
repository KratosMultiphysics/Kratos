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
// We indlude the internal variable interpolation process
#include "custom_processes/nodal_values_interpolation_process.h"
#include "custom_processes/internal_variables_interpolation_process.h"
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

    mFilename = mThisParameters["filename"].GetString();
    mEchoLevel = mThisParameters["echo_level"].GetInt();

    // The framework type
    mFramework = ConvertFramework(mThisParameters["framework"].GetString());

    // The discretization type
    mDiscretization = ConvertDiscretization(mThisParameters["discretization_type"].GetString());

    if ( mDiscretization == DiscretizationOption::ISOSURFACE ){
        mRemoveRegions = mThisParameters["isosurface_parameters"]["remove_regions"].GetBool();
    } else {
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

    /* We print one important information message */
    KRATOS_INFO_IF("MmgProcess", mEchoLevel > 0) << "We clone the first condition and element of each type (we will assume that each sub model part has just one kind of condition, in my opinion it is quite reccomended to create more than one sub model part if you have more than one element or condition)" << std::endl;

    if( mRemoveRegions ){
        // The conditions are re-creted in the process
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
    "//---------------  BEFORE REMESHING   ---------------//" << std::endl <<
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
    "//---------------   AFTER REMESHING   ---------------//" << std::endl <<
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

    // Nodes not belonging to an element are removed
    if(mRemoveRegions) {
        CleanSuperfluousNodes();
    }

    // Save the mesh in an .mdpa format
    const bool save_mdpa_file = mThisParameters["save_mdpa_file"].GetBool();
    if(save_mdpa_file) OutputMdpa();

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
    if (mThisParameters["preserve_flags"].GetBool()) {
        mMmmgUtilities.CreateAuxiliarSubModelPartForFlags(mrThisModelPart);
    }

    // The auxiliar color maps
    ColorsMapType aux_ref_cond, aux_ref_elem;

    // We initialize the mesh data with the given modelpart
    mMmmgUtilities.GenerateMeshDataFromModelPart(mrThisModelPart, mColors, aux_ref_cond, aux_ref_elem, mFramework);

    // Iterate over components
    auto& r_nodes_array = mrThisModelPart.Nodes();

    // We copy the DOF from the first node (after we release, to avoid problem with previous conditions)
    mDofs = r_nodes_array.begin()->GetDofs();
    for (auto it_dof = mDofs.begin(); it_dof != mDofs.end(); ++it_dof)
        it_dof->FreeDof();

    // Generate the maps of reference
    mMmmgUtilities.GenerateReferenceMaps(mrThisModelPart, aux_ref_cond, aux_ref_elem, mpRefCondition, mpRefElement);
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibrary TMMGLibrary>
void MmgProcess<TMMGLibrary>::InitializeSolDataMetric()
{
    // We initialize the solution data with the given modelpart
    mMmmgUtilities.GenerateSolDataFromModelPart(mrThisModelPart);
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibrary TMMGLibrary>
void MmgProcess<TMMGLibrary>::InitializeSolDataDistance()
{
    ////////* SOLUTION FILE for ISOSURFACE*////////
    // Iterate in the nodes
    auto& r_nodes_array = mrThisModelPart.Nodes();
    const auto it_node_begin = r_nodes_array.begin();

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
        auto it_node = it_node_begin + i;

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
    const SizeType buffer_size = mrThisModelPart.NodesBegin()->GetBufferSize();

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
    Model& r_owner_model = mrThisModelPart.GetModel();
    ModelPart& r_old_model_part = r_owner_model.CreateModelPart(mrThisModelPart.Name()+"_Old", mrThisModelPart.GetBufferSize());

    // First we empty the model part
    auto& r_nodes_array = mrThisModelPart.Nodes();

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i)
        (r_nodes_array.begin() + i)->Set(TO_ERASE, true);
    r_old_model_part.AddNodes( mrThisModelPart.NodesBegin(), mrThisModelPart.NodesEnd() );
    mrThisModelPart.RemoveNodesFromAllLevels(TO_ERASE);

    auto& r_conditions_array = mrThisModelPart.Conditions();

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(r_conditions_array.size()); ++i)
        (r_conditions_array.begin() + i)->Set(TO_ERASE, true);
    r_old_model_part.AddConditions( mrThisModelPart.ConditionsBegin(), mrThisModelPart.ConditionsEnd() );
    mrThisModelPart.RemoveConditionsFromAllLevels(TO_ERASE);

    auto& r_elements_array = mrThisModelPart.Elements();

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(r_elements_array.size()); ++i)
        (r_elements_array.begin() + i)->Set(TO_ERASE, true);
    r_old_model_part.AddElements( mrThisModelPart.ElementsBegin(), mrThisModelPart.ElementsEnd() );
    mrThisModelPart.RemoveElementsFromAllLevels(TO_ERASE);

    // Writing the new mesh data on the model part
    mMmmgUtilities.WriteMeshDataToModelPart(mrThisModelPart, mColors, mDofs, mmg_mesh_info, mpRefCondition, mpRefElement);

    // Writing the new solution data on the model part
    mMmmgUtilities.WriteSolDataToModelPart(mrThisModelPart);

    /* After that we reorder nodes, conditions and elements: */
    mMmmgUtilities.ReorderAllIds(mrThisModelPart);

    /* We assign flags and clear the auxiliar model parts created to reassing the flags */
    if (mThisParameters["preserve_flags"].GetBool()) {
        mMmmgUtilities.AssignAndClearAuxiliarSubModelPartForFlags(mrThisModelPart);
    }

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
    if (mThisParameters["initialize_entities"].GetBool()) {
        InitializeElementsAndConditions();
    }

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
    r_owner_model.DeleteModelPart(mrThisModelPart.Name()+"_Old");

    /* We clean conditions with duplicated geometries (this is an error on fluid simulations) */
    if (mFramework == FrameworkEulerLagrange::EULERIAN) {
        ClearConditionsDuplicatedGeometries();
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibrary TMMGLibrary>
void MmgProcess<TMMGLibrary>::InitializeElementsAndConditions()
{
    // Iterate over conditions
    auto& r_conditions_array = mrThisModelPart.Conditions();
    const auto it_cond_begin = r_conditions_array.begin();
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(r_conditions_array.size()); ++i)
        (it_cond_begin + i)->Initialize();

    // Iterate over elements
    auto& r_elements_array = mrThisModelPart.Elements();
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
    const int step = mrThisModelPart.GetProcessInfo()[STEP];
    const std::string file_name = mFilename + "_step=" + std::to_string(step) + (PostOutput ? ".o" : "");

    // Automatically save the mesh
    mMmmgUtilities.OutputMesh(file_name);

    // Automatically save the solution
    mMmmgUtilities.OutputSol(file_name);

    if (mThisParameters["save_colors_files"].GetBool()) {
        // Output the reference files
        mMmmgUtilities.OutputReferenceEntitities(file_name, mpRefCondition, mpRefElement);

        // Writing the colors to a JSON
        AssignUniqueModelPartCollectionTagUtility::WriteTagsToJson(file_name, mColors);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibrary TMMGLibrary>
void MmgProcess<TMMGLibrary>::FreeMemory()
{
    // Free the MMG structures
    mMmmgUtilities.FreeAll();

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
void MmgProcess<TMMGLibrary>::ClearConditionsDuplicatedGeometries()
{
    // Next check that the conditions are oriented accordingly to do so begin by putting all of the conditions in a set
    typedef std::unordered_map<DenseVector<IndexType>, std::vector<IndexType>, KeyHasherRange<DenseVector<IndexType>>, KeyComparorRange<DenseVector<IndexType>> > HashMapType;
    HashMapType faces_map;

    // Iterate over conditions
    auto& r_conditions_array = mrThisModelPart.Conditions();

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
    for (auto& r_face : faces_map) {
        const auto& r_pairs = r_face.second;
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
    Model& r_owner_model = mrThisModelPart.GetModel();
    ModelPart& r_auxiliar_model_part = r_owner_model.CreateModelPart(mrThisModelPart.Name()+"_Auxiliar", mrThisModelPart.GetBufferSize());
    ModelPart& r_copy_old_model_part = r_owner_model.CreateModelPart(mrThisModelPart.Name()+"_Old_Copy", mrThisModelPart.GetBufferSize());

    Properties::Pointer p_prop_1 = r_auxiliar_model_part.pGetProperties(1);
    Properties::Pointer p_prop_2 = r_auxiliar_model_part.pGetProperties(2);

    // We just transfer nodes and elements
    // Current model part
    FastTransferBetweenModelPartsProcess transfer_process_current = FastTransferBetweenModelPartsProcess(r_auxiliar_model_part, mrThisModelPart, FastTransferBetweenModelPartsProcess::EntityTransfered::NODESANDELEMENTS);
    transfer_process_current.Set(MODIFIED); // We replicate, not transfer
    transfer_process_current.Execute();

    // Iterate over first elements
    auto& r_elements_array_1 = r_auxiliar_model_part.Elements();
    const auto it_elem_begin_1 = r_elements_array_1.begin();

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(r_elements_array_1.size()); ++i) {
        auto it_elem = it_elem_begin_1 + i;
        it_elem->SetProperties(p_prop_1);
    }
    // Old model part
    FastTransferBetweenModelPartsProcess transfer_process_old = FastTransferBetweenModelPartsProcess(r_copy_old_model_part, rOldModelPart, FastTransferBetweenModelPartsProcess::EntityTransfered::NODESANDELEMENTS);
    transfer_process_current.Set(MODIFIED); // We replicate, not transfer
    transfer_process_old.Execute();

    // Iterate over second elements
    auto& r_elements_array_2 = r_copy_old_model_part.Elements();
    const auto it_elem_begin_2 = r_elements_array_2.begin();

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(r_elements_array_2.size()); ++i) {
        auto it_elem = it_elem_begin_2 + i;
        it_elem->SetProperties(p_prop_2);
    }

    // Reorder ids to ensure be consecuent
    auto& r_auxiliar_nodes_array = r_auxiliar_model_part.Nodes();
    const SizeType auxiliar_number_of_nodes = (r_auxiliar_nodes_array.end() - 1)->Id();
    auto& r_copy_old_nodes_array = r_copy_old_model_part.Nodes();

    for(IndexType i = 0; i < r_copy_old_nodes_array.size(); ++i) {
        auto it_node = r_copy_old_nodes_array.begin() + i;
        it_node->SetId(auxiliar_number_of_nodes + i + 1);
    }

    // Last transfer
    FastTransferBetweenModelPartsProcess transfer_process_last = FastTransferBetweenModelPartsProcess(r_auxiliar_model_part, r_copy_old_model_part, FastTransferBetweenModelPartsProcess::EntityTransfered::NODESANDELEMENTS);
    transfer_process_last.Set(MODIFIED);
    transfer_process_last.Execute();

    const int step = mrThisModelPart.GetProcessInfo()[STEP];
    const double label = static_cast<double>(step);
    GidIO<> gid_io("BEFORE_AND_AFTER_MMG_MESH_STEP=" + std::to_string(step), GiD_PostBinary, SingleFile, WriteUndeformed,  WriteElementsOnly);

    gid_io.InitializeMesh(label);
    gid_io.WriteMesh(r_auxiliar_model_part.GetMesh());
    gid_io.FinalizeMesh();
    gid_io.InitializeResults(label, r_auxiliar_model_part.GetMesh());

    // Remove auxiliar model parts
    r_owner_model.DeleteModelPart(mrThisModelPart.Name()+"_Auxiliar");
    r_owner_model.DeleteModelPart(mrThisModelPart.Name()+"_Old_Copy");
}

/***********************************************************************************/
/***********************************************************************************/

template<MMGLibrary TMMGLibrary>
void MmgProcess<TMMGLibrary>::CleanSuperfluousNodes()
{
    // Iterate over nodes
    auto& r_nodes_array = mrThisModelPart.Nodes();
    const SizeType initial_num = r_nodes_array.size();

    // Marking all nodes as "superfluous"
    VariableUtils().SetFlag(TO_ERASE, true, r_nodes_array);

    // Iterate over elements
    const auto& r_elements_array = mrThisModelPart.Elements();
    const auto it_elem_begin = r_elements_array.begin();

    // Saving the nodes that belong to an element
    #pragma omp parallel for
    for(int i_elem = 0; i_elem < static_cast<int>(r_elements_array.size()); ++i_elem ){

        const auto it_elem = it_elem_begin + i_elem;
        auto& r_geom = it_elem->GetGeometry();

        for (IndexType i = 0; i < r_geom.size(); ++i){
            r_geom[i].Set(TO_ERASE, false);
        }
    }

    mrThisModelPart.RemoveNodesFromAllLevels(TO_ERASE);
    const SizeType final_num = mrThisModelPart.Nodes().size();
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
        "save_colors_files"                    : false,
        "save_mdpa_file"                       : false,
        "max_number_of_searchs"                : 1000,
        "preserve_flags"                       : true,
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
