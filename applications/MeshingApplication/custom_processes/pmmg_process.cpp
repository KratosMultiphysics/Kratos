// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:         BSD License
//                   license: MeshingApplication/license.txt
//
//  Main authors:    Marc Nunez
//                   Carlos Roig
//                   Vicente Mataix Ferrandiz
//

// System includes
#include <set>

// External includes

// Project includes
#include "custom_processes/pmmg_process.h"
#include "includes/gid_io.h"
#include "includes/model_part_io.h"
#include "mpi/utilities/parallel_fill_communicator.h"


// NOTE: The following contains the license of the PMMG library
/* =============================================================================
**  Copyright (c) Bx INP/Inria/UBordeaux, 2017- .
**
**  parmmg is free software: you can redistribute it and/or modify it
**  under the terms of the GNU Lesser General Public License as published
**  by the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  parmmg is distributed in the hope that it will be useful, but WITHOUT
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

template<PMMGLibrary TPMMGLibrary>
ParMmgProcess<TPMMGLibrary>::ParMmgProcess(
    ModelPart& rThisModelPart,
    Parameters ThisParameters
    ): MmgProcess(rThisModelPart, ThisParameters)
{
    Parameters default_parameters = GetDefaultParameters();
    mThisParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

    mFilename = mThisParameters["filename"].GetString();
    mEchoLevel = mThisParameters["echo_level"].GetInt();

    // The framework type
    mFramework = ConvertFramework(mThisParameters["framework"].GetString());
    KRATOS_ERROR_IF(mFramework != FrameworkEulerLagrange::EULERIAN) << "Only the eulerian framework is currently supported with ParMmg!" << std::endl;
    // The discretization type
    mDiscretization = ConvertDiscretization(mThisParameters["discretization_type"].GetString());
    KRATOS_ERROR_IF(mDiscretization != DiscretizationOption::STANDARD) << "Only the standard discretization is currently supported with ParMmg!" << std::endl;

    mRemoveRegions = false;

    KRATOS_ERROR_IF(mThisParameters["collapse_prisms_elements"].GetBool()) << "Prisms are not currently supported in ParMmg" << std::endl;

    mpRefElement.clear();
    mpRefCondition.clear();
}

/*************************************** EXECUTE ***********************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
void ParMmgProcess<TPMMGLibrary>::Execute()
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

template<PMMGLibrary TPMMGLibrary>
void ParMmgProcess<TPMMGLibrary>::ExecuteInitialize()
{
    KRATOS_TRY;

    /* We print one important information message */
    KRATOS_INFO_IF("ParMmgProcess", mEchoLevel > 0) << "We clone the first condition and element of each type (we will assume that each sub model part has just one kind of condition, it is quite recommended to create more than one sub model part if you have more than one element or condition)" << std::endl;

    /* We restart the PMMG mesh and solution */
    mPMmgUtilities.SetEchoLevel(mEchoLevel);
    mPMmgUtilities.InitMesh();

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
void ParMmgProcess<TPMMGLibrary>::ExecuteBeforeSolutionLoop()
{
    KRATOS_TRY;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
void ParMmgProcess<TPMMGLibrary>::ExecuteInitializeSolutionStep()
{
    KRATOS_TRY;

    const bool save_to_file = mThisParameters["save_external_files"].GetBool();

    /* We print the original model part */
    KRATOS_INFO_IF("", mEchoLevel > 0) <<
    "//---------------------------------------------------//" << std::endl <<
    "//---------------  BEFORE REMESHING   ---------------//" << std::endl <<
    "//---------------------------------------------------//" << std::endl <<
    std::endl << mrThisModelPart << std::endl;

    // We initialize the mesh and solution data
    InitializeMeshData();

    InitializeSolDataMetric();

    // Save to file
    if (save_to_file) SaveSolutionToFile(false);

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

template<PMMGLibrary TPMMGLibrary>
void ParMmgProcess<TPMMGLibrary>::ExecuteFinalizeSolutionStep()
{
    KRATOS_TRY;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
void ParMmgProcess<TPMMGLibrary>::ExecuteBeforeOutputStep()
{
    KRATOS_TRY;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
void ParMmgProcess<TPMMGLibrary>::ExecuteAfterOutputStep()
{
    KRATOS_TRY;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
void ParMmgProcess<TPMMGLibrary>::ExecuteFinalize()
{
    KRATOS_TRY;

    // Getting new metric (not necessary, only for post-process pourposes)
    /* Tensor variable definition */
    const Variable<TensorArrayType>& r_tensor_variable = KratosComponents<Variable<TensorArrayType>>::Get("METRIC_TENSOR_" + std::to_string(Dimension) + "D");

    // Iterate in the nodes
    auto& r_nodes_array = mrThisModelPart.Nodes();
    auto local_to_global_nodes_map = mPMmgUtilities.GetNodalLocalToGlobalMap();
    for(int i = 1; i <= static_cast<int>(r_nodes_array.size()); ++i) {

        auto it_node = mrThisModelPart.pGetNode(local_to_global_nodes_map[i]);

        const bool old_entity = it_node->IsDefined(OLD_ENTITY) ? it_node->Is(OLD_ENTITY) : false;
        if (!old_entity) {
            KRATOS_DEBUG_ERROR_IF_NOT(it_node->Has(r_tensor_variable)) << "METRIC_TENSOR_" + std::to_string(Dimension) + "D  not defined for node " << it_node->Id() << std::endl;

            // We get the metric
            TensorArrayType& r_metric = it_node->GetValue(r_tensor_variable);

            // We set the metric
            mPMmgUtilities.GetMetricTensor(r_metric);
        }
    }

    // We release the memory
    FreeMemory();

    // Call parallel fill communicator
    mrThisModelPart.GetCommunicator().GetDataCommunicator().Barrier();
    ParallelFillCommunicator(mrThisModelPart).Execute();

    // Save the mesh in an .mdpa format
    const bool save_mdpa_file = mThisParameters["save_mdpa_file"].GetBool();
    if(save_mdpa_file) OutputMdpa();

    KRATOS_CATCH("");
}

/************************************* OPERATOR() **********************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
void ParMmgProcess<TPMMGLibrary>::operator()()
{
    Execute();
}

/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
void ParMmgProcess<TPMMGLibrary>::InitializeMeshData()
{
    // We create a list of submodelparts to later reassign flags after remesh
    if (mThisParameters["preserve_flags"].GetBool()) {
        mPMmgUtilities.CreateAuxiliarSubModelPartForFlags(mrThisModelPart);
    }

    // The auxiliar color maps
    ColorsMapType aux_ref_cond, aux_ref_elem;

    // Actually generate mesh data
    mPMmgUtilities.GenerateMeshDataFromModelPart(mrThisModelPart, mColors, aux_ref_cond, aux_ref_elem, mFramework, false);

    // Generate Interface data
    mPMmgUtilities.GenerateParallelInterfaces(mrThisModelPart);

    // We copy the DOF from the first node (after we release, to avoid problem with previous conditions)
    auto& r_nodes_array = mrThisModelPart.Nodes();
    const auto& r_old_dofs = r_nodes_array.begin()->GetDofs();
    for (auto it_dof = r_old_dofs.begin(); it_dof != r_old_dofs.end(); ++it_dof)
        mDofs.push_back(Kratos::make_unique<NodeType::DofType>(**it_dof));
    for (auto it_dof = mDofs.begin(); it_dof != mDofs.end(); ++it_dof)
        (**it_dof).FreeDof();

    // Generate the maps of reference
    mPMmgUtilities.GenerateReferenceMaps(mrThisModelPart, aux_ref_cond, aux_ref_elem, mpRefCondition, mpRefElement);

    SyncMapAcrossRanks(mpRefCondition);
    SyncMapAcrossRanks(mpRefElement);
}

template<PMMGLibrary TPMMGLibrary>
template<typename TPointerType>
void ParMmgProcess<TPMMGLibrary>::SyncMapAcrossRanks(std::unordered_map<IndexType, TPointerType>& rRefEntityMap)
{
    const IndexType size = mrThisModelPart.GetCommunicator().GetDataCommunicator().Size();
    const IndexType rank = mrThisModelPart.GetCommunicator().GetDataCommunicator().Rank();

    if (rank>0) {
        mrThisModelPart.GetCommunicator().GetDataCommunicator().Send(rRefEntityMap, 0);
    }
    else {
        for (IndexType i_rank = 1; i_rank<size; i_rank++) {
            std::unordered_map<IndexType, TPointerType> aux_recv_map;
            mrThisModelPart.GetCommunicator().GetDataCommunicator().Recv(aux_recv_map, i_rank);
            // Filling map from other ranks information
            for (auto& t : aux_recv_map) {
                if (rRefEntityMap.find(t.first) == rRefEntityMap.end()) {
                    rRefEntityMap[t.first] = t.second;
                }
            }
        }
    }

    if (rank==0) {
        for (IndexType i_rank = 1; i_rank < size; i_rank++) {
            mrThisModelPart.GetCommunicator().GetDataCommunicator().Send(rRefEntityMap, i_rank);
        }
    }
    else {
        std::unordered_map<IndexType, TPointerType> aux_recv_map;
        mrThisModelPart.GetCommunicator().GetDataCommunicator().Recv(aux_recv_map, 0);
        rRefEntityMap = aux_recv_map;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
void ParMmgProcess<TPMMGLibrary>::InitializeSolDataMetric()
{
    // We initialize the solution data with the given modelpart
    mPMmgUtilities.GenerateSolDataFromModelPart(mrThisModelPart);
}

/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
void ParMmgProcess<TPMMGLibrary>::ExecuteRemeshing()
{
    // Getting the parameters
    const bool save_to_file = mThisParameters["save_external_files"].GetBool();

    // We initialize some values
    const SizeType step_data_size = mrThisModelPart.GetNodalSolutionStepDataSize();
    const SizeType buffer_size = mrThisModelPart.NodesBegin()->GetBufferSize();

    mThisParameters["step_data_size"].SetInt(step_data_size);
    mThisParameters["buffer_size"].SetInt(buffer_size);

    KRATOS_INFO_IF("ParMmgProcess", mEchoLevel > 0) << "Step data size: " << step_data_size << " Buffer size: " << buffer_size << std::endl;

    ////////* PMMG LIBRARY CALL *////////
    KRATOS_INFO_IF("ParMmgProcess", mEchoLevel > 0) << "////////* PMMG LIBRARY CALL *////////" << std::endl;

    ////////* EMPTY AND BACKUP THE MODEL PART *////////
    Model& r_owner_model = mrThisModelPart.GetModel();
    ModelPart& r_old_model_part = r_owner_model.CreateModelPart(mrThisModelPart.Name()+"_Old", mrThisModelPart.GetBufferSize());

    // Calling the library functions
    mPMmgUtilities.PMMGLibCallMetric(mThisParameters);

    /* Save to file */
     if (save_to_file) SaveSolutionToFile(true);

    // Some information
    PMMGMeshInfo<TPMMGLibrary> mmg_mesh_info;
    mPMmgUtilities.PrintAndGetParMmgMeshInfo(mmg_mesh_info);


    // First we empty the model part
    auto& r_nodes_array = mrThisModelPart.Nodes();
    const auto it_node_begin = r_nodes_array.begin();

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i) {
        auto it_node = it_node_begin + i;

        const bool old_entity = it_node->IsDefined(OLD_ENTITY) ? it_node->Is(OLD_ENTITY) : false;
        if (!old_entity) {
            it_node->Set(TO_ERASE, true);
        }
    }
    r_old_model_part.AddNodes( mrThisModelPart.NodesBegin(), mrThisModelPart.NodesEnd() );
    mrThisModelPart.RemoveNodesFromAllLevels(TO_ERASE);

    auto& r_conditions_array = mrThisModelPart.Conditions();
    const auto it_cond_begin = r_conditions_array.begin();

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(r_conditions_array.size()); ++i) {
        auto it_cond = it_cond_begin + i;

        const bool old_entity = it_cond->IsDefined(OLD_ENTITY) ? it_cond->Is(OLD_ENTITY) : false;
        if (!old_entity) {
            it_cond->Set(TO_ERASE, true);
        }
    }
    r_old_model_part.AddConditions( mrThisModelPart.ConditionsBegin(), mrThisModelPart.ConditionsEnd() );
    mrThisModelPart.RemoveConditionsFromAllLevels(TO_ERASE);

    auto& r_elements_array = mrThisModelPart.Elements();
    const auto it_elem_begin = r_elements_array.begin();

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(r_elements_array.size()); ++i) {
        auto it_elem = it_elem_begin + i;

        const bool old_entity = it_elem->IsDefined(OLD_ENTITY) ? it_elem->Is(OLD_ENTITY) : false;
        if (!old_entity) {
            it_elem->Set(TO_ERASE, true);
        }
    }
    r_old_model_part.AddElements( mrThisModelPart.ElementsBegin(), mrThisModelPart.ElementsEnd() );
    mrThisModelPart.RemoveElementsFromAllLevels(TO_ERASE);


    mPMmgUtilities.WriteMeshDataToModelPart(mrThisModelPart, mColors, mDofs, mmg_mesh_info, mpRefCondition, mpRefElement);

    // Writing the new solution data on the model part
    mPMmgUtilities.WriteSolDataToModelPart(mrThisModelPart);

    /* We assign flags and clear the auxiliar model parts created to reassing the flags */
    if (mThisParameters["preserve_flags"].GetBool()) {
        mPMmgUtilities.AssignAndClearAuxiliarSubModelPartForFlags(mrThisModelPart);
    }

    // We create an auxiliar mesh for debugging purposes
    if (mThisParameters["debug_result_mesh"].GetBool()) {
        CreateDebugPrePostRemeshOutput(r_old_model_part);
    }

    /* We initialize elements and conditions */
    if (mThisParameters["initialize_entities"].GetBool()) {
        InitializeElementsAndConditions();
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

    // We clear the OLD_ENTITY flag
    VariableUtils().ResetFlag(OLD_ENTITY, mrThisModelPart.Nodes());
    VariableUtils().ResetFlag(OLD_ENTITY, mrThisModelPart.Conditions());
    VariableUtils().ResetFlag(OLD_ENTITY, mrThisModelPart.Elements());
}

/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
void ParMmgProcess<TPMMGLibrary>::InitializeElementsAndConditions()
{
    KRATOS_TRY;

    const auto& r_process_info = mrThisModelPart.GetProcessInfo();

    // Iterate over conditions
    auto& r_conditions_array = mrThisModelPart.Conditions();
    const auto it_cond_begin = r_conditions_array.begin();
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(r_conditions_array.size()); ++i)
        (it_cond_begin + i)->Initialize(r_process_info);

    // Iterate over elements
    auto& r_elements_array = mrThisModelPart.Elements();
    const auto it_elem_begin = r_elements_array.begin();
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(r_elements_array.size()); ++i)
        (it_elem_begin + i)->Initialize(r_process_info);

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
void ParMmgProcess<TPMMGLibrary>::SaveSolutionToFile(const bool PostOutput)
{
    /* GET RESULTS */
    const int step = mrThisModelPart.GetProcessInfo()[STEP];
    // const int rank = mrThisModelPart.GetCommunicator().MyPID();
    const std::string file_name = mFilename + "_step=" + std::to_string(step) + (PostOutput ? ".o" : "");
    // const std::string file_name = mFilename + "_rank_" + std::to_string(rank) + "_step=" + std::to_string(step) + (PostOutput ? ".o" : "");

    // Automatically save the mesh
    if (PostOutput) {
        mPMmgUtilities.OutputMesh(file_name);

        // Automatically save the solution
        mPMmgUtilities.OutputSol(file_name);
    }//

    if (mThisParameters["save_colors_files"].GetBool()) {
        // Output the reference files
        mPMmgUtilities.OutputReferenceEntitities(file_name, mpRefCondition, mpRefElement);

        // Writing the colors to a JSON
        AssignUniqueModelPartCollectionTagUtility::WriteTagsToJson(file_name, mColors);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
void ParMmgProcess<TPMMGLibrary>::FreeMemory()
{
    // Free the PMMG structures
    mPMmgUtilities.FreeAll();

    // Free reference std::unordered_map
    mpRefElement.clear();
    mpRefCondition.clear();

    // Clear the colors
    mColors.clear();
}

/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
void ParMmgProcess<TPMMGLibrary>::OutputMdpa()
{
    const auto rank = mrThisModelPart.GetCommunicator().MyPID();
    std::ofstream output_file;
    ModelPartIO model_part_io("output_"+std::to_string(rank), IO::WRITE);
    model_part_io.WriteModelPart(mrThisModelPart);
}


/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
void ParMmgProcess<TPMMGLibrary>::ClearConditionsDuplicatedGeometries()
{
    // Next check that the conditions are oriented accordingly to do so begin by putting all of the conditions in a set
    typedef std::unordered_map<DenseVector<IndexType>, std::vector<IndexType>, KeyHasherRange<DenseVector<IndexType>>, KeyComparorRange<DenseVector<IndexType>> > HashMapType;
    HashMapType faces_map;

    // Iterate over conditions
    auto& r_conditions_array = mrThisModelPart.Conditions();

    // Reset flag
    VariableUtils().ResetFlag(TO_ERASE, r_conditions_array);

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
                KRATOS_INFO_IF("ParMmgProcess", mEchoLevel > 2) << "Condition created ID:\t" << i_id << " will be removed" << std::endl;
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

template<PMMGLibrary TPMMGLibrary>
void ParMmgProcess<TPMMGLibrary>::CreateDebugPrePostRemeshOutput(ModelPart& rOldModelPart)
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
    GidIO<> gid_io("BEFORE_AND_AFTER_PMMG_MESH_STEP=" + std::to_string(step), GiD_PostBinary, SingleFile, WriteUndeformed,  WriteElementsOnly);

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

template<PMMGLibrary TPMMGLibrary>
void ParMmgProcess<TPMMGLibrary>::MarkConditionsSubmodelParts(ModelPart& rModelPart)
{
    // Iterate over submodelparts
    for (auto& r_sub_model_part : rModelPart.SubModelParts()) {
        VariableUtils().SetFlag(MARKER, true, r_sub_model_part.Conditions());
        MarkConditionsSubmodelParts(r_sub_model_part);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
Parameters ParMmgProcess<TPMMGLibrary>::GetDefaultParameters()
{
    Parameters default_parameters = Parameters(R"(
    {
        "filename"                             : "out",
        "discretization_type"                  : "Standard",
        "isosurface_parameters"                :
        {
            "isosurface_variable"              : "DISTANCE",
            "nonhistorical_variable"           : false,
            "remove_internal_regions"          : false
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
            "gradation_value"                     : 1.3,
            "niter"                               : 4,
            "meshSize"                            : 30000,
            "metisRatio"                          : 82,
            "hgradreq"                            : 5.0,
            "APImode"                             : 0
        },
        "collapse_prisms_elements"             : false,
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

template class ParMmgProcess<PMMGLibrary::PMMG3D>;

}// namespace Kratos.
