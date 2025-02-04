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
#include "custom_processes/parmmg/pmmg_process.h"
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
    ): BaseType(&rThisModelPart)
{
    const Parameters default_parameters = GetDefaultParameters();
    ThisParameters.RecursivelyValidateAndAssignDefaults(default_parameters);
    mThisParameters = ThisParameters;

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
    mPMmgUtilities.InitMesh(mrThisModelPart.GetCommunicator().GetDataCommunicator());

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

    /* We print the original model part */
    KRATOS_INFO_IF("", mEchoLevel > 0) <<
    "//---------------------------------------------------//" << std::endl <<
    "//---------------  BEFORE REMESHING   ---------------//" << std::endl <<
    "//---------------------------------------------------//" << std::endl <<
    std::endl << mrThisModelPart << std::endl;

    // We initialize the mesh and solution data
    InitializeMeshData();

    InitializeSolDataMetric();

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
            // We get the metric
            TensorArrayType post_mesh_metric;
            mPMmgUtilities.GetMetricTensor(post_mesh_metric);

            // We set the metric
            it_node->SetValue(r_tensor_variable, post_mesh_metric);
        }
    }

    // We release the memory
    FreeMemory();

    // Call parallel fill communicator
    ParallelFillCommunicator(mrThisModelPart, mrThisModelPart.GetCommunicator().GetDataCommunicator()).Execute();

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
    KRATOS_TRY;

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

    KRATOS_CATCH("");

}

template<PMMGLibrary TPMMGLibrary>
template<typename TPointerType>
void ParMmgProcess<TPMMGLibrary>::SyncMapAcrossRanks(std::unordered_map<IndexType, TPointerType>& rRefEntityMap)
{
    KRATOS_TRY;

    const IndexType size = mrThisModelPart.GetCommunicator().GetDataCommunicator().Size();
    const IndexType rank = mrThisModelPart.GetCommunicator().GetDataCommunicator().Rank();

    KRATOS_ERROR_IF(size < 2) << "Using ParMmg Process in a serial run! The MPI size is: " << size;

    if (rank > 0) {
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

    if (rank == 0) {
        for (IndexType i_rank = 1; i_rank < size; i_rank++) {
            mrThisModelPart.GetCommunicator().GetDataCommunicator().Send(rRefEntityMap, i_rank);
        }
    }
    else {
        std::unordered_map<IndexType, TPointerType> aux_recv_map;
        mrThisModelPart.GetCommunicator().GetDataCommunicator().Recv(aux_recv_map, 0);
        rRefEntityMap = aux_recv_map;
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
void ParMmgProcess<TPMMGLibrary>::InitializeSolDataMetric()
{
    KRATOS_TRY;

    // We initialize the solution data with the given modelpart
    mPMmgUtilities.GenerateSolDataFromModelPart(mrThisModelPart);

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
void ParMmgProcess<TPMMGLibrary>::ExecuteRemeshing()
{
    KRATOS_TRY;

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

    /* Save to file. It is required to call it AFTER writing the model part data */
    if (save_to_file) SaveSolutionToFile(false);

    // Calling the library functions
    mPMmgUtilities.PMMGLibCallMetric(mThisParameters);

    // Some information
    PMMGMeshInfo<TPMMGLibrary> mmg_mesh_info;
    mPMmgUtilities.PrintAndGetParMmgMeshInfo(mmg_mesh_info);


    // First we empty the model part
    auto& r_nodes_array = mrThisModelPart.Nodes();

    block_for_each(r_nodes_array,
        [&](NodeType& rNode) {

        const bool old_entity = rNode.IsDefined(OLD_ENTITY) ? rNode.Is(OLD_ENTITY) : false;
        if (!old_entity) {
            rNode.Set(TO_ERASE, true);
        }
    });
    r_old_model_part.AddNodes( mrThisModelPart.NodesBegin(), mrThisModelPart.NodesEnd() );
    mrThisModelPart.RemoveNodesFromAllLevels(TO_ERASE);

    auto& r_conditions_array = mrThisModelPart.Conditions();

    block_for_each(r_conditions_array,
        [&](Condition& rCondition) {

        const bool old_entity = rCondition.IsDefined(OLD_ENTITY) ? rCondition.Is(OLD_ENTITY) : false;
        if (!old_entity) {
            rCondition.Set(TO_ERASE, true);
        }
    });
    r_old_model_part.AddConditions( mrThisModelPart.ConditionsBegin(), mrThisModelPart.ConditionsEnd() );
    mrThisModelPart.RemoveConditionsFromAllLevels(TO_ERASE);

    auto& r_elements_array = mrThisModelPart.Elements();

    block_for_each(r_elements_array,
        [&](Element& rElement) {

        const bool old_entity = rElement.IsDefined(OLD_ENTITY) ? rElement.Is(OLD_ENTITY) : false;
        if (!old_entity) {
            rElement.Set(TO_ERASE, true);
        }
    });
    r_old_model_part.AddElements( mrThisModelPart.ElementsBegin(), mrThisModelPart.ElementsEnd() );
    mrThisModelPart.RemoveElementsFromAllLevels(TO_ERASE);


    mPMmgUtilities.WriteMeshDataToModelPart(mrThisModelPart, mColors, mDofs, mmg_mesh_info, mpRefCondition, mpRefElement);

    /* Save to file. It is required to call it AFTER writing the model part data */
    if (save_to_file) SaveSolutionToFile(true);

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

    KRATOS_CATCH("");

}

/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
void ParMmgProcess<TPMMGLibrary>::InitializeElementsAndConditions()
{
    KRATOS_TRY;

    BaseType::InitializeElementsAndConditions();

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
void ParMmgProcess<TPMMGLibrary>::SaveSolutionToFile(const bool PostOutput)
{
    KRATOS_TRY;

    /* GET RESULTS */
    const int step = mrThisModelPart.GetProcessInfo()[STEP];
    const int rank = mrThisModelPart.GetCommunicator().MyPID();
    const std::string file_name = mFilename + "_step=" + std::to_string(step) + (PostOutput ? ".o" : "");

    if (PostOutput) {
        mPMmgUtilities.OutputMesh(file_name);

        mPMmgUtilities.OutputSol(file_name);
    }

    if (mThisParameters["save_colors_files"].GetBool()) {
        // Output the reference files
        mPMmgUtilities.OutputReferenceEntitities(file_name+"_"+std::to_string(rank), mpRefCondition, mpRefElement);

        // Writing the colors to a JSON
        AssignUniqueModelPartCollectionTagUtility::WriteTagsToJson(file_name+"_"+std::to_string(rank), mColors);
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
void ParMmgProcess<TPMMGLibrary>::FreeMemory()
{
    KRATOS_TRY;

    // Free the PMMG structures
    mPMmgUtilities.FreeAll();

    // Free reference std::unordered_map
    mpRefElement.clear();
    mpRefCondition.clear();

    // Clear the colors
    mColors.clear();

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
void ParMmgProcess<TPMMGLibrary>::OutputMdpa()
{
    KRATOS_TRY;

    const auto rank = mrThisModelPart.GetCommunicator().MyPID();
    ModelPartIO model_part_io(mFilename+"_"+std::to_string(rank), IO::WRITE);
    model_part_io.WriteModelPart(mrThisModelPart);

    KRATOS_CATCH("");
}


/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
void ParMmgProcess<TPMMGLibrary>::ClearConditionsDuplicatedGeometries()
{
    KRATOS_TRY;

    BaseType::ClearConditionsDuplicatedGeometries();

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
void ParMmgProcess<TPMMGLibrary>::CreateDebugPrePostRemeshOutput(ModelPart& rOldModelPart)
{
    KRATOS_TRY;

    BaseType::CreateDebugPrePostRemeshOutput(rOldModelPart);

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<PMMGLibrary TPMMGLibrary>
const Parameters ParMmgProcess<TPMMGLibrary>::GetDefaultParameters() const
{
    KRATOS_TRY;

    Parameters default_parameters = Parameters(R"(
    {
        "advanced_parameters"    :
        {
            "number_of_iterations"     : 4,
            "mesh_size"                : 30000,
            "metis_ratio"              : 82
        }
    })");

    // Getting base class default parameters
    const Parameters base_default_parameters = BaseType::GetDefaultParameters();
    default_parameters.RecursivelyAddMissingParameters(base_default_parameters);

    return default_parameters;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template class KRATOS_API(MESHING_APPLICATION) ParMmgProcess<PMMGLibrary::PMMG3D>;

}// namespace Kratos.
