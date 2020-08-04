//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Aditya Ghantasala
//
//

// System includes
// External includes

// Project includes
#include "utilities/openmp_utils.h"
#include "custom_processes/sliding_interface_process.h"
#include "utilities/binbased_fast_point_locator_conditions.h"
#include "utilities/geometrical_transformation_utilities.h"
#include "utilities/variable_utils.h"
#include "containers/model.h"
#ifdef KRATOS_USING_MPI
    #include "custom_utilities/gather_modelpart_on_all_ranks.h"
#endif

namespace Kratos
{

SlidingInterfaceProcess::SlidingInterfaceProcess(ModelPart &rMasterModelPart, ModelPart &rSlaveModelPart,
                                Parameters Settings) : Process(Flags()), mrMasterModelPart(rMasterModelPart),
                                mrSlaveModelPart(rSlaveModelPart), mParameters(Settings)
{
    // Initializing
    const Parameters default_parameters = this->GetDefaultParameters();
    mParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

    mSearchMaxResults = mParameters["search_settings"]["max_results"].GetInt();
    mSearchTolerance = mParameters["search_settings"]["tolerance"].GetDouble();
}

SlidingInterfaceProcess::~SlidingInterfaceProcess()
{
}


/**
    * @brief Function initializes the process
    */
void SlidingInterfaceProcess::ExecuteInitialize()
{
}

/**
    * @brief Function initializes the solution step
    */
void SlidingInterfaceProcess::ExecuteInitializeSolutionStep()
{

    KRATOS_TRY;
    mrMasterModelPart.RemoveMasterSlaveConstraintsFromAllLevels(TO_ERASE);
    VariableUtils().SetHistoricalVariableToZero<double>(MESH_VELOCITY_X, mrMasterModelPart.Nodes());
    VariableUtils().SetHistoricalVariableToZero<double>(MESH_VELOCITY_Y, mrMasterModelPart.Nodes());
    VariableUtils().SetHistoricalVariableToZero<double>(MESH_VELOCITY_Z, mrMasterModelPart.Nodes());
    const int domain_size = mrMasterModelPart.GetProcessInfo()[DOMAIN_SIZE];
    // Rotate the master so it goes to the slave
    if (domain_size == 2){
        ApplyConstraintsForSlidingInterface<2>();
    }
    else if (domain_size == 3){
        ApplyConstraintsForSlidingInterface<3>();
    } else {
        KRATOS_ERROR <<"Periodic conditions are designed only for 2 and 3 Dimensional cases ! "<<std::endl;
    }

    KRATOS_CATCH("");
}


void SlidingInterfaceProcess::ExecuteFinalizeSolutionStep()
{
    mrMasterModelPart.RemoveMasterSlaveConstraintsFromAllLevels(TO_ERASE);
    Model& current_model = mrMasterModelPart.GetModel();
    const DataCommunicator &r_comm =
        mrMasterModelPart.GetCommunicator().GetDataCommunicator();
    if (r_comm.IsDistributed())
        current_model.DeleteModelPart("gathered_master");
}



const Parameters SlidingInterfaceProcess::GetDefaultParameters() const
{
    const Parameters default_parameters(R"(
    {
        "variable_names":[],
        "search_settings":{
            "max_results":100000,
            "tolerance": 1E-6
        }
    }  )");
    return default_parameters;
}


/**
    * @brief Function to print the information about this current process
    */
void SlidingInterfaceProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream <<"SlidingInterfaceProcess Process "<<std::endl;
}


ModelPart& SlidingInterfaceProcess::GetSearchModelpart()
{

    const DataCommunicator &r_comm =
        mrMasterModelPart.GetCommunicator().GetDataCommunicator();
    Model& current_model = mrMasterModelPart.GetModel();
    ModelPart &gathered_master = r_comm.IsDistributed() ? current_model.CreateModelPart("gathered_master") : mrMasterModelPart;

    #ifdef KRATOS_USING_MPI
        if (r_comm.IsDistributed())
            GatherModelPartOnAllRanksUtility::GatherModelPartOnAllRanks(mrMasterModelPart, gathered_master);
    #endif

    return gathered_master;
}

template <int TDim>
void SlidingInterfaceProcess::ApplyConstraintsForSlidingInterface()
{
    const double start_apply = OpenMPUtils::GetCurrentTime();
    const int num_vars = mParameters["variable_names"].size();

    ModelPart& r_search_modelpart = GetSearchModelpart();
    BinBasedFastPointLocatorConditions<TDim> bin_based_point_locator(r_search_modelpart);
    bin_based_point_locator.UpdateSearchDatabase();

    const int num_slave_nodes = mrSlaveModelPart.GetCommunicator().LocalMesh().NumberOfNodes();
    const NodeIteratorType it_slave_node_begin = mrSlaveModelPart.GetCommunicator().LocalMesh().NodesBegin();
    IndexType num_slaves_found = 0;

    auto& r_comm = mrMasterModelPart.GetCommunicator();
    const int current_rank = r_comm.MyPID();


    #pragma omp parallel for schedule(guided, 512) reduction( + : num_slaves_found )
    for(int i_node = 0; i_node<num_slave_nodes; ++i_node)
    {
        Condition::Pointer p_host_cond;
        VectorType shape_function_values;
        NodeIteratorType it_slave_node = it_slave_node_begin+i_node;
        const auto& slave_node_coords = it_slave_node->Coordinates();
        // if(r_comm.IsDistributed())
        //     if(it_slave_node->GetSolutionStepValue(PARTITION_INDEX) != current_rank)
        //         continue;

        // Finding the host element for this node
        const bool is_found = bin_based_point_locator.FindPointOnMeshSimplified(slave_node_coords, shape_function_values, p_host_cond, mSearchMaxResults, mSearchTolerance);
        if(is_found)
        {
            ++num_slaves_found;
            for (int j = 0; j < num_vars; j++)
            {
                const std::string var_name = mParameters["variable_names"][j].GetString();
                ConstraintSlaveNodeWithConditionForVariable<TDim>(*it_slave_node, p_host_cond->GetGeometry() , shape_function_values, var_name);
            }
        }
    }
    KRATOS_WARNING_IF("SlidingInterfaceProcess",num_slaves_found != num_slave_nodes)<<"Sliding interface condition cannot be applied for all the nodes.  "<<  num_slaves_found <<std::endl;
    const double end_apply = OpenMPUtils::GetCurrentTime();
    KRATOS_INFO("SlidingInterfaceProcess")<<"Applying sliding interface took : "<<end_apply - start_apply<<" seconds." <<std::endl;
}

template <int TDim>
void SlidingInterfaceProcess::ConstraintSlaveNodeWithConditionForVariable(NodeType& rSlaveNode, const GeometryType& rHostedGeometry, const VectorType& rWeights,const std::string& rVarName )
{
    const VariableType& r_var = KratosComponents<VariableType>::Get(rVarName);
    auto slave_variable_val = rSlaveNode.GetSolutionStepValue(r_var);
    auto n_slave_variable_val = rSlaveNode.GetSolutionStepValue(r_var,1);
    n_slave_variable_val = 0.0;
    slave_variable_val = 0.0;

    // Reference constraint
    const auto& r_clone_constraint = KratosComponents<MasterSlaveConstraint>::Get("LinearMasterSlaveConstraint");

    IndexType master_index = 0;
    for (auto& master_node : rHostedGeometry)
    {
        const double master_weight = rWeights(master_index);
        const auto master_variable_val = master_node.GetSolutionStepValue(r_var);
        const auto n_master_variable_val = master_node.GetSolutionStepValue(r_var,1);
        #pragma omp critical
        {
            slave_variable_val += master_weight*master_variable_val;
            n_slave_variable_val += master_weight*n_master_variable_val;
            int current_constraint_id = mrMasterModelPart.GetRootModelPart().NumberOfMasterSlaveConstraints();
            auto constraint = r_clone_constraint.Create(++current_constraint_id,master_node, r_var, rSlaveNode, r_var, master_weight, 0.0);
            constraint->Set(TO_ERASE);
            mrMasterModelPart.AddMasterSlaveConstraint(constraint);
        }
        master_index++;
    }
}



}