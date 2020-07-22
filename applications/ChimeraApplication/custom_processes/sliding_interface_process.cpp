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
    KRATOS_TRY;
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

/**
    * @brief Function initializes the solution step
    */
void SlidingInterfaceProcess::ExecuteInitializeSolutionStep()
{
}


const Parameters SlidingInterfaceProcess::GetDefaultParameters() const
{
    const Parameters default_parameters(R"(
    {
        "variable_names":["VELOCITY_X", "VELOCITY_Y", "VELOCITY_Z","PRESSURE"],
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


template <int TDim>
void SlidingInterfaceProcess::ApplyConstraintsForSlidingInterface()
{
    const double start_apply = OpenMPUtils::GetCurrentTime();
    const int num_vars = mParameters["variable_names"].size();
    BinBasedFastPointLocatorConditions<TDim> bin_based_point_locator(mrMasterModelPart);
    bin_based_point_locator.UpdateSearchDatabase();

    const int num_slave_nodes = mrSlaveModelPart.NumberOfNodes();
    const NodeIteratorType it_slave_node_begin = mrSlaveModelPart.NodesBegin();

    IndexType num_slaves_found = 0;

    #pragma omp parallel for schedule(guided, 512) reduction( + : num_slaves_found )
    for(int i_node = 0; i_node<num_slave_nodes; ++i_node)
    {
        Condition::Pointer p_host_cond;
        VectorType shape_function_values;
        NodeIteratorType it_slave_node = it_slave_node_begin+i_node;
        const auto& slave_node_coords = it_slave_node->Coordinates();

        // Finding the host element for this node
        const bool is_found = bin_based_point_locator.FindPointOnMeshSimplified(slave_node_coords, shape_function_values, p_host_cond, mSearchMaxResults, mSearchTolerance);
        if(is_found)
        {
            ++num_slaves_found;
            for (int j = 0; j < num_vars; j++)
            {
                const std::string var_name = mParameters["variable_names"][j].GetString();
                // Checking if the variable is a vector variable
                ConstraintSlaveNodeWithConditionForVariable<TDim>(*it_slave_node, p_host_cond->GetGeometry() , shape_function_values, var_name);
            }
        }
    }
    KRATOS_WARNING_IF("SlidingInterfaceProcess",num_slaves_found != mrSlaveModelPart.NumberOfNodes())<<"Sliding interface condition cannot be applied for all the nodes."<<std::endl;
    const double end_apply = OpenMPUtils::GetCurrentTime();
    KRATOS_INFO("SlidingInterfaceProcess")<<"Applying sliding interface took : "<<end_apply - start_apply<<" seconds." <<std::endl;
}

template <int TDim>
void SlidingInterfaceProcess::ConstraintSlaveNodeWithConditionForVariable(NodeType& rSlaveNode, const GeometryType& rHostedGeometry, const VectorType& rWeights,const std::string& rVarName )
{
    const VariableType& r_var = KratosComponents<VariableType>::Get(rVarName);

    // Reference constraint
    const auto& r_clone_constraint = KratosComponents<MasterSlaveConstraint>::Get("LinearMasterSlaveConstraint");

    IndexType master_index = 0;
    for (auto& master_node : rHostedGeometry)
    {
        const double master_weight = rWeights(master_index);
        #pragma omp critical
        {
            int current_constraint_id = ( mrMasterModelPart.GetRootModelPart().MasterSlaveConstraints().end()-1 )->Id();
            auto constraint = r_clone_constraint.Create(++current_constraint_id,master_node, r_var, rSlaveNode, r_var, master_weight, 0.0);
            mrMasterModelPart.AddMasterSlaveConstraint(constraint);
        }
        master_index++;
    }
}



}