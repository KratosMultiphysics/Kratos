import KratosMultiphysics


def CheckVariableStatistics(model_part, input_variable):
    data_comm = KratosMultiphysics.DataCommunicator.GetDefault()
    if data_comm.IsDistributed():
        mean_var = 0.0
        fixed_count = 0
        for node in model_part.GetCommunicator().LocalMesh().Nodes:
            fixed_count += int(node.IsFixed(input_variable))
            mean_var += node.GetSolutionStepValue(input_variable)

        global_mean_var = data_comm.SumAll(mean_var)
        global_nodes = data_comm.SumAll(model_part.GetCommunicator().LocalMesh().NumberOfNodes())
        global_fixed_count = data_comm.SumAll(fixed_count)

        print("Mean", input_variable.Name(), global_mean_var, global_nodes, global_mean_var/global_nodes, global_fixed_count)
    else:
        mean_var = 0.0
        fixed_count = 0
        for node in model_part.Nodes:
            fixed_count += int(node.IsFixed(input_variable))
            mean_var += node.GetSolutionStepValue(input_variable)

        print("Mean", input_variable.Name(), mean_var, model_part.NumberOfNodes(), mean_var/model_part.NumberOfNodes(), fixed_count)% 