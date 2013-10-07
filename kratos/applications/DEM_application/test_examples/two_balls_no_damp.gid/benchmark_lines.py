#including kratos path                                                                        ### BENCHMARK ###
kratos_path = '../../../..'                                                                   ### BENCHMARK ###
kratos_libs_path = '../../../../libs' ##kratos_root/libs                                      ### BENCHMARK ###
kratos_applications_path = '../../../../applications' ##kratos_root/applications              ### BENCHMARK ###
kratos_benchmarking_path = '../../../../benchmarking'                                         ### BENCHMARK ###
import sys                                                                                    ### BENCHMARK ###
sys.path.append(kratos_path)                                                                  ### BENCHMARK ###
sys.path.append(kratos_libs_path)                                                             ### BENCHMARK ###
sys.path.append(kratos_applications_path)                                                     ### BENCHMARK ###
sys.path.append(kratos_benchmarking_path)                                                     ### BENCHMARK ###
from DEM_explicit_solver_var import *                                                         ### BENCHMARK ###
import benchmarking                                                                           ### BENCHMARK ###
def FindNode(node_list,X,Y,Z):                                                                ### BENCHMARK ###
    for node in node_list:                                                                    ### BENCHMARK ###
	if((node.X-X)**2 + (node.Y-Y)**2 + (node.Z-Z)**2 < .000001):                          ### BENCHMARK ###
		return node                                                                   ### BENCHMARK ###
                                                                                              ### BENCHMARK ###
def BenchmarkCheck(time, node1):                                                              ### BENCHMARK ###
    benchmarking.Output(time, "Time")                                                         ### BENCHMARK ###
    benchmarking.Output(node1.GetSolutionStepValue(DISPLACEMENT_Y), "Node Displacement", 0.01)### BENCHMARK ###
    benchmarking.Output(node1.GetSolutionStepValue(VELOCITY_Y), "Node Velocity", 0.01)        ### BENCHMARK ###
###############################################################                               ### BENCHMARK ###
node1 = FindNode(balls_model_part.Nodes , 0.0, 1.0, 0.0)                                      ### BENCHMARK ###
print node1 #there is a memory problem with the string                                        ### BENCHMARK ###
###############################################################                               ### BENCHMARK ###
        BenchmarkCheck(time, node1)                                                           ### BENCHMARK ###
