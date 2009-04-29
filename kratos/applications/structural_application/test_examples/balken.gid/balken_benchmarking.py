##################################################################
##### ekate - Enhanced KRATOS for Advanced Tunnel Enineering #####
##### copyright by CIMNE, Barcelona, Spain                   #####
#####          and Janosch Stascheit for TUNCONSTRUCT        #####
##### all rights reserved                                    #####
##################################################################
#setting the domain size for the problem to be solved
domain_size = 3
##################################################################
##################################################################
## ATTENTION: here the order is important                    #####
##################################################################

#including kratos path
kratos_libs_path = '../../../../libs' ##kratos_root/libs
kratos_applications_path = '../../../../applications' ##kratos_root/applications
kratos_benchmarking_path = '../../../../benchmarking' ##kratos_root/benchmarking
##################################################################
##################################################################
import sys
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)
sys.path.append(kratos_benchmarking_path)
##################################################################
##################################################################
import benchmarking
##################################################################
##################################################################

def BenchmarkCheck(time, node1, node2, node3, node4):
    benchmarking.Output(node1.GetSolutionStepValue(DISPLACEMENT_Z), "Node 1 Displacement_z", 0.00001 )
    benchmarking.Output(node2.GetSolutionStepValue(DISPLACEMENT_Z), "Node 2 Displacement_z", 0.00001 )
    benchmarking.Output(node3.GetSolutionStepValue(DISPLACEMENT_Z), "Node 3 Displacement_z", 0.00001 )
    benchmarking.Output(node4.GetSolutionStepValue(DISPLACEMENT_Z), "Node 4 Displacement_z", 0.00001 )

def AnalyticalResults(time, node1, node2,node3, node4):
    benchmarking.Output(node1.GetSolutionStepValue(DISPLACEMENT_Z), "Node 1 Displacement_z", 0.00001 )
    benchmarking.Output(node2.GetSolutionStepValue(DISPLACEMENT_Z), "Node 2 Displacement_z", 0.00001 )
    benchmarking.Output(node3.GetSolutionStepValue(DISPLACEMENT_Z), "Node 3 Displacement_z", 0.00001 )
    benchmarking.Output(node4.GetSolutionStepValue(DISPLACEMENT_Z), "Node 4 Displacement_z", 0.00001 )


sys.path.append('./balken.gid')
import balken_include
from balken_include import *
# calculate insitu-stress for geology_virgin.gid
model = balken_include.Model('balken','./')
model.InitializeModel()

##################################################################
###  SIMULATION  #################################################
##################################################################

node_1 = model.model_part.Nodes[246]
node_2 = model.model_part.Nodes[247]
node_3 = model.model_part.Nodes[76]
node_4 = model.model_part.Nodes[301]


time = 0.0
model.SetCalculateInSituStress( False )
model.Solve( time, 0, 0, 0, 0 )
model.WriteOutput( time )

if (benchmarking.InBuildReferenceMode()):
	AnalyticalResults(time, node_1, node_2, node_3, node_4)
else:
	BenchmarkCheck(time, node_1, node_2, node_3, node_4)

print "Calculation done"
sys.exit(0)
