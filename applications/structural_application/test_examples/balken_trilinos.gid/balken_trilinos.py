##################################################################
##### ekate - Enhanced KRATOS for Advanced Tunnel Enineering #####
##### copyright by CIMNE, Barcelona, Spain                   #####
#####          and Janosch Stascheit for TUNCONSTRUCT        #####
##### all rights reserved                                    #####
##################################################################
#setting the domain size for the problem to be solved
domain_size = 3

#importing MPI ... for this boost 1.35 or superior is needed
import mpi
print "i am ",mpi.rank , " of ",mpi.size


##################################################################
##################################################################
## ATTENTION: here the order is important                    #####
##################################################################
##################################################################
sys.path.append('./balken_trilinos.gid')
import balken_trilinos_include
from balken_trilinos_include import *
# calculate insitu-stress for geology_virgin.gid
model = balken_trilinos_include.Model('balken_trilinos','./')
model.InitializeModel()

##################################################################
###  SIMULATION  #################################################
##################################################################
time = 0.0
model.SetCalculateInSituStress( False )
print( mpi.rank )
print( model.model_part )
model.Solve( time, 0, 0, 0, 0 )
print( model.model_part )
model.WriteOutput( time )
print "Calculation done"
sys.exit(0)
