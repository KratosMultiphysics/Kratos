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
## including kratos path                                     #####
## ATTENTION: the following lines have to be adapted to      #####
##            match your acrtual configuration               #####
##################################################################
## calypso:
#kratos_root_path = '/home/hurga/kratos_bcn/kratos'
#kratos_root_path = '/home/hurga/kratos_merged/kratos'
## alderaan:
#kratos_root_path = '/localsw/kratos/kratos'
## r2d2
#kratos_root_path = '/home/stasch/kratos_bcn/kratos'
## astra
kratos_root_path = '../../../../'
##setting up paths
kratos_libs_path = kratos_root_path+'/libs' ##kratos_root/libs
kratos_applications_path = kratos_root_path+'/applications' ##kratos_root/applications
##################################################################
##################################################################
import sys
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)

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
