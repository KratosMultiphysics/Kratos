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
## including kratos path                                     #####
## ATTENTION: the following lines have to be adapted to      #####
##            match your acrtual configuration               #####
##################################################################
## kratos test examples:
kratos_root_path = '../../../..'
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
sys.path.append('./balken.gid')
import balken_include
from balken_include import *
# calculate insitu-stress for geology_virgin.gid
model = balken_include.Model('balken','./')
model.InitializeModel()

##################################################################
###  SIMULATION  #################################################
##################################################################
time = 0.0
model.SetCalculateInSituStress( False )
model.Solve( time, 0, 0, 0, 0 )
model.WriteOutput( time )
print "Calculation done"
sys.exit(0)
