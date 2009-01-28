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
kratos_libs_path = '/home/hurga/kratos_local/kratos/libs/' ##kratos_root/libs
#kratos_libs_path = '/home/hurga/kratosR1/libs/' ##kratos_root/libs
kratos_applications_path = '/home/hurga/kratos_local/kratos/applications/' ##kratos_root/applications
#kratos_applications_path = '/home/hurga/kratosR1/applications/' ##kratos_root/applications
##################################################################
##################################################################
import sys
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)

#importing Kratos main library
#from Kratos import *
#kernel = Kernel()   #defining kernel

#importing applications
#import applications_interface
#applications_interface.Import_StructuralApplication = True
#applications_interface.Import_EkateAuxiliaryApplication = True
#applications_interface.ImportApplications(kernel, kratos_applications_path)
#from KratosStructuralApplication import *
#from KratosExternalSolversApplication import *
#from KratosEkateAuxiliaryApplication import *
##################################################################
##################################################################
sys.path.append('./activation_test2.gid')
import activation_test2_include
from activation_test2_include import *
# calculate insitu-stress for geology_virgin.gid
model = activation_test2_include.Model('activation_test2','./')
model.InitializeModel()

##################################################################
###  SIMULATION  #################################################
##################################################################
# =====================
# | USER SCRIPT FOR CALCULATION OF AUSBLAS.GID |
# vvvvvvvvvvvvvvvvvvvvv
# ekate simulation script for activation_test.gid
#setting up activation flags
from_excavation = 0
to_excavation = 1
from_lining = -12
to_lining = -11
Dt =            1
time = 0.0
for excavation_step in range(1,10):
	time = time+Dt
	to_excavation = to_excavation+1
	to_lining = to_lining+1
	model.Solve(time,from_excavation,to_excavation,from_lining,to_lining)
	model.WriteOutput( time )
	print("~~~~~~~~~~~~~~STEP DONE~~~~~~~~~~~~~~")
# END OF USER SCRIPT
