##################################################################
##### KRATOS Simulation script                               #####
##### copyright by CIMNE, Barcelona, Spain                   #####
#####          and Janosch Stascheit for TUNCONSTRUCT        #####
#####          and Giang H. Bui for SD-RUB                   #####
##### all rights reserved                                    #####
##################################################################
##################################################################
## ATTENTION: here the order is important                    #####
##################################################################
## including kratos path                                     #####
## ATTENTION: the following lines have to be adapted to      #####
##            match your acrtual configuration               #####
##################################################################
import sys
import os
kratos_root_path=os.environ['KRATOS_ROOT_PATH']
##setting up paths
kratos_libs_path = kratos_root_path+'libs' ##kratos_root/libs
kratos_applications_path = kratos_root_path+'applications' ##kratos_root/applications
##################################################################
##################################################################
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)

##################################################################
##################################################################
sys.path.append('./rEpLaCeMeNtStRiNg.gid')
*if(strcmp(GenData(Parallel_Execution),"shared")==0)
parallel_type = 'shared'
*elseif(strcmp(GenData(Parallel_Execution),"distributed")==0)
parallel_type = 'distributed'
*endif
if(parallel_type == 'shared'):
    import rEpLaCeMeNtStRiNg_shared_include
    from rEpLaCeMeNtStRiNg_shared_include import **
    model = rEpLaCeMeNtStRiNg_shared_include.Model('rEpLaCeMeNtStRiNg',os.getcwd()+"/")
elif(parallel_type == 'distributed'):
    import rEpLaCeMeNtStRiNg_distributed_include
    from rEpLaCeMeNtStRiNg_distributed_include import **
    model = rEpLaCeMeNtStRiNg_distributed_include.Model('rEpLaCeMeNtStRiNg',os.getcwd()+"/")
model.InitializeModel()

##################################################################
###  SIMULATION  #################################################
##################################################################

