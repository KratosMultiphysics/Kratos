# Import Python libraries
import time
import pickle
try:
    from threadpoolctl import *
except:
    pass

# Import Kratos, XMC, distributed environment
import KratosMultiphysics
from KratosMultiphysics.MultilevelMonteCarloApplication.adaptive_refinement_utilities import AdaptiveRefinement
from xmc.distributedEnvironmentFramework import *


####################################################################################################
############################################ WRAPPERS ##############################################
####################################################################################################


def executeInstanceStochasticAdaptiveRefinementAllAtOnce_Wrapper(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_analysis,time_for_qoi,mapping_flag,adaptive_refinement_jump_to_finest_level,print_to_file,current_contribution):
    if (current_index == 0):
        qoi,time_for_qoi = ExecuteInstanceStochasticAdaptiveRefinementAllAtOnceAuxLev0_Task(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_analysis,time_for_qoi,mapping_flag,adaptive_refinement_jump_to_finest_level,print_to_file,"filename_level_"+str(current_index)+"_contribution_"+str(current_contribution)+"_random_variable_"+str(random_variable[0])+".dat")
    elif (current_index == 1):
        qoi,time_for_qoi = ExecuteInstanceStochasticAdaptiveRefinementAllAtOnceAuxLev1_Task(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_analysis,time_for_qoi,mapping_flag,adaptive_refinement_jump_to_finest_level,print_to_file,"filename_level_"+str(current_index)+"_contribution_"+str(current_contribution)+"_random_variable_"+str(random_variable[0])+".dat")
    elif (current_index == 2):
        qoi,time_for_qoi = ExecuteInstanceStochasticAdaptiveRefinementAllAtOnceAuxLev2_Task(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_analysis,time_for_qoi,mapping_flag,adaptive_refinement_jump_to_finest_level,print_to_file,"filename_level_"+str(current_index)+"_contribution_"+str(current_contribution)+"_random_variable_"+str(random_variable[0])+".dat")
    elif (current_index == 3):
        qoi,time_for_qoi = ExecuteInstanceStochasticAdaptiveRefinementAllAtOnceAuxLev3_Task(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_analysis,time_for_qoi,mapping_flag,adaptive_refinement_jump_to_finest_level,print_to_file,"filename_level_"+str(current_index)+"_contribution_"+str(current_contribution)+"_random_variable_"+str(random_variable[0])+".dat")
    elif (current_index == 4):
        qoi,time_for_qoi = ExecuteInstanceStochasticAdaptiveRefinementAllAtOnceAuxLev4_Task(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_analysis,time_for_qoi,mapping_flag,adaptive_refinement_jump_to_finest_level,print_to_file,"filename_level_"+str(current_index)+"_contribution_"+str(current_contribution)+"_random_variable_"+str(random_variable[0])+".dat")
    elif (current_index == 5):
        qoi,time_for_qoi = ExecuteInstanceStochasticAdaptiveRefinementAllAtOnceAuxLev5_Task(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_analysis,time_for_qoi,mapping_flag,adaptive_refinement_jump_to_finest_level,print_to_file,"filename_level_"+str(current_index)+"_contribution_"+str(current_contribution)+"_random_variable_"+str(random_variable[0])+".dat")
    else:
        raise Exception("Level not supported")
    return qoi,time_for_qoi

def executeInstanceStochasticAdaptiveRefinementMultipleTasks_Wrapper(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_local_index,current_analysis,time_for_qoi,mapping_flag,print_to_file,current_contribution,pickled_mapping_reference_model=None):
    if (current_index == 0):
        qoi,pickled_current_model,time_for_qoi = ExecuteInstanceStochasticAdaptiveRefinementMultipleTasksAuxLev0_Task(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_local_index,current_analysis,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,"filename_level_"+str(current_index)+"_contribution_"+str(current_contribution)+"_random_variable_"+str(random_variable[0])+".dat")
    elif (current_index == 1):
        qoi,pickled_current_model,time_for_qoi = ExecuteInstanceStochasticAdaptiveRefinementMultipleTasksAuxLev1_Task(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_local_index,current_analysis,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,"filename_level_"+str(current_index)+"_contribution_"+str(current_contribution)+"_random_variable_"+str(random_variable[0])+".dat")
    elif (current_index == 2):
        qoi,pickled_current_model,time_for_qoi = ExecuteInstanceStochasticAdaptiveRefinementMultipleTasksAuxLev2_Task(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_local_index,current_analysis,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,"filename_level_"+str(current_index)+"_contribution_"+str(current_contribution)+"_random_variable_"+str(random_variable[0])+".dat")
    elif (current_index == 3):
        qoi,pickled_current_model,time_for_qoi = ExecuteInstanceStochasticAdaptiveRefinementMultipleTasksAuxLev3_Task(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_local_index,current_analysis,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,"filename_level_"+str(current_index)+"_contribution_"+str(current_contribution)+"_random_variable_"+str(random_variable[0])+".dat")
    elif (current_index == 4):
        qoi,pickled_current_model,time_for_qoi = ExecuteInstanceStochasticAdaptiveRefinementMultipleTasksAuxLev4_Task(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_local_index,current_analysis,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,"filename_level_"+str(current_index)+"_contribution_"+str(current_contribution)+"_random_variable_"+str(random_variable[0])+".dat")
    elif (current_index == 5):
        qoi,pickled_current_model,time_for_qoi = ExecuteInstanceStochasticAdaptiveRefinementMultipleTasksAuxLev5_Task(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_local_index,current_analysis,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,"filename_level_"+str(current_index)+"_contribution_"+str(current_contribution)+"_random_variable_"+str(random_variable[0])+".dat")
    else:
        raise Exception("Level not supported")
    return qoi,pickled_current_model,time_for_qoi

def executeInstanceDeterministicAdaptiveRefinement_Wrapper(current_index,pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,current_contribution):
    if (current_index == 0):
        qoi,time_for_qoi = executeInstanceDeterministicAdaptiveRefinementAuxLev0_Task(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,"filename_level_"+str(current_index)+"_contribution_"+str(current_contribution)+"_random_variable_"+str(random_variable[0])+".dat")
    elif (current_index == 1):
        qoi,time_for_qoi = executeInstanceDeterministicAdaptiveRefinementAuxLev1_Task(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,"filename_level_"+str(current_index)+"_contribution_"+str(current_contribution)+"_random_variable_"+str(random_variable[0])+".dat")
    elif (current_index == 2):
        qoi,time_for_qoi = executeInstanceDeterministicAdaptiveRefinementAuxLev2_Task(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,"filename_level_"+str(current_index)+"_contribution_"+str(current_contribution)+"_random_variable_"+str(random_variable[0])+".dat")
    elif (current_index == 3):
        qoi,time_for_qoi = executeInstanceDeterministicAdaptiveRefinementAuxLev3_Task(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,"filename_level_"+str(current_index)+"_contribution_"+str(current_contribution)+"_random_variable_"+str(random_variable[0])+".dat")
    elif (current_index == 4):
        qoi,time_for_qoi = executeInstanceDeterministicAdaptiveRefinementAuxLev4_Task(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,"filename_level_"+str(current_index)+"_contribution_"+str(current_contribution)+"_random_variable_"+str(random_variable[0])+".dat")
    elif (current_index == 5):
        qoi,time_for_qoi = executeInstanceDeterministicAdaptiveRefinementAuxLev5_Task(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,"filename_level_"+str(current_index)+"_contribution_"+str(current_contribution)+"_random_variable_"+str(random_variable[0])+".dat")
    else:
        raise Exception("Level not supported")
    return qoi,time_for_qoi

def executeInstanceReadingFromFile_Wrapper(current_index,pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,current_contribution):
    if (current_index == 0):
        qoi,time_for_qoi = executeInstanceReadingFromFileAuxLev0_Task(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,"filename_level_"+str(current_index)+"_contribution_"+str(current_contribution)+"_random_variable_"+str(random_variable[0])+".dat")
    elif (current_index == 1):
        qoi,time_for_qoi = executeInstanceReadingFromFileAuxLev1_Task(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,"filename_level_"+str(current_index)+"_contribution_"+str(current_contribution)+"_random_variable_"+str(random_variable[0])+".dat")
    elif (current_index == 2):
        qoi,time_for_qoi = executeInstanceReadingFromFileAuxLev2_Task(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,"filename_level_"+str(current_index)+"_contribution_"+str(current_contribution)+"_random_variable_"+str(random_variable[0])+".dat")
    elif (current_index == 3):
        qoi,time_for_qoi = executeInstanceReadingFromFileAuxLev3_Task(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,"filename_level_"+str(current_index)+"_contribution_"+str(current_contribution)+"_random_variable_"+str(random_variable[0])+".dat")
    elif (current_index == 4):
        qoi,time_for_qoi = executeInstanceReadingFromFileAuxLev4_Task(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,"filename_level_"+str(current_index)+"_contribution_"+str(current_contribution)+"_random_variable_"+str(random_variable[0])+".dat")
    elif (current_index == 5):
        qoi,time_for_qoi = executeInstanceReadingFromFileAuxLev5_Task(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,"filename_level_"+str(current_index)+"_contribution_"+str(current_contribution)+"_random_variable_"+str(random_variable[0])+".dat")
    else:
        raise Exception("Level not supported")
    return qoi,time_for_qoi


####################################################################################################
############################################## TASKS ###############################################
####################################################################################################


############################### StochasticAdaptiveRefinementAllAtOnce ##############################

# @ExaquteTask(filename=FILE_OUT,returns=2)
@constraint(ComputingUnits="${computing_units_mlmc_execute_0}")
@ExaquteTask(returns=2)
def ExecuteInstanceStochasticAdaptiveRefinementAllAtOnceAuxLev0_Task(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_analysis,time_for_qoi,mapping_flag,adaptive_refinement_jump_to_finest_level,print_to_file,filename):
    try:
        open_mp_threads = int(os.environ["computing_units_mlmc_execute_0"])
        threadpool_limits(limits=open_mp_threads)
    except:
        open_mp_threads = 1
    pickled_coarsest_model = pickled_coarse_model
    for current_local_index in range(current_index+1):
        if ((adaptive_refinement_jump_to_finest_level is False) or (adaptive_refinement_jump_to_finest_level is True and (current_local_index == 0 or current_local_index == current_index))):
            qoi,pickled_current_model,time_for_qoi = \
                ExecuteInstanceStochasticAdaptiveRefinementAux_Functionality(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_local_index,current_analysis,time_for_qoi,open_mp_threads,mapping_flag,pickled_coarsest_model,print_to_file,filename)
            del(pickled_coarse_model)
            pickled_coarse_model = pickled_current_model
            del(pickled_current_model)
    return qoi,time_for_qoi

# @ExaquteTask(filename=FILE_OUT,returns=2)
@constraint(ComputingUnits="${computing_units_mlmc_execute_1}")
@ExaquteTask(returns=2)
def ExecuteInstanceStochasticAdaptiveRefinementAllAtOnceAuxLev1_Task(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_analysis,time_for_qoi,mapping_flag,adaptive_refinement_jump_to_finest_level,print_to_file,filename):
    try:
        open_mp_threads = int(os.environ["computing_units_mlmc_execute_1"])
        threadpool_limits(limits=open_mp_threads)
    except:
        open_mp_threads = 1
    pickled_coarsest_model = pickled_coarse_model
    for current_local_index in range(current_index+1):
        if ((adaptive_refinement_jump_to_finest_level is False) or (adaptive_refinement_jump_to_finest_level is True and (current_local_index == 0 or current_local_index == current_index))):
            qoi,pickled_current_model,time_for_qoi = \
                ExecuteInstanceStochasticAdaptiveRefinementAux_Functionality(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_local_index,current_analysis,time_for_qoi,open_mp_threads,mapping_flag,pickled_coarsest_model,print_to_file,filename)
            del(pickled_coarse_model)
            pickled_coarse_model = pickled_current_model
            del(pickled_current_model)
    return qoi,time_for_qoi

# @ExaquteTask(filename=FILE_OUT,returns=2)
@constraint(ComputingUnits="${computing_units_mlmc_execute_2}")
@ExaquteTask(returns=2)
def ExecuteInstanceStochasticAdaptiveRefinementAllAtOnceAuxLev2_Task(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_analysis,time_for_qoi,mapping_flag,adaptive_refinement_jump_to_finest_level,print_to_file,filename):
    try:
        open_mp_threads = int(os.environ["computing_units_mlmc_execute_2"])
        threadpool_limits(limits=open_mp_threads)
    except:
        open_mp_threads = 1
    pickled_coarsest_model = pickled_coarse_model
    for current_local_index in range(current_index+1):
        if ((adaptive_refinement_jump_to_finest_level is False) or (adaptive_refinement_jump_to_finest_level is True and (current_local_index == 0 or current_local_index == current_index))):
            qoi,pickled_current_model,time_for_qoi = \
                ExecuteInstanceStochasticAdaptiveRefinementAux_Functionality(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_local_index,current_analysis,time_for_qoi,open_mp_threads,mapping_flag,pickled_coarsest_model,print_to_file,filename)
            del(pickled_coarse_model)
            pickled_coarse_model = pickled_current_model
            del(pickled_current_model)
        else: # not running since we jump from coarsest to finest level
            pass
    return qoi,time_for_qoi

# @ExaquteTask(filename=FILE_OUT,returns=2)
@constraint(ComputingUnits="${computing_units_mlmc_execute_3}")
@ExaquteTask(returns=2)
def ExecuteInstanceStochasticAdaptiveRefinementAllAtOnceAuxLev3_Task(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_analysis,time_for_qoi,mapping_flag,adaptive_refinement_jump_to_finest_level,print_to_file,filename):
    try:
        open_mp_threads = int(os.environ["computing_units_mlmc_execute_3"])
        threadpool_limits(limits=open_mp_threads)
    except:
        open_mp_threads = 1
    pickled_coarsest_model = pickled_coarse_model
    for current_local_index in range(current_index+1):
        if ((adaptive_refinement_jump_to_finest_level is False) or (adaptive_refinement_jump_to_finest_level is True and (current_local_index == 0 or current_local_index == current_index))):
            qoi,pickled_current_model,time_for_qoi = \
                ExecuteInstanceStochasticAdaptiveRefinementAux_Functionality(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_local_index,current_analysis,time_for_qoi,open_mp_threads,mapping_flag,pickled_coarsest_model,print_to_file,filename)
            del(pickled_coarse_model)
            pickled_coarse_model = pickled_current_model
            del(pickled_current_model)
    return qoi,time_for_qoi

# @ExaquteTask(filename=FILE_OUT,returns=2)
@constraint(ComputingUnits="${computing_units_mlmc_execute_4}")
@ExaquteTask(returns=2)
def ExecuteInstanceStochasticAdaptiveRefinementAllAtOnceAuxLev4_Task(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_analysis,time_for_qoi,mapping_flag,adaptive_refinement_jump_to_finest_level,print_to_file,filename):
    try:
        open_mp_threads = int(os.environ["computing_units_mlmc_execute_4"])
        threadpool_limits(limits=open_mp_threads)
    except:
        open_mp_threads = 1
    pickled_coarsest_model = pickled_coarse_model
    for current_local_index in range(current_index+1):
        if ((adaptive_refinement_jump_to_finest_level is False) or (adaptive_refinement_jump_to_finest_level is True and (current_local_index == 0 or current_local_index == current_index))):
            qoi,pickled_current_model,time_for_qoi = \
                ExecuteInstanceStochasticAdaptiveRefinementAux_Functionality(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_local_index,current_analysis,time_for_qoi,open_mp_threads,mapping_flag,pickled_coarsest_model,print_to_file,filename)
            del(pickled_coarse_model)
            pickled_coarse_model = pickled_current_model
            del(pickled_current_model)
    return qoi,time_for_qoi

# @ExaquteTask(filename=FILE_OUT,returns=2)
@constraint(ComputingUnits="${computing_units_mlmc_execute_5}")
@ExaquteTask(returns=2)
def ExecuteInstanceStochasticAdaptiveRefinementAllAtOnceAuxLev5_Task(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_analysis,time_for_qoi,mapping_flag,adaptive_refinement_jump_to_finest_level,print_to_file,filename):
    try:
        open_mp_threads = int(os.environ["computing_units_mlmc_execute_5"])
        threadpool_limits(limits=open_mp_threads)
    except:
        open_mp_threads = 1
    pickled_coarsest_model = pickled_coarse_model
    for current_local_index in range(current_index+1):
        if ((adaptive_refinement_jump_to_finest_level is False) or (adaptive_refinement_jump_to_finest_level is True and (current_local_index == 0 or current_local_index == current_index))):
            qoi,pickled_current_model,time_for_qoi = \
                ExecuteInstanceStochasticAdaptiveRefinementAux_Functionality(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_local_index,current_analysis,time_for_qoi,open_mp_threads,mapping_flag,pickled_coarsest_model,print_to_file,filename)
            del(pickled_coarse_model)
            pickled_coarse_model = pickled_current_model
            del(pickled_current_model)
    return qoi,time_for_qoi


############################# StochasticAdaptiveRefinementMultipleTasks ############################

# @ExaquteTask(filename=FILE_OUT,returns=3)
@constraint(ComputingUnits="${computing_units_mlmc_execute_0}")
@ExaquteTask(returns=3)
def ExecuteInstanceStochasticAdaptiveRefinementMultipleTasksAuxLev0_Task(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_local_index,current_analysis,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,filename):
    try:
        open_mp_threads = int(os.environ["computing_units_mlmc_execute_0"])
        threadpool_limits(limits=open_mp_threads)
    except:
        open_mp_threads = 1
    qoi,pickled_current_model,time_for_qoi = \
        ExecuteInstanceStochasticAdaptiveRefinementAux_Functionality(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_local_index,current_analysis,time_for_qoi,open_mp_threads,mapping_flag,pickled_mapping_reference_model,print_to_file,filename)
    return qoi,pickled_current_model,time_for_qoi

# @ExaquteTask(filename=FILE_OUT,returns=3)
@constraint(ComputingUnits="${computing_units_mlmc_execute_1}")
@ExaquteTask(returns=3)
def ExecuteInstanceStochasticAdaptiveRefinementMultipleTasksAuxLev1_Task(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_local_index,current_analysis,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,filename):
    try:
        open_mp_threads = int(os.environ["computing_units_mlmc_execute_1"])
        threadpool_limits(limits=open_mp_threads)
    except:
        open_mp_threads = 1
    qoi,pickled_current_model,time_for_qoi = \
        ExecuteInstanceStochasticAdaptiveRefinementAux_Functionality(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_local_index,current_analysis,time_for_qoi,open_mp_threads,mapping_flag,pickled_mapping_reference_model,print_to_file,filename)
    return qoi,pickled_current_model,time_for_qoi

# @ExaquteTask(filename=FILE_OUT,returns=3)
@constraint(ComputingUnits="${computing_units_mlmc_execute_2}")
@ExaquteTask(returns=3)
def ExecuteInstanceStochasticAdaptiveRefinementMultipleTasksAuxLev2_Task(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_local_index,current_analysis,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,filename):
    try:
        open_mp_threads = int(os.environ["computing_units_mlmc_execute_2"])
        threadpool_limits(limits=open_mp_threads)
    except:
        open_mp_threads = 1
    qoi,pickled_current_model,time_for_qoi = \
        ExecuteInstanceStochasticAdaptiveRefinementAux_Functionality(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_local_index,current_analysis,time_for_qoi,open_mp_threads,mapping_flag,pickled_mapping_reference_model,print_to_file,filename)
    return qoi,pickled_current_model,time_for_qoi

# @ExaquteTask(filename=FILE_OUT,returns=3)
@constraint(ComputingUnits="${computing_units_mlmc_execute_3}")
@ExaquteTask(returns=3)
def ExecuteInstanceStochasticAdaptiveRefinementMultipleTasksAuxLev3_Task(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_local_index,current_analysis,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,filename):
    try:
        open_mp_threads = int(os.environ["computing_units_mlmc_execute_3"])
        threadpool_limits(limits=open_mp_threads)
    except:
        open_mp_threads = 1
    qoi,pickled_current_model,time_for_qoi = \
        ExecuteInstanceStochasticAdaptiveRefinementAux_Functionality(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_local_index,current_analysis,time_for_qoi,open_mp_threads,mapping_flag,pickled_mapping_reference_model,print_to_file,filename)
    return qoi,pickled_current_model,time_for_qoi

# @ExaquteTask(filename=FILE_OUT,returns=3)
@constraint(ComputingUnits="${computing_units_mlmc_execute_4}")
@ExaquteTask(returns=3)
def ExecuteInstanceStochasticAdaptiveRefinementMultipleTasksAuxLev4_Task(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_local_index,current_analysis,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,filename):
    try:
        open_mp_threads = int(os.environ["computing_units_mlmc_execute_4"])
        threadpool_limits(limits=open_mp_threads)
    except:
        open_mp_threads = 1
    qoi,pickled_current_model,time_for_qoi = \
        ExecuteInstanceStochasticAdaptiveRefinementAux_Functionality(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_local_index,current_analysis,time_for_qoi,open_mp_threads,mapping_flag,pickled_mapping_reference_model,print_to_file,filename)
    return qoi,pickled_current_model,time_for_qoi

# @ExaquteTask(filename=FILE_OUT,returns=3)
@constraint(ComputingUnits="${computing_units_mlmc_execute_5}")
@ExaquteTask(returns=3)
def ExecuteInstanceStochasticAdaptiveRefinementMultipleTasksAuxLev5_Task(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_local_index,current_analysis,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,filename):
    try:
        open_mp_threads = int(os.environ["computing_units_mlmc_execute_5"])
        threadpool_limits(limits=open_mp_threads)
    except:
        open_mp_threads = 1
    qoi,pickled_current_model,time_for_qoi = \
        ExecuteInstanceStochasticAdaptiveRefinementAux_Functionality(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_local_index,current_analysis,time_for_qoi,open_mp_threads,mapping_flag,pickled_mapping_reference_model,print_to_file,filename)
    return qoi,pickled_current_model,time_for_qoi


########################################## DeterministicAdaptiveRefinement ########################################

# @ExaquteTask(filename=FILE_OUT,returns=2)
@constraint(ComputingUnits="${computing_units_mlmc_execute_0}")
@ExaquteTask(returns=2)
def executeInstanceDeterministicAdaptiveRefinementAuxLev0_Task(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,filename):
    try:
        open_mp_threads = int(os.environ["computing_units_mlmc_execute_0"])
        threadpool_limits(limits=open_mp_threads)
    except:
        open_mp_threads = 1
    qoi,time_for_qoi = \
        ExecuteInstanceDeterministicAdaptiveRefinementAux_Functionality(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,filename,open_mp_threads)
    return qoi,time_for_qoi

# @ExaquteTask(filename=FILE_OUT,returns=2)
@constraint(ComputingUnits="${computing_units_mlmc_execute_1}")
@ExaquteTask(returns=2)
def executeInstanceDeterministicAdaptiveRefinementAuxLev1_Task(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,filename):
    try:
        open_mp_threads = int(os.environ["computing_units_mlmc_execute_1"])
        threadpool_limits(limits=open_mp_threads)
    except:
        open_mp_threads = 1
    qoi,time_for_qoi = \
        ExecuteInstanceDeterministicAdaptiveRefinementAux_Functionality(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,filename,open_mp_threads)
    return qoi,time_for_qoi

# @ExaquteTask(filename=FILE_OUT,returns=2)
@constraint(ComputingUnits="${computing_units_mlmc_execute_2}")
@ExaquteTask(returns=2)
def executeInstanceDeterministicAdaptiveRefinementAuxLev2_Task(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,filename):
    try:
        open_mp_threads = int(os.environ["computing_units_mlmc_execute_2"])
        threadpool_limits(limits=open_mp_threads)
    except:
        open_mp_threads = 1
    qoi,time_for_qoi = \
        ExecuteInstanceDeterministicAdaptiveRefinementAux_Functionality(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,filename,open_mp_threads)
    return qoi,time_for_qoi

# @ExaquteTask(filename=FILE_OUT,returns=2)
@constraint(ComputingUnits="${computing_units_mlmc_execute_3}")
@ExaquteTask(returns=2)
def executeInstanceDeterministicAdaptiveRefinementAuxLev3_Task(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,filename):
    try:
        open_mp_threads = int(os.environ["computing_units_mlmc_execute_3"])
        threadpool_limits(limits=open_mp_threads)
    except:
        open_mp_threads = 1
    qoi,time_for_qoi = \
        ExecuteInstanceDeterministicAdaptiveRefinementAux_Functionality(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,filename,open_mp_threads)
    return qoi,time_for_qoi

# @ExaquteTask(filename=FILE_OUT,returns=2)
@constraint(ComputingUnits="${computing_units_mlmc_execute_4}")
@ExaquteTask(returns=2)
def executeInstanceDeterministicAdaptiveRefinementAuxLev4_Task(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,filename):
    try:
        open_mp_threads = int(os.environ["computing_units_mlmc_execute_4"])
        threadpool_limits(limits=open_mp_threads)
    except:
        open_mp_threads = 1
    qoi,time_for_qoi = \
        ExecuteInstanceDeterministicAdaptiveRefinementAux_Functionality(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,filename,open_mp_threads)
    return qoi,time_for_qoi

# @ExaquteTask(filename=FILE_OUT,returns=2)
@constraint(ComputingUnits="${computing_units_mlmc_execute_5}")
@ExaquteTask(returns=2)
def executeInstanceDeterministicAdaptiveRefinementAuxLev5_Task(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,filename):
    try:
        open_mp_threads = int(os.environ["computing_units_mlmc_execute_5"])
        threadpool_limits(limits=open_mp_threads)
    except:
        open_mp_threads = 1
    qoi,time_for_qoi = \
        ExecuteInstanceDeterministicAdaptiveRefinementAux_Functionality(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,filename,open_mp_threads)
    return qoi,time_for_qoi


########################################## ReadingFromFile #########################################

# @ExaquteTask(filename=FILE_OUT,returns=2)
@constraint(ComputingUnits="${computing_units_mlmc_execute_0}")
@ExaquteTask(returns=2)
def executeInstanceReadingFromFileAuxLev0_Task(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,filename):
    try:
        open_mp_threads = int(os.environ["computing_units_mlmc_execute_0"])
        threadpool_limits(limits=open_mp_threads)
    except:
        open_mp_threads = 1
    qoi,time_for_qoi = \
        ExecuteInstanceReadingFromFileAux_Functionality(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,filename,open_mp_threads)
    return qoi,time_for_qoi

# @ExaquteTask(filename=FILE_OUT,returns=2)
@constraint(ComputingUnits="${computing_units_mlmc_execute_1}")
@ExaquteTask(returns=2)
def executeInstanceReadingFromFileAuxLev1_Task(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,filename):
    try:
        open_mp_threads = int(os.environ["computing_units_mlmc_execute_1"])
        threadpool_limits(limits=open_mp_threads)
    except:
        open_mp_threads = 1
    qoi,time_for_qoi = \
        ExecuteInstanceReadingFromFileAux_Functionality(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,filename,open_mp_threads)
    return qoi,time_for_qoi

# @ExaquteTask(filename=FILE_OUT,returns=2)
@constraint(ComputingUnits="${computing_units_mlmc_execute_2}")
@ExaquteTask(returns=2)
def executeInstanceReadingFromFileAuxLev2_Task(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,filename):
    try:
        open_mp_threads = int(os.environ["computing_units_mlmc_execute_2"])
        threadpool_limits(limits=open_mp_threads)
    except:
        open_mp_threads = 1
    qoi,time_for_qoi = \
        ExecuteInstanceReadingFromFileAux_Functionality(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,filename,open_mp_threads)
    return qoi,time_for_qoi

# @ExaquteTask(filename=FILE_OUT,returns=2)
@constraint(ComputingUnits="${computing_units_mlmc_execute_3}")
@ExaquteTask(returns=2)
def executeInstanceReadingFromFileAuxLev3_Task(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,filename):
    try:
        open_mp_threads = int(os.environ["computing_units_mlmc_execute_3"])
        threadpool_limits(limits=open_mp_threads)
    except:
        open_mp_threads = 1
    qoi,time_for_qoi = \
        ExecuteInstanceReadingFromFileAux_Functionality(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,filename,open_mp_threads)
    return qoi,time_for_qoi

# @ExaquteTask(filename=FILE_OUT,returns=2)
@constraint(ComputingUnits="${computing_units_mlmc_execute_4}")
@ExaquteTask(returns=2)
def executeInstanceReadingFromFileAuxLev4_Task(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,filename):
    try:
        open_mp_threads = int(os.environ["computing_units_mlmc_execute_4"])
        threadpool_limits(limits=open_mp_threads)
    except:
        open_mp_threads = 1
    qoi,time_for_qoi = \
        ExecuteInstanceReadingFromFileAux_Functionality(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,filename,open_mp_threads)
    return qoi,time_for_qoi

# @ExaquteTask(filename=FILE_OUT,returns=2)
@constraint(ComputingUnits="${computing_units_mlmc_execute_5}")
@ExaquteTask(returns=2)
def executeInstanceReadingFromFileAuxLev5_Task(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,filename):
    try:
        open_mp_threads = int(os.environ["computing_units_mlmc_execute_5"])
        threadpool_limits(limits=open_mp_threads)
    except:
        open_mp_threads = 1
    qoi,time_for_qoi = \
        ExecuteInstanceReadingFromFileAux_Functionality(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,filename,open_mp_threads)
    return qoi,time_for_qoi


    ############################################# ZeroTask #############################################

@ExaquteTask(returns=2)
def returnZeroQoiAndTime_Task(number_scalar, number_vector, size_vector):
    if number_scalar > 0:
        qoi_scalar = [0.0 for _ in range(number_scalar)]
    if number_vector > 0:
        qoi_vector = [[0.0 for _ in range(size_vector)] for _ in range(number_vector)]
    if number_scalar > 0 and number_vector > 0:
        qoi = qoi_scalar + qoi_vector
    elif number_scalar > 0:
        qoi = qoi_scalar
    elif number_vector > 0:
        qoi = qoi_vector
    time_for_qoi = 0.0
    return qoi,time_for_qoi


####################################################################################################
######################################### FUNCTIONALITIES ##########################################
####################################################################################################


def ExecuteInstanceDeterministicAdaptiveRefinementAux_Functionality(pickled_model,pickled_project_parameters,current_analysis_stage,random_variable,previous_computational_time,mapping_flag,pickled_mapping_reference_model,print_to_file,filename,open_mp_threads):
    """
    Auxiliary method to the solve method of the KratosSolverWrapper class. The problem is solved calling Kratos. To be called if the selected refinement strategy is deterministic_adaptive_refinement.

    Inputs:
    - pickled_model: serialization of the model.
    - pickled_project_parameters: serialization of the parameters.
    - current_analysis_stage: Kratos analysis stage of the problem.
    - random_variable: list containing all random variables.
    - previous_computational_time: time to solution of previous indices, if any.
    - mapping_flag: booleans. Defines if in current_analysis_stage = SimulationScenario class mapping is required.
    - pickled_mapping_reference_model: Kratos Model to which mapping is performed.
    - print_to_file: boolean. If true, the specific filename is passed to the analysis stage.
    - filename: string which defines the name of the file the task will write, if any.
    - open_mp_threads: number of threads we are exploiting to solve current task.

    Outputs:
    - qoi: list containing all the quantities of interest.
    - computational_time: time to solution up to current index.
    """

    start_time = time.time()
    # overwrite the old model serializer with the unpickled one
    model_serializer = pickle.loads(pickled_model)
    current_model = KratosMultiphysics.Model()
    model_serializer.Load("ModelSerialization",current_model)
    del(model_serializer)
    # overwrite the old parameters serializer with the unpickled one
    serialized_project_parameters = pickle.loads(pickled_project_parameters)
    current_project_parameters = KratosMultiphysics.Parameters()
    serialized_project_parameters.Load("ParametersSerialization",current_project_parameters)
    del(serialized_project_parameters)

    # constructor analysis stage
    simulation = current_analysis_stage(current_model,current_project_parameters,random_variable)
    # add filename flag print_to_file is true
    if (print_to_file):
        simulation.filename = filename
    # add flag if current index is maximum index: always True
    simulation.is_current_index_maximum_index = True
    # mapping if in current finest level (always true) and mapping flag is true
    # otherwise standard behavior
    if (mapping_flag is True):
        # unpickle mapping reference model and build Kratos Model object
        serialized_mapping_reference_model = pickle.loads(pickled_mapping_reference_model)
        mapping_reference_model = KratosMultiphysics.Model()
        serialized_mapping_reference_model.Load("ModelSerialization",mapping_reference_model)
        del(serialized_mapping_reference_model)
        # send reference model to analysis stage for mapping and set mapping flag to true
        simulation.mapping_reference_model = mapping_reference_model
        simulation.mapping = True
    simulation.Run()
    # mapping if in current finest level and mapping flag is true
    # otherwise standard qoi evaluation
    if (mapping_flag is True):
        qoi = simulation.MappingAndEvaluateQuantityOfInterest()
    else:
        qoi = simulation.EvaluateQuantityOfInterest()
    del(simulation)
    end_time = time.time()
    computational_time = previous_computational_time + open_mp_threads*(end_time-start_time) # multiply by open mp threads to consider real machine cost
    return qoi,computational_time


def ExecuteInstanceStochasticAdaptiveRefinementAux_Functionality(current_global_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_index,current_analysis_stage,previous_computational_time,open_mp_threads,mapping_flag,pickled_mapping_reference_model,print_to_file,filename):
    """
    Auxiliary method to the solve method of the KratosSolverWrapper class. Firstly, if needed, adaptive refinement is performed. Then the problem is solved calling Kratos. To be called if the selected refinement strategy is stochastic_adaptive_refinement.

    Inputs:
    - current_global_index: index we are interested in and for which we have called the solve method of the solver wrapper.
    - pickled_coarse_model: serialization of the model.
    - pickled_coarse_project_parameters: serialization of the parameters.
    - pickled_custom_metric_refinement_parameters: serialization of metric adaptive refinement parameters.
    - pickled_custom_remesh_refinement_parameters: serialization of remeshing adaptive refinement parameters.
    - random_variable: list containing all random variables.
    - current_index: index we are currently solving at this moment. May be different from current_global_index for "stocharstic_adaptive_refinement" strategy.
    - current_analysis_stage: Kratos analysis stage of the problem.
    - previous_computational_time: time to solution of previous indices, if any.
    - open_mp_threads: number of threads we are exploiting to solve current task.
    - mapping_flag: booleans. Defines if in current_analysis_stage = SimulationScenario class mapping is required.
    - pickled_mapping_reference_model: Kratos Model to which mapping is performed.
    - print_to_file: boolean. If true, the specific filename is passed to the analysis stage.
    - filename: string which defines the name of the file the task will write, if any.

    Outputs:
    - qoi: quantities of interest.
    - pickled_finer_model: serialized adaptively refined model, used at current index.
    - computational_time: time to solution up to current index.
    """

    start_time = time.time()
    # unpickle model and build Kratos Model object
    serialized_model = pickle.loads(pickled_coarse_model)
    current_model = KratosMultiphysics.Model()
    serialized_model.Load("ModelSerialization",current_model)
    del(serialized_model)
    # unpickle parameters and build Kratos Parameters object
    serialized_project_parameters = pickle.loads(pickled_coarse_project_parameters)
    current_project_parameters = KratosMultiphysics.Parameters()
    serialized_project_parameters.Load("ParametersSerialization",current_project_parameters)
    del(serialized_project_parameters)
    # refine if current current_global_index > 0, adaptive refinement based on the solution of previous index
    if (current_index > 0):
        # unpickle metric and remesh refinement parameters and build Kratos Parameters objects
        serialized_custom_metric_refinement_parameters = pickle.loads(pickled_custom_metric_refinement_parameters)
        serialized_custom_remesh_refinement_parameters = pickle.loads(pickled_custom_remesh_refinement_parameters)
        current_custom_metric_refinement_parameters = KratosMultiphysics.Parameters()
        current_custom_remesh_refinement_parameters = KratosMultiphysics.Parameters()
        serialized_custom_metric_refinement_parameters.Load("MetricRefinementParametersSerialization",current_custom_metric_refinement_parameters)
        serialized_custom_remesh_refinement_parameters.Load("RemeshRefinementParametersSerialization",current_custom_remesh_refinement_parameters)
        del(serialized_custom_metric_refinement_parameters,serialized_custom_remesh_refinement_parameters)
        # refine the model Kratos object
        adaptive_refinement_manager = AdaptiveRefinement(current_index,current_model,current_project_parameters,current_custom_metric_refinement_parameters,current_custom_remesh_refinement_parameters)
        refined_model,refined_project_parameters = adaptive_refinement_manager.ComputeAdaptiveRefinement()
        current_model = refined_model
        del(refined_model,refined_project_parameters)
    # constructor analysis stage
    simulation = current_analysis_stage(current_model,current_project_parameters,random_variable)
    # add filename flag print_to_file is true
    if (print_to_file):
        simulation.filename = filename
    # add flag if current index is maximum index
    if (current_index == current_global_index):
        simulation.is_current_index_maximum_index = True
    else:
        simulation.is_current_index_maximum_index = False
    # mapping if in current finest level and mapping flag is true
    # otherwise standard behavior
    if (mapping_flag is True and current_index == current_global_index):
        # unpickle mapping reference model and build Kratos Model object
        serialized_mapping_reference_model = pickle.loads(pickled_mapping_reference_model)
        mapping_reference_model = KratosMultiphysics.Model()
        serialized_mapping_reference_model.Load("ModelSerialization",mapping_reference_model)
        del(serialized_mapping_reference_model)
        # send reference model to analysis stage for mapping and set mapping flag to true
        simulation.mapping_reference_model = mapping_reference_model
        simulation.mapping = True
    simulation.Run()
    # mapping if in current finest level and mapping flag is true
    # otherwise standard qoi evaluation
    if (mapping_flag is True and current_index == current_global_index):
        qoi = simulation.MappingAndEvaluateQuantityOfInterest()
    else:
        qoi = simulation.EvaluateQuantityOfInterest()
    # save model and parameters as MpiSerializer Kratos objects
    serialized_finer_model = KratosMultiphysics.MpiSerializer()
    serialized_finer_model.Save("ModelSerialization",simulation.model)
    # pickle model and parameters
    pickled_finer_model = pickle.dumps(serialized_finer_model, 2) # second argument is the protocol and is NECESSARY (according to pybind11 docs)
    del(simulation)
    end_time = time.time()
    computational_time = previous_computational_time + open_mp_threads*(end_time-start_time) # multiply by open mp threads to consider real machine cost
    return qoi,pickled_finer_model,computational_time


def ExecuteInstanceReadingFromFileAux_Functionality(pickled_model,pickled_project_parameters,current_analysis_stage,random_variable,previous_computational_time,mapping_flag,pickled_mapping_reference_model,print_to_file,filename,open_mp_threads):
    """
    Auxiliary method to the solve method of the KratosSolverWrapper class. The problem is solved calling Kratos. To be called if the selected refinement strategy is reading_from_file.

    Inputs:
    - pickled_model: serialization of the model.
    - pickled_project_parameters: serialization of the parameters.
    - current_analysis_stage: Kratos analysis stage of the problem.
    - random_variable: list containing all random variables.
    - previous_computational_time: time to solution of previous indices, if any.
    - mapping_flag: booleans. Defines if in current_analysis_stage = SimulationScenario class mapping is required.
    - pickled_mapping_reference_model: Kratos Model to which mapping is performed.
    - print_to_file: boolean. If true, the specific filename is passed to the analysis stage.
    - filename: string which defines the name of the file the task will write, if any.
    - open_mp_threads: number of threads we are exploiting to solve current task.

    Outputs:
    - qoi: quantities of interest.
    - computational_time: time to solution up to current index.
    """

    start_time = time.time()
    # unpickle model and build Kratos Model object
    serialized_model = pickle.loads(pickled_model)
    current_model = KratosMultiphysics.Model()
    serialized_model.Load("ModelSerialization",current_model)
    del(serialized_model)
    # unpickle parameters and build Kratos Parameters object
    serialized_project_parameters = pickle.loads(pickled_project_parameters)
    current_project_parameters = KratosMultiphysics.Parameters()
    serialized_project_parameters.Load("ParametersSerialization",current_project_parameters)
    del(serialized_project_parameters)
    # constructor analysis stage
    simulation = current_analysis_stage(current_model,current_project_parameters,random_variable)
    # add filename flag print_to_file is true
    if (print_to_file):
        simulation.filename = filename
    # add flag if current index is maximum index: always True
    simulation.is_current_index_maximum_index = True
    # mapping if in current finest level (always true) and mapping flag is true
    # otherwise standard behavior
    if (mapping_flag is True):
        # unpickle mapping reference model and build Kratos Model object
        serialized_mapping_reference_model = pickle.loads(pickled_mapping_reference_model)
        mapping_reference_model = KratosMultiphysics.Model()
        serialized_mapping_reference_model.Load("ModelSerialization",mapping_reference_model)
        del(serialized_mapping_reference_model)
        # send reference model to analysis stage for mapping and set mapping flag to true
        simulation.mapping_reference_model = mapping_reference_model
        simulation.mapping = True
    simulation.Run()
    # mapping if in current finest level and mapping flag is true
    # otherwise standard qoi evaluation
    if (mapping_flag is True):
        qoi = simulation.MappingAndEvaluateQuantityOfInterest()
    else:
        qoi = simulation.EvaluateQuantityOfInterest()
    del(simulation)
    end_time = time.time()
    computational_time = previous_computational_time + open_mp_threads*(end_time-start_time) # multiply by open mp threads to consider real machine cost
    return qoi,computational_time
