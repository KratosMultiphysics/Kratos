# Import Python libraries
import time
import pickle
import os
try:
    from threadpoolctl import *
except:
    pass

# Import XMC, distributed environment
from xmc.distributedEnvironmentFramework import *
from xmc.classDefs_solverWrapper.methodDefs_KratosSolverWrapper.solve import ExecuteInstanceDeterministicAdaptiveRefinementAux_Functionality,ExecuteInstanceReadingFromFileAux_Functionality,ExecuteInstanceStochasticAdaptiveRefinementAux_Functionality

try:
    computing_procs_mlmc_execute_0 = int(os.environ["computing_procs_mlmc_execute_0"])
except:
    computing_procs_mlmc_execute_0 = 1


####################################################################################################
########################################## SERIALIZATION ###########################################
####################################################################################################

@constraint(computing_units="${computing_units_mlmc_execute_0}")
@mpi(runner="mpirun", processes=computing_procs_mlmc_execute_0)
@ExaquteTask(returns=computing_procs_mlmc_execute_0)
def SerializeMPIModel(pickled_parameters, main_model_part_name, fake_sample_to_serialize, analysis):

    import KratosMultiphysics
    import KratosMultiphysics.mpi as KratosMPI

    serialized_parameters = pickle.loads(pickled_parameters)
    del pickled_parameters
    deserialized_parameters = KratosMultiphysics.Parameters()
    serialized_parameters.Load("ParametersSerialization", deserialized_parameters)

    # prepare the model to serialize
    model = KratosMultiphysics.Model()
    fake_sample = fake_sample_to_serialize
    deserialized_parameters["solver_settings"]["model_import_settings"]["input_type"].SetString("mdpa")

    simulation = analysis(model,deserialized_parameters,fake_sample)
    simulation.Initialize()
    # reset general flags
    simulation.model.GetModelPart(main_model_part_name).ProcessInfo.SetValue(KratosMultiphysics.IS_RESTARTED,True)

    # serialize model
    serialized_model = KratosMultiphysics.MpiSerializer()
    serialized_model.Save("ModelSerialization",simulation.model)
    # self.serialized_model.append(serialized_model)

    # pickle dataserialized_data
    pickled_model = pickle.dumps(serialized_model, 2) # second argument is the protocol and is NECESSARY (according to pybind11 docs)

    return pickled_model


####################################################################################################
############################################ WRAPPERS ##############################################
####################################################################################################

def executeInstanceStochasticAdaptiveRefinementAllAtOnce_Wrapper(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_analysis,time_for_qoi,mapping_flag,adaptive_refinement_jump_to_finest_level,print_to_file,current_contribution):
    if (current_index == 0):
        qoi_and_time_list = ExecuteInstanceStochasticAdaptiveRefinementAllAtOnceAuxLev0_Task(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_analysis,time_for_qoi,mapping_flag,adaptive_refinement_jump_to_finest_level,print_to_file,"filename_level_"+str(current_index)+"_contribution_"+str(current_contribution)+"_random_variable_"+str(random_variable[0])+".dat")
    elif (current_index == 1):
        qoi_and_time_list = ExecuteInstanceStochasticAdaptiveRefinementAllAtOnceAuxLev1_Task(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_analysis,time_for_qoi,mapping_flag,adaptive_refinement_jump_to_finest_level,print_to_file,"filename_level_"+str(current_index)+"_contribution_"+str(current_contribution)+"_random_variable_"+str(random_variable[0])+".dat")
    elif (current_index == 2):
        qoi_and_time_list = ExecuteInstanceStochasticAdaptiveRefinementAllAtOnceAuxLev2_Task(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_analysis,time_for_qoi,mapping_flag,adaptive_refinement_jump_to_finest_level,print_to_file,"filename_level_"+str(current_index)+"_contribution_"+str(current_contribution)+"_random_variable_"+str(random_variable[0])+".dat")
    else:
        raise Exception("Level not supported")
    qoi, time_for_qoi = UnfoldFutureQT(qoi_and_time_list)
    return qoi, time_for_qoi

def executeInstanceStochasticAdaptiveRefinementMultipleTasks_Wrapper(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_local_index,current_analysis,time_for_qoi,mapping_flag,print_to_file,current_contribution,pickled_mapping_reference_model=None):
    if (current_index == 0):
        qoi_pickled_current_model_time_for_qoi_list = ExecuteInstanceStochasticAdaptiveRefinementMultipleTasksAuxLev0_Task(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_local_index,current_analysis,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,"filename_level_"+str(current_index)+"_contribution_"+str(current_contribution)+"_random_variable_"+str(random_variable[0])+".dat")
    elif (current_index == 1):
        qoi_and_time_list = ExecuteInstanceStochasticAdaptiveRefinementMultipleTasksAuxLev1_Task(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_local_index,current_analysis,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,"filename_level_"+str(current_index)+"_contribution_"+str(current_contribution)+"_random_variable_"+str(random_variable[0])+".dat")
    else:
        raise Exception("Level not supported")
    qoi, pickled_current_model, time_for_qoi = UnfoldFutureQMT(qoi_pickled_current_model_time_for_qoi_list)
    return qoi, pickled_current_model, time_for_qoi

def executeInstanceDeterministicAdaptiveRefinement_Wrapper(current_index,pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,current_contribution):
    if (current_index == 0):
        qoi_and_time_list = executeInstanceDeterministicAdaptiveRefinementAuxLev0_Task(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,"filename_level_"+str(current_index)+"_contribution_"+str(current_contribution)+"_random_variable_"+str(random_variable[0])+".dat")
    elif (current_index == 1):
        qoi_and_time_list = executeInstanceDeterministicAdaptiveRefinementAuxLev1_Task(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,"filename_level_"+str(current_index)+"_contribution_"+str(current_contribution)+"_random_variable_"+str(random_variable[0])+".dat")
    else:
        raise Exception("Level not supported")
    qoi, time_for_qoi = UnfoldFutureQT(qoi_and_time_list)
    return qoi, time_for_qoi

def executeInstanceReadingFromFile_Wrapper(current_index,pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,current_contribution):
    if (current_index == 0):
        qoi_and_time_list = executeInstanceReadingFromFileAuxLev0_Task(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,"filename_level_"+str(current_index)+"_contribution_"+str(current_contribution)+"_random_variable_"+str(random_variable[0])+".dat")
    elif (current_index == 1):
        qoi_and_time_list = executeInstanceReadingFromFileAuxLev1_Task(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,"filename_level_"+str(current_index)+"_contribution_"+str(current_contribution)+"_random_variable_"+str(random_variable[0])+".dat")
    else:
        raise Exception("Level not supported")
    qoi, time_for_qoi = UnfoldFutureQT(qoi_and_time_list)
    return qoi, time_for_qoi


####################################################################################################
############################################## TASKS ###############################################
####################################################################################################

@ExaquteTask(qoi_and_time_list={Type: COLLECTION_IN, Depth: 2}, returns=2)
def UnfoldFutureQT(qoi_and_time_list):
    qoi = qoi_and_time_list[0][0] # get first qoi element (all are equal since they are synchronized)
    time_for_qoi = 0.0
    for qoi_and_time  in qoi_and_time_list:
        time_for_qoi += qoi_and_time[1] # sum all times
    return qoi, time_for_qoi

@ExaquteTask(qoi_pickled_current_model_time_for_qoi_list={Type: COLLECTION_IN, Depth: 2}, returns=3)
def UnfoldFutureQMT(qoi_pickled_current_model_time_for_qoi_list):
    qoi = qoi_pickled_current_model_time_for_qoi_list[0][0] # get first qoi element (all are equal since they are synchronized)
    pickled_current_model = qoi_pickled_current_model_time_for_qoi_list[1]
    time_for_qoi = 0.0
    for qoi_pickled_current_model_time_for_qoi  in qoi_pickled_current_model_time_for_qoi_list:
        time_for_qoi += qoi_pickled_current_model_time_for_qoi[-1] # sum all times
    return qoi, pickled_current_model, time_for_qoi

############################### StochasticAdaptiveRefinementAllAtOnce ##############################

# @ExaquteTask(filename=FILE_OUT,pickled_coarse_model=COLLECTION_IN, returns=computing_procs_mlmc_execute_0)
@constraint(computing_units="${computing_units_mlmc_execute_0}")
@mpi(runner="mpirun", processes=computing_procs_mlmc_execute_0, pickled_coarse_model_layout={block_count: computing_procs_mlmc_execute_0, block_length: 1, stride: 1})
@ExaquteTask(pickled_coarse_model=COLLECTION_IN, returns=computing_procs_mlmc_execute_0)
def ExecuteInstanceStochasticAdaptiveRefinementAllAtOnceAuxLev0_Task(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_analysis,time_for_qoi,mapping_flag,adaptive_refinement_jump_to_finest_level,print_to_file,filename):
    # Import Kratos
    import KratosMultiphysics
    import KratosMultiphysics.mpi as KratosMPI
    from KratosMultiphysics.MultilevelMonteCarloApplication.adaptive_refinement_utilities import AdaptiveRefinement

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

############################# StochasticAdaptiveRefinementMultipleTasks ############################

# @ExaquteTask(filename=FILE_OUT,pickled_coarse_model=COLLECTION_IN, returns=computing_procs_mlmc_execute_0)
@constraint(computing_units="${computing_units_mlmc_execute_0}")
@mpi(runner="mpirun", processes=computing_procs_mlmc_execute_0, pickled_coarse_model_layout={block_count: computing_procs_mlmc_execute_0, block_length: 1, stride: 1})
@ExaquteTask(pickled_coarse_model=COLLECTION_IN, returns=computing_procs_mlmc_execute_0)
def ExecuteInstanceStochasticAdaptiveRefinementMultipleTasksAuxLev0_Task(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_local_index,current_analysis,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,filename):
    # Import Kratos
    import KratosMultiphysics
    import KratosMultiphysics.mpi as KratosMPI
    from KratosMultiphysics.MultilevelMonteCarloApplication.adaptive_refinement_utilities import AdaptiveRefinement

    try:
        open_mp_threads = int(os.environ["computing_units_mlmc_execute_0"])
        threadpool_limits(limits=open_mp_threads)
    except:
        open_mp_threads = 1

    qoi,pickled_current_model,time_for_qoi = \
        ExecuteInstanceStochasticAdaptiveRefinementAux_Functionality(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_local_index,current_analysis,time_for_qoi,open_mp_threads,mapping_flag,pickled_mapping_reference_model,print_to_file,filename)
    return qoi,pickled_current_model,time_for_qoi

########################################## DeterministicAdaptiveRefinement ########################################

# @ExaquteTask(filename=FILE_OUT,pickled_model=COLLECTION_IN, returns=computing_procs_mlmc_execute_0)
@constraint(computing_units="${computing_units_mlmc_execute_0}")
@mpi(runner="mpirun", processes=computing_procs_mlmc_execute_0, pickled_model_layout={block_count: computing_procs_mlmc_execute_0, block_length: 1, stride: 1})
@ExaquteTask(pickled_model=COLLECTION_IN, returns=computing_procs_mlmc_execute_0)
def executeInstanceDeterministicAdaptiveRefinementAuxLev0_Task(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,filename):
    # Import Kratos
    import KratosMultiphysics
    import KratosMultiphysics.mpi as KratosMPI
    from KratosMultiphysics.MultilevelMonteCarloApplication.adaptive_refinement_utilities import AdaptiveRefinement

    try:
        open_mp_threads = int(os.environ["computing_units_mlmc_execute_0"])
        threadpool_limits(limits=open_mp_threads)
    except:
        open_mp_threads = 1

    qoi,time_for_qoi = \
        ExecuteInstanceDeterministicAdaptiveRefinementAux_Functionality(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,filename,open_mp_threads)
    return qoi,time_for_qoi

########################################## ReadingFromFile #########################################

# @ExaquteTask(filename=FILE_OUT,pickled_model=COLLECTION_IN, returns=computing_procs_mlmc_execute_0)
@constraint(computing_units="${computing_units_mlmc_execute_0}")
@mpi(runner="mpirun", processes=computing_procs_mlmc_execute_0, pickled_model_layout={block_count: computing_procs_mlmc_execute_0, block_length: 1, stride: 1})
@ExaquteTask(pickled_model=COLLECTION_IN, returns=computing_procs_mlmc_execute_0)
def executeInstanceReadingFromFileAuxLev0_Task(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,filename):
    # Import Kratos
    import KratosMultiphysics
    import KratosMultiphysics.mpi as KratosMPI
    from KratosMultiphysics.MultilevelMonteCarloApplication.adaptive_refinement_utilities import AdaptiveRefinement

    try:
        open_mp_threads = int(os.environ["computing_units_mlmc_execute_0"])
        threadpool_limits(limits=open_mp_threads)
    except:
        open_mp_threads = 1

    qoi,time_for_qoi = \
        ExecuteInstanceReadingFromFileAux_Functionality(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,filename,open_mp_threads)
    return qoi,time_for_qoi
