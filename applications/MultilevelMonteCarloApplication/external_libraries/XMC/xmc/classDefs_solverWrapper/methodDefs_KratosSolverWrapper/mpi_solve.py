# Import Python libraries
import time
import pickle
import os
try:
    from threadpoolctl import *
except:
    pass

# Import Kratos, XMC, distributed environment
from KratosMultiphysics import IsDistributedRun, DataCommunicator
from xmc.classDefs_solverWrapper.methodDefs_KratosSolverWrapper.solve import ExecuteInstanceDeterministicAdaptiveRefinementAux_Functionality,ExecuteInstanceReadingFromFileAux_Functionality,ExecuteInstanceStochasticAdaptiveRefinementAux_Functionality
from exaqute import *

computing_units_mlmc_execute_0 = int(os.getenv("computing_units_mlmc_execute_0", 1))
computing_units_mlmc_execute_1 = int(os.getenv("computing_units_mlmc_execute_1", 1))
computing_units_mlmc_execute_2 = int(os.getenv("computing_units_mlmc_execute_2", 1))

computing_procs_mlmc_execute_0 = int(os.getenv("computing_procs_mlmc_execute_0", 1))
computing_procs_mlmc_execute_1 = int(os.getenv("computing_procs_mlmc_execute_1", 1))
computing_procs_mlmc_execute_2 = int(os.getenv("computing_procs_mlmc_execute_2", 1))

ppn_mlmc_execute_0 = int(os.getenv("ppn_mlmc_execute_0", 1))
ppn_mlmc_execute_1 = int(os.getenv("ppn_mlmc_execute_1", 1))
ppn_mlmc_execute_2 = int(os.getenv("ppn_mlmc_execute_2", 1))

####################################################################################################
############################################ WRAPPERS ##############################################
####################################################################################################

def SerializeMPIModel_Wrapper(pickled_parameters, main_model_part_name, fake_sample_to_serialize, analysis, current_index):
    if current_index == 0:
        pickled_model = SerializeMPIModelAuxLev0_Task(pickled_parameters, main_model_part_name, fake_sample_to_serialize, analysis)
    elif current_index == 1:
        pickled_model = SerializeMPIModelAuxLev1_Task(pickled_parameters, main_model_part_name, fake_sample_to_serialize, analysis)
    elif current_index == 2:
        pickled_model = SerializeMPIModelAuxLev2_Task(pickled_parameters, main_model_part_name, fake_sample_to_serialize, analysis)
    else:
        raise Exception("Level not supported")
    return pickled_model

def SerializeDeterministicAdaptiveRefinementMPIModel_Wrapper(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_analysis,time_for_qoi,adaptive_refinement_jump_to_finest_level):
    if current_index == 0:
        pickled_model = SerializeDeterministicAdaptiveRefinementMPIModelAuxLev0_Task(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_analysis,time_for_qoi,adaptive_refinement_jump_to_finest_level)
    elif current_index == 1:
        pickled_model = SerializeDeterministicAdaptiveRefinementMPIModelAuxLev1_Task(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_analysis,time_for_qoi,adaptive_refinement_jump_to_finest_level)
    elif current_index == 2:
        pickled_model = SerializeDeterministicAdaptiveRefinementMPIModelAuxLev2_Task(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_analysis,time_for_qoi,adaptive_refinement_jump_to_finest_level)
    else:
        raise Exception("Level not supported")
    return pickled_model

def executeInstanceStochasticAdaptiveRefinementAllAtOnce_Wrapper(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_analysis,time_for_qoi,mapping_flag,adaptive_refinement_jump_to_finest_level,print_to_file,current_contribution):
    if (current_index == 0):
        qoi_and_time_list = ExecuteInstanceStochasticAdaptiveRefinementAllAtOnceAuxLev0_Task(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_analysis,time_for_qoi,mapping_flag,adaptive_refinement_jump_to_finest_level,print_to_file,"filename_level_"+str(current_index)+"_contribution_"+str(current_contribution)+"_random_variable_"+str(random_variable[0])+".dat")
    elif (current_index == 1):
        qoi_and_time_list = ExecuteInstanceStochasticAdaptiveRefinementAllAtOnceAuxLev1_Task(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_analysis,time_for_qoi,mapping_flag,adaptive_refinement_jump_to_finest_level,print_to_file,"filename_level_"+str(current_index)+"_contribution_"+str(current_contribution)+"_random_variable_"+str(random_variable[0])+".dat")
    elif (current_index == 2):
        qoi_and_time_list = ExecuteInstanceStochasticAdaptiveRefinementAllAtOnceAuxLev2_Task(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_analysis,time_for_qoi,mapping_flag,adaptive_refinement_jump_to_finest_level,print_to_file,"filename_level_"+str(current_index)+"_contribution_"+str(current_contribution)+"_random_variable_"+str(random_variable[0])+".dat")
    else:
        raise Exception("Level not supported")
    if IsDistributedRun():
        # running with mpirun the whole xmc algorithm
        qoi, time_for_qoi = UnfoldQT(qoi_and_time_list)
    else:
        # running with distributed environment framework, only Kratos tasks are run with mpi
        qoi, time_for_qoi = UnfoldFutureQT(qoi_and_time_list)
    return qoi, time_for_qoi

def executeInstanceStochasticAdaptiveRefinementMultipleTasks_Wrapper(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_local_index,current_analysis,time_for_qoi,mapping_flag,print_to_file,current_contribution,pickled_mapping_reference_model=None):
    if (current_index == 0):
        qoi_pickled_current_model_time_for_qoi_list = ExecuteInstanceStochasticAdaptiveRefinementMultipleTasksAuxLev0_Task(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_local_index,current_analysis,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,"filename_level_"+str(current_index)+"_contribution_"+str(current_contribution)+"_random_variable_"+str(random_variable[0])+".dat")
    else:
        # We cannot run with multiple tasks, since tasks of different levels are normally run with different number of processors,
        # and when running with MPI the model should be pickled with the number of processes of the task.
        # For example, if I want to run with MPI and 4 processes, I need to serialize within an MPI task of 4 processes.
        raise Exception("Level not supported. You should set \"taskAllAtOnce\" to \"true\" to run multi-level algorithms with \"stochastic_adaptive_refinement\" as \"refinement_strategy\".")
    if IsDistributedRun():
        # running with mpirun the whole xmc algorithm
        qoi, pickled_current_model, time_for_qoi = UnfoldQMT(qoi_pickled_current_model_time_for_qoi_list)
    else:
        # running with distributed environment framework, only Kratos tasks are run with mpi
        qoi, pickled_current_model, time_for_qoi = UnfoldFutureQMT(qoi_pickled_current_model_time_for_qoi_list)
    return qoi, pickled_current_model, time_for_qoi

def executeInstanceDeterministicAdaptiveRefinement_Wrapper(current_index,pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,current_contribution):
    if (current_index == 0):
        qoi_and_time_list = executeInstanceDeterministicAdaptiveRefinementAuxLev0_Task(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,"filename_level_"+str(current_index)+"_contribution_"+str(current_contribution)+"_random_variable_"+str(random_variable[0])+".dat")
    elif (current_index == 1):
        qoi_and_time_list = executeInstanceDeterministicAdaptiveRefinementAuxLev1_Task(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,"filename_level_"+str(current_index)+"_contribution_"+str(current_contribution)+"_random_variable_"+str(random_variable[0])+".dat")
    elif (current_index == 2):
        qoi_and_time_list = executeInstanceDeterministicAdaptiveRefinementAuxLev2_Task(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,"filename_level_"+str(current_index)+"_contribution_"+str(current_contribution)+"_random_variable_"+str(random_variable[0])+".dat")
    else:
        raise Exception("Level not supported")
    if IsDistributedRun():
        # running with mpirun the whole xmc algorithm
        qoi, time_for_qoi = UnfoldQT(qoi_and_time_list)
    else:
        # running with distributed environment framework, only Kratos tasks are run with mpi
        qoi, time_for_qoi = UnfoldFutureQT(qoi_and_time_list)
    return qoi, time_for_qoi

def executeInstanceReadingFromFile_Wrapper(current_index,pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,current_contribution):
    if (current_index == 0):
        qoi_and_time_list = executeInstanceReadingFromFileAuxLev0_Task(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,"filename_level_"+str(current_index)+"_contribution_"+str(current_contribution)+"_random_variable_"+str(random_variable[0])+".dat")
    elif (current_index == 1):
        qoi_and_time_list = executeInstanceReadingFromFileAuxLev1_Task(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,"filename_level_"+str(current_index)+"_contribution_"+str(current_contribution)+"_random_variable_"+str(random_variable[0])+".dat")
    elif (current_index == 2):
        qoi_and_time_list = executeInstanceReadingFromFileAuxLev2_Task(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,"filename_level_"+str(current_index)+"_contribution_"+str(current_contribution)+"_random_variable_"+str(random_variable[0])+".dat")
    else:
        raise Exception("Level not supported")
    if IsDistributedRun():
        # running with mpirun the whole xmc algorithm
        qoi, time_for_qoi = UnfoldQT(qoi_and_time_list)
    else:
        # running with distributed environment framework, only Kratos tasks are run with mpi
        qoi, time_for_qoi = UnfoldFutureQT(qoi_and_time_list)
    return qoi, time_for_qoi


####################################################################################################
############################################## TASKS ###############################################
####################################################################################################

@task(keep=True, returns=2)
def UnfoldQT(qoi_and_time_list):
    communicator = DataCommunicator.GetDefault()
    qoi = qoi_and_time_list[0]
    time_for_qoi = communicator.SumAll(qoi_and_time_list[-1])
    return qoi, time_for_qoi

@task(keep=True, returns=3)
def UnfoldQMT(qoi_pickled_current_model_time_for_qoi_list):
    communicator = DataCommunicator.GetDefault()
    qoi = qoi_pickled_current_model_time_for_qoi_list[0]
    pickled_current_model = qoi_pickled_current_model_time_for_qoi_list[1]
    time_for_qoi = communicator.SumAll(qoi_pickled_current_model_time_for_qoi_list[-1])
    return qoi, pickled_current_model, time_for_qoi

@task(keep=True, qoi_and_time_list={Type: COLLECTION_IN, Depth: 2}, returns=2)
def UnfoldFutureQT(qoi_and_time_list):
    qoi = qoi_and_time_list[0][0] # get first qoi element (all are equal since they are synchronized)
    time_for_qoi = 0.0
    for qoi_and_time  in qoi_and_time_list:
        time_for_qoi += qoi_and_time[1] # sum all times
    return qoi, time_for_qoi

@task(keep=True, qoi_pickled_current_model_time_for_qoi_list={Type: COLLECTION_IN, Depth: 2}, returns=3)
def UnfoldFutureQMT(qoi_pickled_current_model_time_for_qoi_list):
    qoi = qoi_pickled_current_model_time_for_qoi_list[0][0] # get first qoi element (all are equal since they are synchronized)
    pickled_current_model = qoi_pickled_current_model_time_for_qoi_list[1]
    time_for_qoi = 0.0
    for qoi_pickled_current_model_time_for_qoi  in qoi_pickled_current_model_time_for_qoi_list:
        time_for_qoi += qoi_pickled_current_model_time_for_qoi[-1] # sum all times
    return qoi, pickled_current_model, time_for_qoi

########################################## Serialization ##########################################

@constraint(computing_units=computing_units_mlmc_execute_0)
@mpi(runner="mpirun", processes=computing_procs_mlmc_execute_0, processes_per_node=ppn_mlmc_execute_0)
@task(keep=True, returns=computing_procs_mlmc_execute_0)
def SerializeMPIModelAuxLev0_Task(pickled_parameters, main_model_part_name, fake_sample_to_serialize, analysis):
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
    # initialize analysis stage
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

@constraint(computing_units=computing_units_mlmc_execute_1)
@mpi(runner="mpirun", processes=computing_procs_mlmc_execute_1, processes_per_node=ppn_mlmc_execute_1)
@task(keep=True, returns=computing_procs_mlmc_execute_1)
def SerializeMPIModelAuxLev1_Task(pickled_parameters, main_model_part_name, fake_sample_to_serialize, analysis):
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
    # initialize analysis stage
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

@constraint(computing_units=computing_units_mlmc_execute_2)
@mpi(runner="mpirun", processes=computing_procs_mlmc_execute_2, processes_per_node=ppn_mlmc_execute_2)
@task(keep=True, returns=computing_procs_mlmc_execute_2)
def SerializeMPIModelAuxLev2_Task(pickled_parameters, main_model_part_name, fake_sample_to_serialize, analysis):
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
    # initialize analysis stage
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

########################################## Serialization DAR ##########################################

@constraint(computing_units=computing_units_mlmc_execute_0)
@mpi(runner="mpirun", processes=computing_procs_mlmc_execute_0, processes_per_node=ppn_mlmc_execute_0, pickled_coarse_model_layout={block_count: computing_procs_mlmc_execute_0, block_length: 1, stride: 1})
@task(keep=True, pickled_coarse_model=COLLECTION_IN, returns=computing_procs_mlmc_execute_0)
def SerializeDeterministicAdaptiveRefinementMPIModelAuxLev0_Task(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_analysis,time_for_qoi,adaptive_refinement_jump_to_finest_level):
    # Import Kratos
    import KratosMultiphysics
    import KratosMultiphysics.mpi as KratosMPI
    from KratosMultiphysics.MultilevelMonteCarloApplication.adaptive_refinement_utilities import AdaptiveRefinement

    try:
        open_mp_threads = computing_units_mlmc_execute_0
        threadpool_limits(limits=open_mp_threads)
    except:
        open_mp_threads = 1

    mapping_flag = False
    print_to_file = False
    filename = ""
    pickled_coarsest_model = pickled_coarse_model
    for current_local_index in range(current_index+1):
        if ((adaptive_refinement_jump_to_finest_level is False) or (adaptive_refinement_jump_to_finest_level is True and (current_local_index == 0 or current_local_index == current_index))):
            qoi,pickled_current_model,time_for_qoi = \
                ExecuteInstanceStochasticAdaptiveRefinementAux_Functionality(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_local_index,current_analysis,time_for_qoi,open_mp_threads,mapping_flag,pickled_coarsest_model,print_to_file,filename)
            del(pickled_coarse_model)
            pickled_coarse_model = pickled_current_model
            del(pickled_current_model)
    return pickled_coarse_model

@constraint(computing_units=computing_units_mlmc_execute_1)
@mpi(runner="mpirun", processes=computing_procs_mlmc_execute_1, processes_per_node=ppn_mlmc_execute_1, pickled_coarse_model_layout={block_count: computing_procs_mlmc_execute_1, block_length: 1, stride: 1})
@task(keep=True, pickled_coarse_model=COLLECTION_IN, returns=computing_procs_mlmc_execute_1)
def SerializeDeterministicAdaptiveRefinementMPIModelAuxLev1_Task(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_analysis,time_for_qoi,adaptive_refinement_jump_to_finest_level):
    # Import Kratos
    import KratosMultiphysics
    import KratosMultiphysics.mpi as KratosMPI
    from KratosMultiphysics.MultilevelMonteCarloApplication.adaptive_refinement_utilities import AdaptiveRefinement

    try:
        open_mp_threads = computing_units_mlmc_execute_1
        threadpool_limits(limits=open_mp_threads)
    except:
        open_mp_threads = 1

    mapping_flag = False
    print_to_file = False
    filename = ""
    pickled_coarsest_model = pickled_coarse_model
    for current_local_index in range(current_index+1):
        if ((adaptive_refinement_jump_to_finest_level is False) or (adaptive_refinement_jump_to_finest_level is True and (current_local_index == 0 or current_local_index == current_index))):
            qoi,pickled_current_model,time_for_qoi = \
                ExecuteInstanceStochasticAdaptiveRefinementAux_Functionality(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_local_index,current_analysis,time_for_qoi,open_mp_threads,mapping_flag,pickled_coarsest_model,print_to_file,filename)
            del(pickled_coarse_model)
            pickled_coarse_model = pickled_current_model
            del(pickled_current_model)
    return pickled_coarse_model

@constraint(computing_units=computing_units_mlmc_execute_2)
@mpi(runner="mpirun", processes=computing_procs_mlmc_execute_2, processes_per_node=ppn_mlmc_execute_2, pickled_coarse_model_layout={block_count: computing_procs_mlmc_execute_2, block_length: 1, stride: 1})
@task(keep=True, pickled_coarse_model=COLLECTION_IN, returns=computing_procs_mlmc_execute_2)
def SerializeDeterministicAdaptiveRefinementMPIModelAuxLev2_Task(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_analysis,time_for_qoi,adaptive_refinement_jump_to_finest_level):
    # Import Kratos
    import KratosMultiphysics
    import KratosMultiphysics.mpi as KratosMPI
    from KratosMultiphysics.MultilevelMonteCarloApplication.adaptive_refinement_utilities import AdaptiveRefinement

    try:
        open_mp_threads = computing_units_mlmc_execute_2
        threadpool_limits(limits=open_mp_threads)
    except:
        open_mp_threads = 1

    mapping_flag = False
    print_to_file = False
    filename = ""
    pickled_coarsest_model = pickled_coarse_model
    for current_local_index in range(current_index+1):
        if ((adaptive_refinement_jump_to_finest_level is False) or (adaptive_refinement_jump_to_finest_level is True and (current_local_index == 0 or current_local_index == current_index))):
            qoi,pickled_current_model,time_for_qoi = \
                ExecuteInstanceStochasticAdaptiveRefinementAux_Functionality(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_local_index,current_analysis,time_for_qoi,open_mp_threads,mapping_flag,pickled_coarsest_model,print_to_file,filename)
            del(pickled_coarse_model)
            pickled_coarse_model = pickled_current_model
            del(pickled_current_model)
    return pickled_coarse_model

############################### StochasticAdaptiveRefinementAllAtOnce ##############################

# @task(keep=True, filename=FILE_OUT, pickled_coarse_model=COLLECTION_IN, returns=computing_procs_mlmc_execute_0)
@constraint(computing_units=computing_units_mlmc_execute_0)
@mpi(runner="mpirun", processes=computing_procs_mlmc_execute_0, processes_per_node=ppn_mlmc_execute_0, pickled_coarse_model_layout={block_count: computing_procs_mlmc_execute_0, block_length: 1, stride: 1})
@task(keep=True, pickled_coarse_model=COLLECTION_IN, returns=computing_procs_mlmc_execute_0)
def ExecuteInstanceStochasticAdaptiveRefinementAllAtOnceAuxLev0_Task(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_analysis,time_for_qoi,mapping_flag,adaptive_refinement_jump_to_finest_level,print_to_file,filename):
    # Import Kratos
    import KratosMultiphysics
    import KratosMultiphysics.mpi as KratosMPI
    from KratosMultiphysics.MultilevelMonteCarloApplication.adaptive_refinement_utilities import AdaptiveRefinement

    try:
        open_mp_threads = computing_units_mlmc_execute_0
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

# @task(keep=True, filename=FILE_OUT, pickled_coarse_model=COLLECTION_IN, returns=computing_procs_mlmc_execute_1)
@constraint(computing_units=computing_units_mlmc_execute_1)
@mpi(runner="mpirun", processes=computing_procs_mlmc_execute_1, processes_per_node=ppn_mlmc_execute_1, pickled_coarse_model_layout={block_count: computing_procs_mlmc_execute_1, block_length: 1, stride: 1})
@task(keep=True, pickled_coarse_model=COLLECTION_IN, returns=computing_procs_mlmc_execute_1)
def ExecuteInstanceStochasticAdaptiveRefinementAllAtOnceAuxLev1_Task(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_analysis,time_for_qoi,mapping_flag,adaptive_refinement_jump_to_finest_level,print_to_file,filename):
    # Import Kratos
    import KratosMultiphysics
    import KratosMultiphysics.mpi as KratosMPI
    from KratosMultiphysics.MultilevelMonteCarloApplication.adaptive_refinement_utilities import AdaptiveRefinement

    try:
        open_mp_threads = computing_units_mlmc_execute_1
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

# @task(keep=True, filename=FILE_OUT, pickled_coarse_model=COLLECTION_IN, returns=computing_procs_mlmc_execute_2)
@constraint(computing_units=computing_units_mlmc_execute_2)
@mpi(runner="mpirun", processes=computing_procs_mlmc_execute_2, processes_per_node=ppn_mlmc_execute_2, pickled_coarse_model_layout={block_count: computing_procs_mlmc_execute_2, block_length: 1, stride: 1})
@task(keep=True, pickled_coarse_model=COLLECTION_IN, returns=computing_procs_mlmc_execute_2)
def ExecuteInstanceStochasticAdaptiveRefinementAllAtOnceAuxLev2_Task(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_analysis,time_for_qoi,mapping_flag,adaptive_refinement_jump_to_finest_level,print_to_file,filename):
    # Import Kratos
    import KratosMultiphysics
    import KratosMultiphysics.mpi as KratosMPI
    from KratosMultiphysics.MultilevelMonteCarloApplication.adaptive_refinement_utilities import AdaptiveRefinement

    try:
        open_mp_threads = computing_units_mlmc_execute_2
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

# @task(keep=True, filename=FILE_OUT,pickled_coarse_model=COLLECTION_IN, returns=computing_procs_mlmc_execute_0)
@constraint(computing_units=computing_units_mlmc_execute_0)
@mpi(runner="mpirun", processes=computing_procs_mlmc_execute_0, processes_per_node=ppn_mlmc_execute_0, pickled_coarse_model_layout={block_count: computing_procs_mlmc_execute_0, block_length: 1, stride: 1})
@task(keep=True, pickled_coarse_model=COLLECTION_IN, returns=computing_procs_mlmc_execute_0)
def ExecuteInstanceStochasticAdaptiveRefinementMultipleTasksAuxLev0_Task(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_local_index,current_analysis,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,filename):
    # Import Kratos
    import KratosMultiphysics
    import KratosMultiphysics.mpi as KratosMPI
    from KratosMultiphysics.MultilevelMonteCarloApplication.adaptive_refinement_utilities import AdaptiveRefinement

    try:
        open_mp_threads = computing_units_mlmc_execute_0
        threadpool_limits(limits=open_mp_threads)
    except:
        open_mp_threads = 1

    qoi,pickled_current_model,time_for_qoi = \
        ExecuteInstanceStochasticAdaptiveRefinementAux_Functionality(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_local_index,current_analysis,time_for_qoi,open_mp_threads,mapping_flag,pickled_mapping_reference_model,print_to_file,filename)
    return qoi,pickled_current_model,time_for_qoi

########################################## DeterministicAdaptiveRefinement ########################################

# @task(keep=True, filename=FILE_OUT,pickled_model=COLLECTION_IN, pickled_mapping_reference_model=COLLECTION_IN, returns=computing_procs_mlmc_execute_0)
@constraint(computing_units=computing_units_mlmc_execute_0)
@mpi(runner="mpirun", processes=computing_procs_mlmc_execute_0, processes_per_node=ppn_mlmc_execute_0, pickled_model_layout={block_count: computing_procs_mlmc_execute_0, block_length: 1, stride: 1}, pickled_mapping_reference_model_layout={block_count: computing_procs_mlmc_execute_0, block_length: 1, stride: 1})
@task(keep=True, pickled_model=COLLECTION_IN, pickled_mapping_reference_model=COLLECTION_IN, returns=computing_procs_mlmc_execute_0)
def executeInstanceDeterministicAdaptiveRefinementAuxLev0_Task(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,filename):
    # Import Kratos
    import KratosMultiphysics
    import KratosMultiphysics.mpi as KratosMPI
    from KratosMultiphysics.MultilevelMonteCarloApplication.adaptive_refinement_utilities import AdaptiveRefinement

    try:
        open_mp_threads = computing_units_mlmc_execute_0
        threadpool_limits(limits=open_mp_threads)
    except:
        open_mp_threads = 1

    qoi,time_for_qoi = \
        ExecuteInstanceDeterministicAdaptiveRefinementAux_Functionality(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,filename,open_mp_threads)
    return qoi,time_for_qoi

# @task(keep=True, filename=FILE_OUT,pickled_model=COLLECTION_IN, pickled_mapping_reference_model=COLLECTION_IN, returns=computing_procs_mlmc_execute_1)
@constraint(computing_units=computing_units_mlmc_execute_1)
@mpi(runner="mpirun", processes=computing_procs_mlmc_execute_1, processes_per_node=ppn_mlmc_execute_1, pickled_model_layout={block_count: computing_procs_mlmc_execute_1, block_length: 1, stride: 1}, pickled_mapping_reference_model_layout={block_count: computing_procs_mlmc_execute_1, block_length: 1, stride: 1})
@task(keep=True, pickled_model=COLLECTION_IN, pickled_mapping_reference_model=COLLECTION_IN, returns=computing_procs_mlmc_execute_1)
def executeInstanceDeterministicAdaptiveRefinementAuxLev1_Task(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,filename):
    # Import Kratos
    import KratosMultiphysics
    import KratosMultiphysics.mpi as KratosMPI
    from KratosMultiphysics.MultilevelMonteCarloApplication.adaptive_refinement_utilities import AdaptiveRefinement

    try:
        open_mp_threads = computing_units_mlmc_execute_1
        threadpool_limits(limits=open_mp_threads)
    except:
        open_mp_threads = 1

    qoi,time_for_qoi = \
        ExecuteInstanceDeterministicAdaptiveRefinementAux_Functionality(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,filename,open_mp_threads)
    return qoi,time_for_qoi

# @task(keep=True, filename=FILE_OUT,pickled_model=COLLECTION_IN, pickled_mapping_reference_model=COLLECTION_IN, returns=computing_procs_mlmc_execute_2)
@constraint(computing_units=computing_units_mlmc_execute_2)
@mpi(runner="mpirun", processes=computing_procs_mlmc_execute_2, processes_per_node=ppn_mlmc_execute_2, pickled_model_layout={block_count: computing_procs_mlmc_execute_2, block_length: 1, stride: 1}, pickled_mapping_reference_model_layout={block_count: computing_procs_mlmc_execute_2, block_length: 1, stride: 1})
@task(keep=True, pickled_model=COLLECTION_IN, pickled_mapping_reference_model=COLLECTION_IN, returns=computing_procs_mlmc_execute_2)
def executeInstanceDeterministicAdaptiveRefinementAuxLev2_Task(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,filename):
    # Import Kratos
    import KratosMultiphysics
    import KratosMultiphysics.mpi as KratosMPI
    from KratosMultiphysics.MultilevelMonteCarloApplication.adaptive_refinement_utilities import AdaptiveRefinement

    try:
        open_mp_threads = computing_units_mlmc_execute_2
        threadpool_limits(limits=open_mp_threads)
    except:
        open_mp_threads = 1

    qoi,time_for_qoi = \
        ExecuteInstanceDeterministicAdaptiveRefinementAux_Functionality(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,filename,open_mp_threads)
    return qoi,time_for_qoi

########################################## ReadingFromFile #########################################

# @task(keep=True, filename=FILE_OUT, pickled_model=COLLECTION_IN, pickled_mapping_reference_model=COLLECTION_IN, returns=computing_procs_mlmc_execute_0)
@constraint(computing_units=computing_units_mlmc_execute_0)
@mpi(runner="mpirun", processes=computing_procs_mlmc_execute_0, processes_per_node=ppn_mlmc_execute_0, pickled_model_layout={block_count: computing_procs_mlmc_execute_0, block_length: 1, stride: 1}, pickled_mapping_reference_model_layout={block_count: computing_procs_mlmc_execute_0, block_length: 1, stride: 1})
@task(keep=True, pickled_model=COLLECTION_IN, pickled_mapping_reference_model=COLLECTION_IN, returns=computing_procs_mlmc_execute_0)
def executeInstanceReadingFromFileAuxLev0_Task(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,filename):
    # Import Kratos
    import KratosMultiphysics
    import KratosMultiphysics.mpi as KratosMPI
    from KratosMultiphysics.MultilevelMonteCarloApplication.adaptive_refinement_utilities import AdaptiveRefinement

    try:
        open_mp_threads = computing_units_mlmc_execute_0
        threadpool_limits(limits=open_mp_threads)
    except:
        open_mp_threads = 1

    qoi,time_for_qoi = \
        ExecuteInstanceReadingFromFileAux_Functionality(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,filename,open_mp_threads)
    return qoi,time_for_qoi

# @task(keep=True, filename=FILE_OUT, pickled_model=COLLECTION_IN, pickled_mapping_reference_model=COLLECTION_IN, returns=computing_procs_mlmc_execute_1)
@constraint(computing_units=computing_units_mlmc_execute_1)
@mpi(runner="mpirun", processes=computing_procs_mlmc_execute_1, processes_per_node=ppn_mlmc_execute_1, pickled_model_layout={block_count: computing_procs_mlmc_execute_1, block_length: 1, stride: 1}, pickled_mapping_reference_model_layout={block_count: computing_procs_mlmc_execute_1, block_length: 1, stride: 1})
@task(keep=True, pickled_model=COLLECTION_IN, pickled_mapping_reference_model=COLLECTION_IN, returns=computing_procs_mlmc_execute_1)
def executeInstanceReadingFromFileAuxLev1_Task(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,filename):
    # Import Kratos
    import KratosMultiphysics
    import KratosMultiphysics.mpi as KratosMPI
    from KratosMultiphysics.MultilevelMonteCarloApplication.adaptive_refinement_utilities import AdaptiveRefinement

    try:
        open_mp_threads = computing_units_mlmc_execute_1
        threadpool_limits(limits=open_mp_threads)
    except:
        open_mp_threads = 1

    qoi,time_for_qoi = \
        ExecuteInstanceReadingFromFileAux_Functionality(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,filename,open_mp_threads)
    return qoi,time_for_qoi

# @task(keep=True, filename=FILE_OUT, pickled_model=COLLECTION_IN, pickled_mapping_reference_model=COLLECTION_IN, returns=computing_procs_mlmc_execute_2)
@constraint(computing_units=computing_units_mlmc_execute_2)
@mpi(runner="mpirun", processes=computing_procs_mlmc_execute_2, processes_per_node=ppn_mlmc_execute_2, pickled_model_layout={block_count: computing_procs_mlmc_execute_2, block_length: 1, stride: 1}, pickled_mapping_reference_model_layout={block_count: computing_procs_mlmc_execute_2, block_length: 1, stride: 1})
@task(keep=True, pickled_model=COLLECTION_IN, pickled_mapping_reference_model=COLLECTION_IN, returns=computing_procs_mlmc_execute_2)
def executeInstanceReadingFromFileAuxLev2_Task(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,filename):
    # Import Kratos
    import KratosMultiphysics
    import KratosMultiphysics.mpi as KratosMPI
    from KratosMultiphysics.MultilevelMonteCarloApplication.adaptive_refinement_utilities import AdaptiveRefinement

    try:
        open_mp_threads = computing_units_mlmc_execute_2
        threadpool_limits(limits=open_mp_threads)
    except:
        open_mp_threads = 1

    qoi,time_for_qoi = \
        ExecuteInstanceReadingFromFileAux_Functionality(pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,filename,open_mp_threads)
    return qoi,time_for_qoi
