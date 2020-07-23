from __future__ import absolute_import, division # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics.MultilevelMonteCarloApplication.tools import ParametersWrapper

# Import packages
import numpy as np
from scipy.stats import norm
import copy
import time
import sys
import os
try:
    from threadpoolctl import *
except:
    pass

# Importing the analysis stage classes of the different problems
from simulation_definition import SimulationScenario

# Import the StatisticalVariable class
from KratosMultiphysics.MultilevelMonteCarloApplication.statistical_variable_utilities import StatisticalVariable

# Import refinement library
from KratosMultiphysics.MultilevelMonteCarloApplication.adaptive_refinement_utilities import AdaptiveRefinement

# Import random variable generator
import KratosMultiphysics.MultilevelMonteCarloApplication.generator_utilities as generator

# Import PyCOMPSs
# from exaqute.ExaquteTaskPyCOMPSs import *   # to execute with runcompss
# from exaqute.ExaquteTaskHyperLoom import *  # to execute with the IT4 scheduler
from exaqute.ExaquteTaskLocal import *      # to execute with python3

# Import cpickle to pickle the serializer
try:
    import cpickle as pickle  # Use cPickle on Python 2.7
except ImportError:
    import pickle


"""
auxiliary function of AddResults of the MultilevelMonteCarlo class
input:  level              : working level
        simulation_results : instances of the MultilevelMonteCarloResults class
output: new_difference_QoI_power_sum_batches : power sums up to level 4 of difference QoI
        new_time_ML_power_sum_batches        : power sums up to level 4 of time MLMC instance
"""
@constraint(ComputingUnits="${computing_units_auxiliar_utilities}")
@ExaquteTask(returns=2,priority=True)
def AddResultsAux_Task(level,*simulation_results):
    def computePowers(value):
        return value, value**2, value**3, value**4
    if (level == 0):
        # each value is inside the relative level list, and only one value per level is computed
        # i.e. results = [[value_level_0],[value_level_1],...]
        difference_QoI_values = list(map(lambda x: x.QoI[level][0], simulation_results))
    else:
        difference_QoI_values = list(map(lambda x: x.QoI[level][0] - x.QoI[level-1][0], simulation_results))
    time_ML_list = list(map(lambda x: x.time_ML[level][0], simulation_results))
    new_difference_QoI_power_sum_batches = np.sum(np.array(list(map(computePowers, difference_QoI_values))), axis=0)
    new_time_ML_power_sum_batches = np.sum(np.array(list(map(computePowers, time_ML_list))), axis=0)
    return new_difference_QoI_power_sum_batches,new_time_ML_power_sum_batches


"""
auxiliary function finalizing the screening phase and the MLMC phase of the MultilevelMonteCarlo class
input:  ConstructorCallback       : an instance of the MultilevelMonteCarlo class
        aux_settings_serialized   : StreamSerializer Kratos class containing the settings of the MultilevelMonteCarlo class
        aux_mesh_parameters       : mesh parameters of the MltilevelMOnteCarlo class
        aux_current_number_levels : current number of levels of the MultilevelMonteCarlo class
        aux_current_iteration     : current iteration of the MultilevelMonteCarlo class
        aux_number_samples        : number of samples of the Multilevel Monte Carlo class
        *args                     : list containing difference_QoI.raw_moment_1, difference_QoI.unbiased_central_moment_2 and time_ML.raw_moment_1
output: auxiliary_MLMC_object.rates_error       : rates error of the MultilevelMonteCarlo class
        auxiliary_MLMC_object.bayesian_variance : bayesian variance of the MultilevelMonteCarlo class
        auxiliary_MLMC_object.mean_mlmc_QoI     : multilevel monte carlo mean estimator of the MultilevelMonteCarlo class
        auxiliary_MLMC_object.total_error       :  total error of the MultilevelMonteCarlo class
        auxiliary_MLMC_object.number_samples    : number of samples of the MultilevelMonteCarlo class
"""
@constraint(ComputingUnits="${computing_units_auxiliar_utilities}")
@ExaquteTask(returns=5,priority=True)
def FinalizePhaseAux_Task(ConstructorCallback,aux_settings_serialized,aux_mesh_parameters,\
aux_current_number_levels,aux_current_iteration,aux_number_samples,*args):
    # retrieve lists
    args = list(args)
    first_occ = args.index("%%%")
    difference_QoI_mean = args[:first_occ]
    second_occ = args[first_occ + 1:].index("%%%")
    difference_QoI_sample_variance = args[first_occ + 1:first_occ + 1 + second_occ]
    time_ML_mean = args[first_occ + 1 + second_occ + 1:]
    # load settings (Kratos Parameters object) from SteramSerializer Kratos object
    aux_settings = KratosMultiphysics.Parameters()
    aux_settings_serialized.Load("ParametersSerialization",aux_settings)
    # create an auxiliary object auxiliary_MLMC_object equivalent to the current one of the problem
    # construct the class with the aux_settings settings passed from outside
    aux_analysis = "auxiliary_analysis"
    aux_project_parameters_path = "auxiliary_project_parameters_path"
    aux_parameters_refinement_path = "auxiliary_parameters_refinement_path"
    auxiliary_MLMC_object = ConstructorCallback(aux_settings,aux_project_parameters_path,aux_parameters_refinement_path,aux_analysis)
    auxiliary_MLMC_object.difference_QoI.h_statistics_1 = difference_QoI_mean
    auxiliary_MLMC_object.difference_QoI.h_statistics_2 = difference_QoI_sample_variance
    auxiliary_MLMC_object.time_ML.h_statistics_1 = time_ML_mean
    auxiliary_MLMC_object.mesh_parameters = aux_mesh_parameters
    auxiliary_MLMC_object.current_number_levels = aux_current_number_levels
    auxiliary_MLMC_object.current_iteration = aux_current_iteration
    auxiliary_MLMC_object.number_samples = aux_number_samples
    # compute the functions needed to finalize the screening phase or the MLMC phase
    if (auxiliary_MLMC_object.current_iteration == 0): # screening phase
        # compute parameters by least square fit
        auxiliary_MLMC_object.ComputeRatesLS()
        # compute Bayesian variance V^c[Y_l]
        auxiliary_MLMC_object.EstimateBayesianVariance(auxiliary_MLMC_object.current_number_levels)
    else: # MLMC phase
        # compute estimator MLMC mean QoI
        auxiliary_MLMC_object.ComputeMeanMLMCQoI()
        # compute parameters by least square fit
        auxiliary_MLMC_object.ComputeRatesLS()
        # compute Bayesian variance V^c[Y_l]
        auxiliary_MLMC_object.EstimateBayesianVariance(auxiliary_MLMC_object.current_number_levels)
        # compute total error of MLMC simulation
        auxiliary_MLMC_object.ComputeTotalErrorMLMC()
    return auxiliary_MLMC_object.rates_error,\
    auxiliary_MLMC_object.bayesian_variance,auxiliary_MLMC_object.mean_mlmc_QoI,\
    auxiliary_MLMC_object.total_error,auxiliary_MLMC_object.number_samples

"""
function wrapper of all the different tasks for the concurrent adaptive refinement all at once strategy
"""
def ExecuteInstanceConcurrentAdaptiveRefinementAllAtOnce_Wrapper(current_MLMC_level,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,sample,current_analysis,mlmc_results):
    if (current_MLMC_level == 0):
        mlmc_results = ExecuteInstanceConcurrentAdaptiveRefinementAllAtOnceAuxLev0_Task(current_MLMC_level,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,sample,current_analysis,mlmc_results)
    elif (current_MLMC_level == 1):
        mlmc_results = ExecuteInstanceConcurrentAdaptiveRefinementAllAtOnceAuxLev1_Task(current_MLMC_level,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,sample,current_analysis,mlmc_results)
    elif (current_MLMC_level == 2):
        mlmc_results = ExecuteInstanceConcurrentAdaptiveRefinementAllAtOnceAuxLev2_Task(current_MLMC_level,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,sample,current_analysis,mlmc_results)
    else:
        raise Exception("Level not supported")
    return mlmc_results

"""
function wrapper of all the different tasks for the concurrent adaptive refinement multiple tasks strategy
"""
def ExecuteInstanceConcurrentAdaptiveRefinementMultipleTasks_Wrapper(current_MLMC_level,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,sample,current_level,current_analysis,mlmc_results):
    if (current_MLMC_level == 0):
        mlmc_results,pickled_current_model = ExecuteInstanceConcurrentAdaptiveRefinementMultipleTasksAuxLev0_Task(current_MLMC_level,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,sample,current_level,current_analysis,mlmc_results)
    elif (current_MLMC_level == 1):
        mlmc_results,pickled_current_model = ExecuteInstanceConcurrentAdaptiveRefinementMultipleTasksAuxLev1_Task(current_MLMC_level,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,sample,current_level,current_analysis,mlmc_results)
    elif (current_MLMC_level == 2):
        mlmc_results,pickled_current_model = ExecuteInstanceConcurrentAdaptiveRefinementMultipleTasksAuxLev2_Task(current_MLMC_level,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,sample,current_level,current_analysis,mlmc_results)
    else:
        raise Exception("Level not supported")
    return mlmc_results,pickled_current_model

"""
function wrapper of all the different tasks for computing adaptive refinement
"""
def ExecuteInstanceOnlyAdaptiveRefinement_Wrapper(pickled_model,pickled_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_level,current_analysis,mlmc_results):
    if (current_level == 0):
        pickled_model,mlmc_results = executeInstanceOnlyAdaptiveRefinementAuxLev0_Task(pickled_model,pickled_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_level,current_analysis,mlmc_results)
    elif (current_level == 1):
        pickled_model,mlmc_results = executeInstanceOnlyAdaptiveRefinementAuxLev1_Task(pickled_model,pickled_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_level,current_analysis,mlmc_results)
    elif (current_level == 2):
        pickled_model,mlmc_results = executeInstanceOnlyAdaptiveRefinementAuxLev2_Task(pickled_model,pickled_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_level,current_analysis,mlmc_results)
    elif (current_level == 3):
        pickled_model,mlmc_results = executeInstanceOnlyAdaptiveRefinementAuxLev3_Task(pickled_model,pickled_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_level,current_analysis,mlmc_results)
    elif (current_level == 4):
        pickled_model,mlmc_results = executeInstanceOnlyAdaptiveRefinementAuxLev4_Task(pickled_model,pickled_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_level,current_analysis,mlmc_results)
    elif (current_level == 5):
        pickled_model,mlmc_results = executeInstanceOnlyAdaptiveRefinementAuxLev5_Task(pickled_model,pickled_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_level,current_analysis,mlmc_results)
    else:
        raise Exception("Level not supported")
    return mlmc_results,pickled_model

"""
function wrapper of all the different tasks for the single adaptive refinement strategy
"""
def ExecuteInstanceSingleRefinement_Wrapper(current_MLMC_level,pickled_model,pickled_project_parameters,random_variable,current_level,current_analysis,mlmc_results):
    if (current_MLMC_level == 0):
        mlmc_results = ExecuteInstanceSingleRefinementAuxLev0_Task(pickled_model,pickled_project_parameters,random_variable,current_level,current_analysis,mlmc_results)
    elif (current_MLMC_level == 1):
        mlmc_results = ExecuteInstanceSingleRefinementAuxLev1_Task(pickled_model,pickled_project_parameters,random_variable,current_level,current_analysis,mlmc_results)
    elif (current_MLMC_level == 2):
        mlmc_results = ExecuteInstanceSingleRefinementAuxLev2_Task(pickled_model,pickled_project_parameters,random_variable,current_level,current_analysis,mlmc_results)
    elif (current_MLMC_level == 3):
        mlmc_results = ExecuteInstanceSingleRefinementAuxLev3_Task(pickled_model,pickled_project_parameters,random_variable,current_level,current_analysis,mlmc_results)
    elif (current_MLMC_level == 4):
        mlmc_results = ExecuteInstanceSingleRefinementAuxLev4_Task(pickled_model,pickled_project_parameters,random_variable,current_level,current_analysis,mlmc_results)
    elif (current_MLMC_level == 5):
        mlmc_results = ExecuteInstanceSingleRefinementAuxLev5_Task(pickled_model,pickled_project_parameters,random_variable,current_level,current_analysis,mlmc_results)
    else:
        raise Exception("Level not supported")
    return mlmc_results


############################### ConcurrentAdaptiveRefinementAllAtOnce ##############################

@constraint(ComputingUnits="${computing_units_mlmc_execute_0}")
@ExaquteTask(returns=1)
def ExecuteInstanceConcurrentAdaptiveRefinementAllAtOnceAuxLev0_Task(current_MLMC_level,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,sample,current_analysis,mlmc_results):
    try:
        open_mp_threads = int(os.environ["computing_units_mlmc_execute_0"])
        threadpool_limits(limits=open_mp_threads)
    except:
        pass
    for current_level in range(current_MLMC_level+1):
        mlmc_results,pickled_current_model = \
            ExecuteInstanceConcurrentAdaptiveRefinementAux_Functionality(current_MLMC_level,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,sample,current_level,current_analysis,mlmc_results)
        del(pickled_coarse_model)
        pickled_coarse_model = pickled_current_model
        del(pickled_current_model)
    return mlmc_results

@constraint(ComputingUnits="${computing_units_mlmc_execute_1}")
@ExaquteTask(returns=1)
def ExecuteInstanceConcurrentAdaptiveRefinementAllAtOnceAuxLev1_Task(current_MLMC_level,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,sample,current_analysis,mlmc_results):
    try:
        open_mp_threads = int(os.environ["computing_units_mlmc_execute_1"])
        threadpool_limits(limits=open_mp_threads)
    except:
        pass
    for current_level in range(current_MLMC_level+1):
        mlmc_results,pickled_current_model = \
            ExecuteInstanceConcurrentAdaptiveRefinementAux_Functionality(current_MLMC_level,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,sample,current_level,current_analysis,mlmc_results)
        del(pickled_coarse_model)
        pickled_coarse_model = pickled_current_model
        del(pickled_current_model)
    return mlmc_results

@constraint(ComputingUnits="${computing_units_mlmc_execute_2}")
@ExaquteTask(returns=1)
def ExecuteInstanceConcurrentAdaptiveRefinementAllAtOnceAuxLev2_Task(current_MLMC_level,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,sample,current_analysis,mlmc_results):
    try:
        open_mp_threads = int(os.environ["computing_units_mlmc_execute_2"])
        threadpool_limits(limits=open_mp_threads)
    except:
        pass
    for current_level in range(current_MLMC_level+1):
        mlmc_results,pickled_current_model = \
            ExecuteInstanceConcurrentAdaptiveRefinementAux_Functionality(current_MLMC_level,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,sample,current_level,current_analysis,mlmc_results)
        del(pickled_coarse_model)
        pickled_coarse_model = pickled_current_model
        del(pickled_current_model)
    return mlmc_results


############################# ConcurrentAdaptiveRefinementMultipleTasks ############################

@constraint(ComputingUnits="${computing_units_mlmc_execute_0}")
@ExaquteTask(returns=2)
def ExecuteInstanceConcurrentAdaptiveRefinementMultipleTasksAuxLev0_Task(current_MLMC_level,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,sample,current_level,current_analysis,mlmc_results):
    try:
        open_mp_threads = int(os.environ["computing_units_mlmc_execute_0"])
        threadpool_limits(limits=open_mp_threads)
    except:
        pass
    mlmc_results,pickled_current_model = \
        ExecuteInstanceConcurrentAdaptiveRefinementAux_Functionality(current_MLMC_level,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,sample,current_level,current_analysis,mlmc_results)
    return mlmc_results,pickled_current_model

@constraint(ComputingUnits="${computing_units_mlmc_execute_1}")
@ExaquteTask(returns=2)
def ExecuteInstanceConcurrentAdaptiveRefinementMultipleTasksAuxLev1_Task(current_MLMC_level,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,sample,current_level,current_analysis,mlmc_results):
    try:
        open_mp_threads = int(os.environ["computing_units_mlmc_execute_1"])
        threadpool_limits(limits=open_mp_threads)
    except:
        pass
    mlmc_results,pickled_current_model = \
        ExecuteInstanceConcurrentAdaptiveRefinementAux_Functionality(current_MLMC_level,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,sample,current_level,current_analysis,mlmc_results)
    return mlmc_results,pickled_current_model

@constraint(ComputingUnits="${computing_units_mlmc_execute_2}")
@ExaquteTask(returns=2)
def ExecuteInstanceConcurrentAdaptiveRefinementMultipleTasksAuxLev2_Task(current_MLMC_level,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,sample,current_level,current_analysis,mlmc_results):
    try:
        open_mp_threads = int(os.environ["computing_units_mlmc_execute_2"])
        threadpool_limits(limits=open_mp_threads)
    except:
        pass
    mlmc_results,pickled_current_model = \
        ExecuteInstanceConcurrentAdaptiveRefinementAux_Functionality(current_MLMC_level,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,sample,current_level,current_analysis,mlmc_results)
    return mlmc_results,pickled_current_model


################################### OnlySingleAdaptiveRefinement ###################################

@constraint(ComputingUnits="${computing_units_mlmc_execute_0}")
@ExaquteTask(returns=2)
def executeInstanceOnlyAdaptiveRefinementAuxLev0_Task(pickled_model,pickled_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_index,current_analysis,mlmc_results):
    try:
        open_mp_threads = int(os.environ["computing_units_mlmc_execute_0"])
        threadpool_limits(limits=open_mp_threads)
    except:
        pass
    pickled_model,mlmc_results = \
        ExecuteInstanceOnlyAdaptiveRefinementAux_Functionality(pickled_model,pickled_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_index,current_analysis,mlmc_results)
    return pickled_model,mlmc_results

@constraint(ComputingUnits="${computing_units_mlmc_execute_1}")
@ExaquteTask(returns=2)
def executeInstanceOnlyAdaptiveRefinementAuxLev1_Task(pickled_model,pickled_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_index,current_analysis,mlmc_results):
    try:
        open_mp_threads = int(os.environ["computing_units_mlmc_execute_1"])
        threadpool_limits(limits=open_mp_threads)
    except:
        pass
    pickled_model,mlmc_results = \
        ExecuteInstanceOnlyAdaptiveRefinementAux_Functionality(pickled_model,pickled_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_index,current_analysis,mlmc_results)
    return pickled_model,mlmc_results

@constraint(ComputingUnits="${computing_units_mlmc_execute_2}")
@ExaquteTask(returns=2)
def executeInstanceOnlyAdaptiveRefinementAuxLev2_Task(pickled_model,pickled_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_index,current_analysis,mlmc_results):
    try:
        open_mp_threads = int(os.environ["computing_units_mlmc_execute_2"])
        threadpool_limits(limits=open_mp_threads)
    except:
        pass
    pickled_model,mlmc_results = \
        ExecuteInstanceOnlyAdaptiveRefinementAux_Functionality(pickled_model,pickled_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_index,current_analysis,mlmc_results)
    return pickled_model,mlmc_results

@constraint(ComputingUnits="${computing_units_mlmc_execute_3}")
@ExaquteTask(returns=2)
def executeInstanceOnlyAdaptiveRefinementAuxLev3_Task(pickled_model,pickled_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_index,current_analysis,mlmc_results):
    try:
        open_mp_threads = int(os.environ["computing_units_mlmc_execute_3"])
        threadpool_limits(limits=open_mp_threads)
    except:
        pass
    pickled_model,mlmc_results = \
        ExecuteInstanceOnlyAdaptiveRefinementAux_Functionality(pickled_model,pickled_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_index,current_analysis,mlmc_results)
    return pickled_model,mlmc_results

@constraint(ComputingUnits="${computing_units_mlmc_execute_4}")
@ExaquteTask(returns=2)
def executeInstanceOnlyAdaptiveRefinementAuxLev4_Task(pickled_model,pickled_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_index,current_analysis,mlmc_results):
    try:
        open_mp_threads = int(os.environ["computing_units_mlmc_execute_4"])
        threadpool_limits(limits=open_mp_threads)
    except:
        pass
    pickled_model,mlmc_results = \
        ExecuteInstanceOnlyAdaptiveRefinementAux_Functionality(pickled_model,pickled_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_index,current_analysis,mlmc_results)
    return pickled_model,mlmc_results

@constraint(ComputingUnits="${computing_units_mlmc_execute_5}")
@ExaquteTask(returns=2)
def executeInstanceOnlyAdaptiveRefinementAuxLev5_Task(pickled_model,pickled_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_index,current_analysis,mlmc_results):
    try:
        open_mp_threads = int(os.environ["computing_units_mlmc_execute_5"])
        threadpool_limits(limits=open_mp_threads)
    except:
        pass
    pickled_model,mlmc_results = \
        ExecuteInstanceOnlyAdaptiveRefinementAux_Functionality(pickled_model,pickled_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_index,current_analysis,mlmc_results)
    return pickled_model,mlmc_results


################################### SingleAdaptiveRefinement ###################################

@constraint(ComputingUnits="${computing_units_mlmc_execute_0}")
@ExaquteTask(returns=1)
def ExecuteInstanceSingleRefinementAuxLev0_Task(pickled_model,pickled_project_parameters,random_variable,current_level,current_analysis,mlmc_results):
    try:
        open_mp_threads = int(os.environ["computing_units_mlmc_execute_0"])
        threadpool_limits(limits=open_mp_threads)
    except:
        pass
    mlmc_results = \
        ExecuteInstanceSingleRefinementAux_Functionality(pickled_model,pickled_project_parameters,random_variable,current_level,current_analysis,mlmc_results)
    return mlmc_results


@constraint(ComputingUnits="${computing_units_mlmc_execute_1}")
@ExaquteTask(returns=1)
def ExecuteInstanceSingleRefinementAuxLev1_Task(pickled_model,pickled_project_parameters,random_variable,current_level,current_analysis,mlmc_results):
    try:
        open_mp_threads = int(os.environ["computing_units_mlmc_execute_1"])
        threadpool_limits(limits=open_mp_threads)
    except:
        pass
    mlmc_results = \
        ExecuteInstanceSingleRefinementAux_Functionality(pickled_model,pickled_project_parameters,random_variable,current_level,current_analysis,mlmc_results)
    return mlmc_results

@constraint(ComputingUnits="${computing_units_mlmc_execute_2}")
@ExaquteTask(returns=1)
def ExecuteInstanceSingleRefinementAuxLev2_Task(pickled_model,pickled_project_parameters,random_variable,current_level,current_analysis,mlmc_results):
    try:
        open_mp_threads = int(os.environ["computing_units_mlmc_execute_2"])
        threadpool_limits(limits=open_mp_threads)
    except:
        pass
    mlmc_results = \
        ExecuteInstanceSingleRefinementAux_Functionality(pickled_model,pickled_project_parameters,random_variable,current_level,current_analysis,mlmc_results)
    return mlmc_results

@constraint(ComputingUnits="${computing_units_mlmc_execute_3}")
@ExaquteTask(returns=1)
def ExecuteInstanceSingleRefinementAuxLev3_Task(pickled_model,pickled_project_parameters,random_variable,current_level,current_analysis,mlmc_results):
    try:
        open_mp_threads = int(os.environ["computing_units_mlmc_execute_3"])
        threadpool_limits(limits=open_mp_threads)
    except:
        pass
    mlmc_results = \
        ExecuteInstanceSingleRefinementAux_Functionality(pickled_model,pickled_project_parameters,random_variable,current_level,current_analysis,mlmc_results)
    return mlmc_results

@constraint(ComputingUnits="${computing_units_mlmc_execute_4}")
@ExaquteTask(returns=1)
def ExecuteInstanceSingleRefinementAuxLev4_Task(pickled_model,pickled_project_parameters,random_variable,current_level,current_analysis,mlmc_results):
    try:
        open_mp_threads = int(os.environ["computing_units_mlmc_execute_4"])
        threadpool_limits(limits=open_mp_threads)
    except:
        pass
    mlmc_results = \
        ExecuteInstanceSingleRefinementAux_Functionality(pickled_model,pickled_project_parameters,random_variable,current_level,current_analysis,mlmc_results)
    return mlmc_results

@constraint(ComputingUnits="${computing_units_mlmc_execute_5}")
@ExaquteTask(returns=1)
def ExecuteInstanceSingleRefinementAuxLev5_Task(pickled_model,pickled_project_parameters,random_variable,current_level,current_analysis,mlmc_results):
    try:
        open_mp_threads = int(os.environ["computing_units_mlmc_execute_5"])
        threadpool_limits(limits=open_mp_threads)
    except:
        pass
    mlmc_results = \
        ExecuteInstanceSingleRefinementAux_Functionality(pickled_model,pickled_project_parameters,random_variable,current_level,current_analysis,mlmc_results)
    return mlmc_results


################################### Functionalities ###################################

"""
function executing an instance with concurrent adaptive refinement strategy
input:  current_MLMC_level :instance MLMC level
        pickled_coarse_model : pickled model of level (current_level-1)
        pickled_coarse_project_parameters : pickled parameters
        pickled_custom_metric_refinement_parameters : pickled metric parameters
        pickled_custom_remesh_refinement_parameters : pickled remeshing parameters
        sample                                      : random variable(s) of current instance
        current_level                               : current working level, \in [0,current_MLMC_level]
        current_analysis_stage                      : analysis stage of the problem
        mlmc_results                                : instance of MultilevelMonteCarloResults class
output: mlmc_results        : instance of MultilevelMonteCarloResults class
        pickled_finer_model : pickled model of level current_level

"""
def ExecuteInstanceConcurrentAdaptiveRefinementAux_Functionality(current_MLMC_level,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,sample,current_level,current_analysis_stage,mlmc_results):
    time_0 = time.time()
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
    time_1 = time.time()
    # start time
    start_MLMC_time = time.time()
    # refine if current current_level > 0, adaptive refinement based on the solution of previous level
    if (current_level > 0):
        time_2 = time.time()
        # unpickle metric and remesh refinement parameters and build Kratos Parameters objects
        serialized_custom_metric_refinement_parameters = pickle.loads(pickled_custom_metric_refinement_parameters)
        serialized_custom_remesh_refinement_parameters = pickle.loads(pickled_custom_remesh_refinement_parameters)
        current_custom_metric_refinement_parameters = KratosMultiphysics.Parameters()
        current_custom_remesh_refinement_parameters = KratosMultiphysics.Parameters()
        serialized_custom_metric_refinement_parameters.Load("MetricRefinementParametersSerialization",current_custom_metric_refinement_parameters)
        serialized_custom_remesh_refinement_parameters.Load("RemeshRefinementParametersSerialization",current_custom_remesh_refinement_parameters)
        del(serialized_custom_metric_refinement_parameters,serialized_custom_remesh_refinement_parameters)
        time_3 = time.time()
        # refine the model Kratos object
        adaptive_refinement_manager = AdaptiveRefinement(current_level,current_model,current_project_parameters,current_custom_metric_refinement_parameters,current_custom_remesh_refinement_parameters)
        refined_model,refined_project_parameters = adaptive_refinement_manager.ComputeAdaptiveRefinement()
        current_model = refined_model
        del(refined_project_parameters,refined_model)
        time_4 = time.time()
        time_5 = time.time()
    time_6 = time.time()
    simulation = current_analysis_stage(current_model,current_project_parameters,sample)
    # print("[CHECK OUTPUT]: sample value:", str(sample))
    sys.stdout.flush()
    simulation.Run()
    QoI = simulation.EvaluateQuantityOfInterest()
    # print("[CHECK OUTPUT]: simulation qoi:", str(QoI))
    sys.stdout.flush()
    time_7 = time.time()
    # save model and parameters as StreamSerializer Kratos objects
    serialized_finer_model = KratosMultiphysics.StreamSerializer()
    serialized_finer_model.Save("ModelSerialization",simulation.model)
    # pickle model and parameters
    pickled_finer_model = pickle.dumps(serialized_finer_model, 2) # second argument is the protocol and is NECESSARY (according to pybind11 docs)
    del(simulation)
    time_8 = time.time()
    end_MLMC_time = time.time()
    # register results of the current level in the MultilevelMonteCarloResults class
    mlmc_results.time_ML[current_level].append(end_MLMC_time-start_MLMC_time) # saving each result in the corresponding list in order to ensure the correctness of the results order and the levels
    mlmc_results.QoI[current_level].append(QoI) # saving each result in the corresponding list in order to ensure the correctness of the results order and the levels
    # post process times of the task
    # print("\n","#"*50," TIMES EXECUTE TASK ","#"*50,"\n")
    deserialization_time = time_1 - time_0
    if (current_level > 0):
        mmg_refinement_time = time_4 - time_3
    refinement_time = time_6 - time_1
    Kratos_run_time = time_7 - time_6
    serialization_time = time_8 - time_7
    total_task_time = time_8 - time_0
    # print("[LEVEL] current level:",current_level)
    # print("[TIMER] total task time:", total_task_time)
    # print("[TIMER] Kratos Run time:",Kratos_run_time)
    # print("[TIMER] Deserialization time:",deserialization_time)
    # print("[TIMER] Serialization time:",serialization_time)
    # print("[TIMER] Refinement time:",refinement_time)
    # if (current_level > 0):
    #     print("[TIMER] mmg refinement time",mmg_refinement_time)
    # print("RATIOs: time of interest / total task time")
    # print("[RATIO] Relative deserialization time:",(deserialization_time)/total_task_time)
    # print("[RATIO] Relative serialization time:",(serialization_time)/total_task_time)
    # print("[RATIO] Relative Kratos run time:",Kratos_run_time/total_task_time)
    # print("[RATIO] Relative refinement time (deserialization + mmg refinement + initialization Kratos)",refinement_time/total_task_time)
    # if (current_level > 0):
    #     print("[RATIO] Relative ONLY mmg refinement time:",mmg_refinement_time/total_task_time)
    # print("\n","#"*50," END TIMES EXECUTE TASK ","#"*50,"\n")
    sys.stdout.flush()
    return mlmc_results,pickled_finer_model

"""
function executing adaptive refinement from a given model
input:  pickled_coarse_model : pickled model of level (current_level-1)
        pickled_coarse_project_parameters : pickled parameters
        pickled_custom_metric_refinement_parameters : pickled metric parameters
        pickled_custom_remesh_refinement_parameters : pickled remeshing parameters
        random_variable                             : random variable(s) of current instance
        current_index                               : current working level
        current_analysis_stage                      : analysis stage of the problem
        mlmc_results                                : instance of MultilevelMonteCarloResults class
output: pickled_finer_model : pickled model of level current_level
        mlmc_results        : instance of MultilevelMonteCarloResults class

"""
def ExecuteInstanceOnlyAdaptiveRefinementAux_Functionality(pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_index,current_analysis_stage,mlmc_results):
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
    simulation = current_analysis_stage(current_model,current_project_parameters,random_variable)
    simulation.Initialize()
    # reset general flags
    main_model_part_name = current_project_parameters["solver_settings"]["model_part_name"].GetString()
    simulation.model.GetModelPart(main_model_part_name).ProcessInfo.SetValue(KratosMultiphysics.IS_RESTARTED,True)
    # save model and parameters as StreamSerializer Kratos objects
    serialized_finer_model = KratosMultiphysics.StreamSerializer()
    serialized_finer_model.Save("ModelSerialization",simulation.model)
    # pickle model and parameters
    pickled_finer_model = pickle.dumps(serialized_finer_model, 2) # second argument is the protocol and is NECESSARY (according to pybind11 docs)
    return pickled_finer_model,mlmc_results

"""
function executing an instance of single adaptive refinement strategy
input:  pickled_coarse_model              : pickled model of level current_level
        pickled_coarse_project_parameters : pickled parameters
        sample                            : random variable(s) of current instance
        current_level                     : current working level
        current_analysis_stage            : analysis stage of the problem
        mlmc_results                      : instance of MultilevelMonteCarloResults class
output: mlmc_results : instance of MultilevelMonteCarloResults class
"""
def ExecuteInstanceSingleRefinementAux_Functionality(pickled_coarse_model,pickled_coarse_project_parameters,sample,current_level,current_analysis_stage,mlmc_results):
    time_0 = time.time()
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
    time_1 = time.time()
    # start time
    start_MLMC_time = time.time()
    time_6 = time.time()
    simulation = current_analysis_stage(current_model,current_project_parameters,sample)
    print("[CHECK OUTPUT]: sample value:", str(sample))
    sys.stdout.flush()
    simulation.Run()
    QoI = simulation.EvaluateQuantityOfInterest()
    print("[CHECK OUTPUT]: simulation qoi:", str(QoI))
    sys.stdout.flush()
    time_7 = time.time()
    del(simulation)
    time_8 = time.time()
    end_MLMC_time = time.time()
    # register results of the current level in the MultilevelMonteCarloResults class
    mlmc_results.QoI[current_level].append(QoI) # saving each result in the corresponding list in order to ensure the correctness of the results order and the levels
    mlmc_results.time_ML[current_level].append(end_MLMC_time-start_MLMC_time) # saving each result in the corresponding list in order to ensure the correctness of the results order and the levels

    # post process times of the task
    print("\n","#"*50," TIMES EXECUTE TASK ","#"*50,"\n")
    deserialization_time = time_1 - time_0
    Kratos_run_time = time_7 - time_6
    total_task_time = time_8 - time_0
    print("[LEVEL] current level:",current_level)
    print("[TIMER] total task time:", total_task_time)
    print("[TIMER] Kratos Run time:",Kratos_run_time)
    print("[TIMER] Deserialization time:",deserialization_time)
    print("RATIOs: time of interest / total task time")
    print("[RATIO] Relative deserialization time:",(deserialization_time)/total_task_time)
    print("[RATIO] Relative Kratos run time:",Kratos_run_time/total_task_time)
    print("\n","#"*50," END TIMES EXECUTE TASK ","#"*50,"\n")

    return mlmc_results


class MultilevelMonteCarlo(object):
    """The base class for the MultilevelMonteCarlo-classes"""
    def __init__(self,custom_parameters_path,project_parameters_path,refinement_parameters_path,custom_analysis):
        """constructor of the MultilevelMonteCarlo-Object
        Keyword arguments:
        self                       : an instance of the class
        custom_parameters_path     : path of settings of the Multilevel Monte Carlo simulation
        project_parameters_path    : path of the project parameters file
        refinement_parameters_path : path of settings of the metric for the refinement and settings of the remeshing
        custom_analysis            : analysis stage of the problem
        """

        # analysis: analysis stage of the current problem
        if (custom_analysis is "auxiliary_analysis"): # needed in FinalizePhaseAux_Task to build an auxiliary MultilevelMonteCarlo class
            pass
        elif (custom_analysis is not None): # standard constructor
            self.SetAnalysis(custom_analysis)
        else:
            raise Exception ("Please provide a Kratos specific application analysis stage for the current problem")
        # project_parameters_path: path to the project parameters json file
        if (project_parameters_path is "auxiliary_project_parameters_path"): # needed in FinalizePhaseAux_Task to build an auxiliary MultilevelMonteCarlo class
            self.to_pickle_model_parameters = False
        elif(project_parameters_path is not None): # standard constructor
            self.to_pickle_model_parameters = True
            self.project_parameters_path = project_parameters_path
        else:
            raise Exception ("Please provide the path of the project parameters json file")
        # custom_settings_metric_refinement: custom settings of the metric refinement
        # custom_settings_remesh_refinement: custom settings of the remeshing
        if (refinement_parameters_path is "auxiliary_parameters_refinement_path"): # needed in FinalizePhaseAux_Task to build an auxiliary MultilevelMonteCarlo class
            self.to_pickle_custom_metric_remesh_refinement_parameters = False
            self.to_pickle_custom_metric_remesh_refinement_parameters = False
        elif (refinement_parameters_path is not None): # standard contructor
            self.to_pickle_custom_metric_remesh_refinement_parameters = True
            self.to_pickle_custom_metric_remesh_refinement_parameters = True
            self.refinement_parameters_path = refinement_parameters_path
            self.SetRefinementParameters()
        else:
            raise Exception ("Please provide the path of the refinement parameters json file")

        # default settings of the Continuation Multilevel Monte Carlo (CMLMC) algorithm
        # run_multilevel_monte_carlo: flag to run or not the algorithm
        # use_lower_levels_samples: flag to store also lower level values of QoI to compute statistics, if avaiable from adaptive refinement --> produces a BIAS in the results
        # adaptive_number_samples: if False, each level is initial_batch_size and keep constant. Add one level per iteration, if possible. If True, automatic computation of number of levels and number of samples per level
        # refinement_strategy: set how to perform refinement for MLMC algorithms. Options: concurrent_adaptive_refinement : refine each time at task level
        #                                                                                  single_refinement              : refine once in SerializeModelParameters()
        # convergence_criteria: convergence criteria exploited to check convergence
        # k0: certainty Parameter 0 rates (confidence in the variance models)
        # k1: certainty Parameter 1 rates (confidence in the weak convergence model)
        # r1: cost increase first iterations CMLMC
        # r2: cost increase final iterations CMLMC
        # initial_tolerance: tolerance iteration 0
        # tolerance: final tolerance
        # tolerance_absolute: safety tolerance (absolute). Useful if expected value (qoi) is zero
        # confidence: confidence on tolerance
        # cphi_confidence: CDF**-1 (confidence) where CDF**-1 is the inverse of the CDF of the standard normal distribution, the default value is computed for confidence = 0.9
        # initial_batch_size: number of samples for screening phase
        # initial_number_batches: number of batches for iteration 0
        # maximum_number_levels: maximum number of levels
        # maximum_number_iterations: maximum number of MLMC iterations
        # minimum_samples_add: minimum number of samples to add if at least one is going to be added
        # splitting_parameter_max : max value parameter splitting discretization and statistical errors
        # splitting_parameter_min : min value parameter splitting discretization and statistical errors
        default_settings = KratosMultiphysics.Parameters("""
        {
            "run_multilevel_monte_carlo"          : true,
            "use_lower_levels_samples"            : false,
            "adaptive_number_samples"             : false,
            "refinement_strategy"                 : "concurrent_adaptive_refinement",
            "convergence_criteria"                : "total_error_stopping_rule",
            "tasks_refinement_sample_all_at_once" : false,
            "k0"                                  : 0.1,
            "k1"                                  : 0.1,
            "r1"                                  : 1.25,
            "r2"                                  : 1.15,
            "initial_tolerance"                   : -1,
            "tolerance"                           : 0.1,
            "tolerance_absolute"                  : 1e-6,
            "confidence"                          : 0.9,
            "cphi_confidence"                     : 1.28155156554,
            "initial_batch_size"                  : [25,25,25],
            "initial_number_batches"              : 1,
            "maximum_number_levels"               : 4,
            "maximum_number_iterations"           : 5,
            "minimum_samples_add_level"           : 6.0,
            "splitting_parameter_max"             : 0.9,
            "splitting_parameter_min"             : 0.1
        }
        """)
        if (type(custom_parameters_path) == KratosMultiphysics.Parameters):
            self.settings = custom_parameters_path # needed in FinalizePhaseAux_Task to build an auxiliary MultilevelMonteCarlo class, then in FinalizePhaseAux_Task we set the true problem parameters
        else: # standard contructor
            self.custom_parameters_path = custom_parameters_path
            self.SetXMCParameters()
            # compute cphi = CDF**-1 (confidence)
            self.settings.AddEmptyValue("cphi_confidence")
            if (self.settings["confidence"].GetDouble()==1.0):
                self.settings["confidence"].SetDouble(0.999) # reduce confidence to not get +inf for cphi_confidence
            self.settings["cphi_confidence"].SetDouble(norm.ppf(self.settings["confidence"].GetDouble()))

        # validate and assign default parameters
        self.settings.ValidateAndAssignDefaults(default_settings)
        # distinguish between MLMC and CMLMC
        if (self.settings["initial_tolerance"].GetDouble() == -1): # MLMC
            self.settings["initial_tolerance"].SetDouble(self.settings["tolerance"].GetDouble()) # set initial tolerance equal to tolerance: does not change as we run more iterations
        else: # CMLMC: as long as iteration_counter grows, we decrease the tolerance from initial_tolerance to tolerance
            if (self.settings["initial_tolerance"].GetDouble() < self.settings["tolerance"].GetDouble()):
                raise Exception ("Set initial_tolerance greater or equal than tolerance")
        # warning if use_lower_levels_samples is True
        if (self.settings["use_lower_levels_samples"].GetBool()):
            print("\n ######## WARNING: use_lower_levels_samples set to True --> introducing a BIAS in the problem ########\n")
        # current_number_levels: number of levels of current iteration
        self.current_number_levels = len(self.settings["initial_batch_size"].GetVector())-1

        # batches_number_samples: total number of samples, organized in level and batches
        # batches_number_samples = [ [ [level0_batch1] [level1_batch1] .. ] [ [level0_batch2] [level1_batch2] ..] .. ]
        # for MLMC: batches_number_samples = [ [ level0_batch1 level1_batch1 .. ] [ level0_batch2 level1_batch2 .. ] [ level0_batch3 level1_batch3 .. ] .. ]
        self.batches_number_samples = []
        # number_samples: total number of samples used for the computation of global power sums
        # number_samples = [total_level0 total_level1 ..]
        self.number_samples = []
        # running_number_samples: total number of samples running
        self.running_number_samples = []
        # batches_launched: boolean true or false if batch launched or not
        self.batches_launched = []
        # batches_execution_finished: boolean true or false if batch finished or not
        self.batches_execution_finished = []
        # batches_analysis_finished: boolean true or false if statistical analysis of batch is finished or not
        self.batches_analysis_finished = []
        # batches_convergence_finished: boolean true or false if convergence computation of batch is finished or not
        self.batches_convergence_finished = []
        # # batch_size: number of iterations of each epoch
        self.batch_size = []
        # current_convergence_batch: current batch for which convergence is computed
        self.current_convergence_batch = None

        # rates_error: dictionary containing the values of the parameters
        # calpha: coefficient of the function maximizing bias
        # alpha: exponent of the function maximizing bias
        # cbeta: coefficient of the function maximizing statistical error
        # beta: exponent of the function maximizing statistical error
        # cgamma: coefficient of the function maximizing cost
        # gamma: exponent of the function maximizing cost
        self.rates_error = {"calpha":None, "alpha":None, "cbeta":None, "beta":None, "cgamma":None, "gamma":None}
        # mesh_parameters: reciprocal of minimal mesh size
        self.mesh_parameters = []
        # size_mesh: minimal mesh size
        self.mesh_sizes = []
        # bayesian_variance: Bayesian variance
        self.bayesian_variance = []
        # number_iterations: theoretical number of iterations the MLMC algorithm will perform
        self.number_iterations_iE = None
        # convergence: boolean variable defining if MLMC algorithm is converged
        self.convergence = False
        # current_iteration: current iteration of MLMC algorithm
        self.iteration_counter = 0
        # tolerance_i: tolerance of i^th-iteration of MLMC algorithm
        self.tolerance_i = None
        # theta_i: splitting parameter \in (0,1), this affects bias and statistical error in the computation of the total error
        self.theta_i = None
        # mean_mlmc_QoI: MLMC estimator for the mean value of the Quantity of Interest
        self.mean_mlmc_QoI = None
        # total_error: total error of MLMC algorithm, the sum of bias and statistical error is an overestmation of the real total error
        #               total_error := \abs(E^MLMC[QoI] - E[QoI])
        self.total_error = None

        # difference_QoI: Quantity of Interest of the considered problem organized per levels
        #                 difference_QoI.values := Y_l = QoI_M_l - Q_M_l-1
        self.difference_QoI = StatisticalVariable() # list containing Y_{l}^{i} = Q_{m_l} - Q_{m_{l-1}}
        self.difference_QoI.InitializeLists(self.settings["maximum_number_levels"].GetInt()+1,self.settings["initial_number_batches"].GetInt())
        self.difference_QoI.type = "scalar"
        # time_ML: time to perform a single MLMC simulation (i.e. one value of difference_QoI.values) organized per levels
        self.time_ML = StatisticalVariable() # list containing the time to compute the level=l simulations
        self.time_ML.InitializeLists(self.settings["maximum_number_levels"].GetInt()+1,self.settings["initial_number_batches"].GetInt())

        ########################################################################
        # observation: levels start from level 0                               #
        #              length arrays and lists starts from 1                   #
        # thus there is a difference of 1 between length lists and levels      #
        # e.g. self.current_level = len(self.number_samples) - 1               #
        #      or                                                              #
        #      self.current_level = len(self.difference_QoI.values) - 1        #
        ########################################################################

        # pickled_model: serialization of model Kratos object of the problem
        self.pickled_model = None
        self.serialized_model = None
        # pickled_project_parameters: serialization of project parameters Kratos object of the problem
        self.pickled_project_parameters = None
        self.serialized_project_parameters = None
        # pickled_custom_settings_metric_refinement: serialization of the metric refinement custom settings
        self.pickled_custom_metric_refinement_parameters = None
        self.serialized_custom_metric_refinement_parameters = None
        # pickled_custom_settings_remesh_refinement: serialization of the remesh custom settings
        self.pickled_custom_remesh_refinement_parameters = None
        self.serialized_custom_remesh_refinement_parameters = None
        # construct the pickled model and pickled project parameters of the problem
        self.is_project_parameters_pickled = False
        self.is_model_pickled = False
        # self.SerializeModelParameters()
        # construct the pickled custom settings metric refinement and the picled custom settings remesh of the problem
        self.is_custom_settings_metric_refinement_pickled = False
        self.is_custom_settings_remesh_refinement_pickled = False
        # self.SerializeRefinementParameters()

    """
    function executing the MLMC algorithm
    input:  self : an instance of the class
    """
    def Run(self):
        if (self.settings["run_multilevel_monte_carlo"].GetBool()):
            self.SerializeRefinementParameters()
            self.SerializeModelParameters()
            self.InitializeScreeningPhase()
            self.LaunchEpoch()
            self.FinalizeScreeningPhase()
            self.ScreeningInfoScreeningPhase()
            self.ComputeMeanMLMCQoI()
            while self.convergence is not True:
                self.InitializeMLMCPhase()
                self.ScreeningInfoInitializeMLMCPhase()
                self.LaunchEpoch()
                self.FinalizeMLMCPhase()
                self.ScreeningInfoFinalizeMLMCPhase()
        else:
            print("\n","#"*50,"Not running Multilevel Monte Carlo algorithm","#"*50)
            pass

    """
    function running one Multilevel Monte Carlo epoch
    input:  self : an instance of the class
    """
    def LaunchEpoch(self):
        for batch in range (len(self.batches_number_samples)):
            if (self.batches_launched[batch] is not True):
                self.batches_launched[batch] = True
                self.objects_to_delete.append([])
                # batch_results = []
                # for level in (range(len(self.batches_number_samples[batch]))):
                for level in reversed(range(len(self.batches_number_samples[batch]))):
                    batch_results = []
                    for instance in range (self.batches_number_samples[batch][level]):
                        self.running_number_samples[level] = self.running_number_samples[level] + 1
                        if (self.settings["refinement_strategy"].GetString() == "concurrent_adaptive_refinement"):
                            batch_results.append(self.ExecuteInstanceConcurrentAdaptiveRefinement(level))
                        elif (self.settings["refinement_strategy"].GetString() == "single_refinement"):
                            batch_results.append(self.ExecuteInstanceSingleRefinement(level))
                        else:
                            raise Exception ("Specify refinement_strategy: concurrent_adaptive_refinement or single_refinement")
                        self.running_number_samples[level] = self.running_number_samples[level] - 1
                    self.AddResults(batch_results,batch)
                self.batches_execution_finished[batch] = True

    """
    function executing an instance of the MLMC algorithm, i.e. a single MC simulation and eventually the refinement (that occurs before the simulation run), with concurrent adaptive refinement strategy
    input:  self : an instance of the class
    output: mlmc_results       : instance of the MultilevelMonteCarloResult class
            current_MLMC_level : level of the current MLMC simulation
    """
    def ExecuteInstanceConcurrentAdaptiveRefinement(self,current_MLMC_level):
        # local variables
        pickled_coarse_model = self.pickled_model[0] # first component is relative to original model part
        pickled_coarse_project_parameters = self.pickled_project_parameters[0] # first component is relative to original project parameters
        pickled_custom_metric_refinement_parameters = self.pickled_custom_metric_refinement_parameters
        pickled_custom_remesh_refinement_parameters = self.pickled_custom_remesh_refinement_parameters
        current_analysis = self.analysis
        different_tasks = not self.settings["tasks_refinement_sample_all_at_once"].GetBool()
        # generate the sample
        sample = generator.GenerateSample(self.problem_name)
        # initialize the MultilevelMonteCarloResults class
        mlmc_results = MultilevelMonteCarloResults(current_MLMC_level)
        if (current_MLMC_level == 0):
            current_level = 0
            mlmc_results,pickled_current_model = \
                ExecuteInstanceConcurrentAdaptiveRefinementMultipleTasks_Wrapper(current_MLMC_level,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,sample,current_level,current_analysis,mlmc_results)
                # ExecuteInstanceConcurrentAdaptiveRefinementAux_Task(current_MLMC_level,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,sample,current_level,current_analysis,mlmc_results)
        else:
            if (different_tasks is True):
                for current_level in range(current_MLMC_level+1):
                    mlmc_results,pickled_current_model = \
                        ExecuteInstanceConcurrentAdaptiveRefinementMultipleTasks_Wrapper(current_MLMC_level,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,sample,current_level,current_analysis,mlmc_results)
                        # ExecuteInstanceConcurrentAdaptiveRefinementAux_Task(current_MLMC_level,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,sample,current_level,current_analysis,mlmc_results)
                    if (current_level > 0):
                        self.objects_to_delete[-1].append(pickled_coarse_model)
                    if (current_level > 0 and current_level == current_MLMC_level):
                        self.objects_to_delete[-1].append(pickled_current_model)
                    del(pickled_coarse_model)
                    pickled_coarse_model = pickled_current_model
                    del(pickled_current_model)
            elif (different_tasks is False):
                mlmc_results = \
                    ExecuteInstanceConcurrentAdaptiveRefinementAllAtOnce_Wrapper(current_MLMC_level,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,sample,current_analysis,mlmc_results)
            else:
                raise Exception ("boolean variable different task is not a boolean, instead is equal to",different_tasks)
        return mlmc_results,current_MLMC_level

    """
    function executing an instance of the MLMC algorithm, i.e. a single MC simulation and eventually the refinement (that occurs before the simulation run), with single adaptive refinement strategy
    input:  self : an instance of the class
    output: mlmc_results       : instance of the MultilevelMonteCarloResult class
            current_MLMC_level : level of the current MLMC simulation
    """
    def ExecuteInstanceSingleRefinement(self,current_MLMC_level):
        current_analysis = self.analysis
        # generate the sample
        sample = generator.GenerateSample(self.problem_name)
        # initialize the MultilevelMonteCarloResults class
        mlmc_results = MultilevelMonteCarloResults(current_MLMC_level)
        if (current_MLMC_level == 0):
            current_level = 0
            mlmc_results = \
                ExecuteInstanceSingleRefinement_Wrapper(current_MLMC_level,self.pickled_model[current_level],self.pickled_project_parameters[0],sample,current_level,current_analysis,mlmc_results)
        else:
            for current_level in range(current_MLMC_level-1,current_MLMC_level+1):
                mlmc_results = \
                    ExecuteInstanceSingleRefinement_Wrapper(current_MLMC_level,self.pickled_model[current_level],self.pickled_project_parameters[0],sample,current_level,current_analysis,mlmc_results)
        mlmc_results = get_value_from_remote(mlmc_results)
        return mlmc_results,current_MLMC_level

    """
    function initializing the screening phase (MLMC phase iteration zero)
    input:  self : an instance of the class
    """
    def InitializeScreeningPhase(self):
        if (self.iteration_counter == 0):
            self.batch_size = [self.settings["initial_batch_size"][level].GetInt() for level in range (self.current_number_levels+1)]
            self.batches_number_samples = [[self.settings["initial_batch_size"][level].GetInt() for level in range (self.current_number_levels+1)] for _ in range (self.settings["initial_number_batches"].GetInt())]
            self.number_samples = [0 for _ in range (self.settings["maximum_number_levels"].GetInt()+1)]
            self.running_number_samples = [0 for _ in range (self.settings["maximum_number_levels"].GetInt()+1)]
            self.batches_launched = [False for _ in range (self.settings["initial_number_batches"].GetInt())]
            self.batches_execution_finished = [False for _ in range (self.settings["initial_number_batches"].GetInt())]
            self.batches_analysis_finished = [False for _ in range (self.settings["initial_number_batches"].GetInt())]
            self.batches_convergence_finished = [False for _ in range (self.settings["initial_number_batches"].GetInt())]
            self.objects_to_delete = []
            # estimate the mesh sizes
            self.ComputeMeshParameters()

    """
    function serializing and pickling the model and the project parameters of the problem
    the serialization-pickling process is the following:
    i)   from Model/Parameters Kratos object to StreamSerializer Kratos object
    ii)  from StreamSerializer Kratos object to pickle string
    iii) from pickle string to StreamSerializer Kratos object
    iv)  from StreamSerializer Kratos object to Model/Parameters Kratos object
    requires: self.project_parameters_path : path of the Project Parameters file
    builds: self.pickled_model              : pickled model
            self.pickled_project_parameters : pickled project parameters
    input:  self : an instance of the class
    """
    def SerializeModelParameters(self):
        if (self.to_pickle_model_parameters):
            self.serialized_model = []
            self.serialized_project_parameters = []
            self.pickled_model = []
            self.pickled_project_parameters = []

            if (self.settings["refinement_strategy"].GetString() == "concurrent_adaptive_refinement"):
                self.SerializeModelParametersConcurrentAdaptiveRefinement()
            elif (self.settings["refinement_strategy"].GetString() == "single_refinement"):
                self.SerializeModelParametersSingleRefinement()
            else:
                raise Exception ("Specify refinement_strategy: concurrent_adaptive_refinement or single_refinement")
            self.is_project_parameters_pickled = True
            self.is_model_pickled = True
            print("\n","#"*50," SERIALIZATION MODEL AND PROJECT PARAMETERS COMPLETED ","#"*50,"\n")

    def SerializeModelParametersConcurrentAdaptiveRefinement(self):
        with open(self.project_parameters_path,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        # create wrapper instance to modify current project parameters
        self.wrapper = ParametersWrapper(parameters)
        # save problem name
        self.problem_name = parameters["problem_data"]["problem_name"].GetString()
        # serialize parmeters (to avoid adding new data dependent on the application)
        parameters = self.wrapper.SetModelImportSettingsInputType("use_input_model_part")
        serialized_project_parameters = KratosMultiphysics.StreamSerializer()
        serialized_project_parameters.Save("ParametersSerialization",parameters)
        self.serialized_project_parameters.append(serialized_project_parameters)
        # reset to read the model part
        parameters = self.wrapper.SetModelImportSettingsInputType("mdpa")
        # prepare the model to serialize
        model = KratosMultiphysics.Model()
        fake_sample = generator.GenerateSample(self.problem_name) # only used to serialize
        simulation = self.analysis(model,parameters,fake_sample)
        simulation.Initialize()
        # reset general flags
        main_model_part_name = self.wrapper.GetModelPartName()
        simulation.model.GetModelPart(main_model_part_name).ProcessInfo.SetValue(KratosMultiphysics.IS_RESTARTED,True)
        # serialize model
        serialized_model = KratosMultiphysics.StreamSerializer()
        serialized_model.Save("ModelSerialization",simulation.model)
        self.serialized_model.append(serialized_model)
        # pickle dataserialized_data
        pickled_model = pickle.dumps(serialized_model, 2) # second argument is the protocol and is NECESSARY (according to pybind11 docs)
        pickled_project_parameters = pickle.dumps(serialized_project_parameters, 2)
        self.pickled_model.append(pickled_model)
        self.pickled_project_parameters.append(pickled_project_parameters)

    def SerializeModelParametersSingleRefinement(self):
        self.SerializeModelParametersConcurrentAdaptiveRefinement() # to prepare parameters and model part of coarsest level
        number_levels_to_serialize = self.settings["maximum_number_levels"].GetInt()
        # same routine of ExecuteInstanceConcurrentAdaptiveRefinemnt() to build models and parameters, but here we save models and parameters
        pickled_coarse_model = self.pickled_model[0]
        pickled_coarse_project_parameters = self.pickled_project_parameters[0]
        pickled_custom_metric_refinement_parameters = self.pickled_custom_metric_refinement_parameters
        pickled_custom_remesh_refinement_parameters = self.pickled_custom_remesh_refinement_parameters
        current_analysis = self.analysis
        # generate the sample and prepare mlmc result class
        fake_sample = generator.GenerateSample(self.problem_name)
        fake_mlmc_results = MultilevelMonteCarloResults(number_levels_to_serialize)
        for serialize_up_to_level in range (number_levels_to_serialize+1):
            if (serialize_up_to_level > 0):
                for current_level in range(serialize_up_to_level+1):
                    if (current_level < serialize_up_to_level):
                        fake_mlmc_results,pickled_current_model = \
                            ExecuteInstanceConcurrentAdaptiveRefinementMultipleTasks_Wrapper(serialize_up_to_level,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,fake_sample,current_level,current_analysis,fake_mlmc_results)
                            # ExecuteInstanceConcurrentAdaptiveRefinementAux_Task(serialize_up_to_level,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,fake_sample,current_level,current_analysis,fake_mlmc_results)
                        del(pickled_coarse_model)
                        pickled_coarse_model = pickled_current_model
                    elif (current_level == serialize_up_to_level):
                        fake_mlmc_results,pickled_current_model = ExecuteInstanceOnlyAdaptiveRefinement_Wrapper(pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,fake_sample,current_level,current_analysis,fake_mlmc_results)
                    # save if current level > 0
                    if (current_level == serialize_up_to_level and current_level > 0):
                        # save pickled and serialized model and parameters
                        self.pickled_model.append(pickled_current_model)
                        # use next line only if not using pickle
                        # self.serialized_model.append(pickle.loads(pickled_current_model))
                    del(pickled_current_model)

    """
    function reading the XMC parameters passed from json file
    input:  self : an instance of the class
    """
    def SetXMCParameters(self):
        with open(self.custom_parameters_path,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        self.settings = parameters["multilevel_monte_carlo"]

    """
    function reading the refinement parameters passed from json file
    input:  self : an instance of the class
    """
    def SetRefinementParameters(self):
        with open(self.refinement_parameters_path,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        self.custom_metric_refinement_parameters = parameters["hessian_metric"]
        self.custom_remesh_refinement_parameters = parameters["refinement_mmg"]

    """
    function serializing and pickling the custom setting metric and remeshing for the refinement
    requires: self.custom_metric_refinement_parameters
              self.custom_remesh_refinement_parameters
    builds: self.pickled_custom_metric_refinement_parameters
            self.pickled_custom_remesh_refinement_parameters
    input:  self : an instance of the class
    """
    def SerializeRefinementParameters(self):
        if (self.to_pickle_custom_metric_remesh_refinement_parameters):
            metric_refinement_parameters = self.custom_metric_refinement_parameters
            remeshing_refinement_parameters = self.custom_remesh_refinement_parameters
            # save parameters as StreamSerializer Kratos objects
            serialized_metric_refinement_parameters = KratosMultiphysics.StreamSerializer()
            serialized_metric_refinement_parameters.Save("MetricRefinementParametersSerialization",metric_refinement_parameters)
            serialized_remesh_refinement_parameters = KratosMultiphysics.StreamSerializer()
            serialized_remesh_refinement_parameters.Save("RemeshRefinementParametersSerialization",remeshing_refinement_parameters)
            # pickle parameters
            pickled_metric_refinement_parameters = pickle.dumps(serialized_metric_refinement_parameters, 2) # second argument is the protocol and is NECESSARY (according to pybind11 docs)
            pickled_remesh_refinement_parameters = pickle.dumps(serialized_remesh_refinement_parameters, 2) # second argument is the protocol and is NECESSARY (according to pybind11 docs)
            self.pickled_custom_metric_refinement_parameters = pickled_metric_refinement_parameters
            self.pickled_custom_remesh_refinement_parameters = pickled_remesh_refinement_parameters
            print("\n","#"*50," SERIALIZATION REFINEMENT PARAMETERS COMPLETED ","#"*50,"\n")

    """
    function setting the Kratos analysis stage
    input:  self                       : an instance of the class
            application_analysis_stage : analysis stage of the problem
    """
    def SetAnalysis(self,application_analysis_stage):
        self.analysis = application_analysis_stage

    """
    function returning the Kratos analysis stage of the problem
    input:  self : an instance of the class
    output: self.analysis : analysis stage of the problem
    """
    def GetAnalysis(self):
        if (self.analysis is not None):
            return self.analysis
        else:
            print("Please provide a Kratos specific application analysis stage for the current problem.")

    """
    function finalizing the screening phase of MLMC algorithm
    usage: designed to be called ONCE, AFTER the screening phase
    input:  self : an instance of the class
    """
    def FinalizeScreeningPhase(self):
        # compute number of iterations of CMLMC algorithm
        self.ComputeNumberIterationsMLMC()
        # update power sums batches
        for batch in range (len(self.batches_number_samples)): # i.e. total number of batches
            if (self.batches_execution_finished[batch] is True and self.batches_analysis_finished[batch] is not True): # consider batches completed and not already analysed
                for level in range (len(self.batches_number_samples[batch])):
                    self.difference_QoI.UpdateBatchesPassPowerSum(level,batch)
                    self.time_ML.UpdateBatchesPassPowerSum(level,batch)
                self.batches_analysis_finished[batch] = True
        continue_iterating = True
        for batch in range (len(self.batches_number_samples)):
            if (self.batches_execution_finished[batch] is True and self.batches_analysis_finished[batch] is True and self.batches_convergence_finished[batch] is not True and continue_iterating): # consider batches completed, analysed and for which convergence has not been computed
                continue_iterating = False
                # update working convergence batch
                self.current_convergence_batch = batch
                for level in range (len(self.batches_number_samples[batch])):
                    # update global power sums from batches power sums
                    self.difference_QoI.UpdateGlobalPowerSums(level,batch)
                    self.time_ML.UpdateGlobalPowerSums(level,batch)
                    # update number of samples used to compute global power sums
                    self.number_samples[level] = self.number_samples[level] + self.batches_number_samples[batch][level]
                    self.difference_QoI.ComputeHStatistics(level)
                    self.time_ML.ComputeHStatistics(level)
                self.batches_convergence_finished[batch] = True
                """
                the functions executed in FinalizePhaseAux_Task are the followings:
                self.ComputeRatesLS()
                self.EstimateBayesianVariance(self.current_number_levels)
                """
                # store lists in synchro_element to execute FinalizePhaseAux_Task
                synchro_elements = [x for x in self.difference_QoI.h_statistics_1]
                synchro_elements.extend(["%%%"])
                synchro_elements.extend(self.difference_QoI.h_statistics_2)
                synchro_elements.extend(["%%%"])
                synchro_elements.extend(self.time_ML.h_statistics_1)
                # create a StreamSerializer Kratos object containing the problem settings
                serial_settings = KratosMultiphysics.StreamSerializer()
                serial_settings.Save("ParametersSerialization", self.settings)
                # compute remaining MLMC finalize process operations
                # observation: we are passing self.settings and we will exploit it to construct the class
                self.rates_error,self.bayesian_variance,self.mean_mlmc_QoI,\
                self.total_error,self.number_samples\
                = FinalizePhaseAux_Task(MultilevelMonteCarlo,
                serial_settings,self.mesh_parameters,self.current_number_levels,\
                self.iteration_counter,self.number_samples,*synchro_elements)
                # synchronization point needed to compute the other functions of MLMC algorithm
                # put as in the end as possible the synchronization point
                self.difference_QoI.h_statistics_1 = get_value_from_remote(self.difference_QoI.h_statistics_1)
                self.difference_QoI.h_statistics_2 = get_value_from_remote(self.difference_QoI.h_statistics_2)
                self.time_ML.h_statistics_1 = get_value_from_remote(self.time_ML.h_statistics_1)
                self.total_error = get_value_from_remote(self.total_error)
                self.rates_error = get_value_from_remote(self.rates_error)
                self.bayesian_variance = get_value_from_remote(self.bayesian_variance)
                self.mean_mlmc_QoI = get_value_from_remote(self.mean_mlmc_QoI)
                self.number_samples = get_value_from_remote(self.number_samples)
                # TODO: this is a workaround of a compss bug
                # delete models and parameters no more used
                for elem in self.objects_to_delete[batch]:
                    delete_object(elem)
                # update iteration counter
                self.iteration_counter = self.iteration_counter + 1
                break
        if (self.iteration_counter >= self.settings["maximum_number_iterations"].GetInt()):
            self.convergence = True

    """
    function performing all the required operations BEFORE the MLMC solution step
    input:  self : an instance of the class
    """
    def InitializeMLMCPhase(self):
        # add new batches in case convergence = False
        if (self.convergence is not True):
            # compute tolerance
            self.ComputeTolerancei()
            # estimate batches to append and batch size
            self.UpdateBatches()
        for batch in range (len(self.batches_launched)):
            if (self.batches_launched[batch] is False):
                self.batches_number_samples.append(self.batch_size)
                # append execution, analysis and convergence False booleans
                self.batches_execution_finished.append(False)
                self.batches_analysis_finished.append(False)
                self.batches_convergence_finished.append(False)
                # append new batch lists to self.QoI.values and self.QoI.power_sum_batches_*
                maximum_number_levels = self.settings["maximum_number_levels"].GetInt()+1
                self.difference_QoI.values.append([[] for _ in range (maximum_number_levels+1)])
                self.difference_QoI.power_sum_batches_1.append([[] for _ in range (maximum_number_levels+1)])
                self.difference_QoI.power_sum_batches_2.append([[] for _ in range (maximum_number_levels+1)])
                self.difference_QoI.power_sum_batches_3.append([[] for _ in range (maximum_number_levels+1)])
                self.difference_QoI.power_sum_batches_4.append([[] for _ in range (maximum_number_levels+1)])
                self.difference_QoI.batches_number_samples.append([0 for _ in range (maximum_number_levels+1)])
                self.time_ML.values.append([[] for _ in range (maximum_number_levels+1)])
                self.time_ML.power_sum_batches_1.append([[] for _ in range (maximum_number_levels+1)])
                self.time_ML.power_sum_batches_2.append([[] for _ in range (maximum_number_levels+1)])
                self.time_ML.power_sum_batches_3.append([[] for _ in range (maximum_number_levels+1)])
                self.time_ML.power_sum_batches_4.append([[] for _ in range (maximum_number_levels+1)])
                self.time_ML.batches_number_samples.append([0 for _ in range (maximum_number_levels+1)])

    """
    function updating number of batches and batch size
    input:  self : an instance of the class
    """
    def UpdateBatches(self):
        # set here number of batches to append
        if (len(self.batches_number_samples) >= self.settings["maximum_number_iterations"].GetInt()):
            new_number_batches = 0
        else:
            new_number_batches = 1
        if (self.settings["adaptive_number_samples"].GetBool() is True):
            # compute optimal number of levels
            self.ComputeLevels()
            # compute theta splitting parameter according to the current_number_levels and tolerance_i
            self.ComputeTheta(self.current_number_levels)
            # compute number of samples according to bayesian_variance and theta_i parameters
            # TODO: rename UpdateBatchSize this function
            self.ComputeNumberSamples()
            self.batch_size = [min(self.batch_size[level],2*self.number_samples[level]) for level in range (self.current_number_levels+1)]
        elif (self.settings["adaptive_number_samples"].GetBool() is False):
            #  add one level per time
            self.current_number_levels = min(self.current_number_levels+1,self.settings["maximum_number_levels"].GetInt())
            # compute theta splitting parameter according to the current_number_levels and tolerance_i
            self.ComputeTheta(self.current_number_levels)
            # add batch size per level: multiply by 1 current number samples and take the maximum between this integer and default batch size
            # self.batch_size = [max(self.number_samples[level],self.settings["initial_batch_size"][level].GetInt()) for level in range (self.current_number_levels+1)]
            self.batch_size = [self.settings["initial_batch_size"][level].GetInt() for level in range (self.current_number_levels+1)]
        else:
            raise Exception ("Assign True or False to adaptive_number_samples setting.")
        for new_batch in range (new_number_batches):
            self.batches_launched.append(False)

    """
    function performing all the required operations AFTER the MLMC solution step
    input:  self : an instance of the class
    """
    def FinalizeMLMCPhase(self):
        # update power sums batches
        for batch in range (len(self.batches_number_samples)): # i.e. total number of batches
            if (self.batches_execution_finished[batch] is True and self.batches_analysis_finished[batch] is not True): # consider batches completed and not already analysed
                for level in range (len(self.batches_number_samples[batch])):
                    self.difference_QoI.UpdateBatchesPassPowerSum(level,batch)
                    self.time_ML.UpdateBatchesPassPowerSum(level,batch)
                self.batches_analysis_finished[batch] = True
        continue_iterating = True
        for batch in range (len(self.batches_number_samples)):
            if (self.batches_execution_finished[batch] is True and self.batches_analysis_finished[batch] is True and self.batches_convergence_finished[batch] is not True and continue_iterating): # consider batches completed, analysed and for which convergence has not been computed
                continue_iterating = False
                # update working convergence batch
                self.current_convergence_batch = batch
                for level in range (len(self.batches_number_samples[batch])):
                    # update global power sums from batches power sums
                    self.difference_QoI.UpdateGlobalPowerSums(level,batch)
                    self.time_ML.UpdateGlobalPowerSums(level,batch)
                    # update number of samples used to compute global power sums
                    self.number_samples[level] = self.number_samples[level] + self.batches_number_samples[batch][level]
                    self.difference_QoI.ComputeHStatistics(level)
                    self.time_ML.ComputeHStatistics(level)
                self.batches_convergence_finished[batch] = True
                """
                the functions executed in FinalizePhaseAux_Task are the followings:
                self.ComputeMeanMLMCQoI()
                self.ComputeRatesLS()
                self.EstimateBayesianVariance(self.current_number_levels)
                self.ComputeTotalErrorMLMC()
                """
                # store lists in synchro_element to execute FinalizePhaseAux_Task
                synchro_elements = [x for x in self.difference_QoI.h_statistics_1]
                synchro_elements.extend(["%%%"])
                synchro_elements.extend(self.difference_QoI.h_statistics_2)
                synchro_elements.extend(["%%%"])
                synchro_elements.extend(self.time_ML.h_statistics_1)
                # create a StreamSerializer Kratos object containing the problem settings
                serial_settings = KratosMultiphysics.StreamSerializer()
                serial_settings.Save("ParametersSerialization", self.settings)
                # compute remaining MLMC finalize process operations
                # observation: we are passing self.settings and we will exploit it to construct the class
                self.rates_error,self.bayesian_variance,self.mean_mlmc_QoI,\
                self.total_error,self.number_samples\
                = FinalizePhaseAux_Task(MultilevelMonteCarlo,
                serial_settings,self.mesh_parameters,self.current_number_levels,\
                self.iteration_counter,self.number_samples,*synchro_elements)
                # synchronization point needed to compute the other functions of MLMC algorithm
                # put as in the end as possible the synchronization point
                self.difference_QoI.h_statistics_1 = get_value_from_remote(self.difference_QoI.h_statistics_1)
                self.difference_QoI.h_statistics_2 = get_value_from_remote(self.difference_QoI.h_statistics_2)
                self.time_ML.h_statistics_1 = get_value_from_remote(self.time_ML.h_statistics_1)
                self.total_error = get_value_from_remote(self.total_error)
                self.rates_error = get_value_from_remote(self.rates_error)
                self.bayesian_variance = get_value_from_remote(self.bayesian_variance)
                self.mean_mlmc_QoI = get_value_from_remote(self.mean_mlmc_QoI)
                self.number_samples = get_value_from_remote(self.number_samples)
                # check convergence
                # convergence reached if: total_error < tolerance
                # if not update current_iteration
                if (self.settings["convergence_criteria"].GetString() == "total_error_stopping_rule"):
                    if (self.total_error < self.settings["tolerance"].GetDouble()):
                        self.convergence = True
                elif (self.settings["convergence_criteria"].GetString() == "relative_total_error_stopping_rule"):
                    factor =  np.abs(self.settings["tolerance"].GetDouble() * self.mean_mlmc_QoI) + self.settings["tolerance_absolute"].GetDouble()
                    if (self.total_error < np.abs(factor)):
                        self.convergence = True
                # TODO: this is a workaround of a compss bug
                # delete models and parameters no more used
                for elem in self.objects_to_delete[batch]:
                    delete_object(elem)
                # update iteration counter
                self.iteration_counter = self.iteration_counter + 1
                break
        if (self.iteration_counter >= self.settings["maximum_number_iterations"].GetInt()):
            self.convergence = True


    """
    function printing informations about screening phase
    input:  self : an instance of the class
    """
    def ScreeningInfoScreeningPhase(self):
        print("\n","#"*50," SCREENING PHASE ","#"*50,"\n")
        print("current convergence batch =",self.current_convergence_batch)
        # print("values computed of QoI = ",self.difference_QoI.values)
        # print("values computed time_ML",self.time_ML.values)
        print("current batches",self.batches_number_samples)
        print("current number of samples = ",self.number_samples)
        print("mean and variance difference_QoI = ",self.difference_QoI.h_statistics_1,self.difference_QoI.h_statistics_2)
        print("mean time_ML",self.time_ML.h_statistics_1)
        print("rates coefficient = ",self.rates_error)
        print("estimated Bayesian variance = ",self.bayesian_variance)
        print("minimum number of MLMC iterations = ",self.number_iterations_iE)

    """
    function printing informations about initializing MLMC phase
    input:  self : an instance of the class
    """
    def ScreeningInfoInitializeMLMCPhase(self):
        print("\n","#"*50," MLMC iter =  ",self.iteration_counter,"#"*50,"\n")
        print("current tolerance = ",self.tolerance_i)
        print("updated estimated Bayesian Variance initialize phase = ",self.bayesian_variance)
        print("current number of levels = ",self.current_number_levels)
        print("current splitting parameter = ",self.theta_i)
        print("batch size = ",self.batch_size)
        print("batches prepared = ",self.batches_number_samples)

    """
    function printing informations about finalizing MLMC phase
    input:  self : an instance of the class
    """
    def ScreeningInfoFinalizeMLMCPhase(self):
        # print("values computed of QoI = ",self.difference_QoI.values)
        # print("values computed time_ML",self.time_ML.values)
        print("current convergence batch =",self.current_convergence_batch)
        print("current batches",self.batches_number_samples)
        print("check number samples of batch statistical variable class",self.difference_QoI.batches_number_samples)
        print("current number of samples",self.number_samples)
        print("mean and variance difference_QoI = ",self.difference_QoI.h_statistics_1,self.difference_QoI.h_statistics_2)
        print("mean time_ML",self.time_ML.h_statistics_1)
        print("rates coefficient = ",self.rates_error)
        print("estimated Bayesian variance = ",self.bayesian_variance)
        print("multilevel monte carlo mean estimator = ",self.mean_mlmc_QoI)
        print("total_error = bias + statistical error = ",self.total_error)
        print("relative total error = total_error / MLMC expected value = ",self.total_error/self.mean_mlmc_QoI)
        import sys
        sys.stdout.flush()

    """
    function adding QoI and MLMC time values to the corresponding level
    input:  self : an instance of the class
            simulation_results : tuple=(instance of MultilevelMonteCarloResults class, working level)
            batch_number       : number of working batch
            batch_size         : compute add result for with this size
    """
    def AddResults(self,simulation_results,batch_number,mini_batch_size=30):
        # store MLMC level of working batch in a list
        simulation_levels = list(map(lambda x: x[1], simulation_results))
        # check all the levels are equal (we are adding results level by level)
        if (all(lev==simulation_levels[0] for lev in simulation_levels)):
            current_level = simulation_levels[0]
        else:
            raise Exception ("All levels should be the same, since the hierarchy is: values > big batch > levels > small batch")
        # store MultilevelMonteCarloResult class of working branch in a list
        simulation_results = list(map(lambda x: x[0], simulation_results))
        # check if storing lower levels samples or not and prepare for loop
        if (self.settings["use_lower_levels_samples"].GetBool()):
            start_level = 0
        else:
            start_level = np.maximum(0,current_level)
        number_samples_batches_level = [0 for _ in range (current_level+1)]
        while (len(simulation_results) >= 1):
            new_simulations_results = simulation_results[mini_batch_size:]
            current_simulations = simulation_results[:mini_batch_size]
            for lev in range (start_level,current_level+1):
                number_samples_batches_level[lev] += len(current_simulations)
                # call to add results task
                # TODO: temporarily use this name, but it is already the power sum!!
                difference_QoI_list,time_ML_list = AddResultsAux_Task(lev,*current_simulations)
                self.difference_QoI.values[batch_number][lev].append(difference_QoI_list)
                self.time_ML.values[batch_number][lev].append(time_ML_list)
                # update number of samples of the batch if needed
                if (lev != current_level):
                    self.batches_number_samples[batch_number][lev] = self.batches_number_samples[batch_number][lev] + len(current_simulations)
            # prepare new list
            simulation_results = new_simulations_results
        for lev in range (start_level,current_level+1):
            self.difference_QoI.batches_number_samples[batch_number][lev] = number_samples_batches_level[lev]
            self.time_ML.batches_number_samples[batch_number][lev] = number_samples_batches_level[lev]

    """
    function giving as output the mesh discretization parameter
    the mesh parameter is the reciprocal of the minimum mesh size of the grid
    input:  self : an instance of the class
    """
    def ComputeMeshParameters(self):
        # unpickle and unserialize model and build Kratos Model object
        serialized_model = pickle.loads(self.pickled_model[0])
        current_model = KratosMultiphysics.Model()
        serialized_model.Load("ModelSerialization",current_model)
        # unpickle and unserialize parameters and build Kratos Parameters object
        serialized_project_parameters = pickle.loads(self.pickled_project_parameters[0])
        current_project_parameters = KratosMultiphysics.Parameters()
        serialized_project_parameters.Load("ParametersSerialization",current_project_parameters)
        # unpickle and unserialize metric refinement parameters and build Kratos Parameters objects
        serialized_custom_metric_refinement_parameters = pickle.loads(self.pickled_custom_metric_refinement_parameters)
        current_custom_metric_refinement_parameters = KratosMultiphysics.Parameters()
        serialized_custom_metric_refinement_parameters.Load("MetricRefinementParametersSerialization",current_custom_metric_refinement_parameters)
        for level in range(0,self.settings["maximum_number_levels"].GetInt()+1):
            adaptive_refinement_manager = AdaptiveRefinement(level,current_model,current_project_parameters,current_custom_metric_refinement_parameters,None)
            adaptive_refinement_manager.EstimateMeshSizeCurrentLevel()
            h_current_level = adaptive_refinement_manager.mesh_size
            mesh_parameter_current_level = h_current_level**(-1)
            self.mesh_sizes.append(h_current_level)
            self.mesh_parameters.append(mesh_parameter_current_level)

    """
    function computing the problem parameters P=[calpha,alpha,cbeta,beta,cgamma,gamma] using least squares fit
    we consider level > 0 to compute calpha,alpha,cbeta,beta for robustness reasons [see [1] for details]
    input:  self : an instance of the class
    """
    def ComputeRatesLS(self):
        if (self.difference_QoI.h_statistics_1[-1] == []):
            current_number_levels = self.difference_QoI.h_statistics_1.index([]) # index() method finds the given element in a list and returns its position (starts from 1, not from 0)
        else:
            current_number_levels = self.settings["maximum_number_levels"].GetInt()+1
        bias_ratesLS = np.abs(self.difference_QoI.h_statistics_1[:current_number_levels])
        variance_ratesLS = self.difference_QoI.h_statistics_2[:current_number_levels]
        cost_ML_ratesLS = self.time_ML.h_statistics_1[:current_number_levels]
        mesh_param_ratesLS = self.mesh_parameters[:current_number_levels]
        # mean - alpha
        # linear fit
        pa = np.polyfit(np.log2(mesh_param_ratesLS[1::]),np.log2(bias_ratesLS[1::]),1)
        alpha   = -pa[0]
        C1      = 2**pa[1]
        # variance - beta
        # linear fit
        pb          = np.polyfit(np.log2(mesh_param_ratesLS[1::]),np.log2(variance_ratesLS[1::]),1)
        beta        = -pb[0]
        C2          = 2**pb[1]
        # cost of computation - gamma
        # linear fit
        pg          = np.polyfit(np.log2(mesh_param_ratesLS),np.log2(cost_ML_ratesLS),1)
        gamma       = pg[0]
        C3          = 2**pg[1]
        # update the rates error dictionary
        self.rates_error["calpha"] = C1
        self.rates_error["alpha"] = alpha
        self.rates_error["cbeta"] = C2
        self.rates_error["beta"] = beta
        self.rates_error["cgamma"] = C3
        self.rates_error["gamma"] = gamma
        del(bias_ratesLS,variance_ratesLS,cost_ML_ratesLS,mesh_param_ratesLS,C1,C2,C3,alpha,beta,gamma)

    """
    function performing the Bayesian update of the variance
    using samples generated on all levels in order to locally improve the estimation of variance[difference_QoI]
    input:  self   : an instance of the class
            levels : number of levels for which we want to estimate the Bayesian variance
    """
    def EstimateBayesianVariance(self,levels): # need to keep levels because in ComputeLevels I use the maximum number of levels
        # use local variables
        k0 = self.settings["k0"].GetDouble()
        k1 = self.settings["k1"].GetDouble()
        calfa = self.rates_error["calpha"]
        alfa  = self.rates_error["alpha"]
        cbeta = self.rates_error["cbeta"]
        beta  = self.rates_error["beta"]
        mesh_param = self.mesh_parameters
        # use local variables, in order to not modify the global variables
        if (self.difference_QoI.h_statistics_1[-1] == []):
            current_number_levels = self.difference_QoI.h_statistics_1.index([])
        else:
            current_number_levels = self.settings["maximum_number_levels"].GetInt()+1
        mean_local = copy.copy(self.difference_QoI.h_statistics_1[:current_number_levels])
        variance_local = copy.copy(self.difference_QoI.h_statistics_2[:current_number_levels])
        nsam_local = copy.copy(self.number_samples)
        # append null values to evaluate the Bayesian variance for all levels
        if len(mean_local) < (levels+1):
            for _ in range (0,(levels+1)-len(mean_local)):
                mean_local.append(0.0)
        if len(variance_local) < (levels+1):
            for _ in range (0,(levels+1)-len(variance_local)):
                variance_local.append(0.0)
        if len(nsam_local) < (levels+1):
            for _ in range ((levels+1)-len(nsam_local)):
                nsam_local.append(0)
        # estimate the Bayesian variance
        bayesian_variance = []
        for level in range (0, (levels+1)):
            mu = calfa*mesh_param[level]**(-alfa)
            lam = (1/cbeta)*mesh_param[level]**(beta)
            G1_l = 0.5 + np.multiply(k1,lam) + np.divide(nsam_local[level],2.0)
            G2_l = k1 + (nsam_local[level]-1)*0.5*variance_local[level] + k0*nsam_local[level]*((mean_local[level]-mu)**2)/(2.0*(k0+nsam_local[level]))
            bayesian_variance.append(np.divide(G2_l,G1_l-0.5))
        self.bayesian_variance = bayesian_variance
        del(bayesian_variance)

    """
    function computing the minimum number of iterations of MLMC algorithm, from [1]:
    iteration < number_iterations_iE : iterations needed to obtain increasingly accurate estimates of the
                                       problem dependent parameters P=[calpha,alpha,cbeta,beta,cgamma,gamma]
    iteration > number_iterations_iE : iterations preventing redundant computations due to fluctuations
                                       in the estimate of P=[calpha,alpha,cbeta,beta,cgamma,gamma]
                                       by solving the problem for a slightly smaller tolerance than the desired one
    input:  self : an instance of the class
    """
    def ComputeNumberIterationsMLMC(self):
        tolF = self.settings["tolerance"].GetDouble()
        tol0 = self.settings["initial_tolerance"].GetDouble()
        r2 = self.settings["r2"].GetDouble()
        r1 = self.settings["r1"].GetDouble()
        self.number_iterations_iE = np.floor((-np.log(tolF)+np.log(r2)+np.log(tol0))/(np.log(r1)))

    """
    function computing the tolerance for the i^th iteration
    input:  self : an instance of the class
    """
    def ComputeTolerancei(self):
        tol0 = self.settings["initial_tolerance"].GetDouble()
        tolF = self.settings["tolerance"].GetDouble()
        if (tol0 <= tolF):
            self.tolerance_i = tolF
        elif(tol0 > tolF):
            r1 = self.settings["r1"].GetDouble()
            current_iter = self.iteration_counter
            tol = tol0 / (r1**(current_iter))
            self.tolerance_i = tol

    """
    function computing the number of levels for i^th iteration of the algorithm
    input:  self :   an instance of the class
    """
    def ComputeLevels(self):
        tol = self.tolerance_i
        Wmin   = 1e10 # high cost to compare with all the level costs (needs to be a high value)
        current_number_levels = self.current_number_levels
        cgamma = self.rates_error["cgamma"]
        gamma  = self.rates_error["gamma"]
        calpha = self.rates_error["calpha"]
        alpha = self.rates_error["alpha"]
        cphi = self.settings["cphi_confidence"].GetDouble()
        mesh_parameters_all_levels = self.mesh_parameters
        Lmax = self.settings["maximum_number_levels"].GetInt()
        Lmin = self.current_number_levels
        if len(self.bayesian_variance) < (Lmax+1):
            self.EstimateBayesianVariance(Lmax)
        bayesian_variance = self.bayesian_variance
        model_cost = np.multiply(cgamma,np.power(mesh_parameters_all_levels,gamma)) # model_cost has length = Lmax + 1
        # observation: not mandatory to increase the number of levels, we may continue using the number of levels of the previous MLMC iteration
        for lev in range(Lmin, Lmax+1):
            # theta_i = 1.0 - (calpha * (mesh_parameters_all_levels[lev])**(-alpha))/tol # use the ComputeTheta function to evaluate theta for level lev
            self.ComputeTheta(lev) # this function bounds theta with a maximum and a minimum value that we set in settings
            theta_i = self.theta_i
            if (theta_i > 0.0) and (theta_i < 1.0):
                coeff2 = np.sum(np.sqrt(np.multiply(model_cost[0:lev+1],bayesian_variance[0:lev+1])))
                coeff2 = coeff2**2.0
                coeff1 = (cphi/(theta_i*tol))**2.0 # formula in case QoI is scalar
            else:
                raise Exception ("The splitting parameter theta_i assumed a value outside the range (0,1): ",self.theta_i)
            Wtot = coeff1 * coeff2
            # print("print level and correspondent cost",lev,Wtot)
            # change number of levels if the cost condition is satisfied
            if Wtot < Wmin:
                Wmin = Wtot
                current_number_levels = lev
        # add maximum one level per time
        if current_number_levels > Lmin:
            current_number_levels = Lmin + 1
        # update length of Bayesian variance with respect to the current levels
        bayesian_variance = bayesian_variance[0:current_number_levels+1]
        # update variables
        self.bayesian_variance = bayesian_variance
        self.current_number_levels = current_number_levels
        del(tol,Wmin,current_number_levels,cgamma,gamma,calpha,alpha,cphi,mesh_parameters_all_levels,Lmax,Lmin)

    """
    function computing the splitting parameter theta \in (0,1)
    input:  self  : an instance of the class
            level : working level
    """
    def ComputeTheta(self,level):
        calpha = self.rates_error["calpha"]
        alpha = self.rates_error["alpha"]
        tol = self.tolerance_i
        mesh_param = self.mesh_parameters[level]
        self.theta_i = 1.0 - (calpha * (mesh_param)**(-alpha))/tol
        if (self.theta_i > self.settings["splitting_parameter_max"].GetDouble()):
            self.theta_i = self.settings["splitting_parameter_max"].GetDouble()
        elif (self.theta_i < self.settings["splitting_parameter_min"].GetDouble()):
            self.theta_i = self.settings["splitting_parameter_min"].GetDouble()
        del(calpha,alpha,tol,mesh_param)

    """
    function computing the updated number of samples
    input:  self : an instance of the class
    """
    def ComputeNumberSamples(self):
        current_number_levels = self.current_number_levels
        bayesian_variance = self.bayesian_variance
        min_samples_add = np.multiply(np.ones(current_number_levels+1),self.settings["minimum_samples_add_level"].GetDouble())
        cgamma = self.rates_error["cgamma"]
        gamma  = self.rates_error["gamma"]
        cphi = self.settings["cphi_confidence"].GetDouble()
        mesh_parameters_current_levels = self.mesh_parameters[0:current_number_levels+1]
        theta = self.theta_i
        tol = self.tolerance_i
        # TODO: now considering self.number_samples, i.e. only the samples used up to this moment for the statistics
        # in future may be that we consider here all the samples, even if not analysed yet
        nsamples = [0 for _ in range (self.settings["maximum_number_levels"].GetInt()+1)]
        for batch in range (len(self.batches_number_samples)):
            for level in range (len(self.batches_number_samples[batch])):
                nsamples[level] = nsamples[level] + self.batches_number_samples[batch][level]
        nsam = copy.copy(nsamples)
        # compute optimal number of samples and store previous number of samples
        coeff1 = (cphi/(theta*tol))**2.0
        model_cost = np.multiply(cgamma,np.power(mesh_parameters_current_levels,gamma))
        coeff2 = np.sqrt(np.divide(bayesian_variance,model_cost))
        coeff3 = np.sum(np.sqrt(np.multiply(model_cost,bayesian_variance)))
        opt_number_samples = np.multiply(coeff1*coeff3,coeff2)
        opt_number_samples = 1.25 * opt_number_samples # safely increasing number of samples since running in parallel
        # print("optimal number of samples computed = ",opt_number_samples)
        if (len(opt_number_samples) != current_number_levels+1):
            raise Exception ("length optimal number of samples and current optimal number of level not coherent")
        for lev in range (current_number_levels+1):
            opt_number_samples[lev] = np.ceil(opt_number_samples[lev])
            opt_number_samples[lev] = opt_number_samples[lev].astype(int)
        if len(nsam) < len(opt_number_samples):
            for _ in range (len(opt_number_samples)-len(nsam)):
                nsam.append(0)
        # copy local string and avoid reference since we modify it
        previous_number_samples = nsam[:]
        # evaluate difference between new and previous number of samples
        diff_nsam = []
        for lev in range(current_number_levels+1):
            diff_nsam.append(opt_number_samples[lev] - nsam[lev])
            # set here that if the optimal number of samples is smaller than
            # the previous number of samples, keep the previous number of samples,
            # i.e. do not decrease the number of samples w.r.t. previous iteration
            if diff_nsam[lev] <= 0.:
                diff_nsam[lev] = 0.
                opt_number_samples[lev] = nsam[lev]
            # set here to have an addition of minimum min_samples_add per level
            if (diff_nsam[lev] > 0.) and (diff_nsam[lev] < min_samples_add[lev]):
                diff_nsam[lev] = min_samples_add[lev]
                opt_number_samples[lev] = nsam[lev] + diff_nsam[lev]
            nsam[lev] = opt_number_samples[lev]
        # convert current number of samples and difference number of samples to integers
        for lev in range (current_number_levels+1):
            diff_nsam[lev] = int(diff_nsam[lev])
            nsam[lev] = int(nsam[lev])
        # update variables and delete local variables
        self.batch_size = diff_nsam
        del(current_number_levels,bayesian_variance,min_samples_add,cgamma,gamma,cphi,mesh_parameters_current_levels,\
        theta,tol,nsam,opt_number_samples,diff_nsam)

    """
    function computing the mlmc estimator for the mean of the Quantity of Interest
    input:  self : an instance of the class
    """
    def ComputeMeanMLMCQoI(self):
        if (self.difference_QoI.h_statistics_1[-1] == []):
            current_number_levels = self.difference_QoI.h_statistics_1.index([]) # index() method finds the given element in a list and returns its position (starts from 1, not from 0)
        else:
            current_number_levels = self.settings["maximum_number_levels"].GetInt()+1
        self.mean_mlmc_QoI = np.sum(self.difference_QoI.h_statistics_1[:current_number_levels])

    """
    function computing the total error:
    total_error = bias contrinution + statistical error contribution
    input:  self : an instance of the class
    """
    def ComputeTotalErrorMLMC(self):
        if (self.difference_QoI.h_statistics_1[-1] == []):
            current_number_levels = self.difference_QoI.h_statistics_1.index([])-1 # index() method finds the given element in a list and returns its position (starts from 1, not from 0)
        else:
            current_number_levels = self.settings["maximum_number_levels"].GetInt()

        self.difference_QoI.bias_error = np.abs(self.difference_QoI.h_statistics_1[current_number_levels])
        variance_from_bayesian = np.zeros(np.size(self.number_samples))
        for lev in range(current_number_levels+1):
            # variance_from_bayesian[lev] = self.bayesian_variance[lev]/self.number_samples[lev]
            variance_from_bayesian[lev] = self.difference_QoI.h_statistics_2[lev]/self.number_samples[lev]
        self.difference_QoI.statistical_error = self.settings["cphi_confidence"].GetDouble() * np.sqrt(np.sum(variance_from_bayesian))
        total_error = self.difference_QoI.bias_error + self.difference_QoI.statistical_error
        self.total_error = np.abs(total_error)


class MultilevelMonteCarloResults(object):
    """The base class for the MultilevelMonteCarloResults-classes"""
    def __init__(self,number_levels):
        """constructor of the MultilevelMonteCarloResults-Object
        Keyword arguments:
        self: an instance of the class
        """
        # Quantity of Interest
        self.QoI = [[] for _ in range (number_levels+1)]
        # time cost
        self.time_ML = [[] for _ in range (number_levels+1)]
        # level of QoI and time_ML
        self.finer_level = number_levels
