from __future__ import absolute_import, division # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import numpy as np
import KratosMultiphysics
import time
import copy

# Import refinement library
import adaptive_refinement_utilities as refinement

# Import cpickle to pickle the serializer
try:
    import cpickle as pickle  # Use cPickle on Python 2.7
except ImportError:
    import pickle

# Import exaqute
from exaqute.ExaquteTaskPyCOMPSs import *   # to exequte with pycompss
# from exaqute.ExaquteTaskHyperLoom import *  # to exequte with the IT4 scheduler
# from exaqute.ExaquteTaskLocal import *      # to execute with python3
'''
get_value_from_remote is the equivalent of compss_wait_on: a synchronization point
in future, when everything is integrated with the it4i team, importing exaqute.ExaquteTaskHyperLoom you can launch your code with their scheduler instead of BSC
'''


'''
This utility contains all the compss tasks to perform the Continuation Multilevel Monte Carlo (CMLMC) algorithm
'''


'''
auxiliary function of UpdateOnepassMeanVariance of the StatisticalVariable class
this function is needed since in compss we do operations among future objects,
and we need to handle the singular future values
'''
@ExaquteTask(returns=4)
def UpdateOnepassMeanVarianceAux(sample, old_mean, old_M2, nsamples):
    nsamples = nsamples + 1
    if nsamples == 1:
        new_mean = sample
        new_M2 = np.zeros(np.size(sample))
        new_M2 = np.asscalar(new_M2) # do so to have a list of scalars, and not a list of arrays of one element
        new_sample_variance = np.zeros(np.size(sample))
        new_sample_variance = np.asscalar(new_sample_variance) # do so to have a list of scalars, and not a list of arrays of one element
    else:
        delta = np.subtract(sample, old_mean)
        new_mean = old_mean + np.divide(delta,nsamples)
        new_M2 = old_M2 + delta*np.subtract(sample,new_mean)
        new_sample_variance = np.divide(new_M2,np.subtract(nsamples,1))
    return new_mean, new_M2, new_sample_variance, nsamples


'''
auxiliary function of AddResults of the MultilevelMonteCarlo class
this function is needed since in compss we do operations among future objects,
and we need to handle the singular future values
'''
@ExaquteTask(returns=2)
def AddResultsAux(simulation_results,QoI_values_level,timeML_values_level):
    difference_QoI_value = simulation_results["QoI_finer_level"] - simulation_results["QoI_coarser_level"]
    return difference_QoI_value, simulation_results["total_MLMC_time"]


'''
auxiliary function finalizing the screening phase and the MLMC phase in the MultilevelMonteCarlo class
'''
def FinalizePhaseTask(ConstructorCallback,difference_QoI_mean,difference_QoI_sample_variance,time_ML_mean,aux_settings,aux_mesh_parameters,\
aux_current_number_levels,aux_current_iteration,aux_number_samples):
    '''create an auxiliary object equivalent to the current one of the problem'''
    auxiliary_settings = KratosMultiphysics.Parameters("""{ }""")
    auxiliary_MLMC_object = ConstructorCallback(auxiliary_settings)
    auxiliary_MLMC_object.settings = aux_settings
    auxiliary_MLMC_object.difference_QoI.mean = difference_QoI_mean
    auxiliary_MLMC_object.difference_QoI.sample_variance = difference_QoI_sample_variance
    auxiliary_MLMC_object.time_ML.mean = time_ML_mean
    auxiliary_MLMC_object.mesh_parameters = aux_mesh_parameters
    auxiliary_MLMC_object.current_number_levels = aux_current_number_levels
    auxiliary_MLMC_object.current_iteration = aux_current_iteration
    auxiliary_MLMC_object.number_samples = aux_number_samples
    '''compute the functions needed to finalize the screening phase or the MLMC phase'''
    if (auxiliary_MLMC_object.current_iteration == 0):
        '''compute parameters by least square fit'''
        auxiliary_MLMC_object.ComputeRatesLS()
        '''compute Bayesian variance V^c[Y_l]'''
        auxiliary_MLMC_object.EstimateBayesianVariance(auxiliary_MLMC_object.current_number_levels)
    else:
        '''compute estimatior MLMC mean QoI'''
        auxiliary_MLMC_object.ComputeMeanMLMCQoI()
        '''compute parameters by least square fit'''
        auxiliary_MLMC_object.ComputeRatesLS()
        '''compute Bayesian variance V^c[Y_l]'''
        auxiliary_MLMC_object.EstimateBayesianVariance(auxiliary_MLMC_object.current_number_levels)
        '''compute total error of the MLMC simulation'''
        auxiliary_MLMC_object.ComputeTotalErrorMLMC()
    return auxiliary_MLMC_object.rates_error["calpha"],auxiliary_MLMC_object.rates_error["alpha"],\
    auxiliary_MLMC_object.rates_error["cbeta"],auxiliary_MLMC_object.rates_error["beta"],\
    auxiliary_MLMC_object.rates_error["cgamma"],auxiliary_MLMC_object.rates_error["gamma"],\
    auxiliary_MLMC_object.BayesianVariance,auxiliary_MLMC_object.mean_mlmc_QoI,\
    auxiliary_MLMC_object.difference_QoI.bias_error,auxiliary_MLMC_object.difference_QoI.statistical_error,\
    auxiliary_MLMC_object.TErr,auxiliary_MLMC_object.number_samples,auxiliary_MLMC_object.BayesianVariance