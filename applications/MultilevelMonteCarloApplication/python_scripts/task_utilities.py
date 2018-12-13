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
function evaluationg the QoI and the cost of simulation, computing the mesh of level finest_level
refining recursively from the coarsest mesh
input:
        finest_level              : current Multilevel MOnte Carlo level we are solving
        pickled_coarse_model      : pickled model
        pickled_coarse_parameters : pickled parameters
        size_meshes               : mesh sizes for all levels
output:
        results_simulation : QoI_finer_level   : QoI of level fisest_level
                             QoI_coarser_level : QoI of level finest_level - 1
                             finer_level       : finest level
                             coarser_level     : finest_level - 1
                             total_MLMC_time   : execution time
'''
@ExaquteTask(returns=1)
def ExecuteMultilevelMonteCarloAnalisys(MultilevelMonteCarloAnalysisCallback,GenerateSampleCallback,EvaluateQuantityOfInterestCallback,finest_level,pickled_coarse_model,pickled_coarse_parameters,size_meshes):
    '''overwrite the old model serializer with the unpickled one'''
    model_serializer = pickle.loads(pickled_coarse_model)
    current_model = KratosMultiphysics.Model()
    model_serializer.Load("ModelSerialization",current_model)
    del(model_serializer)
    '''overwrite the old parameters serializer with the unpickled one'''
    serialized_parameters = pickle.loads(pickled_coarse_parameters)
    current_parameters = KratosMultiphysics.Parameters()
    serialized_parameters.Load("ParametersSerialization",current_parameters)
    del(serialized_parameters)
    '''generate the sample'''
    sample = GenerateSampleCallback()
    QoI = []
    start_MLMC_time = time.time()
    if(finest_level == 0):
        QoI.append(0.0)
        simulation = MultilevelMonteCarloAnalysisCallback(current_model,current_parameters,sample)
        simulation.Run()
        QoI.append(EvaluateQuantityOfInterestCallback(simulation))
        del(simulation)
    else:
        for lev in range(finest_level+1):
            simulation = MultilevelMonteCarloAnalysisCallback(current_model,current_parameters,sample)
            simulation.Run()
            QoI.append(EvaluateQuantityOfInterestCallback(simulation))
            '''refine if level < finest level exploiting the solution just computed'''
            if (lev < finest_level):
                '''refine the model Kratos object'''
                model_refined = refinement.compute_refinement_hessian_metric(simulation,size_meshes[lev+1],size_meshes[lev])
                '''initialize the model Kratos object'''
                simulation = MultilevelMonteCarloAnalysisCallback(model_refined,current_parameters,sample)
                simulation.Initialize()
                '''update model Kratos object'''
                current_model = simulation.model
            del(simulation)
    end_MLMC_time = time.time()
    '''prepare results of the simulation'''
    results_simulation = {
        "QoI_finer_level":QoI[-1],\
        "QoI_coarser_level":QoI[-2],\
        "finer_level":finest_level,\
        "coarser_level":np.maximum(0,finest_level-1),\
        "total_MLMC_time":end_MLMC_time - start_MLMC_time}
    return results_simulation


'''
function serializing and pickling the model and the parameters of the problem
the idea is the following:
i)   from Model/Parameters Kratos object to StreamSerializer Kratos object
ii)  from StreamSerializer Kratos object to pickle string
iii) from pickle string to StreamSerializer Kratos object
iv)  from StreamSerializer Kratos object to Model/Parameters Kratos object
input:
        parameter_file_name   : path of the Project Parameters file
output:
        serialized_model      : model serialized
        serialized_parameters : project parameters serialized
'''
@ExaquteTask(parameter_file_name=FILE_IN,returns=2)
def SerializeModelParameters(MultilevelMonteCarloAnalysisCallback,parameter_file_name):
    with open(parameter_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    local_parameters = parameters
    model = KratosMultiphysics.Model()
    # local_parameters["solver_settings"]["model_import_settings"]["input_filename"].SetString(model_part_file_name[:-5])
    fake_sample = 1.0
    '''initialize'''
    simulation = MultilevelMonteCarloAnalysisCallback(model,local_parameters,fake_sample)
    simulation.Initialize()
    '''save and pickle model and parameters'''
    serialized_model = KratosMultiphysics.StreamSerializer()
    serialized_model.Save("ModelSerialization",simulation.model)
    serialized_parameters = KratosMultiphysics.StreamSerializer()
    serialized_parameters.Save("ParametersSerialization",simulation.project_parameters)
    # pickle dataserialized_data
    pickled_model = pickle.dumps(serialized_model, 2) # second argument is the protocol and is NECESSARY (according to pybind11 docs)
    pickled_parameters = pickle.dumps(serialized_parameters, 2) # second argument is the protocol and is NECESSARY (according to pybind11 docs)
    print("\n############## Serialization completed ##############\n")
    return pickled_model, pickled_parameters


'''
function executing the problem
input:
        model       : serialization of the model
        parameters  : serialization of the Project Parameters
output:
        QoI                   : Quantity of Interest
        serialized_model      : model serialized
        serialized_parameters : parameters serialized
'''
@ExaquteTask(returns=2)
def ExecuteRefinement(MultilevelMonteCarloAnalysisCallback,EvaluateQuantityOfInterestCallback,pickled_model_coarse, pickled_parameters, min_size, max_size):
    fake_sample = 1.0
    '''overwrite the old model serializer with the unpickled one'''
    model_serializer_coarse = pickle.loads(pickled_model_coarse)
    model_coarse = KratosMultiphysics.Model()
    model_serializer_coarse.Load("ModelSerialization",model_coarse)
    del(model_serializer_coarse)
    '''overwrite the old parameters serializer with the unpickled one'''
    serialized_parameters = pickle.loads(pickled_parameters)
    parameters_refinement = KratosMultiphysics.Parameters()
    serialized_parameters.Load("ParametersSerialization",parameters_refinement)
    del(serialized_parameters)
    simulation_coarse = MultilevelMonteCarloAnalysisCallback(model_coarse,parameters_refinement,fake_sample)
    simulation_coarse.Run()
    QoI =  EvaluateQuantityOfInterestCallback(simulation_coarse)
    '''refine'''
    model_refined = refinement.compute_refinement_hessian_metric(simulation_coarse,min_size,max_size)
    '''initialize'''
    simulation = MultilevelMonteCarloAnalysisCallback(model_refined,parameters_refinement,fake_sample)
    simulation.Initialize()
    '''serialize model and pickle it'''
    serialized_model = KratosMultiphysics.StreamSerializer()
    serialized_model.Save("ModelSerialization",simulation.model)
    pickled_model_refined = pickle.dumps(serialized_model, 2)
    return QoI,pickled_model_refined


'''
function computing the relative error between the Multilevel Monte Carlo expected value and the exact expected value
input :
        AveragedMeanQoI       : Multilevel Monte Carlo expected value
        ExactExpectedValueQoI : exact expected value
output :
        relative_error        : relative error
'''
@ExaquteTask(returns=1)
def CompareMean(AveragedMeanQoI,ExactExpectedValueQoI):
    relative_error = abs((AveragedMeanQoI - ExactExpectedValueQoI)/ExactExpectedValueQoI)
    return relative_error



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