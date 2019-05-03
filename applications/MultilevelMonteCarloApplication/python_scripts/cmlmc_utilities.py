from __future__ import absolute_import, division # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import packages
import numpy as np
from scipy.stats import norm
import copy
import time

# Importing the analysis stage classes of the different problems
from simulation_definition import SimulationScenario

# Import the StatisticalVariable class
from KratosMultiphysics.MultilevelMonteCarloApplication.statistical_variable_utilities import StatisticalVariable

# Import refinement library
from KratosMultiphysics.MultilevelMonteCarloApplication.adaptive_refinement_utilities import AdaptiveRefinement

# Import random variable generator
import KratosMultiphysics.MultilevelMonteCarloApplication.generator_utilities as generator

# Import PyCOMPSs
from exaqute.ExaquteTaskPyCOMPSs import *   # to execute with runcompss
# from exaqute.ExaquteTaskHyperLoom import *  # to execute with the IT4 scheduler
# from exaqute.ExaquteTaskLocal import *      # to execute with python3
"""
get_value_from_remote is the equivalent of compss_wait_on: a synchronization point
in future, when everything is integrated with the it4i team, importing exaqute.ExaquteTaskHyperLoom you can launch your code with their scheduler instead of BSC
"""

# Import cpickle to pickle the serializer
try:
    import cpickle as pickle  # Use cPickle on Python 2.7
except ImportError:
    import pickle

"""
This utility contains the functions to perform the Continuation Multilevel Monte Carlo (CMLMC) algorithm
References:
M. Pisaroni, F. Nobile, P. Leyland; A Continuation Multi Level Monte Carlo (C-MLMC) method for uncertainty quantification in compressible inviscid aerodynamics; Computer Methods in Applied Mechanics and Engineering, vol 326, pp 20-50, 2017. DOI : 10.1016/j.cma.2017.07.030.
M. Pisaroni, S. Krumscheid, F. Nobile; Quantifying uncertain system outputs via the multilevel Monte Carlo method - Part I: Central moment estimation ;  available as MATHICSE technical report no. 23.2017
"""


# TODO: insert distinction between scalar and field Quantity of Interests
# TODO: create function updating lists
# TODO: run cmlmc without refinement, i.e. adaptive_refinement setting


"""
auxiliary function of AddResults of the MultilevelMonteCarlo class
input:  simulation_results: an instance of the multilevel monte carlo result class
        level:              working level
output: difference_QoI_value:              difference QoIvalue to be added
        simulation_results.time_ML[level]: time value to be added
"""
@ExaquteTask(returns=2)
def AddResultsAux_Task(simulation_results,level):
    if (level == 0):
        # each value is inside the relative level list, and only one value per level is computed
        # i.e. results = [[value_level_0],[value_level_1],...]
        difference_QoI_value = simulation_results.QoI[level][0]
    else:
        difference_QoI_value = simulation_results.QoI[level][0] - simulation_results.QoI[level-1][0]
    return difference_QoI_value,simulation_results.time_ML[level]


"""
auxiliary function finalizing the screening phase and the MLMC phase of the MultilevelMonteCarlo class
input:  ConstructorCallback:       an instance of the MultilevelMonteCarlo class
        aux_settings_serialized:   StreamSerializer Kratos class containing the settings of the MultilevelMonteCarlo class
        aux_mesh_parameters:       mesh parameters of the MltilevelMOnteCarlo class
        aux_current_number_levels: current number of levels of the MultilevelMonteCarlo class
        aux_current_iteration:     current iteration of the MultilevelMonteCarlo class
        aux_number_samples:        number of samples of the Multilevel Monte Carlo class
        *args:                     list containing difference_QoI.raw_moment_1, difference_QoI.unbiased_central_moment_2 and time_ML.raw_moment_1
output: auxiliary_MLMC_object.rates_error:       rates error of the MultilevelMonteCarlo class
        auxiliary_MLMC_object.bayesian_variance: bayesian variance of the MultilevelMonteCarlo class
        auxiliary_MLMC_object.mean_mlmc_QoI:     multilevel monte carlo mean estimator of the MultilevelMonteCarlo class
        auxiliary_MLMC_object.total_error:       total error of the MultilevelMonteCarlo class
        auxiliary_MLMC_object.number_samples:    number of samples of the MultilevelMonteCarlo class
"""
@ExaquteTask(returns=5)
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
    auxiliary_MLMC_object.difference_QoI.raw_moment_1 = difference_QoI_mean
    auxiliary_MLMC_object.difference_QoI.unbiased_central_moment_2 = difference_QoI_sample_variance
    auxiliary_MLMC_object.time_ML.raw_moment_1 = time_ML_mean
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
function evaluating the QoI and the cost of the simulation, computing the mesh of current level (current_MLMC_level)
refining recursively from the coarsest mesh
input:
        current_MLMC_level:                current Multilevel Monte Carlo level we are solving
        pickled_coarse_model:              pickled model
        pickled_coarse_project_parameters: pickled parameters
        mesh_sizes:                        mesh sizes for all levels
output:
        mlmc_results: QoI:                list of QoI for all levels computed in the current simulation
                      finer_level:        finest level
                      time_ML:            list of MLMC time for all levels computed in the current simulation
        pickled_finer_model:              pickled finer model
        pickled_finer_project_parameters: pickled finer parameters
"""
@ExaquteTask(returns=3)
def ExecuteInstanceAux_Task(current_MLMC_level,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,mesh_sizes,recursive_maximal_size,sample,current_level,current_analysis_stage,mlmc_results):
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
    # start time
    start_MLMC_time = time.time()
    # refine if current current_level > 0, adaptive refinement based on the solution of previous level
    if (current_level > 0):
        # unpickle metric and remesh refinement parameters and build Kratos Parameters objects
        serialized_custom_metric_refinement_parameters = pickle.loads(pickled_custom_metric_refinement_parameters)
        serialized_custom_remesh_refinement_parameters = pickle.loads(pickled_custom_remesh_refinement_parameters)
        current_custom_metric_refinement_parameters = KratosMultiphysics.Parameters()
        current_custom_remesh_refinement_parameters = KratosMultiphysics.Parameters()
        serialized_custom_metric_refinement_parameters.Load("MetricRefinementParametersSerialization",current_custom_metric_refinement_parameters)
        serialized_custom_remesh_refinement_parameters.Load("RemeshRefinementParametersSerialization",current_custom_remesh_refinement_parameters)
        del(serialized_custom_metric_refinement_parameters,serialized_custom_remesh_refinement_parameters)
        # set minimal and maximal mesh sizes
        minimal_mesh_size_level = mesh_sizes[current_level]
        if (recursive_maximal_size == "True"):
            maximal_mesh_size_level = mesh_sizes[current_level-1]
        else:
            maximal_mesh_size_level = 10.0
        # refine the model Kratos object
        adaptive_refinement_manager = AdaptiveRefinement(current_model,current_project_parameters,current_custom_metric_refinement_parameters,current_custom_remesh_refinement_parameters,minimal_mesh_size_level,maximal_mesh_size_level)
        refined_model,refined_project_parameters = adaptive_refinement_manager.ComputeAdaptiveRefinement()
        # initialize the model Kratos object
        simulation = current_analysis_stage(refined_model,refined_project_parameters,sample)
        simulation.Initialize()
        # update model Kratos object
        current_model = simulation.model
        current_project_parameters = simulation.project_parameters
        del(simulation)
    simulation = current_analysis_stage(current_model,current_project_parameters,sample)
    simulation.Run()
    QoI = simulation.EvaluateQuantityOfInterest()
    # save model and parameters as StreamSerializer Kratos objects
    serialized_finer_model = KratosMultiphysics.StreamSerializer()
    serialized_finer_model.Save("ModelSerialization",simulation.model)
    serialized_finer_project_parameters = KratosMultiphysics.StreamSerializer()
    serialized_finer_project_parameters.Save("ParametersSerialization",simulation.project_parameters)
    # pickle model and parameters
    pickled_finer_model = pickle.dumps(serialized_finer_model, 2) # second argument is the protocol and is NECESSARY (according to pybind11 docs)
    pickled_finer_project_parameters = pickle.dumps(serialized_finer_project_parameters, 2) # second argument is the protocol and is NECESSARY (according to pybind11 docs)
    del(simulation)
    end_MLMC_time = time.time()
    # register results of the current level in the MultilevelMonteCarloResults class
    mlmc_results.time_ML[current_level].append(end_MLMC_time-start_MLMC_time) # saving each result in the corresponding list in order to ensure the correctness of the results order and the levels
    mlmc_results.QoI[current_level].append(QoI) # saving each result in the corresponding list in order to ensure the correctness of the results order and the levels
    return mlmc_results,pickled_finer_model,pickled_finer_project_parameters


class MultilevelMonteCarlo(object):
    """The base class for the MultilevelMonteCarlo-classes"""
    def __init__(self,custom_parameters_path,project_parameters_path,refinement_parameters_path,custom_analysis):
        """constructor of the MultilevelMonteCarlo-Object
        Keyword arguments:
        self:                              an instance of the class
        custom_parameters_path:            path of settings of the Multilevel Monte Carlo simulation
        project_parameters_path:           path of the project parameters file
        refinement_parameters_path:        path of settings of the metric for the refinement and settings of the remeshing
        custom_analysis:                   analysis stage of the problem
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
        # adaptive_refinement: if True refine from mesh level 0 each time a different random variable is used
        # store_lower_levels_samples: if True, store also lower level values of QoI to compute statistics, if avaiable from adaptive refinement --> produces a BIAS in the results
        # initial_tolerance: tolerance iter 0
        # tolerance: tolerance final
        # confidence: confidence on tolerance
        # cphi_confidence: CDF**-1 (confidence) where CDF**-1 is the inverse of the CDF of the standard normal distribution, the default value is computed for confidence = 0.9
        # number_samples_screening: number of samples for screening phase
        # number_levels_screening: number of levels for screening phase
        # maximum_number_levels: maximum number of levels
        # mesh_refinement_coefficient: coefficient of mesh refinement
        # initial_mesh_size: size of coarsest/initial mesh
        # recursive_maximal_size : set as maximal mesh size the minimal mesh size of the previous coarser level, otherwise maximal mesh size default equal to 10.0 (directly set into ExecuteInstanceAux_Task)
        # minimum_samples_add: minimum number of samples to add if at least one is going to be added
        # k0: certainty Parameter 0 rates (confidence in the variance models)
        # k1: certainty Parameter 1 rates (confidence in the weak convergence model)
        # r1: cost increase first iterations CMLMC
        # r2: cost increase final iterations CMLMC
        # splitting_parameter_max : max value parameter splitting discretization and statistical errors
        # splitting_parameter_min : min value parameter splitting discretization and statistical errors
        default_settings = KratosMultiphysics.Parameters("""
        {
            "run_multilevel_monte_carlo" : "True",
            "adaptive_refinement"        : "True",
            "store_lower_levels_samples" : "False",
            "k0" : 0.1,
            "k1" : 0.1,
            "r1" : 1.25,
            "r2" : 1.15,
            "initial_tolerance" : 0.25,
            "tolerance" : 0.1,
            "confidence" : 0.9,
            "cphi_confidence" : 1.28155156554,
            "number_samples_screening" : 25,
            "levels_screening" : 2,
            "maximum_number_levels" : 4,
            "mesh_refinement_coefficient" : 2,
            "recursive_maximal_size" : "False",
            "initial_mesh_size" : 0.5,
            "minimum_samples_add_level" : 6.0,
            "splitting_parameter_max" : 0.9,
            "splitting_parameter_min" : 0.1
        }
        """)
        if (type(custom_parameters_path) == KratosMultiphysics.Parameters):
            self.settings = custom_parameters_path # needed in FinalizePhaseAux_Task to build an auxiliary MultilevelMonteCarlo class, then in FinalizePhaseAux_Task we set the true problem parameters
        else: # standard contructor
            self.custom_parameters_path = custom_parameters_path
            self.SetXMCParameters()
            # warning if initial_mesh_size parameter not set by the user
            if not (self.settings.Has("initial_mesh_size")):
                print("\n ######## WARNING: initial_mesh_size parameter not set --> using defalut value 0.5 as initial minimum size ########\n")
            # compute cphi = CDF**-1 (confidence)
            self.settings.AddEmptyValue("cphi_confidence")
            if (self.settings["confidence"].GetDouble()==1.0):
                self.settings["confidence"].SetDouble(0.999) # reduce confidence to not get +inf for cphi_confidence
            self.settings["cphi_confidence"].SetDouble(norm.ppf(self.settings["confidence"].GetDouble()))
        # validate and assign default parameters
        self.settings.ValidateAndAssignDefaults(default_settings)
        # warning if store_lower_levels_samples is True
        if (self.settings["store_lower_levels_samples"].GetString() == "True"):
            print("\n ######## WARNING: store_lower_levels_samples set to True --> introducing a BIAS in the problem ########\n")
        # current_number_levels: number of levels of current iteration
        self.current_number_levels = self.settings["levels_screening"].GetInt()
        # previous_number_levels: number of levels of previous iteration
        self.previous_number_levels = None
        # current_level: current level of work
        self.current_level = 0
        # number_samples: total number of samples of current iteration
        self.number_samples = [self.settings["number_samples_screening"].GetInt() for _ in range (self.settings["levels_screening"].GetInt()+1)]
        # difference_number_samples: difference between number of samples of current and previous iterations
        self.difference_number_samples = [self.settings["number_samples_screening"].GetInt() for _ in range (self.settings["levels_screening"].GetInt()+1)]
        # previous_number_samples: total number of samples of previous iteration
        self.previous_number_samples = None
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
        self.current_iteration = 0
        # tolerance_i: tolerance of i^th-iteration of MLMC algorithm
        self.tolerance_i = None
        # theta_i: splitting parameter \in (0,1), this affects bias and statistical error in the computation of the total error
        self.theta_i = None
        # mean_mlmc_QoI: MLMC estimator for the mean value of the Quantity of Interest
        self.mean_mlmc_QoI = None
        # total_error: total error of MLMC algorithm, the sum of bias and statistical error is an overestmation of the real total error
        #               total_error := \abs(E^MLMC[QoI] - E[QoI])
        self.total_error = None
        # compute mesh parameter for each mesh
        self.ComputeMeshParameters()

        # difference_QoI: Quantity of Interest of the considered problem organized per levels
        #                 difference_QoI.values := Y_l = QoI_M_l - Q_M_l-1
        self.difference_QoI = StatisticalVariable(self.current_number_levels)
        self.difference_QoI.values = [[] for _ in range (self.settings["levels_screening"].GetInt()+1)] # list containing Y_{l}^{i} = Q_{m_l} - Q_{m_{l-1}}
        self.difference_QoI.type = "scalar"
        # time_ML: time to perform a single MLMC simulation (i.e. one value of difference_QoI.values) organized per levels
        self.time_ML = StatisticalVariable(self.current_number_levels)
        self.time_ML.values = [[] for _ in range (self.settings["levels_screening"].GetInt()+1)] # list containing the time to compute the level=l simulations

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
        # pickled_project_parameters: serialization of project parameters Kratos object of the problem
        self.pickled_project_parameters = None
        # pickled_custom_settings_metric_refinement: serialization of the metric refinement custom settings
        self.pickled_custom_metric_refinement_parameters = None
        # pickled_custom_settings_remesh_refinement: serialization of the remesh custom settings
        self.pickled_custom_remesh_refinement_parameters = None
        # construct the pickled model and pickled project parameters of the problem
        self.is_project_parameters_pickled = False
        self.is_model_pickled = False
        self.SerializeModelParameters()
        # construct the pickled custom settings metric refinement and the picled custom settings remesh of the problem
        self.is_custom_settings_metric_refinement_pickled = False
        self.is_custom_settings_remesh_refinement_pickled = False
        self.SerializeRefinementParameters()
        # construct the pickled custom settings metric refinement and the picled custom settings remesh of the problem
        self.is_custom_settings_metric_refinement_pickled = False
        self.is_custom_settings_remesh_refinement_pickled = False
        self.SerializeRefinementParameters()

    """
    function executing the CMLMC algorithm
    input:  self: an instance of the class
    """
    def Run(self):
        if (self.settings["run_multilevel_monte_carlo"].GetString() == "True"):
            # start screening phase
            self.LaunchEpoch()
            # finalize screening phase
            self.FinalizeScreeningPhase()
            self.ScreeningInfoScreeningPhase()
            # start MLMC phase
            while self.convergence is not True:
                # initialize MLMC phase
                self.InitializeMLMCPhase()
                self.ScreeningInfoInitializeMLMCPhase()
                # MLMC execution phase
                self.LaunchEpoch()
                # finalize MLMC phase
                self.FinalizeMLMCPhase()
                self.ScreeningInfoFinalizeMLMCPhase()
        else:
            print("\n","#"*50,"Not running Multilevel Monte Carlo algorithm","#"*50)
            pass

    """
    function running one Monte Carlo epoch
    input:  self: an instance of the class
    """
    def LaunchEpoch(self):
        for level in range(self.current_number_levels+1):
            for _ in range (self.difference_number_samples[level]):
                self.AddResults(self.ExecuteInstance(level))

    """
    function executing an instance of the CMLMC algorithm, i.e. a single MC simulation and eventually the refinement (that occurs before the simulation run)
    input:  self: an instance of the class
    output: mlmc_results:       instance of the MultilevelMonteCarloResult class
            current_MLMC_level: level of the current MLMC simulation
    """
    def ExecuteInstance(self,current_MLMC_level):
        # local variables
        pickled_coarse_model = self.pickled_model
        pickled_coarse_project_parameters = self.pickled_project_parameters
        mesh_sizes = self.mesh_sizes
        recursive_maximal_size = self.settings["recursive_maximal_size"].GetString() # check if setting or not maximal mesh size
        pickled_custom_metric_refinement_parameters = self.pickled_custom_metric_refinement_parameters
        pickled_custom_remesh_refinement_parameters = self.pickled_custom_remesh_refinement_parameters
        current_analysis = self.analysis
        # generate the sample
        sample = generator.GenerateSample()
        # initialize the MultilevelMonteCarloResults class
        mlmc_results = MultilevelMonteCarloResults(current_MLMC_level)
        if (current_MLMC_level == 0):
            current_level = 0
            mlmc_results,pickled_current_model,pickled_current_project_parameters = \
                ExecuteInstanceAux_Task(current_MLMC_level,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,mesh_sizes,recursive_maximal_size,sample,current_level,current_analysis,mlmc_results)
        else:
            for current_level in range(current_MLMC_level+1):
                mlmc_results,pickled_current_model,pickled_current_project_parameters = \
                    ExecuteInstanceAux_Task(current_MLMC_level,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,mesh_sizes,recursive_maximal_size,sample,current_level,current_analysis,mlmc_results)
                del(pickled_coarse_model,pickled_coarse_project_parameters)
                pickled_coarse_model = pickled_current_model
                pickled_coarse_project_parameters = pickled_current_project_parameters
                del(pickled_current_model,pickled_current_project_parameters)
        return mlmc_results,current_MLMC_level

    """
    function serializing and pickling the model and the project parameters of the problem
    the serialization-pickling process is the following:
    i)   from Model/Parameters Kratos object to StreamSerializer Kratos object
    ii)  from StreamSerializer Kratos object to pickle string
    iii) from pickle string to StreamSerializer Kratos object
    iv)  from StreamSerializer Kratos object to Model/Parameters Kratos object
    requires: self.project_parameters_path: path of the Project Parameters file
    builds: self.pickled_model:              pickled model
            self.pickled_project_parameters: pickled project parameters
    input:  self: an instance of the class
    """
    def SerializeModelParameters(self):
        if (self.to_pickle_model_parameters):
            with open(self.project_parameters_path,'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())
            model = KratosMultiphysics.Model()
            fake_sample = generator.GenerateSample() # only used to serialize
            simulation = self.analysis(model,parameters,fake_sample)
            simulation.Initialize()
            serialized_model = KratosMultiphysics.StreamSerializer()
            serialized_model.Save("ModelSerialization",simulation.model)
            serialized_project_parameters = KratosMultiphysics.StreamSerializer()
            serialized_project_parameters.Save("ParametersSerialization",simulation.project_parameters)
            # pickle dataserialized_data
            pickled_model = pickle.dumps(serialized_model, 2) # second argument is the protocol and is NECESSARY (according to pybind11 docs)
            pickled_project_parameters = pickle.dumps(serialized_project_parameters, 2)
            self.pickled_model = pickled_model
            self.pickled_project_parameters = pickled_project_parameters
            self.is_project_parameters_pickled = True
            self.is_model_pickled = True
            print("\n","#"*50," SERIALIZATION MODEL AND PROJECT PARAMETERS COMPLETED ","#"*50,"\n")

    """
    function reading the XMC parameters passed from json file
    input:  self: an instance of the class
    """
    def SetXMCParameters(self):
        with open(self.custom_parameters_path,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        self.settings = parameters["multilevel_monte_carlo"]

    """
    function reading the refinement parameters passed from json file
    input:  self: an instance of the class
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
    input:  self: an instance of the class
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
    input:  self:                       an instance of the class
            application_analysis_stage: analysis stage of the problem
    """
    def SetAnalysis(self,application_analysis_stage):
        self.analysis = application_analysis_stage

    """
    function returning the Kratos analysis stage of the problem
    input:  self: an instance of the class
    output: self.analysis: analysis stage of the problem
    """
    def GetAnalysis(self):
        if (self.analysis is not None):
            return self.analysis
        else:
            print("Please provide a Kratos specific application analysis stage for the current problem.")

    """
    function finalizing the screening phase of MLMC algorithm
    usage: designed to be called ONCE, AFTER the screening phase
    input:  self: an instance of the class
    """
    def FinalizeScreeningPhase(self):
        # prepare lists
        self.difference_QoI.raw_moment_1 = [[] for _ in range (self.settings["levels_screening"].GetInt()+1)]
        self.difference_QoI.unbiased_central_moment_2 = [[] for _ in range (self.settings["levels_screening"].GetInt()+1)]
        self.difference_QoI.central_moment_1 = [[] for _ in range (self.settings["levels_screening"].GetInt()+1)]
        self.difference_QoI.central_moment_2 = [[] for _ in range (self.settings["levels_screening"].GetInt()+1)]
        self.difference_QoI.central_moment_3 = [[] for _ in range (self.settings["levels_screening"].GetInt()+1)]
        self.difference_QoI.central_moment_4 = [[] for _ in range (self.settings["levels_screening"].GetInt()+1)]
        self.time_ML.raw_moment_1 = [[] for _ in range (self.settings["levels_screening"].GetInt()+1)]
        self.time_ML.unbiased_central_moment_2 = [[] for _ in range (self.settings["levels_screening"].GetInt()+1)]
        self.time_ML.central_moment_1 = [[] for _ in range (self.settings["levels_screening"].GetInt()+1)]
        self.time_ML.central_moment_2 = [[] for _ in range (self.settings["levels_screening"].GetInt()+1)]
        self.time_ML.central_moment_3 = [[] for _ in range (self.settings["levels_screening"].GetInt()+1)]
        self.time_ML.central_moment_4 = [[] for _ in range (self.settings["levels_screening"].GetInt()+1)]
        # compute mean, sample variance and second moment for difference QoI and time ML
        for level in range (self.current_number_levels+1):
            for i_sample in range(self.number_samples[level]):
                self.difference_QoI.UpdateOnePassCentralMoments(level,i_sample)
                self.time_ML.UpdateOnePassCentralMoments(level,i_sample)
        # compute i_E, number of iterations of Multilevel Monte Carlo algorithm
        self.ComputeNumberIterationsMLMC()
        """
        the functions executed in FinalizePhaseAux_Task are the followings:
        self.ComputeRatesLS()
        self.EstimateBayesianVariance(self.current_number_levels)
        """
        # store lists in synchro_element to execute FinalizePhaseAux_Task
        synchro_elements = [x for x in self.difference_QoI.raw_moment_1]
        synchro_elements.extend(["%%%"])
        synchro_elements.extend(self.difference_QoI.unbiased_central_moment_2)
        synchro_elements.extend(["%%%"])
        synchro_elements.extend(self.time_ML.raw_moment_1)
        # create a StreamSerializer Kratos object containing the problem settings
        serial_settings = KratosMultiphysics.StreamSerializer()
        serial_settings.Save("ParametersSerialization", self.settings)
        # compute remaining MLMC finalize process operations
        # observation: we are passing self.settings and we will exploit it to construct the class
        self.rates_error,self.bayesian_variance,self.mean_mlmc_QoI,\
        self.total_error,self.number_samples\
        = FinalizePhaseAux_Task(MultilevelMonteCarlo,
        serial_settings,self.mesh_parameters,self.current_number_levels,\
        self.current_iteration,self.number_samples,*synchro_elements)
        # synchronization point needed to compute the other functions of MLMC algorithm
        # put as in the end as possible the synchronization point
        self.difference_QoI.raw_moment_1 = get_value_from_remote(self.difference_QoI.raw_moment_1)
        self.difference_QoI.unbiased_central_moment_2 = get_value_from_remote(self.difference_QoI.unbiased_central_moment_2)
        self.time_ML.raw_moment_1 = get_value_from_remote(self.time_ML.raw_moment_1)
        self.total_error = get_value_from_remote(self.total_error)
        self.rates_error = get_value_from_remote(self.rates_error)
        self.bayesian_variance = get_value_from_remote(self.bayesian_variance)
        self.mean_mlmc_QoI = get_value_from_remote(self.mean_mlmc_QoI)
        self.number_samples = get_value_from_remote(self.number_samples)
        # start first iteration, we enter in the MLMC algorithm
        self.current_iteration = 1

    """
    function performing all the required operations BEFORE the MLMC solution step
    input:  self: an instance of the class
    """
    def InitializeMLMCPhase(self):
        # compute tolerance
        self.ComputeTolerancei()
        # compute optimal number of levels
        self.ComputeLevels()
        # compute theta splitting parameter according to the current_number_levels and tolerance_i
        self.ComputeTheta(self.current_number_levels)
        # compute number of samples according to bayesian_variance and theta_i parameters
        self.ComputeNumberSamples()
        # prepare lists
        for _ in range (self.current_number_levels - self.previous_number_levels): # append a list for the new level
            self.difference_QoI.values.append([])
            self.time_ML.values.append([])

    """
    function performing all the required operations AFTER the MLMC solution step
    input:  self: an instance of the class
    """
    def FinalizeMLMCPhase(self):
        # prepare lists
        for _ in range (self.current_number_levels - self.previous_number_levels): # append a list for the new level
            self.difference_QoI.raw_moment_1.append([])
            self.difference_QoI.unbiased_central_moment_2.append([])
            self.difference_QoI.central_moment_1.append([])
            self.difference_QoI.central_moment_2.append([])
            self.difference_QoI.central_moment_3.append([])
            self.difference_QoI.central_moment_4.append([])
            self.difference_QoI.number_samples.append(0)
            self.time_ML.raw_moment_1.append([])
            self.time_ML.unbiased_central_moment_2.append([])
            self.time_ML.central_moment_1.append([])
            self.time_ML.central_moment_2.append([])
            self.time_ML.central_moment_3.append([])
            self.time_ML.central_moment_4.append([])
            self.time_ML.number_samples.append(0)
        # compute mean, sample variance and second moment for difference QoI and time ML
        for level in range (self.current_number_levels+1):
            for i_sample in range(self.previous_number_samples[level],self.number_samples[level]):
                self.difference_QoI.UpdateOnePassCentralMoments(level,i_sample)
                self.time_ML.UpdateOnePassCentralMoments(level,i_sample)
        # update number of levels
        self.previous_number_levels = self.current_number_levels
        """
        the functions executed in FinalizePhaseAux_Task are the followings:
        self.ComputeMeanMLMCQoI()
        self.ComputeRatesLS()
        self.EstimateBayesianVariance(self.current_number_levels)
        self.ComputeTotalErrorMLMC()
        """
        # store lists in synchro_element to execute FinalizePhaseAux_Task
        synchro_elements = [x for x in self.difference_QoI.raw_moment_1]
        synchro_elements.extend(["%%%"])
        synchro_elements.extend(self.difference_QoI.unbiased_central_moment_2)
        synchro_elements.extend(["%%%"])
        synchro_elements.extend(self.time_ML.raw_moment_1)
        # create a StreamSerializer Kratos object containing the problem settings
        serial_settings = KratosMultiphysics.StreamSerializer()
        serial_settings.Save("ParametersSerialization", self.settings)
        # compute remaining MLMC finalize process operations
        # observation: we are passing self.settings and we will exploit it to construct the class
        self.rates_error,self.bayesian_variance,self.mean_mlmc_QoI,\
        self.total_error,self.number_samples\
        = FinalizePhaseAux_Task(MultilevelMonteCarlo,
        serial_settings,self.mesh_parameters,self.current_number_levels,\
        self.current_iteration,self.number_samples,*synchro_elements)
        # synchronization point needed to compute the other functions of MLMC algorithm
        # put as in the end as possible the synchronization point
        self.difference_QoI.raw_moment_1 = get_value_from_remote(self.difference_QoI.raw_moment_1)
        self.difference_QoI.unbiased_central_moment_2 = get_value_from_remote(self.difference_QoI.unbiased_central_moment_2)
        self.time_ML.raw_moment_1 = get_value_from_remote(self.time_ML.raw_moment_1)
        self.total_error = get_value_from_remote(self.total_error)
        self.rates_error = get_value_from_remote(self.rates_error)
        self.bayesian_variance = get_value_from_remote(self.bayesian_variance)
        self.mean_mlmc_QoI = get_value_from_remote(self.mean_mlmc_QoI)
        self.number_samples = get_value_from_remote(self.number_samples)
        # check convergence
        # convergence reached if: i) current_iteration >= number_iterations_iE
        #                        ii) total_error < tolerance_i
        # if not update current_iteration
        if (self.total_error < self.settings["tolerance"].GetDouble()):
            self.convergence = True
        else:
            self.current_iteration = self.current_iteration + 1

    """
    function printing informations about screening phase
    input:  self: an instance of the class
    """
    def ScreeningInfoScreeningPhase(self):
        print("\n","#"*50," SCREENING PHASE ","#"*50,"\n")
        # print("values computed of QoI = ",self.difference_QoI.values)
        # print("values computed time_ML",self.time_ML.values)
        print("mean and variance difference_QoI = ",self.difference_QoI.raw_moment_1,self.difference_QoI.unbiased_central_moment_2)
        print("mean time_ML",self.time_ML.raw_moment_1)
        print("rates coefficient = ",self.rates_error)
        print("estimated Bayesian variance = ",self.bayesian_variance)
        print("estimated number of CMLMC iterations = ",self.number_iterations_iE)

    """
    function printing informations about initializing MLMC phase
    input:  self: an instance of the class
    """
    def ScreeningInfoInitializeMLMCPhase(self):
        print("\n","#"*50," MLMC iter =  ",self.current_iteration,"#"*50,"\n")
        print("current tolerance = ",self.tolerance_i)
        print("updated estimated Bayesian Variance initialize phase = ",self.bayesian_variance)
        print("current number of levels = ",self.current_number_levels)
        print("previous number of levels = ",self.previous_number_levels)
        print("current splitting parameter = ",self.theta_i)
        print("difference number of samples = ",self.difference_number_samples)
        print("previous number of samples = ",self.previous_number_samples)

    """
    function printing informations about finalizing MLMC phase
    input:  self: an instance of the class
    """
    def ScreeningInfoFinalizeMLMCPhase(self):
        # print("values computed of QoI = ",self.difference_QoI.values)
        # print("values computed time_ML",self.time_ML.values)
        print("current number of samples",self.number_samples)
        print("mean and variance difference_QoI = ",self.difference_QoI.raw_moment_1,self.difference_QoI.unbiased_central_moment_2)
        print("mean time_ML",self.time_ML.raw_moment_1)
        print("rates coefficient = ",self.rates_error)
        print("estimated Bayesian variance = ",self.bayesian_variance)
        print("multilevel monte carlo mean estimator = ",self.mean_mlmc_QoI)
        print("total_error = bias + statistical error = ",self.total_error)

    """
    function adding QoI and MLMC time values to the corresponding level
    input:  self:               an instance of the class
            simulation_results: tuple=(instance of MultilevelMonteCarloResults class, working level)
    """
    def AddResults(self,simulation_results):
        simulation_results_class = simulation_results[0]
        level = simulation_results[1]
        if (self.settings["store_lower_levels_samples"].GetString() == "True"):
            start_level = 0
        else:
            start_level = np.maximum(0,level)
        for lev in range (start_level,level+1):
            difference_QoI_value, time_ML_value = AddResultsAux_Task(simulation_results_class,lev)
            # update values of difference QoI and time ML per each level
            self.difference_QoI.values[lev] = np.append(self.difference_QoI.values[lev], difference_QoI_value)
            self.time_ML.values[lev] = np.append(self.time_ML.values[lev],time_ML_value)
            # update number of samples
            if (lev != level):
                self.number_samples[lev] = self.number_samples[lev] + 1

    """
    function giving as output the mesh discretization parameter
    the mesh parameter is the reciprocal of the minimum mesh size of the grid
    h_lev=h_0*M^(-lev)
    input:  self: an instance of the class
    """
    def ComputeMeshParameters(self):
        h0 = self.settings["initial_mesh_size"].GetDouble()
        M  = self.settings["mesh_refinement_coefficient"].GetInt()
        for level in range(self.settings["maximum_number_levels"].GetInt()+1):
            h_current_level = h0 * M**(-level)
            mesh_parameter_current_level = h_current_level**(-1)
            self.mesh_sizes.append(h_current_level)
            self.mesh_parameters.append(mesh_parameter_current_level)

    """
    function computing the problem parameters P=[calpha,alpha,cbeta,beta,cgamma,gamma] using least squares fit
    we consider level > 0 to compute calpha,alpha,cbeta,beta for robustness reasons [see [1] for details]
    input:  self: an instance of the class
    """
    def ComputeRatesLS(self):
        bias_ratesLS = np.abs(self.difference_QoI.raw_moment_1)
        variance_ratesLS = self.difference_QoI.unbiased_central_moment_2
        cost_ML_ratesLS = self.time_ML.raw_moment_1
        mesh_param_ratesLS = self.mesh_parameters[0:self.current_number_levels+1]
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
    input:  self:   an instance of the class
            levels: number of levels for which we want to estimate the Bayesian variance
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
        mean_local = copy.copy(self.difference_QoI.raw_moment_1)
        variance_local = copy.copy(self.difference_QoI.unbiased_central_moment_2)
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
    function computing the minimum number of iterations of MLMC algorithm
    iteration < number_iterations_iE : iterations needed to obtain increasingly accurate estimates of the
                                       problem dependent parameters P=[calpha,alpha,cbeta,beta,cgamma,gamma]
    iteration > number_iterations_iE : iterations preventing redundant computations due to fluctuations
                                       in the estimate of P=[calpha,alpha,cbeta,beta,cgamma,gamma]
                                       by solving the problem for a slightly smaller tolerance than the desired one
    input:  self: an instance of the class
    """
    def ComputeNumberIterationsMLMC(self):
        tolF = self.settings["tolerance"].GetDouble()
        tol0 = self.settings["initial_tolerance"].GetDouble()
        r2 = self.settings["r2"].GetDouble()
        r1 = self.settings["r1"].GetDouble()
        self.number_iterations_iE = np.floor((-np.log(tolF)+np.log(r2)+np.log(tol0))/(np.log(r1)))

    """
    function computing the tolerance for the i^th iteration
    input:  self: an instance of the class
    """
    def ComputeTolerancei(self):
        tol0 = self.settings["initial_tolerance"].GetDouble()
        r1 = self.settings["r1"].GetDouble()
        current_iter = self.current_iteration
        tol = tol0 / (r1**(current_iter))
        self.tolerance_i = tol

    """
    function computing the number of levels for i^th iteration of the algorithm
    input:  self:   an instance of the class
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
        self.previous_number_levels = Lmin
        del(tol,Wmin,current_number_levels,cgamma,gamma,calpha,alpha,cphi,mesh_parameters_all_levels,Lmax,Lmin)

    """
    function computing the splitting parameter theta \in (0,1)
    input:  self:  an instance of the class
            level: working level
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
    input:  self: an instance of the class
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
        nsam = self.number_samples
        # compute optimal number of samples and store previous number of samples
        coeff1 = (cphi/(theta*tol))**2.0
        model_cost = np.multiply(cgamma,np.power(mesh_parameters_current_levels,gamma))
        coeff2 = np.sqrt(np.divide(bayesian_variance,model_cost))
        coeff3 = np.sum(np.sqrt(np.multiply(model_cost,bayesian_variance)))
        opt_number_samples = np.multiply(coeff1*coeff3,coeff2)
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
        self.number_samples = nsam
        self.difference_number_samples = diff_nsam
        self.previous_number_samples = previous_number_samples
        del(current_number_levels,bayesian_variance,min_samples_add,cgamma,gamma,cphi,mesh_parameters_current_levels,\
        theta,tol,nsam,opt_number_samples,previous_number_samples,diff_nsam)

    """
    function computing the mlmc estimator for the mean of the Quantity of Interest
    input:  self: an instance of the class
    """
    def ComputeMeanMLMCQoI(self):
        self.mean_mlmc_QoI = np.sum(self.difference_QoI.raw_moment_1)

    """
    function computing the total error:
    total_error = bias contrinution + statistical error contribution
    input:  self: an instance of the class
    """
    def ComputeTotalErrorMLMC(self):
        self.difference_QoI.bias_error = np.abs(self.difference_QoI.raw_moment_1[self.current_number_levels])
        variance_from_bayesian = np.zeros(np.size(self.number_samples))
        for lev in range(self.current_number_levels+1):
            variance_from_bayesian[lev] = self.bayesian_variance[lev]/self.number_samples[lev]
        self.difference_QoI.statistical_error = self.settings["cphi_confidence"].GetDouble() * np.sqrt(np.sum(variance_from_bayesian))
        total_error = self.difference_QoI.bias_error + self.difference_QoI.statistical_error
        self.total_error = total_error


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