from __future__ import absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.ConvectionDiffusionApplication as KratosConvDiff
import KratosMultiphysics.MultilevelMonteCarloApplication as KratosMLMC

# Avoid printing of Kratos informations
KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING) # avoid printing of Kratos things

# Importing the base class
from analysis_stage import AnalysisStage

# Import packages
import numpy as np
import time

# Import Continuation Multilevel Monte Carlo library
import test_cmlmc_utilities as mlmc

# Import refinement library
import adaptive_refinement_utilities as refinement

# Import cpickle to pickle the serializer
try:
    import cpickle as pickle  # Use cPickle on Python 2.7
except ImportError:
    import pickle


'''Adapt the following class depending on the problem, deriving the MultilevelMonteCarloAnalysis class from the problem of interest'''

'''
This Analysis Stage implementation solves the elliptic PDE in (0,1)^2 with zero Dirichlet boundary conditions
-lapl(u) = xi*f,    f= -432*x*(x-1)*y*(y-1)
                    f= -432*(x**2+y**2-x-y)
where xi is a Beta(2,6) random variable, and computes statistic of the QoI
Q = int_(0,1)^2 u(x,y)dxdy
'''


class MultilevelMonteCarloAnalysis(AnalysisStage):
    '''Main analysis stage for MultilevelMonte Carlo simulations'''
    def __init__(self,input_model,input_parameters,sample):
        self.sample = sample
        super(MultilevelMonteCarloAnalysis,self).__init__(input_model,input_parameters)
        self._GetSolver().main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)

    def _CreateSolver(self):
        import convection_diffusion_stationary_solver
        return convection_diffusion_stationary_solver.CreateSolver(self.model,self.project_parameters["solver_settings"])

    def _GetSimulationName(self):
        return "Multilevel Monte Carlo Analysis"

    '''Introduce here the stochasticity in the right hand side defining the forcing function and apply the stochastic contribute'''
    def ModifyInitialProperties(self):
        for node in self.model.GetModelPart("MLMCLaplacianModelPart").Nodes:
            coord_x = node.X
            coord_y = node.Y
            # forcing = -432.0 * coord_x * (coord_x - 1) * coord_y * (coord_y - 1)
            forcing = -432.0 * (coord_x**2 + coord_y**2 - coord_x - coord_y) # this forcing presents an analytical solution
            node.SetSolutionStepValue(KratosMultiphysics.HEAT_FLUX,forcing*self.sample)


'''
function generating the random sample
'''
def GenerateSample():
    alpha = 2.0
    beta = 6.0
    number_samples = 1
    sample = np.random.beta(alpha,beta,number_samples)
    return sample


'''
function evaluating the QoI of the problem: int_{domain} TEMPERATURE(x,y) dx dy
right now we are using the midpoint rule to evaluate the integral: improve!
'''
def EvaluateQuantityOfInterest(simulation):
    '''TODO: in future once is ready: declare and initialize NODAL_AREA as non historical variable:
    this way I only store the values in the current nodes and not also in the buffer'''
    # KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.NODAL_AREA, 0.0, simulation._GetSolver().main_model_part.Nodes)
    KratosMultiphysics.CalculateNodalAreaProcess(simulation._GetSolver().main_model_part,2).Execute()
    Q = 0.0
    for node in simulation._GetSolver().main_model_part.Nodes:
        Q = Q + (node.GetSolutionStepValue(KratosMultiphysics.NODAL_AREA)*node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE))
        # print(node.Id,"NODAL AREA = ",node.GetSolutionStepValue(KratosMultiphysics.NODAL_AREA),"NODAL SOLUTION = ",node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE),"CURRENT Q = ",Q)
    return Q


'''
function called in the main returning a future object (the result class) and an integer (the finer level)
input:
        finest_level              : current Multilevel MOnte Carlo level we are solving
        pickled_coarse_model      : pickled model
        pickled_coarse_parameters : pickled parameters
        size_meshes               : mesh sizes for all levels
output:
        MultilevelMonteCarloResults class : class of the simulation results
        finest_level                      : level of the current MLMC simulation
'''
def ExecuteMultilevelMonteCarloAnalisys(finest_level,pickled_coarse_model,pickled_coarse_parameters,size_meshes,pickled_settings_metric_refinement,pickled_settings_remesh_refinement):
    return ExecuteMultilevelMonteCarloAnalisys_Task(finest_level,pickled_coarse_model,pickled_coarse_parameters,size_meshes,pickled_settings_metric_refinement,pickled_settings_remesh_refinement),finest_level


'''
function evaluating the QoI and the cost of simulation, computing the mesh of level finest_level
refining recursively from the coarsest mesh
input:
        finest_level              : current Multilevel MOnte Carlo level we are solving
        pickled_coarse_model      : pickled model
        pickled_coarse_parameters : pickled parameters
        size_meshes               : mesh sizes for all levels
output:
        mlmc_results_class : QoI         : list of QoI for all levels computed in the current simulation
                             finer_level : finest level
                             time_ML     : list of MLMC time for all levels computed in the current simulation
'''
def ExecuteMultilevelMonteCarloAnalisys_Task(finest_level,pickled_coarse_model,pickled_coarse_parameters,size_meshes,pickled_settings_metric_refinement,pickled_settings_remesh_refinement):
    '''unpickle model and build Kratos Model object'''
    model_serializer = pickle.loads(pickled_coarse_model)
    current_model = KratosMultiphysics.Model()
    model_serializer.Load("ModelSerialization",current_model)
    del(model_serializer)
    '''unpickle parameters and build Kratos Parameters object'''
    serialized_parameters = pickle.loads(pickled_coarse_parameters)
    current_parameters = KratosMultiphysics.Parameters()
    serialized_parameters.Load("ParametersSerialization",current_parameters)
    del(serialized_parameters)
    '''unpickle metric and remesh refinement parameters and build Kratos Parameters objects'''
    settings_metric_refinement_serializer = pickle.loads(pickled_settings_metric_refinement)
    settings_remesh_refinement_serializer = pickle.loads(pickled_settings_remesh_refinement)
    current_settings_metric_refinement = KratosMultiphysics.Parameters()
    current_settings_remesh_refinement = KratosMultiphysics.Parameters()
    settings_metric_refinement_serializer.Load("MetricRefinementParametersSerialization",current_settings_metric_refinement)
    settings_remesh_refinement_serializer.Load("RemeshRefinementParametersSerialization",current_settings_remesh_refinement)
    del(settings_metric_refinement_serializer,settings_remesh_refinement_serializer)
    '''generate the sample'''
    sample = GenerateSample()
    '''initialize the MultilevelMonteCarloResults class and prepare the results'''
    mlmc_results_class = mlmc.MultilevelMonteCarloResults()
    QoI = []
    start_MLMC_time = time.time()
    end_MLMC_time = []
    if(finest_level == 0):
        simulation = MultilevelMonteCarloAnalysis(current_model,current_parameters,sample)
        simulation.Run()
        QoI.append(EvaluateQuantityOfInterest(simulation))
        del(simulation)
        end_MLMC_time.append(time.time())
    else:
        for lev in range(finest_level+1):
            simulation = MultilevelMonteCarloAnalysis(current_model,current_parameters,sample)
            simulation.Run()
            QoI.append(EvaluateQuantityOfInterest(simulation))
            end_MLMC_time.append(time.time())
            '''refine if level < finest level exploiting the solution just computed'''
            if (lev < finest_level):
                '''refine the model Kratos object'''
                model_refined = refinement.compute_refinement_hessian_metric(simulation,size_meshes[lev+1],size_meshes[lev],current_settings_metric_refinement,current_settings_remesh_refinement)
                '''initialize the model Kratos object'''
                simulation = MultilevelMonteCarloAnalysis(model_refined,current_parameters,sample)
                simulation.Initialize()
                '''update model Kratos object'''
                current_model = simulation.model
            del(simulation)
    '''prepare results of the simulation in the MultilevelMonteCarloResults class'''
    mlmc_results_class.finer_level = finest_level
    for lev in range(finest_level+1):
        mlmc_results_class.time_ML.append(end_MLMC_time[lev]-start_MLMC_time)
        mlmc_results_class.QoI.append(QoI[lev])
    return mlmc_results_class


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
        pickled_model      : model serializaton
        pickled_parameters : project parameters serialization
'''
def SerializeModelParameters_Task(parameter_file_name):
    with open(parameter_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    local_parameters = parameters
    model = KratosMultiphysics.Model()
    fake_sample = GenerateSample()
    '''initialize'''
    simulation = MultilevelMonteCarloAnalysis(model,local_parameters,fake_sample)
    simulation.Initialize()
    '''save model and parameters as StreamSerializer Kratos objects'''
    serialized_model = KratosMultiphysics.StreamSerializer()
    serialized_model.Save("ModelSerialization",simulation.model)
    serialized_parameters = KratosMultiphysics.StreamSerializer()
    serialized_parameters.Save("ParametersSerialization",simulation.project_parameters)
    '''pickle model and parameters'''
    pickled_model = pickle.dumps(serialized_model, 2) # second argument is the protocol and is NECESSARY (according to pybind11 docs)
    pickled_parameters = pickle.dumps(serialized_parameters, 2) # second argument is the protocol and is NECESSARY (according to pybind11 docs)
    print("\n","#"*50," SERIALIZATION MODEL AND PARAMETERS COMPLETED ","#"*50,"\n")
    return pickled_model,pickled_parameters

'''
function serializing and pickling the metric and remeshing parameters of the problem
the idea is the following:
i)   from Parameters Kratos object to StreamSerializer Kratos object
ii)  from StreamSerializer Kratos object to pickle string
iii) from pickle string to StreamSerializer Kratos object
iv)  from StreamSerializer Kratos object to Parameters Kratos object
input:
        metric_refinement_parameters    : Kratos Parameters object
        remeshing_refinement_parameters : Kratos Parameters object
output:
        pickled_metric_refinement_parameters    : project parameters serialization
        pickled_remeshing_refinement_parameters : project parameters serialization
'''
def SerializeRefinementParameters(metric_refinement_parameters,remeshing_refinement_parameters):
    '''save parameters as StreamSerializer Kratos objects'''
    serialized_metric_refinement_parameters = KratosMultiphysics.StreamSerializer()
    serialized_metric_refinement_parameters.Save("MetricRefinementParametersSerialization",metric_refinement_parameters)
    serialized_remesh_refinement_parameters = KratosMultiphysics.StreamSerializer()
    serialized_remesh_refinement_parameters.Save("RemeshRefinementParametersSerialization",remeshing_refinement_parameters)
    '''pickle parameters'''
    pickled_metric_refinement_parameters = pickle.dumps(serialized_metric_refinement_parameters, 2) # second argument is the protocol and is NECESSARY (according to pybind11 docs)
    pickled_remesh_refinement_parameters = pickle.dumps(serialized_remesh_refinement_parameters, 2) # second argument is the protocol and is NECESSARY (according to pybind11 docs)
    print("\n","#"*50," SERIALIZATION REFINEMENT PARAMETERS COMPLETED ","#"*50,"\n")
    return pickled_metric_refinement_parameters,pickled_remesh_refinement_parameters

'''
function executing the refinement of the problem
input:
        pickled_model_coarse : serialization of the model with coarser model part
        pickled_parameters   : serialization of the Project Parameters
        min_size             : minimum size of the refined model part
        max_size             : maximum size of the refined mesh
output:
        QoI                   : Quantity of Interest
        pickled_model_refined : serialization of the model with refined model part
'''
def ExecuteRefinement_Task(pickled_model_coarse, pickled_parameters, min_size, max_size):
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
    simulation_coarse = MultilevelMonteCarloAnalysis(model_coarse,parameters_refinement,fake_sample)
    simulation_coarse.Run()
    QoI =  EvaluateQuantityOfInterest(simulation_coarse)
    '''refine'''
    model_refined = refinement.compute_refinement_hessian_metric(simulation_coarse,min_size,max_size)
    '''initialize'''
    simulation = MultilevelMonteCarloAnalysis(model_refined,parameters_refinement,fake_sample)
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
def CompareMean_Task(AveragedMeanQoI,ExactExpectedValueQoI):
    relative_error = abs((AveragedMeanQoI - ExactExpectedValueQoI)/ExactExpectedValueQoI)
    return relative_error