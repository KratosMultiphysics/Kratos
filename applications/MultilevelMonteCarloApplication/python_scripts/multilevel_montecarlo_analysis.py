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
import cmlmc_utilities as mlmc

# Import refinement library
import adaptive_refinement_utilities as refinement

# Import cpickle to pickle the serializer
try:
    import cpickle as pickle  # Use cPickle on Python 2.7
except ImportError:
    import pickle

# Import exaqute
from exaqute.ExaquteTaskPyCOMPSs import *   # to execute with pycompss
# from exaqute.ExaquteTaskHyperLoom import *  # to execute with the IT4 scheduler
# from exaqute.ExaquteTaskLocal import *      # to execute with python3
'''
get_value_from_remote is the equivalent of compss_wait_on: a synchronization point
in future, when everything is integrated with the it4i team, importing exaqute.ExaquteTaskHyperLoom you can launch your code with their scheduler instead of BSC
'''


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
            forcing = -432.0 * (coord_x**2 + coord_y**2 - coord_x - coord_y)
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
    QoI = 0.0
    for node in simulation._GetSolver().main_model_part.Nodes:
        QoI = QoI + (node.GetSolutionStepValue(KratosMultiphysics.NODAL_AREA)*node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE))
    return QoI


'''
function called in the main returning a future object (the result class) and an integer (the finer level)
input:
        current_MLMC_level        : current Multilevel Monte Carlo level we are solving
        pickled_coarse_model      : pickled model
        pickled_coarse_parameters : pickled parameters
        size_meshes               : mesh sizes for all levels
output:
        MultilevelMonteCarloResults class : class of the simulation results
        current_MLMC_level                      : level of the current MLMC simulation
'''
def ExecuteMultilevelMonteCarloAnalisys(current_MLMC_level,pickled_coarse_model,pickled_coarse_parameters,size_meshes,pickled_settings_metric_refinement,pickled_settings_remesh_refinement):
    '''generate the sample'''
    sample = GenerateSample()
    '''initialize the MultilevelMonteCarloResults class'''
    mlmc_results_class = mlmc.MultilevelMonteCarloResults(current_MLMC_level)
    if (current_MLMC_level == 0):
        mlmc_results_class,pickled_current_model,pickled_current_parameters = ExecuteMultilevelMonteCarloAnalisys_Task(current_MLMC_level,pickled_coarse_model,pickled_coarse_parameters,size_meshes,pickled_settings_metric_refinement,pickled_settings_remesh_refinement,sample,current_MLMC_level,mlmc_results_class)
    else:
        for level in range(current_MLMC_level+1):
            mlmc_results_class,pickled_current_model,pickled_current_parameters = ExecuteMultilevelMonteCarloAnalisys_Task(current_MLMC_level,pickled_coarse_model,pickled_coarse_parameters,size_meshes,pickled_settings_metric_refinement,pickled_settings_remesh_refinement,sample,level,mlmc_results_class)
            del(pickled_coarse_model,pickled_coarse_parameters)
            pickled_coarse_model = pickled_current_model
            pickled_coarse_parameters = pickled_current_parameters
            del(pickled_current_model,pickled_current_parameters)
    return mlmc_results_class,current_MLMC_level


'''
function evaluating the QoI and the cost of simulation, computing the mesh of level current_MLMC_level
refining recursively from the coarsest mesh
input:
        current_MLMC_level        : current Multilevel Monte Carlo level we are solving
        pickled_coarse_model      : pickled model
        pickled_coarse_parameters : pickled parameters
        size_meshes               : mesh sizes for all levels
output:
        mlmc_results_class : QoI         : list of QoI for all levels computed in the current simulation
                             finer_level : finest level
                             time_ML     : list of MLMC time for all levels computed in the current simulation
'''
@ExaquteTask(returns=3)
def ExecuteMultilevelMonteCarloAnalisys_Task(current_MLMC_level,pickled_coarse_model,pickled_coarse_parameters,size_meshes,pickled_settings_metric_refinement,pickled_settings_remesh_refinement,sample,current_level,mlmc_results_class):
    '''unpickle model and build Kratos Model object'''
    serialized_model = pickle.loads(pickled_coarse_model)
    current_model = KratosMultiphysics.Model()
    serialized_model.Load("ModelSerialization",current_model)
    del(serialized_model)
    '''unpickle parameters and build Kratos Parameters object'''
    serialized_parameters = pickle.loads(pickled_coarse_parameters)
    current_parameters = KratosMultiphysics.Parameters()
    serialized_parameters.Load("ParametersSerialization",current_parameters)
    del(serialized_parameters)
    '''start time'''
    start_MLMC_time = time.time()
    '''refine if current current_level > 0, adaptive refinement based on the solution of previous level'''
    if (current_level > 0):
        '''unpickle metric and remesh refinement parameters and build Kratos Parameters objects'''
        settings_metric_refinement_serializer = pickle.loads(pickled_settings_metric_refinement)
        settings_remesh_refinement_serializer = pickle.loads(pickled_settings_remesh_refinement)
        current_settings_metric_refinement = KratosMultiphysics.Parameters()
        current_settings_remesh_refinement = KratosMultiphysics.Parameters()
        settings_metric_refinement_serializer.Load("MetricRefinementParametersSerialization",current_settings_metric_refinement)
        settings_remesh_refinement_serializer.Load("RemeshRefinementParametersSerialization",current_settings_remesh_refinement)
        del(settings_metric_refinement_serializer,settings_remesh_refinement_serializer)
        '''refine the model Kratos object'''
        refined_model,refined_parameters = refinement.compute_refinement_hessian_metric(current_model,current_parameters,size_meshes[current_level],size_meshes[current_level-1],current_settings_metric_refinement,current_settings_remesh_refinement)
        '''initialize the model Kratos object'''
        simulation = MultilevelMonteCarloAnalysis(refined_model,refined_parameters,sample)
        simulation.Initialize()
        '''update model Kratos object'''
        current_model = simulation.model
        current_parameters = simulation.project_parameters
        del(simulation)
    simulation = MultilevelMonteCarloAnalysis(current_model,current_parameters,sample)
    simulation.Run()
    QoI = EvaluateQuantityOfInterest(simulation)
    '''save model and parameters as StreamSerializer Kratos objects'''
    serialized_finer_model = KratosMultiphysics.StreamSerializer()
    serialized_finer_model.Save("ModelSerialization",simulation.model)
    serialized_finer_parameters = KratosMultiphysics.StreamSerializer()
    serialized_finer_parameters.Save("ParametersSerialization",simulation.project_parameters)
    '''pickle model and parameters'''
    pickled_finer_model = pickle.dumps(serialized_finer_model, 2) # second argument is the protocol and is NECESSARY (according to pybind11 docs)
    pickled_finer_parameters = pickle.dumps(serialized_finer_parameters, 2) # second argument is the protocol and is NECESSARY (according to pybind11 docs)
    del(simulation)
    end_MLMC_time = time.time()
    '''register results of the current level in the MultilevelMonteCarloResults class'''
    mlmc_results_class.time_ML[current_level].append(end_MLMC_time-start_MLMC_time) # saving each result in the corresponding list in order to ensure the correctness of the results order and the levels
    mlmc_results_class.QoI[current_level].append(QoI) # saving each result in the corresponding list in order to ensure the correctness of the results order and the levels
    return mlmc_results_class,pickled_finer_model,pickled_finer_parameters


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
@ExaquteTask(parameter_file_name=FILE_IN,returns=2)
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
#TODO: to execute this function in a compss task, metric_refinement_parameters and remeshing_refinement_parameters should be read from a file
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


if __name__ == '__main__':

    '''set the ProjectParameters.json path'''
    parameter_file_name = "../tests/PoissonSquareTest/parameters_poisson_coarse.json"
    '''create a serialization of the model and of the project parameters'''
    pickled_model,pickled_parameters = SerializeModelParameters_Task(parameter_file_name)
    '''customize setting parameters of the ML simulation'''
    settings_ML_simulation = KratosMultiphysics.Parameters("""
    {
        "tol0"                            : 0.25,
        "tolF"                            : 0.1,
        "cphi"                            : 1.0,
        "number_samples_screening"        : 25,
        "Lscreening"                      : 2,
        "Lmax"                            : 4,
        "initial_mesh_size"               : 0.5
    }
    """)
    '''customize setting parameters of the metric of the adaptive refinement utility'''
    settings_metric_refinement = KratosMultiphysics.Parameters("""
        {
            "hessian_strategy_parameters"           :{
                    "metric_variable"               : ["TEMPERATURE"],
                    "estimate_interpolation_error"  : false,
                    "interpolation_error"           : 0.004
            },
            "anisotropy_remeshing"                  : true,
            "anisotropy_parameters":{
                "reference_variable_name"           : "TEMPERATURE",
                "hmin_over_hmax_anisotropic_ratio"  : 0.15,
                "boundary_layer_max_distance"       : 1.0,
                "interpolation"                     : "Linear"
            },
            "local_gradient_variable"               : "TEMPERATURE"
        }
    """)
    '''customize setting parameters of the remesh of the adaptive refinement utility'''
    settings_remesh_refinement = KratosMultiphysics.Parameters("""
        {
            "echo_level"                            : 0
        }
    """)
    pickled_settings_metric_refinement,pickled_settings_remesh_refinement = SerializeRefinementParameters(settings_metric_refinement,settings_remesh_refinement)

    '''contruct MultilevelMonteCarlo class'''
    mlmc_class = mlmc.MultilevelMonteCarlo(settings_ML_simulation)
    ''''start screening phase'''
    for level in range(mlmc_class.current_number_levels+1):
        for instance in range (mlmc_class.difference_number_samples[level]):
            mlmc_class.AddResults(ExecuteMultilevelMonteCarloAnalisys(level,pickled_model,pickled_parameters,mlmc_class.sizes_mesh,pickled_settings_metric_refinement,pickled_settings_remesh_refinement))
    '''finalize screening phase'''
    mlmc_class.FinalizeScreeningPhase()
    mlmc_class.ScreeningInfoScreeningPhase()
    '''start MLMC phase'''
    while mlmc_class.convergence is not True:
        '''initialize MLMC phase'''
        mlmc_class.InitializeMLMCPhase()
        mlmc_class.ScreeningInfoInitializeMLMCPhase()
        '''MLMC execution phase'''
        for level in range (mlmc_class.current_number_levels+1):
            for instance in range (mlmc_class.difference_number_samples[level]):
                mlmc_class.AddResults(ExecuteMultilevelMonteCarloAnalisys(level,pickled_model,pickled_parameters,mlmc_class.sizes_mesh,pickled_settings_metric_refinement,pickled_settings_remesh_refinement))
        '''finalize MLMC phase'''
        mlmc_class.FinalizeMLMCPhase()
        mlmc_class.ScreeningInfoFinalizeMLMCPhase()

    print("\niterations = ",mlmc_class.current_iteration,\
    "total error TErr computed = ",mlmc_class.TErr,"mean MLMC QoI = ",mlmc_class.mean_mlmc_QoI)

    '''### OBSERVATIONS ###

    compss: set absolute path when launching with compss

    MmgProcess: need to use conditions in the model part to preserve the boundary conditions in the refinement process
    The submodelpart: Subpart_Boundary contains only nodes and no geometries (conditions/elements).
    It is not guaranteed that the submodelpart will be preserved.
    PLEASE: Add some "dummy" conditions to the submodelpart to preserve it'''