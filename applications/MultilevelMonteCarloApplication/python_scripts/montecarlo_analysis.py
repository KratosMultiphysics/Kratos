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

# Import Monte Carlo library
import mc_utilities as mc

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


'''Adapt the following class depending on the problem, deriving the MonteCarloAnalysis class from the problem of interest'''

'''
This Analysis Stage implementation solves the elliptic PDE in (0,1)^2 with zero Dirichlet boundary conditions
-lapl(u) = xi*f,    f= -432*x*(x-1)*y*(y-1)
                    f= -432*(x**2+y**2-x-y)
where xi is a Beta(2,6) random variable, and computes statistic of the QoI
Q = int_(0,1)^2 u(x,y)dxdy
'''
class MonteCarloAnalysis(AnalysisStage):
    '''Main analysis stage for Monte Carlo simulations'''
    def __init__(self,input_model,input_parameters,sample):
        self.sample = sample
        super(MonteCarloAnalysis,self).__init__(input_model,input_parameters)
        self._GetSolver().main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)

    def _CreateSolver(self):
        import convection_diffusion_stationary_solver
        return convection_diffusion_stationary_solver.CreateSolver(self.model,self.project_parameters["solver_settings"])

    def _GetSimulationName(self):
        return "Monte Carlo Analysis"

    '''Introduce here the stochasticity in the right hand side defining the forcing function and apply the stochastic contribute'''
    def ModifyInitialProperties(self):
        for node in self.model.GetModelPart("MLMCLaplacianModelPart").Nodes:
            coord_x = node.X
            coord_y = node.Y
            # forcing = -432.0 * coord_x * (coord_x - 1) * coord_y * (coord_y - 1)
            forcing = -432.0 * (coord_x**2 + coord_y**2 - coord_x - coord_y) # this forcing presents the below commented analytical solution
            node.SetSolutionStepValue(KratosMultiphysics.HEAT_FLUX,forcing*self.sample)


##################################################
######## END OF CLASS MONTECARLOANALYSIS #########
##################################################


'''
function generating the random sample
here the sample has a beta distribution with parameters alpha = 2.0 and beta = 6.0
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
    """here we evaluate the QoI of the problem: int_{domain} SOLUTION(x,y) dx dy
    we use the midpoint rule to evaluate the integral"""
    KratosMultiphysics.CalculateNodalAreaProcess(simulation._GetSolver().main_model_part,2).Execute()
    Q = 0.0
    for node in simulation._GetSolver().main_model_part.Nodes:
        Q = Q + (node.GetSolutionStepValue(KratosMultiphysics.NODAL_AREA)*node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE))
    return Q


'''
function called in the main returning a future object (the result class) and an integer (the finer level)
input:
        pickled_coarse_model      : pickled model
        pickled_coarse_parameters : pickled parameters
output:
        MonteCarloResults class   : class of the simulation results
        current_MC_level          : level of the current MLMC simulation
'''
def ExecuteMonteCarloAnalysis(pickled_model, pickled_parameters):
    current_MC_level = 0 # MC has only level 0
    return (ExecuteMonteCarloAnalysis_Task(pickled_model, pickled_parameters),current_MC_level)


'''
function executing the problem
input:
        model       : serialization of the model
        parameters  : serialization of the Project Parameters
output:
        QoI         : Quantity of Interest
'''
@ExaquteTask(returns=1)
def ExecuteMonteCarloAnalysis_Task(pickled_model, pickled_parameters):
    '''overwrite the old model serializer with the unpickled one'''
    model_serializer = pickle.loads(pickled_model)
    current_model = KratosMultiphysics.Model()
    model_serializer.Load("ModelSerialization",current_model)
    del(model_serializer)
    '''overwrite the old parameters serializer with the unpickled one'''
    serialized_parameters = pickle.loads(pickled_parameters)
    current_parameters = KratosMultiphysics.Parameters()
    serialized_parameters.Load("ParametersSerialization",current_parameters)
    del(serialized_parameters)
    '''initialize the MonteCarloResults class'''
    current_level = 0 # always 0 for MC
    mc_results_class = mc.MonteCarloResults(current_level)
    sample = GenerateSample()
    simulation = MonteCarloAnalysis(current_model,current_parameters,sample)
    simulation.Run()
    QoI = EvaluateQuantityOfInterest(simulation)
    mc_results_class.QoI[current_level].append(QoI) # saving results in the corresponding list, for MC only list of level 0
    return mc_results_class


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
    # local_parameters["solver_settings"]["model_import_settings"]["input_filename"].SetString(model_part_file_name[:-5])
    fake_sample = GenerateSample()
    simulation = MonteCarloAnalysis(model,local_parameters,fake_sample)
    simulation.Initialize()
    serialized_model = KratosMultiphysics.StreamSerializer()
    serialized_model.Save("ModelSerialization",simulation.model)
    serialized_parameters = KratosMultiphysics.StreamSerializer()
    serialized_parameters.Save("ParametersSerialization",simulation.project_parameters)
    # pickle dataserialized_data
    pickled_model = pickle.dumps(serialized_model, 2) # second argument is the protocol and is NECESSARY (according to pybind11 docs)
    pickled_parameters = pickle.dumps(serialized_parameters, 2)
    print("\n","#"*50," SERIALIZATION COMPLETED ","#"*50,"\n")
    return pickled_model,pickled_parameters


if __name__ == '__main__':

    '''set the ProjectParameters.json path'''
    parameter_file_name = "../tests/PoissonSquareTest/parameters_poisson_finer.json"
    '''create a serialization of the model and of the project parameters'''
    pickled_model,pickled_parameters = SerializeModelParameters_Task(parameter_file_name)
    '''customize setting parameters of the ML simulation'''
    settings_MC_simulation = KratosMultiphysics.Parameters("""
    {
        "tolerance" : 0.1,
        "cphi" : 5e-1,
        "batch_size" : 20,
        "convergence_criteria" : "MC_higher_moments_sequential_stopping_rule"
    }
    """)
    '''contruct MonteCarlo class'''
    mc_class = mc.MonteCarlo(settings_MC_simulation)
    '''start MC algorithm'''
    while mc_class.convergence is not True:
        mc_class.InitializeMCPhase()
        mc_class.ScreeningInfoInitializeMCPhase()
        for instance in range (mc_class.difference_number_samples[0]):
            mc_class.AddResults(ExecuteMonteCarloAnalysis(pickled_model,pickled_parameters))
        mc_class.FinalizeMCPhase()
        mc_class.ScreeningInfoFinalizeMCPhase()

    mc_class.QoI.mean = get_value_from_remote(mc_class.QoI.mean)
    print("\nMC mean = ",mc_class.QoI.mean)


    ''' The below part evaluates the relative L2 error between the numerical solution SOLUTION(x,y,sample) and the analytical solution, also dependent on sample.
    Analytical solution available in case FORCING = sample * -432.0 * (coord_x**2 + coord_y**2 - coord_x - coord_y)'''
    # model_serializer = pickle.loads(pickled_model)
    # current_model = KratosMultiphysics.Model()
    # model_serializer.Load("ModelSerialization",current_model)
    # del(model_serializer)
    # serialized_parameters = pickle.loads(pickled_parameters)
    # current_parameters = KratosMultiphysics.Parameters()
    # serialized_parameters.Load("ParametersSerialization",current_parameters)
    # del(serialized_parameters)
    # sample = 1.0
    # simulation = MonteCarloAnalysis(current_model,current_parameters,sample)
    # simulation.Run()
    # KratosMultiphysics.CalculateNodalAreaProcess(simulation._GetSolver().main_model_part,2).Execute()
    # error = 0.0
    # L2norm_analyticalsolution = 0.0
    # for node in simulation._GetSolver().main_model_part.Nodes:
    #     local_error = ((node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE) - (432.0*simulation.sample*node.X*node.Y*(1-node.X)*(1-node.Y)*0.5))**2) * node.GetSolutionStepValue(KratosMultiphysics.NODAL_AREA)
    #     error = error + local_error
    #     local_analyticalsolution = (432.0*simulation.sample*node.X*node.Y*(1-node.X)*(1-node.Y)*0.5)**2 * node.GetSolutionStepValue(KratosMultiphysics.NODAL_AREA)
    #     L2norm_analyticalsolution = L2norm_analyticalsolution + local_analyticalsolution
    # error = np.sqrt(error)
    # L2norm_analyticalsolution = np.sqrt(L2norm_analyticalsolution)
    # print("L2 relative error = ", error/L2norm_analyticalsolution)