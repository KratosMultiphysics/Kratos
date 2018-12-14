from __future__ import absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import numpy as np
import time

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.ConvectionDiffusionApplication as KratosConvDiff
import KratosMultiphysics.MultilevelMonteCarloApplication as KratosMLMC

# Avoid printing of Kratos informations
KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING) # avoid printing of Kratos things

# Importing the base class
from analysis_stage import AnalysisStage

# Import exaqute
from exaqute.ExaquteTaskPyCOMPSs import *   # to exequte with pycompss
# from exaqute.ExaquteTaskHyperLoom import *  # to exequte with the IT4 scheduler
# from exaqute.ExaquteTaskLocal import *      # to execute with python3
'''
get_value_from_remote is the equivalent of compss_wait_on: a synchronization point
in future, when everything is integrated with the it4i team, importing exaqute.ExaquteTaskHyperLoom you can launch your code with their scheduler instead of BSC
'''

# Import Continuation Multilevel Monte Carlo library
import cmlmc_utilities as mlmc

# Import refinement library
import adaptive_refinement_utilities as refinement

# Import cpickle to pickle the serializer
try:
    import cpickle as pickle  # Use cPickle on Python 2.7
except ImportError:
    import pickle


'''Adapt the following class depending on the problem, deriving the MultilevelMonteCarloAnalysis class from the problem of interest'''

'''This Analysis Stage implementation solves the elliptic PDE in (0,1)^2 with zero Dirichlet boundary conditions
-lapl(u) = xi*f,    f= -432*x*(x-1)*y*(y-1)
                    f= -432*(x**2+y**2-x-y)
where xi is a Beta(2,6) random variable, and computes statistic of the QoI
Q = int_(0,1)^2 u(x,y)dxdy
more details in Section 5.2 of [PKN17]

References:
[PKN17] M. Pisaroni; S. Krumscheid; F. Nobile : Quantifying uncertain system outputs via the multilevel Monte Carlo method - Part I: Central moment estimation; MATHICSE technical report no. 23.2017.
'''


class MultilevelMonteCarloAnalysis(AnalysisStage):
    '''Main analysis stage for MultilevelMonte Carlo simulations'''
    def __init__(self,input_model,input_parameters,sample):
        self.sample = sample
        super(MultilevelMonteCarloAnalysis,self).__init__(input_model,input_parameters)
        self._GetSolver().main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)

    def _CreateSolver(self):
        import convection_diffusion_stationary_solver
        solver = convection_diffusion_stationary_solver.CreateSolver(self.model,self.project_parameters["solver_settings"])
        self.LaplacianSolver = solver
        return self.LaplacianSolver
    
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



###########################################################
######## END OF CLASS MULTILEVELMONTECARLOANALYSIS ########
###########################################################


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
    '''TODO: in future once is ready: declare and initialize NODAL_AREA as non historical variable:
    this way I only store the values in the current nodes and not also in the buffer'''
    # KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.NODAL_AREA, 0.0, simulation._GetSolver().main_model_part.Nodes)
    KratosMultiphysics.CalculateNodalAreaProcess(simulation._GetSolver().main_model_part,2).Execute()
    Q = 0.0
    for node in simulation._GetSolver().main_model_part.Nodes:
        Q = Q + (node.GetSolutionStepValue(KratosMultiphysics.NODAL_AREA)*node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE))
        # print(node.Id,"NODAL AREA = ",node.GetSolutionStepValue(KratosMultiphysics.NODAL_AREA),"NODAL SOLUTION = ",node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE),"CURRENT Q = ",Q)
    return Q


# '''this function returns a future object (the result class) and a master integer (the finer level)
# this would avoid to the user to see lev in the call in the main as an input for the AddResults function'''
# def ExecuteMultilevelMonteCarloAnalisys(finest_level,pickled_coarse_model,pickled_coarse_parameters,size_meshes):
#     return ExecuteMultilevelMonteCarloAnalisysAux(finest_level,pickled_coarse_model,pickled_coarse_parameters,size_meshes),finest_level


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
def ExecuteMultilevelMonteCarloAnalisys_Task(finest_level,pickled_coarse_model,pickled_coarse_parameters,size_meshes):
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
                model_refined = refinement.compute_refinement_hessian_metric(simulation,size_meshes[lev+1],size_meshes[lev])
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
        serialized_model      : model serialized
        serialized_parameters : project parameters serialized
'''
@ExaquteTask(parameter_file_name=FILE_IN,returns=2)
def SerializeModelParameters_Task(parameter_file_name):
    with open(parameter_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    local_parameters = parameters
    model = KratosMultiphysics.Model()
    fake_sample = 1.0
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
    print("\n############## Serialization completed ##############\n")
    return pickled_model, pickled_parameters


'''
function executing the refinement of the problem
input:
        model       : serialization of the model
        parameters  : serialization of the Project Parameters
output:
        QoI                   : Quantity of Interest
        serialized_model      : model serialized
        serialized_parameters : parameters serialized
'''
@ExaquteTask(returns=2)
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
@ExaquteTask(returns=1)
def CompareMean_Task(AveragedMeanQoI,ExactExpectedValueQoI):
    relative_error = abs((AveragedMeanQoI - ExactExpectedValueQoI)/ExactExpectedValueQoI)
    return relative_error


if __name__ == '__main__':

    '''set the ProjectParameters.json path'''
    parameter_file_name = "/home/kratos105b/Kratos/applications/MultilevelMonteCarloApplication/tests/MeshCoarse8Nodes/ProjectParameters.json"
    '''create a serialization of the model and of the project parameters'''
    pickled_model,pickled_parameters = SerializeModelParameters_Task(parameter_file_name)
    '''customize setting parameters of the ML simulation'''
    settings_ML_simulation = KratosMultiphysics.Parameters("""
    {
        "tol0"                            : 0.25,
        "tolF"                            : 0.1,
        "cphi"                            : 1.0,
        "number_samples_screening"        : 15,
        "Lscreening"                      : 2,
        "Lmax"                            : 4,
        "initial_mesh_size"               : 0.5
    }
    """)
    '''contruct MultilevelMonteCarlo class'''
    mlmc_class = mlmc.MultilevelMonteCarlo(settings_ML_simulation)
    ''''start screening phase'''
    for lev in range(mlmc_class.current_number_levels+1):
        for instance in range (mlmc_class.number_samples[lev]):
            mlmc_class.AddResults(ExecuteMultilevelMonteCarloAnalisys_Task(lev,pickled_model,pickled_parameters,mlmc_class.sizes_mesh),lev)
    '''finalize screening phase'''
    mlmc_class.FinalizeScreeningPhase()
    mlmc_class.ScreeningInfoScreeningPhase()
    '''start MLMC phase'''
    while mlmc_class.convergence is not True:
        '''initialize MLMC phase'''
        mlmc_class.InitializeMLMCPhase()
        mlmc_class.ScreeningInfoInitializeMLMCPhase()
        '''MLMC execution phase'''
        for lev in range (mlmc_class.current_number_levels+1):
            for instance in range (mlmc_class.difference_number_samples[lev]):
                mlmc_class.AddResults(ExecuteMultilevelMonteCarloAnalisys_Task(lev,pickled_model,pickled_parameters,mlmc_class.sizes_mesh),lev)
        '''finalize MLMC phase'''
        mlmc_class.FinalizeMLMCPhase()
        mlmc_class.ScreeningInfoFinalizeMLMCPhase()

    print("\niterations = ",mlmc_class.current_iteration,\
    "total error TErr computed = ",mlmc_class.TErr,"mean MLMC QoI = ",mlmc_class.mean_mlmc_QoI)

    '''### OBSERVATIONS ###

    compss: between different tasks you don't need compss_wait_on/get_value_from_remote, it's pycompss who handles everything automatically
    if the output of a task is given directly to the input of an other, pycompss handles everything

    MmgProcess: need to use conditions in the model part to preserve the boundary conditions in the refinement process
    The submodelpart: Subpart_Boundary contains only nodes and no geometries (conditions/elements).
    It is not guaranteed that the submodelpart will be preserved.
    PLEASE: Add some "dummy" conditions to the submodelpart to preserve it'''