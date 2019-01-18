from __future__ import absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import numpy as np
import time

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.CompressiblePotentialFlowApplication as KratosCompressFlow
import KratosMultiphysics.MultilevelMonteCarloApplication as KratosMLMC

# Avoid printing of Kratos informations
KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING) # avoid printing of Kratos things

# Importing the base class
# from analysis_stage import AnalysisStage

# Importing derived classes
from potential_flow_analysis import PotentialFlowAnalysis

# Import exaqute
# from exaqute.ExaquteTaskPyCOMPSs import *   # to execute with pycompss
# from exaqute.ExaquteTaskHyperLoom import *  # to execute with the IT4 scheduler
from exaqute.ExaquteTaskLocal import *      # to execute with python3
'''
get_value_from_remote is the equivalent of compss_wait_on: a synchronization point
in future, when everything is integrated with the it4i team, importing exaqute.ExaquteTaskHyperLoom you can launch your code with their scheduler instead of BSC
'''

# Import Continuation Multilevel Monte Carlo library
import cmlmc_utilities as mlmc

# Import refinement library
import compressible_adaptive_refinement_utilities as refinement

# Import cpickle to pickle the serializer
try:
    import cpickle as pickle  # Use cPickle on Python 2.7
except ImportError:
    import pickle

class MultilevelMonteCarloAnalysis(PotentialFlowAnalysis):
    '''Main analysis stage for MultilevelMonte Carlo simulations'''
    def __init__(self,input_model,input_parameters,sample):
        self.sample = sample
        super(MultilevelMonteCarloAnalysis,self).__init__(input_model,input_parameters)
        # self._GetSolver().main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)

    def _GetSimulationName(self):
        return "Multilevel Monte Carlo Analysis"

    def ModifyInitialProperties(self):
        '''Introduce here the stochasticity in the Mach number and the angle of attack'''
        Mach = self.sample[0]
        a_infinity = 340 # [m/s] velocity of sound at infinity
        alpha =  self.sample[1]
        v_norm = Mach * a_infinity
        velocity = [v_norm*np.cos(alpha),v_norm*np.sin(alpha),0]
        boundary_processes = self.project_parameters["processes"]["boundary_conditions_process_list"]
        problem_name=self.project_parameters["solver_settings"]["model_import_settings"]["input_filename"].GetString()
        for i in range(0,boundary_processes.size()):
            python_module = boundary_processes[i]["python_module"].GetString()
            if python_module == "apply_far_field_process":
                self.project_parameters["processes"]["boundary_conditions_process_list"][i]["Parameters"]["velocity_infinity"].SetVector(velocity)
        # self.project_parameters["output_processes"]["gid_output"][0]["Parameters"]["output_name"].SetString(problem_name+'_M'+str(Mach)+'_A'+str(alpha))


############################################################
######## END OF CLASS MULTILEVELMONTECARLOANALYSIS #########
############################################################


'''
function generating the random sample
'''
def GenerateSample():
    sample = []
    mean_Mach = 0.3
    std_deviation_Mach = 0.01
    number_samples = 1
    sample.append(np.random.normal(mean_Mach,std_deviation_Mach,number_samples))
    mean_angle_attack = 0.0 # [rad] = 0 [degrees] airfoil already has 5 degrees
    std_deviation_angle_attack = np.deg2rad(0.1)
    sample.append(np.random.normal(mean_angle_attack,std_deviation_angle_attack,number_samples))
    print("Mach number = ",sample[0],"angle of attack = ",sample[1])
    if sample[0] >= 1.0 or sample[0] <= 0.0 :
        raise Exception ("stochastic Mach number computed > 1 or < 0")
    return sample


'''
function evaluating the QoI of the problem
'''
def EvaluateQuantityOfInterest(simulation):
    """here we evaluate the QoI of the problem: the lift coefficient"""
    Q = simulation._GetSolver().main_model_part.GetValue(KratosMultiphysics.FRICTION_COEFFICIENT)
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
def ExecuteMultilevelMonteCarloAnalisys(finest_level,pickled_coarse_model,pickled_coarse_parameters,size_meshes):
    return ExecuteMultilevelMonteCarloAnalisys_Task(finest_level,pickled_coarse_model,pickled_coarse_parameters,size_meshes),finest_level


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
                model_refined = refinement.compute_refinement_hessian_metric(simulation,size_meshes[lev+1],1.0)
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
@ExaquteTask(parameter_file_name=FILE_IN,returns=2)
def SerializeModelParameters_Task(parameter_file_name):
    with open(parameter_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    local_parameters = parameters
    model = KratosMultiphysics.Model()
    fake_sample = [0.3,0.0]
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
    print("\n","#"*50," SERIALIZATION COMPLETED ","#"*50,"\n")
    return pickled_model, pickled_parameters


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


if __name__ == '__main__':

    '''set the ProjectParameters.json path'''
    parameter_file_name = "../tests/CompressiblePotentialFlowTest/parameters_potential_naca_mesh0.json"
    '''create a serialization of the model and of the project parameters'''
    pickled_model,pickled_parameters = SerializeModelParameters_Task(parameter_file_name)
    '''customize setting parameters of the ML simulation'''
    settings_ML_simulation = KratosMultiphysics.Parameters("""
    {
        "tol0"                            : 0.25,
        "tolF"                            : 0.1,
        "cphi"                            : 1.0,
        "number_samples_screening"        : 2,
        "Lscreening"                      : 2,
        "Lmax"                            : 4,
        "initial_mesh_size"               : 0.05
    }
    """)
    '''contruct MultilevelMonteCarlo class'''
    mlmc_class = mlmc.MultilevelMonteCarlo(settings_ML_simulation)
    ''''start screening phase'''
    for lev in range(mlmc_class.current_number_levels+1):
        for instance in range (mlmc_class.number_samples[lev]):
            mlmc_class.AddResults(ExecuteMultilevelMonteCarloAnalisys(lev,pickled_model,pickled_parameters,mlmc_class.sizes_mesh))
    '''finalize screening phase'''
    mlmc_class.FinalizeScreeningPhase()
    mlmc_class.ScreeningInfoScreeningPhase()
    '''start MLMC phase'''
    # while mlmc_class.convergence is not True:
    #     '''initialize MLMC phase'''
    #     mlmc_class.InitializeMLMCPhase()
    #     mlmc_class.ScreeningInfoInitializeMLMCPhase()
    #     '''MLMC execution phase'''
    #     for lev in range (mlmc_class.current_number_levels+1):
    #         for instance in range (mlmc_class.difference_number_samples[lev]):
    #             mlmc_class.AddResults(ExecuteMultilevelMonteCarloAnalisys(lev,pickled_model,pickled_parameters,mlmc_class.sizes_mesh))
    #     '''finalize MLMC phase'''
    #     mlmc_class.FinalizeMLMCPhase()
    #     mlmc_class.ScreeningInfoFinalizeMLMCPhase()

    # print("\niterations = ",mlmc_class.current_iteration,\
    # "total error TErr computed = ",mlmc_class.TErr,"mean MLMC QoI = ",mlmc_class.mean_mlmc_QoI)

    '''### OBSERVATIONS ###

    compss: between different tasks you don't need compss_wait_on/get_value_from_remote, it's pycompss who handles everything automatically
    if the output of a task is given directly to the input of an other, pycompss handles everything

    compss: set absolute path when launching with compss

    MmgProcess: need to use conditions in the model part to preserve the boundary conditions in the refinement process
    The submodelpart: Subpart_Boundary contains only nodes and no geometries (conditions/elements).
    It is not guaranteed that the submodelpart will be preserved.
    PLEASE: Add some "dummy" conditions to the submodelpart to preserve it'''