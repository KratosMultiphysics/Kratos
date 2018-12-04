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

# Import pycompss
# from pycompss.api.task import task
# from pycompss.api.api import compss_wait_on
# from pycompss.api.parameter import *

# Import exaqute
# from exaqute.ExaquteTaskPyCOMPSs import *   # to exequte with pycompss
# from exaqute.ExaquteTaskHyperLoom import *  # to exequte with the IT4 scheduler
# from exaqute.ExaquteTaskLocal import *      # to execute with python3
# get_value_from_remote is the equivalent of compss_wait_on
# in the future, when everything is integrated with the it4i team, putting exaqute.ExaquteTaskHyperLoom you can launch your code with their scheduler instead of BSC

# Import Continuation Multilevel Monte Carlo library
import cmlmc_utilities_CLASS as mlmc

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
        self._GetSolver().main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_H)

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
def GenerateBetaSample(alpha,beta):
    number_samples = 1
    sample = np.random.beta(alpha,beta,number_samples)
    return sample


'''
function evaluating the QoI of the problem: int_{domain} TEMPERATURE(x,y) dx dy
right now we are using the midpoint rule to evaluate the integral: improve!
'''
def EvaluateQuantityOfInterest(simulation):
    KratosMultiphysics.CalculateNodalAreaProcess(simulation._GetSolver().main_model_part,2).Execute()
    Q = 0.0
    for node in simulation._GetSolver().main_model_part.Nodes:
        Q = Q + (node.GetSolutionStepValue(KratosMultiphysics.NODAL_AREA)*node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE))
        #print("NODAL AREA = ",node.GetSolutionStepValue(KratosMultiphysics.NODAL_AREA),"NODAL SOLUTION = ",node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE),"CURRENT Q = ",Q)
    return Q


# '''
# function executing the problem
# input:
#         model       : serialization of the model
#         parameters  : serialization of the Project Parameters
#         sample      : stochastic random variable
# output:
#         QoI         : Quantity of Interest
# '''
# # @task(parameter_file_name=FILE_IN, returns=1)
# def execution_task(parameter_file_name, sample):
#     with open(parameter_file_name,'r') as parameter_file:
#         parameters = KratosMultiphysics.Parameters(parameter_file.read())
#     local_parameters = parameters # in case there are more parameters file, we rename them
#     model = KratosMultiphysics.Model()
#     # local_parameters["solver_settings"]["model_import_settings"]["input_filename"].SetString(model_part_file_name[:-5])
#     simulation = MultilevelMonteCarloAnalysis(model,local_parameters,sample)
#     simulation.Run()
#     QoI = EvaluateQuantityOfInterest(simulation)
#     return QoI
    

# '''
# function executing the problem for sample = 1.0
# input:
#         model       : serialization of the model
#         parameters  : serialization of the Project Parameters
# output:
#         ExactExpectedValueQoI : Quantity of Interest for sample = 1.0
# OBSERVATION: here we multiply by 0.25 because it is the mean value of beta(2,6)
# '''
# # @ExaquteTask(returns=1)
# def exact_execution_task(model, parameters):
#     sample = 1.0
#     simulation = MultilevelMonteCarloAnalysis(model,parameters,sample)
#     simulation.Run()
#     QoI =  EvaluateQuantityOfInterest(simulation)
#     ExactExpectedValueQoI = 0.25 * EvaluateQuantityOfInterest(simulation)
#     return ExactExpectedValueQoI


# '''
# function serializing the model and the parameters of the problem
# input:
#         parameter_file_name   : path of the Project Parameters file
# output:
#         serialized_model      : model serialized
#         serialized_parameters : project parameters serialized
# '''
# # @ExaquteTask(parameter_file_name=FILE_IN,returns=1)
# def serialize_model_projectparameters(parameter_file_name):
#     with open(parameter_file_name,'r') as parameter_file:
#         parameters = KratosMultiphysics.Parameters(parameter_file.read())
#     local_parameters = parameters
#     model = KratosMultiphysics.Model()      
#     # local_parameters["solver_settings"]["model_import_settings"]["input_filename"].SetString(model_part_file_name[:-5])
#     fake_sample = 1.0

#     simulation = MultilevelMonteCarloAnalysis(model,local_parameters,fake_sample)
#     simulation.Initialize()

#     serialized_model = KratosMultiphysics.StreamSerializer()
#     serialized_model.Save("ModelSerialization",simulation.model)
#     serialized_parameters = KratosMultiphysics.StreamSerializer()
#     serialized_parameters.Save("ParametersSerialization",simulation.project_parameters)
#     objects_serialized = [serialized_model, serialized_parameters]
#     return objects_serialized


# '''
# function executing the problem
# input:
#         model       : serialization of the model
#         parameters  : serialization of the Project Parameters
# output:
#         QoI         : Quantity of Interest
# '''
# # @ExaquteTask(returns=1)
# def execution_task(pickled_model, pickled_parameters):
#     '''overwrite the old model serializer with the unpickled one'''
#     model_serializer = pickle.loads(pickled_model)
#     current_model = KratosMultiphysics.Model()
#     model_serializer.Load("ModelSerialization",current_model)
#     del(model_serializer)

#     '''overwrite the old parameters serializer with the unpickled one'''
#     serialized_parameters = pickle.loads(pickled_parameters)
#     current_parameters = KratosMultiphysics.Parameters()
#     serialized_parameters.Load("ParametersSerialization",current_parameters)
#     del(serialized_parameters)

#     sample = GenerateBetaSample()
#     simulation = MonteCarloAnalysis(current_model,current_parameters,sample)
#     simulation.Run()
#     QoI = EvaluateQuantityOfInterest(simulation)
#     return QoI


# '''
# function executing the problem for sample = 1.0
# input:
#         model       : serialization of the model
#         parameters  : serialization of the Project Parameters
# output:
#         ExactExpectedValueQoI : Quantity of Interest for sample = 1.0
# OBSERVATION: here we multiply by 0.25 because it is the mean value of beta(2,6)
# '''
# # @ExaquteTask(returns=1)
# def exact_execution_task(pickled_model, pickled_parameters):
#     '''overwrite the old model serializer with the unpickled one'''
#     model_serializer = pickle.loads(pickled_model)
#     current_model = KratosMultiphysics.Model()
#     model_serializer.Load("ModelSerialization",current_model)
#     del(model_serializer)

#     '''overwrite the old parameters serializer with the unpickled one'''
#     serialized_parameters = pickle.loads(pickled_parameters)
#     current_parameters = KratosMultiphysics.Parameters()
#     serialized_parameters.Load("ParametersSerialization",current_parameters)
#     del(serialized_parameters)

#     sample = 1.0
#     simulation = MonteCarloAnalysis(current_model,current_parameters,sample)
#     simulation.Run()
#     ExactExpectedValueQoI = 0.25 * EvaluateQuantityOfInterest(simulation)
#     return ExactExpectedValueQoI


# '''
# function executing the problem
# input:
#         model       : serialization of the model
#         parameters  : serialization of the Project Parameters
# output:
#         QoI                   : Quantity of Interest
#         serialized_model      : model serialized
#         serialized_parameters : parameters serialized
# '''
# # @ExaquteTask(returns=2)
# def execution_task_refinement(pickled_model_coarse, pickled_parameters, min_size, max_size):
#     sample = GenerateBetaSample()

#     '''overwrite the old model serializer with the unpickled one'''
#     model_serializer_coarse = pickle.loads(pickled_model_coarse)
#     model_coarse = KratosMultiphysics.Model()
#     model_serializer_coarse.Load("ModelSerialization",model_coarse)
#     del(model_serializer_coarse)

#     '''overwrite the old parameters serializer with the unpickled one'''
#     serialized_parameters = pickle.loads(pickled_parameters)
#     parameters_refinement = KratosMultiphysics.Parameters()
#     serialized_parameters.Load("ParametersSerialization",parameters_refinement)
#     del(serialized_parameters)

#     simulation_coarse = MonteCarloAnalysis(model_coarse,parameters_refinement,sample)
#     simulation_coarse.Run()
#     QoI =  EvaluateQuantityOfInterest(simulation_coarse)

#     '''refine here'''
#     model_refined = refinement.compute_refinement_from_analysisstage_object(simulation_coarse,min_size,max_size)
    
#     '''initialize'''
#     simulation = MonteCarloAnalysis(model_refined,parameters_refinement,sample)
#     simulation.Initialize()

#     '''serialize model and pickle it'''
#     serialized_model = KratosMultiphysics.StreamSerializer()
#     serialized_model.Save("ModelSerialization",simulation.model)
#     pickled_model_refined = pickle.dumps(serialized_model, 2)

#     return QoI,pickled_model_refined


##with refinement
def ExecuteMultilevelMonteCarloAnalisys(finest_level,pickled_coarse_model,pickled_coarse_parameters):
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
    sample = GenerateBetaSample(2.0,6.0)
    QoI = []
    start_MLMC_time = time.time()
    if(finest_level == 0):
        QoI.append(0.0)
        simulation = MultilevelMonteCarloAnalysis(current_model,current_parameters,sample)
        simulation.Run()
        QoI.append(EvaluateQuantityOfInterest(simulation))
        del(simulation)
    else:
        for lev in range(finest_level+1):
            simulation = MultilevelMonteCarloAnalysis(current_model,current_parameters,sample)
            simulation.Run()
            QoI.append(EvaluateQuantityOfInterest(simulation))
            '''refine if level < finest level'''
            if (lev < finest_level):
                print("CURRENT LEVEL, MUST ARRIVE TO FINEST LEVEL - 1",lev,finest_level)
                '''refine the model Kratos object'''
                model_refined = refinement.compute_refinement_from_analysisstage_object(simulation,0.3,0.4)
                '''initialize the model Kratos object'''
                simulation = MultilevelMonteCarloAnalysis(model_refined,current_parameters,sample)
                simulation.Initialize()
                '''update model Kratos object'''
                current_model = simulation.model
            del(simulation)
    end_MLMC_time = time.time()
    '''prepare results of the simulation'''
    print(QoI)
    results_simulation = {
        "QoI_finer_level":QoI[finest_level+1],\
        "QoI_coarser_level":QoI[finest_level],\
        "finer_level":finest_level,\
        "coarser_level":finest_level-1,\
        "total_MLMC_time":end_MLMC_time - start_MLMC_time}
    return(results_simulation)




'''
function serializing the model and the parameters of the problem
the idea is the following:
first from Model/Parameters Kratos object to StreamSerializer Kratos object
second from StreamSerializer Kratos object to pickle string
third from pickle string to StreamSerializer Kratos object
fourth from StreamSerializer Kratos object to Model/Parameters Kratos object
input:
        parameter_file_name   : path of the Project Parameters file
output:
        serialized_model      : model serialized
        serialized_parameters : project parameters serialized
'''
# @ExaquteTask(parameter_file_name=FILE_IN,returns=2)
def serialize_model_projectparameters(parameter_file_name):
    with open(parameter_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    local_parameters = parameters
    model = KratosMultiphysics.Model()
    # local_parameters["solver_settings"]["model_import_settings"]["input_filename"].SetString(model_part_file_name[:-5])
    fake_sample = 1.0

    simulation = MultilevelMonteCarloAnalysis(model,local_parameters,fake_sample)
    simulation.Initialize()

    serialized_model = KratosMultiphysics.StreamSerializer()
    serialized_model.Save("ModelSerialization",simulation.model)

    serialized_parameters = KratosMultiphysics.StreamSerializer()
    serialized_parameters.Save("ParametersSerialization",simulation.project_parameters)

    # pickle dataserialized_data
    pickled_model = pickle.dumps(serialized_model, 2)
    pickled_parameters = pickle.dumps(serialized_parameters, 2)

    return pickled_model, pickled_parameters



'''
function computing the relative error between the Multilevel Monte Carlo expected value and the exact expected value
input :
        AveragedMeanQoI       : Multilevel Monte Carlo expected value
        ExactExpectedValueQoI : exact expected value
output :
        relative_error        : relative error
'''
# @ExaquteTask(returns=1)
def compare_mean(AveragedMeanQoI,ExactExpectedValueQoI):
    relative_error = abs((AveragedMeanQoI - ExactExpectedValueQoI)/ExactExpectedValueQoI)
    return relative_error


if __name__ == '__main__':

    '''set the ProjectParameters.json path in the parameter_file_name list'''
    # parameter_file_name = "/home/kratos105b/Kratos/applications/MultilevelMonteCarloApplication/tests/Level0/ProjectParameters.json"

    # parameter_file_name =[]
    # parameter_file_name.append("../tests/Level0/ProjectParameters.json")
    # with open(parameter_file_name[0],'r') as parameter_file:
    #     parameters = KratosMultiphysics.Parameters(parameter_file.read())
    # local_parameters_0 = parameters # in case there are more parameters file, we rename them
    # parameter_file_name.append("../tests/Level1/ProjectParameters.json")
    # with open(parameter_file_name[1],'r') as parameter_file:
    #     parameters = KratosMultiphysics.Parameters(parameter_file.read())
    # local_parameters_1 = parameters # in case there are more parameters file, we rename them
    # parameter_file_name.append("../tests/Level2/ProjectParameters.json")
    # with open(parameter_file_name[2],'r') as parameter_file:
    #     parameters = KratosMultiphysics.Parameters(parameter_file.read())
    # local_parameters_2 = parameters # in case there are more parameters file, we rename them
    # parameter_file_name.append("../tests/Level3/ProjectParameters.json")
    # with open(parameter_file_name[3],'r') as parameter_file:
    #     parameters = KratosMultiphysics.Parameters(parameter_file.read())
    # local_parameters_3 = parameters # in case there are more parameters file, we rename them
    # parameter_file_name.append("../tests/Level4/ProjectParameters.json")
    # with open(parameter_file_name[4],'r') as parameter_file:
    #     parameters = KratosMultiphysics.Parameters(parameter_file.read())
    # local_parameters_4 = parameters # in case there are more parameters file, we rename them
    parameter_file_name = "/home/riccardo/Kratos/applications/MultilevelMonteCarloApplication/tests/MeshCoarse8Nodes/ProjectParameters.json"
    '''create a serialization of the model and of the project parameters'''
    pickled_model,pickled_parameters = serialize_model_projectparameters(parameter_file_name)
    print("\n############## Serialization completed ##############\n")

    '''define setting parameters of the ML simulation'''
    k0   = 0.1 # Certainty Parameter 0 rates
    k1   = 0.1 # Certainty Parameter 1 rates
    r1   = 1.25 # Cost increase first iterations C-MLMC
    r2   = 1.15 # Cost increase final iterations C-MLMC
    tol0 = 0.25 # Tolerance iter 0
    tolF = 0.1 # Tolerance final
    cphi = 1.0 # Confidence on tolerance
    N0   = 25 # Number of samples for iter 0
    L0   = 2 # Number of levels for iter 0
    Lmax = 4
    M = 4 # mesh refinement coefficient
    initial_mesh_size = 0.25
    settings_ML_simulation = [k0,k1,r1,r2,tol0,tolF,cphi,N0,L0,Lmax,M,initial_mesh_size]

    '''contruct MultilevelMonteCarlo class'''
    mlmc_class = mlmc.MultilevelMonteCarlo(settings_ML_simulation)

    print(mlmc_class.difference_QoI.values,mlmc_class.time_ML.values)
    
    results = ExecuteMultilevelMonteCarloAnalisys(3,pickled_model,pickled_parameters)
    print(results["QoI_finer_level"])

    stop

    print("\n ######## SCREENING PHASE ######## \n")

    for level in range(mlmc_class.current_number_levels+1):
        for instance in range (mlmc_class.number_samples[level]):
            sample = GenerateBetaSample(2.0,6.0) # generate a random variable with beta pdf, alpha = 2.0 and beta = 6.0
            run_results = []
            start_time_ML = time.time() # I can insert this operation in "GenerateBetaSample", or better to create a new function?

            if level == 0: # evaluating QoI in the coarsest grid
                run_results.append(execution_task(parameter_file_name[level], sample)) # append to run_results QoI for the coarsest grid
                time_MLi = time.time() - start_time_ML # create a new function?
                # difference_QoI[level].append(run_results[-1])
                mlmc_class.difference_QoI.values[level] = np.append(mlmc_class.difference_QoI.values[level],run_results[-1]) # with list[-1] we read the last element of the list
                # time_ML[level].append(time_MLi)
                mlmc_class.time_ML.values[level] = np.append(mlmc_class.time_ML.values[level],time_MLi)
                
            else:
                for cycle_level in range (level-1,level+1):
                    run_results.append(execution_task(parameter_file_name[cycle_level], sample))
                  
                time_MLi = time.time() - start_time_ML
                # difference_QoI[level].append(run_results[-1] - run_results[-2])
                mlmc_class.difference_QoI.values[level] = np.append(mlmc_class.difference_QoI.values[level],run_results[-1] - run_results[-2])
                # time_ML[level].append(time_MLi)
                mlmc_class.time_ML.values[level] = np.append(mlmc_class.time_ML.values[level],time_MLi)
    
    # print("values",mlmc_class.difference_QoI.values,mlmc_class.time_ML.values)
    mlmc_class.FinalizeScreeningPhase()
    print("mean,variance QoI",mlmc_class.difference_QoI.mean,mlmc_class.difference_QoI.sample_variance)
    print("mean,variance time",mlmc_class.time_ML.mean,mlmc_class.time_ML.sample_variance)
    print(mlmc_class.mesh_parameters)
    print(mlmc_class.rates_error)
    print(mlmc_class.BayesianVariance)
    print(mlmc_class.number_iterations_iE)

    while mlmc_class.convergence is not True:
        print("\n ######## CMLMC iter = ",mlmc_class.current_iteration,"######## \n")
        
        mlmc_class.InitializeMLMCPhase()
        print(mlmc_class.tolerance_i)

        print (mlmc_class.BayesianVariance)
        print (mlmc_class.current_number_levels)
        print (mlmc_class.previous_number_levels)

        print(mlmc_class.theta_i)

        print(mlmc_class.number_samples,mlmc_class.difference_number_samples,mlmc_class.previous_number_samples)

        '''run hierarchy using optimal levels and optimal number of samples
        I use the "old" difference_QoI, and append the new values for the added samples,
        or append at the end for the new level'''

        for level in range (mlmc_class.current_number_levels+1):
            for instance in range (mlmc_class.difference_number_samples[level]):
                sample = GenerateBetaSample(2.0,6.0)
                run_results = []
                start_time_ML = time.time()
                if level == 0: # evaluating QoI in the coarsest grid
                    run_results.append(execution_task(parameter_file_name[level], sample)) # append to run_results QoI for the coarsest grid
                    time_MLi = time.time() - start_time_ML
                    # difference_QoI[level].append(run_results[-1]) # with list[-1] we read the last element of the list
                    # time_ML[level].append(time_MLi)
                    mlmc_class.difference_QoI.values[level] = np.append(mlmc_class.difference_QoI.values[level],run_results[-1])
                    mlmc_class.time_ML.values[level] = np.append(mlmc_class.time_ML.values[level],time_MLi)
                else:
                    for cycle_level in range (level-1,level+1):
                        run_results.append(execution_task(parameter_file_name[cycle_level], sample))

                    time_MLi = time.time() - start_time_ML
                    # difference_QoI[level].append(run_results[-1] - run_results[-2])
                    # time_ML[level].append(time_MLi)
                    mlmc_class.difference_QoI.values[level] = np.append(mlmc_class.difference_QoI.values[level],run_results[-1] - run_results[-2])
                    mlmc_class.time_ML.values[level] = np.append(mlmc_class.time_ML.values[level],time_MLi)

        # print("values",mlmc_class.difference_QoI.values,mlmc_class.time_ML.values)
        mlmc_class.FinalizeMLMCPhase()
        print("mean,variance QoI",mlmc_class.difference_QoI.mean,mlmc_class.difference_QoI.sample_variance)
        print("mean,variance time",mlmc_class.time_ML.mean,mlmc_class.time_ML.sample_variance)

        print(mlmc_class.mean_mlmc_QoI)
        print(mlmc_class.rates_error)
        print(mlmc_class.BayesianVariance)

        print(mlmc_class.TErr)
        

    print("\niterations = ",mlmc_class.current_iteration,\
    "total error TErr computed = ",mlmc_class.TErr,"mean MLMC QoI = ",mlmc_class.mean_mlmc_QoI)

    '''### OBSERVATION ###
    between different tasks you don't need compss_wait_on, it's pycompss who handles everything automatically
    if the output of a task is given directly to the input of an other, pycompss handles everything'''
