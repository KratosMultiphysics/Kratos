from __future__ import absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import numpy as np
import time

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.ConvectionDiffusionApplication as KratosConvDiff
import KratosMultiphysics.MonteCarloApplication as KratosMC

# Avoid printing of Kratos informations
KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING) # avoid printing of Kratos things

# Importing the base class
from analysis_stage import AnalysisStage

# Import pycompss
from pycompss.api.task import task
from pycompss.api.api import compss_wait_on
from pycompss.api.parameter import *

# Import Monte Carlo library
import mc as mc


'''Adapt the following class depending on the problem, deriving the MonteCarloAnalysis class from the problem of interest'''

'''This Analysis Stage implementation solves the elliptic PDE in (0,1)^2 with zero Dirichlet boundary conditions
-lapl(u) = xi*f,    f= -432*x*(x-1)*y*(y-1)
                    f= -432*(x**2+y**2-x-y)
where xi is a Beta(2,6) random variable, and computes statistic of the QoI
Q = int_(0,1)^2 u(x,y)dxdy
more details in Section 5.2 of [PKN17]

References:
[PKN17] M. Pisaroni; S. Krumscheid; F. Nobile : Quantifying uncertain system outputs via the multilevel Monte Carlo method - Part I: Central moment estimation; MATHICSE technical report no. 23.2017.
'''
class MonteCarloAnalysis(AnalysisStage):
    '''Main script for Monte Carlo simulations'''

    def __init__(self,model,parameters,sample):
        self.sample = sample
        super(MonteCarloAnalysis,self).__init__(model,parameters)
        self._GetSolver().main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)

    def _CreateSolver(self):
        import convection_diffusion_stationary_solver
        solver = convection_diffusion_stationary_solver.CreateSolver(self.model,self.project_parameters["solver_settings"])
        self.LaplacianSolver = solver
        return self.LaplacianSolver
    
    def _GetSimulationName(self):
        return "Monte Carlo Analysis"

    '''Introduce here the stochasticity in the right hand side defining the forcing function and apply the stochastic contribute'''
    def ModifyInitialProperties(self):
        for node in self.model.GetModelPart("MLMCLaplacianModelPart").Nodes:
            coord_x = node.X
            coord_y = node.Y
            # forcing = -432.0 * coord_x * (coord_x - 1) * coord_y * (coord_y - 1)
            forcing = -432.0 * (coord_x**2 + coord_y**2 - coord_x - coord_y) # this forcing presents an analytical solution
            node.SetSolutionStepValue(KratosMultiphysics.HEAT_FLUX,forcing*self.sample)


##################################################
######## END OF CLASS MONTECARLOANALYSIS #########
##################################################


'''
function generating the random sample
here the sample has a beta distribution with parameters alpha = 2.0 and beta = 6.0
'''
def GenerateBetaSample():
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
        #print("NODAL AREA = ",node.GetSolutionStepValue(KratosMultiphysics.NODAL_AREA),"NODAL SOLUTION = ",node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE),"CURRENT Q = ",Q)
    return Q


'''
function executing the problem
input:
        model_part_file_name : path of the model part file (still to implement how to import in efficient way in a loop where I have different model part files and different ProjectParameters files, thus for now read model part name from the ProjectParameters.json file)
        parameter_file_name  : path of the Project Parameters file
        sample               : stochastic random variable
output:
        QoI                  : Quantity of Interest
'''
@task(model_part_file_name=FILE_IN, parameter_file_name=FILE_IN, returns=1)
def execution_task(model_part_file_name, parameter_file_name):
    with open(parameter_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    local_parameters = parameters # in case there are more parameters file, we rename them
    model = KratosMultiphysics.Model()
    local_parameters["solver_settings"]["model_import_settings"]["input_filename"].SetString(model_part_file_name[:-5])
    sample = GenerateBetaSample()
    simulation = MonteCarloAnalysis(model,local_parameters,sample)
    simulation.Run()
    QoI =  EvaluateQuantityOfInterest(simulation)
    return QoI


'''
function executing the problem for sample = 1.0
input:
        model_part_file_name  : path of the model part file
        parameter_file_name   : path of the Project Parameters file
output:
        ExactExpectedValueQoI : Quantity of Interest for sample = 1.0
OBSERVATION: here we multiply by 0.25 because it is the mean value of beta(2,6)
'''
@task(model_part_file_name=FILE_IN, parameter_file_name=FILE_IN,returns=1)
def exact_execution_task(model_part_file_name, parameter_file_name):
    with open(parameter_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    local_parameters = parameters # in case there are more parameters file, we rename them
    model = KratosMultiphysics.Model()      
    local_parameters["solver_settings"]["model_import_settings"]["input_filename"].SetString(model_part_file_name[:-5])
    sample = 1.0
    simulation = MonteCarloAnalysis(model,local_parameters, sample)
    simulation.Run() 
    QoI =  EvaluateQuantityOfInterest(simulation)
    ExactExpectedValueQoI = 0.25 * EvaluateQuantityOfInterest(simulation)
    return ExactExpectedValueQoI



'''
function computing the relative error between the Multilevel Monte Carlo expected value and the exact expected value
input :
        AveragedMeanQoI       : Multilevel Monte Carlo expected value
        ExactExpectedValueQoI : exact expected value
output :
        relative_error        : relative error
'''
@task(returns=1)
def compare_mean(AveragedMeanQoI,ExactExpectedValueQoI):
    relative_error = abs((AveragedMeanQoI - ExactExpectedValueQoI)/ExactExpectedValueQoI)
    return relative_error


if __name__ == '__main__':
    from sys import argv

    if len(argv) > 2:
        err_msg = 'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default parameter file (assumed to be called "ProjectParameters.json"):\n'
        err_msg += '    "python montecarlo_analysis.py"\n'
        err_msg += '- With custom parameter file:\n'
        err_msg += '    "python3 montecarlo_analysis.py <my-parameter-file>.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        parameter_file_name = argv[1]
    else: # using default name
        parameter_file_name = "../tests/Level0/ProjectParameters.json"

    with open(parameter_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    local_parameters = parameters # in case there are more parameters file, we rename them
    
    number_samples = 10
    Qlist = []

    '''evaluate the exact expected value of Q (sample = 1.0)'''
    ExactExpectedValueQoI = exact_execution_task(local_parameters["solver_settings"]["model_import_settings"]["input_filename"].GetString() + ".mdpa", parameter_file_name)
    
    for instance in range (0,number_samples):
        Qlist.append(execution_task(local_parameters["solver_settings"]["model_import_settings"]["input_filename"].GetString() + ".mdpa", parameter_file_name))
        

    '''Compute mean, second moment and sample variance'''
    MC_mean = 0.0
    MC_second_moment = 0.0
    for i in range (0,number_samples):
        nsam = i+1
        MC_mean, MC_second_moment, MC_variance = mc.update_onepass_M(Qlist[i], MC_mean, MC_second_moment, nsam)
    '''Evaluation of the relative error between the computed mean value and the expected value of the QoI'''
    relative_error = compare_mean(MC_mean,ExactExpectedValueQoI)
    # print("Values QoI:",Qlist)
    MC_mean = compss_wait_on(MC_mean)
    ExactExpectedValueQoI = compss_wait_on(ExactExpectedValueQoI)
    relative_error = compss_wait_on(relative_error)
    print("\nMC mean = ",MC_mean,"exact mean = ",ExactExpectedValueQoI)
    print("relative error: ",relative_error)


    ''' The below part evaluates the relative L2 error between the numerical solution SOLUTION(x,y,sample) and the analytical solution, also dependent on sample.
    Analytical solution available in case FORCING = sample * -432.0 * (coord_x**2 + coord_y**2 - coord_x - coord_y)'''
    # model = KratosMultiphysics.Model()
    # sample = 1.0
    # simulation = MonteCarloAnalysis(model,local_parameters,sample)
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
   
