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

class MonteCarloAnalysis(AnalysisStage):
    """Main script for Monte Carlo simulations using the pure_diffusion solver"""

    
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


    def ApplyBoundaryConditions(self):
        super(MonteCarloAnalysis,self).ApplyBoundaryConditions()
        '''define the forcing function'''
        for node in self.model.GetModelPart("MCLaplacianModelPart").Nodes:
            coord_x = node.X
            coord_y = node.Y
            # forcing = -432.0 * coord_x * (coord_x - 1) * coord_y * (coord_y - 1)
            forcing = -432.0 * (coord_x**2 + coord_y**2 - coord_x - coord_y)
            node.SetSolutionStepValue(KratosMultiphysics.HEAT_FLUX,forcing*self.sample)
            
        
        '''need to use SetSolutionStepValue and not SetValue (equivalent for GetSolutionStepValue and GetValue) because in custom_elements and custom_conditions I use this function
        Set/GetValue uses less memory and is always used for element and conditions
        Set/GetSolutionStepValue uses more memory, and allows to memorize values stored in the nodes in case of buffer > 1 (e.g. sol_{i-1}, sol_{i-2} remain stored there to evaluate sol_{i})
        if I use SetValue, I need to use GetValue (the same for the other)
        # node.SetValue(Poisson.FORCING,node.GetSolutionStepValue(Poisson.FORCING)*sample)'''


    
##################################################
######## END OF CLASS MONTECARLOANALYSIS #########
##################################################



# def GenerateStochasticContribute(parameters):
#     number_samples = 1 # i.e. generate one stochastic variable per time

#     if (parameters["problem_data"].Has("stochastic_pdf")):
#         stochastic_pdf = parameters["problem_data"]["stochastic_pdf"]
#         print(stochastic_pdf)
#     else:
#         raise Exception('Please provide the "stochastic_pdf" parameter in the .json file')

#     if stochastic_pdf.Has("normal_distribution"):
#         if stochastic_pdf["normal_distribution"].Has("mean"):
#             mu = stochastic_pdf["normal_distribution"]["mean"].GetDouble()
#         else:
#             raise Exception('Please define the "mean" for the normal distribution in the .json file')            
#         if stochastic_pdf["normal_distribution"].Has("variance"):
#             sigma = stochastic_pdf["normal_distribution"]["variance"].GetDouble()
#         else:
#             raise Exception('Please define the "variance" for the normal distribution in the .json file')
#         sample = np.random.normal(mu,sigma,number_samples)

#     elif stochastic_pdf.Has("beta_distribution"):
#         if stochastic_pdf["beta_distribution"].Has("alpha"):
#             alpha = stochastic_pdf["beta_distribution"]["alpha"].GetDouble()
#         else:
#             raise Exception('Please define the "alpha" for the beta distribution in the .json file')            
#         if stochastic_pdf["beta_distribution"].Has("beta"):
#             beta = stochastic_pdf["beta_distribution"]["beta"].GetDouble()
#         else:
#             raise Exception('Please define the "beta" for the beta distribution in the .json file')
#         sample = np.random.beta(alpha,beta,number_samples)

#     else:
#         raise Exception('Please provide "normal_distribution" or "beta_distribution" in the .json file; at the moment only "normal_distribution" and "beta_distribution" are implemented')
#     return sample


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
The following commented execution_task function is analogous to the uncommented execution_task,
the difference is what the function returns, and we need only the Quantity of Interest (QoI)
'''
#@task(model_part_file_name=FILE_IN, parameter_file_name=FILE_IN, returns=4)
# def execution_task(model_part_file_name, parameter_file_name):
#     with open(parameter_file_name,'r') as parameter_file:
#         parameters = KratosMultiphysics.Parameters(parameter_file.read())
#     local_parameters = parameters # in case there are more parameters file, we rename them
#     model = KratosMultiphysics.Model()      

#     local_parameters["solver_settings"]["model_import_settings"]["input_filename"].SetString(model_part_file_name[:-5])
#     alpha = 2.0
#     beta = 6.0
#     sample = GenerateBetaSample(alpha,beta)
#     simulation = MonteCarloAnalysis(model,local_parameters, sample)
#     simulation.Run() 
#     QoI =  EvaluateQuantityOfInterest(simulation)
#     KratosMultiphysics.CalculateNodalAreaProcess(simulation._GetSolver().main_model_part,2).Execute()
#     node_information = []
#     for node in simulation._GetSolver().main_model_part.Nodes:
#         node_information.append([node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE), node.GetSolutionStepValue(KratosMultiphysics.NODAL_AREA), node.X, node.Y])
#     return 1, QoI, node_information, simulation.sample


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
The following is the function implemented with Ramon to compute the mean
Right now we are using the update_onepass_M function, computed in mc.py
'''
# @task(returns=2)
# def mean_task(*args):
#     sum_val = 0
#     amount_elems = 0
#     for i in range(int(len(args) / 2)):
#         sum_val = sum_val + (args[2 * i] * args[2 * i + 1])
#         amount_elems = amount_elems + args[2 * i]
#     return amount_elems, sum_val / amount_elems


'''
function executing the problem for sample = 1.0
input:
        model_part_file_name  : path of the model part file
        parameter_file_name   : path of the Project Parameters file
output:
        ExactExpectedValueQoI : Quantity of Interest for sample = 1.0
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
        parameter_file_name = "/home/kratos105b/Kratos/applications/MonteCarloApplication/tests/Level1/ProjectParameters.json"

    with open(parameter_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    local_parameters = parameters # in case there are more parameters file, we rename them
    
    number_samples = 10
    Qlist = []
    run_results = []

    '''evaluate the exact expected value of Q (sample = 1.0)'''
    ExactExpectedValueQoI = exact_execution_task(local_parameters["solver_settings"]["model_import_settings"]["input_filename"].GetString() + ".mdpa", parameter_file_name)
    
    for instance in range (0,number_samples):
        Qlist.append(execution_task(local_parameters["solver_settings"]["model_import_settings"]["input_filename"].GetString() + ".mdpa", parameter_file_name))
        
    
    '''The following lines compute the mean using the mean_task function'''
    # for elem in Qlist:
        #curr_elem = compss_wait_on(elem[0])
        # Qlist_aux.append(1)
        # Qlist_aux.append(elem)

    # chunk_size = 4
    '''use this second while loop if you want to append at the end of Qlist the mean value you evaluate'''
    # while len(Qlist_aux) > 2:
    #     newQlist = Qlist_aux[2*chunk_size:]
    #     card, val = mean_task(*Qlist_aux[0:2*chunk_size])
    #     newQlist.append(card)
    #     newQlist.append(val)
    #     Qlist_aux = newQlist
    '''use this second while loop if you want to append at the beginning of Qlist the mean value you evaluate'''
    # while len(Qlist_aux) > 2:
    #     card, val = mean_task(*Qlist_aux[0:2*chunk_size])
    #     newQlist = []
    #     newQlist.append(card)
    #     newQlist.append(val)
    #     newQlist.extend(Qlist_aux[2*chunk_size:])
    #     Qlist_aux = newQlist

    # print("Size of Qlist_aux: " + str(len(Qlist_aux)) + " " + str(Qlist_aux))
    # mean = Qlist_aux[1]

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


    """ The below part evaluates the relative L2 error between the numerical solution SOLUTION(x,y,sample) and the analytical solution, also dependent on sample.
    Analytical solution available in case FORCING = sample * -432.0 * (coord_x**2 + coord_y**2 - coord_x - coord_y)"""
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
    # print("")
    
