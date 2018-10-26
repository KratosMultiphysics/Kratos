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

class MonteCarloAnalysis(AnalysisStage):
    """Main script for Monte Carlo simulations using the pure_diffusion solver"""

    
    def __init__(self,model,parameters,sample):
        self.sample = sample
        super(MonteCarloAnalysis,self).__init__(model,parameters)


    def _CreateSolver(self):
        import pure_diffusion_solver
        solver = pure_diffusion_solver.CreateSolver(self.model,self.project_parameters["solver_settings"])
        self.PureDiffusionSolver = solver
        return self.PureDiffusionSolver

    
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


def GenerateBetaSample(alpha,beta):
    number_samples = 1
    sample = np.random.beta(alpha,beta,number_samples)
    return sample


def EvaluateQuantityOfInterest(simulation):
    """here we evaluate the QoI of the problem: int_{domain} SOLUTION(x,y) dx dy
    we use the midpoint rule to evaluate the integral"""
    KratosMultiphysics.CalculateNodalAreaProcess(simulation._GetSolver().main_model_part,2).Execute()
    Q = 0.0
    for node in simulation._GetSolver().main_model_part.Nodes:
        Q = Q + (node.GetSolutionStepValue(KratosMultiphysics.NODAL_AREA)*node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE))
        #print("NODAL AREA = ",node.GetSolutionStepValue(KratosMultiphysics.NODAL_AREA),"NODAL SOLUTION = ",node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE),"CURRENT Q = ",Q)
    return Q


def update_onepass_M(sample, old_mean, old_M2, nsam):
    '''update mean and second moment values
    M_{2,n} = sum_{i=1}^{n} (x_i - mean(x)_n)^2
    M_{2,n} = M_{2,n-1} + (x_n - mean(x)_{n-1}) * (x_n - mean(x)_{n})
    s_n^2 = M_{2,n} / (n-1)'''
    delta = np.subtract(sample, old_mean)
    if nsam == 1:
        new_mean = sample
        new_M2 = np.zeros(np.size(sample))
        new_M2 = np.asscalar(new_M2)
        '''do so to have a list of scalars, and not a list of arrays of one element'''
    else:
        new_mean = old_mean + np.divide(delta,nsam)
        new_M2 = old_M2 + delta*np.subtract(sample,new_mean)
    return new_mean, new_M2


def compute_sample_variance_from_M2(M2,nsam):
    sample_variance = np.divide(M2,np.subtract(nsam,1))
    return sample_variance


@task(model_part_file_name=FILE_IN, parameter_file_name=FILE_IN, returns=4)
def execution_task(model_part_file_name, parameter_file_name):
    with open(parameter_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    local_parameters = parameters # in case there are more parameters file, we rename them
    model = KratosMultiphysics.Model()      

    local_parameters["solver_settings"]["model_import_settings"]["input_filename"].SetString(model_part_file_name[:-5])
    alpha = 2.0
    beta = 6.0
    sample = GenerateBetaSample(alpha,beta)
    simulation = MonteCarloAnalysis(model,local_parameters, sample)
    simulation.Run() 
    QoI =  EvaluateQuantityOfInterest(simulation)
    KratosMultiphysics.CalculateNodalAreaProcess(simulation._GetSolver().main_model_part,2).Execute()
    node_information = []
    for node in simulation._GetSolver().main_model_part.Nodes:
        node_information.append([node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE), node.GetSolutionStepValue(KratosMultiphysics.NODAL_AREA), node.X, node.Y])
    return 1, QoI, node_information, simulation.sample


# def execution_task(parameter_file_name, sample):
#     with open(parameter_file_name,'r') as parameter_file:
#         parameters = KratosMultiphysics.Parameters(parameter_file.read())
#     local_parameters = parameters # in case there are more parameters file, we rename them
#     model = KratosMultiphysics.Model()
#     simulation = MultilevelMonteCarloAnalysis(model,local_parameters,sample)
#     simulation.Run()
#     QoI =  EvaluateQuantityOfInterest(simulation)
# return QoI


# @task(returns=2)
def mean_task(*args):
    sum_val = 0
    amount_elems = 0
    for i in range(int(len(args) / 2)):
        sum_val = sum_val + (args[2 * i] * args[2 * i + 1])
        amount_elems = amount_elems + args[2 * i]
    return amount_elems, sum_val / amount_elems
    
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
        parameter_file_name = "/home/kratos105b/Kratos/applications/MonteCarloApplication/tests/Level0/ProjectParameters.json"

    with open(parameter_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    local_parameters = parameters # in case there are more parameters file, we rename them
    '''Read the number of cycles of the Monte Carlo algorithm from the .json file'''
    instances = local_parameters["problem_data"]["number_samples"].GetInt()
    
    instances = 5
    Qlist = []
    run_results = []

    '''evaluate the exact expected value of Q (sample = 1.0)'''
    ExactExpectedValueQoI = exact_execution_task(local_parameters["solver_settings"]["model_import_settings"]["input_filename"].GetString() + ".mdpa", parameter_file_name)
    
    for instance in range (0,instances):
#        model = KratosMultiphysics.Model()      
#        sample = GenerateStochasticContribute(local_parameters)
#        simulation = MonteCarloAnalysis(model,local_parameters,sample)
#        simulation.Run()
#        QoI =  EvaluateQuantityOfInterest(simulation)
#        Qlist.append(QoI)
        run_results.append(execution_task(local_parameters["solver_settings"]["model_import_settings"]["input_filename"].GetString() + ".mdpa", parameter_file_name))
        
    
    for elem in run_results:
        #curr_elem = compss_wait_on(elem[0])
        Qlist.append(elem[0])
        Qlist.append(elem[1])

    chunk_size = 4
    ## use this second while loop if you want to append at the end of Qlist the mean value you evaluate
    while len(Qlist) > 2:
        newQlist = Qlist[2*chunk_size:]
        card, val = mean_task(*Qlist[0:2*chunk_size])
        newQlist.append(card)
        newQlist.append(val)
        Qlist = newQlist
    ## use this second while loop if you want to append at the beginning of Qlist the mean value you evaluate
    # while len(Qlist) > 2:
    #     card, val = mean_task(*Qlist[0:2*chunk_size])
    #     newQlist = []
    #     newQlist.append(card)
    #     newQlist.append(val)
    #     newQlist.extend(Qlist[2*chunk_size:])
    #     Qlist = newQlist
    # variance = CalcVariance(Qlist,mean)
    # print("Variance value QoI:",variance)
    print("Size of Qlist: " + str(len(Qlist)) + " " + str(Qlist))
    mean = Qlist[1]
    
    ## evaluate now the exact expected value of Q (sample = 1.0)
    # model = KratosMultiphysics.Model()
    # sample = 1.0
    # simulation = MonteCarloAnalysis(model,local_parameters, sample)
    # simulation.Run()
    # ExactExpectedValueQoI = 0.25 * EvaluateQuantityOfInterest(simulation)
    
    print("")
    print("!!!Evaluation of the relative error between the computed mean value and the expected value of the QoI!!!")
    relative_error = compare_mean(mean,ExactExpectedValueQoI)
    # print("Values QoI:",Qlist)
    mean = compss_wait_on(mean)
    ExactExpectedValueQoI = compss_wait_on(ExactExpectedValueQoI)
    relative_error = compss_wait_on(relative_error)
    print("stochastic mean = ",mean,"exact mean = ",ExactExpectedValueQoI)
    print("relative error: ",relative_error)
    print("")


    """ The below part evaluates the relative L2 error between the numerical solution SOLUTION(x,y,sample) and the analytical solution, also dependent on sample.
    Analytical solution available in case FORCING = sample * -432.0 * (coord_x**2 + coord_y**2 - coord_x - coord_y)"""
    # print("!!!Evaluation of the L2 relative error between the numerical and the analytical solutions!!!")
    # KratosMultiphysics.CalculateNodalAreaProcess(simulation._GetSolver().main_model_part,2).Execute()
    # error = 0.0
    # L2norm_analyticalsolution = 0.0
    # for node in simulation._GetSolver().main_model_part.Nodes:
    #     local_error = ((node.GetSolutionStepValue(Poisson.SOLUTION) - (432.0*simulation.sample*node.X*node.Y*(1-node.X)*(1-node.Y)*0.5))**2) * node.GetSolutionStepValue(KratosMultiphysics.NODAL_AREA)
    #     error = error + local_error
    #     local_analyticalsolution = (432.0*simulation.sample*node.X*node.Y*(1-node.X)*(1-node.Y)*0.5)**2 * node.GetSolutionStepValue(KratosMultiphysics.NODAL_AREA)
    #     L2norm_analyticalsolution = L2norm_analyticalsolution + local_analyticalsolution
    # error = np.sqrt(error)
    # L2norm_analyticalsolution = np.sqrt(L2norm_analyticalsolution)
    # print("L2 relative error = ", error/L2norm_analyticalsolution)
    # print("")
    
