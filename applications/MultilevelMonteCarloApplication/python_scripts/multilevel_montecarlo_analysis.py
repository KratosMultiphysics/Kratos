from __future__ import absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import numpy as np
import time

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.ConvectionDiffusionApplication as KratosConvDiff
import KratosMultiphysics.MultilevelMonteCarloApplication as Poisson

# Avoid printing of Kratos informations
KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING) # avoid printing of Kratos things

# Importing the base class
from analysis_stage import AnalysisStage

# Import pycompss
from pycompss.api.task import task
from pycompss.api.api import compss_wait_on
from pycompss.api.parameter import *

# Import Continuation Multilevel Monte Carlo library
import cmlmc as mlmc

class MultilevelMonteCarloAnalysis(AnalysisStage):
    '''Main script for MultilevelMonte Carlo simulations using the pure_diffusion solver'''


    def __init__(self,model,parameters,sample):
        self.sample = sample
        super(MultilevelMonteCarloAnalysis,self).__init__(model,parameters)
        self._GetSolver().main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
            
    def _CreateSolver(self):
        import convection_diffusion_stationary_solver
        solver = convection_diffusion_stationary_solver.CreateSolver(self.model,self.project_parameters["solver_settings"])
        self.LaplacianSolver = solver
        return self.LaplacianSolver

    
    def _GetSimulationName(self):
        return "Multilevel Monte Carlo Analysis"

    
    def ApplyBoundaryConditions(self):
        super(MultilevelMonteCarloAnalysis,self).ApplyBoundaryConditions()
        ## define the forcing function
        for node in self.model.GetModelPart("MLMCLaplacianModelPart").Nodes:
            coord_x = node.X
            coord_y = node.Y
            # forcing = -432.0 * coord_x * (coord_x - 1) * coord_y * (coord_y - 1)
            forcing = -432.0 * (coord_x**2 + coord_y**2 - coord_x - coord_y)
            node.SetSolutionStepValue(KratosMultiphysics.HEAT_FLUX,forcing*self.sample)
        # print("SAMPLE = ",self.sample)
            '''need to use SetSolutionStepValue and not SetValue (equivalent for GetSolutionStepValue and GetValue) because in custom_elements and custom_conditions I use this function
            Set/GetValue uses less memory and is always used for element and conditions
            Set/GetSolutionStepValue uses more memory, and allows to memorize values stored in the nodes in case of buffer > 1 (e.g. sol_{i-1}, sol_{i-2} remain stored there to evaluate sol_{i})
            if I use SetValue, I need to use GetValue (the same for the other)
            node.SetValue(Poisson.FORCING,node.GetSolutionStepValue(Poisson.FORCING)*sample)'''
            

    
    
###########################################################
######## END OF CLASS MULTILEVELMONTECARLOANALYSIS ########
###########################################################



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


'''
function executing the problem
input:
        model_part_file_name : path of the model part file (still to implement how to import in efficient way in a loop where I have different model part files and different ProjectParameters files, thus for now read model part name from the ProjectParameters.json file)
        parameter_file_name  : path of the Project Parameters file
        sample               : stochastic random variable
output:
        QoI                  : Quantity of Interest
still to implement how to import in efficient way in a loop where I have different model part files and different ProjectParameters files, thus for now read model part name from the ProjectParameters.json file
'''
# @task(model_part_file_name=FILE_IN, parameter_file_name=FILE_IN, returns=1)
@task(parameter_file_name=FILE_IN, returns=1)
def execution_task(parameter_file_name, sample):
    with open(parameter_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    local_parameters = parameters # in case there are more parameters file, we rename them
    model = KratosMultiphysics.Model()
    # local_parameters["solver_settings"]["model_import_settings"]["input_filename"].SetString(model_part_file_name[:-5])
    simulation = MultilevelMonteCarloAnalysis(model,local_parameters,sample)
    simulation.Run()
    QoI = EvaluateQuantityOfInterest(simulation)
    return QoI
    

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
    simulation = MultilevelMonteCarloAnalysis(model,local_parameters, sample)
    simulation.Run()
    ExactExpectedValueQoI = 0.25 * EvaluateQuantityOfInterest(simulation)
    # return simulation,ExactExpectedValueQoI
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
    # from sys import argv
    
    # if len(argv) > 2:
    #     err_msg = 'Too many input arguments!\n'
    #     err_msg += 'Use this script in the following way:\n'
    #     err_msg += '- With default parameter file (assumed to be called "ProjectParameters.json"):\n'
    #     err_msg += '    "python montecarlo_laplacian_analysis.py"\n'
    #     err_msg += '- With custom parameter file:\n'
    #     err_msg += '    "python montecarlo_laplacian_analysis.py <my-parameter-file>.json"\n'
    #     raise Exception(err_msg)

    # if len(argv) == 2: # ProjectParameters is being passed from outside
    #     parameter_file_name = argv[1]
    # else: # using default name
    #     parameter_file_name = "/home/kratos105b/Kratos/applications/MultilevelMonteCarloApplication/tests/Level0/ProjectParameters.json"

    '''
    set the ProjectParameters.json path in the parameter_file_name list
    '''
    parameter_file_name =[]
    parameter_file_name.append("/home/kratos105b/Kratos/applications/MultilevelMonteCarloApplication/tests/Level0/ProjectParameters.json")
    with open(parameter_file_name[0],'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    local_parameters_0 = parameters # in case there are more parameters file, we rename them
    parameter_file_name.append("/home/kratos105b/Kratos/applications/MultilevelMonteCarloApplication/tests/Level1/ProjectParameters.json")
    with open(parameter_file_name[1],'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    local_parameters_1 = parameters # in case there are more parameters file, we rename them
    parameter_file_name.append("/home/kratos105b/Kratos/applications/MultilevelMonteCarloApplication/tests/Level2/ProjectParameters.json")
    with open(parameter_file_name[2],'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    local_parameters_2 = parameters # in case there are more parameters file, we rename them
    parameter_file_name.append("/home/kratos105b/Kratos/applications/MultilevelMonteCarloApplication/tests/Level3/ProjectParameters.json")
    with open(parameter_file_name[3],'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    local_parameters_3 = parameters # in case there are more parameters file, we rename them
    parameter_file_name.append("/home/kratos105b/Kratos/applications/MultilevelMonteCarloApplication/tests/Level4/ProjectParameters.json")
    with open(parameter_file_name[4],'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    local_parameters_4 = parameters # in case there are more parameters file, we rename them
    L_max = len(parameter_file_name) - 1
    print("Maximum number of levels = ",L_max)

    '''
    evaluate the exact expected value of Q (sample = 1.0)
    need to change both local_parameters_LEVEL and parameter_file_name[LEVEL] to compute for level LEVEL
    '''
    ExactExpectedValueQoI = exact_execution_task(local_parameters_2["solver_settings"]["model_import_settings"]["input_filename"].GetString() + ".mdpa", parameter_file_name[2])
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
    # print("\n L2 relative error = ", error/L2norm_analyticalsolution,"\n")
    
    '''define setting parameters of the ML simulation'''
    settings_ML_simulation = [0.1, 0.1, 1.25, 1.15, 0.25, 0.1, 1.0, 10, 2]
    '''
    k0   = settings_ML_simulation[0] # Certainty Parameter 0 rates
    k1   = settings_ML_simulation[1] # Certainty Parameter 1 rates
    r1   = settings_ML_simulation[2] # Cost increase first iterations C-MLMC
    r2   = settings_ML_simulation[3] # Cost increase final iterations C-MLMC
    tol0 = settings_ML_simulation[4] # Tolerance iter 0
    tolF = settings_ML_simulation[5] # Tolerance final
    cphi = settings_ML_simulation[6] # Confidence on tolerance
    N0   = settings_ML_simulation[7] # Number of samples for iter 0
    L0   = settings_ML_simulation[8] # Number of levels for iter 0
    '''

    difference_QoI = [] # list containing Y_{l}^{i} = Q_{m_l} - Q_{m_{l-1}}
    time_ML = []        # list containing the time to compute the level=l simulations
    L_screening = settings_ML_simulation[8] # recall the levels start from zero
    number_samples = []
    for lev in range(0,L_screening+1):
        number_samples.append(settings_ML_simulation[7])

    if (L_screening+1) > len(difference_QoI):
        for i in range (0,(L_screening+1)-len(difference_QoI)):
            difference_QoI.append([]) # append a list in Y_l for every level
            time_ML.append([])
    print("\n ######## SCREENING PHASE ######## \n")
    
    for level in range (0,(L_screening+1)):
        for instance in range (0,number_samples[level]):
            sample = GenerateBetaSample(2.0,6.0) # generate a random variable with beta pdf, alpha = 2.0 and beta = 6.0
            run_results = []
            start_time_ML = time.time() # I can insert this operation in "GenerateBetaSample", or better to create a new function?

            if level == 0: # evaluating QoI in the coarsest grid
                run_results.append(execution_task(parameter_file_name[level], sample)) # append to run_results QoI for the coarsest grid
                time_MLi = time.time() - start_time_ML # create a new function?
                # difference_QoI[level].append(run_results[-1])
                difference_QoI[level] = np.append(difference_QoI[level],run_results[-1]) # with list[-1] we read the last element of the list
                # time_ML[level].append(time_MLi)
                time_ML[level] = np.append(time_ML[level],time_MLi)
                
            else:
                for cycle_level in range (0,level+1):
                    run_results.append(execution_task(parameter_file_name[cycle_level], sample))
                  
                time_MLi = time.time() - start_time_ML
                # difference_QoI[level].append(run_results[-1] - run_results[-2])
                difference_QoI[level] = np.append(difference_QoI[level],run_results[-1] - run_results[-2])
                # time_ML[level].append(time_MLi)
                time_ML[level] = np.append(time_ML[level],time_MLi)
       
    '''compute {E^(MC)[Y_l]} = 1/N * sum_{i=1}^{N} Y_l(sample(i))
    compute {V^(MC)[Y_l]}'''
    mean_difference_QoI = []
    if len(mean_difference_QoI) < (L_screening +1):
        for i in range (0,(L_screening+1)-len(mean_difference_QoI)):
            mean_difference_QoI.append([]) # append a list in E^(MC)[Y_l] for every level
    variance_difference_QoI = []
    if len(variance_difference_QoI) < (L_screening +1):
        for i in range (0,(L_screening+1)-len(variance_difference_QoI)):
            variance_difference_QoI.append([]) # append a list in Var^(MC)[Y_l] for every level
    second_moment_difference_QoI = []
    if len(second_moment_difference_QoI) < (L_screening +1):
        for i in range (0,(L_screening+1)-len(second_moment_difference_QoI)):
            second_moment_difference_QoI.append([]) # append a list in Var^(MC)[Y_l] for every level
    
    for level in range (0,L_screening+1):
        for i in range(0,number_samples[level]):
            nsam = i+1
            mean_difference_QoI[level],second_moment_difference_QoI[level],variance_difference_QoI[level] = mlmc.update_onepass_M(difference_QoI[level][i],mean_difference_QoI[level],second_moment_difference_QoI[level],nsam)
    # print("list Y_l",difference_QoI)
    print("mean Y_l",mean_difference_QoI)
    print("sample variance Y_l",variance_difference_QoI)

    '''compute {E^(MC)[time_ML_l]} = 1/N * sum_{i=1}^{N} time_ML_l(sample(i))
    compute {V^(MC)[time_ML_l]}'''
    mean_time_ML = []
    if len(mean_time_ML) < (L_screening +1):
        for i in range (0,(L_screening+1)-len(mean_time_ML)):
            mean_time_ML.append([]) # append a list in E^(MC)[Y_l] for every level
    variance_time_ML = []
    if len(variance_time_ML) < (L_screening +1):
        for i in range (0,(L_screening+1)-len(variance_time_ML)):
            variance_time_ML.append([]) # append a list in Var^(MC)[Y_l] for every level
    second_moment_time_ML = []
    if len(second_moment_time_ML) < (L_screening +1):
        for i in range (0,(L_screening+1)-len(second_moment_time_ML)):
            second_moment_time_ML.append([])

    for level in range (0,L_screening+1):    
        for i in range(0,number_samples[level]):
            nsam = i+1
            mean_time_ML[level],second_moment_time_ML[level],variance_time_ML[level] = mlmc.update_onepass_M(time_ML[level][i],mean_time_ML[level],second_moment_time_ML[level],nsam)
    # print("list time ML",time_ML)
    print("mean time ML",mean_time_ML)
    print("sample variance time ML",variance_time_ML)

    '''compute nDoF: number degrees of freedom for each mesh'''
    nDoF = []
    for level in range (0,L_max + 1):
        nDoF.append(mlmc.Nf_law(level))
    
    '''compute parameters by least square fit to estimate Bayesian VAR'''
    ratesLS = mlmc.compute_ratesLS(mean_difference_QoI,variance_difference_QoI,mean_time_ML,nDoF[0:L_screening+1])
    # print("rates LS computed through least square fit = ",ratesLS)

    '''compute Bayesian VAR V^c[Y_l]'''
    BayesianVariance = mlmc.EstimateBayesianVariance(mean_difference_QoI,variance_difference_QoI,settings_ML_simulation,ratesLS,nDoF,number_samples,L_screening)
    print("Bayesian Variance estimated = ", BayesianVariance)
    
    '''compute i_E, number of iterations'''
    # iE_cmlmc = np.floor((-np.log(tolF)+np.log(r2)+np.log(tol0))/(np.log(r1)))
    iE_cmlmc = mlmc.compute_iE_cmlmc(settings_ML_simulation)
    print("\nnumber of iterations we are going to perform for CMLMC = ",iE_cmlmc)

    convergence = False
    iter_MLMC = 1
    L_old = L_screening

    while convergence is not True:
        print("\n ######## CMLMC iter = ",iter_MLMC,"######## \n")
        '''Compute Tolerance for the iteration i
        eventually, we may still run the algorithm for a few more iterations wrt iE_cmlmc'''
        tol_i = mlmc.compute_tolerance_i(settings_ML_simulation,iE_cmlmc,iter_MLMC)
        
        '''Compute Optimal Number of Levels L_i'''
        # print("Bayesian variance before computing optimal number of levels",BayesianVariance)
        L_opt, BayesianVariance, L_old = mlmc.compute_levels(tol_i,number_samples,ratesLS,nDoF,BayesianVariance,
                                                        mean_difference_QoI,variance_difference_QoI,settings_ML_simulation,L_max,L_old)
        print("optimal number of levels for current iteration = ",L_opt,"previous number of levels = ",L_old)
        print("Bayesian variance after new optimal number of levels = ",BayesianVariance)
        
        '''#################################################################
        # observation: levels start from level 0                           #
        #              length arrays starts from 1                         #
        # then we have a difference of 1 between length arrays and levels  #
        # current_level = len(number_samples) - 1                          #
        # or                                                               #
        # current_level = len(difference_value) - 1                        #
        ####################################################################'''

        '''compute new theta splitting,
        since we updated the number of levels we have a new M_l
        theta_i is a scalar'''
        theta_i = mlmc.theta_model(ratesLS,tol_i,nDoF[L_opt])
        if not((theta_i > 0.0) or (theta_i < 1.0)):
            raise Exception ("The splitting parameter theta_i assumed a value outside the range (0,1)")

        '''compute number of samples according to bayesian variance and theta splitting parameters'''
        number_samples, difference_number_samples, previous_number_samples = mlmc.compute_number_samples(L_opt,BayesianVariance,ratesLS,theta_i,tol_i,nDoF,number_samples,settings_ML_simulation)
        # difference_number_samples = [2,2,2,2]
        # number_samples = [4,4,4,2]

        '''run hierarchy using optimal levels and optimal number of samples
        I use the "old" difference_QoI, and append the new values for the added samples,
        or append at the end for the new level'''
        if (L_opt+1) > len(difference_QoI):
            for i in range (0,(L_opt+1)-len(difference_QoI)):
                difference_QoI.append([]) # append a list in Y_l for the new level
        if (L_opt+1) > len(time_ML):
            for i in range (0,(L_opt+1)-len(time_ML)):
                time_ML.append([]) # append a list in time_ML for the new level
        for level in range (0,L_opt+1):
            for instance in range (0,difference_number_samples[level]):
                sample = GenerateBetaSample(2.0,6.0)
                run_results = []
                start_time_ML = time.time()
                if level == 0: # evaluating QoI in the coarsest grid
                    run_results.append(execution_task(parameter_file_name[level], sample)) # append to run_results QoI for the coarsest grid
                    time_MLi = time.time() - start_time_ML
                    # difference_QoI[level].append(run_results[-1]) # with list[-1] we read the last element of the list
                    # time_ML[level].append(time_MLi)
                    difference_QoI[level] = np.append(difference_QoI[level],run_results[-1])
                    time_ML[level] = np.append(time_ML[level],time_MLi)
                else:
                    for cycle_level in range (0,level+1):
                        run_results.append(execution_task(parameter_file_name[cycle_level], sample))
                    time_MLi = time.time() - start_time_ML
                    # difference_QoI[level].append(run_results[-1] - run_results[-2])
                    # time_ML[level].append(time_MLi)
                    difference_QoI[level] = np.append(difference_QoI[level],run_results[-1] - run_results[-2])
                    time_ML[level] = np.append(time_ML[level],time_MLi)

        # print("iteration",iter_MLMC,"Y_l",difference_QoI)
        # print("iteration",iter_MLMC,"time ML",time_ML)
    
        '''compute mean, second moment, sample variance for Y_l'''
        mean_difference_QoI = mean_difference_QoI[0:L_old+1]
        variance_difference_QoI = variance_difference_QoI[0:L_old+1]
        second_moment_difference_QoI = second_moment_difference_QoI[0:L_old+1]
        if len(mean_difference_QoI) < (L_opt+1):
            for i in range (0,(L_opt+1)-len(mean_difference_QoI)):
                mean_difference_QoI.append([]) # append a list in E^(MC)[Y_l] for every level
        if len(variance_difference_QoI) < (L_opt+1):
            for i in range (0,(L_opt+1)-len(variance_difference_QoI)):
                variance_difference_QoI.append([]) # append a list in Var^(MC)[Y_l] for every level
        if len(second_moment_difference_QoI) < (L_opt+1):
            for i in range (0,(L_opt+1)-len(second_moment_difference_QoI)):
                second_moment_difference_QoI.append([])

        for level in range (0,L_opt+1):
            for i in range(0,difference_number_samples[level]):
                nsam = previous_number_samples[level] + (i+1)
                mean_difference_QoI[level],second_moment_difference_QoI[level],variance_difference_QoI[level] = mlmc.update_onepass_M(
                    difference_QoI[level][previous_number_samples[level]+i],mean_difference_QoI[level],second_moment_difference_QoI[level],nsam)
        
        print("updated mean Y_l",mean_difference_QoI)
        print("updated sample variance Y_l",variance_difference_QoI)

        '''now compute mean, second moment, sample variance for time ML'''
        mean_time_ML = mean_time_ML[0:L_old+1]
        variance_time_ML = variance_time_ML[0:L_old+1]
        second_moment_time_ML = second_moment_time_ML[0:L_old+1]
        if len(mean_time_ML) < (L_opt+1):
            for i in range (0,(L_opt+1)-len(mean_time_ML)):
                mean_time_ML.append([]) # append a list in E^(MC)[time_ML] for every level
        if len(variance_time_ML) < (L_opt+1):
            for i in range (0,(L_opt+1)-len(variance_time_ML)):
                variance_time_ML.append([]) # append a list in Var^(MC)[time_ML] for every level
        if len(second_moment_time_ML) < (L_opt+1):
            for i in range (0,(L_opt+1)-len(second_moment_time_ML)):
                second_moment_time_ML.append([])

        for level in range (0,L_opt+1):
            for i in range(0,difference_number_samples[level]):
                nsam = previous_number_samples[level] + (i+1)
                mean_time_ML[level],second_moment_time_ML[level],variance_time_ML[level] = mlmc.update_onepass_M(
                    time_ML[level][previous_number_samples[level]+i],mean_time_ML[level],second_moment_time_ML[level],nsam)
        
        print("updated mean time ML",mean_time_ML)
        print("updated sample variance time ML",variance_time_ML)

        '''compute E^MLMC [QoI]'''
        mean_mlmc_QoI = mlmc.compute_mean_mlmc_QoI(mean_difference_QoI)
        print("updated mean MLMC QoI = ",mean_mlmc_QoI)

        '''estimation problem parameters for Bayesian updates
        compute parameters by least square fit to estimate Bayesian VAR'''
        ratesLS = mlmc.compute_ratesLS(mean_difference_QoI,variance_difference_QoI,mean_time_ML,nDoF[0:L_opt+1])

        '''compute Bayesian VAR V^c[Y_l]
        I compute from zero, since I have new ratesLS
        also Nobile computes each time from zero the Bayesian variance in "var_estim"'''
        BayesianVariance = mlmc.EstimateBayesianVariance(mean_difference_QoI,variance_difference_QoI,settings_ML_simulation,ratesLS,nDoF,number_samples,L_opt)
        print("updated Bayesian Variance estimated = ", BayesianVariance)

        '''compute total error of the MLMC simulation'''
        TErr = mlmc.compute_total_error_MLMC(mean_difference_QoI,number_samples,L_opt,BayesianVariance,settings_ML_simulation)
        print("total error TErr of current iteration",TErr)

        L_old = L_opt

        '''in [PNL17] go out of the cycle if: i) iter >= iE_cmlmc
                                             ii) TErr < tolerance_iter'''
        if iter_MLMC >= iE_cmlmc:
            if (TErr < tol_i):
                convergence = True
        else:
            iter_MLMC = iter_MLMC + 1

        
    relative_error = compare_mean(mean_mlmc_QoI,ExactExpectedValueQoI)
    # mean_mlmc_QoI = compss_wait_on(mean_mlmc_QoI)
    # ExactExpectedValueQoI = compss_wait_on(ExactExpectedValueQoI)
    # relative_error = compss_wait_on(relative_error)
    print("\niterations = ",iter_MLMC,"total error TErr computed = ",TErr,"mean MLMC QoI = ",mean_mlmc_QoI,"exact mean = ",ExactExpectedValueQoI)
    print("relative error: ",relative_error,"\n")

    '''### OBSERVATION ###
    between different tasks you don't need compss_wait_on, it's pycompss who handles everything automatically
    if the output of a task is given directly to the input of an other, pycompss handles everything'''
