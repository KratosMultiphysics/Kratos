import json, os, math, time
import numpy as np
import KratosMultiphysics
from KratosMultiphysics import Parameters, Logger
import KratosMultiphysics.CompressiblePotentialFlowApplication as KCPFApp
import KratosMultiphysics.ShapeOptimizationApplication as KSO
from KratosMultiphysics.response_functions.response_function_interface import ResponseFunctionInterface
import KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_analysis as potential_flow_analysis
import KratosMultiphysics.MappingApplication
from scipy.interpolate import splev
from KratosMultiphysics.gid_output_process import GiDOutputProcess

# Import Kratos, XMC, PyCOMPSs API
import KratosMultiphysics.MultilevelMonteCarloApplication
import xmc
import xmc.methodDefs_momentEstimator.computeCentralMoments as mdccm
from exaqute import get_value_from_remote

def _GetModelPart(model, solver_settings):
    model_part_name = solver_settings["model_part_name"].GetString()
    if not model.HasModelPart(model_part_name):
        model_part = model.CreateModelPart(model_part_name, 2)
        domain_size = solver_settings["domain_size"].GetInt()
        if domain_size < 0:
            raise Exception('Please specify a "domain_size" >= 0!')
        model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)
    else:
        model_part = model.GetModelPart(model_part_name)

    return model_part

class AdjointResponseFunction(ResponseFunctionInterface):

    def __init__(self, identifier, response_settings, model):
        default_parameters = KratosMultiphysics.Parameters( """
            {
                "response_type": "stochastic_adjoint_lift_potential_jump",
                "risk_measure": "expected_value",
                "primal_settings": "",
                "adjoint_settings": "",
                "xmc_settings": "",
                "design_surface_sub_model_part_name": "",
                "auxiliary_mdpa_path": "auxiliary_mdpa",
                "primal_data_transfer_with_python": true,
                "output_dict_results_file_name": "",
                "output_pressure_file_path": ""
            }  """ )
        response_settings.ValidateAndAssignDefaults(default_parameters)

        self.identifier = identifier
        self.response_settings = response_settings

        if not response_settings["primal_settings"].GetString() == "":
            self.primal_settings = response_settings["primal_settings"].GetString()
        else:
            raise Exception("Please set the path to the primal parameters in \"primal_settings\"")

        if not response_settings["adjoint_settings"].GetString() == "":
            self.adjoint_settings = response_settings["adjoint_settings"].GetString()
        else:
            raise Exception("Please set the path to the adjoint parameters in \"adjoint_settings\"")

        if not response_settings["xmc_settings"].GetString() == "":
            self.xmc_settings_path = response_settings["xmc_settings"].GetString()
        else:
            raise Exception("Please set the path to the XMC parameters in \"xmc_settings\"")

        if not response_settings["design_surface_sub_model_part_name"].GetString() == "":
            self.design_surface_sub_model_part_name = response_settings["design_surface_sub_model_part_name"].GetString()
        else:
            raise Exception("Please set the name of the design surface submodelpart in \"design_surface_sub_model_part_name\"")

        self.auxiliary_mdpa_path = response_settings["auxiliary_mdpa_path"].GetString()
        self.risk_measure = response_settings["risk_measure"].GetString()

        if response_settings.Has("output_dict_results_file_name"):
            self.output_dict_results_file_name = response_settings["output_dict_results_file_name"].GetString()
            self.results_dict = {}
        else:
            self.output_dict_results_file_name = ""

        if response_settings.Has("output_pressure_file_path"):
            self.output_pressure_file_path = response_settings["output_pressure_file_path"].GetString()
        else:
            self.output_pressure_file_path = ""
        # Create the primal solver
        with open(self.response_settings["primal_settings"].GetString(),'r') as parameter_file:
            primal_parameters = Parameters( parameter_file.read() )

        primal_parameters = _CheckParameters(primal_parameters)
        if primal_parameters.Has("adjoint_parameters_path"):
            primal_parameters["adjoint_parameters_path"].SetString(self.response_settings["adjoint_settings"].GetString())
        else:
            primal_parameters.AddString("adjoint_parameters_path", self.response_settings["adjoint_settings"].GetString())
        if primal_parameters.Has("design_surface_sub_model_part_name"):
            primal_parameters["design_surface_sub_model_part_name"].SetString(self.design_surface_sub_model_part_name)
        else:
            primal_parameters.AddString("design_surface_sub_model_part_name", self.design_surface_sub_model_part_name)
        if primal_parameters.Has("auxiliary_mdpa_path"):
            primal_parameters["auxiliary_mdpa_path"].SetString(self.auxiliary_mdpa_path)
        else:
            primal_parameters.AddString("auxiliary_mdpa_path", self.auxiliary_mdpa_path)
        open(self.response_settings["primal_settings"].GetString(), 'w').write(primal_parameters.PrettyPrintJsonString())

        # Store current design
        self.current_model_part = _GetModelPart(model, primal_parameters["solver_settings"])

    def Initialize(self):

        if not self.output_pressure_file_path == "" and not os.path.exists(self.output_pressure_file_path):
            os.makedirs(self.output_pressure_file_path)

    def InitializeSolutionStep(self):
        self.current_model_part.RemoveSubModelPart("fluid_computational_model_part")
        self.step = self.current_model_part.ProcessInfo[KratosMultiphysics.STEP]
        KratosMultiphysics.ModelPartIO(self.auxiliary_mdpa_path, KratosMultiphysics.IO.WRITE | KratosMultiphysics.IO.MESH_ONLY | KratosMultiphysics.IO.SKIP_TIMER).WriteModelPart( self.current_model_part)

        initial_time = time.time()
        if not self.risk_measure == "cvar":
            self.InitDevelopXMC()
        else:
            self.InitCVARXMC()

    def InitDevelopXMC(self):

        initial_time = time.time()
        self._RunXMC()
        elapsed_time = time.time() - initial_time

        if self.risk_measure == "expected_value":
            order = 1 ; is_central = False
        elif self.risk_measure == "variance":
            order = 2 ; is_central = True

        if not self.output_dict_results_file_name == "":
            self.results_dict[self.step] = {}

        # save lift coefficient
        qoi_counter = 0
        estimator_container = [] # here we append the estimator for each index/level
        error_container = [] # here we append the variance of the estimator for each index/level
        n_samples_container = []
        for index in range (len(self.xmc_analysis.monteCarloSampler.indices)):
            self.xmc_analysis.monteCarloSampler.indices[index].qoiEstimator[qoi_counter] = get_value_from_remote(self.xmc_analysis.monteCarloSampler.indices[index].qoiEstimator[qoi_counter])
            estimator_container.append(float(get_value_from_remote(self.xmc_analysis.monteCarloSampler.indices[index].qoiEstimator[qoi_counter].value(order=order, isCentral=is_central))))
            error_container.append(float(get_value_from_remote(self.xmc_analysis.monteCarloSampler.indices[index].qoiEstimator[qoi_counter].value(order=order, isCentral=is_central, isErrorEstimationRequested=True)[1])))
            n_samples_container.append(int(get_value_from_remote(self.xmc_analysis.monteCarloSampler.indices[index].qoiEstimator[qoi_counter]._sampleCounter)))
        qoi_counter += 1
        # linearly sum estimators: this summation operation is valid for expected value and central moments
        # we refer to equation 4 of Krumscheid, S., Nobile, F., & Pisaroni, M. (2020). Quantifying uncertain system outputs via the multilevel Monte Carlo method — Part I: Central moment estimation. Journal of Computational Physics. https://doi.org/10.1016/j.jcp.2020.109466
        self._value = sum(estimator_container)
        # compute statistical error as in section 2.2 of
        # Pisaroni, M., Nobile, F., & Leyland, P. (2017). A Continuation Multi Level Monte Carlo (C-MLMC) method for uncertainty quantification in compressible inviscid aerodynamics. Computer Methods in Applied Mechanics and Engineering, 326, 20–50. https://doi.org/10.1016/j.cma.2017.07.030
        statistical_error = math.sqrt(sum(error_container))

        if not self.output_dict_results_file_name == "":
            self.results_dict[self.step]["run_time"]=elapsed_time
            self.results_dict[self.step]["number_of_samples"]=n_samples_container
            self.results_dict[self.step]["lift_coefficient"]={}
            self.results_dict[self.step]["lift_coefficient"]["risk_measure"]=self.risk_measure
            self.results_dict[self.step]["lift_coefficient"]["value"]=self._value
            self.results_dict[self.step]["lift_coefficient"]["statistical_error"]=statistical_error

        # save pressure coefficient
        pressure_dict = {}
        member = 0
        for node in self.current_model_part.GetSubModelPart(self.design_surface_sub_model_part_name).Nodes:
            estimator_container = [] # here we append contribution for each index/level
            variance_container = [] # here we append contribution for each index/level
            for index in range (len(self.xmc_analysis.monteCarloSampler.indices)):
                self.xmc_analysis.monteCarloSampler.indices[index].qoiEstimator[qoi_counter] = get_value_from_remote(self.xmc_analysis.monteCarloSampler.indices[index].qoiEstimator[qoi_counter])
                estimator_container.append(float(get_value_from_remote(self.xmc_analysis.monteCarloSampler.indices[index].qoiEstimator[qoi_counter].multiValue(order=order, component = member, isCentral=is_central))))
                variance_container.append(float(get_value_from_remote(self.xmc_analysis.monteCarloSampler.indices[index].qoiEstimator[qoi_counter].multiValue(order=2, component = member, isCentral=True))))
            pressure_coefficient = sum(estimator_container) # sum raw/central moment estimations on different indeces/levels
            variance_pressure_coefficient = sum(variance_container) # sum raw/central moment estimations on different indeces/levels
            member += 1
            pressure_dict[node.Id] = {}
            pressure_dict[node.Id]["coordinates"] = [node.X, node.Y, node.Z]
            pressure_dict[node.Id]["pressure_coefficient"] = pressure_coefficient
            pressure_dict[node.Id]["variance_pressure_coefficient"] = variance_pressure_coefficient
            node.SetValue(KratosMultiphysics.PRESSURE_COEFFICIENT, pressure_coefficient)
        qoi_counter += 1
        if not self.output_pressure_file_path == "":
            with open(self.output_pressure_file_path+"/pressure_"+str(self.step)+".json", 'w') as fp:
                json.dump(pressure_dict, fp,indent=4, sort_keys=True)

        # save shape sensitivity
        member = 0
        for node in self.current_model_part.GetSubModelPart(self.design_surface_sub_model_part_name).Nodes:
            shape_sensitivity = KratosMultiphysics.Vector(3, 0.0)
            for idim in range(3):
                estimator_container = [] # here we append contribution for each index/level
                for index in range (len(self.xmc_analysis.monteCarloSampler.indices)):
                    self.xmc_analysis.monteCarloSampler.indices[index].qoiEstimator[qoi_counter] = get_value_from_remote(self.xmc_analysis.monteCarloSampler.indices[index].qoiEstimator[qoi_counter])
                    estimator_container.append(float(get_value_from_remote(self.xmc_analysis.monteCarloSampler.indices[index].qoiEstimator[qoi_counter].multiValue(order=order, component = member, isCentral=is_central))))
                shape_sensitivity[idim] = sum(estimator_container) # sum raw/central moment estimations on different indeces/levels
                member += 1

            node.SetValue(KratosMultiphysics.SHAPE_SENSITIVITY, shape_sensitivity)


    def InitCVARXMC(self):

        self._RunCVaRXMC()
        print ("FINISHED XMC!")

        self.splineKnotsPsi = self.xmc_analysis.splineKnotsPsi
        self.splineCoeffsPsi = self.xmc_analysis.splineCoeffsPsi

        i_spline=0
        # print(self.argmin)
        ini_time = time.time()
        for node in self.current_model_part.GetSubModelPart(self.design_surface_sub_model_part_name).Nodes:
            # Hack becuase of how scipy splines works. 10 interpolation points hard coded and never changes
            # You will need to change the numbers to the left and right of *knots to the limit of the CVaR function
            shape_sensitivity = KratosMultiphysics.Vector(3, 0.0)
            if i_spline < len(self.splineKnotsPsi):
                for i_dim in range(2):
                    coeffs_psi = [*self.splineCoeffsPsi[i_spline+i_dim],0.0,0.0,0.0,0.0]
                    knots_psi = [-0.9,-0.9,-0.9,*self.splineKnotsPsi[i_spline+i_dim],-0.5,-0.5,-0.5]
                    shape_sensitivity[i_dim] = splev(self.argmin,(knots_psi, coeffs_psi,3),der=1)

                i_spline += 2
                node.SetValue(KratosMultiphysics.SHAPE_SENSITIVITY, shape_sensitivity)
            else:
                shape_sensitivity[1]=1.0
        print("time spent computing gradients from spline", time.time()-ini_time)

    def CalculateValue(self):
        pass

    def CalculateGradient(self):
        pass

    def GetValue(self):
        return self._value

    def GetNodalGradient(self, variable):
        if variable != KratosMultiphysics.SHAPE_SENSITIVITY:
            raise RuntimeError("GetNodalGradient: No gradient for {}!".format(variable.Name))

        gradient = {node.Id : node.GetValue(variable) for node in self.current_model_part.GetSubModelPart(self.design_surface_sub_model_part_name).Nodes}

        return gradient

    def Finalize(self):
        if not self.output_dict_results_file_name == "":
            with open(self.output_dict_results_file_name, 'w') as fp:
                json.dump(self.results_dict, fp,indent=4, sort_keys=True)

    def _GetLabel(self):
        type_labels = {
            "stochastic_adjoint_lift_potential_jump" : "StochasticLiftPotentialJump"
        }
        response_type = self.response_settings["response_type"].GetString()
        return "Adjoint" + type_labels[response_type]  +"Response"

    def _RunXMC(self):
        # read parameters
        with open(self.xmc_settings_path,'r') as parameter_file:
                parameters = json.load(parameter_file)

        # SolverWrapper
        parameters["solverWrapperInputDictionary"]["qoiEstimator"] = parameters["monteCarloIndexInputDictionary"]["qoiEstimator"]

        # SampleGenerator
        samplerInputDictionary = parameters["samplerInputDictionary"]
        samplerInputDictionary['randomGeneratorInputDictionary'] = parameters["randomGeneratorInputDictionary"]
        samplerInputDictionary['solverWrapperInputDictionary'] = parameters["solverWrapperInputDictionary"]

        # MonteCarloIndex
        monteCarloIndexInputDictionary = parameters["monteCarloIndexInputDictionary"]
        monteCarloIndexInputDictionary["samplerInputDictionary"] = samplerInputDictionary

        # MonoCriterion
        criteriaArray = []
        criteriaInputs = []
        for monoCriterion in (parameters["monoCriteriaInpuctDictionary"]):
            criteriaArray.append(xmc.monoCriterion.MonoCriterion(\
                parameters["monoCriteriaInpuctDictionary"][monoCriterion]["criteria"],\
                parameters["monoCriteriaInpuctDictionary"][monoCriterion]["tolerance"]))
            criteriaInputs.append([parameters["monoCriteriaInpuctDictionary"][monoCriterion]["input"]])

        # MultiCriterion
        multiCriterionInputDictionary=parameters["multiCriterionInputDictionary"]
        multiCriterionInputDictionary["criteria"] = criteriaArray
        multiCriterionInputDictionary["inputsForCriterion"] = criteriaInputs
        criterion = xmc.multiCriterion.MultiCriterion(**multiCriterionInputDictionary)

        # ErrorEstimator
        errorEstimator = xmc.errorEstimator.ErrorEstimator(**parameters["errorEstimatorInputDictionary"])

        # HierarchyOptimiser
        hierarchyCostOptimiser = xmc.hierarchyOptimiser.HierarchyOptimiser(**parameters["hierarchyOptimiserInputDictionary"])

        # EstimationAssembler
        assemblers = []
        if "expectationAssembler" in parameters["estimationAssemblerInputDictionary"].keys():
            expectationAssembler = xmc.estimationAssembler.EstimationAssembler(**parameters["estimationAssemblerInputDictionary"]["expectationAssembler"])
            assemblers.append(expectationAssembler)
        if "discretizationErrorAssembler" in parameters["estimationAssemblerInputDictionary"].keys():
            discretizationErrorAssembler = xmc.estimationAssembler.EstimationAssembler(**parameters["estimationAssemblerInputDictionary"]["discretizationErrorAssembler"])
            assemblers.append(discretizationErrorAssembler)
        if "varianceAssembler" in parameters["estimationAssemblerInputDictionary"].keys():
            varianceAssembler = xmc.estimationAssembler.EstimationAssembler(**parameters["estimationAssemblerInputDictionary"]["varianceAssembler"])
            assemblers.append(varianceAssembler)

        # MonteCarloSampler
        monteCarloSamplerInputDictionary = parameters["monteCarloSamplerInputDictionary"]
        monteCarloSamplerInputDictionary["indexConstructorDictionary"] = monteCarloIndexInputDictionary
        monteCarloSamplerInputDictionary["assemblers"] = assemblers
        monteCarloSamplerInputDictionary["errorEstimators"] = [errorEstimator]
        mcSampler = xmc.monteCarloSampler.MonteCarloSampler(**monteCarloSamplerInputDictionary)

        # XMCAlgorithm
        XMCAlgorithmInputDictionary = parameters["XMCAlgorithmInputDictionary"]
        XMCAlgorithmInputDictionary["monteCarloSampler"] = mcSampler
        XMCAlgorithmInputDictionary["hierarchyOptimiser"] = hierarchyCostOptimiser
        XMCAlgorithmInputDictionary["stoppingCriterion"] = criterion

        self.xmc_analysis = xmc.XMCAlgorithm(**XMCAlgorithmInputDictionary)

        if (parameters["solverWrapperInputDictionary"]["asynchronous"] is True):
            self.xmc_analysis.runAsynchronousXMC()
        else:
            self.xmc_analysis.runXMC()
    def _RunCVaRXMC(self):

        # read parameters
        with open(self.xmc_settings_path,'r') as parameter_file:
                parameters = json.load(parameter_file)

        significance = parameters["significance"]
        errorTolerance = parameters["errorTolerance"]
        tolRefinementRatio = 1.5
        Nref = int(parameters["Nref"])
        # iterationBounds = [Nref,Nref+3]
        iterationBounds = [0,0]
        toleranceSequence = [errorTolerance*(tolRefinementRatio**(Nref-i-1)) for i in range(Nref)]
        print(toleranceSequence)
        toleranceSplitting = [.001, .8] # fractions on interpolation and bias errors, respectively
        if "optiParameters" in parameters.keys():
            optiParameters = parameters["optiParameters"]
        else:
            optiParameters = self.current_model_part.NumberOfNodes()*2
        # optiParameters = 2
        optiWeights = [1.0]*(optiParameters+1)# THIS SHOULD BE LEN() == optiparameters
        derivationOrder = 1
        indexSpace = [0,3]
        parameterPointsSpace = [10,10**3+1]
        parameterPoints = np.linspace(-0.9,-0.5,num=parameterPointsSpace[0]).tolist() #TRY to set this interval as small as possible.
        sampleNumberSpace = [10,10**5]
        initialHierarchy = parameters["initialHierarchy"]
        for i, _ in enumerate(initialHierarchy):
            initialHierarchy[i].append(parameterPoints)
        sampleNumberSpace = list(np.ceil(sampleNumberSpace).astype(int))

        # RandomGeneratorWrapper
        randomGeneratorInputDict=parameters["randomGeneratorInputDictionary"]
        # SolverWrapper
        # Parameters are [c_mean,r_mean,c_std,r_std,c_time,r_time],
        # marc: READING FROM JSON
        solverWrapperInputDict = parameters["solverWrapperInputDictionary"]
        solverWrapperInputDict["outputBatchSize"] = optiParameters+1

        # SampleGenerator
        # marc: READING FROM JSON
        samplerInputDict =  parameters["samplerInputDictionary"]
        samplerInputDict["randomGeneratorInputDict"] = randomGeneratorInputDict
        samplerInputDict["solverWrapperInputDict"] = solverWrapperInputDict
                # SolverWrapper

        # Indexwise Estimators
        # Two for the function Phi, two for each Psi_i, so (optiParameters+1)*2
        # Ordering is [ParamEst(Phi), KDEst(Phi), ParamEst(Psi_1), KDEst(Psi_1),...,ParamEst(Psi_N), KDEst(Psi_N)]
        qoiEstPaths = [None]*(optiParameters+1)*2
        qoiEstDicts = [None]*(optiParameters+1)*2
        iParamEstFunction = 0
        iKDEstFunction = 1

        qoiEstDicts[iParamEstFunction] = {'evaluator':'xmc.methodDefs_parametricMomentEstimator.evaluator.allInOneEvaluator',
                                        'evaluationArguments':[significance],
                                        'parameterPoints':parameterPoints,
                                        'derivationOrder':derivationOrder}
        qoiEstPaths[iParamEstFunction] = 'xmc.parametricMomentEstimator.ParametricExpectationEstimator'
        qoiEstDicts[iKDEstFunction] = {'bandwidthMethod':'scott',
                                        'analyticFunction':'xmc.methodDefs_kernelDensityEstimator.analyticFunction.allInOne',
                                        'analyticFunctionParameters':significance,
                                        'evaluationPoints':parameterPoints}
        qoiEstPaths[iKDEstFunction] = 'xmc.kernelDensityEstimator.KernelDensityEstimator'

        for i in range(optiParameters):
            qoiEstDicts[2+2*i] = {'evaluator':'xmc.methodDefs_parametricMomentEstimator.evaluator.allInOneEvaluator',
                                        'evaluationArguments':[significance],
                                        'parameterPoints':parameterPoints,
                                        'derivationOrder':derivationOrder}
            qoiEstPaths[2+2*i] = 'xmc.parametricSensitivityEstimator.ParametricSensitivityEstimator'

            qoiEstDicts[3+2*i] = {'bandwidthMethod':'scott',
                            'analyticFunction':'xmc.methodDefs_kernelDensityEstimator.analyticFunction.allInOneJustPositive',
                            'analyticFunctionParameters':significance,
                            'evaluationPoints':parameterPoints}
            qoiEstPaths[3+2*i] = 'xmc.kernelDensityEstimator.KernelDensityEstimator'

        qoiForEstimatorsDict = [[0],[0]]
        for i in range(optiParameters):
            qoiForEstimatorsDict.append([0,i+1])
            qoiForEstimatorsDict.append([0,i+1])

        # # MonteCarloIndex Constructor
        # mci_inputDict = parameters["monteCarloIndexInputDictionary"]
        # mci_inputDict["samplerInputDict"] = samplerInputDict
        # mci_inputDict["qoiForEstimators"] = qoiForEstimatorsDict
        # mci_inputDict["qoiEstimator"] = qoiEstPaths
        # mci_inputDict["qoiEstimatorInputDict"] = qoiEstDicts

            # MonteCarloIndex Constructor
        mci_inputDict = {'indexValue': None,
                     'sampler':'xmc.sampleGenerator.SampleGenerator',
                     'samplerInputDict':samplerInputDict,
                     #'qoiForEstimators':[[0]]*2,
                     'qoiForEstimators':qoiForEstimatorsDict,
                     'qoiEstimator':qoiEstPaths,
                     'qoiEstimatorInputDict':qoiEstDicts,
                     'costEstimator':'xmc.momentEstimator.MomentEstimator',
                     'costEstimatorInputDict':{'order':1},
                     'areSamplesRecycledi':False}

        solverWrapperInputDict["qoiEstimator"] = mci_inputDict["qoiEstimator"]

        ################################# RUN TIME GENERATED ENTITIES END HERE #######################

        # MonoCriterion
        criteriaArray = [None]*3
        criteriaInputs = [None]*3
        iTotalErrCrit = 0
        #iInterpErrCrit = 0
        #iBiasErrCrit = 1
        #iStatErrCrit = 2
        criteriaArray[iTotalErrCrit] = xmc.monoCriterion.MonoCriterion(
            'xmc.methodDefs_monoCriterion.criterionFunctions.isLowerThanOrEqualTo',
            errorTolerance)
        criteriaInputs[iTotalErrCrit] = ['error'+str(iTotalErrCrit)]

        # min iterations criterion
        criteriaArray.append(xmc.monoCriterion.MonoCriterion(
            'xmc.methodDefs_monoCriterion.criterionFunctions.isGreaterThanOrEqualTo',
            iterationBounds[0]))
        criteriaInputs.append(['algorithmCost'])
        # max iterations criterion
        criteriaArray.append(xmc.monoCriterion.MonoCriterion(
            'xmc.methodDefs_monoCriterion.criterionFunctions.isGreaterThanOrEqualTo',
            iterationBounds[1]))
        criteriaInputs.append(['algorithmCost'])

        # MultiCriterion
        criterionDict = {'criteria':criteriaArray,
                        'inputsForCriterion':criteriaInputs,
    #                     'toleranceToSplit':errorTolerance,
    #                     'splitCriteria':[iInterpErrCrit,iBiasErrCrit,iStatErrCrit],
                        'interpreter':'xmc.methodDefs_multiCriterion.interpreter.interpretAsMultipleRequiredConvergencesAndIterationBounds',
                        'flag':'xmc.methodDefs_multiCriterion.flag.plainFlag'}
        criterion = xmc.multiCriterion.MultiCriterion(**criterionDict)

        # ModelEstimator
        genericModelDict = {'valueForParameters':
                            'xmc.methodDefs_modelEstimator.valueForParameters.geometricModel',
                            'updater': 'xmc.methodDefs_modelEstimator.update.updatePredictorGeometric_Task'}

        qoiPredictors = [None]*(optiParameters+1)*2
        estForPred = [[None]]*(optiParameters+1)*2
        iBiasPredictorFunction = 0
        iVariancePredictorFunction = 1

        qoiPredictors[iBiasPredictorFunction] = xmc.modelEstimator.ModelEstimator(**genericModelDict)
        estForPred[iBiasPredictorFunction] = [iParamEstFunction,[13],xmc.tools.returnInput]
        qoiPredictors[iVariancePredictorFunction] = xmc.modelEstimator.ModelEstimator(**genericModelDict)
        estForPred[iVariancePredictorFunction] = [iParamEstFunction,[3],xmc.tools.returnInput]

        for i in range(optiParameters):
            qoiPredictors[2+2*i] = xmc.modelEstimator.ModelEstimator(**genericModelDict)
            estForPred[2+2*i] = [2+2*i,[13],xmc.tools.returnInput]
            qoiPredictors[3+2*i] = xmc.modelEstimator.ModelEstimator(**genericModelDict)
            estForPred[3+2*i] = [2+2*i,[3],xmc.tools.returnInput]

        costPredictor = xmc.modelEstimator.ModelEstimator(**genericModelDict)
        costEstForPred = [1,True,False]

        # HierarchyOptimiser
        inputDict = {'indexSpace': indexSpace,
                    'sampleNumberSpace':sampleNumberSpace,
                    'parameterPointsSpace':parameterPointsSpace,
                    'tolerances':toleranceSequence,
                    'toleranceRefinementFactor':1.1,
                    'toleranceSplittingBounds': [[s]*2 for s in toleranceSplitting],
                    'derivationOrder':derivationOrder,
                    'sensParameters':optiParameters,
    #                 'optimalIndexSet': 'xmc.methodDefs_hierarchyOptimiser.optimalIndexSet.incrementLevelsByOne',
    #                 'optimalSampleNumbers': 'xmc.methodDefs_hierarchyOptimiser.optimalSampleNumbers.multiLevelDoubleAllSamples',
                    'optimalIndexSet': 'xmc.methodDefs_hierarchyOptimiser.optimalIndexSet.derivativeInterpolationMSEPosterioriWeighted',
                    'optimalSampleNumbers': 'xmc.methodDefs_hierarchyOptimiser.optimalSampleNumbers.derivativeInterpolationMSEBootstrapWeighted',
                    'optimalParameterPoints':'xmc.methodDefs_hierarchyOptimiser.optimalParameterPoints.derivativeInterpolationMSEWeighted',
                    'optimalToleranceSplitting':'xmc.methodDefs_hierarchyOptimiser.optimalToleranceSplitting.returnCurrent',
                    'defaultHierarchy': initialHierarchy,
                    'isVarianceBlended':False
        }
        hierarchyCostOptimiser = xmc.hierarchyOptimiser.HierarchyOptimiser(**inputDict)

        # EstimationAssembler
        mcsAssemblers = [None]*(3+3*optiParameters)
        estForAssblr = [None]*(3+3*optiParameters)
        iCVaRAssembler = 0
        iDerivativeNormAssemblerFunction = 1
        iBiasAssemblerFunction = 2
        iDerivativeNormAssemblerSens = list(range(3,3+2*optiParameters,2))
        iBiasAssemblerSens = list(range(4,3+2*optiParameters,2))
        iSensAssembler = list(range(3+2*optiParameters,3+3*optiParameters))

        # Assemble the actual VaR and CVaR from Phi
        pointwiseAssemblerSchematics = {'constructor':'xmc.estimationAssembler.EstimationAssembler',
                                        'assembleEstimation':'xmc.methodDefs_estimationAssembler.assembleEstimation.assembleValue_Task'}
        interpolatorSchematics = {'constructor':'xmc.interpolator.SplineInterpolator',
                                'domain':[parameterPoints[0],parameterPoints[-1]]}
        cvarAssemblerInputDict = {'interpolator':interpolatorSchematics,
                                'pointwiseAssembler':pointwiseAssemblerSchematics,
                                'post-processing':['argmin','min']}
        mcsAssemblers[iCVaRAssembler] = xmc.parametricEstimationInterpolator.ParametricEstimationInterpolator(**cvarAssemblerInputDict)
        estForAssblr[iCVaRAssembler] = [[iParamEstFunction,[0]]]

        # norm of 4th deriveative and bias for Phi
        mcsAssemblers[iDerivativeNormAssemblerFunction] = xmc.estimationAssembler.EstimationAssembler(
            assembleEstimation=
            'xmc.methodDefs_estimationAssembler.assembleEstimation.normOfDerivative',
            parameters=[4,0,1])
        estForAssblr[iDerivativeNormAssemblerFunction] = [[iKDEstFunction,[2,10**2+6]]]
        mcsAssemblers[iBiasAssemblerFunction] = xmc.estimationAssembler.EstimationAssembler(
            assembleEstimation=
            'xmc.methodDefs_estimationAssembler.assembleEstimation.assembleBias_Task')
        estForAssblr[iBiasAssemblerFunction] = [[iParamEstFunction,[13]]]

        # norm of 4th derivative and bias for Psi_i
        for i in range(optiParameters):
            mcsAssemblers[3+2*i] = xmc.estimationAssembler.EstimationAssembler(
                assembleEstimation=
                'xmc.methodDefs_estimationAssembler.assembleEstimation.normOfDerivative',
                parameters=[4,0,1])
            estForAssblr[3+2*i] = [[3+2*i,[3,10**2+6]]]

            mcsAssemblers[4+2*i] = xmc.estimationAssembler.EstimationAssembler(
                assembleEstimation=
                'xmc.methodDefs_estimationAssembler.assembleEstimation.assembleBias_Task')
            estForAssblr[4+2*i] = [[2+2*i,[13]]]

        for i in range(len(iSensAssembler)):
            # Assemble the sensitivity functions Psi_i
            pointwiseAssemblerSchematics = {'constructor':'xmc.estimationAssembler.EstimationAssembler',
                                            'assembleEstimation':'xmc.methodDefs_estimationAssembler.assembleEstimation.assembleValue_Task'}
            interpolatorSchematics = {'constructor':'xmc.interpolator.SplineInterpolator',
                                    'domain':[parameterPoints[0],parameterPoints[-1]]}
            cvarAssemblerInputDict = {'interpolator':interpolatorSchematics,
                                    'pointwiseAssembler':pointwiseAssemblerSchematics,
                                    'post-processing':['argmin','min']}

            j = iSensAssembler[i]
            mcsAssemblers[j] = xmc.parametricEstimationInterpolator.ParametricEstimationInterpolator(**cvarAssemblerInputDict)
            estForAssblr[j] = [[2+2*i,[0]]]
        # ErrorEstimator
        mcsErrEst = [None]*(1+4+optiParameters*4)
        errorOrders = [None]*(1+4+optiParameters*4)

        iInterpErrPhi = 0
        iBiasErrPhi = 1
        iStatErrPhi = 2
        iTotalErrPhi = 3

        iInterpErrPsi = list(range(4,4+optiParameters*4,4))
        iBiasErrPsi = list(range(5,4+optiParameters*4,4))
        iStatErrPsi = list(range(6,4+optiParameters*4,4))
        iTotalErrPsi = list(range(7,4+optiParameters*4,4))

        iTotalErr = 4*(optiParameters+1)

        ############ PHI1 ERROR ###################
        mcsErrEst[iInterpErrPhi] = xmc.errorEstimator.ErrorEstimatorQ(
            error='xmc.methodDefs_errorEstimator.errorEstimation.derivativeUniformInterpolationMSE',
            parameters=1)
        errorOrders[iInterpErrPhi] = {'manifest':'order',
                                'interpolatorConstants':{'method':'interpolatorConstants',
                                                            'arguments':[iCVaRAssembler,1]},
                                'numberOfInterpolationPoints':{'method':'numberOfParameterPoints',
                                                                'arguments':[0]},
                                'derivativeNorm':{'method':'estimation',
                                                    'arguments':[iDerivativeNormAssemblerFunction]},
                                'interpolatorDomain':{'method':'interpolatorDomain','arguments':[iCVaRAssembler]}}

        mcsErrEst[iBiasErrPhi] = xmc.errorEstimator.ErrorEstimatorQ(
            error='xmc.methodDefs_errorEstimator.errorEstimation.uniformlyInterpolatedDerivativeBiasMSEPosteriori',
            parameters=1)
        errorOrders[iBiasErrPhi] = {'manifest':'order',
                                'bias':{'method':'estimation','arguments':[iBiasAssemblerFunction]},
                                'biasParameters':{'method':'predictorParameters',
                                                    'arguments':[iBiasPredictorFunction]}}

        mcsErrEst[iStatErrPhi] = xmc.errorEstimator.ErrorEstimatorQ(
            error='xmc.methodDefs_errorEstimator.errorEstimation.uniformlyInterpolatedDerivativeVarianceMSEBootstrap',
            parameters=1)
        errorOrders[iStatErrPhi] = {'manifest':'order',
                                'bootstrapSamples':{'method':'indexEstimation','arguments':[iParamEstFunction,[9]]},
                                'parameterPoints':{'method':'parameterPoints','arguments':[iParamEstFunction]}}

        mcsErrEst[iTotalErrPhi] = xmc.errorEstimator.ErrorEstimatorQ(
            error='xmc.methodDefs_errorEstimator.errorEstimation.derivativeUniformTotalError',
            parameters=1)

        errorOrders[iTotalErrPhi] = {'manifest':'order',
                                'errors':{'method':'errorEstimation','arguments':[[iInterpErrPhi,iBiasErrPhi,iStatErrPhi]]}}

        ############ PSI1 ERROR ###################
        for i in range(optiParameters):
            mcsErrEst[iInterpErrPsi[i]] = xmc.errorEstimator.ErrorEstimatorQ(
                error='xmc.methodDefs_errorEstimator.errorEstimation.derivativeUniformInterpolationMSE',
                parameters=1)
            errorOrders[iInterpErrPsi[i]] = {'manifest':'order',
                                    'interpolatorConstants':{'method':'interpolatorConstants',
                                                                'arguments':[iCVaRAssembler,1]},
                                    'numberOfInterpolationPoints':{'method':'numberOfParameterPoints',
                                                                    'arguments':[0]},
                                    'derivativeNorm':{'method':'estimation',
                                                        'arguments':[3+2*i]},
                                    'interpolatorDomain':{'method':'interpolatorDomain','arguments':[iCVaRAssembler]}}

            mcsErrEst[iBiasErrPsi[i]] = xmc.errorEstimator.ErrorEstimatorQ(
                error='xmc.methodDefs_errorEstimator.errorEstimation.uniformlyInterpolatedDerivativeBiasMSEPosteriori',
                parameters=1)
            errorOrders[iBiasErrPsi[i]] = {'manifest':'order',
                                    'bias':{'method':'estimation','arguments':[4+2*i]},
                                    'biasParameters':{'method':'predictorParameters',
                                                        'arguments':[0]}}

            mcsErrEst[iStatErrPsi[i]] = xmc.errorEstimator.ErrorEstimatorQ(
                error='xmc.methodDefs_errorEstimator.errorEstimation.uniformlyInterpolatedDerivativeVarianceMSEBootstrap',
                parameters=derivationOrder)
            errorOrders[iStatErrPsi[i]] = {'manifest':'order',
                                    'bootstrapSamples':{'method':'indexEstimation','arguments':[2+2*i,[9]]},
                                    'parameterPoints':{'method':'parameterPoints','arguments':[0]}}

            mcsErrEst[iTotalErrPsi[i]] = xmc.errorEstimator.ErrorEstimatorQ(
                error='xmc.methodDefs_errorEstimator.errorEstimation.derivativeUniformTotalError',
                parameters=1)

            errorOrders[iTotalErrPsi[i]] = {'manifest':'order',
                                    'errors':{'method':'errorEstimation','arguments':[[iInterpErrPsi[i],iBiasErrPsi[i],iStatErrPsi[i]]]}}

        ############ GRADIENT ERROR ###################
        mcsErrEst[iTotalErr] = xmc.errorEstimator.ErrorEstimatorQ(
            error='xmc.methodDefs_errorEstimator.errorEstimation.derivativeUniformTotalErrorWeighted',
            parameters=optiWeights)

        errorOrders[iTotalErr] = {'manifest':'order',
                                'errors':{'method':'errorEstimation','arguments':[[iTotalErrPhi, *iTotalErrPsi]]}}

        # MonteCarloSampler
        samplerInputDict = {'indices': [],
                            'indexConstructor':'xmc.monteCarloIndex.MonteCarloIndex',
                            'indexConstructorDictionary':mci_inputDict,
                            'assemblers': mcsAssemblers,
                            'estimatorsForAssembler': estForAssblr,
                            'qoiPredictor': qoiPredictors,
                            'estimatorsForPredictor': estForPred,
                            'costEstimatorsForPredictor': [1, True, False],
                            'costPredictor': costPredictor,
                            'errorEstimators': mcsErrEst,
                            'assemblersForError': [[]]*len(mcsErrEst), # Deprecated
                            'errorDataOrders':errorOrders,
                            'isCostUpdated':True,
        }
        mcSampler = xmc.monteCarloSampler.MonteCarloSampler(**samplerInputDict)

        # XMCAlgorithm
        folder = os.getcwd()+'/output_function_sensitivity_potential_flow_new_cmlmc_'+str(errorTolerance)
        print(folder)

        indexEstimationsForHierarchy = [[0,[13]], [0,[3]]]
        for i in range(optiParameters):
            indexEstimationsForHierarchy.append([2+2*i,[13]])
            indexEstimationsForHierarchy.append([2+2*i,[3]])

        for i in range(optiParameters+1):
            indexEstimationsForHierarchy.append([2*i,[9]])

        errForCrit = [None]*1
        errForCrit[0] = iTotalErr
        algoInputDict = {
            'monteCarloSampler': mcSampler,
            'hierarchyOptimiser': hierarchyCostOptimiser,
            'stoppingCriterion': criterion,
            'errorsForStoppingCriterion': errForCrit,
            'predictorsForHierarchy': list(range(2*(optiParameters+1))),# order important
            'indexEstimationsForHierarchy': indexEstimationsForHierarchy,# order important
            'globalEstimationsForHierarchy': list(range(1,3+2*optiParameters,2)),
            'errorParametersForHierarchy': [len(mcsErrEst)-1],
            'costEstimatorForHierarchy': [1,True,False],
            'assemblersForHierarchy':[iCVaRAssembler],
            'tolerancesForHierarchy':[-1],
            'outputFolderPath':folder,
            'isDataDumped': True,
            'toleranceSplitting':toleranceSplitting
        }
        self.xmc_analysis = xmc.XMCAlgorithm(**algoInputDict)
        self.xmc_analysis.runXMC()

        self._value = self.xmc_analysis.estimation(0)[0]['min'] #armgmin is quantile, min is cvar.
        self.argmin = self.xmc_analysis.estimation(0)[0]['argmin'] #armgmin is quantile, min is cvar.

        print("OBTAINED VALUE FROM XMC", self._value)

    def _GetAdjointParameters(self):
        with open(self.response_settings["adjoint_settings"].GetString(),'r') as parameter_file:
            adjoint_parameters = Parameters( parameter_file.read() )

        return adjoint_parameters

class SimulationScenario(potential_flow_analysis.PotentialFlowAnalysis):
    def __init__(self,input_model,input_parameters,sample):
        self.sample = sample
        self.mapping = False
        self.adjoint_parameters_path =input_parameters["adjoint_parameters_path"].GetString()
        self.design_surface_sub_model_part_name = input_parameters["design_surface_sub_model_part_name"].GetString()
        self.main_model_part_name = input_parameters["solver_settings"]["model_part_name"].GetString()
        self.auxiliary_mdpa_path = input_parameters["auxiliary_mdpa_path"].GetString()

        super().__init__(input_model,input_parameters)

    def Finalize(self):

        super().Finalize()
        aoa = self.project_parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].GetDouble()
        mach = self.project_parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].GetDouble()
        self.primal_model_part = self._GetSolver().main_model_part
        nodal_velocity_process = KCPFApp.ComputeNodalValueProcess(self.primal_model_part, ["VELOCITY"])
        nodal_velocity_process.Execute()

        # Store mesh to solve with adjoint after remeshing
        self.primal_model_part.RemoveSubModelPart("fluid_computational_model_part")
        self.primal_model_part.RemoveSubModelPart("wake_sub_model_part")
        auxiliary_mdpa_path = self.auxiliary_mdpa_path+"_"+str(self.sample[0])+"_"+str(math.floor(time.time()*100000))[6:]
        if not os.path.exists(auxiliary_mdpa_path):
            KratosMultiphysics.ModelPartIO(auxiliary_mdpa_path, KratosMultiphysics.IO.WRITE | KratosMultiphysics.IO.MESH_ONLY | KratosMultiphysics.IO.SKIP_TIMER).WriteModelPart(self.primal_model_part)

        with open(self.adjoint_parameters_path,'r') as parameter_file:
            adjoint_parameters = KratosMultiphysics.Parameters( parameter_file.read() )
        # Create the adjoint solver
        adjoint_parameters = _CheckParameters(adjoint_parameters)
        adjoint_model = KratosMultiphysics.Model()

        adjoint_parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(mach)
        adjoint_parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(aoa)
        adjoint_parameters["solver_settings"]["model_import_settings"]["input_filename"].SetString(auxiliary_mdpa_path)
        self.adjoint_analysis = potential_flow_analysis.PotentialFlowAnalysis(adjoint_model, adjoint_parameters)

        self.primal_state_variables = [KCPFApp.VELOCITY_POTENTIAL, KCPFApp.AUXILIARY_VELOCITY_POTENTIAL]

        self.adjoint_analysis.Initialize()
        self.adjoint_model_part = self.adjoint_analysis._GetSolver().main_model_part

        # synchronize the modelparts
        self._SynchronizeAdjointFromPrimal()

        self.adjoint_analysis.RunSolutionLoop()
        self.adjoint_analysis.Finalize()
        self.response_function = self.adjoint_analysis._GetSolver()._GetResponseFunction()

    def ModifyInitialProperties(self):
        """
        Method introducing the stochasticity in the right hand side. Mach number and angle of attack are random varaibles.
        """
        mach = abs(self.sample[1])
        alpha = self.sample[2]
        self.project_parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(0.2+mach*0.001)
        self.project_parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(0.0)
        super().ModifyInitialProperties()


    def EvaluateQuantityOfInterest(self):
        """
        Method evaluating the QoI of the problem: lift coefficient.
        """
        qoi_list = [self.response_function.CalculateValue(self.primal_model_part)]
        print("StochasticAdjointResponse", " Lift Coefficient: ",qoi_list[0])

        pressure_coefficient = []
        nodal_value_process = KCPFApp.ComputeNodalValueProcess(self.adjoint_analysis._GetSolver().main_model_part, ["PRESSURE_COEFFICIENT"])
        nodal_value_process.Execute()
        if (self.mapping is not True):
            for node in self.adjoint_analysis._GetSolver().main_model_part.GetSubModelPart(self.design_surface_sub_model_part_name).Nodes:
                this_pressure = node.GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT)
                pressure_coefficient.append(this_pressure)

        elif (self.mapping is True):
            for node in self.mapping_reference_model.GetModelPart(self.main_model_part_name).GetSubModelPart(self.design_surface_sub_model_part_name).Nodes:
                this_pressure = node.GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT)
                pressure_coefficient.append(this_pressure)
            # Fill the rest of the list to match SHAPE_SENSITIVITY data structure length
            pressure_coefficient.extend([0.0]*self.mapping_reference_model.GetModelPart(self.main_model_part_name).GetSubModelPart(self.design_surface_sub_model_part_name).NumberOfNodes()*2)
        qoi_list.append(pressure_coefficient)

        shape_sensitivity = []
        if (self.mapping is not True):
            for node in self.adjoint_analysis._GetSolver().main_model_part.GetSubModelPart(self.design_surface_sub_model_part_name).Nodes:
                this_shape = node.GetSolutionStepValue(KratosMultiphysics.SHAPE_SENSITIVITY)
                shape_sensitivity.extend(this_shape)

        elif (self.mapping is True):
            for node in self.mapping_reference_model.GetModelPart(self.main_model_part_name).GetSubModelPart(self.design_surface_sub_model_part_name).Nodes:
                this_shape = node.GetValue(KratosMultiphysics.SHAPE_SENSITIVITY)
                shape_sensitivity.extend(this_shape)
        qoi_list.append(shape_sensitivity)
        Logger.PrintInfo("StochasticAdjointResponse", "Total number of QoI:",len(qoi_list))
        return qoi_list

    def MappingAndEvaluateQuantityOfInterest(self):

        nodal_value_process = KCPFApp.ComputeNodalValueProcess(self.adjoint_analysis._GetSolver().main_model_part, ["PRESSURE_COEFFICIENT"])
        nodal_value_process.Execute()

        KratosMultiphysics.VariableUtils().SetNonHistoricalVariableToZero(KratosMultiphysics.PRESSURE_COEFFICIENT, self.mapping_reference_model.GetModelPart(self.main_model_part_name).Nodes)
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariableToZero(KratosMultiphysics.SHAPE_SENSITIVITY, self.mapping_reference_model.GetModelPart(self.main_model_part_name).Nodes)

        # map from current model part of interest to reference model part
        mapping_parameters = KratosMultiphysics.Parameters("""{
            "mapper_type": "nearest_element",
            "echo_level" : 0
            }""")
        mapping_parameters.AddString("interface_submodel_part_origin", self.design_surface_sub_model_part_name)
        mapping_parameters.AddString("interface_submodel_part_destination", self.design_surface_sub_model_part_name)
        mapper = KratosMultiphysics.MappingApplication.MapperFactory.CreateMapper(self.adjoint_analysis._GetSolver().main_model_part,self.mapping_reference_model.GetModelPart(self.main_model_part_name),mapping_parameters)
        mapper.Map(KratosMultiphysics.PRESSURE_COEFFICIENT, \
            KratosMultiphysics.PRESSURE_COEFFICIENT,        \
            KratosMultiphysics.MappingApplication.Mapper.FROM_NON_HISTORICAL |     \
            KratosMultiphysics.MappingApplication.Mapper.TO_NON_HISTORICAL)
        mapper.Map(KratosMultiphysics.SHAPE_SENSITIVITY, \
            KratosMultiphysics.SHAPE_SENSITIVITY,
            KratosMultiphysics.MappingApplication.Mapper.TO_NON_HISTORICAL)
        # evaluate qoi
        qoi_list = self.EvaluateQuantityOfInterest()
        return qoi_list

    def _SynchronizeAdjointFromPrimal(self):

        if len(self.primal_model_part.Nodes) != len(self.adjoint_model_part.Nodes):
            raise RuntimeError("_SynchronizeAdjointFromPrimal: Model parts have a different number of nodes!")

        # TODO this should happen automatically
        for primal_node, adjoint_node in zip(self.primal_model_part.Nodes, self.adjoint_model_part.Nodes):
            adjoint_node.X0 = primal_node.X0
            adjoint_node.Y0 = primal_node.Y0
            adjoint_node.Z0 = primal_node.Z0
            adjoint_node.X = primal_node.X
            adjoint_node.Y = primal_node.Y
            adjoint_node.Z = primal_node.Z

        variable_utils = KratosMultiphysics.VariableUtils()
        for variable in self.primal_state_variables:
            variable_utils.CopyModelPartNodalVar(variable, self.primal_model_part, self.adjoint_model_part, 0)

class EmbeddedSimulationScenario(potential_flow_analysis.PotentialFlowAnalysis):
    def __init__(self,input_model,input_parameters,sample):
        self.sample = sample
        self.mapping = False
        self.adjoint_parameters_path =input_parameters["adjoint_parameters_path"].GetString()
        self.design_surface_sub_model_part_name = input_parameters["design_surface_sub_model_part_name"].GetString()
        self.main_model_part_name = input_parameters["solver_settings"]["model_part_name"].GetString()
        self.auxiliary_mdpa_path = input_parameters["auxiliary_mdpa_path"].GetString()
        self.model = input_model
        super().__init__(self.model,input_parameters)

    def Initialize(self):
        if self.project_parameters["solver_settings"]["model_import_settings"]["input_type"].GetString()=="use_input_model_part":
            self.project_parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["remeshing_flag"].SetBool(False)
            # if not self.project_parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"].Has("perform_moving"):
            #     self.project_parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"].AddBool("perform_moving", False)
            # else:
            #     self.project_parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["perform_moving"].SetBool(False)
            premain_model_part = self.model.GetModelPart("model")
            KratosMultiphysics.ModelPartIO("failed_preprimal_mdpa", KratosMultiphysics.IO.WRITE | KratosMultiphysics.IO.MESH_ONLY | KratosMultiphysics.IO.SKIP_TIMER).WriteModelPart(premain_model_part)
            self.model.DeleteModelPart("skin")
        super().Initialize()

    def Finalize(self):

        super().Finalize()
        aoa = self.project_parameters["processes"]["boundary_conditions_process_list"][self.i_boundary]["Parameters"]["angle_of_attack"].GetDouble()
        mach = self.project_parameters["processes"]["boundary_conditions_process_list"][self.i_boundary]["Parameters"]["mach_infinity"].GetDouble()
        self.primal_model_part = self._GetSolver().main_model_part
        self.skin_model_part = self.model.GetModelPart(self.design_surface_sub_model_part_name)

        nodal_velocity_process = KCPFApp.ComputeNodalValueProcess(self.primal_model_part, ["VELOCITY", "PRESSURE_COEFFICIENT"])
        nodal_velocity_process.Execute()

        # Store mesh to solve with adjoint after remeshing
        self.primal_model_part.RemoveSubModelPart("fluid_computational_model_part")
        self.primal_model_part.RemoveSubModelPart("wake_sub_model_part")
        if not self.primal_model_part.HasProperties(0):
            self.primal_model_part.AddProperties(KratosMultiphysics.Properties(0))
        if not self.primal_model_part.HasProperties(1):
            self.primal_model_part.AddProperties(KratosMultiphysics.Properties(1))

        auxiliary_mdpa_path = self.auxiliary_mdpa_path+"_"+str(self.sample[0])+"_"+str(math.floor(time.time()*100000))[6:]
        if not os.path.exists(auxiliary_mdpa_path):
            KratosMultiphysics.ModelPartIO(auxiliary_mdpa_path, KratosMultiphysics.IO.WRITE | KratosMultiphysics.IO.MESH_ONLY | KratosMultiphysics.IO.SKIP_TIMER).WriteModelPart(self.primal_model_part)


        with open(self.adjoint_parameters_path,'r') as parameter_file:
            adjoint_parameters = KratosMultiphysics.Parameters( parameter_file.read() )
        # Create the adjoint solver
        adjoint_parameters = _CheckParameters(adjoint_parameters)
        adjoint_model = KratosMultiphysics.Model()

        adjoint_parameters["processes"]["boundary_conditions_process_list"][self.i_boundary]["Parameters"]["mach_infinity"].SetDouble(mach)
        adjoint_parameters["processes"]["boundary_conditions_process_list"][self.i_boundary]["Parameters"]["angle_of_attack"].SetDouble(aoa)
        adjoint_parameters["solver_settings"]["formulation"]["element_type"].SetString(self.project_parameters["solver_settings"]["formulation"]["element_type"].GetString())
        adjoint_parameters["solver_settings"]["model_import_settings"]["input_filename"].SetString(auxiliary_mdpa_path)
        if not adjoint_parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"].Has("perform_moving"):
            adjoint_parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"].AddBool("perform_moving", False)
        else:
            adjoint_parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["perform_moving"].SetBool(False)

        self.this_skin_model_part = adjoint_model.CreateModelPart("skin")
        self.this_skin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.this_skin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)
        self.this_skin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL_SENSITIVITY)
        self.this_skin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SHAPE_SENSITIVITY)
        self.this_skin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_VAUX)
        self.this_skin_model_part.AddNodalSolutionStepVariable(KSO.DF1DX_MAPPED)
        skin_model_part_name = self.project_parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["skin_model_part_name"].GetString()
        KratosMultiphysics.ModelPartIO(skin_model_part_name).ReadModelPart(self.this_skin_model_part)
        self._SynchronizeSkinToCurrent(self.this_skin_model_part)

        self.adjoint_analysis = potential_flow_analysis.PotentialFlowAnalysis(adjoint_model, adjoint_parameters)

        self.primal_state_variables = [KCPFApp.VELOCITY_POTENTIAL, KCPFApp.AUXILIARY_VELOCITY_POTENTIAL]

        self.adjoint_analysis.Initialize()
        self.adjoint_model_part = self.adjoint_analysis._GetSolver().main_model_part

        # synchronize the modelparts
        self._SynchronizeAdjointFromPrimal()

        self.adjoint_analysis.RunSolutionLoop()
        self.adjoint_analysis.Finalize()
        self.response_function = self.adjoint_analysis._GetSolver()._GetResponseFunction()
        self.response_function.InitializeSolutionStep()

        # from KratosMultiphysics.gid_output_process import GiDOutputProcess
        # gid_output = GiDOutputProcess(
        #         self.adjoint_model_part,
        #         "adjoint_"+str(self.sample[0]),
        #         KratosMultiphysics.Parameters("""
        #             {
        #                 "result_file_configuration" : {
        #                     "gidpost_flags": {
        #                         "GiDPostMode": "GiD_PostBinary",
        #                         "MultiFileFlag": "SingleFile"
        #                     },
        #                     "nodal_results"       : ["VELOCITY_POTENTIAL","AUXILIARY_VELOCITY_POTENTIAL","GEOMETRY_DISTANCE"],
        #                     "nodal_nonhistorical_results": ["VELOCITY","METRIC_TENSOR_2D","TEMPERATURE","DISTANCE","TRAILING_EDGE"],
        #                     "gauss_point_results" : ["WAKE","KUTTA", "PRESSURE_COEFFICIENT"],
        #                     "nodal_flags_results": [],
        #                     "elemental_conditional_flags_results": ["TO_SPLIT","THERMAL","STRUCTURE"]
        #                 }
        #             }
        #             """)
        #         )
        # gid_output.ExecuteInitialize()
        # gid_output.ExecuteBeforeSolutionLoop()
        # gid_output.ExecuteInitializeSolutionStep()
        # gid_output.PrintOutput()
        # gid_output.ExecuteFinalizeSolutionStep()
        # gid_output.ExecuteFinalize()
        # stop

    def ModifyInitialProperties(self):
        """
        Method introducing the stochasticity in the right hand side. Mach number and angle of attack are random varaibles.
        """
        mach = abs(self.sample[1])
        alpha = self.sample[2]
        # print("MACH ALPHA", mach, alpha)
        for i_boundary, boundary_process in enumerate(self.project_parameters["processes"]["boundary_conditions_process_list"]):
            if boundary_process["python_module"].GetString() == "apply_far_field_process":
                boundary_process["Parameters"]["mach_infinity"].SetDouble(mach)
                boundary_process["Parameters"]["angle_of_attack"].SetDouble(alpha)
                self.i_boundary = i_boundary
                break
        super().ModifyInitialProperties()


    def EvaluateQuantityOfInterest(self):
        """
        Method evaluating the QoI of the problem: lift coefficient.
        """
        qoi_list = [self.response_function.CalculateValue(self.primal_model_part)]
        print("StochasticAdjointResponse", " Lift Coefficient: ",qoi_list[0], "Number of nodes", self.primal_model_part.NumberOfNodes())

        model_part = self.adjoint_model_part
        skin_model_part = self.this_skin_model_part

        KratosMultiphysics.ParallelDistanceCalculator2D().CalculateDistances(model_part,
            KratosMultiphysics.CompressiblePotentialFlowApplication.GEOMETRY_DISTANCE,
            KratosMultiphysics.NODAL_AREA,
            10,
            2.0)

        local_gradient = KratosMultiphysics.ComputeNodalGradientProcess2D(model_part,
        KratosMultiphysics.CompressiblePotentialFlowApplication.GEOMETRY_DISTANCE,
        KratosMultiphysics.DISTANCE_GRADIENT,
        KratosMultiphysics.NODAL_AREA)
        local_gradient.Execute()

        nodal_value_process = KCPFApp.ComputeNodalValueProcess(self.adjoint_analysis._GetSolver().main_model_part, ["PRESSURE_COEFFICIENT"])
        nodal_value_process.Execute()

        mapping_parameters = KratosMultiphysics.Parameters("""{
            "mapper_type": "nearest_element",
            "interface_submodel_part_origin": "Parts_Parts_Auto1",
            "search_radius" : 0.005,
            "echo_level" : 0
            }""")
        mapper = KratosMultiphysics.MappingApplication.MapperFactory.CreateMapper(model_part, skin_model_part,mapping_parameters)
        mapper.Map(KratosMultiphysics.NORMAL_SENSITIVITY,KratosMultiphysics.NORMAL_SENSITIVITY)
        mapper.Map(KratosMultiphysics.DISTANCE_GRADIENT,KratosMultiphysics.DISTANCE_GRADIENT)

        mapping_parameters = KratosMultiphysics.Parameters("""{
            "mapper_type": "nearest_element",
            "interface_submodel_part_origin": "Parts_Parts_Auto1",
            "search_radius" : 0.005,
            "echo_level" : 0
            }""")
        primal_mapper = KratosMultiphysics.MappingApplication.MapperFactory.CreateMapper(self.primal_model_part, skin_model_part,mapping_parameters)
        primal_mapper.Map(KratosMultiphysics.PRESSURE_COEFFICIENT, \
            KratosMultiphysics.PRESSURE_COEFFICIENT,        \
            KratosMultiphysics.MappingApplication.Mapper.FROM_NON_HISTORICAL |     \
            KratosMultiphysics.MappingApplication.Mapper.TO_NON_HISTORICAL)

        pressure_coefficient = []

        if (self.mapping is not True):
            for node in skin_model_part.Nodes:
                this_pressure = node.GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT)
                pressure_coefficient.append(this_pressure)
            pressure_coefficient.extend([0.0]*skin_model_part.NumberOfNodes()*2)

        elif (self.mapping is True):
            raise(Exception("XMC mapping is NOT needed in embedded, as the skin stays the same"))
        qoi_list.append(pressure_coefficient)

        shape_sensitivity = []
        if (self.mapping is not True):
            for node in skin_model_part.Nodes:
                distance_gradient=node.GetSolutionStepValue(KratosMultiphysics.DISTANCE_GRADIENT)
                sensitivity=node.GetSolutionStepValue(KratosMultiphysics.NORMAL_SENSITIVITY)
                this_shape_sensitivity =[-1*sensitivity*i for i in distance_gradient]
                shape_sensitivity.extend(this_shape_sensitivity)
        elif (self.mapping is True):
            raise(Exception("XMC mapping is NOT needed in embedded, as the skin stays the same"))

        qoi_list.append(shape_sensitivity)
        Logger.PrintInfo("StochasticAdjointResponse", "Total number of QoI:",len(qoi_list))
        return qoi_list

    def MappingAndEvaluateQuantityOfInterest(self):
        raise(Exception("XMC mapping is NOT needed in embedded, as the skin stays the same"))

    def _SynchronizeAdjointFromPrimal(self):

        if len(self.primal_model_part.Nodes) != len(self.adjoint_model_part.Nodes):
            raise RuntimeError("_SynchronizeAdjointFromPrimal: Model parts have a different number of nodes!")

        # TODO this should happen automatically
        for primal_node, adjoint_node in zip(self.primal_model_part.Nodes, self.adjoint_model_part.Nodes):
            adjoint_node.X0 = primal_node.X0
            adjoint_node.Y0 = primal_node.Y0
            adjoint_node.Z0 = primal_node.Z0
            adjoint_node.X = primal_node.X
            adjoint_node.Y = primal_node.Y
            adjoint_node.Z = primal_node.Z

        variable_utils = KratosMultiphysics.VariableUtils()
        for variable in self.primal_state_variables:
            variable_utils.CopyModelPartNodalVar(variable, self.primal_model_part, self.adjoint_model_part, 0)

    def _SynchronizeSkinToCurrent(self, this_skin_model_part):
        Logger.PrintInfo("EmbeddedSimulationScenario", "Synchronizing skin coordinates")

        self.current_model_part = self.skin_model_part
        if len(self.current_model_part.Nodes) != len(this_skin_model_part.Nodes):
            raise RuntimeError("_SynchronizeAdjointFromPrimal: Model parts have a different number of nodes!")

        # TODO this should happen automatically
        for primal_node, adjoint_node in zip(self.current_model_part.Nodes, this_skin_model_part.Nodes):
            adjoint_node.X0 = primal_node.X0
            adjoint_node.Y0 = primal_node.Y0
            adjoint_node.Z0 = primal_node.Z0
            adjoint_node.X = primal_node.X
            adjoint_node.Y = primal_node.Y
            adjoint_node.Z = primal_node.Z

class EmbeddedCVaRSimulationScenario(potential_flow_analysis.PotentialFlowAnalysis):
    def __init__(self,input_model,input_parameters,sample):
        self.sample = sample
        self.mapping = False
        self.adjoint_parameters_path =input_parameters["adjoint_parameters_path"].GetString()
        self.design_surface_sub_model_part_name = input_parameters["design_surface_sub_model_part_name"].GetString()
        self.main_model_part_name = input_parameters["solver_settings"]["model_part_name"].GetString()
        self.auxiliary_mdpa_path = input_parameters["auxiliary_mdpa_path"].GetString()
        self.model = input_model
        super().__init__(self.model,input_parameters)

    def Initialize(self):
        if self.project_parameters["solver_settings"]["model_import_settings"]["input_type"].GetString()=="use_input_model_part":
            self.project_parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["remeshing_flag"].SetBool(False)
            if not self.project_parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"].Has("perform_moving"):
                self.project_parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"].AddBool("perform_moving", False)
            else:
                self.project_parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["perform_moving"].SetBool(False)

        super().Initialize()

        if self.project_parameters["solver_settings"]["model_import_settings"]["input_type"].GetString()=="use_input_model_part":
            copy_model_part = self.model.GetModelPart(self.main_model_part_name)
            # copy_model_part.RemoveSubModelPart("fluid_computational_model_part")
            if not copy_model_part.HasProperties(0):
                copy_model_part.AddProperties(KratosMultiphysics.Properties(0))
            if not copy_model_part.HasProperties(1):
                copy_model_part.AddProperties(KratosMultiphysics.Properties(1))
            self.auxiliary_mdpa_path_new = self.auxiliary_mdpa_path+"_"+str(self.sample[0])+"_"+str(math.floor(time.time()*100000))[6:]
            if not os.path.exists(self.auxiliary_mdpa_path_new):
                KratosMultiphysics.ModelPartIO(self.auxiliary_mdpa_path_new, KratosMultiphysics.IO.WRITE | KratosMultiphysics.IO.MESH_ONLY | KratosMultiphysics.IO.SKIP_TIMER).WriteModelPart(copy_model_part)
            #time.sleep(1)
            print("SOLVING mdpa:", self.auxiliary_mdpa_path_new)
    def Finalize(self):

        super().Finalize()
        aoa = self.project_parameters["processes"]["boundary_conditions_process_list"][self.i_boundary]["Parameters"]["angle_of_attack"].GetDouble()
        mach = self.project_parameters["processes"]["boundary_conditions_process_list"][self.i_boundary]["Parameters"]["mach_infinity"].GetDouble()
        self.primal_model_part = self._GetSolver().main_model_part
        self.skin_model_part = self.model.GetModelPart(self.design_surface_sub_model_part_name)
        nodal_velocity_process = KCPFApp.ComputeNodalValueProcess(self.primal_model_part, ["VELOCITY", "PRESSURE_COEFFICIENT"])
        nodal_velocity_process.Execute()

        for elem in self.primal_model_part.Elements:
            elem.Set(KratosMultiphysics.ACTIVE, True)
        # gid_output = GiDOutputProcess(
        #         self.primal_model_part,
        #         "gid_output/primal_"+str(self.sample[0])+"_"+str(self.primal_model_part.NumberOfNodes()),
        #         KratosMultiphysics.Parameters("""
        #             {
        #                 "result_file_configuration" : {
        #                     "gidpost_flags": {
        #                         "GiDPostMode": "GiD_PostBinary",
        #                         "MultiFileFlag": "SingleFile"
        #                     },
        #                     "gauss_point_results" : ["VELOCITY", "PRESSURE_COEFFICIENT"]
        #                 }
        #             }
        #             """)
        #         )
        # gid_output.ExecuteInitialize()
        # gid_output.ExecuteBeforeSolutionLoop()
        # gid_output.ExecuteInitializeSolutionStep()
        # gid_output.PrintOutput()
        # gid_output.ExecuteFinalizeSolutionStep()
        # gid_output.ExecuteFinalize()

        # Store mesh to solve with adjoint after remeshing
        self.primal_model_part.RemoveSubModelPart("fluid_computational_model_part")
        self.primal_model_part.RemoveSubModelPart("wake_sub_model_part")
        if not self.primal_model_part.HasProperties(0):
            self.primal_model_part.AddProperties(KratosMultiphysics.Properties(0))
        if not self.primal_model_part.HasProperties(1):
            self.primal_model_part.AddProperties(KratosMultiphysics.Properties(1))

        auxiliary_mdpa_path = self.auxiliary_mdpa_path+"_"+str(self.sample[0])+"_"+str(math.floor(time.time()*100000))[6:]
        if not os.path.exists(auxiliary_mdpa_path):
            KratosMultiphysics.ModelPartIO(auxiliary_mdpa_path, KratosMultiphysics.IO.WRITE | KratosMultiphysics.IO.MESH_ONLY | KratosMultiphysics.IO.SKIP_TIMER).WriteModelPart(self.primal_model_part)
        #time.sleep(1)

        with open(self.adjoint_parameters_path,'r') as parameter_file:
            adjoint_parameters = KratosMultiphysics.Parameters( parameter_file.read() )
        # Create the adjoint solver
        adjoint_parameters = _CheckParameters(adjoint_parameters)
        adjoint_model = KratosMultiphysics.Model()

        adjoint_parameters["processes"]["boundary_conditions_process_list"][self.i_boundary]["Parameters"]["mach_infinity"].SetDouble(mach)
        adjoint_parameters["processes"]["boundary_conditions_process_list"][self.i_boundary]["Parameters"]["angle_of_attack"].SetDouble(aoa)
        adjoint_parameters["solver_settings"]["formulation"]["element_type"].SetString(self.project_parameters["solver_settings"]["formulation"]["element_type"].GetString())
        adjoint_parameters["solver_settings"]["model_import_settings"]["input_filename"].SetString(auxiliary_mdpa_path)
        self.this_skin_model_part = adjoint_model.CreateModelPart("skin")
        self.this_skin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.this_skin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)
        self.this_skin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL_SENSITIVITY)
        self.this_skin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SHAPE_SENSITIVITY)
        self.this_skin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_VAUX)
        self.this_skin_model_part.AddNodalSolutionStepVariable(KSO.DF1DX_MAPPED)
        skin_model_part_name = self.project_parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["skin_model_part_name"].GetString()
        KratosMultiphysics.ModelPartIO(skin_model_part_name).ReadModelPart(self.this_skin_model_part)
        self._SynchronizeSkinToCurrent(self.this_skin_model_part)

        self.adjoint_analysis = potential_flow_analysis.PotentialFlowAnalysis(adjoint_model, adjoint_parameters)

        self.primal_state_variables = [KCPFApp.VELOCITY_POTENTIAL, KCPFApp.AUXILIARY_VELOCITY_POTENTIAL]

        self.adjoint_analysis.Initialize()
        self.adjoint_model_part = self.adjoint_analysis._GetSolver().main_model_part

        # synchronize the modelparts
        self._SynchronizeAdjointFromPrimal()

        self.adjoint_analysis.RunSolutionLoop()
        self.adjoint_analysis.Finalize()
        self.response_function = self.adjoint_analysis._GetSolver()._GetResponseFunction()
        self.response_function.InitializeSolutionStep()


    def ModifyInitialProperties(self):
        """
        Method introducing the stochasticity in the right hand side. Mach number and angle of attack are random varaibles.
        """
        mach = abs(self.sample[1])
        alpha = self.sample[2]
        KratosMultiphysics.Logger.PrintInfo("MACH ALPHA", mach, alpha)
        for i_boundary, boundary_process in enumerate(self.project_parameters["processes"]["boundary_conditions_process_list"]):
            if boundary_process["python_module"].GetString() == "apply_far_field_process":
                boundary_process["Parameters"]["mach_infinity"].SetDouble(mach)
                boundary_process["Parameters"]["angle_of_attack"].SetDouble(alpha)
                self.i_boundary = i_boundary
                break
        super().ModifyInitialProperties()


    def EvaluateQuantityOfInterest(self):
        """
        Method evaluating the QoI of the problem: lift coefficient.
        """
        # REPORTING NEGATIVE OF LIFT TO MINIMIZE THE NEGATIVE -> maximize lift
        qoi_list = [-1*self.response_function.CalculateValue(self.primal_model_part)]
        print("StochasticAdjointResponse", " Lift Coefficient: ",qoi_list[0], "Number of nodes", self.primal_model_part.NumberOfNodes())

        model_part = self.adjoint_model_part
        skin_model_part = self.this_skin_model_part

        KratosMultiphysics.ParallelDistanceCalculator2D().CalculateDistances(model_part,
            KratosMultiphysics.CompressiblePotentialFlowApplication.GEOMETRY_DISTANCE,
            KratosMultiphysics.NODAL_AREA,
            10,
            2.0)

        local_gradient = KratosMultiphysics.ComputeNodalGradientProcess2D(model_part,
        KratosMultiphysics.CompressiblePotentialFlowApplication.GEOMETRY_DISTANCE,
        KratosMultiphysics.DISTANCE_GRADIENT,
        KratosMultiphysics.NODAL_AREA)
        local_gradient.Execute()

        nodal_value_process = KCPFApp.ComputeNodalValueProcess(self.adjoint_analysis._GetSolver().main_model_part, ["PRESSURE_COEFFICIENT"])
        nodal_value_process.Execute()

        find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(self.adjoint_analysis._GetSolver().main_model_part)
        find_nodal_h.Execute()

        KratosMultiphysics.CompressiblePotentialFlowApplication.PotentialFlowUtilities.ScaleSensitivity(self.adjoint_analysis._GetSolver().main_model_part)

        mapping_parameters = KratosMultiphysics.Parameters("""{
            "mapper_type": "nearest_element",
            "interface_submodel_part_origin": "Parts_Parts_Auto1",
            "search_radius" : 0.005,
            "echo_level" : 0
            }""")
        mapper = KratosMultiphysics.MappingApplication.MapperFactory.CreateMapper(model_part, skin_model_part,mapping_parameters)
        mapper.Map(KratosMultiphysics.NORMAL_SENSITIVITY,KratosMultiphysics.NORMAL_SENSITIVITY)
        mapper.Map(KratosMultiphysics.DISTANCE_GRADIENT,KratosMultiphysics.DISTANCE_GRADIENT)

        mapping_parameters = KratosMultiphysics.Parameters("""{
            "mapper_type": "nearest_element",
            "interface_submodel_part_origin": "Parts_Parts_Auto1",
            "search_radius" : 0.005,
            "echo_level" : 0
            }""")
        primal_mapper = KratosMultiphysics.MappingApplication.MapperFactory.CreateMapper(self.primal_model_part, skin_model_part,mapping_parameters)
        primal_mapper.Map(KratosMultiphysics.PRESSURE_COEFFICIENT, \
            KratosMultiphysics.PRESSURE_COEFFICIENT,        \
            KratosMultiphysics.MappingApplication.Mapper.FROM_NON_HISTORICAL |     \
            KratosMultiphysics.MappingApplication.Mapper.TO_NON_HISTORICAL)

        # pressure_coefficient = []
        # nodal_value_process = KCPFApp.ComputeNodalValueProcess(self.adjoint_analysis._GetSolver().main_model_part, ["PRESSURE_COEFFICIENT"])
        # nodal_value_process.Execute()
        # if (self.mapping is not True):
            # for node in skin_model_part.Nodes:
            #     this_pressure = node.GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT)
            #     pressure_coefficient.append(this_pressure)
            # pressure_coefficient.extend([0.0]*skin_model_part.NumberOfNodes()*2)

        # elif (self.mapping is True):
        #     raise(Exception("XMC mapping is NOT needed in embedded, as the skin stays the same"))
        # qoi_list.append(pressure_coefficient)

        if (self.mapping is not True):
            for node in skin_model_part.Nodes:
                distance_gradient=node.GetSolutionStepValue(KratosMultiphysics.DISTANCE_GRADIENT)
                sensitivity=node.GetSolutionStepValue(KratosMultiphysics.NORMAL_SENSITIVITY)
                # one negative due to level set and one negative as to minimize the negative
                this_shape_sensitivity =[-1*-1*sensitivity*i for i in distance_gradient]
                qoi_list.extend(this_shape_sensitivity[0:2])
        elif (self.mapping is True):
            raise(Exception("XMC mapping is NOT needed in embedded, as the skin stays the same"))

        Logger.PrintInfo("StochasticAdjointResponse", "Total number of QoI:",len(qoi_list))
        return qoi_list

    def MappingAndEvaluateQuantityOfInterest(self):
        raise(Exception("XMC mapping is NOT needed in embedded, as the skin stays the same"))

    def _SynchronizeAdjointFromPrimal(self):

        if len(self.primal_model_part.Nodes) != len(self.adjoint_model_part.Nodes):
            raise RuntimeError("_SynchronizeAdjointFromPrimal: Model parts have a different number of nodes!")

        # TODO this should happen automatically
        for primal_node, adjoint_node in zip(self.primal_model_part.Nodes, self.adjoint_model_part.Nodes):
            adjoint_node.X0 = primal_node.X0
            adjoint_node.Y0 = primal_node.Y0
            adjoint_node.Z0 = primal_node.Z0
            adjoint_node.X = primal_node.X
            adjoint_node.Y = primal_node.Y
            adjoint_node.Z = primal_node.Z

        variable_utils = KratosMultiphysics.VariableUtils()
        for variable in self.primal_state_variables:
            variable_utils.CopyModelPartNodalVar(variable, self.primal_model_part, self.adjoint_model_part, 0)

    def _SynchronizeSkinToCurrent(self, this_skin_model_part):
        Logger.PrintInfo("EmbeddedSimulationScenario", "Synchronizing skin coordinates")

        self.current_model_part = self.skin_model_part
        if len(self.current_model_part.Nodes) != len(this_skin_model_part.Nodes):
            raise RuntimeError("_SynchronizeAdjointFromPrimal: Model parts have a different number of nodes!")

        # TODO this should happen automatically
        for primal_node, adjoint_node in zip(self.current_model_part.Nodes, this_skin_model_part.Nodes):
            adjoint_node.X0 = primal_node.X0
            adjoint_node.Y0 = primal_node.Y0
            adjoint_node.Z0 = primal_node.Z0
            adjoint_node.X = primal_node.X
            adjoint_node.Y = primal_node.Y
            adjoint_node.Z = primal_node.Z

def _CheckParameters(parameters):
    if not parameters["solver_settings"].Has("reform_dofs_at_each_step") or not parameters["solver_settings"]["reform_dofs_at_each_step"].GetBool():
        if not parameters["solver_settings"].Has("reform_dofs_at_each_step"):
            parameters["solver_settings"].AddEmptyValue("reform_dofs_at_each_step")
        parameters["solver_settings"]["reform_dofs_at_each_step"].SetBool(True)
        wrn_msg = 'This solver requires the setting reform the dofs at each step in optimization.'
        wrn_msg += 'The solver setting has been set to True'
    for subproc_keys, subproc_values in parameters["processes"].items():
        for process  in subproc_values:
            if "wake" in process["python_module"].GetString() and not "embedded" in process["python_module"].GetString():
                if not process["Parameters"].Has("compute_wake_at_each_step") or not process["Parameters"]["compute_wake_at_each_step"].GetBool():
                    if not process["Parameters"].Has("compute_wake_at_each_step"):
                        process["Parameters"].AddEmptyValue("compute_wake_at_each_step")
                process["Parameters"]["compute_wake_at_each_step"].SetBool(True)
    return parameters
