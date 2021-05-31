import KratosMultiphysics
from KratosMultiphysics import Parameters, Logger
import KratosMultiphysics.CompressiblePotentialFlowApplication as KCPFApp
from KratosMultiphysics.response_functions.response_function_interface import ResponseFunctionInterface
import KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_analysis as potential_flow_analysis
import time as timer

# Import Kratos, XMC, PyCOMPSs API
import KratosMultiphysics.MultilevelMonteCarloApplication
import xmc
import xmc.methodDefs_momentEstimator.computeCentralMoments as mdccm
from exaqute import get_value_from_remote
import json

def _GetModelPart(model, solver_settings):
    #TODO can be removed once model is fully available
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
        self.identifier = identifier
        self.response_settings = response_settings
        # Create the primal solver
        with open(self.response_settings["primal_settings"].GetString(),'r') as parameter_file:
            primal_parameters = Parameters( parameter_file.read() )

        primal_parameters = self._CheckParameters(primal_parameters)

        self.primal_model_part = _GetModelPart(model, primal_parameters["solver_settings"])

        self.primal_analysis = potential_flow_analysis.PotentialFlowAnalysis(model, primal_parameters)

        self.primal_data_transfer_with_python = self.response_settings["primal_data_transfer_with_python"].GetBool()

        # Create the adjoint solver
        adjoint_parameters = self._CheckParameters(self._GetAdjointParameters())
        adjoint_model = KratosMultiphysics.Model()
        self.adjoint_model_part = _GetModelPart(adjoint_model, adjoint_parameters["solver_settings"])

        self.adjoint_analysis = potential_flow_analysis.PotentialFlowAnalysis(adjoint_model, adjoint_parameters)

        self.primal_state_variables = [KCPFApp.VELOCITY_POTENTIAL, KCPFApp.AUXILIARY_VELOCITY_POTENTIAL]

    def Initialize(self):
        self.primal_analysis.Initialize()
        self.adjoint_analysis.Initialize()

    def InitializeSolutionStep(self):
        self.primal_model_part.RemoveSubModelPart("fluid_computational_model_part")
        this_path = "/media/kratos105a/datos/branch_cases/20210511_bodyfitted_optimization_stochastic"
        KratosMultiphysics.ModelPartIO(this_path+'/current_design', KratosMultiphysics.IO.WRITE | KratosMultiphysics.IO.MESH_ONLY).WriteModelPart( self.primal_model_part)

        self._RunXMC()

        # save lift coefficient
        for qoi_counter in range (0,1):
            for index in range (len(self.xmc_analysis.monteCarloSampler.indices)):
                self.xmc_analysis.monteCarloSampler.indices[index].qoiEstimator[qoi_counter] = get_value_from_remote(self.xmc_analysis.monteCarloSampler.indices[index].qoiEstimator[qoi_counter])
                sample_counter = self.xmc_analysis.monteCarloSampler.indices[index].qoiEstimator[qoi_counter]._sampleCounter
                S1 = float(get_value_from_remote(self.xmc_analysis.monteCarloSampler.indices[index].qoiEstimator[qoi_counter].powerSums[0][0]))
                S2 = float(get_value_from_remote(self.xmc_analysis.monteCarloSampler.indices[index].qoiEstimator[qoi_counter].powerSums[1][0]))
                h1 = float(get_value_from_remote(mdccm.computeCentralMomentsOrderOneDimensionZero(S1,sample_counter)))
                h2 = float(get_value_from_remote(mdccm.computeCentralMomentsOrderTwoDimensionZero(S1,S2,sample_counter)))
                self._value = h1

        # save shape sensitivity
        qoi_counter = 1
        member = 0
        for node in self.adjoint_model_part.GetSubModelPart("Body2D_Body").Nodes:
            # for index in range (len(self.xmc_analysis.monteCarloSampler.indices)):
            shape_sensitivity = KratosMultiphysics.Vector(3, 0.0)
            for idim in range(3):
                self.xmc_analysis.monteCarloSampler.indices[index].qoiEstimator[qoi_counter] = get_value_from_remote(self.xmc_analysis.monteCarloSampler.indices[index].qoiEstimator[qoi_counter])
                sample_counter = self.xmc_analysis.monteCarloSampler.indices[index].qoiEstimator[qoi_counter]._sampleCounter
                S1 = float(get_value_from_remote(self.xmc_analysis.monteCarloSampler.indices[index].qoiEstimator[qoi_counter]._powerSums["1"][member]))
                S2 = float(get_value_from_remote(self.xmc_analysis.monteCarloSampler.indices[index].qoiEstimator[qoi_counter]._powerSums["2"][member]))
                h1 = float(get_value_from_remote(mdccm.computeCentralMomentsOrderOneDimensionZero(S1,sample_counter)))
                h2 = float(get_value_from_remote(mdccm.computeCentralMomentsOrderTwoDimensionZero(S1,S2,sample_counter)))
                shape_sensitivity[idim] = h1
                member += 1

            node.SetSolutionStepValue(KratosMultiphysics.SHAPE_SENSITIVITY, shape_sensitivity)

    def CalculateValue(self):
        pass

    def CalculateGradient(self):
        pass

    def GetValue(self):
        return self._value

    def GetNodalGradient(self, variable):
        if variable != KratosMultiphysics.SHAPE_SENSITIVITY:
            raise RuntimeError("GetNodalGradient: No gradient for {}!".format(variable.Name))

        gradient = {node.Id : node.GetSolutionStepValue(variable) for node in self.adjoint_model_part.Nodes}

        return gradient

    def Finalize(self):
        pass

    def _GetLabel(self):
        type_labels = {
            "stochastic_adjoint_lift_potential_jump" : "StochasticLiftPotentialJump"
        }
        response_type = self.response_settings["response_type"].GetString()
        return "Adjoint" + type_labels[response_type]  +"Response"

    def _RunXMC(self):
        parametersPath = "/media/kratos105a/datos/branch_cases/20210511_bodyfitted_optimization_stochastic/parameters_xmc_asynchronous_mc_potentialFlow.json"

        # read parameters
        with open(parametersPath,'r') as parameter_file:
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
        statErrorEstimator = xmc.errorEstimator.ErrorEstimator(**parameters["errorEstimatorInputDictionary"])

        # HierarchyOptimiser
        hierarchyCostOptimiser = xmc.hierarchyOptimiser.HierarchyOptimiser(**parameters["hierarchyOptimiserInputDictionary"])

        # EstimationAssembler
        if "expectationAssembler" in parameters["estimationAssemblerInputDictionary"].keys():
            expectationAssembler = xmc.estimationAssembler.EstimationAssembler(**parameters["estimationAssemblerInputDictionary"]["expectationAssembler"])
        if "varianceAssembler" in parameters["estimationAssemblerInputDictionary"].keys():
            varianceAssembler = xmc.estimationAssembler.EstimationAssembler(**parameters["estimationAssemblerInputDictionary"]["varianceAssembler"])

        # MonteCarloSampler
        monteCarloSamplerInputDictionary = parameters["monteCarloSamplerInputDictionary"]
        monteCarloSamplerInputDictionary["indexConstructorDictionary"] = monteCarloIndexInputDictionary
        monteCarloSamplerInputDictionary["assemblers"] =  [expectationAssembler,varianceAssembler]
        monteCarloSamplerInputDictionary["errorEstimators"] = [statErrorEstimator]
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

    def _CheckParameters(self, parameters):
        if not parameters["solver_settings"].Has("reform_dofs_at_each_step") or not parameters["solver_settings"]["reform_dofs_at_each_step"].GetBool():
            if not parameters["solver_settings"].Has("reform_dofs_at_each_step"):
                parameters["solver_settings"].AddEmptyValue("reform_dofs_at_each_step")
            parameters["solver_settings"]["reform_dofs_at_each_step"].SetBool(True)
            wrn_msg = 'This solver requires the setting reform the dofs at each step in optimization.'
            wrn_msg += 'The solver setting has been set to True'
            Logger.PrintWarning(self._GetLabel(), wrn_msg)
        for subproc_keys, subproc_values in parameters["processes"].items():
            for process  in subproc_values:
                if "wake" in process["python_module"].GetString():
                    if not process["Parameters"].Has("compute_wake_at_each_step") or not process["Parameters"]["compute_wake_at_each_step"].GetBool():
                        if not process["Parameters"].Has("compute_wake_at_each_step"):
                            process["Parameters"].AddEmptyValue("compute_wake_at_each_step")
                    process["Parameters"]["compute_wake_at_each_step"].SetBool(True)
        return parameters


    def _GetAdjointParameters(self):
        with open(self.response_settings["adjoint_settings"].GetString(),'r') as parameter_file:
            adjoint_parameters = Parameters( parameter_file.read() )

        return adjoint_parameters