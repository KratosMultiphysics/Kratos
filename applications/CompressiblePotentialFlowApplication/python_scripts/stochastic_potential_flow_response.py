import KratosMultiphysics
from KratosMultiphysics import Parameters, Logger
import KratosMultiphysics.CompressiblePotentialFlowApplication as KCPFApp
from KratosMultiphysics.response_functions.response_function_interface import ResponseFunctionInterface
import KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_analysis as potential_flow_analysis

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
        self.xmc_settings_path = response_settings["xmc_settings"].GetString()
        self.design_surface_sub_model_part_name = response_settings["design_surface_sub_model_part_name"].GetString()
        self.auxiliary_mdpa_path = response_settings["auxiliary_mdpa_path"].GetString()
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
        open(self.response_settings["primal_settings"].GetString(), 'w').write(primal_parameters.PrettyPrintJsonString())

        self.primal_model_part = _GetModelPart(model, primal_parameters["solver_settings"])

        self.primal_analysis = potential_flow_analysis.PotentialFlowAnalysis(model, primal_parameters)

        # Create the adjoint solver
        adjoint_parameters = _CheckParameters(self._GetAdjointParameters())
        adjoint_model = KratosMultiphysics.Model()
        self.adjoint_model_part = _GetModelPart(adjoint_model, adjoint_parameters["solver_settings"])

        self.adjoint_analysis = potential_flow_analysis.PotentialFlowAnalysis(adjoint_model, adjoint_parameters)

        self.primal_state_variables = [KCPFApp.VELOCITY_POTENTIAL, KCPFApp.AUXILIARY_VELOCITY_POTENTIAL]

    def Initialize(self):
        self.primal_analysis.Initialize()
        self.adjoint_analysis.Initialize()

    def InitializeSolutionStep(self):
        self.primal_model_part.RemoveSubModelPart("fluid_computational_model_part")
        KratosMultiphysics.ModelPartIO(self.auxiliary_mdpa_path, KratosMultiphysics.IO.WRITE | KratosMultiphysics.IO.MESH_ONLY).WriteModelPart( self.primal_model_part)

        self._RunXMC()

        # TODO: here we should retrieve the risk measure we want to use for optimization
        #       at the moment we are using the expected value estimation
        risk_measure = "expected_value"
        if risk_measure == "expected_value":
            order = 1 ; is_central = False
        elif risk_measure == "variance":
            order = 2 ; is_central = True

        # save lift coefficient
        qoi_counter = 0
        estimator_container = [] # here we append contribution for each index/level
        for index in range (len(self.xmc_analysis.monteCarloSampler.indices)):
            self.xmc_analysis.monteCarloSampler.indices[index].qoiEstimator[qoi_counter] = get_value_from_remote(self.xmc_analysis.monteCarloSampler.indices[index].qoiEstimator[qoi_counter])
            estimator_container.append(float(get_value_from_remote(self.xmc_analysis.monteCarloSampler.indices[index].qoiEstimator[qoi_counter].value(order=order, isCentral=is_central))))

        # linearly sum estimators: this summation operation is valid for expected value and central moments
        # we refer to equation 4 of Krumscheid, S., Nobile, F., & Pisaroni, M. (2020). Quantifying uncertain system outputs via the multilevel Monte Carlo method — Part I: Central moment estimation. Journal of Computational Physics. https://doi.org/10.1016/j.jcp.2020.109466
        self._value = sum(estimator_container)

        # save shape sensitivity
        qoi_counter = 1
        member = 0
        for node in self.adjoint_model_part.GetSubModelPart(self.design_surface_sub_model_part_name).Nodes:
            shape_sensitivity = KratosMultiphysics.Vector(3, 0.0)
            for idim in range(3):
                estimator_container = [] # here we append contribution for each index/level
                for index in range (len(self.xmc_analysis.monteCarloSampler.indices)):
                    self.xmc_analysis.monteCarloSampler.indices[index].qoiEstimator[qoi_counter] = get_value_from_remote(self.xmc_analysis.monteCarloSampler.indices[index].qoiEstimator[qoi_counter])
                    estimator_container.append(float(get_value_from_remote(self.xmc_analysis.monteCarloSampler.indices[index].qoiEstimator[qoi_counter].multiValue(order=order, component = member, isCentral=is_central))))
                shape_sensitivity[idim] = sum(estimator_container) # sum raw/central moment estimations on different indeces/levels
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
        super(SimulationScenario,self).__init__(input_model,input_parameters)

    def Finalize(self):

        super().Finalize()
        aoa = self.project_parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].GetDouble()
        mach = self.project_parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].GetDouble()
        self.primal_model_part = self._GetSolver().main_model_part


        with open(self.adjoint_parameters_path,'r') as parameter_file:
            adjoint_parameters = KratosMultiphysics.Parameters( parameter_file.read() )
        # Create the adjoint solver
        adjoint_parameters = _CheckParameters(adjoint_parameters)
        adjoint_model = KratosMultiphysics.Model()

        adjoint_parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(mach)
        adjoint_parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(aoa)
        self.adjoint_analysis = potential_flow_analysis.PotentialFlowAnalysis(adjoint_model, adjoint_parameters)

        self.primal_state_variables = [KCPFApp.VELOCITY_POTENTIAL, KCPFApp.AUXILIARY_VELOCITY_POTENTIAL]

        self.adjoint_analysis.Initialize()
        self.adjoint_model_part = self.adjoint_analysis._GetSolver().main_model_part

        self._SynchronizeAdjointFromPrimal()

        self.response_function = self.adjoint_analysis._GetSolver()._GetResponseFunction()
        self.response_function.InitializeSolutionStep()
        # synchronize the modelparts
        self.adjoint_analysis.RunSolutionLoop()

    """
    function introducing the stochasticity in the right hand side
    input:  self: an instance of the class
    """
    def ModifyInitialProperties(self):
        '''Introduce here the stochasticity in the Mach number and the angle of attack'''
        mach = self.sample[1]
        alpha =  self.sample[2]
        self.project_parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(mach)
        self.project_parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(alpha)
        super(SimulationScenario,self).ModifyInitialProperties()


    """
    function evaluating the QoI of the problem: lift coefficient
    input:  self: an instance of the class
    """
    def EvaluateQuantityOfInterest(self):
        qoi_list = [self.response_function.CalculateValue(self.primal_model_part)]
        Logger.PrintInfo("StochasticAdjointResponse", " Lift Coefficient: ",qoi_list[0])

        shape_sensitivity = []
        if (self.mapping is not True):
            for node in self.adjoint_analysis._GetSolver().main_model_part.GetSubModelPart(self.design_surface_sub_model_part_name).Nodes:
                this_shape = node.GetSolutionStepValue(KratosMultiphysics.SHAPE_SENSITIVITY)
                shape_sensitivity.extend(this_shape)

        elif (self.mapping is True):
            raise(Exception(("Mapping not contemplated yet for shape sensitivity")))
        qoi_list.append(shape_sensitivity)
        Logger.PrintInfo("StochasticAdjointResponse", "Total number of QoI:",len(qoi_list))
        Logger.PrintInfo("StochasticAdjointResponse", "Total number of shape:",len(qoi_list[1]))
        Logger.PrintInfo("StochasticAdjointResponse", "Total number of exampleshape:",qoi_list[1][0])
        return qoi_list

    """
    function mapping the pressure field on reference model
    input:  self: an instance of the class
    """
    def MappingAndEvaluateQuantityOfInterest(self):
        Logger.PrintInfo("StochasticAdjointResponse",  "Start Mapping")
        # map from current model part of interest to reference model part
        mapping_parameters = KratosMultiphysics.Parameters("""{
            "mapper_type": "nearest_element",
            "echo_level" : 0
            }""")
        mapping_parameters.AddString("interface_sub_model_part_origin", self.design_surface_sub_model_part_name)
        mapping_parameters.AddString("interface_submodel_part_destination", self.design_surface_sub_model_part_name)
        mapper = KratosMultiphysics.MappingApplication.MapperFactory.CreateMapper(self._GetSolver().main_model_part,self.mapping_reference_model.GetModelPart("model"),mapping_parameters)
        mapper.Map(KratosMultiphysics.PRESSURE_COEFFICIENT, \
            KratosMultiphysics.PRESSURE_COEFFICIENT,        \
            KratosMultiphysics.MappingApplication.Mapper.FROM_NON_HISTORICAL |     \
            KratosMultiphysics.MappingApplication.Mapper.TO_NON_HISTORICAL)
        mapper.Map(KratosMultiphysics.SHAPE_SENSITIVITY, KratosMultiphysics.SHAPE_SENSITIVITY)
        Logger.PrintInfo("StochasticAdjointResponse", "End Mapping")
        # evaluate qoi
        Logger.PrintInfo("StochasticAdjointResponse", "Start evaluating QoI")
        qoi_list = self.EvaluateQuantityOfInterest()
        Logger.PrintInfo("StochasticAdjointResponse", "End evaluating QoI")
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

def _CheckParameters(parameters):
    if not parameters["solver_settings"].Has("reform_dofs_at_each_step") or not parameters["solver_settings"]["reform_dofs_at_each_step"].GetBool():
        if not parameters["solver_settings"].Has("reform_dofs_at_each_step"):
            parameters["solver_settings"].AddEmptyValue("reform_dofs_at_each_step")
        parameters["solver_settings"]["reform_dofs_at_each_step"].SetBool(True)
        wrn_msg = 'This solver requires the setting reform the dofs at each step in optimization.'
        wrn_msg += 'The solver setting has been set to True'
    for subproc_keys, subproc_values in parameters["processes"].items():
        for process  in subproc_values:
            if "wake" in process["python_module"].GetString():
                if not process["Parameters"].Has("compute_wake_at_each_step") or not process["Parameters"]["compute_wake_at_each_step"].GetBool():
                    if not process["Parameters"].Has("compute_wake_at_each_step"):
                        process["Parameters"].AddEmptyValue("compute_wake_at_each_step")
                process["Parameters"]["compute_wake_at_each_step"].SetBool(True)
    return parameters