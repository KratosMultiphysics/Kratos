import json, os, math, time
import KratosMultiphysics
from KratosMultiphysics import Parameters, Logger
import KratosMultiphysics.CompressiblePotentialFlowApplication as KCPFApp
from KratosMultiphysics.response_functions.response_function_interface import ResponseFunctionInterface
import KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_analysis as potential_flow_analysis
import KratosMultiphysics.MappingApplication

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
        KratosMultiphysics.ModelPartIO(self.auxiliary_mdpa_path, KratosMultiphysics.IO.WRITE | KratosMultiphysics.IO.MESH_ONLY).WriteModelPart( self.current_model_part)

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

    def CalculateValue(self):
        pass

    def CalculateGradient(self):
        pass

    def GetValue(self):
        return self._value

    def GetNodalGradient(self, variable):
        if variable != KratosMultiphysics.SHAPE_SENSITIVITY:
            raise RuntimeError("GetNodalGradient: No gradient for {}!".format(variable.Name))

        gradient = {node.Id : node.GetValue(variable) for node in self.current_model_part.Nodes}

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
        KratosMultiphysics.ModelPartIO(self.auxiliary_mdpa_path+"_"+str(self.sample[0]), KratosMultiphysics.IO.WRITE | KratosMultiphysics.IO.MESH_ONLY).WriteModelPart(self.primal_model_part)

        with open(self.adjoint_parameters_path,'r') as parameter_file:
            adjoint_parameters = KratosMultiphysics.Parameters( parameter_file.read() )
        # Create the adjoint solver
        adjoint_parameters = _CheckParameters(adjoint_parameters)
        adjoint_model = KratosMultiphysics.Model()

        adjoint_parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(mach)
        adjoint_parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(aoa)
        adjoint_parameters["solver_settings"]["model_import_settings"]["input_filename"].SetString(self.auxiliary_mdpa_path+"_"+str(self.sample[0]))
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
        self.project_parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["mach_infinity"].SetDouble(mach)
        self.project_parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["angle_of_attack"].SetDouble(alpha)
        super().ModifyInitialProperties()


    def EvaluateQuantityOfInterest(self):
        """
        Method evaluating the QoI of the problem: lift coefficient.
        """
        qoi_list = [self.response_function.CalculateValue(self.primal_model_part)]
        Logger.PrintInfo("StochasticAdjointResponse", " Lift Coefficient: ",qoi_list[0])

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
