import KratosMultiphysics
from KratosMultiphysics import Parameters, Logger
import KratosMultiphysics.CompressiblePotentialFlowApplication as KCPFApp
import KratosMultiphysics.ShapeOptimizationApplication as KSO
from KratosMultiphysics.response_functions.response_function_interface import ResponseFunctionInterface
import KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_analysis as potential_flow_analysis
import time as timer

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

# ==============================================================================
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
        # Run the primal analysis.
        Logger.PrintInfo(self._GetLabel(), "Starting primal analysis for response:", self.identifier)
        Logger.PrintInfo("\n> Starting primal analysis for response:", self.identifier)
        startTime = timer.time()
        if not self.primal_analysis.time < self.primal_analysis.end_time:
            self.primal_analysis.end_time += 1
        self.primal_analysis.RunSolutionLoop()
        Logger.PrintInfo(self._GetLabel(), "Time needed for solving the primal analysis = ",round(timer.time() - startTime,2),"s")


    def CalculateValue(self):
        startTime = timer.time()
        self._GetResponseFunctionUtility().InitializeSolutionStep()
        value = self._GetResponseFunctionUtility().CalculateValue(self.primal_model_part)
        Logger.PrintInfo(self._GetLabel(), "Time needed for calculating the response value = ",round(timer.time() - startTime,2),"s")

        self._value = value

    def CalculateGradient(self):
        # Solving Adjoint
        # synchronize the modelparts
        self._SynchronizeAdjointFromPrimal()
        startTime = timer.time()
        Logger.PrintInfo("\n> Starting adjoint analysis for response:", self.identifier)
        if not self.adjoint_analysis.time < self.adjoint_analysis.end_time:
            self.adjoint_analysis.end_time += 1
        self.adjoint_analysis.RunSolutionLoop()
        Logger.PrintInfo("> Time needed for solving the adjoint analysis = ",round(timer.time() - startTime,2),"s")

    def GetValue(self):
        return self._value

    def GetNodalGradient(self, variable):
        if variable != KratosMultiphysics.SHAPE_SENSITIVITY:
            raise RuntimeError("GetNodalGradient: No gradient for {}!".format(variable.Name))

        gradient = {node.Id : node.GetSolutionStepValue(variable) for node in self.adjoint_model_part.Nodes}

        return gradient

    def Finalize(self):
        self.primal_analysis.Finalize()
        self.adjoint_analysis.Finalize()

    def _GetResponseFunctionUtility(self):
        return self.adjoint_analysis._GetSolver()._GetResponseFunction()

    def _SynchronizeAdjointFromPrimal(self):
        Logger.PrintInfo(self._GetLabel(), "Synchronize primal and adjoint modelpart for response:", self.identifier)

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

        # Put primal solution on adjoint model
        if self.primal_data_transfer_with_python:
            Logger.PrintInfo(self._GetLabel(), "Transfer primal state to adjoint model part.")
            variable_utils = KratosMultiphysics.VariableUtils()
            for variable in self.primal_state_variables:
                variable_utils.CopyModelPartNodalVar(variable, self.primal_model_part, self.adjoint_model_part, 0)


    def _GetAdjointParameters(self):
        with open(self.response_settings["adjoint_settings"].GetString(),'r') as parameter_file:
            adjoint_parameters = Parameters( parameter_file.read() )

        return adjoint_parameters

    def _GetLabel(self):
        type_labels = {
            "adjoint_lift_potential_jump" : "LiftPotentialJump"
        }
        response_type = self.response_settings["response_type"].GetString()
        return "Adjoint" + type_labels[response_type]  +"Response"

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
# ==============================================================================
class EmbeddedAdjointResponseFunction(ResponseFunctionInterface):
    def __init__(self, identifier, response_settings, model):
        self.identifier = identifier
        self.response_settings = response_settings
        self.skin_model_part = model.GetModelPart("MainModelPart")
        self.primal_state_variables = [KCPFApp.VELOCITY_POTENTIAL, KCPFApp.AUXILIARY_VELOCITY_POTENTIAL]
        self.primal_data_transfer_with_python = True

    def InitializeSolutionStep(self):
        with open(self.response_settings["primal_settings"].GetString(),'r') as parameter_file:
            primal_parameters = Parameters( parameter_file.read() )

        primal_parameters = self._CheckParameters(primal_parameters)
        case = primal_parameters["solver_settings"]["formulation"]["element_type"].GetString()

        primal_model = KratosMultiphysics.Model()
        this_skin_model_part = primal_model.CreateModelPart("skin")
        this_skin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        this_skin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)
        this_skin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL_SENSITIVITY)
        this_skin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SHAPE_SENSITIVITY)
        this_skin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_VAUX)
        this_skin_model_part.AddNodalSolutionStepVariable(KSO.DF1DX_MAPPED)
        skin_model_part_name = primal_parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["skin_model_part_name"].GetString()
        KratosMultiphysics.ModelPartIO(skin_model_part_name).ReadModelPart(this_skin_model_part)
        self._SynchronizeSkinToCurrent(this_skin_model_part)
        KratosMultiphysics.ModelPartIO("current_skin_model_part", KratosMultiphysics.IO.WRITE | KratosMultiphysics.IO.MESH_ONLY).WriteModelPart(this_skin_model_part)

        startTime = timer.time()
        Logger.PrintInfo(self._GetLabel(), "Starting primal analysis for response:", self.identifier)
        self.primal_analysis = potential_flow_analysis.PotentialFlowAnalysis(primal_model, primal_parameters)
        self.primal_analysis.Run()
        self.primal_model_part =  self.primal_analysis._GetSolver().main_model_part
        self.primal_model_part.RemoveSubModelPart("fluid_computational_model_part")
        self.primal_model_part.RemoveSubModelPart("trailing_edge_sub_model_part")

        if not self.primal_model_part.HasProperties(0):
            self.primal_model_part.AddProperties(KratosMultiphysics.Properties(0))
        if not self.primal_model_part.HasProperties(1):
            self.primal_model_part.AddProperties(KratosMultiphysics.Properties(1))
        KratosMultiphysics.ModelPartIO("remeshed_mdpa", KratosMultiphysics.IO.WRITE | KratosMultiphysics.IO.MESH_ONLY).WriteModelPart(self.primal_model_part)
        Logger.PrintInfo(self._GetLabel(), "Time needed for solving the primal analysis = ",round(timer.time() - startTime,2),"s")


        Logger.PrintInfo("\n> Starting adjoint analysis for response:", self.identifier)
        adjoint_model = KratosMultiphysics.Model()
        adjoint_parameters = self._CheckParameters(self._GetAdjointParameters())

        adjoint_parameters["solver_settings"]["formulation"]["penalty_coefficient"].SetDouble(0.0)
        adjoint_parameters["solver_settings"]["formulation"]["stabilization_factor"].SetDouble(0.0)
        adjoint_parameters["solver_settings"]["formulation"]["element_type"].SetString(case)
        self.this_skin_model_part = adjoint_model.CreateModelPart("skin")
        self.this_skin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.this_skin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)
        self.this_skin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL_SENSITIVITY)
        self.this_skin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SHAPE_SENSITIVITY)
        self.this_skin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_VAUX)
        self.this_skin_model_part.AddNodalSolutionStepVariable(KSO.DF1DX_MAPPED)
        skin_model_part_name = primal_parameters["processes"]["boundary_conditions_process_list"][0]["Parameters"]["skin_model_part_name"].GetString()
        KratosMultiphysics.ModelPartIO(skin_model_part_name).ReadModelPart(self.this_skin_model_part)
        self._SynchronizeSkinToCurrent(self.this_skin_model_part)
        self.adjoint_analysis = potential_flow_analysis.PotentialFlowAnalysis(adjoint_model, adjoint_parameters)
        self.adjoint_model_part = self.adjoint_analysis._GetSolver().main_model_part

        self.adjoint_analysis.Initialize()
        self._SynchronizeAdjointFromPrimal()
        self.adjoint_analysis.RunSolutionLoop()
        self.adjoint_analysis.Finalize()

    def CalculateValue(self):
        startTime = timer.time()
        self._GetResponseFunctionUtility().InitializeSolutionStep()
        value = self._GetResponseFunctionUtility().CalculateValue(self.primal_model_part)
        print("LIFT VALUE", value)
        self._value = value

    def CalculateGradient(self):
        pass

    def GetValue(self):
        return self._value

    def GetNodalGradient(self, variable):
        if variable != KratosMultiphysics.SHAPE_SENSITIVITY:
            raise RuntimeError("GetNodalGradient: No gradient for {}!".format(variable.Name))
        gradient = {}

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

        mapping_parameters = KratosMultiphysics.Parameters("""{
            "mapper_type": "nearest_element",
            "interface_submodel_part_origin": "Parts_Parts_Auto1",
            "search_radius" : 0.005,
            "echo_level" : 0
            }""")
        mapper = KratosMultiphysics.MappingApplication.MapperFactory.CreateMapper(model_part, skin_model_part,mapping_parameters)
        mapper.Map(KratosMultiphysics.NORMAL_SENSITIVITY,KratosMultiphysics.NORMAL_SENSITIVITY)
        mapper.Map(KratosMultiphysics.DISTANCE_GRADIENT,KratosMultiphysics.DISTANCE_GRADIENT)

        for node in skin_model_part.Nodes:
            distance_gradient=node.GetSolutionStepValue(KratosMultiphysics.DISTANCE_GRADIENT)

            sensitivity=node.GetSolutionStepValue(KratosMultiphysics.NORMAL_SENSITIVITY)
            shape_sensitivity =[-1*sensitivity*i for i in distance_gradient]
            node.SetSolutionStepValue(KratosMultiphysics.SHAPE_SENSITIVITY, shape_sensitivity)
            gradient[node.Id] = shape_sensitivity

        return gradient

    def _GetResponseFunctionUtility(self):
        return self.adjoint_analysis._GetSolver()._GetResponseFunction()

    def _SynchronizeSkinToCurrent(self, this_skin_model_part):
        Logger.PrintInfo(self._GetLabel(), "Synchronize skin coordinates for response:", self.identifier)

        print(len(self.skin_model_part.Nodes) , len(this_skin_model_part.Nodes) )

        if len(self.skin_model_part.Nodes) != len(this_skin_model_part.Nodes):
            raise RuntimeError("_SynchronizeAdjointFromPrimal: Model parts have a different number of nodes!")

        # TODO this should happen automatically
        for primal_node, adjoint_node in zip(self.skin_model_part.Nodes, this_skin_model_part.Nodes):
            adjoint_node.X0 = primal_node.X0
            adjoint_node.Y0 = primal_node.Y0
            adjoint_node.Z0 = primal_node.Z0
            adjoint_node.X = primal_node.X
            adjoint_node.Y = primal_node.Y
            adjoint_node.Z = primal_node.Z

    def _SynchronizeAdjointFromPrimal(self):
        Logger.PrintInfo(self._GetLabel(), "Synchronize primal and adjoint modelpart for response:", self.identifier)

        if len(self.primal_model_part.Nodes) != len(self.adjoint_model_part.Nodes):
            raise RuntimeError("_SynchronizeAdjointFromPrimal: Model parts have a different number of nodes!")

        # Put primal solution on adjoint model
        if self.primal_data_transfer_with_python:
            Logger.PrintInfo(self._GetLabel(), "Transfer primal state to adjoint model part.")
            variable_utils = KratosMultiphysics.VariableUtils()
            for variable in self.primal_state_variables:
                variable_utils.CopyModelPartNodalVar(variable, self.primal_model_part, self.adjoint_model_part, 0)

    def _GetAdjointParameters(self):
        adjoint_settings = self.response_settings["adjoint_settings"].GetString()
        with open(self.response_settings["adjoint_settings"].GetString(),'r') as parameter_file:
            adjoint_parameters = Parameters( parameter_file.read() )

        return adjoint_parameters

    def _GetLabel(self):
        type_labels = {
            "embedded_adjoint_lift_potential_jump" : "EmbeddedLiftPotentialJump"
        }
        response_type = self.response_settings["response_type"].GetString()
        return "Adjoint" + type_labels[response_type]  +"Response"

    def _CheckParameters(self, parameters):
        # if not parameters["solver_settings"].Has("reform_dofs_at_each_step") or not parameters["solver_settings"]["reform_dofs_at_each_step"].GetBool():
        #     if not parameters["solver_settings"].Has("reform_dofs_at_each_step"):
        #         parameters["solver_settings"].AddEmptyValue("reform_dofs_at_each_step")
        #     parameters["solver_settings"]["reform_dofs_at_each_step"].SetBool(True)
        #     wrn_msg = 'This solver requires the setting reform the dofs at each step in optimization.'
        #     wrn_msg += 'The solver setting has been set to True'
        #     Logger.PrintWarning(self._GetLabel(), wrn_msg)
        return parameters