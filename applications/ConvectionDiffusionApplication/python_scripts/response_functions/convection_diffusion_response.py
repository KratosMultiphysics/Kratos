"""This module contains the available response functions and their base class"""
import time as timer

import KratosMultiphysics as KM
from KratosMultiphysics import Parameters, Logger
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication # TODO
import KratosMultiphysics.ConvectionDiffusionApplication as ConvectionDiffusionApplication
from KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_analysis import ConvectionDiffusionAnalysis


class ResponseFunctionBase(object):
    """The base class for response functions. This is a copy from StructuralMechanicsApplication - move to core?
    """

    def RunCalculation(self, calculate_gradient):
        self.Initialize()
        self.InitializeSolutionStep()
        self.CalculateValue()
        if calculate_gradient:
            self.CalculateGradient()
        self.FinalizeSolutionStep()
        self.Finalize()

    def Initialize(self):
        pass

    def InitializeSolutionStep(self):
        pass

    def CalculateValue(self):
        raise NotImplementedError("CalculateValue needs to be implemented by the derived class")

    def CalculateGradient(self):
        raise NotImplementedError("CalculateGradient needs to be implemented by the derived class")

    def FinalizeSolutionStep(self):
        pass

    def Finalize(self):
        pass

    def GetValue(self):
        raise NotImplementedError("GetValue needs to be implemented by the derived class")

    def GetShapeGradient(self):
        raise NotImplementedError("GetShapeGradient needs to be implemented by the derived class")


class AdjointResponseFunction(ResponseFunctionBase):
    """Linear adjoint response function.
    - runs the primal analysis
    - primal results are transferred to adjoint model part via python
    - uses primal results to calculate value
    - uses primal results to calculate gradient by running the adjoint analysis

    Attributes
    ----------
    primal_analysis : Primal analysis object of the response function
    adjoint_analysis : Adjoint analysis object of the response function
    """
    def __init__(self, identifier, response_settings, model):
        self.identifier = identifier
        self.response_settings = response_settings

        # Create the primal solver
        with open(self.response_settings["primal_settings"].GetString(),'r') as parameter_file:
            primal_parameters = Parameters( parameter_file.read() )

        self.primal_analysis = ConvectionDiffusionAnalysis(model, primal_parameters)
        self.primal_model_part = model.GetModelPart(primal_parameters["solver_settings"]["model_part_name"].GetString())

        # Create the adjoint solver
        adjoint_parameters = self._GetAdjointParameters()

        self.adjoint_analysis = ConvectionDiffusionAnalysis(model, adjoint_parameters)
        self.adjoint_model_part = model.GetModelPart(adjoint_parameters["solver_settings"]["model_part_name"].GetString())

        self.primal_state_variables = [
            KM.CONDUCTIVITY,
            KM.TEMPERATURE,
            KM.HEAT_FLUX,
            KM.FACE_HEAT_FLUX
        ]

    def Initialize(self):
        self.primal_analysis.Initialize()
        self.adjoint_analysis.Initialize()

    def InitializeSolutionStep(self):
        # Run the primal analysis.
        Logger.PrintInfo(self._GetLabel(), "Starting primal analysis for response:", self.identifier)
        startTime = timer.time()
        if not self.primal_analysis.time < self.primal_analysis.end_time:
            self.primal_analysis.end_time += 1
        self.primal_analysis.RunSolutionLoop()
        Logger.PrintInfo(self._GetLabel(), "Time needed for solving the primal analysis = ",round(timer.time() - startTime,2),"s")

    def CalculateValue(self):
        startTime = timer.time()
        value = self._GetResponseFunctionUtility().CalculateValue(self.primal_model_part)
        Logger.PrintInfo(self._GetLabel(), "Time needed for calculating the response value = ",round(timer.time() - startTime,2),"s")

        self.primal_model_part.ProcessInfo[StructuralMechanicsApplication.RESPONSE_VALUE] = value

    def CalculateGradient(self):
        # synchronize the modelparts
        self._SynchronizeAdjointFromPrimal()
        startTime = timer.time()
        Logger.PrintInfo(self._GetLabel(), "Starting adjoint analysis for response:", self.identifier)
        if not self.adjoint_analysis.time < self.adjoint_analysis.end_time:
            self.adjoint_analysis.end_time += 1
        self.adjoint_analysis.RunSolutionLoop()
        Logger.PrintInfo(self._GetLabel(), "Time needed for solving the adjoint analysis = ",round(timer.time() - startTime,2),"s")

    def GetValue(self):
        return self.primal_model_part.ProcessInfo[StructuralMechanicsApplication.RESPONSE_VALUE]

    def GetShapeGradient(self):
        gradient = {}
        for node in self.adjoint_model_part.Nodes:
            gradient[node.Id] = node.GetSolutionStepValue(KM.SHAPE_SENSITIVITY)
        return gradient

    def Finalize(self):
        self.primal_analysis.Finalize()
        self.adjoint_analysis.Finalize()

    def _GetResponseFunctionUtility(self):
        return self.adjoint_analysis._GetSolver().response_function

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
        Logger.PrintInfo(self._GetLabel(), "Transfer primal state to adjoint model part.")
        variable_utils = KM.VariableUtils()
        for variable in self.primal_state_variables:
            variable_utils.CopyModelPartNodalVar(variable, self.primal_model_part, self.adjoint_model_part, 0)

    def _GetAdjointParameters(self):
        adjoint_settings = self.response_settings["adjoint_settings"].GetString()

        if adjoint_settings == "auto":
            raise RuntimeError("auto adjoint settings not yet implemented.")
            # Logger.PrintInfo(self._GetLabel(), "Automatic set up adjoint parameters for response:", self.identifier)

            # with open(self.response_settings["primal_settings"].GetString(),'r') as parameter_file:
            #     primal_parameters = Parameters( parameter_file.read() )

            # # check that HDF5 process is not there
            # if primal_parameters["processes"].Has("list_other_processes"):
            #     for i in range(0,primal_parameters["processes"]["list_other_processes"].size()):
            #         process = primal_parameters["processes"]["list_other_processes"][i]
            #         raise Exception("Auto setup of adjoint parameters does not support {} in list_other_processes".format(process["python_module"].GetString()))

            # # clone primal settings as base for adjoint
            # adjoint_parameters = primal_parameters.Clone()

            # # analysis settings
            # solver_settings = adjoint_parameters["solver_settings"]
            # primal_solver_type = solver_settings["solver_type"].GetString()
            # if primal_solver_type != "static":
            #     raise Exception("Auto setup of adjoint parameters does not support {} solver_type. Only available for 'static'".format(primal_solver_type))
            # solver_settings["solver_type"].SetString("adjoint_"+primal_solver_type)

            # if not solver_settings.Has("compute_reactions"):
            #     solver_settings.AddEmptyValue("compute_reactions")
            # solver_settings["compute_reactions"].SetBool(False)

            # if not solver_settings.Has("move_mesh_flag"):
            #     solver_settings.AddEmptyValue("move_mesh_flag")
            # solver_settings["move_mesh_flag"].SetBool(False)

            # if solver_settings.Has("scheme_settings"):
            #     depr_msg = '\nDEPRECATION-WARNING: "scheme_settings" is deprecated, please remove it from your json parameters.\n'
            #     Logger.PrintWarning(__name__, depr_msg)
            #     solver_settings.RemoveValue("scheme_settings")

            # if solver_settings["model_import_settings"]["input_type"].GetString() == "use_input_model_part":
            #     solver_settings["model_import_settings"]["input_type"].SetString("mdpa")
            #     if solver_settings["model_import_settings"].Has("input_filename"):
            #         file_name = solver_settings["model_import_settings"]["input_filename"].GetString()
            #     else:
            #         Logger.PrintWarning(self._GetLabel(), "Automatic adjoint settings creator assumes the model_part_name as input_filename.")
            #         solver_settings["model_import_settings"].AddEmptyValue("input_filename")
            #         file_name = solver_settings["model_part_name"].GetString()
            #     solver_settings["model_import_settings"]["input_filename"].SetString(file_name)

            # # Dirichlet conditions: change variables
            # for i in range(0,primal_parameters["processes"]["constraints_process_list"].size()):
            #     process = adjoint_parameters["processes"]["constraints_process_list"][i]
            #     variable_name = process["Parameters"]["variable_name"].GetString()
            #     process["Parameters"]["variable_name"].SetString("ADJOINT_"+variable_name)

            # # Neumann conditions - do not modify to read the same load values as in primal:

            # # Output process:
            # # TODO how to add the output process? How find out about the variables?
            # if adjoint_parameters.Has("output_processes"):
            #     Logger.PrintInfo(self._GetLabel(), "Output process is removed for adjoint analysis. To enable it define adjoint_parameters yourself.")
            #     adjoint_parameters.RemoveValue("output_processes")

            # # sensitivity settings
            # adjoint_parameters["solver_settings"].AddValue("sensitivity_settings", self.response_settings["sensitivity_settings"])

            # # response settings
            # adjoint_parameters["solver_settings"].AddValue("response_function_settings", self.response_settings)

        else: # adjoint parameters file is explicitely given - do not change it.
            with open(self.response_settings["adjoint_settings"].GetString(),'r') as parameter_file:
                adjoint_parameters = Parameters( parameter_file.read() )

        return adjoint_parameters

    def _GetLabel(self):
        type_labels = {
            "point_temperature" : "LocalTemperatureAverageResponseFunction"
        }
        response_type = self.response_settings["response_type"].GetString()
        return "Adjoint" + type_labels[response_type] + "Response"
