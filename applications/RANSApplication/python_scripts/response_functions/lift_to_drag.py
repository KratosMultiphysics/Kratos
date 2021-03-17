"""
This module contains an interface to the available response functions
"""
import time as timer
import shutil
from pathlib import Path
from math import isclose, sqrt

import KratosMultiphysics as Kratos
from KratosMultiphysics.response_functions.response_function_interface import ResponseFunctionInterface

from KratosMultiphysics.RANSApplication.rans_analysis import RANSAnalysis
from KratosMultiphysics.RANSApplication.adjoint_rans_analysis import AdjointRANSAnalysis

try:
    import KratosMultiphysics.MeshingApplication as KMA
    from KratosMultiphysics.MeshingApplication.mmg_process import MmgProcess as automatic_remeshing_process
    if not hasattr(KMA, "MmgProcess2D"):
        automatic_remeshing_process = None
        automatic_remeshing_error_msg = "MeshingApplication is not compiled with '-DINCLUDE_MMG=ON'"
except ImportError as err:
    automatic_remeshing_process = None
    automatic_remeshing_error_msg = str(err)

class LiftToDrag(ResponseFunctionInterface):
    def __init__(self, identifier, response_settings, model):
        self.identifier = identifier
        self.response_settings = response_settings

        default_parameters = Kratos.Parameters("""
        {
            "response_type"       : "lift_to_drag_response",
            "problem_setup_folder": "PLEASE_SPECIFY_PROBLEM_SETUP_FOLDER",
            "problem_setup_files" : {
                "primal_project_parameters_file"      : "PLEASE_SPECIFY_PRIMAL_PROJECT_PARAMETERS_FILE",
                "lift_adjoint_project_parameters_file": "PLEASE_SPECIFY_LIFT_ADJOINT_PROJECT_PARAMETERS_FILE",
                "drag_adjoint_project_parameters_file": "PLEASE_SPECIFY_DRAG_ADJOINT_PROJECT_PARAMETERS_FILE"
            },
            "use_automatic_mesh_refinement"     : false,
            "automatic_mesh_refinement_settings": {
                "moving_average_number_of_iterations"         : 10,
                "moving_average_difference_relative_tolerance": 1e-3,
                "moving_average_difference_absolute_tolerance": 1e-5,
                "radius_of_influence"                         : 0.0,
                "influence_model_part_name"                   : "PLEASE_SPECIFY_INFLUENCE_MODEL_PART",
                "refinement_data_settings": {
                    "refinement_data_from": "lift_adjoint_solution",
                    "list_of_variables": [
                        "VELOCITY",
                        "PRESSURE",
                        "ADJOINT_FLUID_SCALAR_1",
                        "ADJOINT_FLUID_VECTOR_1"
                    ]
                },
                "mmg_process_settings": {}
            }
        }
        """)

        self.model = model

        self.response_settings.ValidateAndAssignDefaults(default_parameters)

        self.problem_setup_folder = Path(self.response_settings["problem_setup_folder"].GetString()).absolute()
        if (not self.problem_setup_folder.is_dir()):
            raise Exception("The provided \"problem_setup_folder\" is not a directory [ \"problem_setup_folder\" = {:s} ].".format(str(self.problem_setup_folder)))

        self.problem_setup_file_settings = self.response_settings["problem_setup_files"]

        # checks for lift response settings
        lift_adjoint_project_parameters_file = self.problem_setup_folder / self.problem_setup_file_settings["lift_adjoint_project_parameters_file"].GetString()
        with open(lift_adjoint_project_parameters_file, "r") as file_input:
            lift_adjoint_settings = Kratos.Parameters(file_input.read())

        lift_adjoint_response_function_settings = lift_adjoint_settings["solver_settings"]["response_function_settings"]
        lift_response_type = lift_adjoint_response_function_settings["response_type"].GetString()
        if (lift_response_type != "drag"):
            raise Exception("Lift to drag needs \"response_type\" in \"response_function_settings\" of \"lift_adjoint_project_parameters_file\" file to be drag. [ \"lift_adjoint_project_parameters_file\" = {:s}, \"response_type\" = {:s} ].".format(
                lift_adjoint_project_parameters_file, lift_response_type))

        # check for drag response settings
        drag_adjoint_project_parameters_file = self.problem_setup_folder / self.problem_setup_file_settings["drag_adjoint_project_parameters_file"].GetString()
        with open(drag_adjoint_project_parameters_file, "r") as file_input:
            drag_adjoint_settings = Kratos.Parameters(file_input.read())

        drag_adjoint_response_function_settings = drag_adjoint_settings["solver_settings"]["response_function_settings"]
        drag_response_type = drag_adjoint_response_function_settings["response_type"].GetString()
        if (drag_response_type != "drag"):
            raise Exception("Lift to drag needs \"response_type\" in \"response_function_settings\" of \"drag_adjoint_project_parameters_file\" file to be drag. [ \"drag_adjoint_project_parameters_file\" = {:s}, \"response_type\" = {:s} ].".format(
                drag_adjoint_project_parameters_file, drag_response_type))

        self.lift_model_part_name = lift_adjoint_response_function_settings["custom_settings"]["structure_model_part_name"].GetString()
        self.lift_direction = lift_adjoint_response_function_settings["custom_settings"]["drag_direction"].GetVector()

        self.drag_model_part_name = drag_adjoint_response_function_settings["custom_settings"]["structure_model_part_name"].GetString()
        self.drag_direction = drag_adjoint_response_function_settings["custom_settings"]["drag_direction"].GetVector()

        # check for model part name settings
        primal_project_parameters_file = self.problem_setup_folder / self.problem_setup_file_settings["primal_project_parameters_file"].GetString()
        with open(primal_project_parameters_file, "r") as file_input:
            primal_settings = Kratos.Parameters(file_input.read())

        self.mdpa_name = primal_settings["solver_settings"]["model_import_settings"]["input_filename"].GetString()

        self.main_model_part_name = primal_settings["solver_settings"]["model_part_name"].GetString()
        self.lift_model_part_name =  self.main_model_part_name + "." + self.lift_model_part_name
        self.drag_model_part_name =  self.main_model_part_name + "." + self.drag_model_part_name

        self.use_automatic_mesh_refinement = self.response_settings["use_automatic_mesh_refinement"].GetBool()
        if (self.use_automatic_mesh_refinement):
            if (automatic_remeshing_process is None):
                raise RuntimeError("Automatic mesh refinement requires to import MeshingApplication. Importing failed with following error msg.\n\t" + automatic_remeshing_error_msg)

            automatic_mesh_refinement_settings = self.response_settings["automatic_mesh_refinement_settings"]
            self.current_objective_moving_average = 0.0
            self.previous_objective_moving_average = 0.0
            self.moving_average_difference_relative_tolerance = automatic_mesh_refinement_settings["moving_average_difference_relative_tolerance"].GetDouble()
            self.moving_average_difference_absolute_tolerance = automatic_mesh_refinement_settings["moving_average_difference_absolute_tolerance"].GetDouble()
            self.moving_average_number_of_iterations = automatic_mesh_refinement_settings["moving_average_number_of_iterations"].GetInt()
            self.read_refinement_data_from = automatic_mesh_refinement_settings["refinement_data_settings"]["refinement_data_from"].GetString()
            self.read_variables = automatic_mesh_refinement_settings["refinement_data_settings"]["list_of_variables"].GetStringArray()
            self.radius_of_influence = automatic_mesh_refinement_settings["radius_of_influence"].GetDouble()
            self.influence_model_part_name = automatic_mesh_refinement_settings["influence_model_part_name"].GetString()
            self.remeshing_process = automatic_remeshing_process(model, automatic_mesh_refinement_settings["mmg_process_settings"])
            self.objective_values = []

    def Initialize(self):
        if (self.use_automatic_mesh_refinement):
            self.remeshing_process.ExecuteInitialize()

    def InitializeIteration(self, optimization_model_part, optimizationIteration):
        if (self.use_automatic_mesh_refinement and len(self.objective_values) == self.moving_average_number_of_iterations):
            Kratos.Logger.PrintInfo(self._GetLabel(), "Initiating automatic mesh refinement")

            self.current_objective_moving_average = 0.0
            for v in self.objective_values:
                self.current_objective_moving_average += v
            self.current_objective_moving_average /= self.moving_average_number_of_iterations

            if (isclose(
                    self.current_objective_moving_average,
                    self.previous_objective_moving_average,
                    self.moving_average_difference_relative_tolerance,
                    self.moving_average_difference_absolute_tolerance)):

                startTime = timer.time()

                if (self.read_refinement_data_from == "lift_adjoint_solution"):
                    self.input_model_part = self.primal_simulation.__GetSolver().GetComputingModelPart()
                elif (self.read_refinement_data_from == "drag_adjoint_solution"):
                    self.input_model_part = self.lift_adjoint_simulation.__GetSolver().GetComputingModelPart()
                elif (self.read_refinement_data_from == "primal_solution"):
                    self.input_model_part = self.drag_adjoint_simulation.__GetSolver().GetComputingModelPart()
                else:
                    raise Exception("Unsupported \"refinement_data_from\" requested. Supported types are \n\t lift_adjoint_solution\n\t drag_adjoint_solution \n\t primal_solution")

                geometric_center = [0.0, 0.0, 0.0]
                for node in self.model[self.influcence_model_part_name].Nodes:
                    geometric_center[0] += node.X
                    geometric_center[1] += node.Y
                    geometric_center[2] += node.Z

                geometric_center[0] /= len(self.influcence_model_part.Nodes)
                geometric_center[1] /= len(self.influcence_model_part.Nodes)
                geometric_center[2] /= len(self.influcence_model_part.Nodes)

                for node in optimization_model_part.Nodes:
                    distance = sqrt((node.X - geometric_center[0]) ** 2 + (node.Y - geometric_center[1]) ** 2 + (node.Z - geometric_center[2]) ** 2)
                    node.Set(Kratos.BLOCKED, distance > self.radius_of_influence)

                for variable_name in self.read_variables:
                    variable = Kratos.KratosGlobals.GetVariable(variable_name)
                    Kratos.VariableUtils().CopyModelPartNodalVarToNonHistoricalVar(variable, variable, self.input_model_part, optimization_model_part)

                self.remeshing_process.ExecuteInitializeSolutionStep()
                self.remeshing_process.ExecuteFinalizeSolutionStep()

                self.objective_values.clear()

                Kratos.Logger.PrintInfo(self._GetLabel(), "Time needed for automatic mesh refinement = ",round(timer.time() - startTime,2),"s")

            self.previous_objective_moving_average = self.current_objective_moving_average

    def UpdateDesign(self, updated_model_part, variable):
        self.updated_model_part = updated_model_part

    def InitializeSolutionStep(self):
        self.lift = None
        self.drag = None

        startTime = timer.time()

        # copy data to the working directory
        LiftToDrag._RecursiveCopy(str(self.problem_setup_folder), ".")

        # write the new shape mdpa
        Kratos.ModelPartIO(self.mdpa_name, Kratos.IO.WRITE | Kratos.IO.MESH_ONLY).WriteModelPart(self.updated_model_part)

        Kratos.Logger.PrintInfo(self._GetLabel(), "Time needed to copy data = ",round(timer.time() - startTime,2),"s")

    def CalculateValue(self):
        startTime = timer.time()

        # open primal settings
        primal_settings_file_name = self.problem_setup_file_settings["primal_project_parameters_file"].GetString()
        primal_filled_settings_file = Path(primal_settings_file_name[:primal_settings_file_name.rfind(".")] + "_final.json")

        if (not primal_filled_settings_file.is_file()):
            with open(primal_settings_file_name, "r") as file_input:
                primal_parameters = Kratos.Parameters(file_input.read())

            # add response function outputs
            LiftToDrag._AddResponseFunctionOutput(primal_parameters, self.lift_model_part_name)
            LiftToDrag._AddResponseFunctionOutput(primal_parameters, self.drag_model_part_name)

            # set the loggers
            file_logger = Kratos.FileLoggerOutput("primal_evaluation.log")
            default_severity = Kratos.Logger.GetDefaultOutput().GetSeverity()
            Kratos.Logger.GetDefaultOutput().SetSeverity(Kratos.Logger.Severity.WARNING)
            Kratos.Logger.AddOutput(file_logger)

            # run the primal analysis
            self.primal_model = Kratos.Model()
            self.primal_simulation = RANSAnalysis(self.primal_model, primal_parameters)
            self.primal_simulation.Run()

            with open(str(primal_filled_settings_file), "w") as file_output:
                file_output.write(primal_parameters.PrettyPrintJsonString())

            # flush the primal output
            Kratos.Logger.Flush()
            Kratos.Logger.GetDefaultOutput().SetSeverity(default_severity)
            Kratos.Logger.RemoveOutput(file_logger)
        else:
            # following is done to initialize default settings
            with open(str(primal_filled_settings_file), "r") as file_input:
                primal_parameters = Kratos.Parameters(file_input.read())
            Kratos.Logger.PrintInfo(self._GetLabel(), "Found existing completed primal evaluation.")

        # calculate lift and drag
        self.lift = LiftToDrag._CalculateTimeAveragedDrag(primal_parameters, self.lift_model_part_name, self.lift_direction)
        self.drag = LiftToDrag._CalculateTimeAveragedDrag(primal_parameters, self.drag_model_part_name, self.drag_direction)

        if (self.use_automatic_mesh_refinement):
            self.objective_values.append(self.GetValue())
            if (len(self.objective_values) > self.moving_average_number_of_iterations):
                del self.objective_values[0]

        Kratos.Logger.PrintInfo(self._GetLabel(), "Time needed for calculating the response value = ",round(timer.time() - startTime,2),"s")

    def CalculateGradient(self):
        # solve adjoint lift problem
        lift_start_time = timer.time()

        self.lift_adjoint_model = Kratos.Model()
        self.lift_adjoint_simulation = LiftToDrag._SolveAdjointProblem(
                self.lift_adjoint_model,
                self.problem_setup_file_settings["lift_adjoint_project_parameters_file"].GetString(),
                "lift_adjoint_evaluation.log")

        Kratos.Logger.PrintInfo(self._GetLabel(), "Time needed for solving the lift adjoint analysis = ",round(timer.time() - lift_start_time,2),"s")

        # solve adjoint drag problem
        drag_start_time = timer.time()

        self.drag_adjoint_model = Kratos.Model()
        self.drag_adjoint_simulation = LiftToDrag._SolveAdjointProblem(
                self.drag_adjoint_model,
                self.problem_setup_file_settings["drag_adjoint_project_parameters_file"].GetString(),
                "drag_adjoint_evaluation.log")

        Kratos.Logger.PrintInfo(self._GetLabel(), "Time needed for solving the drag adjoint analysis = ",round(timer.time() - drag_start_time,2),"s")
        Kratos.Logger.PrintInfo(self._GetLabel(), "Time needed for solving the total adjoint analysis = ",round(timer.time() - lift_start_time,2),"s")

    def GetValue(self):
        if (self.lift is None or self.drag is None):
            raise RuntimeError("GetValue: Lift or Drag is not computed. Please use CalculateValue method first.")
        return self.lift / self.drag

    def GetNodalGradient(self, variable):
        if (self.lift is None or self.drag is None):
            raise RuntimeError("GetNodalGradient: Lift or Drag is not computed. Please use CalculateValue method first.")
        if variable != Kratos.SHAPE_SENSITIVITY:
            raise RuntimeError("GetNodalGradient: No gradient for {}!".format(variable.Name))

        gradient = {}
        lift_adjoint_model_part = self.lift_adjoint_model[self.main_model_part_name]
        drag_adjoint_model_part = self.drag_adjoint_model[self.main_model_part_name]
        for lift_node, drag_node in zip(lift_adjoint_model_part.Nodes, drag_adjoint_model_part.Nodes):
            node_id = lift_node.Id
            if (lift_node.Id != drag_node.Id):
                raise RuntimeError("Node mismatch")
            gradient[node_id] = \
                lift_node.GetSolutionStepValue(variable) / self.drag - \
                self.lift * drag_node.GetSolutionStepValue(variable) / (self.drag ** 2)
        return gradient

    def IsEvaluatedInFolder(self):
        return True

    @staticmethod
    def _GetLabel():
        return "AdjointLiftToDrag"

    @staticmethod
    def _GetResponseFunctionOutputProcess(kratos_parameters, model_part_name):
        auxiliar_process_list = kratos_parameters["processes"]["auxiliar_process_list"]
        for process_settings in auxiliar_process_list:
            if (
                process_settings.Has("python_module") and process_settings["python_module"].GetString() == "compute_body_fitted_drag_process" and
                process_settings.Has("kratos_module") and process_settings["kratos_module"].GetString() == "KratosMultiphysics.FluidDynamicsApplication" and
                process_settings["Parameters"].Has("model_part_name") and process_settings["Parameters"]["model_part_name"].GetString() == model_part_name
                ):
                return process_settings

        return None

    @staticmethod
    def _SolveAdjointProblem(model, adjoint_parameters_file_name, log_file_name):
        with open(adjoint_parameters_file_name, "r") as file_input:
            adjoint_parameters = Kratos.Parameters(file_input.read())

        # set the lift loggers
        file_logger = Kratos.FileLoggerOutput(log_file_name)
        default_severity = Kratos.Logger.GetDefaultOutput().GetSeverity()
        Kratos.Logger.GetDefaultOutput().SetSeverity(Kratos.Logger.Severity.WARNING)
        Kratos.Logger.AddOutput(file_logger)

        adjoint_simulation = AdjointRANSAnalysis(model, adjoint_parameters)
        adjoint_simulation.Run()

        Kratos.Logger.Flush()

        Kratos.Logger.GetDefaultOutput().SetSeverity(default_severity)
        Kratos.Logger.RemoveOutput(file_logger)

        return adjoint_simulation

    @staticmethod
    def _CalculateTimeAveragedDrag(kratos_parameters, model_part_name, direction):
        output_process = LiftToDrag._GetResponseFunctionOutputProcess(kratos_parameters, model_part_name)
        if (output_process is not None):
            output_file_name = output_process["Parameters"]["output_file_settings"]["file_name"].GetString()
            time_steps, reactions = LiftToDrag._ReadDrag(output_file_name)
            total_drag = 0.0
            for reaction in reversed(reactions):
                total_drag += reaction[0] * direction[0] + reaction[1] * direction[1] + reaction[2] * direction[2]
                if (kratos_parameters["solver_settings"]["time_scheme_settings"]["scheme_type"].GetString() == "steady"):
                    break
            if len(time_steps) > 1:
                delta_time = time_steps[1] - time_steps[0]
                total_drag *= delta_time
            return total_drag
        else:
            raise RuntimeError("No \"compute_body_fitted_drag_process\" found in auxiliar_process_list.")

    @staticmethod
    def _ReadDrag(file_name):
        with open(file_name, "r") as file_input:
            lines = file_input.readlines()
        time_steps = []
        reaction = []
        for line in lines:
            line = line.strip()
            if len(line) == 0 or line[0] == '#':
                continue
            time_step_data = [float(v) for v in line.split()]
            time, fx, fy, fz = time_step_data
            time_steps.append(time)
            reaction.append([fx, fy, fz])
        return time_steps, reaction

    @staticmethod
    def _AddResponseFunctionOutput(kratos_parameters, model_part_name):
        response_parameters = Kratos.Parameters(R'''
        {
            "python_module": "compute_body_fitted_drag_process",
            "kratos_module": "KratosMultiphysics.FluidDynamicsApplication",
            "process_name": "ComputeBodyFittedDragProcess",
            "Parameters": {
                "model_part_name": "<model_part_name>",
                "interval": [
                    0.0,
                    1e30
                ],
                "write_drag_output_file": true,
                "print_drag_to_screen": false,
                "print_format": "22.15e"
            }
        }
        ''')

        if (LiftToDrag._GetResponseFunctionOutputProcess(kratos_parameters, model_part_name) is None):
            response_parameters["Parameters"]["model_part_name"].SetString(model_part_name)
            kratos_parameters["processes"]["auxiliar_process_list"].Append(response_parameters)

    @staticmethod
    def _RecursiveCopy(src, dest):
        src_path = Path(src)
        dest_path = Path(dest)
        for item in src_path.iterdir():
            if (item.is_file()):
                shutil.copy(str(item), dest)
            elif (item.is_dir()):
                new_path = dest_path / item.relative_to(src)
                new_path.mkdir(exist_ok=True)
                LiftToDrag._RecursiveCopy(str(item), str(new_path))
