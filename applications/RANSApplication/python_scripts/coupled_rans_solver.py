from pathlib import Path

# Importing the Kratos Library
import KratosMultiphysics as Kratos
from KratosMultiphysics import IsDistributedRun
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable
from KratosMultiphysics.python_solver import PythonSolver

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.RANSApplication as KratosRANS

# Import application specific modules
from KratosMultiphysics.RANSApplication.formulations import Factory as FormulationFactory
from KratosMultiphysics.RANSApplication.formulations.utilities import InitializeWallLawProperties
from KratosMultiphysics.RANSApplication.formulations.utilities import GetTimeDerivativeVariablesRecursively
from KratosMultiphysics.RANSApplication.formulations.utilities import AddFileLoggerOutput
from KratosMultiphysics.RANSApplication.formulations.utilities import RemoveFileLoggerOutput
from KratosMultiphysics.RANSApplication import RansVariableUtilities
from KratosMultiphysics.FluidDynamicsApplication.check_and_prepare_model_process_fluid import CheckAndPrepareModelProcess

# case specific imports
if (IsDistributedRun() and CheckIfApplicationsAvailable("TrilinosApplication")):
    import KratosMultiphysics.TrilinosApplication as KratosTrilinos     # MPI solvers
    from KratosMultiphysics.mpi.distributed_import_model_part_utility import DistributedImportModelPartUtility
elif (IsDistributedRun()):
    raise Exception("Distributed run requires TrilinosApplication")

class CoupledRANSSolver(PythonSolver):
    def __init__(self, model, custom_settings):
        """RANS solver to be used with RANSFormulations

        This solver creates PythonSolver based on the RANSFormulations specified in custom_settings.

        Args:
            model (Kratos.Model): Model to be used in the solver.
            custom_settings (Kratos.Parameters): Settings to be used in the solver.
        """

        self._validate_settings_in_baseclass = True  # To be removed eventually
        super().__init__(model, custom_settings)

        model_part_name = self.settings["model_part_name"].GetString()
        if model_part_name == "":
            raise Exception(
                'Please provide the model part name as the "model_part_name" (string) parameter!'
            )

        if self.model.HasModelPart(model_part_name):
            self.main_model_part = self.model.GetModelPart(model_part_name)
        else:
            self.main_model_part = self.model.CreateModelPart(model_part_name)

        self.domain_size = self.settings["domain_size"].GetInt()
        if self.domain_size == -1:
            raise Exception(
                'Please provide the domain size as the "domain_size" (int) parameter!'
            )
        self.main_model_part.ProcessInfo.SetValue(Kratos.DOMAIN_SIZE,
                                                  self.domain_size)

        self.formulation = FormulationFactory(self.main_model_part,
                                self.settings["formulation_settings"])

        self.formulation.SetConstants(self.settings["constants"])
        self.formulation.SetIsPeriodic(self.settings["consider_periodic_conditions"].GetBool())

        self.is_periodic = self.formulation.IsPeriodic()

        self.formulation.SetTimeSchemeSettings(self.settings["time_scheme_settings"])
        self.formulation.SetWallFunctionSettings(self.settings["wall_function_settings"])
        scheme_type = self.settings["time_scheme_settings"]["scheme_type"].GetString()
        if (scheme_type == "steady"):
            self.is_steady = True
        else:
            self.is_steady = False

        self.main_model_part.ProcessInfo[KratosRANS.RANS_IS_STEADY] = self.is_steady

        self.is_converged = False
        self.min_buffer_size = self.formulation.GetMinimumBufferSize()
        self.move_mesh = self.settings["move_mesh"].GetBool()
        self.echo_level = self.settings["echo_level"].GetInt()

        if (self.settings["time_scheme_settings"].Has("ramp_up_interval")):
            self.ramp_up_interval = self.settings["time_scheme_settings"]["ramp_up_interval"].GetVector()
            self.ramp_up_interval = sorted(self.ramp_up_interval)
            if (len(self.ramp_up_interval) != 2):
                raise Exception("Ramp up interval should only have two values indicating ramp up start time and end time.")
        else:
            self.ramp_up_interval = None

        self.automatic_mesh_refinement_step = 0
        self.perform_automatic_mesh_refinement = False

        Kratos.Logger.PrintInfo(self.__class__.__name__,
                                            "Solver construction finished.")

    @classmethod
    def GetDefaultParameters(cls):
        ##settings string in json format
        default_settings = Kratos.Parameters("""
        {
            "solver_type": "CoupledRANS",
            "model_part_name": "FluidModelPart",
            "domain_size": -1,
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name",
                "reorder": false
            },
            "material_import_settings": {
                "materials_filename": ""
            },
            "consider_periodic_conditions": false,
            "formulation_settings": {},
            "wall_function_settings": {},
            "echo_level": 0,
            "volume_model_part_name": "volume_model_part",
            "skin_parts"   : [""],
            "no_skin_parts": [""],
            "assign_neighbour_elements_to_conditions": true,
            "move_mesh": false,
            "time_scheme_settings":{
                "scheme_type": "steady",
                "ramp_up_interval": [-1.0, -1.0]
            },
            "time_stepping": {
                "automatic_time_step" : false,
                "CFL_number"          : 1,
                "minimum_delta_time"  : 1e-4,
                "maximum_delta_time"  : 0.01,
                "time_step"           : 0.0
            },
            "constants": {},
            "adaptive_mesh_refinement_based_on_response_function": false,
            "adaptive_mesh_refinement_based_on_response_function_settings": {
                "primal_problem_project_parameters_file_name" : "PLEASE_SPECIFY_PRIMAL_PROBLEM_PROJECT_PARAMETERS_FILE_NAME",
                "adjoint_problem_project_parameters_file_name": "PLEASE_SPECIFY_ADJOINT_PROBLEM_PROJECT_PARAMETERS_FILE_NAME",
                "time_range"                                  : [0.0, 1e+30],
                "step_interval"                               : 10,
                "output_model_part_name"                      : "adapted_mesh",
                "output_log_files_folder_name"                : "logs",
                "mmg_process_parameters": {
                    "strategy"               : "hessian",
                    "automatic_remesh"       : false,
                    "enforce_current"        : false,
                    "maximal_size"           : 0.5,
                    "minimal_size"           : 1e-3,
                    "hessian_strategy_parameters":{
                        "metric_variable": ["VELOCITY_X", "VELOCITY_Y", "PRESSURE"],
                        "non_historical_metric_variable": [false, false, false],
                        "interpolation_error": 1e-5,
                        "use_response_function_interpolation_error": true,
                        "response_function_interpolation_variable_index": [0, 1, 2]
                    },
                    "model_part_name" : "FluidModelPart",
                    "step_frequency"  : 1,
                    "force_min"       : true,
                    "force_max"       : true,
                    "echo_level"      : 3
                }
            }
        }""")

        default_settings.AddMissingParameters(super().GetDefaultParameters())
        return default_settings

    def AddVariables(self):
        self.formulation.AddVariables()

        if self.is_periodic:
            self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.PATCH_INDEX)

        if (IsDistributedRun()):
            self.main_model_part.AddNodalSolutionStepVariable(Kratos.PARTITION_INDEX)

        if (self.formulation.ElementHasNodalProperties()):
            self.main_model_part.AddNodalSolutionStepVariable(Kratos.DENSITY)
            self.main_model_part.AddNodalSolutionStepVariable(Kratos.VISCOSITY)

        Kratos.Logger.PrintInfo(
            self.__class__.__name__, "Solver variables added correctly.")

    def AddDofs(self):
        self.formulation.AddDofs()

        Kratos.Logger.PrintInfo(self.__class__.__name__,
                                            "Solver dofs added correctly.")

    def ImportModelPart(self):
        if (IsDistributedRun()):
            ## Construct the MPI import model part utility
            self.distributed_model_part_importer = DistributedImportModelPartUtility(self.main_model_part, self.settings)
            ## Execute the Metis partitioning and reading
            self.distributed_model_part_importer.ImportModelPart()
        else:
            # we can use the default implementation in the base class
            self._ImportModelPart(self.main_model_part,
                                  self.settings["model_import_settings"])

    def PrepareModelPart(self):
        if not self.main_model_part.ProcessInfo[
                Kratos.IS_RESTARTED]:
            ## Set fluid properties from materials json file
            materials_imported = self._SetPhysicalProperties()
            if not materials_imported:
                Kratos.Logger.PrintWarning(
                    self.__class__.__name__,
                    "Material properties have not been imported. Check \'material_import_settings\' in your ProjectParameters.json."
                )
            ## remove fluid_computational_model_part if it exists. It will be there if adaptive mesh refinement is used.
            if (self.main_model_part.HasSubModelPart("fluid_computational_model_part")):
                self.main_model_part.RemoveSubModelPart("fluid_computational_model_part")
            ## Executes the check and prepare model process
            self._ExecuteCheckAndPrepare()
            ## Set buffer size
            self.main_model_part.SetBufferSize(self.min_buffer_size)

        if (IsDistributedRun()):
            self.distributed_model_part_importer.CreateCommunicators()

        self.formulation.PrepareModelPart()

        Kratos.Logger.PrintInfo(self.__class__.__name__,
                                            "Model reading finished.")

    def ExportModelPart(self):
        ## Model part writing
        name_out_file = self.settings["model_import_settings"][
            "input_filename"].GetString() + ".out"
        Kratos.ModelPartIO(
            name_out_file,
            Kratos.IO.WRITE).WriteModelPart(self.main_model_part)

        Kratos.Logger.PrintInfo(self.__class__.__name__,
                                            "Model export finished.")

    def GetMinimumBufferSize(self):
        return self.min_buffer_size

    def Initialize(self):
        if (IsDistributedRun()):
            self.EpetraComm = KratosTrilinos.CreateCommunicator()
            self.formulation.SetCommunicator(self.EpetraComm)
        else:
            self.formulation.SetCommunicator(None)

        self.main_model_part.ProcessInfo[Kratos.STEP] = 0

        # If needed, create the estimate time step utility
        if (self.settings["time_stepping"]["automatic_time_step"].GetBool()):
            self.EstimateDeltaTimeUtility = self._GetAutomaticTimeSteppingUtility()

        RansVariableUtilities.AssignBoundaryFlagsToGeometries(self.main_model_part)
        self.formulation.Initialize()

        if (self.settings["adaptive_mesh_refinement_based_on_response_function"].GetBool()):
            number_of_solving_variables = len(self.formulation.GetSolvingVariables())
            Kratos.VariableUtils().SetNonHistoricalVariable(Kratos.RESPONSE_FUNCTION_INTERPOLATION_ERROR, Kratos.Vector(number_of_solving_variables, 0.0), self.main_model_part.Nodes)

        Kratos.Logger.PrintInfo(self.__class__.__name__, self.formulation.GetInfo())

        Kratos.Logger.PrintInfo(self.__class__.__name__,
                                            "Solver initialization finished.")

    def AdvanceInTime(self, current_time):
        dt = self._ComputeDeltaTime()
        new_time = current_time + dt

        self.main_model_part.CloneTimeStep(new_time)
        self.main_model_part.ProcessInfo[Kratos.STEP] += 1

        return new_time

    def InitializeSolutionStep(self):
        if (self.perform_automatic_mesh_refinement):
            # perform adaptive remeshing
            Kratos.Logger.PrintInfo(self.__class__.__name__, "Remeshing based on response function error analysis...")
            from KratosMultiphysics.MeshingApplication.mmg_process import MmgProcess
            remeshing_process = MmgProcess(self.main_model_part.GetModel(), self.settings["adaptive_mesh_refinement_based_on_response_function_settings"]["mmg_process_parameters"])
            remeshing_process.ExecuteInitialize()
            remeshing_process.ExecuteInitializeSolutionStep()
            remeshing_process.ExecuteFinalizeSolutionStep()
            remeshing_process.ExecuteFinalize()

            self.perform_automatic_mesh_refinement = False
            Kratos.Logger.PrintInfo(self.__class__.__name__, "Finished adaptive mesh refinement.")

        self.formulation.InitializeSolutionStep()

    def Predict(self):
        self.formulation.Predict()

    def SolveSolutionStep(self):
        self.formulation.SolveCouplingStep()
        self.is_converged = self.formulation.IsConverged()

        if not self.is_converged and not self.IsSteadySimulation() and self.echo_level > -1:
            msg = "Fluid solver did not converge for step " + str(self.main_model_part.ProcessInfo[Kratos.STEP]) + "\n"
            msg += "corresponding to time " + str(self.main_model_part.ProcessInfo[Kratos.TIME]) + "\n"
            Kratos.Logger.PrintWarning(self.__class__.__name__, msg)
        return self.is_converged

    def FinalizeSolutionStep(self):
        self.formulation.FinalizeSolutionStep()

        if (self.settings["adaptive_mesh_refinement_based_on_response_function"].GetBool()):
            adaptive_mesh_refinement_settings = self.settings["adaptive_mesh_refinement_based_on_response_function_settings"]
            automatic_mesh_refinement_time_range = adaptive_mesh_refinement_settings["time_range"].GetVector()

            current_time = self.main_model_part.ProcessInfo[Kratos.TIME]
            if (current_time >= automatic_mesh_refinement_time_range[0] and current_time <= automatic_mesh_refinement_time_range[1]):
                if (self.automatic_mesh_refinement_step % adaptive_mesh_refinement_settings["step_interval"].GetInt() == 0):
                    self._ExecuteAdaptiveMeshRefinementBasedOnResponseFunction(current_time)
                self.automatic_mesh_refinement_step += 1

    def Check(self):
        self.formulation.Check()

    def Clear(self):
        self.formulation.Clear()

    def IsSteadySimulation(self):
        return self.is_steady

    def IsConverged(self):
        if (self.is_steady):
            current_time = self.main_model_part.ProcessInfo[Kratos.TIME]
            if (self.ramp_up_interval is not None):
                if (current_time >= self.ramp_up_interval[0] and current_time <= self.ramp_up_interval[1]):
                    if (self.is_converged):
                        Kratos.Logger.PrintInfo(self.__class__.__name__, "Continuing steady simulation because it is still within ramp up interval.")
                    return False
                else:
                    return self.is_converged
            else:
                return self.is_converged
        else:
            return False

    def GetComputingModelPart(self):
        if not self.main_model_part.HasSubModelPart(
                "fluid_computational_model_part"):
            raise Exception("The ComputingModelPart was not created yet!")
        return self.main_model_part.GetSubModelPart(
            "fluid_computational_model_part")

    def _ComputeDeltaTime(self):
        # Automatic time step computation according to user defined CFL number
        if (self.settings["time_stepping"]["automatic_time_step"].GetBool()):
            delta_time = self.EstimateDeltaTimeUtility.EstimateDt()
        # User-defined delta time
        else:
            delta_time = self.settings["time_stepping"]["time_step"].GetDouble(
            )

        return delta_time

    def _GetAutomaticTimeSteppingUtility(self):
        estimate_delta_time_utility = KratosCFD.EstimateDtUtility(
            self.GetComputingModelPart(),
            self.settings["time_stepping"])

        return estimate_delta_time_utility

    def _ExecuteCheckAndPrepare(self):
        ## Check that the input read has the shape we like
        prepare_model_part_settings = Kratos.Parameters("{}")
        prepare_model_part_settings.AddValue(
            "volume_model_part_name", self.settings["volume_model_part_name"])
        prepare_model_part_settings.AddValue("skin_parts",
                                             self.settings["skin_parts"])
        if (self.settings.Has("assign_neighbour_elements_to_conditions")):
            prepare_model_part_settings.AddValue(
                "assign_neighbour_elements_to_conditions",
                self.settings["assign_neighbour_elements_to_conditions"])
        else:
            warn_msg = "\"assign_neighbour_elements_to_conditions\" should be added to defaults of " + self.__class__.__name__
            Kratos.Logger.PrintWarning(
                '\n\x1b[1;31mDEPRECATION-WARNING\x1b[0m', warn_msg)

        CheckAndPrepareModelProcess(self.main_model_part,
                                    prepare_model_part_settings).Execute()

    def _SetPhysicalProperties(self):
        # Check if the fluid properties are provided using a .json file
        materials_filename = self.settings["material_import_settings"][
            "materials_filename"].GetString()
        if (materials_filename != ""):
            # Add constitutive laws and material properties from json file to model parts.
            material_settings = Kratos.Parameters(
                """{"Parameters": {"materials_filename": ""}} """)
            material_settings["Parameters"]["materials_filename"].SetString(
                materials_filename)
            Kratos.ReadMaterialsUtility(material_settings,
                                                    self.model)
            # add wall law properties
            InitializeWallLawProperties(self.model)

            # initialize constitutive laws
            RansVariableUtilities.SetElementConstitutiveLaws(self.main_model_part.Elements)

            materials_imported = True
        else:
            materials_imported = False

        # If the element uses nodal material properties, transfer them to the nodes
        if self.formulation.ElementHasNodalProperties():
            self._SetNodalProperties()

        return materials_imported

    def _SetNodalProperties(self):
        # Get density and dynamic viscostity from the properties of the first element
        for el in self.main_model_part.Elements:
            rho = el.Properties.GetValue(Kratos.DENSITY)
            mu = el.Properties.GetValue(Kratos.DYNAMIC_VISCOSITY)

            if (rho <= 0.0):
                raise Exception ("DENSITY is not properly set in material properties.")

            nu = mu / rho

            Kratos.VariableUtils().SetVariable(Kratos.DENSITY, rho, self.main_model_part.Nodes)
            Kratos.VariableUtils().SetVariable(Kratos.VISCOSITY, nu, self.main_model_part.Nodes)

            break

    def _ExecuteAdaptiveMeshRefinementBasedOnResponseFunction(self, current_time):
        if (not self.is_steady):
            Kratos.Logger.PrintInfo(self.__class__.__name__, "Performing adaptive mesh refinement...")

            adaptive_mesh_refinement_based_on_response_function_settings = self.settings["adaptive_mesh_refinement_based_on_response_function_settings"]

            log_path = Path(adaptive_mesh_refinement_based_on_response_function_settings["output_log_files_folder_name"].GetString())
            log_path.mkdir(parents=True, exist_ok=True)
            current_step = self.main_model_part.ProcessInfo[Kratos.STEP]
            primal_log_file_name = str(log_path / "primal_error_analysis_evaluation_step_{:d}.log".format(current_step))
            adjoint_log_file_name = str(log_path / "adjoint_error_analysis_evaluation_step_{:d}.log".format(current_step))

            # write the current mesh so primal and adjoint problems can read it properly
            Kratos.ModelPartIO(adaptive_mesh_refinement_based_on_response_function_settings["output_model_part_name"].GetString(), Kratos.IO.WRITE | Kratos.IO.MESH_ONLY).WriteModelPart(self.main_model_part)

            # set the starting point of the simulation to current step
            dt = self._ComputeDeltaTime()
            time_steps = adaptive_mesh_refinement_based_on_response_function_settings["step_interval"].GetInt()
            start_time = current_time
            end_time = current_time + dt * time_steps

            # compute the list of variables to be passed as initial conditions for primal
            # and adjoint problem
            list_of_solving_variables = []
            for var in self.formulation.GetSolvingVariables():
                time_derivative_vars_list = GetTimeDerivativeVariablesRecursively(var)
                for time_derivative_var in time_derivative_vars_list:
                    list_of_solving_variables.append(time_derivative_var)

            total_variables_list = []
            for var in list_of_solving_variables:
                total_variables_list.append((var.Name(), True, var.Name(), True))
            total_variables_list.append(("RELAXED_ACCELERATION", False, "RELAXED_ACCELERATION", False))

            # open primal parameters
            with open(adaptive_mesh_refinement_based_on_response_function_settings["primal_problem_project_parameters_file_name"].GetString(),'r') as parameter_file:
                primal_parameters = Kratos.Parameters(parameter_file.read())

            primal_parameters["problem_data"]["start_time"].SetDouble(start_time)
            primal_parameters["problem_data"]["end_time"].SetDouble(end_time)

            # solve the primal problem
            Kratos.Logger.PrintInfo(self.__class__.__name__, "Solving primal problem for error analysis...")
            default_severity, file_logger = AddFileLoggerOutput(primal_log_file_name)
            from KratosMultiphysics.RANSApplication.rans_analysis import RANSAnalysis
            primal_model = Kratos.Model()
            primal_simulation = RANSAnalysis(primal_model, primal_parameters)
            primal_simulation.Initialize()

            # copy initialization data to primal solver
            KratosRANS.RansVariableDataTransferProcess(
                self.main_model_part.GetModel(),
                primal_model,
                self.main_model_part.FullName(),
                primal_simulation._GetSolver().GetComputingModelPart().FullName(),
                ["execute"],
                total_variables_list,
                self.echo_level).Execute()

            primal_simulation.RunSolutionLoop()
            primal_simulation.Finalize()
            RemoveFileLoggerOutput(default_severity, file_logger)

            # open adjoint problem
            with open(adaptive_mesh_refinement_based_on_response_function_settings["adjoint_problem_project_parameters_file_name"].GetString(),'r') as parameter_file:
                adjoint_parameters = Kratos.Parameters(parameter_file.read())

            adjoint_parameters["problem_data"]["start_time"].SetDouble(start_time + dt)
            adjoint_parameters["problem_data"]["end_time"].SetDouble(end_time)

            # run adjoint problem
            # importing is done here to avoid circular dependency
            Kratos.Logger.PrintInfo(self.__class__.__name__, "Solving adjoint problem for error analysis...")
            default_severity, file_logger = AddFileLoggerOutput(adjoint_log_file_name)
            from KratosMultiphysics.RANSApplication.adjoint_rans_analysis import AdjointRANSAnalysis
            adjoint_model = Kratos.Model()
            adjoint_simulation = AdjointRANSAnalysis(adjoint_model, adjoint_parameters)
            adjoint_simulation.Run()
            RemoveFileLoggerOutput(default_severity, file_logger)

            # copy interpolation error data from adjoint model.
            KratosRANS.RansVariableDataTransferProcess(
                adjoint_model,
                self.main_model_part.GetModel(),
                adjoint_simulation._GetSolver().GetComputingModelPart().FullName(),
                self.main_model_part.FullName(),
                ["execute"],
                [("RESPONSE_FUNCTION_INTERPOLATION_ERROR", False, "RESPONSE_FUNCTION_INTERPOLATION_ERROR", False)],
                self.echo_level).Execute()

            self.perform_automatic_mesh_refinement = True
            Kratos.Logger.PrintInfo(self.__class__.__name__, "Computed response based interpolation errors for adaptive mesh refinement.")

    def Finalize(self):
        self.formulation.Finalize()
