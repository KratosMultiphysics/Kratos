from pathlib import Path
from shutil import rmtree

# Importing the Kratos Library
import KratosMultiphysics as Kratos
from KratosMultiphysics import IsDistributedRun
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable
from KratosMultiphysics.kratos_utilities import IssueDeprecationWarning
from KratosMultiphysics.python_solver import PythonSolver
from KratosMultiphysics.process_factory import KratosProcessFactory

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.RANSApplication as KratosRANS
import KratosMultiphysics.MappingApplication as KratosMapping

# Import application specific modules
from KratosMultiphysics.RANSApplication.formulations import Factory as FormulationFactory
from KratosMultiphysics.RANSApplication.formulations.utilities import InitializeWallLawProperties
from KratosMultiphysics.RANSApplication.formulations.utilities import ExecutionScope
from KratosMultiphysics.RANSApplication.formulations.utilities import SolveProblem
from KratosMultiphysics.RANSApplication.formulations.utilities import GetConvergenceInfo
from KratosMultiphysics.RANSApplication import RansVariableUtilities
from KratosMultiphysics.FluidDynamicsApplication.check_and_prepare_model_process_fluid import CheckAndPrepareModelProcess

from KratosMultiphysics.RANSApplication.post_processing_utilities import GetHDF5File
from KratosMultiphysics.RANSApplication.post_processing_utilities import OutputNodalResultsToHDF5
from KratosMultiphysics.RANSApplication.post_processing_utilities import InputNodalResultsFromHDF5

# case specific imports
if (IsDistributedRun() and CheckIfApplicationsAvailable("TrilinosApplication")):
    import KratosMultiphysics.TrilinosApplication as KratosTrilinos     # MPI solvers
    from KratosMultiphysics.mpi.distributed_import_model_part_utility import DistributedImportModelPartUtility
elif (IsDistributedRun()):
    raise Exception("Distributed run requires TrilinosApplication")

class MeshInterpolationMethod:
    def __init__(self, model_part, settings):
        self.model_part = model_part
        self.model = self.model_part.GetModel()
        self.settings = settings
        self.interpolation_settings = settings["mesh_interpolation_settings"]

    def Initialize(self):
        pass

    def InitializeSolutionStep(self):
        pass

    def Interpolate(self, delta_time):
        pass

    def FinalizeSolutionStep(self):
        pass


class MappedMeshInterpolationMethod(MeshInterpolationMethod):
    def __init__(self, model_part, settings):
        super().__init__(model_part, settings)

        default_settings = Kratos.Parameters("""
        {
            "interpolation_method": "mapped",
            "mapper_settings": {
                "mapper_type": "nearest_element",
                "echo_level" : 1,
                "search_radius": 0.1
            },
            "post_processing_processes": []
        }""")

        self.interpolation_settings.ValidateAndAssignDefaults(default_settings)
        self.list_of_post_processing_processes = []

        Kratos.Logger.PrintInfo(self.__class__.__name__, "Created mapped mesh interpolation for {:s}".format(self.model_part.FullName()))

    def Initialize(self):
        # create list of post processing processes
        factory = KratosProcessFactory(self.model)
        self.post_processing_processes = factory.ConstructListOfProcesses(self.interpolation_settings["post_processing_processes"])
        for process in self.post_processing_processes:
            self.list_of_post_processing_processes.append(process)

    def Interpolate(self, delta_time):
        destination_model_part = self.model_part
        model = destination_model_part.GetModel()
        source_model_part = model[destination_model_part.FullName() + "_Old"]

        # remove conditions from the source model_part because it cannot have
        # elements and conditions both for mapping
        Kratos.VariableUtils().SetFlag(Kratos.TO_ERASE, True, source_model_part.Conditions)
        source_model_part.RemoveConditionsFromAllLevels(Kratos.TO_ERASE)

        # create mapper
        mapper = KratosMapping.MapperFactory.CreateMapper(source_model_part, destination_model_part, self.interpolation_settings["mapper_settings"].Clone())

        solution_step_variable_names_list = KratosRANS.RansVariableUtilities.GetSolutionstepVariableNamesList(destination_model_part)
        variable_transfer_list = []

        # map current solution step data
        for solution_step_variable_name in solution_step_variable_names_list:
            solution_step_variable = Kratos.KratosGlobals.GetVariable(solution_step_variable_name)
            mapper.Map(solution_step_variable, solution_step_variable)
            variable_transfer_list.append([
                solution_step_variable_name,
                True,
                0,
                solution_step_variable_name,
                True,
                1
            ])


        # copy current mapped solution step data to previous solution step
        # since at this point both of the solution steps are cloned
        # It is just the start of the the solve step.
        KratosRANS.RansVariableDataTransferProcess(
            model,
            destination_model_part.FullName(),
            destination_model_part.FullName(),
            ["execute"],
            variable_transfer_list,
            2
        ).Execute()

        # delete the old model part and mapper model part
        model.DeleteModelPart(source_model_part.FullName())

        # now correct boundary values
        for process in self.list_of_post_processing_processes:
            process.Execute()

class SmoothMappedMeshInterpolationMethod(MeshInterpolationMethod):
    def __init__(self, model_part, settings):
        super().__init__(model_part, settings)

        default_settings = Kratos.Parameters("""
        {
            "interpolation_method": "smooth_mapped",
            "primal_transient_smoothing_project_parameters_file_name": "PLEASE_SPECIFY_PROJECT_PARAMETERS_FILE_NAME",
            "number_of_smoothing_step": 10,
            "echo_level": 0,
            "initialization_hdf5_file": "initialization.h5",
            "clear_intermediate_files": true,
            "mapper_settings": {
                "mapper_type": "nearest_element",
                "echo_level" : 1,
                "search_radius": 0.1
            },
            "post_processing_processes": []
        }""")

        self.interpolation_settings.ValidateAndAssignDefaults(default_settings)
        self.step_data_list = []
        self.number_of_smoothing_steps = self.interpolation_settings["number_of_smoothing_step"].GetInt()
        self.output_path = Path(self.settings["output_path"].GetString()) / "smooth_mapped_interpolation"
        self.list_of_post_processing_processes = []
        self.echo_level = self.interpolation_settings["echo_level"].GetInt()

        primal_transient_smoothing_project_parameters_file_name = self.interpolation_settings["primal_transient_smoothing_project_parameters_file_name"].GetString()
        with open(primal_transient_smoothing_project_parameters_file_name, "r") as file_input:
            self.primal_parameters = Kratos.Parameters(file_input.read())

        Kratos.Logger.PrintInfo(self.__class__.__name__, "Created smoothed mapped mesh interpolation for {:s}".format(self.model_part.FullName()))

    def Initialize(self):
        # create list of post processing processes
        factory = KratosProcessFactory(self.model)
        self.post_processing_processes = factory.ConstructListOfProcesses(self.interpolation_settings["post_processing_processes"])
        for process in self.post_processing_processes:
            self.list_of_post_processing_processes.append(process)

    def Interpolate(self, delta_time):
        if (len(self.step_data_list) == 0):
            raise Exception("No completed steps found.")

        process_info = self.model_part.ProcessInfo
        current_time = process_info[Kratos.TIME]

        # calculate starting point for smoothing simulation
        required_start_time = max(current_time - delta_time * (self.number_of_smoothing_steps+1), 0.0)
        end_time = current_time - delta_time

        start_index = 0
        start_time = self.step_data_list[0][0]
        start_path = self.step_data_list[0][2]
        for index, step_data in enumerate(self.step_data_list):
            if start_time >= required_start_time:
                break
            start_index = index
            start_time = step_data[0]
            start_path = step_data[2]

        primal_parameters = self.primal_parameters.Clone()

        # set start time and end time
        primal_parameters["problem_data"]["start_time"].SetDouble(start_time)
        primal_parameters["problem_data"]["end_time"].SetDouble(end_time)
        primal_parameters["solver_settings"]["time_stepping"]["time_step"].SetDouble(delta_time)

        with ExecutionScope(start_path):
            if (self.echo_level > 1):
                Kratos.Logger.PrintInfo(self.__class__.__name__, "Initializing data for smmothed mapped mesh interpolation from {:s}".format(start_path))

            source_model_part = self.model[self.model_part.FullName() + "_Old"]
            h5_file = GetHDF5File("initialization_old.h5", "read_only")
            InputNodalResultsFromHDF5(source_model_part, h5_file, ["ALL_VARIABLES_FROM_FILE"])
            del h5_file

            # remove conditions from the source model_part because it cannot have
            # elements and conditions both for mapping
            Kratos.VariableUtils().SetFlag(Kratos.TO_ERASE, True, source_model_part.Conditions)
            source_model_part.RemoveConditionsFromAllLevels(Kratos.TO_ERASE)

            # create mapper
            mapper = KratosMapping.MapperFactory.CreateMapper(source_model_part, self.model_part, self.interpolation_settings["mapper_settings"].Clone())

            solution_step_variable_names_list = KratosRANS.RansVariableUtilities.GetSolutionstepVariableNamesList(self.model_part)
            variable_transfer_list = []

            # map current solution step data
            for solution_step_variable_name in solution_step_variable_names_list:
                solution_step_variable = Kratos.KratosGlobals.GetVariable(solution_step_variable_name)
                mapper.Map(solution_step_variable, solution_step_variable)
                variable_transfer_list.append([
                    solution_step_variable_name,
                    True,
                    0,
                    solution_step_variable_name,
                    True,
                    1
                ])

            process_info[Kratos.TIME] = start_time

            # now correct boundary values
            for process in self.list_of_post_processing_processes:
                process.Execute()

            # write existing model part
            Kratos.ModelPartIO(primal_parameters["solver_settings"]["model_import_settings"]["input_filename"].GetString(), Kratos.IO.WRITE | Kratos.IO.MESH_ONLY).WriteModelPart(self.model_part)

            # write initialization hdf5 output
            h5_file = GetHDF5File(self.interpolation_settings["initialization_hdf5_file"].GetString(), "truncate")
            OutputNodalResultsToHDF5(self.model_part, h5_file, ["ALL_VARIABLES_FROM_VARIABLES_LIST"])

            from KratosMultiphysics.RANSApplication.rans_analysis import RANSAnalysis
            primal_model, primal_simulation = SolveProblem(RANSAnalysis, primal_parameters, "smoothing")

            # now transfer values from the smoothing simulation to current model part
            KratosRANS.RansVariableDataTransferProcess(
                primal_model,
                self.model,
                primal_simulation._GetSolver().main_model_part.FullName(),
                self.model_part.FullName(),
                ["execute"],
                variable_transfer_list,
                2
            ).Execute()

            del self.step_data_list[start_index]

        # delete the old model part and mapper model part
        self.model.DeleteModelPart(source_model_part.FullName())

        if (self.interpolation_settings["clear_intermediate_files"].GetBool()):
            for step_data in self.step_data_list:
                rmtree(step_data[2])
            self.step_data_list = []

        process_info[Kratos.TIME] = current_time

        if (self.echo_level > 1):
            Kratos.Logger.PrintInfo(self.__class__.__name__, "Computed smoothed interpolation using time step {:f}".format(start_time))

    def FinalizeSolutionStep(self):
        process_info = self.model_part.ProcessInfo
        current_step = process_info[Kratos.STEP]
        current_time = process_info[Kratos.TIME]

        current_output_path = Path(self.output_path / "step_{:d}".format(current_step))
        current_output_path.mkdir(exist_ok=True, parents=True)

        h5_file = GetHDF5File(str(current_output_path / "initialization_old.h5"), "truncate")
        OutputNodalResultsToHDF5(self.model_part, h5_file, ["ALL_VARIABLES_FROM_VARIABLES_LIST"])
        Kratos.ModelPartIO(str(current_output_path / "mesh_old"), Kratos.IO.WRITE | Kratos.IO.MESH_ONLY).WriteModelPart(self.model_part)

        self.step_data_list.append([current_time, current_step, str(current_output_path)])

class CoupledRANSSolver(PythonSolver):
    def __init__(self, model, custom_settings):
        """RANS solver to be used with RANSFormulations

        This solver creates PythonSolver based on the RANSFormulations specified in custom_settings.

        Args:
            model (Kratos.Model): Model to be used in the solver.
            custom_settings (Kratos.Parameters): Settings to be used in the solver.
        """

        self._validate_settings_in_baseclass = True  # To be removed eventually
        self.BackwardCompatibilityHelper(custom_settings)
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
                                self.settings["formulation_settings"],
                                self.deprecated_settings)

        self.formulation.SetTimeSchemeSettings(self.settings["time_scheme_settings"])
        self.formulation.SetConstants(self.settings["constants"])
        self.formulation.SetIsPeriodic(self.settings["consider_periodic_conditions"].GetBool())

        self.is_periodic = self.formulation.IsPeriodic()

        self.formulation.SetTimeSchemeSettings(self.settings["time_scheme_settings"])
        self.formulation.SetWallFunctionSettings()
        scheme_type = self.settings["time_scheme_settings"]["scheme_type"].GetString()
        if (scheme_type == "steady"):
            self.is_steady = True
            default_variable_tolerance_data = Kratos.Parameters("""{
                "name": "PLEASE_SPECIFY_SCALAR_VARIABLE_NAME",
                "relative_tolerance": 1e-5,
                "absolute_tolerance": 1e-7
            }""")
            self.steady_variable_tolerance_data = []
            for variable_tolerance_data in self.settings["steady_coupling_tolerances"]:
                variable_tolerance_data.ValidateAndAssignDefaults(default_variable_tolerance_data)
                variable = Kratos.KratosGlobals.GetVariable(variable_tolerance_data["name"].GetString())
                if variable not in self.formulation.GetSolvingVariables():
                    Kratos.Logger.PrintWarning(self.__class__.__name__, "Variable {:s} not found in solving variables list hence ignoring specified steady state tolerances.".format(variable.Name()))
                else:
                    self.steady_variable_tolerance_data.append([
                        variable,
                        variable_tolerance_data["relative_tolerance"].GetDouble(),
                        variable_tolerance_data["absolute_tolerance"].GetDouble()])
        else:
            self.is_steady = False

        if (self.settings["time_scheme_settings"].Has("ramp_up_interval")):
            self.ramp_up_interval = self.settings["time_scheme_settings"]["ramp_up_interval"].GetVector()
            self.ramp_up_interval = sorted(self.ramp_up_interval)
            if (len(self.ramp_up_interval) != 2):
                raise Exception("Ramp up interval should only have two values indicating ramp up start time and end time.")
        else:
            self.ramp_up_interval = None

        self.main_model_part.ProcessInfo[KratosRANS.RANS_IS_STEADY] = self.is_steady

        self.is_converged = False
        self.min_buffer_size = self.formulation.GetMinimumBufferSize()
        self.main_model_part.SetBufferSize(self.min_buffer_size)
        self.move_mesh = self.settings["move_mesh"].GetBool()
        self.echo_level = self.settings["echo_level"].GetInt()

        if (self.settings["time_scheme_settings"].Has("ramp_up_interval")):
            self.ramp_up_interval = self.settings["time_scheme_settings"]["ramp_up_interval"].GetVector()
            self.ramp_up_interval = sorted(self.ramp_up_interval)
            if (len(self.ramp_up_interval) != 2):
                raise Exception("Ramp up interval should only have two values indicating ramp up start time and end time.")
        else:
            self.ramp_up_interval = None

        self.perform_adaptive_mesh_refinement = False
        if (custom_settings["adaptive_mesh_refinement_based_on_response_function"].GetBool()):
            self.perform_adaptive_mesh_refinement = True
            self.adaptive_mesh_refinement_based_on_response_function_settings = self.settings["adaptive_mesh_refinement_based_on_response_function_settings"]
            self.adaptive_mesh_refinement_interval_utility = Kratos.IntervalUtility(self.adaptive_mesh_refinement_based_on_response_function_settings)
            self.adaptive_mesh_refinement_based_on_response_function_settings.ValidateAndAssignDefaults(self.GetDefaultParameters()["adaptive_mesh_refinement_based_on_response_function_settings"])

            # create the folder structure
            self.amr_output_path = Path(self.adaptive_mesh_refinement_based_on_response_function_settings["output_path"].GetString())

            self.automatic_mesh_refinement_step = 0
            self.reform_dofs_dict = {}
            self.mesh_indication_variable = Kratos.KratosGlobals.GetVariable(self.adaptive_mesh_refinement_based_on_response_function_settings["mesh_change_indication_variable_name"].GetString())

            if not self.adaptive_mesh_refinement_based_on_response_function_settings["mesh_interpolation_settings"].Has("interpolation_method"):
                self.adaptive_mesh_refinement_based_on_response_function_settings["mesh_interpolation_settings"].AddEmptyValue("interpolation_method")
                self.adaptive_mesh_refinement_based_on_response_function_settings["mesh_interpolation_settings"]["interpolation_method"].SetString("none")

            mesh_interpolation_method_name = self.adaptive_mesh_refinement_based_on_response_function_settings["mesh_interpolation_settings"]["interpolation_method"].GetString()
            if (mesh_interpolation_method_name == "none"):
                self.mesh_interpolation_method = MeshInterpolationMethod(self.main_model_part, self.adaptive_mesh_refinement_based_on_response_function_settings)
            elif (mesh_interpolation_method_name == "smooth_mapped"):
                self.mesh_interpolation_method = SmoothMappedMeshInterpolationMethod(self.main_model_part, self.adaptive_mesh_refinement_based_on_response_function_settings)
            elif (mesh_interpolation_method_name == "mapped"):
                self.mesh_interpolation_method = MappedMeshInterpolationMethod(self.main_model_part, self.adaptive_mesh_refinement_based_on_response_function_settings)
            else:
                raise Exception(
                    r"""Unsupported mesh interpolation method requested. Supported method names are:
                            "none",
                            "smooth_mapped",
                            "mapped"
                    """)


            Kratos.Logger.PrintInfo(self.__class__.__name__,
                                            "Solver prepared for transient adaptive mesh refinement.")

        self.is_initialized = False
        Kratos.Logger.PrintInfo(self.__class__.__name__,
                                            "Solver construction finished.")

    def BackwardCompatibilityHelper(self, settings):
        self.deprecated_settings = {}
        if settings.Has("wall_function_settings"):
            IssueDeprecationWarning(self.__class__.__name__, "Using global \"wall_function_settings\". Please define leaf formulation specific settings.")
            self.deprecated_settings["wall_function_settings"] = settings["wall_function_settings"].Clone()
            settings.RemoveValue("wall_function_settings")

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
            "steady_coupling_tolerances": [
                {
                    "name": "REACTION_X",
                    "relative_tolerance":1e-5,
                    "absolute_tolerance":1e-7
                }
            ],
            "adaptive_mesh_refinement_based_on_response_function": false,
            "adaptive_mesh_refinement_based_on_response_function_settings": {
                "interval"                               : [0.0, 1e+30],
                "execute_every_nth_step"                 : 50,
                "mesh_change_indication_variable_name"   : "AUX_MESH_VAR",
                "output_path"                            : "adaptive_mesh_refinement",
                "re_calculate_time_step_after_refinement": true,
                "time_step_re_calculation_settings": {
                    "desired_maximum_cfl_number": 10.0,
                    "minimum_time_step"         : 1e-5,
                    "maximum_time_step"         : 0.01
                },
                "response_function_interpolation_error_computation_settings": {},
                "mmg_mesh_refinement_process_parameters": {
                    "strategy"               : "hessian",
                    "automatic_remesh"       : false,
                    "enforce_current"        : false,
                    "maximal_size"           : 0.5,
                    "minimal_size"           : 1e-3,
                    "hessian_strategy_parameters":{
                        "metric_variable": ["TIME_AVERAGED_VELOCITY_X", "TIME_AVERAGED_VELOCITY_Y", "TIME_AVERAGED_PRESSURE"],
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
                },
                "apply_mesh_optimization": true,
                "mmg_mesh_optimization_process_parameters": {
                    "strategy"        : "optimization",
                    "model_part_name" : "FluidModelPart",
                    "step_frequency"  : 1,
                    "automatic_remesh": true,
                    "automatic_remesh_parameters": {
                        "automatic_remesh_type": "Ratio",
                        "min_size_ratio": 1e-4,
                        "max_size_ratio": 1e+4,
                        "refer_type"    : "Mean"
                    },
                    "echo_level": 0,
                    "force_min" : true,
                    "force_max" : true
                },
                "mesh_interpolation_settings": {
                    "interpolation_method": "none"
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

        if not self.is_initialized:
            self.main_model_part.ProcessInfo[Kratos.STEP] = 0

        # If needed, create the estimate time step utility
        if (self.settings["time_stepping"]["automatic_time_step"].GetBool()):
            self.EstimateDeltaTimeUtility = self._GetAutomaticTimeSteppingUtility()

        RansVariableUtilities.AssignBoundaryFlagsToGeometries(self.main_model_part)
        self.formulation.Initialize()

        if (self.perform_adaptive_mesh_refinement):
            number_of_solving_variables = len(self.formulation.GetSolvingVariables())
            Kratos.VariableUtils().SetNonHistoricalVariable(Kratos.RESPONSE_FUNCTION_INTERPOLATION_ERROR, Kratos.Vector(number_of_solving_variables, 0.0), self.main_model_part.Nodes)

            # communicate to application to write the modified meshes
            self.main_model_part.ProcessInfo.SetValue(self.mesh_indication_variable, 1.0)

            if not self.is_initialized:
                self.mesh_interpolation_method.Initialize()

        if not self.is_initialized:
            Kratos.Logger.PrintInfo(self.__class__.__name__, self.formulation.GetInfo())

            Kratos.Logger.PrintInfo(self.__class__.__name__, "Solver initialization finished.")

        self.is_initialized = True

    def AdvanceInTime(self, current_time):
        dt = self._ComputeDeltaTime()
        new_time = current_time + dt

        self.main_model_part.CloneTimeStep(new_time)
        self.main_model_part.ProcessInfo[Kratos.STEP] += 1

        return new_time

    def InitializeSolutionStep(self):
        self.formulation.InitializeSolutionStep()

        if (self.perform_adaptive_mesh_refinement):
            self.mesh_interpolation_method.InitializeSolutionStep()

    def Predict(self):
        self.formulation.Predict()

    def SolveSolutionStep(self):
        self.formulation.SolveCouplingStep()
        self.is_converged = self.formulation.IsConverged()

        if not self.is_converged and not self.IsSteadySimulation() and self.echo_level > -1:
            msg = "Fluid solver did not converge for step " + str(self.main_model_part.ProcessInfo[Kratos.STEP]) + "\n"
            msg += "corresponding to time " + str(self.main_model_part.ProcessInfo[Kratos.TIME]) + "\n"
            Kratos.Logger.PrintWarning(self.__class__.__name__, msg)

        if self.is_converged and self.IsSteadySimulation():
            non_converged_variables = []
            for variable, relative_tolerance, absolute_tolerance in self.steady_variable_tolerance_data:
                relative_error, absolute_error = RansVariableUtilities.CalculateTransientVariableConvergence(
                    self.main_model_part,
                    variable)
                info = GetConvergenceInfo(
                    variable,
                    relative_error,
                    relative_tolerance,
                    absolute_error,
                    absolute_tolerance)
                Kratos.Logger.PrintInfo(self.__class__.__name__, info)

                current_variable_convergence = (relative_error <= relative_tolerance or absolute_error <= absolute_tolerance)
                if not current_variable_convergence:
                    non_converged_variables.append(variable.Name())

                self.is_converged = self.is_converged and current_variable_convergence

            if len(non_converged_variables) > 0:
                Kratos.Logger.PrintInfo(self.__class__.__name__, "Following variables are not converged to steady state tolerances hence continuing psedo time stepping: \n\t{:s}".format("\n\t".join(non_converged_variables)))

        return self.is_converged

    def FinalizeSolutionStep(self):
        self.formulation.FinalizeSolutionStep()

        if (self.perform_adaptive_mesh_refinement):
            self.mesh_interpolation_method.FinalizeSolutionStep()

    def Finalize(self):
        self.formulation.Finalize()

    def Check(self):
        self.formulation.Check()

    def Clear(self):
        self.formulation.Clear()

    def IsSteadySimulation(self):
        return self.is_steady

    def IsConverged(self):
        if (self.IsSteadySimulation()):
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

    def CheckAndExecuteAdaptiveMeshRefinement(self):
        if (self.perform_adaptive_mesh_refinement):
            if (self.adaptive_mesh_refinement_interval_utility.IsInInterval(self.main_model_part.ProcessInfo[Kratos.TIME])):
                if (self.automatic_mesh_refinement_step % self.adaptive_mesh_refinement_based_on_response_function_settings["execute_every_nth_step"].GetInt() == 0):
                    self.formulation.ComputeTransientResponseFunctionInterpolationError(
                        self.adaptive_mesh_refinement_based_on_response_function_settings["response_function_interpolation_error_computation_settings"].Clone(),
                        self.amr_output_path,
                        self.settings["model_import_settings"])

                    # this will interpolate all the nodal variables
                    self._PerformAdaptiveMeshRefinement()

                    # communicate to applications to write the modified mesh
                    self.main_model_part.ProcessInfo[self.mesh_indication_variable] = 1.0
                else:
                    # communicate to applications to not to write the modified mesh
                    self.main_model_part.ProcessInfo[self.mesh_indication_variable] = 0.0

                self.automatic_mesh_refinement_step += 1
                return (self.main_model_part.ProcessInfo[self.mesh_indication_variable] == 1.0)
            else:
                self.main_model_part.ProcessInfo[self.mesh_indication_variable] = 0.0
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

    def _PerformAdaptiveMeshRefinement(self):
        adaptive_mesh_refinement_settings = self.settings["adaptive_mesh_refinement_based_on_response_function_settings"]
        # perform adaptive remeshing
        Kratos.Logger.PrintInfo(self.__class__.__name__, "Performing adaptive mesh refinement based on response function error analysis...")
        from KratosMultiphysics.MeshingApplication.mmg_process import MmgProcess
        remeshing_process = MmgProcess(self.main_model_part.GetModel(), adaptive_mesh_refinement_settings["mmg_mesh_refinement_process_parameters"].Clone())
        remeshing_process.ExecuteInitialize()
        remeshing_process.ExecuteInitializeSolutionStep()
        remeshing_process.ExecuteFinalizeSolutionStep()
        remeshing_process.ExecuteFinalize()

        if adaptive_mesh_refinement_settings["apply_mesh_optimization"].GetBool():
            Kratos.Logger.PrintInfo(self.__class__.__name__, "Optimizing refined mesh...")
            remeshing_process = MmgProcess(self.main_model_part.GetModel(), adaptive_mesh_refinement_settings["mmg_mesh_optimization_process_parameters"].Clone())
            remeshing_process.ExecuteInitialize()
            remeshing_process.ExecuteInitializeSolutionStep()
            remeshing_process.ExecuteFinalizeSolutionStep()
            remeshing_process.ExecuteFinalize()

        # re-calculate time step if required
        if adaptive_mesh_refinement_settings["re_calculate_time_step_after_refinement"].GetBool():
            time_step_re_calculation_settings = adaptive_mesh_refinement_settings["time_step_re_calculation_settings"]
            time_step_re_calculation_settings.ValidateAndAssignDefaults(self.GetDefaultParameters()["adaptive_mesh_refinement_based_on_response_function_settings"]["time_step_re_calculation_settings"])

            # now calculate the max cfl number in the domain
            KratosCFD.FluidCharacteristicNumbersUtilities.CalculateLocalCFL(self.main_model_part)
            import KratosMultiphysics.StatisticsApplication as KratosStats
            max_cfl, _ = KratosStats.SpatialMethods.NonHistorical.Elements.NormMethods.Max(self.main_model_part, Kratos.CFL_NUMBER, "value")

            # now calculate the new time step
            if (max_cfl > 0.0):
                desired_cfl_number_after_refinement = time_step_re_calculation_settings["desired_maximum_cfl_number"].GetDouble()
                new_time_step = min(time_step_re_calculation_settings["maximum_time_step"].GetDouble(), max(self._ComputeDeltaTime() * desired_cfl_number_after_refinement / max_cfl, time_step_re_calculation_settings["minimum_time_step"].GetDouble()))
                self.settings["time_stepping"]["time_step"].SetDouble(new_time_step)
                Kratos.Logger.PrintInfo(self.__class__.__name__, "Estimated max cfl number of refined mesh is {:f}, therefore time step is changed to {:f} to adhere to desired cfl number of {:f}.".format(max_cfl, new_time_step, desired_cfl_number_after_refinement))

        self.mesh_interpolation_method.Interpolate(self._ComputeDeltaTime())

        # now remove all the formulation model parts, because they have to be created
        # from the refined mesh again
        self.formulation.DeleteModelPartsRecursively()

        Kratos.Logger.PrintInfo(self.__class__.__name__, "Finished adaptive mesh refinement.")