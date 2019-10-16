from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
## Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.mpi as KratosMPI

# Import applications
import KratosMultiphysics.TrilinosApplication as TrilinosApplication
from KratosMultiphysics.TrilinosApplication import trilinos_linear_solver_factory
import KratosMultiphysics.FluidDynamicsApplication as FluidDynamicsApplication
from KratosMultiphysics.FluidDynamicsApplication.adjoint_vmsmonolithic_solver import AdjointVMSMonolithicSolver

from KratosMultiphysics.mpi.distributed_import_model_part_utility import DistributedImportModelPartUtility

def CreateSolver(main_model_part, custom_settings):
    return AdjointVMSMonolithicMPISolver(main_model_part, custom_settings)

class AdjointVMSMonolithicMPISolver(AdjointVMSMonolithicSolver):

    @classmethod
    def GetDefaultSettings(cls):
        # default settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "trilinos_adjoint_vmsmonolithic_solver",
            "scheme_settings" : {
                "scheme_type": "bossak"
            },
            "response_function_settings" : {
                "response_type" : "drag"
            },
            "sensitivity_settings" : {},
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name"
            },
            "material_import_settings": {
                "materials_filename": ""
            },
            "linear_solver_settings" : {
                "solver_type" : "multi_level"
            },
            "volume_model_part_name" : "volume_model_part",
            "skin_parts": [""],
            "dynamic_tau": 0.0,
            "oss_switch": 0,
            "echo_level": 0,
            "time_stepping"               : {
                "automatic_time_step" : false,
                "time_step"           : -0.1
            },
            "domain_size": -1,
            "model_part_name": "",
            "time_stepping": {
                "automatic_time_step" : false,
                "time_step"           : -0.1
            }
        }""")

        default_settings.AddMissingParameters(super(AdjointVMSMonolithicMPISolver, cls).GetDefaultSettings())
        return default_settings

    def __init__(self, model, custom_settings):
        self._validate_settings_in_baseclass=True # To be removed eventually
        super(AdjointVMSMonolithicSolver, self).__init__(model, custom_settings)

        self.element_name = "VMSAdjointElement"
        if self.settings["domain_size"].GetInt() == 2:
            self.condition_name = "LineCondition"
        elif self.settings["domain_size"].GetInt() == 3:
            self.condition_name = "SurfaceCondition"
        self.min_buffer_size = 2
        self.element_has_nodal_properties = True

        # construct the linear solver
        self.trilinos_linear_solver = trilinos_linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Construction of AdjointVMSMonolithicMPISolver finished.")

    def AddVariables(self):
        ## Add variables from the base class
        super(self.__class__, self).AddVariables()

        ## Add specific MPI variables
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Variables for the AdjointVMSMonolithicMPISolver added correctly in each processor.")

    def ImportModelPart(self):
        ## Construct the distributed import model part utility
        self.distributed_model_part_importer = DistributedImportModelPartUtility(self.main_model_part, self.settings)
        ## Execute the Metis partitioning and reading
        self.distributed_model_part_importer.ExecutePartitioningAndReading()

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "TrilinosNavierStokesSolverMonolithic","MPI model reading finished.")

    def PrepareModelPart(self):
        super(self.__class__,self).PrepareModelPart()
        ## Construct the MPI communicators
        self.distributed_model_part_importer.CreateCommunicators()

    def AddDofs(self):
        super(self.__class__, self).AddDofs()

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "DOFs for the AdjointVMSMonolithicMPISolver added correctly in all processors.")

    def Initialize(self):

        ## Construct the communicator
        self.EpetraCommunicator = TrilinosApplication.CreateCommunicator()

        ## Get the computing model part
        self.computing_model_part = self.GetComputingModelPart()

        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        if self.settings["response_function_settings"]["response_type"].GetString() == "drag":
            if (domain_size == 2):
                self.response_function = FluidDynamicsApplication.DragResponseFunction2D(self.settings["response_function_settings"]["custom_settings"], self.main_model_part)
            elif (domain_size == 3):
                self.response_function = FluidDynamicsApplication.DragResponseFunction3D(self.settings["response_function_settings"]["custom_settings"], self.main_model_part)
            else:
                raise Exception("Invalid DOMAIN_SIZE: " + str(domain_size))
        else:
            raise Exception("invalid response_type: " + self.settings["response_function_settings"]["response_type"].GetString())

        self.sensitivity_builder = KratosMultiphysics.SensitivityBuilder(self.settings["sensitivity_settings"], self.main_model_part, self.response_function)

        if self.settings["scheme_settings"]["scheme_type"].GetString() == "bossak":
            self.time_scheme = TrilinosApplication.TrilinosResidualBasedAdjointBossakScheme(self.settings["scheme_settings"], self.response_function)
        elif self.settings["scheme_settings"]["scheme_type"].GetString() == "steady":
            self.time_scheme = TrilinosApplication.TrilinosResidualBasedAdjointSteadyScheme(self.response_function)
        else:
            raise Exception("invalid scheme_type: " + self.settings["scheme_settings"]["scheme_type"].GetString())

        if self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 3:
            guess_row_size = 20*4
        elif self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
            guess_row_size = 10*3

        self.builder_and_solver = TrilinosApplication.TrilinosBlockBuilderAndSolver(self.EpetraCommunicator,
                                                                               guess_row_size,
                                                                               self.trilinos_linear_solver)

        self.solver = TrilinosApplication.TrilinosLinearStrategy(self.main_model_part,
                                                                 self.time_scheme,
                                                                 self.trilinos_linear_solver,
                                                                 self.builder_and_solver,
                                                                 False,
                                                                 False,
                                                                 False,
                                                                 False)

        (self.solver).SetEchoLevel(self.settings["echo_level"].GetInt())

        (self.solver).Initialize()
        (self.response_function).Initialize()
        (self.sensitivity_builder).Initialize()

        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DYNAMIC_TAU, self.settings["dynamic_tau"].GetDouble())
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.OSS_SWITCH, self.settings["oss_switch"].GetInt())

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__,"Monolithic MPI solver initialization finished.")
