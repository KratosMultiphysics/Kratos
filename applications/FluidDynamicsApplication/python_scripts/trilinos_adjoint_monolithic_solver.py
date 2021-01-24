## Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.mpi as KratosMPI

# Import applications
import KratosMultiphysics.TrilinosApplication as TrilinosApplication
import KratosMultiphysics.FluidDynamicsApplication.TrilinosExtension as KratosCFD
from KratosMultiphysics.TrilinosApplication import trilinos_linear_solver_factory
from KratosMultiphysics.FluidDynamicsApplication.adjoint_monolithic_solver import AdjointMonolithicSolver

from KratosMultiphysics.mpi.distributed_import_model_part_utility import DistributedImportModelPartUtility

def CreateSolver(main_model_part, custom_settings):
    return AdjointMonolithicMPISolver(main_model_part, custom_settings)

class AdjointMonolithicMPISolver(AdjointMonolithicSolver):

    @classmethod
    def GetDefaultParameters(cls):
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
                "solver_type" : "amgcl"
            },
            "volume_model_part_name" : "volume_model_part",
            "skin_parts": [""],
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
            },
            "consider_periodic_conditions": false,
            "assign_neighbour_elements_to_conditions": true,
            "formulation": {
                "element_type": "vms"
            }
        }""")

        default_settings.AddMissingParameters(super(AdjointMonolithicMPISolver, cls).GetDefaultParameters())
        return default_settings

    def __init__(self, model, custom_settings):
        self._validate_settings_in_baseclass=True # To be removed eventually
        super().__init__(model, custom_settings)

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Construction of AdjointMonolithicMPISolver finished.")

    def AddVariables(self):
        ## Add variables from the base class
        super(self.__class__, self).AddVariables()

        ## Add specific MPI variables
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Variables for the AdjointMonolithicMPISolver added correctly in each processor.")

    def ImportModelPart(self):
        ## Construct the distributed import model part utility
        self.distributed_model_part_importer = DistributedImportModelPartUtility(self.main_model_part, self.settings)
        ## Execute the Metis partitioning and reading
        self.distributed_model_part_importer.ExecutePartitioningAndReading()

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "MPI model reading finished.")

    def PrepareModelPart(self):
        super().PrepareModelPart()
        ## Construct the MPI communicators
        self.distributed_model_part_importer.CreateCommunicators()

    def _GetEpetraCommunicator(self):
        if not hasattr(self, '_epetra_communicator'):
            self._epetra_communicator = TrilinosApplication.CreateCommunicator()
        return self._epetra_communicator

    def _CreateScheme(self):
        response_function = self.GetResponseFunction()
        scheme_type = self.settings["scheme_settings"]["scheme_type"].GetString()
        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        if scheme_type == "bossak":
            if domain_size == 2:
                scheme = KratosCFD.TrilinosVelocityBossakAdjointScheme2D(self.settings["scheme_settings"], response_function)
            elif domain_size == 3:
                scheme = KratosCFD.TrilinosVelocityBossakAdjointScheme3D(self.settings["scheme_settings"], response_function)
            else:
                raise Exception("Invalid DOMAIN_SIZE: " + str(domain_size))
        elif scheme_type == "steady":
            if domain_size == 2:
                scheme = KratosCFD.TrilinosSimpleSteadyAdjointScheme2D(response_function)
            elif domain_size == 3:
                scheme = KratosCFD.TrilinosSimpleSteadyAdjointScheme3D(response_function)
            else:
                raise Exception("Invalid DOMAIN_SIZE: " + str(domain_size))
        else:
            raise Exception("Invalid scheme_type: " + scheme_type)
        return scheme

    def _CreateLinearSolver(self):
        linear_solver_configuration = self.settings["linear_solver_settings"]
        return trilinos_linear_solver_factory.ConstructSolver(linear_solver_configuration)

    def _CreateBuilderAndSolver(self):
        # Set the guess_row_size (guess about the number of zero entries) for the Trilinos builder and solver
        domain_size = self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        if domain_size == 3:
            guess_row_size = 20*4
        else:
            guess_row_size = 10*3
        # Construct the Trilinos builder and solver
        trilinos_linear_solver = self._GetLinearSolver()
        epetra_communicator = self._GetEpetraCommunicator()
        if self.settings["consider_periodic_conditions"].GetBool():
            builder_and_solver = TrilinosApplication.TrilinosBlockBuilderAndSolverPeriodic(
                epetra_communicator,
                guess_row_size,
                trilinos_linear_solver,
                KratosFluid.PATCH_INDEX)
        else:
            builder_and_solver = TrilinosApplication.TrilinosBlockBuilderAndSolver(
                epetra_communicator,
                guess_row_size,
                trilinos_linear_solver)
        return builder_and_solver

    def _CreateSolutionStrategy(self):
        computing_model_part = self.GetComputingModelPart()
        time_scheme = self._GetScheme()
        linear_solver = self._GetLinearSolver()
        builder_and_solver = self._GetBuilderAndSolver()
        calculate_reaction_flag = False
        reform_dof_set_at_each_step = False
        calculate_norm_dx_flag = False
        move_mesh_flag = False
        return TrilinosApplication.TrilinosLinearStrategy(
            computing_model_part,
            time_scheme,
            linear_solver,
            builder_and_solver,
            calculate_reaction_flag,
            reform_dof_set_at_each_step,
            calculate_norm_dx_flag,
            move_mesh_flag)
