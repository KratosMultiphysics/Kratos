# Importing the Kratos Library
import KratosMultiphysics

# Importing MPI extensions to Kratos
from KratosMultiphysics.mpi.distributed_import_model_part_utility import DistributedImportModelPartUtility

# Import applications
import KratosMultiphysics.TrilinosApplication as TrilinosApplication
from KratosMultiphysics.TrilinosApplication import trilinos_linear_solver_factory

# Import base class file
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_solver import MechanicalSolver

# Other imports
from KratosMultiphysics.StructuralMechanicsApplication import trilinos_convergence_criteria_factory as convergence_criteria_factory

def CreateSolver(model, custom_settings):
    return TrilinosMechanicalSolver(model, custom_settings)

class TrilinosMechanicalSolver(MechanicalSolver):
    """The base class for trilinos structural mechanics solver.

    See structural_mechanics_solver.py for more information.
    """
    def __init__(self, model, custom_settings):
        # Construct the base solver.
        super().__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[TrilinosMechanicalSolver]:: ", "Construction finished")

    @classmethod
    def GetDefaultParameters(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "multi_point_constraints_used": false,
            "linear_solver_settings" : {
                "solver_type" : "amesos",
                "amesos_solver_type" : "Amesos_Klu"
            }
        }""")
        this_defaults.AddMissingParameters(super().GetDefaultParameters())
        return this_defaults

    def AddVariables(self):
        super().AddVariables()
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)
        KratosMultiphysics.Logger.PrintInfo("::[TrilinosMechanicalSolver]:: ", "Variables ADDED")

    def ImportModelPart(self):
        KratosMultiphysics.Logger.PrintInfo("::[TrilinosMechanicalSolver]:: ", "Importing model part.")
        self.distributed_model_part_importer = DistributedImportModelPartUtility(self.main_model_part, self.settings)
        self.distributed_model_part_importer.ImportModelPart()
        KratosMultiphysics.Logger.PrintInfo("::[TrilinosMechanicalSolver]:: ", "Finished importing model part.")

    def PrepareModelPart(self):
        super().PrepareModelPart()
        # Construct the mpi-communicator
        self.distributed_model_part_importer.CreateCommunicators()
        KratosMultiphysics.Logger.PrintInfo("::[TrilinosMechanicalSolver]::", "ModelPart prepared for Solver.")

    def Finalize(self):
        super().Finalize()
        self._GetSolutionStrategy().Clear() # needed for proper finalization of MPI

    #### Specific internal functions ####

    def _GetEpetraCommunicator(self):
        if not hasattr(self, '_epetra_communicator'):
            self._epetra_communicator = self._create_epetra_communicator()
        return self._epetra_communicator

    #### Private functions ####

    def _create_epetra_communicator(self):
        return TrilinosApplication.CreateEpetraCommunicator(self.main_model_part.GetCommunicator().GetDataCommunicator())

    def _CreateConvergenceCriterion(self):
        convergence_criterion = convergence_criteria_factory.convergence_criterion(self._get_convergence_criterion_settings())
        return convergence_criterion.mechanical_convergence_criterion

    def _CreateLinearSolver(self):
        return trilinos_linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

    def _CreateBuilderAndSolver(self):
        if self.settings["multi_point_constraints_used"].GetBool():
            raise Exception("MPCs not yet implemented in MPI")

        if (self.GetComputingModelPart().NumberOfMasterSlaveConstraints() > 0):
            KratosMultiphysics.Logger.PrintWarning("Constraints are not yet implemented in MPI and will therefore not be considered!")

        linear_solver = self._GetLinearSolver()
        epetra_communicator = self._GetEpetraCommunicator()
        if(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
            guess_row_size = 15
        else:
            guess_row_size = 45

        builder_and_solver_type: str = self.settings["builder_and_solver_settings"]["type"].GetString()
        if builder_and_solver_type == "block":
            builder_and_solver = TrilinosApplication.TrilinosBlockBuilderAndSolver(epetra_communicator,
                                                                                   guess_row_size,
                                                                                   linear_solver)
        elif builder_and_solver_type == "elimination":
            builder_and_solver = TrilinosApplication.TrilinosEliminationBuilderAndSolver(epetra_communicator,
                                                                                         guess_row_size,
                                                                                         linear_solver)
        else:
            raise ValueError(f"Unsupported builder and solver type '{builder_and_solver_type}'. Options are 'block', 'elimination'.")

        return builder_and_solver

    def _create_linear_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        mechanical_scheme = self._GetScheme()
        linear_solver = self._GetLinearSolver()
        builder_and_solver = self._GetBuilderAndSolver()
        return TrilinosApplication.TrilinosLinearStrategy(computing_model_part,
                                                          mechanical_scheme,
                                                          builder_and_solver,
                                                          self.settings["compute_reactions"].GetBool(),
                                                          self.settings["reform_dofs_at_each_step"].GetBool(),
                                                          False,
                                                          self.settings["move_mesh_flag"].GetBool())

    def _create_newton_raphson_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        solution_scheme = self._GetScheme()
        linear_solver = self._GetLinearSolver()
        convergence_criterion = self._GetConvergenceCriterion()
        builder_and_solver = self._GetBuilderAndSolver()
        return TrilinosApplication.TrilinosNewtonRaphsonStrategy(computing_model_part,
                                                                 solution_scheme,
                                                                 convergence_criterion,
                                                                 builder_and_solver,
                                                                 self.settings["max_iteration"].GetInt(),
                                                                 self.settings["compute_reactions"].GetBool(),
                                                                 self.settings["reform_dofs_at_each_step"].GetBool(),
                                                                 self.settings["move_mesh_flag"].GetBool())
