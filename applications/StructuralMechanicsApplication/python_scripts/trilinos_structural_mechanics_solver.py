from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.mpi as KratosMPI

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("StructuralMechanicsApplication","TrilinosApplication")

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.TrilinosApplication as TrilinosApplication

# Import base class file
import structural_mechanics_solver


def CreateSolver(model, custom_settings):
    return TrilinosMechanicalSolver(model, custom_settings)


class TrilinosMechanicalSolver(structural_mechanics_solver.MechanicalSolver):
    """The base class for trilinos structural mechanics solver.

    See structural_mechanics_solver.py for more information.
    """
    def __init__(self, model, custom_settings):
        if not custom_settings.Has("linear_solver_settings"): # Override defaults in the base class.
            linear_solver_settings = KratosMultiphysics.Parameters("""{
                "solver_type" : "AmesosSolver",
                "amesos_solver_type" : "Amesos_Klu"
            }""")
            custom_settings.AddValue("linear_solver_settings", linear_solver_settings)

        # Construct the base solver.
        super(TrilinosMechanicalSolver, self).__init__(model, custom_settings)
        self.print_on_rank_zero("::[TrilinosMechanicalSolver]:: ", "Construction finished")

    def AddVariables(self):
        if not self.is_restarted():
            super(TrilinosMechanicalSolver, self).AddVariables()
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)
            self.print_on_rank_zero("::[TrilinosMechanicalSolver]:: ", "Variables ADDED")

    def ReadModelPart(self):
        self.print_on_rank_zero("::[TrilinosMechanicalSolver]:: ", "Reading model part.")
        if self.is_restarted():
            self.get_restart_utility().LoadRestart()
        elif(self.settings["model_import_settings"]["input_type"].GetString() == "mdpa"):
            # Construct the Trilinos import model part utility.
            import trilinos_import_model_part_utility
            self.trilinos_model_part_importer = trilinos_import_model_part_utility.TrilinosImportModelPartUtility(self.main_model_part, self.settings)
            # Execute the Metis partitioning and reading.
            self.trilinos_model_part_importer.ExecutePartitioningAndReading()
        else:
            raise Exception("Other model part input options are not yet implemented.")
        self.print_on_rank_zero("::[TrilinosMechanicalSolver]:: ", "Finished reading model part.")

    def PrepareModel(self):
        if not self.is_restarted():
            super(TrilinosMechanicalSolver, self).PrepareModelPart()
            # Construct the communicators
            self.trilinos_model_part_importer.CreateCommunicators()
        KratosMultiphysics.Logger.PrintInfo("::[TrilinosMechanicalSolver]::", "ModelPart prepared for Solver.")

    #### Specific internal functions ####

    def get_epetra_communicator(self):
        if not hasattr(self, '_epetra_communicator'):
            self._epetra_communicator = self._create_epetra_communicator()
        return self._epetra_communicator

    def print_on_rank_zero(self, *args):
        KratosMPI.mpi.world.barrier()
        if KratosMPI.mpi.rank == 0:
            KratosMultiphysics.Logger.PrintInfo(" ".join(map(str,args)))

    #### Private functions ####

    def _create_epetra_communicator(self):
        return TrilinosApplication.CreateCommunicator()

    def _create_convergence_criterion(self):
        import trilinos_convergence_criteria_factory as convergence_criteria_factory
        convergence_criterion = convergence_criteria_factory.convergence_criterion(self._get_convergence_criterion_settings())
        return convergence_criterion.mechanical_convergence_criterion

    def _create_linear_solver(self):
        import trilinos_linear_solver_factory # TODO: Is new_trilinos_linear_solver_factory or trilinos_linear_solver_factory?
        linear_solver = trilinos_linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])
        return linear_solver

    def _create_builder_and_solver(self):
        if self.settings["multi_point_constraints_used"].GetBool():
            raise Exception("MPCs not yet implemented in MPI")

        linear_solver = self.get_linear_solver()
        epetra_communicator = self.get_epetra_communicator()
        if(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
            guess_row_size = 15
        else:
            guess_row_size = 45
        if(self.settings["block_builder"].GetBool() == True):
            builder_and_solver = TrilinosApplication.TrilinosBlockBuilderAndSolver(epetra_communicator,
                                                                                   guess_row_size,
                                                                                   linear_solver)
        else:
            builder_and_solver = TrilinosApplication.TrilinosEliminationBuilderAndSolver(epetra_communicator,
                                                                                         guess_row_size,
                                                                                         linear_solver)
        return builder_and_solver

    def _create_linear_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        mechanical_scheme = self.get_solution_scheme()
        linear_solver = self.get_linear_solver()
        builder_and_solver = self.get_builder_and_solver()
        return TrilinosApplication.TrilinosLinearStrategy(computing_model_part,
                                                          mechanical_scheme,
                                                          linear_solver,
                                                          builder_and_solver,
                                                          self.settings["compute_reactions"].GetBool(),
                                                          self.settings["reform_dofs_at_each_step"].GetBool(),
                                                          False,
                                                          self.settings["move_mesh_flag"].GetBool())

    def _create_newton_raphson_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        solution_scheme = self.get_solution_scheme()
        linear_solver = self.get_linear_solver()
        convergence_criterion = self.get_convergence_criterion()
        builder_and_solver = self.get_builder_and_solver()
        return TrilinosApplication.TrilinosNewtonRaphsonStrategy(computing_model_part,
                                                                 solution_scheme,
                                                                 linear_solver,
                                                                 convergence_criterion,
                                                                 builder_and_solver,
                                                                 self.settings["max_iteration"].GetInt(),
                                                                 self.settings["compute_reactions"].GetBool(),
                                                                 self.settings["reform_dofs_at_each_step"].GetBool(),
                                                                 self.settings["move_mesh_flag"].GetBool())

    def _create_restart_utility(self):
        """Create the restart utility."""
        import trilinos_restart_utility as restart_utility
        rest_utility = restart_utility.RestartUtility(self.main_model_part,
                                                      self._get_restart_settings())
        return rest_utility
