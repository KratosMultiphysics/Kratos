from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.mpi as mpi
import KratosMultiphysics.TrilinosApplication as TrilinosApplication
import KratosMultiphysics.MetisApplication as MetisApplication
import structural_mechanics_solver

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()


def CreateSolver(main_model_part, custom_settings):
    return TrilinosMechanicalSolver(main_model_part, custom_settings)


class TrilinosMechanicalSolver(structural_mechanics_solver.MechanicalSolver):
    """The base class for trilinos structural mechanics solver.

    See structural_mechanics_solver.py for more information.
    """
    def __init__(self, main_model_part, custom_settings):
        if not custom_settings.Has("linear_solver_settings"): # Override defaults in the base class.
            linear_solver_settings = KratosMultiphysics.Parameters("""{
                "solver_type" : "Klu",
                "scaling" : false
            }""")
            custom_settings.AddValue("linear_solver_settings", linear_solver_settings)

        # Construct the base solver.
        super(TrilinosMechanicalSolver, self).__init__(main_model_part, custom_settings)
        print("::[TrilinosMechanicalSolver]:: Construction finished")

    def AddVariables(self):
        super(TrilinosMechanicalSolver, self).AddVariables()
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)
        print("::[TrilinosMechanicalSolver]:: Variables ADDED")

    def ImportModelPart(self):
        print("::[TrilinosMechanicalSolver]:: Importing model part.")
        # Construct the Trilinos import model part utility.
        import trilinos_import_model_part_utility
        TrilinosModelPartImporter = trilinos_import_model_part_utility.TrilinosImportModelPartUtility(self.main_model_part, self.settings)
        # Execute the Metis partitioning and reading.
        TrilinosModelPartImporter.ExecutePartitioningAndReading()
        # Call the base class execute after reading (check and prepare model process and set the constitutive law).
        super(TrilinosMechanicalSolver, self)._execute_after_reading()
        # Call the base class set and fill buffer size.
        super(TrilinosMechanicalSolver, self)._set_and_fill_buffer()
        # Construct the communicators
        TrilinosModelPartImporter.CreateCommunicators()
        print ("::[TrilinosMechanicalSolver]:: Finished importing model part.")

    #### Specific internal functions ####

    def get_epetra_communicator(self):
        if not hasattr(self, '_epetra_communicator'):
            self._epetra_communicator = self._create_epetra_communicator()
        return self._epetra_communicator

    #### Private functions ####

    def _create_epetra_communicator(self):
        return TrilinosApplication.CreateCommunicator()

    def _create_convergence_criterion(self):
        # Create an auxiliary Kratos parameters object to store the convergence settings.
        conv_params = KratosMultiphysics.Parameters("{}")
        conv_params.AddValue("convergence_criterion",self.settings["convergence_criterion"])
        conv_params.AddValue("rotation_dofs",self.settings["rotation_dofs"])
        conv_params.AddValue("echo_level",self.settings["echo_level"])
        conv_params.AddValue("component_wise",self.settings["component_wise"])
        conv_params.AddValue("displacement_relative_tolerance",self.settings["displacement_relative_tolerance"])
        conv_params.AddValue("displacement_absolute_tolerance",self.settings["displacement_absolute_tolerance"])
        conv_params.AddValue("residual_relative_tolerance",self.settings["residual_relative_tolerance"])
        conv_params.AddValue("residual_absolute_tolerance",self.settings["residual_absolute_tolerance"])
        import trilinos_convergence_criteria_factory
        convergence_criterion = trilinos_convergence_criteria_factory.convergence_criterion(conv_params)
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
