from __future__ import print_function, absolute_import, division # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.mpi as KratosMPI

# Import applications
import KratosMultiphysics.TrilinosApplication as TrilinosApplication
import KratosMultiphysics.MetisApplication as MetisApplication
import KratosMultiphysics.PoromechanicsApplication as KratosPoro
from KratosMultiphysics.TrilinosApplication import trilinos_linear_solver_factory

# Import base class file
from KratosMultiphysics.PoromechanicsApplication.poromechanics_U_Pw_solver import UPwSolver

def CreateSolver(model, custom_settings):
    return MPIUPwSolver(model, custom_settings)

class MPIUPwSolver(UPwSolver):

    def __init__(self, model, custom_settings):
        super(MPIUPwSolver,self).__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("MPIUPwSolver: ", "Construction of MPI UPwSolver finished.")

    def AddVariables(self):

        super(MPIUPwSolver, self).AddVariables()

        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)

    def ImportModelPart(self):
        # Construct the import model part utility
        from KratosMultiphysics.mpi.distributed_import_model_part_utility import DistributedImportModelPartUtility
        self.distributed_model_part_importer = DistributedImportModelPartUtility(self.main_model_part, self.settings)

        ## Execute the Metis partitioning and reading
        self.distributed_model_part_importer.ImportModelPart()

    def PrepareModelPart(self):
        super(MPIUPwSolver, self).PrepareModelPart()

        # Set ProcessInfo variables
        self.main_model_part.ProcessInfo.SetValue(KratosPoro.NODAL_SMOOTHING, False)

        # Construct the communicators
        self.distributed_model_part_importer.CreateCommunicators()

        KratosMultiphysics.Logger.PrintInfo("MPIUPwSolver: ", "Model reading finished.")

    def Initialize(self):
        self.computing_model_part = self.GetComputingModelPart()

        # Fill the previous steps of the buffer with the initial conditions
        self._FillBuffer()

        # Construct the communicator
        self.EpetraCommunicator = TrilinosApplication.CreateCommunicator()

        # Construct the linear solver
        self.linear_solver = self._ConstructLinearSolver()

        # Builder and solver creation
        builder_and_solver = self._ConstructBuilderAndSolver(self.settings["block_builder"].GetBool())

        # Solution scheme creation
        self.scheme = self._ConstructScheme(self.settings["scheme_type"].GetString(),
                                         self.settings["solution_type"].GetString())

        # Get the convergence criterion
        self.convergence_criterion = self._ConstructConvergenceCriterion(self.settings["convergence_criterion"].GetString())

        # Solver creation
        self.solver = self._ConstructSolver(builder_and_solver,
                                            self.settings["strategy_type"].GetString())

        # Set echo_level
        self.SetEchoLevel(self.settings["echo_level"].GetInt())

        # Initialize Strategy
        if self.settings["clear_storage"].GetBool():
            self.Clear()

        self.solver.Initialize()

        # Check if everything is assigned correctly
        self.Check()

        KratosMultiphysics.Logger.PrintInfo("MPIUPwSolver: ", "Solver initialization finished.")


    #### Specific internal functions ####

    def _ConstructLinearSolver(self):
        return trilinos_linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

    def _ConstructBuilderAndSolver(self, block_builder):

        if(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
            guess_row_size = 15
        else:
            guess_row_size = 45

        if(block_builder):
            builder_and_solver = TrilinosApplication.TrilinosBlockBuilderAndSolver(self.EpetraCommunicator,
                                                                                   guess_row_size,
                                                                                   self.linear_solver
                                                                                   )
        else:
            builder_and_solver = TrilinosApplication.TrilinosEliminationBuilderAndSolver(self.EpetraCommunicator,
                                                                                         guess_row_size,
                                                                                         self.linear_solver
                                                                                         )

        return builder_and_solver

    def _ConstructScheme(self, scheme_type, solution_type):

        if(scheme_type == "Newmark"):
            beta = self.settings["newmark_beta"].GetDouble()
            gamma = self.settings["newmark_gamma"].GetDouble()
            theta = self.settings["newmark_theta"].GetDouble()
            rayleigh_m = self.settings["rayleigh_m"].GetDouble()
            rayleigh_k = self.settings["rayleigh_k"].GetDouble()
            self.main_model_part.ProcessInfo.SetValue(KratosStructural.RAYLEIGH_ALPHA,rayleigh_m)
            self.main_model_part.ProcessInfo.SetValue(KratosStructural.RAYLEIGH_BETA,rayleigh_k)
            if(solution_type == "quasi_static"):
                if(rayleigh_m<1.0e-20 and rayleigh_k<1.0e-20):
                    scheme = KratosPoro.TrilinosNewmarkQuasistaticUPwScheme(beta,gamma,theta)
                else:
                    scheme = KratosPoro.TrilinosNewmarkQuasistaticDampedUPwScheme(beta,gamma,theta)
            else:
                scheme = KratosPoro.TrilinosNewmarkDynamicUPwScheme(beta,gamma,theta)
        else:
            raise Exception("Apart from Newmark, other scheme_type are not available.")

        return scheme

    def _ConstructConvergenceCriterion(self, convergence_criterion):

        D_RT = self.settings["displacement_relative_tolerance"].GetDouble()
        D_AT = self.settings["displacement_absolute_tolerance"].GetDouble()
        R_RT = self.settings["residual_relative_tolerance"].GetDouble()
        R_AT = self.settings["residual_absolute_tolerance"].GetDouble()
        echo_level = self.settings["echo_level"].GetInt()

        if(convergence_criterion == "Displacement_criterion"):
            convergence_criterion = TrilinosApplication.TrilinosDisplacementCriteria(D_RT, D_AT, self.EpetraCommunicator)
            convergence_criterion.SetEchoLevel(echo_level)
        elif(convergence_criterion == "Residual_criterion"):
            convergence_criterion = TrilinosApplication.TrilinosResidualCriteria(R_RT, R_AT)
            convergence_criterion.SetEchoLevel(echo_level)
        elif(convergence_criterion == "And_criterion"):
            Displacement = TrilinosApplication.TrilinosDisplacementCriteria(D_RT, D_AT, self.EpetraCommunicator)
            Displacement.SetEchoLevel(echo_level)
            Residual = TrilinosApplication.TrilinosResidualCriteria(R_RT, R_AT)
            Residual.SetEchoLevel(echo_level)
            convergence_criterion = TrilinosApplication.TrilinosAndCriteria(Residual, Displacement)
        elif(convergence_criterion == "Or_criterion"):
            Displacement = TrilinosApplication.TrilinosDisplacementCriteria(D_RT, D_AT, self.EpetraCommunicator)
            Displacement.SetEchoLevel(echo_level)
            Residual = TrilinosApplication.TrilinosResidualCriteria(R_RT, R_AT)
            Residual.SetEchoLevel(echo_level)
            convergence_criterion = TrilinosApplication.TrilinosOrCriteria(Residual, Displacement)

        return convergence_criterion

    def _ConstructSolver(self, builder_and_solver, strategy_type):

        max_iters = self.settings["max_iteration"].GetInt()
        compute_reactions = self.settings["compute_reactions"].GetBool()
        reform_step_dofs = self.settings["reform_dofs_at_each_step"].GetBool()
        move_mesh_flag = self.settings["move_mesh_flag"].GetBool()

        if strategy_type == "newton_raphson":
            self.main_model_part.ProcessInfo.SetValue(KratosPoro.IS_CONVERGED, True)
            solving_strategy = TrilinosApplication.TrilinosNewtonRaphsonStrategy(self.main_model_part,
                                                                       self.scheme,
                                                                       self.linear_solver,
                                                                       self.convergence_criterion,
                                                                       builder_and_solver,
                                                                       max_iters,
                                                                       compute_reactions,
                                                                       reform_step_dofs,
                                                                       move_mesh_flag
                                                                       )
        else:
            raise Exception("Apart from newton_raphson, other strategy_type are not available.")

        return solving_strategy