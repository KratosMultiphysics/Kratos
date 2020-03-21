from __future__ import print_function, absolute_import, division # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.mpi as KratosMPI

# Import applications
import KratosMultiphysics.TrilinosApplication as KratosTrilinos
import KratosMultiphysics.StructuralMechanicsApplication as KratosStructural
import KratosMultiphysics.PoromechanicsApplication as KratosPoro
from KratosMultiphysics.PoromechanicsApplication import TrilinosExtension as TrilinosPoro

# Import base class file
from KratosMultiphysics.PoromechanicsApplication.poromechanics_U_Pw_solver import UPwSolver
from KratosMultiphysics.mpi.distributed_import_model_part_utility import DistributedImportModelPartUtility

def CreateSolver(model, custom_settings):
    return TrilinosUPwSolver(model, custom_settings)

class TrilinosUPwSolver(UPwSolver):

    def __init__(self, model, custom_settings):
        super(TrilinosUPwSolver,self).__init__(model, custom_settings)

        KratosMultiphysics.Logger.PrintInfo("TrilinosUPwSolver: ", "Construction of MPI UPwSolver finished.")

    def AddVariables(self):

        super(TrilinosUPwSolver, self).AddVariables()

        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)

    def ImportModelPart(self):
        # Construct the import model part utility
        self.distributed_model_part_importer = DistributedImportModelPartUtility(self.main_model_part, self.settings)
        ## Execute the Metis partitioning and reading
        self.distributed_model_part_importer.ImportModelPart()

    def PrepareModelPart(self):
        super(TrilinosUPwSolver, self).PrepareModelPart()

        # Set ProcessInfo variables
        # self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME,
        #                                           self.settings["start_time"].GetDouble())
        # self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME,
        #                                           self.settings["time_step"].GetDouble())
        # self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.STEP, 0)
        # self.main_model_part.ProcessInfo.SetValue(KratosPoro.TIME_UNIT_CONVERTER, 1.0)
        self.main_model_part.ProcessInfo.SetValue(KratosPoro.NODAL_SMOOTHING, False)

        # if not self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
        #     materials_imported = self.import_constitutive_laws()
        #     ## Set buffer size
        #     self._SetBufferSize()
        #     if materials_imported:
        #         KratosMultiphysics.Logger.PrintInfo("TrilinosUPwSolver", "Constitutive law was successfully imported via json.")
        #     else:
        #         KratosMultiphysics.Logger.PrintInfo("TrilinosUPwSolver", "Constitutive law was not successfully imported.")


        # if not self.model.HasModelPart(self.settings["model_part_name"].GetString()):
        #     self.model.AddModelPart(self.main_model_part)

        # Construct the communicators
        self.distributed_model_part_importer.CreateCommunicators()

        # stop

        KratosMultiphysics.Logger.PrintInfo("TrilinosUPwSolver: ", "Model reading finished.")

    def Initialize(self):
        # Construct the communicator
        self.EpetraCommunicator = KratosTrilinos.CreateCommunicator()

        ## Get the computing model part
        self.computing_model_part = self.GetComputingModelPart()

        # Fill the previous steps of the buffer with the initial conditions
        self._FillBuffer()

        # Construct the linear solver
        self.linear_solver = self._ConstructLinearSolver()

        # Builder and solver creation
        self.builder_and_solver = self._ConstructBuilderAndSolver(self.settings["block_builder"].GetBool())

        # Solution scheme creation
        self.scheme = self._ConstructScheme(self.settings["scheme_type"].GetString(),
                                         self.settings["solution_type"].GetString())

        # Get the convergence criterion
        self.convergence_criterion = self._ConstructConvergenceCriterion(self.settings["convergence_criterion"].GetString())

        # Solver creation
        self.solver = self._ConstructSolver(self.settings["strategy_type"].GetString())

        # Set echo_level
        self.SetEchoLevel(self.settings["echo_level"].GetInt())

        # Initialize Strategy
        if self.settings["clear_storage"].GetBool():
            self.Clear()

        self.solver.Initialize()

        # Check if everything is assigned correctly
        self.Check()

        KratosMultiphysics.Logger.PrintInfo("TrilinosUPwSolver: ", "Solver initialization finished.")

    # def GetComputingModelPart(self):
    #     return self.main_model_part

    #### Specific internal functions ####

    def _ConstructLinearSolver(self):
        from KratosMultiphysics.TrilinosApplication import trilinos_linear_solver_factory
        return trilinos_linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

    def _ConstructBuilderAndSolver(self, block_builder):

        if(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
            guess_row_size = 10*3
        else:
            guess_row_size = 20*4

        if(block_builder):
            builder_and_solver = KratosTrilinos.TrilinosBlockBuilderAndSolver(self.EpetraCommunicator,
                                                                                   guess_row_size,
                                                                                   self.linear_solver
                                                                                   )
        else:
            builder_and_solver = KratosTrilinos.TrilinosEliminationBuilderAndSolver(self.EpetraCommunicator,
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
                    scheme = TrilinosPoro.TrilinosNewmarkQuasistaticUPwScheme(beta,gamma,theta)
                else:
                    scheme = TrilinosPoro.TrilinosNewmarkQuasistaticDampedUPwScheme(beta,gamma,theta)
            else:
                scheme = TrilinosPoro.TrilinosNewmarkDynamicUPwScheme(beta,gamma,theta)
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
            convergence_criterion = KratosTrilinos.TrilinosDisplacementCriteria(D_RT, D_AT)
            convergence_criterion.SetEchoLevel(echo_level)
        elif(convergence_criterion == "Residual_criterion"):
            convergence_criterion = KratosTrilinos.TrilinosResidualCriteria(R_RT, R_AT)
            convergence_criterion.SetEchoLevel(echo_level)
        elif(convergence_criterion == "And_criterion"):
            Displacement = KratosTrilinos.TrilinosDisplacementCriteria(D_RT, D_AT)
            Displacement.SetEchoLevel(echo_level)
            Residual = KratosTrilinos.TrilinosResidualCriteria(R_RT, R_AT)
            Residual.SetEchoLevel(echo_level)
            convergence_criterion = KratosTrilinos.TrilinosAndCriteria(Residual, Displacement)
        elif(convergence_criterion == "Or_criterion"):
            Displacement = KratosTrilinos.TrilinosDisplacementCriteria(D_RT, D_AT)
            Displacement.SetEchoLevel(echo_level)
            Residual = KratosTrilinos.TrilinosResidualCriteria(R_RT, R_AT)
            Residual.SetEchoLevel(echo_level)
            convergence_criterion = KratosTrilinos.TrilinosOrCriteria(Residual, Displacement)

        return convergence_criterion

    def _ConstructSolver(self, strategy_type):

        max_iters = self.settings["max_iteration"].GetInt()
        compute_reactions = self.settings["compute_reactions"].GetBool()
        reform_step_dofs = self.settings["reform_dofs_at_each_step"].GetBool()
        move_mesh_flag = self.settings["move_mesh_flag"].GetBool()

        if strategy_type == "newton_raphson":
            self.main_model_part.ProcessInfo.SetValue(KratosPoro.IS_CONVERGED, True)
            solving_strategy = KratosTrilinos.TrilinosNewtonRaphsonStrategy(self.main_model_part,
                                                                       self.scheme,
                                                                       self.linear_solver,
                                                                       self.convergence_criterion,
                                                                       self.builder_and_solver,
                                                                       max_iters,
                                                                       compute_reactions,
                                                                       reform_step_dofs,
                                                                       move_mesh_flag
                                                                       )
        else:
            raise Exception("Apart from newton_raphson, other strategy_type are not available.")

        return solving_strategy