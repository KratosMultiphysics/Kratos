# Importing the Kratos Library
import KratosMultiphysics

# Importing MPI extensions to Kratos
from KratosMultiphysics.mpi.distributed_import_model_part_utility import DistributedImportModelPartUtility

# Import applications and dependencies
import KratosMultiphysics.ParticleMechanicsApplication as KratosParticle
import KratosMultiphysics.TrilinosApplication as TrilinosApplication
from KratosMultiphysics.TrilinosApplication import trilinos_linear_solver_factory
from KratosMultiphysics.ParticleMechanicsApplication import TrilinosExtension as TrilinosParticle

# Importing the base class
from KratosMultiphysics.ParticleMechanicsApplication.mpm_solver import MPMSolver

def CreateSolver(model, custom_settings):
    return TrilinosMPMSolver(model, custom_settings)

class TrilinosMPMSolver(MPMSolver):

    def __init__(self, model, custom_settings):
        super().__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[TrilinosMPMSolver]:: ", "Construction is finished.")

    def AddVariables(self):
        super().AddVariables()
        self.grid_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)
        KratosMultiphysics.Logger.PrintInfo("::[TrilinosMPMSolver]:: ", "Variables ADDED")

    def ImportModelPart(self):
        KratosMultiphysics.Logger.PrintInfo("::[TrilinosMPMSolver]:: ", "Importing model part.")
        local_import_settings = KratosMultiphysics.Parameters("""{
            "model_import_settings" : {},
            "echo_level"            : 0
            }
        """)
        local_import_settings["model_import_settings"].AddMissingParameters(self.settings["grid_model_import_settings"])
        local_import_settings["echo_level"].SetInt(self.settings["echo_level"].GetInt())

        # Reading the background grid model part
        if(self.settings["grid_model_import_settings"]["input_type"].GetString() == "mdpa"):
            self.distributed_model_part_importer = DistributedImportModelPartUtility(self.grid_model_part, local_import_settings)
            self.distributed_model_part_importer.ImportModelPart()
        else:
            raise Exception("Other input options are not implemented yet.")
        KratosMultiphysics.Logger.PrintInfo("::[TrilinosMPMSolver]:: ", "Finished importing the background grid model part.")

        # Reading the initial mesh model part
        if(self.settings["model_import_settings"]["input_type"].GetString() == "mdpa"):
            self._ImportModelPart(self.initial_mesh_model_part, self.settings["model_import_settings"])
        else:
            raise Exception("Other input options are not implemented yet.")
        KratosMultiphysics.Logger.PrintInfo("::[TrilinosMPMSolver]:: ", "Finished importing the initial mesh model part.")

    def PrepareModelPart(self):
        super().PrepareModelPart()
        # Construct the mpi-communicator
        self.distributed_model_part_importer.CreateCommunicators()
        KratosMultiphysics.Logger.PrintInfo("::[TrilinosMPMSolver]::", "ModelPart prepared for Solver.")
        # Clear elements and non-grid conditions from communciator
        TrilinosParticle.MPM_MPI_Utilities.ClearLocalElementsFromCommunicator(self.grid_model_part)
        # Copy mpi-communicator to material_point_model_part
        TrilinosParticle.MPM_MPI_Utilities.SetMPICommunicator(self.grid_model_part, self.material_point_model_part)
        KratosMultiphysics.Logger.PrintInfo("::[TrilinosMPMSolver]::", "Communicator prepared for Solver.")

    def _SearchElement(self):
        searching_alg_type = self.settings["element_search_settings"]["search_algorithm_type"].GetString()
        max_number_of_search_results = self.settings["element_search_settings"]["max_number_of_results"].GetInt()
        searching_tolerance          = self.settings["element_search_settings"]["searching_tolerance"].GetDouble()
        if (searching_alg_type == "bin_based"):
            KratosParticle.SearchElementMPI(self.grid_model_part, self.material_point_model_part, max_number_of_search_results, searching_tolerance)
        else:
            err_msg  = "The requested searching algorithm \"" + searching_alg_type
            err_msg += "\" is not available for ParticleMechanicsApplication!\n"
            err_msg += "Available options are: \"bin_based\""
            raise Exception(err_msg)

    def Finalize(self):
        super().Finalize()
        self._GetSolutionStrategy().Clear() # needed for proper finalization of MPI

    def _GetEpetraCommunicator(self):
        if not hasattr(self, '_epetra_communicator'):
            self._epetra_communicator = self._CreateEpetraCommunicator()
        return self._epetra_communicator

    def _CreateEpetraCommunicator(self):
        return TrilinosApplication.CreateCommunicator()

    def _CreateConvergenceCriteria(self):
        convergence_criterion_parameters = self._GetConvergenceCriteriaSettings()
        if (convergence_criterion_parameters["convergence_criterion"].GetString() == "residual_criterion"):
            R_RT = convergence_criterion_parameters["residual_relative_tolerance"].GetDouble()
            R_AT = convergence_criterion_parameters["residual_absolute_tolerance"].GetDouble()
            convergence_criterion = TrilinosApplication.TrilinosResidualCriteria(R_RT, R_AT)
            convergence_criterion.SetEchoLevel(convergence_criterion_parameters["echo_level"].GetInt())
        elif (convergence_criterion_parameters["convergence_criterion"].GetString() == "displacement_criterion"):
            D_RT = convergence_criterion_parameters["displacement_relative_tolerance"].GetDouble()
            D_AT = convergence_criterion_parameters["displacement_absolute_tolerance"].GetDouble()
            convergence_criterion = TrilinosApplication.TrilinosDisplacementCriteria(D_RT, D_AT)
            convergence_criterion.SetEchoLevel(convergence_criterion_parameters["echo_level"].GetInt())
        else:
            err_msg  = "The requested convergence criteria \"" + convergence_criterion_parameters["convergence_criterion"].GetString()
            err_msg += "\" is not supported for ParticleMechanicsApplication!\n"
            err_msg += "Available options are: \"residual_criterion\" or \"displacement_criterion\""
            raise Exception(err_msg)

        return convergence_criterion

    def _CreateLinearSolver(self):
        print(self.settings)
        return trilinos_linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

    def _CreateBuilderAndSolver(self):
        linear_solver = self._GetLinearSolver()
        epetra_communicator = self._GetEpetraCommunicator()
        if(self.material_point_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
            guess_row_size = 15
        else:
            guess_row_size = 45
        if(self.settings["block_builder"].GetBool()):
            return TrilinosParticle.TrilinosMPMBlockBuilderAndSolver(epetra_communicator,guess_row_size,linear_solver)
        else:
            raise Exception("Elimination Builder is not yet implementde for MPI!")

    def _CreateLinearStrategy(self):
        computing_model_part = self.GetComputingModelPart()
        mechanical_scheme = self._GetSolutionScheme()
        linear_solver = self._GetLinearSolver()
        builder_and_solver = self._GetBuilderAndSolver()
        reform_dofs_at_each_step = False
        calc_norm_dx_flag = False
        return TrilinosApplication.TrilinosLinearStrategy(computing_model_part,
                                                          mechanical_scheme,
                                                          linear_solver,
                                                          builder_and_solver,
                                                          self.settings["compute_reactions"].GetBool(),
                                                          reform_dofs_at_each_step,
                                                          calc_norm_dx_flag,
                                                          self.settings["move_mesh_flag"].GetBool())

    def _CreateNewtonRaphsonStrategy(self):
        computing_model_part = self.GetComputingModelPart()
        solution_scheme = self._GetSolutionScheme()
        convergence_criterion = self._GetConvergenceCriteria()
        builder_and_solver = self._GetBuilderAndSolver()
        reform_dofs_at_each_step = False
        return TrilinosParticle.TrilinosMPMResidualBasedNewtonRaphsonStrategy(computing_model_part,
                                                                              solution_scheme,
                                                                              convergence_criterion,
                                                                              builder_and_solver,
                                                                              self.settings["max_iteration"].GetInt(),
                                                                              self.settings["compute_reactions"].GetBool(),
                                                                              reform_dofs_at_each_step,
                                                                              self.settings["move_mesh_flag"].GetBool())
