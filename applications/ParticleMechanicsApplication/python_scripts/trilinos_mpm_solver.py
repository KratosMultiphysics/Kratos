from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Importing MPI extensions to Kratos
from KratosMultiphysics.mpi.distributed_import_model_part_utility import DistributedImportModelPartUtility

# Import applications and dependencies
import KratosMultiphysics.ParticleMechanicsApplication as KratosParticle
import KratosMultiphysics.TrilinosApplication as TrilinosApplication
from KratosMultiphysics.TrilinosApplication import trilinos_linear_solver_factory
from KratosMultiphysics.ParticleMechanicsApplication import TrilinosExtension as TrilinosExtension

#Other imports (Exchange this!!!)
from KratosMultiphysics.StructuralMechanicsApplication import trilinos_convergence_criteria_factory as convergence_criteria_factory
# Importing the base class
from KratosMultiphysics.ParticleMechanicsApplication.mpm_solver import MPMSolver

def CreateSolver(model, custom_settings):
    return TrilinosMPMSolver(model, custom_settings)

class TrilinosMPMSolver(MPMSolver):
    def __init__(self, model, custom_settings):
        # Set defaults and validate custom settings in the base class.
        # Construct the base solver.
        super(TrilinosMPMSolver, self).__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[TrilinosMPMSolver]:: ", "Construction is finished.")

    def AddVariables(self):
        super(TrilinosMPMSolver, self).AddVariables()
        self.grid_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)
        KratosMultiphysics.Logger.PrintInfo("::[TrilinosMPMSolver]:: ", "Variables ADDED")

    def ImportModelPart(self):
        KratosMultiphysics.Logger.PrintInfo("::[TrilinosMPMSolver]:: ", "Importing model part.")
        grid_model_part_settings = self.settings["grid_model_import_settings"]
        mdpa_file_name = grid_model_part_settings["input_filename"].GetString()
        if(self.settings["grid_model_import_settings"]["input_type"].GetString() == "mdpa"):
            self.distributed_model_part_importer = DistributedImportModelPartUtility(self.grid_model_part, self.settings)
            self.distributed_model_part_importer.ImportModelPart()
        else:
            raise Exception("Other input options are not implemented yet.")
        KratosMultiphysics.Logger.PrintInfo("::[TrilinosMPMSolver]:: ", "Finished importing model part.")
        # reading the model part of the material point
        if(self.settings["model_import_settings"]["input_type"].GetString() == "mdpa"):
            self._ImportModelPart(self.initial_mesh_model_part, self.settings["model_import_settings"])
        else:
            raise Exception("Other input options are not implemented yet.")

    def PrepareModelPart(self):
        super(TrilinosMPMSolver, self).PrepareModelPart()
        # Construct the mpi-communicator
        self.distributed_model_part_importer.CreateCommunicators()
        KratosMultiphysics.Logger.PrintInfo("::[TrilinosMPMSolver]::", "ModelPart prepared for Solver.")
        # Copy mpi-communicator to material_point_model_part
        KratosParticle.CopyMPICommunicator(self.grid_model_part, self.material_point_model_part)

    def Finalize(self):
        super(TrilinosMPMSolver, self).Finalize()
        self._GetSolutionStrategy().Clear() # needed for proper finalization of MPI


#### Specific internal functions ####
    def _GetEpetraCommunicator(self):
        if not hasattr(self, '_epetra_communicator'):
            self._epetra_communicator = self._CreateEpetraCommunicator()
        return self._epetra_communicator

#### Protected functions ####
    def _CreateEpetraCommunicator(self):
        return TrilinosApplication.CreateCommunicator()

    def _CreateConvergenceCriteria(self):
        convergence_criterion = convergence_criteria_factory.convergence_criterion(self._GetConvergenceCriterionSettings())
        return convergence_criterion.mechanical_convergence_criterion

    def _CreateLinearSolver(self):
        return trilinos_linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

    def _GetConvergenceCriterionSettings(self):
        # Create an auxiliary Kratos parameters object to store the convergence settings.
        conv_params = KratosMultiphysics.Parameters("{}")
        conv_params.AddValue("convergence_criterion",self.settings["convergence_criterion"])
        conv_params.AddValue("echo_level",self.settings["echo_level"])
        conv_params.AddValue("displacement_relative_tolerance",self.settings["displacement_relative_tolerance"])
        conv_params.AddValue("displacement_absolute_tolerance",self.settings["displacement_absolute_tolerance"])
        conv_params.AddValue("residual_relative_tolerance",self.settings["residual_relative_tolerance"])
        conv_params.AddValue("residual_absolute_tolerance",self.settings["residual_absolute_tolerance"])

        return conv_params

    def _CreateBuilderAndSolver(self):
        linear_solver = self._GetLinearSolver()
        epetra_communicator = self._GetEpetraCommunicator()
        if(self.material_point_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
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
        linear_solver = self._GetLinearSolver()
        convergence_criterion = self._GetConvergenceCriteria()
        builder_and_solver = self._GetBuilderAndSolver()
        reform_dofs_at_each_step = False
        return TrilinosExtension.TrilinosMPMResidualBasedNewtonRaphsonStrategy( computing_model_part,
                                                                                solution_scheme,
                                                                                linear_solver,
                                                                                convergence_criterion,
                                                                                builder_and_solver,
                                                                                self.settings["max_iteration"].GetInt(),
                                                                                self.settings["compute_reactions"].GetBool(),
                                                                                reform_dofs_at_each_step,
                                                                                self.settings["move_mesh_flag"].GetBool())


