from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.ChimeraApplication as KratosChimera
# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.TrilinosApplication as KratosTrilinos     # MPI solvers
from KratosMultiphysics.ChimeraApplication import TrilinosExtension as TrilinosChimera
from KratosMultiphysics.FluidDynamicsApplication import TrilinosExtension as TrilinosFluid
from KratosMultiphysics.ChimeraApplication import chimera_setup_utils

# Import base class file
from KratosMultiphysics.FluidDynamicsApplication.trilinos_navier_stokes_solver_vmsmonolithic import TrilinosNavierStokesSolverMonolithic

def CreateSolver(main_model_part, custom_settings):
    return TrilinosChimeraNavierStokesSolverMonolithic(main_model_part, custom_settings)

class TrilinosChimeraNavierStokesSolverMonolithic(TrilinosNavierStokesSolverMonolithic):
    def __init__(self, model, custom_settings):
        # self.chimera_settings = custom_settings["chimera_settings"].Clone()
        # custom_settings.RemoveValue("chimera_settings")
        [self.chimera_settings, self.chimera_internal_parts, custom_settings] = chimera_setup_utils.SeparateAndValidateChimeraSettings(custom_settings)

        super(TrilinosChimeraNavierStokesSolverMonolithic,self).__init__(model,custom_settings)
        KratosMultiphysics.Logger.PrintInfo("TrilinosChimeraNavierStokesSolverMonolithic", "Construction of NavierStokesSolverMonolithic finished.")

    def AddVariables(self):
        super(TrilinosChimeraNavierStokesSolverMonolithic,self).AddVariables()
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosChimera.CHIMERA_DISTANCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosChimera.ROTATIONAL_ANGLE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosChimera.ROTATIONAL_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosChimera.ROTATION_MESH_DISPLACEMENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosChimera.ROTATION_MESH_VELOCITY)

        KratosMultiphysics.Logger.PrintInfo("TrilinosChimeraNavierStokesSolverMonolithic", "Fluid solver variables added correctly.")

    def ImportModelPart(self):
        if(self.settings["model_import_settings"]["input_type"].GetString() == "chimera"):
            chimera_mp_import_settings = []
            for chimera_part_levels in self.chimera_settings["chimera_parts"]:
                for chimera_part_settings in chimera_part_levels:
                    chimera_mp_import_settings.append( chimera_part_settings["model_import_settings"].Clone() )

            material_file_name = self.settings["material_import_settings"]["materials_filename"].GetString()
            import KratosMultiphysics.ChimeraApplication.chimera_modelpart_import as chim_mp_imp
            chim_mp_imp.ImportChimeraModelparts(self.main_model_part, chimera_mp_import_settings, material_file=material_file_name, parallel_type="MPI")
            KratosMultiphysics.Logger.PrintInfo("NavierStokesSolverMonolithicChimera", " Import of all chimera modelparts completed.")
        else:# we can use the default implementation in the base class
            super(TrilinosChimeraNavierStokesSolverMonolithic,self).ImportModelPart()

    def PrepareModelPart(self):
        if( not self.settings["model_import_settings"]["input_type"].GetString() == "chimera"):
            # Call the base solver to do the PrepareModelPart
            # Note that his also calls the PrepareModelPart of the turbulence model
            super(TrilinosChimeraNavierStokesSolverMonolithic, self).PrepareModelPart()
        else :
            super(TrilinosNavierStokesSolverMonolithic,self).PrepareModelPart()


    def Initialize(self):
        # Call the base solver to create the solution strategy
        super(TrilinosChimeraNavierStokesSolverMonolithic,self).Initialize()

        # Chimera utilities initialization
        self.chimera_process = chimera_setup_utils.GetApplyChimeraProcess(
            self.model,
            self.chimera_settings,
            self.settings)
        chimera_setup_utils.SetChimeraInternalPartsFlag(
            self.model,
            self.chimera_internal_parts)

    def GetComputingModelPart(self):
        return self.main_model_part

    def InitializeSolutionStep(self):
        self.chimera_process.ExecuteInitializeSolutionStep()
        super(TrilinosChimeraNavierStokesSolverMonolithic,self).InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        super(TrilinosChimeraNavierStokesSolverMonolithic,self).FinalizeSolutionStep()
        ## Depending on the setting this will clear the created constraints
        self.chimera_process.ExecuteFinalizeSolutionStep()

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
        ## Construct the Trilinos builder and solver
        if self.settings["consider_periodic_conditions"].GetBool() == True:
            KratosMultiphysics.Logger.PrintInfo("NavierStokesSolverMonolithicForChimera Periodic conditions are not implemented in this case .")
            raise NotImplementedError
        else:
            # TODO: This should be trilinos version
            builder_and_solver =  TrilinosChimera.TrilinosChimeraBlockBuilderAndSolverType(epetra_communicator,guess_row_size,trilinos_linear_solver)
            return builder_and_solver