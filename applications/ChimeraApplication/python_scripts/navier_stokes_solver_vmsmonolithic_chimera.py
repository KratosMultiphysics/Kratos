
# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.ChimeraApplication as KratosChimera
from KratosMultiphysics.ChimeraApplication import chimera_setup_utils

# Import base class file
from KratosMultiphysics.FluidDynamicsApplication.navier_stokes_monolithic_solver import NavierStokesMonolithicSolver

def CreateSolver(main_model_part, custom_settings):
    return NavierStokesSolverMonolithicChimera(main_model_part, custom_settings)

class NavierStokesSolverMonolithicChimera(NavierStokesMonolithicSolver):
    def __init__(self, model, custom_settings):
        [self.chimera_settings, self.chimera_internal_parts, custom_settings] = chimera_setup_utils.SeparateAndValidateChimeraSettings(custom_settings)
        super(NavierStokesSolverMonolithicChimera,self).__init__(model,custom_settings)
        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Construction of NavierStokesSolverMonolithicChimera finished.")

    def AddVariables(self):
        super(NavierStokesSolverMonolithicChimera,self).AddVariables()
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosChimera.CHIMERA_DISTANCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosChimera.ROTATIONAL_ANGLE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosChimera.ROTATIONAL_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosChimera.ROTATION_MESH_DISPLACEMENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosChimera.ROTATION_MESH_VELOCITY)

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Fluid chimera solver variables added correctly.")

    def ImportModelPart(self):
        if(self.settings["model_import_settings"]["input_type"].GetString() == "chimera"):
            chimera_mp_import_settings = []
            for chimera_part_levels in self.chimera_settings["chimera_parts"]:
                for chimera_part_settings in chimera_part_levels:
                    chimera_mp_import_settings.append( chimera_part_settings["model_import_settings"].Clone() )

            material_file_name = self.settings["material_import_settings"]["materials_filename"].GetString()
            import KratosMultiphysics.ChimeraApplication.chimera_modelpart_import as chim_mp_imp
            chim_mp_imp.ImportChimeraModelparts(self.main_model_part, chimera_mp_import_settings, material_file=material_file_name, parallel_type="OpenMP")
            KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, " Import of all chimera modelparts completed.")
        else:# we can use the default implementation in the base class
            super(NavierStokesSolverMonolithicChimera,self).ImportModelPart()

    def Initialize(self):
        # Call the base solver to create the solution strategy
        super(NavierStokesSolverMonolithicChimera,self).Initialize()

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
        super(NavierStokesSolverMonolithicChimera,self).InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        super(NavierStokesSolverMonolithicChimera,self).FinalizeSolutionStep()
        ## Depending on the setting this will clear the created constraints
        self.chimera_process.ExecuteFinalizeSolutionStep()

    def _CreateBuilderAndSolver(self):
        linear_solver = self._GetLinearSolver()
        if self.settings["consider_periodic_conditions"].GetBool():
            KratosMultiphysics.Logger.PrintInfo("NavierStokesSolverMonolithicForChimera Periodic conditions are not implemented in this case .")
            raise NotImplementedError
        else:
            builder_and_solver = KratosChimera.ResidualBasedBlockBuilderAndSolverWithConstraintsForChimera(linear_solver)
        return builder_and_solver
