from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.ChimeraApplication as KratosChimera
# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.TrilinosApplication as KratosTrilinos     # MPI solvers
from KratosMultiphysics.ChimeraApplication import TrilinosExtension as TrilinosChimera
from KratosMultiphysics.FluidDynamicsApplication.TrilinosExtension import *
from KratosMultiphysics.FluidDynamicsApplication import TrilinosExtension as TrilinosFluid
from KratosMultiphysics.ChimeraApplication import chimera_setup_utils

# Import base class file
from KratosMultiphysics.FluidDynamicsApplication.trilinos_navier_stokes_solver_fractionalstep import TrilinosNavierStokesSolverFractionalStep

def CreateSolver(main_model_part, custom_settings):
    return TrilinosChimeraNavierStokesSolverFractionalStep(main_model_part, custom_settings)

class TrilinosChimeraNavierStokesSolverFractionalStep(TrilinosNavierStokesSolverFractionalStep):
    def __init__(self, model, custom_settings):
        # self.chimera_settings = custom_settings["chimera_settings"].Clone()
        # custom_settings.RemoveValue("chimera_settings")
        [self.chimera_settings, self.chimera_internal_parts, custom_settings] = chimera_setup_utils.SeparateAndValidateChimeraSettings(custom_settings)

        super(TrilinosChimeraNavierStokesSolverFractionalStep,self).__init__(model,custom_settings)
        KratosMultiphysics.Logger.PrintInfo("TrilinosChimeraNavierStokesSolverFractionalStep", "Construction of NavierStokesSolverMonolithic finished.")

    def AddVariables(self):
        super(TrilinosChimeraNavierStokesSolverFractionalStep,self).AddVariables()
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosChimera.CHIMERA_DISTANCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosChimera.ROTATIONAL_ANGLE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosChimera.ROTATIONAL_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosChimera.ROTATION_MESH_DISPLACEMENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosChimera.ROTATION_MESH_VELOCITY)

        KratosMultiphysics.Logger.PrintInfo("TrilinosChimeraNavierStokesSolverFractionalStep", "Fluid solver variables added correctly.")

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
            super(TrilinosChimeraNavierStokesSolverFractionalStep,self).ImportModelPart()

    def PrepareModelPart(self):
        if( not self.settings["model_import_settings"]["input_type"].GetString() == "chimera"):
            # Call the base solver to do the PrepareModelPart
            # Note that his also calls the PrepareModelPart of the turbulence model
            super(TrilinosChimeraNavierStokesSolverFractionalStep, self).PrepareModelPart()
        else :
            super(TrilinosNavierStokesSolverMonolithic,self).PrepareModelPart()


    def Initialize(self):
        # Call the base solver to create the solution strategy
        super(TrilinosChimeraNavierStokesSolverFractionalStep,self).Initialize()

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
        super(TrilinosChimeraNavierStokesSolverFractionalStep,self).InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        super(TrilinosChimeraNavierStokesSolverFractionalStep,self).FinalizeSolutionStep()
        ## Depending on the setting this will clear the created constraints
        self.chimera_process.ExecuteFinalizeSolutionStep()

    def _CreateBuilderAndSolver(self):
        pass

    def _CreateSolutionStrategy(self):
        computing_model_part = self.GetComputingModelPart()
        epetra_communicator = self._GetEpetraCommunicator()
        domain_size = computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

        # Create the pressure and velocity linear solvers
        # Note that linear_solvers is a tuple. The first item is the pressure
        # linear solver. The second item is the velocity linear solver.
        linear_solvers = self._GetLinearSolver()

        # Create the fractional step settings instance
        # TODO: next part would be much cleaner if we passed directly the parameters to the c++
        if self.settings["consider_periodic_conditions"].GetBool():
            raise Exception("Periodic conditions are not implemented in Trilinos FractionalStep.")
        else:
            fractional_step_settings = TrilinosChimera.TrilinosFractionalStepSettings(
                epetra_communicator,
                computing_model_part,
                domain_size,
                self.settings["time_order"].GetInt(),
                self.settings["use_slip_conditions"].GetBool(),
                self.settings["move_mesh_flag"].GetBool(),
                self.settings["reform_dofs_at_each_step"].GetBool())

        # Set the strategy echo level
        fractional_step_settings.SetEchoLevel(self.settings["echo_level"].GetInt())

        # Set the velocity and pressure fractional step strategy settings
        fractional_step_settings.SetStrategy(TrilinosFluid.TrilinosStrategyLabel.Pressure,
            linear_solvers[0],
            self.settings["pressure_tolerance"].GetDouble(),
            self.settings["maximum_pressure_iterations"].GetInt())

        fractional_step_settings.SetStrategy(TrilinosFluid.TrilinosStrategyLabel.Velocity,
            linear_solvers[1],
            self.settings["velocity_tolerance"].GetDouble(),
            self.settings["maximum_velocity_iterations"].GetInt())

        # Create the fractional step strategy
        if self.settings["consider_periodic_conditions"].GetBool():
            raise Exception("Periodic conditions are not implemented in Trilinos FractionalStep.")
        else:
            solution_strategy = TrilinosFluid.TrilinosFractionalStepStrategy(
                computing_model_part,
                fractional_step_settings,
                self.settings["predictor_corrector"].GetBool(),
                self.settings["compute_reactions"].GetBool())

        return solution_strategy