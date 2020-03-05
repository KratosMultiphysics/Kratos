from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.ChimeraApplication as KratosChimera
from KratosMultiphysics.ChimeraApplication import chimera_setup_utils
# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

# Import base class file
from KratosMultiphysics.FluidDynamicsApplication.navier_stokes_solver_vmsmonolithic import NavierStokesSolverMonolithic

def CreateSolver(main_model_part, custom_settings):
    return NavierStokesSolverMonolithicChimera(main_model_part, custom_settings)

class NavierStokesSolverMonolithicChimera(NavierStokesSolverMonolithic):
    def __init__(self, model, custom_settings):
        [self.chimera_settings, self.chimera_internal_parts, custom_settings] = chimera_setup_utils.SeparateAndValidateChimeraSettings(custom_settings)
        super(NavierStokesSolverMonolithicChimera,self).__init__(model,custom_settings)
        KratosMultiphysics.Logger.PrintInfo("NavierStokesSolverMonolithicChimera", "Construction of NavierStokesSolverMonolithic finished.")

    def AddVariables(self):
        super(NavierStokesSolverMonolithicChimera,self).AddVariables()
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosChimera.CHIMERA_DISTANCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosChimera.ROTATIONAL_ANGLE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosChimera.ROTATIONAL_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosChimera.ROTATION_MESH_DISPLACEMENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosChimera.ROTATION_MESH_VELOCITY)

        KratosMultiphysics.Logger.PrintInfo("NavierStokesSolverMonolithicChimera", "Fluid solver variables added correctly.")

    def ImportModelPart(self):
        if(self.settings["model_import_settings"]["input_type"].GetString() == "chimera"):
            chimera_mp_import_settings = []
            for chimera_part_levels in self.chimera_settings["chimera_parts"]:
                for chimera_part_settings in chimera_part_levels:
                    chimera_mp_import_settings.append( chimera_part_settings["model_import_settings"].Clone() )

            material_file_name = self.settings["material_import_settings"]["materials_filename"].GetString()
            import KratosMultiphysics.ChimeraApplication.chimera_modelpart_import as chim_mp_imp
            chim_mp_imp.ImportChimeraModelparts(self.main_model_part, chimera_mp_import_settings, material_file=material_file_name, parallel_type="OpenMP")
            KratosMultiphysics.Logger.PrintInfo("NavierStokesSolverMonolithicChimera", " Import of all chimera modelparts completed.")
        else:# we can use the default implementation in the base class
            super(NavierStokesSolverMonolithicChimera,self).ImportModelPart()

    def Initialize(self):
        self.chimera_process = chimera_setup_utils.GetApplyChimeraProcess(self.model, self.chimera_settings, self.settings)
        self.computing_model_part = self.main_model_part
        # If needed, create the estimate time step utility
        if (self.settings["time_stepping"]["automatic_time_step"].GetBool()):
            self.EstimateDeltaTimeUtility = self._get_automatic_time_stepping_utility()

        # Creating the solution strategy
        self.conv_criteria = KratosCFD.VelPrCriteria(self.settings["relative_velocity_tolerance"].GetDouble(),
                                                     self.settings["absolute_velocity_tolerance"].GetDouble(),
                                                     self.settings["relative_pressure_tolerance"].GetDouble(),
                                                     self.settings["absolute_pressure_tolerance"].GetDouble())

        (self.conv_criteria).SetEchoLevel(self.settings["echo_level"].GetInt())

       # Creating the time integration scheme
        if (self.element_integrates_in_time):
            # "Fake" scheme for those cases in where the element manages the time integration
            # It is required to perform the nodal update once the current time step is solved
            self.time_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticSchemeSlip(
                self.computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE],
                self.computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]+1)
            # In case the BDF2 scheme is used inside the element, the BDF time discretization utility is required to update the BDF coefficients
            if (self.settings["time_scheme"].GetString() == "bdf2"):
                time_order = 2
                self.time_discretization = KratosMultiphysics.TimeDiscretization.BDF(time_order)
            else:
                err_msg = "Requested elemental time scheme \"" + self.settings["time_scheme"].GetString()+ "\" is not available.\n"
                err_msg += "Available options are: \"bdf2\""
                raise Exception(err_msg)
        else:
            if not hasattr(self, "_turbulence_model_solver"):
                # Bossak time integration scheme
                if self.settings["time_scheme"].GetString() == "bossak":
                    if self.settings["consider_periodic_conditions"].GetBool() == True:
                        self.time_scheme = KratosCFD.ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent(
                            self.settings["alpha"].GetDouble(),
                            self.computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE],
                            KratosCFD.PATCH_INDEX)
                    else:
                        self.time_scheme = KratosCFD.ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent(
                            self.settings["alpha"].GetDouble(),
                            self.settings["move_mesh_strategy"].GetInt(),
                            self.computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])
                # BDF2 time integration scheme
                elif self.settings["time_scheme"].GetString() == "bdf2":
                    self.time_scheme = KratosCFD.GearScheme()
                # Time scheme for steady state fluid solver
                elif self.settings["time_scheme"].GetString() == "steady":
                    self.time_scheme = KratosCFD.ResidualBasedSimpleSteadyScheme(
                            self.settings["velocity_relaxation"].GetDouble(),
                            self.settings["pressure_relaxation"].GetDouble(),
                            self.computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])
                else:
                    err_msg = "Requested time scheme " + self.settings["time_scheme"].GetString() + " is not available.\n"
                    err_msg += "Available options are: \"bossak\", \"bdf2\" and \"steady\""
                    raise Exception(err_msg)
            else:
                KratosMultiphysics.Logger.PrintInfo("NavierStokesSolverMonolithicForChimera turbulent solver is not possible.")
                raise NotImplementedError

        if self.settings["consider_periodic_conditions"].GetBool() == True:
            KratosMultiphysics.Logger.PrintInfo("NavierStokesSolverMonolithicForChimera Periodic conditions are not implemented in this case .")
            raise NotImplementedError
        else:
            builder_and_solver = KratosChimera.ResidualBasedBlockBuilderAndSolverWithConstraintsForChimera(self.linear_solver)

        self.solver = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(self.computing_model_part,
                                                                            self.time_scheme,
                                                                            self.linear_solver,
                                                                            self.conv_criteria,
                                                                            builder_and_solver,
                                                                            self.settings["maximum_iterations"].GetInt(),
                                                                            self.settings["compute_reactions"].GetBool(),
                                                                            self.settings["reform_dofs_at_each_step"].GetBool(),
                                                                            self.settings["move_mesh_flag"].GetBool())

        (self.solver).SetEchoLevel(self.settings["echo_level"].GetInt())

        self.formulation.SetProcessInfo(self.computing_model_part)

        (self.solver).Initialize()

        self.solver.Check()

        KratosMultiphysics.Logger.PrintInfo("NavierStokesSolverMonolithicChimera", "Solver initialization finished.")

        chimera_setup_utils.SetChimeraInternalPartsFlag(self.model, self.chimera_internal_parts)


    def InitializeSolutionStep(self):
        self.chimera_process.ExecuteInitializeSolutionStep()
        super(NavierStokesSolverMonolithicChimera,self).InitializeSolutionStep()


    def FinalizeSolutionStep(self):
        super(NavierStokesSolverMonolithicChimera,self).FinalizeSolutionStep()
        ## Depending on the setting this will clear the created constraints
        self.chimera_process.ExecuteFinalizeSolutionStep()
