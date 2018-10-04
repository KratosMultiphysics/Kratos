from __future__ import absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.mpi as KratosMPI                          # MPI-python interface

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("FluidDynamicsApplication","MetisApplication","TrilinosApplication")

# Import applications
import KratosMultiphysics.MetisApplication as KratosMetis           # Partitioning
import KratosMultiphysics.TrilinosApplication as KratosTrilinos     # MPI solvers
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid   # Fluid dynamics application

# Import base class file
import navier_stokes_solver_fractionalstep
import trilinos_import_model_part_utility

def CreateSolver(model, custom_settings):
    return TrilinosNavierStokesSolverFractionalStep(model, custom_settings)

class TrilinosNavierStokesSolverFractionalStep(navier_stokes_solver_fractionalstep.NavierStokesSolverFractionalStep):

    def _ValidateSettings(self,settings):
        ## Default settings string in Json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "FractionalStep",
            "model_part_name": "",
            "domain_size": -1,
            "model_import_settings": {
                    "input_type": "mdpa",
                    "input_filename": "unknown_name"
            },
            "predictor_corrector": false,
            "maximum_velocity_iterations": 3,
            "maximum_pressure_iterations": 3,
            "velocity_tolerance": 1e-3,
            "pressure_tolerance": 1e-2,
            "dynamic_tau": 0.01,
            "oss_switch": 0,
            "echo_level": 0,
            "consider_periodic_conditions": false,
            "time_order": 2,
            "compute_reactions": false,
            "reform_dofs_at_each_step": false,
            "pressure_linear_solver_settings": {
                "solver_type"                        : "MultiLevelSolver",
                "max_iteration"                      : 200,
                "tolerance"                          : 1e-6,
                "symmetric"                          : true,
                "scaling"                            : true,
                "reform_preconditioner_at_each_step" : true,
                "verbosity"                          : 0
            },
            "velocity_linear_solver_settings": {
                "solver_type"                        : "MultiLevelSolver",
                "max_iteration"                      : 200,
                "tolerance"                          : 1e-6,
                "symmetric"                          : false,
                "scaling"                            : true,
                "reform_preconditioner_at_each_step" : true,
                "verbosity"                          : 0
            },
            "volume_model_part_name" : "volume_model_part",
            "skin_parts":[""],
            "no_skin_parts":[""],
            "time_stepping": {
                "automatic_time_step" : false,
                "CFL_number"          : 1,
                "minimum_delta_time"  : 1e-4,
                "maximum_delta_time"  : 0.01
            },
            "move_mesh_flag": false,
            "use_slip_conditions": true
        }""")

        settings.ValidateAndAssignDefaults(default_settings)
        return settings

    def __init__(self, model, custom_settings):
        # Note: deliberately calling the constructor of the base python solver (the parent of my parent)
        super(navier_stokes_solver_fractionalstep.NavierStokesSolverFractionalStep,self).__init__(model,custom_settings)

        self.element_name = "FractionalStep"
        self.condition_name = "WallCondition"
        self.min_buffer_size = 3

        self._is_printing_rank = (KratosMPI.mpi.rank == 0)

        ## Construct the linear solvers
        import trilinos_linear_solver_factory
        self.pressure_linear_solver = trilinos_linear_solver_factory.ConstructSolver(self.settings["pressure_linear_solver_settings"])
        self.velocity_linear_solver = trilinos_linear_solver_factory.ConstructSolver(self.settings["velocity_linear_solver_settings"])

        self.compute_reactions = self.settings["compute_reactions"].GetBool()

        if self._IsPrintingRank():
            #TODO: CHANGE THIS ONCE THE MPI LOGGER IS IMPLEMENTED
            KratosMultiphysics.Logger.PrintInfo("TrilinosNavierStokesSolverFractionalStep","Construction of TrilinosNavierStokesSolverFractionalStep solver finished.")


    def AddVariables(self):
        ## Add variables from the base class
        super(TrilinosNavierStokesSolverFractionalStep, self).AddVariables()

        ## Add specific MPI variables
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)
        KratosMPI.mpi.world.barrier()

        if self._IsPrintingRank():
            #TODO: CHANGE THIS ONCE THE MPI LOGGER IS IMPLEMENTED
            KratosMultiphysics.Logger.PrintInfo("TrilinosNavierStokesSolverFractionalStep","variables for the trilinos fractional step solver added correctly")


    def ImportModelPart(self):

        ## Construct the Trilinos import model part utility
        self.trilinos_model_part_importer = trilinos_import_model_part_utility.TrilinosImportModelPartUtility(self.main_model_part, self.settings)
        ## Execute the Metis partitioning and reading
        self.trilinos_model_part_importer.ImportModelPart()

        if self._IsPrintingRank():
            #TODO: CHANGE THIS ONCE THE MPI LOGGER IS IMPLEMENTED
            KratosMultiphysics.Logger.PrintInfo("TrilinosNavierStokesSolverFractionalStep","MPI model reading finished.")

    def PrepareModelPart(self):
        super(TrilinosNavierStokesSolverFractionalStep,self).PrepareModelPart()
        ## Construct Trilinos the communicators
        self.trilinos_model_part_importer.CreateCommunicators()


    def AddDofs(self):
        ## Base class DOFs addition
        super(TrilinosNavierStokesSolverFractionalStep, self).AddDofs()
        KratosMPI.mpi.world.barrier()

        if self._IsPrintingRank():
            #TODO: CHANGE THIS ONCE THE MPI LOGGER IS IMPLEMENTED
            KratosMultiphysics.Logger.PrintInfo("TrilinosNavierStokesSolverFractionalStep","DOFs for the VMS Trilinos fluid solver added correctly in all processors.")


    def Initialize(self):
        ## Construct the communicator
        self.EpetraComm = KratosTrilinos.CreateCommunicator()

        ## Get the computing model part
        self.computing_model_part = self.GetComputingModelPart()

        ## If needed, create the estimate time step utility
        if (self.settings["time_stepping"]["automatic_time_step"].GetBool()):
            self.EstimateDeltaTimeUtility = self._GetAutomaticTimeSteppingUtility()

        #TODO: next part would be much cleaner if we passed directly the parameters to the c++
        if self.settings["consider_periodic_conditions"] == True:
            self.solver_settings = KratosTrilinos.TrilinosFractionalStepSettingsPeriodic(self.EpetraComm,
                                                                                        self.computing_model_part,
                                                                                        self.computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE],
                                                                                        self.settings["time_order"].GetInt(),
                                                                                        self.settings["use_slip_conditions"].GetBool(),
                                                                                        self.settings["move_mesh_flag"].GetBool(),
                                                                                        self.settings["reform_dofs_at_each_step]"].GetBool(),
                                                                                        KratosFluid.PATCH_INDEX)

        else:
            self.solver_settings = KratosTrilinos.TrilinosFractionalStepSettings(self.EpetraComm,
                                                                                 self.computing_model_part,
                                                                                 self.computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE],
                                                                                 self.settings["time_order"].GetInt(),
                                                                                 self.settings["use_slip_conditions"].GetBool(),
                                                                                 self.settings["move_mesh_flag"].GetBool(),
                                                                                 self.settings["reform_dofs_at_each_step"].GetBool())

        self.solver_settings.SetEchoLevel(self.settings["echo_level"].GetInt())

        self.solver_settings.SetStrategy(KratosTrilinos.TrilinosStrategyLabel.Velocity,
                                         self.velocity_linear_solver,
                                         self.settings["velocity_tolerance"].GetDouble(),
                                         self.settings["maximum_velocity_iterations"].GetInt())

        self.solver_settings.SetStrategy(KratosTrilinos.TrilinosStrategyLabel.Pressure,
                                         self.pressure_linear_solver,
                                         self.settings["pressure_tolerance"].GetDouble(),
                                         self.settings["maximum_pressure_iterations"].GetInt())

        self.solver = KratosTrilinos.TrilinosFSStrategy(self.computing_model_part,
                                                        self.solver_settings,
                                                        self.settings["predictor_corrector"].GetBool(),
                                                        KratosFluid.PATCH_INDEX)

        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DYNAMIC_TAU, self.settings["dynamic_tau"].GetDouble())
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.OSS_SWITCH, self.settings["oss_switch"].GetInt())

        (self.solver).Initialize()
        (self.solver).Check()

        if self._IsPrintingRank():
            #TODO: CHANGE THIS ONCE THE MPI LOGGER IS IMPLEMENTED
            KratosMultiphysics.Logger.PrintInfo("TrilinosNavierStokesSolverFractionalStep","Initialization TrilinosNavierStokesSolverFractionalStep finished")

