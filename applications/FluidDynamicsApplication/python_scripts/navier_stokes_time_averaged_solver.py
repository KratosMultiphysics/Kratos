from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("FluidDynamicsApplication")

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

# Import base class file
from fluid_solver import FluidSolver

def CreateSolver(model, custom_settings):
    return NavierStokesTimeAveragedMonolithicSolver(model, custom_settings)

class NavierStokesTimeAveragedMonolithicSolver(FluidSolver):

    def _ValidateSettings(self, settings):
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "TimeAveraged",
            "model_part_name": "MainModelPart",
            "domain_size": -1,
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name"
            },
            "maximum_iterations": 7,
            "predictor_corrector": true,
            "dynamic_tau": 1.0,
            "echo_level": 0,
            "time_order": 2,
            "compute_reactions": false,
            "reform_dofs_at_each_step": false,
            "relative_velocity_tolerance": 1e-3,
            "absolute_velocity_tolerance": 1e-5,
            "relative_pressure_tolerance": 1e-3,
            "absolute_pressure_tolerance": 1e-5,
            "linear_solver_settings"       : {
                "solver_type"         : "BICGSTABSolver",
                "max_iteration"       : 5000,
                "tolerance"           : 1e-7,
                "preconditioner_type" : "DiagonalPreconditioner",
                "scaling"             : false
            },
            "volume_model_part_name" : "volume_model_part",
            "skin_parts": [""],
            "no_skin_parts":[""],
            "time_stepping"                : {
                "automatic_time_step" : true,
                "CFL_number"          : 1,
                "minimum_delta_time"  : 1e-2,
                "maximum_delta_time"  : 1.0
            },
            "time_averaging_acceleration"        :{
            "considered_time"           : 100.0,
            "minimum_delta_time"        : 0.5,
            "maximum_delta_time"        : 20.0,
            "minimum_CFL"               : 1.0,
            "maximum_CFL"               : 20.0,
            "start_acceleration_time"   : 100,
            "end_acceleration_time"     : 5000.0,
            "end_time"                  : 10000.0
            },
            "move_mesh_flag": false,
            "use_slip_conditions": true
        }""")

        settings.ValidateAndAssignDefaults(default_settings)
        return settings

    def __init__(self, model, custom_settings):
        super(NavierStokesTimeAveragedMonolithicSolver,self).__init__(model,custom_settings)

        self.element_name = "TimeAveragedNavierStokes"
        self.condition_name = "TimeAveragedNavierStokesWallCondition"
        self.min_buffer_size = 5
        
        # Get domain size
        self.domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

        # There is only a single rank in OpenMP, we always print
        self._is_printing_rank = True

        ## Construct the linear solver
        import linear_solver_factory
        self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

        KratosMultiphysics.Logger.PrintInfo("NavierStokesTimeAveragedMonolithicSolver", "Construction of NavierStokesTimeAveragedMonolithicSolver finished.")


    def AddVariables(self):
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY) # TODO: Remove this once the "old" embedded elements get the density from the properties (or once we delete them)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DYNAMIC_VISCOSITY) # TODO: Remove this once the "old" embedded elements get the density from the properties (or once we delete them)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.TIME_AVERAGED_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.TIME_AVERAGED_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.IS_STRUCTURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.BODY_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_H)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_WATER_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.EXTERNAL_PRESSURE)

        KratosMultiphysics.Logger.PrintInfo("NavierStokesTimeAveragedMonolithicSolver", "Fluid solver variables added correctly.")


    def AddDofs(self):
        KratosMultiphysics.VariableUtils().AddDof(KratosCFD.TIME_AVERAGED_VELOCITY_X,  self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosCFD.TIME_AVERAGED_VELOCITY_Y,  self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosCFD.TIME_AVERAGED_VELOCITY_Z,  self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosCFD.TIME_AVERAGED_PRESSURE,    self.main_model_part)

        if self._IsPrintingRank():
            KratosMultiphysics.Logger.PrintInfo("NavierStokesTimeAveragedMonolithicSolver", "Fluid solver DOFs added correctly.")


    def ImportModelPart(self):
        super(NavierStokesTimeAveragedMonolithicSolver, self).ImportModelPart()


    def PrepareModelPart(self):
        super(NavierStokesTimeAveragedMonolithicSolver, self).PrepareModelPart()
        if not self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            ## Sets DENSITY, DYNAMIC_VISCOSITY and SOUND_VELOCITY
            self._set_physical_properties()
            ## Sets the constitutive law
            self._set_constitutive_law()
            ## Sets averaging time length
            self._set_averaging_time_length()


    def Initialize(self):
        self.computing_model_part = self.GetComputingModelPart()

        # If needed, create the estimate time step utility
        if (self.settings["time_stepping"]["automatic_time_step"].GetBool()):
            self.EstimateDeltaTimeUtility = self._GetAutomaticTimeSteppingUtility()

        # Creating the solution strategy
        self.conv_criteria = KratosCFD.VelPrCriteria(self.settings["relative_velocity_tolerance"].GetDouble(),
                                                     self.settings["absolute_velocity_tolerance"].GetDouble(),
                                                     self.settings["relative_pressure_tolerance"].GetDouble(),
                                                     self.settings["absolute_pressure_tolerance"].GetDouble())

        (self.conv_criteria).SetEchoLevel(self.settings["echo_level"].GetInt())

        self.bdf_process = KratosMultiphysics.ComputeBDFCoefficientsProcess(self.computing_model_part,
                                                                            self.settings["time_order"].GetInt())

        time_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticSchemeSlip(self.domain_size,   # Domain size (2,3)
                                                                                        self.domain_size+1) # DOFs (3,4)

        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(self.linear_solver)

        self.solver = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(self.computing_model_part,
                                                                            time_scheme,
                                                                            self.linear_solver,
                                                                            self.conv_criteria,
                                                                            builder_and_solver,
                                                                            self.settings["maximum_iterations"].GetInt(),
                                                                            self.settings["compute_reactions"].GetBool(),
                                                                            self.settings["reform_dofs_at_each_step"].GetBool(),
                                                                            self.settings["move_mesh_flag"].GetBool())

        (self.solver).SetEchoLevel(self.settings["echo_level"].GetInt())

        (self.solver).Initialize() # Initialize the solver. Otherwise the constitutive law is not initializated.
        (self.solver).Check()
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DYNAMIC_TAU, self.settings["dynamic_tau"].GetDouble())

        # Compute the fluid domain NODAL_AREA values (required as weight for steady state estimation)
        KratosMultiphysics.CalculateNodalAreaProcess(self.main_model_part, 
                                                     self.domain_size).Execute()

        # parameters for sample length calculation
        self.initial_averaging_time_length = self.settings["time_averaging_acceleration"]["considered_time"].GetDouble()
        self.averaging_time_length = self.initial_averaging_time_length
        self.restart_time = self.initial_averaging_time_length 
        self.end_time = self.settings["time_averaging_acceleration"]["end_time"].GetDouble()
        # parameters for dt accleration
        self.start_acceleration_time = self.settings["time_averaging_acceleration"]["start_acceleration_time"].GetDouble()
        self.end_acceleration_time = self.settings["time_averaging_acceleration"]["end_acceleration_time"].GetDouble()
        
        KratosMultiphysics.Logger.PrintInfo("NavierStokesTimeAveragedMonolithicSolver", "Solver initialization finished.")


    def AdvanceInTime(self, current_time):
        # dt = self._ComputeDeltaTime()        
        self._check_steady_state()
        
        dt = self._compute_increased_delta_time(current_time)
        new_time = current_time + dt
        self.main_model_part.CloneTimeStep(new_time)
        self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1
        
        self._compute_averaging_time_length(new_time, dt)
        self._set_averaging_time_length(self.averaging_time_length)

        return new_time


    def Solve(self):
        self.InitializeSolutionStep()
        self.Predict()
        self.SolveSolutionStep()
        self.FinalizeSolutionStep()


    def InitializeSolutionStep(self):
        if self._TimeBufferIsInitialized():
            (self.bdf_process).Execute()
            (self.solver).InitializeSolutionStep()


    def _set_physical_properties(self):
        ## Set the SOUND_VELOCITY value (wave velocity)
        if self.main_model_part.Properties[1].Has(KratosMultiphysics.SOUND_VELOCITY):
            self.main_model_part.ProcessInfo[KratosMultiphysics.SOUND_VELOCITY] = self.main_model_part.Properties[1][KratosMultiphysics.SOUND_VELOCITY]
        else:
            # If the wave velocity is not defined take a large enough value to consider the fluid as incompressible
            default_sound_velocity = 1e+12
            self.main_model_part.ProcessInfo[KratosMultiphysics.SOUND_VELOCITY] = default_sound_velocity

        # Transfer density and (dynamic) viscostity to the nodes
        for el in self.main_model_part.Elements:
            rho = el.Properties.GetValue(KratosMultiphysics.DENSITY)
            if rho <= 0.0:
                raise Exception("DENSITY set to {0} in Properties {1}, positive number expected.".format(rho,el.Properties.Id))
            dyn_viscosity = el.Properties.GetValue(KratosMultiphysics.DYNAMIC_VISCOSITY)
            if dyn_viscosity <= 0.0:
                raise Exception("DYNAMIC_VISCOSITY set to {0} in Properties {1}, positive number expected.".format(dyn_viscosity,el.Properties.Id))
            break

        # TODO: Remove this once the "old" embedded elements get the density from the properties (or once we delete them)
        KratosMultiphysics.VariableUtils().SetScalarVar(KratosMultiphysics.DENSITY, rho, self.main_model_part.Nodes)
        KratosMultiphysics.VariableUtils().SetScalarVar(KratosMultiphysics.DYNAMIC_VISCOSITY, dyn_viscosity, self.main_model_part.Nodes)


    def _compute_increased_delta_time(self, current_time):
        dt_min = self.settings["time_averaging_acceleration"]["minimum_delta_time"].GetDouble()
        dt_max = self.settings["time_averaging_acceleration"]["maximum_delta_time"].GetDouble()
        CFL_min =  self.settings["time_averaging_acceleration"]["minimum_CFL"].GetDouble()
        CFL_max = self.settings["time_averaging_acceleration"]["maximum_CFL"].GetDouble()

        if ( current_time > self.start_acceleration_time): 
            if (current_time <= self.end_acceleration_time):
                CFL = CFL_min + current_time / self.end_acceleration_time * (CFL_max - CFL_min)
            else:
                CFL = CFL_max
        else: CFL = CFL_min
        
        if(self.domain_size == 3):
            EstimateDeltaTimeUtility = KratosCFD.EstimateDtUtility3D(self.computing_model_part, CFL, dt_min, dt_max, True)
        elif(self.domain_size == 2):
            EstimateDeltaTimeUtility = KratosCFD.EstimateDtUtility2D(self.computing_model_part, CFL, dt_min, dt_max, True)
        
        new_dt = EstimateDeltaTimeUtility.EstimateDt()

        print("New dt is: ", new_dt)
        return new_dt


    def _compute_averaging_time_length(self, current_time, dt):
        if (current_time > self.initial_averaging_time_length):
            if ( current_time > self.restart_time):
                self.averaging_time_length = self.initial_averaging_time_length
                self.restart_time += self.restart_time
            else: 
                self.averaging_time_length += dt
            print("Averaging time length set to ", self.averaging_time_length, ", droping previous time infomation")


    def _set_averaging_time_length(self, new_averaging_time_length=0.0):
        if new_averaging_time_length==0.0:
            averaging_time_length = self.settings["time_averaging_acceleration"]["considered_time"].GetDouble()
        else:
            averaging_time_length = new_averaging_time_length
        self.main_model_part.ProcessInfo.SetValue(KratosCFD.AVERAGING_TIME_LENGTH, averaging_time_length)


    def _check_steady_state(self):
        SteadyStateIndicatorUtility = KratosCFD.SteadyStateIndicatorUtility(self.computing_model_part)
        SteadyStateIndicatorUtility.EstimateQuantityChangesInTime()
        change_in_velocity = SteadyStateIndicatorUtility.GetVelocityChange()
        change_in_pressure = SteadyStateIndicatorUtility.GetPressureChange()
        print("Change in velocity in percentage: " + str(change_in_velocity))
        print("Change in pressure in percentage: " + str(change_in_pressure))


    def _set_constitutive_law(self):
        ## Construct the constitutive law needed for the embedded element
        if(self.domain_size == 3):
            self.main_model_part.Properties[1][KratosMultiphysics.CONSTITUTIVE_LAW] = KratosCFD.Newtonian3DLaw()
        elif(self.domain_size == 2):
            self.main_model_part.Properties[1][KratosMultiphysics.CONSTITUTIVE_LAW] = KratosCFD.Newtonian2DLaw()
