from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("FluidDynamicsApplication")

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

# Import base class file
from fluid_solver import FluidSolver

from numpy import log

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
                "restart" : {
                    "consider_restart"              : false,
                    "considered_time"               : 1.0,
                    "initial_restart_time"          : 0.0,
                    "end_time"                      : 20
                },
                "acceleration" : {
                    "type"                          : "exponential",
                    "exponential_factor"            : 1.05,
                    "log_factor"                    : 1.0,
                    "minimum_delta_time"            : 0.05,
                    "maximum_delta_time"            : 1.0,
                    "minimum_CFL"                   : 1.0,
                    "maximum_CFL"                   : 40.0,
                    "start_acceleration_time"       : 1.0,
                    "end_acceleration_time"         : 10,
                    "maximum_iteration_number"      : 8,
                    "minimum_iteration_number"      : 3
                }
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

        self.model = model

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
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.BODY_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_H)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ADVPROJ)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DIVPROJ)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_WATER_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.Y_WALL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.EXTERNAL_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.Q_VALUE)

        KratosMultiphysics.Logger.PrintInfo("NavierStokesTimeAveragedMonolithicSolver", "Fluid solver variables added correctly.")


    def AddDofs(self):
        KratosMultiphysics.VariableUtils().AddDof(KratosCFD.TIME_AVERAGED_VELOCITY_X, KratosMultiphysics.REACTION_X,  self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosCFD.TIME_AVERAGED_VELOCITY_Y, KratosMultiphysics.REACTION_Y, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosCFD.TIME_AVERAGED_VELOCITY_Z, KratosMultiphysics.REACTION_Z, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosCFD.TIME_AVERAGED_PRESSURE, KratosMultiphysics.REACTION_WATER_PRESSURE,  self.main_model_part)

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
            

    def Initialize(self):

        self.computing_model_part = self.GetComputingModelPart()

        # Creating the solution strategy
        self.conv_criteria = KratosCFD.VelPrCriteria(self.settings["relative_velocity_tolerance"].GetDouble(),
                                                     self.settings["absolute_velocity_tolerance"].GetDouble(),
                                                     self.settings["relative_pressure_tolerance"].GetDouble(),
                                                     self.settings["absolute_pressure_tolerance"].GetDouble())

        (self.conv_criteria).SetEchoLevel(self.settings["echo_level"].GetInt())

        # TODO -> move bdf coefficients calculation to ComputeBDFCoefficientsProcess
        # ComputeBDFCoefficientsProcess is only suitable for constant time steps -> bdf coefficients are now calculated within the
        # time_averaged_navier_stokes element 
        # self.bdf_process = KratosMultiphysics.ComputeBDFCoefficientsProcess(self.computing_model_part,
        #                                                                     self.settings["time_order"].GetInt())

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

        #(self.solver).Initialize() # Initialize the solver. Otherwise the constitutive law is not initializated.
        self._set_constitutive_law()
        #(self.solver).Check()
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DYNAMIC_TAU, self.settings["dynamic_tau"].GetDouble())

        # Compute the fluid domain NODAL_AREA values (required as weight for steady state estimation)
        KratosMultiphysics.CalculateNodalAreaProcess(self.main_model_part, self.domain_size).Execute()
        
        self.accelerated = self.settings["time_averaging_acceleration"]["acceleration"]["consider_acceleration"].GetBool()
        if self.accelerated == True:
            self.delta_time_acceleration_type = self.settings["time_averaging_acceleration"]["acceleration"]["type"].GetString()
            # parameters for dt accleration
            self.start_acceleration_time = self.settings["time_averaging_acceleration"]["acceleration"]["start_acceleration_time"].GetDouble()
            self.end_acceleration_time = self.settings["time_averaging_acceleration"]["acceleration"]["end_acceleration_time"].GetDouble()
            self.min_dt = self.settings["time_averaging_acceleration"]["acceleration"]["minimum_delta_time"].GetDouble()
            self.max_dt = self.settings["time_averaging_acceleration"]["acceleration"]["maximum_delta_time"].GetDouble()
            self.CFL_min =  self.settings["time_averaging_acceleration"]["acceleration"]["minimum_CFL"].GetDouble()        
            self.max_it_number =  self.settings["time_averaging_acceleration"]["acceleration"]["maximum_iteration_number"].GetDouble()
            self.min_it_number =  self.settings["time_averaging_acceleration"]["acceleration"]["minimum_iteration_number"].GetDouble()

            if self.delta_time_acceleration_type == "exponential":
                self.delta_time_acceleration_exp_factor = self.settings["time_averaging_acceleration"]["acceleration"]["exponential_factor"].GetDouble()
            if self.delta_time_acceleration_type == "log":
                self.delta_time_acceleration_log_factor = self.settings["time_averaging_acceleration"]["acceleration"]["log_factor"].GetDouble()
            if self.delta_time_acceleration_type == "linear":
                self.CFL_max = self.settings["time_averaging_acceleration"]["acceleration"]["maximum_CFL"].GetDouble()
        else:
            if self.settings["time_stepping"]["automatic_time_step"].GetBool() == True:
                self.automatic_time_step = True
                self.CFL_min =  self.settings["time_stepping"]["CFL_number"].GetDouble()
                self.min_dt = self.settings["time_stepping"]["minimum_delta_time"].GetDouble()
                self.max_dt =  self.settings["time_stepping"]["maximum_delta_time"].GetDouble()
                if(self.domain_size == 3):
                    self.EstimateDeltaTimeUtility = KratosCFD.EstimateDtUtility3D(self.computing_model_part, self.CFL_min, self.min_dt, self.max_dt, True)
                elif(self.domain_size == 2):
                    self.EstimateDeltaTimeUtility = KratosCFD.EstimateDtUtility2D(self.computing_model_part, self.CFL_min, self.min_dt, self.max_dt, True)
            else:
                self.automatic_time_step = False
                self.CFL_min = 1.0
                self.min_dt = self.settings["time_stepping"]["time_step"].GetDouble()
            self.delta_time_acceleration_type = "None"

        self.restart = self.settings["time_averaging_acceleration"]["restart"]["consider_restart"].GetBool()
        self.end_time = self.settings["time_averaging_acceleration"]["restart"]["end_time"].GetDouble()
        
        if self.restart == True:
            # parameters for time averaging acceleration
            self.initial_averaging_time_length = self.settings["time_averaging_acceleration"]["restart"]["considered_time"].GetDouble()
            self.averaging_time_length = self.initial_averaging_time_length
            self.restart_time = self.initial_averaging_time_length 
            self.initial_restart = self.settings["time_averaging_acceleration"]["restart"]["initial_restart_time"].GetDouble()
            # initial ratio
            self.initial_ratio = self.initial_averaging_time_length/self.min_dt
            self.restarted = False
        else:
            self.averaging_time_length = self.end_time
        # sets initial averaging time length
        self._set_averaging_time_length(self.averaging_time_length)
        # Initializing STEP
        self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] = 0
        self.main_model_part.ProcessInfo[KratosMultiphysics.TIME] = 0
        self.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME] = self.min_dt


        KratosMultiphysics.Logger.PrintInfo("NavierStokesTimeAveragedMonolithicSolver", "Solver initialization finished.")
        KratosMultiphysics.Logger.PrintInfo("NavierStokesTimeAveragedMonolithicSolver", "Using", self.delta_time_acceleration_type, "acceleration")


    def AdvanceInTime(self, previous_time):   
        # check if already converged, when yes, stop progressing in time
        self._check_steady_state()
        old_dt = self.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
        
        step = self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] + 1
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.STEP, step)

        if step < self.min_buffer_size:
            new_dt = self.min_dt
        else:
            if self.accelerated and previous_time > self.start_acceleration_time:
                if previous_time <= self.end_acceleration_time:
                    if self.delta_time_acceleration_type == "exponential": 
                        print("old dt:", old_dt )
                        new_dt = self._compute_increased_exponential_delta_time(previous_time, old_dt)
                    elif self.delta_time_acceleration_type == "linear":
                        new_dt = self._compute_increased_delta_time_with_CFL(previous_time)
                    elif self.delta_time_acceleration_type == "log":    
                        new_dt = self._compute_increased_logarithmic_delta_time(previous_time, step)
                    else:
                        new_dt = self.EstimateDeltaTimeUtility.EstimateDt()
        
                    if new_dt > self.max_dt:
                            new_dt = self.max_dt
                    if new_dt < self.min_dt:
                        new_dt = self.min_dt

                else:
                    old_dt = self.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
                    new_dt = old_dt
            else:
                if self.automatic_time_step:
                    new_dt = self.EstimateDeltaTimeUtility.EstimateDt()
                else:
                    new_dt = self.min_dt
        print(new_dt)
        # clone the time steps to save place for new variables new_time t and new_dt dt
        self.main_model_part.CloneTimeStep(previous_time)
        # tn+1 = tn + dtn
        new_time = previous_time + old_dt
        #print("new_time ", new_time)
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, new_time)
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, new_dt)

        if (self.restart == True):
            self._compute_averaging_time_length(new_time, new_dt)
            self._set_averaging_time_length(self.averaging_time_length)
        
        print("New Dt: ", new_dt)
        return new_time


    def InitializeSolutionStep(self):
        if self._TimeBufferIsInitialized():
            #(self.bdf_process).Execute()
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


    def _compute_initial_ratio(self):
        if(self.domain_size == 3):
            EstimateDeltaTimeUtility = KratosCFD.EstimateDtUtility3D(self.computing_model_part, self.CFL_min, self.min_dt, self.max_dt, True)
        elif(self.domain_size == 2):
            EstimateDeltaTimeUtility = KratosCFD.EstimateDtUtility2D(self.computing_model_part, self.CFL_min, self.min_dt, self.max_dt, True)      
        self.CFL_min_dt = EstimateDeltaTimeUtility.EstimateDt()
        self.initial_ratio = self.initial_averaging_time_length/self.CFL_min_dt


    def _compute_increased_exponential_delta_time(self, current_time, old_dt):

        print("Nr Nonlinear Iterations: " + str(self.main_model_part.ProcessInfo[KratosMultiphysics.NL_ITERATION_NUMBER]))
        # if few iteration steps needed -> increase time step
        if self.main_model_part.ProcessInfo[KratosMultiphysics.NL_ITERATION_NUMBER] <= self.min_it_number:
            new_dt = old_dt*self.delta_time_acceleration_exp_factor 
        # if a lot of iteration steps needed -> increase time step
        elif self.main_model_part.ProcessInfo[KratosMultiphysics.NL_ITERATION_NUMBER] >= self.max_it_number:
            new_dt = old_dt/self.delta_time_acceleration_exp_factor 
        else:
            new_dt = old_dt

        return new_dt


    def _compute_increased_logarithmic_delta_time(self, current_time, step):
        step -= self.step_before_accleration
        step += 4
        new_dt = self.delta_time_acceleration_log_factor*self.min_dt*(log(step+1)-log(step))
        if new_dt > self.max_dt:
            new_dt = self.max_dt

        return new_dt


    def _compute_increased_delta_time_with_CFL(self, current_time):

        if ( current_time > self.start_acceleration_time): 
            if (current_time <= self.end_acceleration_time):
                CFL = self.CFL_min + current_time / self.end_acceleration_time * (self.CFL_max - self.CFL_min)
            else:
                CFL = self.CFL_max
        else: CFL = self.CFL_min
        
        if(self.domain_size == 3):
            EstimateDeltaTimeUtility = KratosCFD.EstimateDtUtility3D(self.computing_model_part, CFL, self.min_dt, self.max_dt, True)
        elif(self.domain_size == 2):
            EstimateDeltaTimeUtility = KratosCFD.EstimateDtUtility2D(self.computing_model_part, CFL, self.min_dt, self.max_dt, True)
        
        new_dt = EstimateDeltaTimeUtility.EstimateDt()

        print("New dt is: ", new_dt)
        return new_dt


    def _compute_averaging_time_length(self, current_time, dt):
        if self.restart == True:
            if ( current_time > self.restart_time):
                # RESTART -> reset averaging time length to the initial averaging time length
                self.averaging_time_length = self.initial_averaging_time_length
                self.restart_time += self.averaging_time_length
                self._restart()
                print("### RESTART at " + str(current_time) + " ###")
                print("New Restart Time: " + str(self.restart_time))
            else: 
                self.averaging_time_length += dt
                if (self.restarted == False):
                    self._initial_restart(current_time) 
        else:
            self.averaging_time_length += dt
            print("Averaging time length set to " + str(self.averaging_time_length) )           
        print("=================================================================")



    def _compute_averaging_time_length_with_const_ratio(self, current_time, dt):
        if ( current_time > self.restart_time):
            if (self.restarted == False):
                # Initializing CFL_min_dt
                self._compute_initial_ratio()
                self.restarted = True
            # RESTART -> reset averaging time length to the initial averaging time length
            self.averaging_time_length = self.initial_ratio * dt
            if (self.averaging_time_length > current_time):
                self.averaging_time_length = current_time - dt
            self.restart_time += self.averaging_time_length
            self._restart()
            print("### RESTART at " + str(current_time) + " ###")
            print("New Restart Time: " + str(self.restart_time))
        else: 
            self.averaging_time_length += dt
            if (self.restarted == False):
                self._initial_restart(current_time) 
        print("Averaging time length set to ", self.averaging_time_length)           
        print("=================================================================")


    def _initial_restart(self, current_time):
        if current_time > self.initial_restart:
            print("INITIAL RESTART: DROPPING PREVIOUS INFORMATION")
            self._restart()
            self.averaging_time_length = 0.0
            self.restarted = True


    def _restart(self):
        # At restart, all the previous averaged velocity and pressure should be set to the current one
        self.main_model_part.CloneSolutionStep()
        self.main_model_part.CloneSolutionStep()
        self.main_model_part.CloneSolutionStep()

    
    def _compute_averaging_time_length_with_ghost_time(self, current_time, dt):
        if ( current_time > self.restart_time):
            # RESTART -> reset averaging time length to the initial averaging time length
            self.averaging_time_length = 2*self.restart_time 
            self.restart_time += self.initial_averaging_time_length
            self._restart()
            print("### RESTART at " + str(current_time) + " ###")
            print("New Restart Time: " + str(self.restart_time))
        else:
            if (self.restarted == False):
                self._initial_restart(current_time)
            self.averaging_time_length += dt
        print("Averaging time length set to ", self.averaging_time_length)           
        print("=================================================================")


    def _set_averaging_time_length(self, new_averaging_time_length=0.0):
        if new_averaging_time_length==0.0:
            averaging_time_length = self.settings["time_averaging_acceleration"]["restart"]["considered_time"].GetDouble()
        else:
            averaging_time_length = new_averaging_time_length
        self.main_model_part.ProcessInfo.SetValue(KratosCFD.AVERAGING_TIME_LENGTH, averaging_time_length)


    def _check_steady_state(self):
        SteadyStateIndicatorUtility = KratosCFD.SteadyStateIndicatorUtility(self.main_model_part)
        SteadyStateIndicatorUtility.EstimateTimeAveragedQuantityChangesInTime()
        change_in_velocity = SteadyStateIndicatorUtility.GetTimeAveragedVelocityChange()
        change_in_pressure = SteadyStateIndicatorUtility.GetTimeAveragedPressureChange()
        print("Change in velocity in percentage: " + str(change_in_velocity))
        print("Change in pressure in percentage: " + str(change_in_pressure))


    def _set_constitutive_law(self):
        ## Construct the constitutive law needed for the embedded element
        if(self.domain_size == 3):
            self.main_model_part.Properties[1][KratosMultiphysics.CONSTITUTIVE_LAW] = KratosCFD.Newtonian3DLaw()
        elif(self.domain_size == 2):
            self.main_model_part.Properties[1][KratosMultiphysics.CONSTITUTIVE_LAW] = KratosCFD.Newtonian2DLaw()

