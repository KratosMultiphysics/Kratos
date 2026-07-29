import KratosMultiphysics
import KratosMultiphysics.MPMApplication as KratosMPM
import numpy as np


def Factory(settings: KratosMultiphysics.Parameters, model: KratosMultiphysics.Model) -> KratosMultiphysics.Process:
    if not isinstance(model, KratosMultiphysics.Model):
        raise Exception("expected input shall be a Model object")
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return MPMInletProcess(model, settings["Parameters"])

class MPMInletProcess(KratosMultiphysics.Process):
    def __init__(self, Model: KratosMultiphysics.Model, settings: KratosMultiphysics.Parameters) -> None:
        KratosMultiphysics.Process.__init__(self)

        # "modulus" must be a number
        if settings.Has("modulus") and not settings["modulus"].IsNumber():
            raise Exception('Parameter "modulus" must be a number')

        # "direction" must be a vector with three components
        if settings.Has("direction"):
            if not settings["direction"].IsVector() or settings["direction"].GetVector().Size() != 3:
                raise Exception('Parameter "direction" must be a vector of length 3')

        settings.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.inlet_model_part_name              = settings["inlet_model_part_name"].GetString()
        if settings["inlet_bottom_nodes_model_part_name"].GetString().startswith('Background_Grid'):
            self.inlet_bottom_nodes_model_part_name = settings["inlet_bottom_nodes_model_part_name"].GetString().split(".",1)[1]
        else:
            self.inlet_bottom_nodes_model_part_name = settings["inlet_bottom_nodes_model_part_name"].GetString()

        self.inlet_model_part                   = Model[self.inlet_model_part_name]
        self.inlet_bottom_nodes_model_part      = Model[settings["inlet_bottom_nodes_model_part_name"].GetString()]
        self.mpm_model_part                     = Model[settings["mpm_model_part_name"].GetString()]
        self.grid_model_part                    = Model["Background_Grid"]

        self.assigned_velocity                  = np.array(settings["velocity"].GetVector())
        self._IsVelocityNumberOrFunction()

        self.density                            = settings["density"].GetDouble()
        self.inlet_perpendicular_distance       = settings["inlet_perpendicular_distance"].GetDouble()
        self.inlet_normal_distance              = settings["inlet_normal_distance"].GetDouble()
        self.is_fixed                           = settings["is_fixed"].GetBool()
        self.gravity                            = settings["gravity"].GetVector()
        self.unit_gravity                       = self.gravity / np.linalg.norm(self.gravity)


        if settings.Has("height") and settings.Has("free_surface_elevation"):
            raise Exception('MPMInletProcess: Both height and free_surface_elevation exist in parameter. Choose one!')

        if settings.Has("height"):
            self.is_calculate_height_from_surface_elevation = False
            self.height                 = settings["height"].GetDouble()



        if settings.Has("free_surface_elevation"):
            self.is_calculate_height_from_surface_elevation = True
            self.free_surface_elevation = settings["free_surface_elevation"].GetDouble()


        self.interval = KratosMultiphysics.IntervalUtility(settings)
        self.variable_utils = KratosMultiphysics.VariableUtils()
        self.water_depth_variable = KratosMultiphysics.KratosGlobals.GetVariable("WATER_DEPTH")
        self.free_surface_elevation_variable = KratosMultiphysics.KratosGlobals.GetVariable("MPM_FREE_SURFACE_ELEVATION")

        # self.mpm_model_part.CreateSubModelPart(self.inlet_bottom_nodes_model_part_name)
        self.accumulated_displacement = np.array([0.0, 0.0, 0.0])

        self.delta_time = self.grid_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
        print("self.delta_time= ", self.delta_time)




    @staticmethod
    def GetDefaultParameters():
        return KratosMultiphysics.Parameters("""
            {
                "inlet_model_part_name"             : "please_specify_model_part_name",
                "inlet_bottom_nodes_model_part_name": "please_specify_model_part_name",
                "grid_model_part_name"              : "please_specify_model_part_name",
                "mpm_model_part_name"               : "MPM_Material",
                "height"                            : 0.0,
                "density"                           : 1.0,
                "inlet_perpendicular_distance"      : 1.0,
                "inlet_normal_distance"             : 1.0,
                "velocity"                          : [0.0, 0.0, 0.0],
                "is_fixed"                          : true,
                "gravity"                           : [0.0, -9.81, 0.0],
                "interval"                          : [0.0, 1e30]
            }
            """)
                # "free_surface_elevation"            : 0.0,

    # def ExecuteInitialize(self):

    def ExecuteBeforeSolutionLoop(self):
        # Assign height and velocity to the nodes
        # if self.value_is_numeric:
        # self.grid_model_part.AddNodalSolutionStepVariable(self.water_depth_variable)
        self.variable_utils.SetVariable(KratosMultiphysics.DISPLACEMENT, self.assigned_velocity * self.delta_time, self.inlet_model_part.Nodes)
        self.variable_utils.SetVariable(KratosMultiphysics.VELOCITY, self.assigned_velocity, self.inlet_model_part.Nodes)

        if self.is_calculate_height_from_surface_elevation:
            self.variable_utils.SetNonHistoricalVariable(self.free_surface_elevation_variable, self.free_surface_elevation, self.inlet_bottom_nodes_model_part.Nodes)
            KratosMPM.ComputeHeightFromFreeSurface(self.inlet_bottom_nodes_model_part, self.unit_gravity)

        else:
            for node in self.inlet_bottom_nodes_model_part.Nodes:
                node.SetSolutionStepValue(KratosMPM.WATER_DEPTH, self.height)
                # node.SetSolutionStepValue(self.water_depth_variable, self.height)
            # self.variable_utils.SetNonHistoricalVariable(self.water_depth_variable, self.height, self.inlet_bottom_nodes_model_part.Nodes)

        KratosMPM.GenerateMaterialPointAtInlet(self.grid_model_part, self.mpm_model_part, self.inlet_bottom_nodes_model_part, False, self.inlet_normal_distance, self.gravity, self.unit_gravity, self.inlet_perpendicular_distance)

    def ExecuteInitializeSolutionStep(self):
        current_time = self.inlet_model_part.ProcessInfo[KratosMultiphysics.TIME]
        self.delta_time = self.mpm_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]


        if self.interval.IsInInterval(current_time):
            # Assign current velocities and displacement at the active nodes of the inlet
            print("debug: is in interval")
            self.step_is_active = True


            if self.is_fixed == True:
                # self.variable_utils.ApplyFixity(KratosMultiphysics.VELOCITY_X    , self.is_fixed, self.inlet_model_part.Nodes)
                self.variable_utils.ApplyFixity(KratosMultiphysics.DISPLACEMENT_X, self.is_fixed, self.inlet_model_part.Nodes)

            # impose displacement and velocity to the nodes
            # if self.velocity_value_is_numeric:
            # self.variable_utils.SetVariable(KratosMultiphysics.VELOCITY    , self.assigned_velocity                  , self.inlet_model_part.Nodes)
            self.variable_utils.SetVariable(KratosMultiphysics.DISPLACEMENT, self.assigned_velocity * self.delta_time, self.inlet_model_part.Nodes)
            for node in self.inlet_model_part.Nodes:
                print(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT))

            # elif self.value_is_function:
            #     if self.aux_function.DependsOnSpace() == False: #depends on time only
            #         self.value = self.aux_function.CallFunction(0.0,0.0,0.0,current_time,0.0,0.0,0.0)
            #         self.variable_utils.SetVariable(self.variable, self.value, self.inlet_model_part.Nodes)
            #     else: #most general case - space varying function (possibly also time varying)
            #         self.cpp_apply_function_utility.ApplyFunction(self.variable, current_time)

            if (self._IsInletStep()):
                # Generate material points at the inlet
                KratosMPM.GenerateMaterialPointAtInlet(self.grid_model_part, self.mpm_model_part, self.inlet_bottom_nodes_model_part, False, np.linalg.norm(self.accumulated_displacement), self.gravity, self.unit_gravity, self.inlet_perpendicular_distance)
                self.accumulated_displacement = np.array([0.0, 0.0, 0.0])
            # else:
                # if not inlet step, accumulate the volume and mass to be assigned to during the inlet step


    def ExecuteFinalizeSolutionStep(self):
        """This method is executed in order to finalize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        if self.step_is_active:
            # Here we free all of the nodes in the mesh
            fixity_status  = False
            self.variable_utils.ApplyFixity(KratosMultiphysics.DISPLACEMENT_X, fixity_status, self.inlet_model_part.Nodes)

        self.step_is_active = False

    def _IsInletStep(self):
        # Note: Should only be used if velocity does not vary in space
        # Determine if this is the step of which material points will be generated at the inlet
        # It is determined by the spacing between material point given by the user
        # If the accumulated velocity * dt is bigger than the spacing, returns true
        # else, accumulate velocity * dt
        print(self.accumulated_displacement)
        print(self.assigned_velocity)
        print(self.assigned_velocity * self.delta_time)
        # TODO:the accumulated_displacement should be calculated using surface normal velocity instead!
        self.accumulated_displacement += self.assigned_velocity * self.delta_time

        print('calculating accumulated disp: ', np.linalg.norm(self.accumulated_displacement))
        print('threshold: ', self.inlet_normal_distance)
        if (np.linalg.norm(self.accumulated_displacement) >= self.inlet_normal_distance):
            print("debug: is_inlet_step? true")
            print("Generating material points at inlet")
            return True
        else:
            print("debug: is_inlet_step? false")
            return False


    def _IsVelocityNumberOrFunction(self):
        # if self.assigned_velocity.IsDouble():
        self.velocity_value_is_numeric = True

        # elif self.assigned_velocity.IsString():
        #     raise Exception('MPMInletProcess: space/time varying inlet is not implemented yet!')

        # else:
        #     raise Exception('MPMInletProcess: Please provide vector for the assigned velocity')



        #     self.value_is_function = True
        #     self.function_string = value.GetString()
        #     self.aux_function = KratosMultiphysics.GenericFunctionUtility(self.function_string, settings["local_axes"])

        #     if self.aux_function.DependsOnSpace():
        #         self.cpp_apply_function_utility = KratosMultiphysics.ApplyFunctionToNodesUtility(self.model_part.Nodes, self.aux_function )