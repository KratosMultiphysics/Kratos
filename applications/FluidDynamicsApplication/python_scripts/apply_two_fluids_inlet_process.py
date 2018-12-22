import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
import math

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyTwoFluidsInletProcess(Model, settings["Parameters"])


class ApplyTwoFluidsInletProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings):
        KratosMultiphysics.Process.__init__(self)

        # settings for inlet with interface between fluids and separate velocities
        default_settings = KratosMultiphysics.Parameters("""
        {
            "mesh_id"                   : 0,
            "model_part_name"           : "",
            "variable_name"             : "VELOCITY",
            "constrained"               : true,
            "direction"                 : [1.0,0.0,0.0],
            "interval"                  : [0.0,"End"],
            "two_fluid_settings" : {
                "modulus_air"               : 0.0,
                "modulus_water"             : 0.0,
                "interface_normal"          : [0.0,1.0,0.0],
                "point_on_interface"        : [0.0,0.25,0.0],
                "inlet_transition_radius"   : 0.05
            }
        }
        """)

        # Trick: allow "modulus" and "direction" to be a double or a string value (otherwise the ValidateAndAssignDefaults might fail).
        if (settings.Has("modulus_air")):
            if (settings["two_fluid_settings"]["modulus_air"].IsString()):
                default_settings["two_fluid_settings"]["modulus_air"].SetString("0.0")

        if (settings.Has("modulus_water")):
            if (settings["two_fluid_settings"]["modulus_water"].IsString()):
                default_settings["two_fluid_settings"]["modulus_water"].SetString("0.0")

        if (settings.Has("direction")):
            if (settings["direction"].IsString()):
                default_settings["direction"].SetString("automatic_inwards_normal")

        # compare against the appropriate default settings (one or two fluids)
        settings.ValidateAndAssignDefaults(default_settings)

        # Check input data
        if (settings["model_part_name"].GetString() == ""):
            raise Exception("Empty inlet model part name string. Set a valid model part name.")
        elif (settings["variable_name"].GetString() != "VELOCITY"):
            raise Exception("Inlet variable_name is not VELOCITY.")
        else:
            # checking for the air
            if (settings["two_fluid_settings"]["modulus_air"].IsDouble()):
                if (settings["two_fluid_settings"]["modulus_air"].GetDouble == 0.0):
                    KratosMultiphysics.Logger.PrintInfo("ApplyTwoFluidsInletProcess", "Inlet velocity for the air was set to 0.0" )
                self.modulus_air = settings["two_fluid_settings"]["modulus_air"].GetDouble()
            elif (settings["two_fluid_settings"]["modulus_air"].IsString()):
                if (settings["two_fluid_settings"]["modulus_air"].GetString == ""):
                    raise Exception("Inlet function string for the air velocity at the inlet is empty.")
                self.modulus_air = settings["two_fluid_settings"]["modulus_air"].GetString()
            # checking for the water
            if (settings["two_fluid_settings"]["modulus_water"].IsDouble()):
                if (settings["two_fluid_settings"]["modulus_water"].GetDouble == 0.0):
                    KratosMultiphysics.Logger.PrintInfo("ApplyTwoFluidsInletProcess", "Inlet velocity for the water was set to 0.0 - this may not be desired" )
                self.modulus_water = settings["two_fluid_settings"]["modulus_water"].GetDouble()
            elif (settings["two_fluid_settings"]["modulus_water"].IsString()):
                if (settings["two_fluid_settings"]["modulus_water"].GetString == ""):
                    raise Exception("Inlet function string for the water velocity at the inlet is empty.")
                self.modulus_water = settings["two_fluid_settings"]["modulus_water"].GetString()

        self.inlet_model_part = Model[settings["model_part_name"].GetString()]
        self.complete_model = self.inlet_model_part.GetRootModelPart()

        self.neighbour_search = KratosMultiphysics.FindNodalNeighboursProcess(self.complete_model, 10, 10)
        self.neighbour_search.Execute()

        # Adding the C++ constructor
        self.two_fluid_inlet_process = KratosFluid.TwoFluidsInletProcess(self.inlet_model_part, settings)

        self.point_on_interface = settings["two_fluid_settings"]["point_on_interface"].GetVector()
        self.interface_normal = settings["two_fluid_settings"]["interface_normal"].GetVector()

        # removing the additional entries to restore the format
        settings.RemoveValue("two_fluid_settings")
        settings.AddEmptyValue("modulus")
        import assign_vector_by_direction_process

        # adapt the (base) settings for the water values
        if ( isinstance(self.modulus_water, str) ):
            settings["modulus"].SetString(self.modulus_water)
        elif ( isinstance(self.modulus_water, float) ):
            settings["modulus"].SetDouble(self.modulus_water)
        # Sub model part "water_inlet" is defined and filled in KratosFluid.TwoFluidsInletProcess
        settings["model_part_name"].SetString("water_inlet")

        if ( self.inlet_model_part.GetSubModelPart("water_inlet").NumberOfNodes() > 0):
            self.aux_process_water = assign_vector_by_direction_process.AssignVectorByDirectionProcess(Model, settings)

        # adapt the (base) settings for the air values
        if ( isinstance(self.modulus_air, str) ):
            settings["modulus"].SetString(self.modulus_air)
        elif ( isinstance(self.modulus_air, float) ):
            settings["modulus"].SetDouble(self.modulus_air)
        # Sub model part "air_inlet" is defined and filled in KratosFluid.TwoFluidsInletProcess
        settings["model_part_name"].SetString("air_inlet")

        if ( self.inlet_model_part.GetSubModelPart("air_inlet").NumberOfNodes() > 0):
            self.aux_process_air = assign_vector_by_direction_process.AssignVectorByDirectionProcess(Model, settings)


    def ExecuteInitializeSolutionStep(self):
        # Call the base process ExecuteInitializeSolutionStep()
        if ( self.inlet_model_part.GetSubModelPart("water_inlet").NumberOfNodes() > 0):
            self.aux_process_water.ExecuteInitializeSolutionStep()
        if ( self.inlet_model_part.GetSubModelPart("air_inlet").NumberOfNodes() > 0):
            self.aux_process_air.ExecuteInitializeSolutionStep()

        # Not sure if appropriate here...
        # PROBLEM: Could distort the physical properties
        # BUT: It could make sense to stabilize cases with CFL >> 1
        self.two_fluid_inlet_process.SmoothDistanceField()


    def ExecuteFinalizeSolutionStep(self):
        # Call the base process ExecuteFinalizeSolutionStep()
        if ( self.inlet_model_part.GetSubModelPart("water_inlet").NumberOfNodes() > 0):
            self.aux_process_water.ExecuteInitializeSolutionStep()
        if ( self.inlet_model_part.GetSubModelPart("air_inlet").NumberOfNodes() > 0):
            self.aux_process_air.ExecuteInitializeSolutionStep()

        self.two_fluid_inlet_process.SmoothDistanceField()