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


        if Model.HasModelPart( settings["model_part_name"].GetString() ):
            self.inlet_model_part = Model[settings["model_part_name"].GetString()]
        else:
            raise Exception("The model does not contain a part with 'model_part_name' = " + settings["model_part_name"].GetString())


        self.neighbour_search = KratosMultiphysics.FindNodalNeighboursProcess(self.complete_model, 10, 10)
        self.neighbour_search.Execute()

        # TODO: Insert the constrcutor of the C++ process
        self.two_fluid_inlet_process = KratosFluid.TwoFluidsInletProcess(self.inlet_model_part, settings)

        # removing the additional entries to restore the format
        settings.RemoveValue("two_fluid_settings")
        import assign_vector_by_direction_process

        # adapt the (base) settings for the water values
        if ( self.modulus_water.IsString() ):
            settings["modulus"].SetString(self.modulus_water)
        elif ( self.modulus_water.IsDouble() ):
            settings["modulus"].SetDouble(self.modulus_water)
        settings["model_part_name"].SetString("water_inlet")
        self.aux_process_water = assign_vector_by_direction_process.AssignVectorByDirectionProcess(Model, settings)

        # adapt the (base) settings for the air values
        if ( self.modulus_air.IsString() ):
            settings["modulus"].SetString(self.modulus_air)
        elif ( self.modulus_air.IsDouble() ):
            settings["modulus"].SetDouble(self.modulus_air)
        settings["model_part_name"].SetString("air_inlet")
        self.aux_process_air = assign_vector_by_direction_process.AssignVectorByDirectionProcess(Model, settings)


    def ExecuteInitializeSolutionStep(self):
        # Call the base process ExecuteInitializeSolutionStep()
        self.aux_process_water.ExecuteInitializeSolutionStep()
        self.aux_process_air.ExecuteInitializeSolutionStep()


    def ExecuteFinalizeSolutionStep(self):
        # Call the base process ExecuteFinalizeSolutionStep()
        self.aux_process_water.ExecuteInitializeSolutionStep()
        self.aux_process_air.ExecuteInitializeSolutionStep()

        # re-initial the distance field and blending the inlet field and the field inside the domain

        for node in self.main_model_part.Nodes:
            if ( node.GetValue(KratosFluid.AUX_DISTANCE) > 0.0 ):
                weighting_factor_inlet_field = node.GetValue(KratosFluid.AUX_DISTANCE)
                weighting_factor_domain_field = 1.0 - weighting_factor_inlet_field

                distance_vector = [ node.X - self.point_on_interface[0], node.Y - self.point_on_interface[1], node.Z - self.point_on_interface[2] ]
                inlet_distance_field =  self.interface_normal[0]*distance_vector[0] + self.interface_normal[1]*distance_vector[1] + self.interface_normal[2]*distance_vector[2]

                domain_distance_field = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
                value = weighting_factor_inlet_field * inlet_distance_field + weighting_factor_domain_field * domain_distance_field
                node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, value)