import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
import math

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyInletProcess(Model, settings["Parameters"])


class ApplyInletProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings):
        KratosMultiphysics.Process.__init__(self)

        self.preliminray_two_fluid_flag = True

        if (self.preliminray_two_fluid_flag == False):
            # settings for regular inlet
            default_settings = KratosMultiphysics.Parameters("""
            {
                "mesh_id"         : 0,
                "model_part_name" : "",
                "variable_name"   : "VELOCITY",
                "modulus"         : 0.0,
                "constrained"     : true,
                "direction"       : [1.0,0.0,0.0],
                "interval"        : [0.0,"End"]
            }
            """)

        elif (self.preliminray_two_fluid_flag == True):
            # settings for inlet with interface between fluids
            default_settings = KratosMultiphysics.Parameters("""
            {
                "mesh_id"                   : 0,
                "model_part_name"           : "",
                "variable_name"             : "VELOCITY",
                "modulus"                   : 0.0,
                "constrained"               : true,
                "direction"                 : [1.0,0.0,0.0],
                "interval"                  : [0.0,"End"],
                "two_fluid_settings" : {
                    "modulus_air"               : 0.2,
                    "modulus_water"             : 1.0,
                    "interface_normal"          : [0.0,1.0,0.0],
                    "point_on_interface"        : [0.0,0.25,0.0],
                    "inlet_radius"              : 0.05
                }
            }
            """)

        # Trick: allow "modulus" and "direction" to be a double or a string value (otherwise the ValidateAndAssignDefaults might fail)
        if (settings.Has("modulus")):
            if (settings["modulus"].IsString()):
                default_settings["modulus"].SetString("0.0")

        if (settings.Has("direction")):
            if (settings["direction"].IsString()):
                default_settings["direction"].SetString("automatic_inwards_normal")

        settings.ValidateAndAssignDefaults(default_settings)

        # Check input data
        if (settings["model_part_name"].GetString() == ""):
            raise Exception("Empty inlet model part name string. Set a valid model part name.")
        elif (settings["variable_name"].GetString() != "VELOCITY"):
            raise Exception("Inlet variable_name is not VELOCITY.")
        else:
            if (settings["modulus"].IsDouble()):
                if (settings["modulus"].GetDouble == 0.0):
                    raise Exception("Inlet scalar value equal to 0.")
            elif (settings["modulus"].IsString()):
                if (settings["modulus"].GetString == ""):
                    raise Exception("Inlet function sting is empty.")

        # Set the INLET flag in the inlet model part nodes and conditions
        self.inlet_model_part = Model[settings["model_part_name"].GetString()]

        if ( self.preliminray_two_fluid_flag == True ):
            self.inlet_radius = settings["two_fluid_settings"]["inlet_radius"].GetDouble()
            self.interface_normal = settings["two_fluid_settings"]["interface_normal"].GetVector()
            self.point_on_interface = settings["two_fluid_settings"]["point_on_interface"].GetVector()

            # needed for the inlet distance calculation
            self.main_model_part = self.inlet_model_part.GetRootModelPart()
            print( self.main_model_part.Nodes )
            self.neighbour_search = KratosMultiphysics.FindNodalNeighboursProcess(self.main_model_part, 10, 10)
            self.dimension = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
            self.neighbour_search.Execute()

        for node in self.inlet_model_part.Nodes:
            node.Set(KratosMultiphysics.INLET, True)
        for condition in self.inlet_model_part.Conditions:
            condition.Set(KratosMultiphysics.INLET, True)

        # Compute distance-from-inlet function
        if( self.preliminray_two_fluid_flag == True):

            self.distance_calculator = KratosMultiphysics.BodyDistanceCalculationUtils()

            # mind that the variable is not endowed with a hsitory in the buffer - thus this function is needed
            KratosMultiphysics.VariableUtils().SetNonHistoricalVariable( KratosMultiphysics.IS_VISITED, 0.0, self.main_model_part.Nodes)
            KratosMultiphysics.VariableUtils().SetNonHistoricalVariable( KratosMultiphysics.IS_VISITED, 1.0, self.inlet_model_part.Nodes)
            KratosMultiphysics.VariableUtils().SetNonHistoricalVariable( KratosFluid.DISTANCE_FROM_INLET, -self.inlet_radius, self.inlet_model_part.Nodes)

            if ( self.dimension == 2 ):
                # new version of the function is used (!!!)
                self.distance_calculator.CalculateDistances2D(self.main_model_part.Elements, KratosFluid.DISTANCE_FROM_INLET, 0.0, False)
            elif ( self.dimension == 3 ):
                # new version of the function is used (!!!)
                self.distance_calculator.CalculateDistances3D(self.main_model_part.Elements, KratosFluid.DISTANCE_FROM_INLET, 0.0, False)

            for node in self.main_model_part.Nodes:
                # is a scaling operation possible with VariableUtils()
                scaling_factor = - 1.0 / self.inlet_radius
                node.SetValue( KratosFluid.DISTANCE_FROM_INLET, scaling_factor * node.GetValue(KratosFluid.DISTANCE_FROM_INLET) )

            for node in self.main_model_part.Nodes:
                # print( node.GetSolutionStepValue(KratosMultiphysics.IS_VISITED) )
                # print( node.GetSolutionStepValue(KratosFluid.DISTANCE_FROM_INLET) )
                print( node.GetValue(KratosFluid.DISTANCE_FROM_INLET) )

            # input( "Paused ..." )

        import assign_vector_by_direction_process
        # Construct the base process AssignVectorByDirectionProcess
        if ( self.preliminray_two_fluid_flag == True) :
            base_settings = settings.Clone()
            base_settings.RemoveValue("two_fluid_settings")
            self.aux_process = assign_vector_by_direction_process.AssignVectorByDirectionProcess(Model, base_settings)
        else:
            self.aux_process = assign_vector_by_direction_process.AssignVectorByDirectionProcess(Model, settings)


    def ExecuteInitializeSolutionStep(self):
        # Call the base process ExecuteInitializeSolutionStep()
        self.aux_process.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        # Call the base process ExecuteFinalizeSolutionStep()
        self.aux_process.ExecuteFinalizeSolutionStep()

        for node in self.main_model_part.Nodes:

            weighting_factor_inlet_field = node.GetValue(KratosFluid.DISTANCE_FROM_INLET)
            weighting_factor_domain_field = 1.0 - weighting_factor_inlet_field

            vector_length = math.sqrt(  self.interface_normal[0]*self.interface_normal[0] +
                                        self.interface_normal[1]*self.interface_normal[1] +
                                        self.interface_normal[2]*self.interface_normal[2] )

            self.interface_normal = [ self.interface_normal[0] / vector_length, self.interface_normal[1] / vector_length, self.interface_normal[2] / vector_length ]

            distance_vector = [ node.X - self.point_on_interface[0], node.Y - self.point_on_interface[1], node.Z - self.point_on_interface[2] ]

            inlet_distance_field =  self.interface_normal[0]*distance_vector[0] + self.interface_normal[1]*distance_vector[1] + self.interface_normal[2]*distance_vector[2]

            domain_distance_field = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
            value = weighting_factor_inlet_field * inlet_distance_field + weighting_factor_domain_field * domain_distance_field
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, value)
