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

        # introduction of flag TWO_FLUIDS should be discussed first
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
            # settings for inlet with interface between fluids and separate velocities
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
                    "modulus_air"               : 0.1,
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

        # compare against the appropriate default settings (one or two fluids)
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
        for node in self.inlet_model_part.Nodes:
            node.Set(KratosMultiphysics.INLET, True)
        for condition in self.inlet_model_part.Conditions:
            condition.Set(KratosMultiphysics.INLET, True)


        if ( self.preliminray_two_fluid_flag == True ):
            self.inlet_radius = settings["two_fluid_settings"]["inlet_radius"].GetDouble()
            self.interface_normal = settings["two_fluid_settings"]["interface_normal"].GetVector()
            self.point_on_interface = settings["two_fluid_settings"]["point_on_interface"].GetVector()
            self.modulus_water = settings["two_fluid_settings"]["modulus_water"].GetDouble()
            self.modulus_air = settings["two_fluid_settings"]["modulus_air"].GetDouble()

            # needed to calculate the distance field from the inlet that is used to blend the distance fields
            self.main_model_part = self.inlet_model_part.GetRootModelPart()
            self.neighbour_search = KratosMultiphysics.FindNodalNeighboursProcess(self.main_model_part, 10, 10)
            self.dimension = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
            self.neighbour_search.Execute()

            self.distance_calculator = KratosMultiphysics.BodyDistanceCalculationUtils()

            # set the inlet nodes as the origin for the distance-from-inlet calculation
            KratosMultiphysics.VariableUtils().SetNonHistoricalVariable( KratosMultiphysics.IS_VISITED, 0.0, self.main_model_part.Nodes)
            KratosMultiphysics.VariableUtils().SetNonHistoricalVariable( KratosMultiphysics.IS_VISITED, 1.0, self.inlet_model_part.Nodes)

            # temporally storing the distance field as an older version of itself (it can be assured that nothing is over-written at the start)
            for node in self.main_model_part.Nodes:
                node.SetSolutionStepValue( KratosMultiphysics.DISTANCE, 2, node.GetSolutionStepValue(KratosMultiphysics.DISTANCE, 0) )

            # working on the distance-from-inlet field that will be finally stored as a non-historical variable
            KratosMultiphysics.VariableUtils().SetScalarVar( KratosMultiphysics.DISTANCE, -self.inlet_radius, self.inlet_model_part.Nodes)
            if ( self.dimension == 2 ):
                self.distance_calculator.CalculateDistances2D(self.main_model_part.Elements, KratosMultiphysics.DISTANCE, 0.0, True)
            elif ( self.dimension == 3 ):
                self.distance_calculator.CalculateDistances3D(self.main_model_part.Elements, KratosMultiphysics.DISTANCE, 0.0, True)
            # scaling to achieve a value of exactly one at the inlet
            scaling_factor = - 1.0 / self.inlet_radius
            for node in self.main_model_part.Nodes:
                node.SetSolutionStepValue( KratosMultiphysics.DISTANCE, 0, scaling_factor * node.GetSolutionStepValue(KratosMultiphysics.DISTANCE, 0) )

            # copying the computed distance-from-inlet field to a non-historical variable
            KratosMultiphysics.VariableUtils().SaveScalarVar( KratosMultiphysics.DISTANCE, KratosFluid.DISTANCE_FROM_INLET, self.main_model_part.Nodes )

            # restoring the original distance field from its stored version
            for node in self.main_model_part.Nodes:
                node.SetSolutionStepValue( KratosMultiphysics.DISTANCE, 0, node.GetSolutionStepValue(KratosMultiphysics.DISTANCE, 2) )

            # splitting the inlet model part to assign different velocities for water and air
            self.water_inlet_name = "water_inlet"
            self.air_inlet_name = "air_inlet"
            self.inlet_model_part.CreateSubModelPart(self.water_inlet_name)
            self.inlet_model_part.CreateSubModelPart(self.air_inlet_name)

            self.water_inlet_sub_model_part = self.inlet_model_part.GetSubModelPart(self.water_inlet_name)
            self.air_inlet_sub_model_part = self.inlet_model_part.GetSubModelPart(self.air_inlet_name)

            normal_length = math.sqrt( self.interface_normal[0]*self.interface_normal[0] +
                                        self.interface_normal[1]*self.interface_normal[1] +
                                        self.interface_normal[2]*self.interface_normal[2] )
            self.interface_normal = [ self.interface_normal[0] / normal_length, self.interface_normal[1] / normal_length, self.interface_normal[2] / normal_length ]

            # separating nodes by means of their distance value at the inlet
            list_of_water_nodes = set()
            list_of_air_nodes = set()
            for node in self.inlet_model_part.Nodes:
                distance_vector = [ node.X - self.point_on_interface[0], node.Y - self.point_on_interface[1], node.Z - self.point_on_interface[2] ]
                inlet_distance_field =  self.interface_normal[0]*distance_vector[0] + self.interface_normal[1]*distance_vector[1] + self.interface_normal[2]*distance_vector[2]
                if ( inlet_distance_field < 0.0 ):
                    list_of_water_nodes.add(node.Id)
                else:
                    list_of_air_nodes.add(node.Id)
            # assigning nodes to the respective sub model parts
            self.water_inlet_sub_model_part.AddNodes( list(list_of_water_nodes) )
            self.air_inlet_sub_model_part.AddNodes( list(list_of_air_nodes) )

            # removing the additional entries to restore the format
            settings.RemoveValue("two_fluid_settings")
            import assign_vector_by_direction_process

            # adapt the (base) settings for the water values
            settings["modulus"].SetDouble(self.modulus_water)
            settings["model_part_name"].SetString(self.water_inlet_name)
            self.aux_process_water = assign_vector_by_direction_process.AssignVectorByDirectionProcess(Model, settings)
            # adapt the (base) settings for the air values
            settings["modulus"].SetDouble(self.modulus_air)
            settings["model_part_name"].SetString(self.air_inlet_name)
            self.aux_process_air = assign_vector_by_direction_process.AssignVectorByDirectionProcess(Model, settings)

        else:

            import assign_vector_by_direction_process
            self.aux_process = assign_vector_by_direction_process.AssignVectorByDirectionProcess(Model, settings)


    def ExecuteInitializeSolutionStep(self):
        # Call the base process ExecuteInitializeSolutionStep()
        if ( self.preliminray_two_fluid_flag == True ):
            self.aux_process_water.ExecuteInitializeSolutionStep()
            self.aux_process_air.ExecuteInitializeSolutionStep()
        else:
            self.aux_process.ExecuteInitializeSolutionStep()


    def ExecuteFinalizeSolutionStep(self):
        # Call the base process ExecuteFinalizeSolutionStep()
        if ( self.preliminray_two_fluid_flag == True ):

            # executing both processes to set the respective velocity fields
            self.aux_process_water.ExecuteInitializeSolutionStep()
            self.aux_process_air.ExecuteInitializeSolutionStep()
            # re-initial the distance field and blending the inlet field and the field inside the domain

            for node in self.main_model_part.Nodes:
                if ( node.GetValue(KratosFluid.DISTANCE_FROM_INLET) > 0.0 ):
                    weighting_factor_inlet_field = node.GetValue(KratosFluid.DISTANCE_FROM_INLET)
                    weighting_factor_domain_field = 1.0 - weighting_factor_inlet_field

                    distance_vector = [ node.X - self.point_on_interface[0], node.Y - self.point_on_interface[1], node.Z - self.point_on_interface[2] ]
                    inlet_distance_field =  self.interface_normal[0]*distance_vector[0] + self.interface_normal[1]*distance_vector[1] + self.interface_normal[2]*distance_vector[2]

                    domain_distance_field = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
                    value = weighting_factor_inlet_field * inlet_distance_field + weighting_factor_domain_field * domain_distance_field
                    node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, value)

        else:
            self.aux_process.ExecuteInitializeSolutionStep()