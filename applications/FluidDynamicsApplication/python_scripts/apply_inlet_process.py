import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
import numpy as np

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
                "interface_normal"          : [0.0,1.0,0.0],
                "point_on_interface"        : [0.0,0.25,0.0],
                "inlet_radius"              : 0.05
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
        self.inlet_radius = settings["inlet_radius"].GetDouble()
        self.interface_normal = settings["interface_normal"].GetVector()
        self.point_on_interface = settings["point_on_interface"].GetVector()

        if ( self.preliminray_two_fluid_flag == True ):
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

            for node in self.main_model_part.Nodes:
                node.SetValue(KratosMultiphysics.IS_VISITED, 0.0)

            for node in self.inlet_model_part.Nodes:
                node.SetValue(KratosMultiphysics.IS_VISITED, 1.0)
                node.SetSolutionStepValue(KratosFluid.DISTANCE_FROM_INLET, 0, -self.inlet_radius)

            if ( self.dimension == 2 ):
                self.distance_calculator.CalculateDistances2D(self.main_model_part.Elements, KratosFluid.DISTANCE_FROM_INLET, 0.0)
            elif ( self.dimension == 3 ):
                self.distance_calculator.CalculateDistances3D(self.main_model_part.Elements, KratosFluid.DISTANCE_FROM_INLET, 0.0)

            for node in self.main_model_part.Nodes:
                scaling_factor = - 1.0 / self.inlet_radius
                node.SetSolutionStepValue( KratosFluid.DISTANCE_FROM_INLET, 0, scaling_factor * node.GetSolutionStepValue(KratosFluid.DISTANCE_FROM_INLET) )

        import assign_vector_by_direction_process
        # Construct the base process AssignVectorByDirectionProcess
        if ( self.preliminray_two_fluid_flag == True) :
            base_settings = settings.Clone()
            base_settings.RemoveValue("interface_normal")
            base_settings.RemoveValue("point_on_interface")
            base_settings.RemoveValue("inlet_radius")
            self.aux_process = assign_vector_by_direction_process.AssignVectorByDirectionProcess(Model, base_settings)
        else:
            self.aux_process = assign_vector_by_direction_process.AssignVectorByDirectionProcess(Model, settings)


    def ExecuteInitializeSolutionStep(self):
        # Call the base process ExecuteInitializeSolutionStep()
        self.aux_process.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        # Call the base process ExecuteFinalizeSolutionStep()
        self.aux_process.ExecuteFinalizeSolutionStep()

        print( self.inlet_radius )

        for node in self.main_model_part.Nodes:

            weighting_factor_inlet_field = node.GetSolutionStepValue(KratosFluid.DISTANCE_FROM_INLET)
            weighting_factor_domain_field = 1.0 - weighting_factor_inlet_field

            start_point = np.array( [ self.point_on_interface[0], self.point_on_interface[1], self.point_on_interface[2] ] )
            normal = np.array( [ self.interface_normal[0], self.interface_normal[1], self.interface_normal[2] ] )
            normal = normal / np.linalg.norm(normal)
            node_position = np.array( [node.X, node.Y, node.Z] )

            inlet_distance_field = np.dot( node_position - start_point, normal )
            domain_distance_field = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)

            value = weighting_factor_inlet_field * inlet_distance_field + weighting_factor_domain_field * domain_distance_field

            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, value)
