import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

def Factory(settings, Model):
    if( not isinstance(settings, KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyTwoFluidsInletProcess(Model, settings["Parameters"])


class ApplyTwoFluidsInletProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings):
        KratosMultiphysics.Process.__init__(self)

        # settings for inlet with interface between fluids and separate velocities
        default_settings = KratosMultiphysics.Parameters("""
        {
            "model_part_name" : "",
            "interface_settings": {
                "interface_normal": [0.0,1.0,0.0],
                "point_on_interface": [0.0,0.25,0.0],
                "inlet_transition_radius": 0.05
            },
            "fluid_1_settings" : {
                "model_part_name": "fluid_1_inlet",
                "variable_name": "VELOCITY",
                "modulus": 0.1,
                "constrained": true,
                "direction": [1.0,0.0,0.0],
                "interval": [0.0,"End"]
            },
            "fluid_2_settings" : {
                "model_part_name": "fluid_2_inlet",
                "variable_name": "VELOCITY",
                "modulus": 0.0,
                "constrained": true,
                "direction": [1.0,0.0,0.0],
                "interval": [0.0,"End"]
            }
        }
        """)

        # Trick: allow "modulus" and "direction" to be a double or a string value (otherwise the ValidateAndAssignDefaults might fail).
        # checking (and adapting) type of the modulus for both fluids
        if (settings["fluid_1_settings"].Has("modulus")):
            if (settings["fluid_1_settings"]["modulus"].IsString()):
                default_settings["fluid_1_settings"]["modulus"].SetString("0.0")

        if (settings["fluid_2_settings"].Has("modulus")):
            if (settings["fluid_2_settings"]["modulus"].IsString()):
                default_settings["fluid_2_settings"]["modulus"].SetString("0.0")

        # checking (and adapting) type of direction for both fluids
        if (settings["fluid_1_settings"].Has("direction")):
            if (settings["fluid_1_settings"]["direction"].IsString()):
                default_settings["fluid_1_settings"]["direction"].SetString("automatic_inwards_normal")
        if (settings["fluid_2_settings"].Has("direction")):
            if (settings["fluid_2_settings"]["direction"].IsString()):
                default_settings["fluid_2_settings"]["direction"].SetString("automatic_inwards_normal")

        # compare against the appropriate default settings
        settings.ValidateAndAssignDefaults(default_settings)

        # checking for empty model part name
        if (settings["model_part_name"].GetString() == ""):
            raise Exception("ApplyTwoFluidsInletProcess: A value (string) for the entry 'model_part_name' must be given in the parameters of the process.")

        # Check input data for only allowed value
        if (settings["fluid_1_settings"]["model_part_name"].GetString() != "fluid_1_inlet"):
            raise Exception("Entry 'model_part_name' in 'fluid_1_settings' must be 'fluid_1_inlet'. Other settings are not admissible.")
        elif (settings["fluid_2_settings"]["model_part_name"].GetString() != "fluid_2_inlet"):
            raise Exception("Entry 'model_part_name' in 'fluid_2_settings' must be 'fluid_2_inlet'. Other settings are not admissible.")
        elif (settings["fluid_1_settings"]["variable_name"].GetString() != "VELOCITY"):
            raise Exception("Inlet 'variable_name' in 'fluid_1_settings' is not VELOCITY. This is not admissible.")
        elif (settings["fluid_2_settings"]["variable_name"].GetString() != "VELOCITY"):
            raise Exception("Inlet 'variable_name' in 'fluid_2_settings' is not VELOCITY. This is not admissible.")
        else:

            # checking for empty strings
            if (settings["fluid_1_settings"]["modulus"].IsString()):
                if (settings["fluid_1_settings"]["modulus"].GetString() == ""):
                    raise Exception("Inlet function string for 'modulus' in 'fluid_1_settings' at the inlet is empty.")
            if (settings["fluid_2_settings"]["modulus"].IsString()):
                if (settings["fluid_2_settings"]["modulus"].GetString() == ""):
                    raise Exception("Inlet function string for 'modulus' in 'fluid_2_settings' at the inlet is empty.")

            # checking for different directions
            if ( (settings["fluid_1_settings"]["direction"].IsString() and settings["fluid_2_settings"]["direction"].IsString() )
                or (settings["fluid_1_settings"]["direction"].IsVector() and settings["fluid_2_settings"]["direction"].IsVector() ) ):

                if ( settings["fluid_1_settings"]["direction"].IsString() and settings["fluid_2_settings"]["direction"].IsString() ):
                    if ( settings["fluid_1_settings"]["direction"].GetString() != settings["fluid_2_settings"]["direction"].GetString() ):
                        raise Exception("The values (string) of the entry 'direction' differ between the two fluids. This is not admissible.")
                if ( settings["fluid_1_settings"]["direction"].IsVector() and settings["fluid_2_settings"]["direction"].IsVector() ):
                    v1 = settings["fluid_1_settings"]["direction"].GetVector()
                    v2 = settings["fluid_2_settings"]["direction"].GetVector()
                    if ( abs(v1[0] - v2[0]) > 1.0e-7 or abs(v1[1] - v2[1]) > 1.0e-7 or abs(v1[2] - v2[2]) > 1.0e-7 ):
                        raise Exception("The values (vector) of the entry 'direction' differ between the two fluids. This is not admissible.")
            else:
                raise Exception("The values (different type) of the entry 'direction' differ between the two fluids. This is not admissible.")

        self.inlet_model_part = Model[settings["model_part_name"].GetString()]
        self.complete_model = self.inlet_model_part.GetRootModelPart()

        self.neighbour_search = KratosMultiphysics.FindNodalNeighboursProcess(self.complete_model)
        self.neighbour_search.Execute()

        self.variational_distance_process = self.set_variational_distance_process()
        self.two_fluid_inlet_process = KratosFluid.TwoFluidsInletProcess( self.inlet_model_part, settings["interface_settings"], self.variational_distance_process )

        import KratosMultiphysics.assign_vector_by_direction_process as assign_vector_by_direction_process

        if ( self.inlet_model_part.GetSubModelPart("fluid_1_inlet").NumberOfNodes() > 0):
            self.aux_process_fluid_1 = assign_vector_by_direction_process.AssignVectorByDirectionProcess(Model, settings["fluid_1_settings"])

        if ( self.inlet_model_part.GetSubModelPart("fluid_2_inlet").NumberOfNodes() > 0):
            self.aux_process_fluid_2 = assign_vector_by_direction_process.AssignVectorByDirectionProcess(Model, settings["fluid_2_settings"])


    def ExecuteInitializeSolutionStep(self):
        # Call the base process ExecuteInitializeSolutionStep()
        if ( self.inlet_model_part.GetSubModelPart("fluid_1_inlet").NumberOfNodes() > 0):
            self.aux_process_fluid_1.ExecuteInitializeSolutionStep()
        if ( self.inlet_model_part.GetSubModelPart("fluid_2_inlet").NumberOfNodes() > 0):
            self.aux_process_fluid_2.ExecuteInitializeSolutionStep()

        self.two_fluid_inlet_process.SmoothDistanceField()


    def set_variational_distance_process(self):
        # Construct the variational distance calculation process
        maximum_iterations = 5

        ### for serial and OpenMP
        serial_settings = KratosMultiphysics.Parameters("""
        {
            "linear_solver_settings"   : {
                "solver_type" : "amgcl"
            }
        }
        """)
        import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory
        self.linear_solver = linear_solver_factory.ConstructSolver(serial_settings["linear_solver_settings"])

        if self.complete_model.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
            variational_distance_process = KratosMultiphysics.VariationalDistanceCalculationProcess2D(
                self.complete_model,
                self.linear_solver,
                maximum_iterations,
                KratosMultiphysics.VariationalDistanceCalculationProcess2D.CALCULATE_EXACT_DISTANCES_TO_PLANE,
                "distance_from_inlet_2D")
        else:
            variational_distance_process = KratosMultiphysics.VariationalDistanceCalculationProcess3D(
                self.complete_model,
                self.linear_solver,
                maximum_iterations,
                KratosMultiphysics.VariationalDistanceCalculationProcess3D.CALCULATE_EXACT_DISTANCES_TO_PLANE,
                "distance_from_inlet_3D")

        return variational_distance_process