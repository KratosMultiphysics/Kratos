# Importing the Kratos Library
import KratosMultiphysics as KM

import KratosMultiphysics.ContactStructuralMechanicsApplication as CSMA

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return MeshTyingProcess(Model, settings["Parameters"])

import KratosMultiphysics.ContactStructuralMechanicsApplication.search_base_process as search_base_process

class MeshTyingProcess(search_base_process.SearchBaseProcess):
    """This class is used in order to compute the a mortar mesh tying formulation

    This class constructs the model parts containing the mesh tying conditions and
    initializes parameters and variables related with the mesh tying. The class creates
    search utilities to be used to create the tying pairs

    Only the member variables listed below should be accessed directly.

    Public member variables:
    Model -- the model used to construct the process.
    settings -- Kratos parameters containing solver settings.
    """

    def __init__(self, Model, settings):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        Model -- the model part used to construct the process.
        settings -- Kratos parameters containing solver settings.
        """

        # NOTE: Due to recursive check "search_model_part" and "assume_master_slave" requires to pre-define configurations, if more that 10 pairs of contact are required, just add. I assume nobody needs that much
        # Settings string in json format
        default_parameters = KM.Parameters("""
        {
            "help"                         : "This class is used in order to compute the a mortar mesh tying formulation. This class constructs the model parts containing the mesh tying conditions and initializes parameters and variables related with the mesh tying. The class creates search utilities to be used to create the tying pairs",
            "model_part_name"              : "Structure",
            "mesh_tying_model_part"        : {"0":[],"1":[],"2":[],"3":[],"4":[],"5":[],"6":[],"7":[],"8":[],"9":[]},
            "assume_master_slave"          : {"0":[],"1":[],"2":[],"3":[],"4":[],"5":[],"6":[],"7":[],"8":[],"9":[]},
            "mesh_tying_property_ids"      : {"0": 0,"1": 0,"2": 0,"3": 0,"4": 0,"5": 0,"6": 0,"7": 0,"8": 0,"9": 0},
            "interval"                     : [0.0,"End"],
            "variable_name"                : "DISPLACEMENT",
            "consider_static_condensation" : false,
            "zero_tolerance_factor"        : 1.0,
            "integration_order"            : 2,
            "consider_tessellation"        : true,
            "normal_check_proportion"      : 0.1,
            "search_parameters"            : {
                "type_search"                 : "in_radius_with_obb",
                "search_factor"               : 3.5,
                "active_check_factor"         : 0.01,
                "max_number_results"          : 1000,
                "bucket_size"                 : 4,
                "dynamic_search"              : false,
                "database_step_update"        : 999999999,
                "debug_mode"                  : false,
                "check_gap"                   : "check_mapping",
                "octree_search_parameters"    : {
                    "bounding_box_factor"             : 0.1,
                    "debug_obb"                       : false,
                    "OBB_intersection_type"           : "SeparatingAxisTheorem",
                    "lower_bounding_box_coefficient"  : 0.0,
                    "higher_bounding_box_coefficient" : 1.0
                }
            },
            "scale_factor_parameters"      : {
                "manual_scale_factor"         : false,
                "stiffness_factor"            : 1.0,
                "scale_factor"                : 1.0e0
            }
        }
        """)

        # Overwrite the default settings with user-provided parameters
        self.mesh_tying_settings = settings
        self.mesh_tying_settings.ValidateAndAssignDefaults(default_parameters)

        # We transfer the parameters to the base class
        base_process_settings = KM.Parameters("""{}""")
        base_process_settings.AddValue("search_model_part", self.mesh_tying_settings["mesh_tying_model_part"])
        base_process_settings.AddValue("search_property_ids", self.mesh_tying_settings["mesh_tying_property_ids"])
        parameter_list = ["model_part_name", "assume_master_slave", "interval", "zero_tolerance_factor", "integration_order", "consider_tessellation", "normal_check_proportion", "search_parameters"]
        base_process_settings.CopyValuesFromExistingParameters(self.mesh_tying_settings, parameter_list)

        # Construct the base process.
        super().__init__(Model, base_process_settings)

        # Mesh tying configurations
        # Determine if the variable is components or scalar
        self.variable_name = self.mesh_tying_settings["variable_name"].GetString()
        if KM.KratosGlobals.HasVariable(self.variable_name):
            var_type = KM.KratosGlobals.GetVariableType(self.variable_name)
            if var_type == "Array":
                self.type_variable = "Components"
            elif var_type == "Double":
                self.type_variable = "Scalar"
            else:
                raise Exception("Variable " + self.variable_name + " not compatible")
        else:
            raise Exception("Variable " + self.variable_name + " not registered")
        
        # If we consider static condensation
        self.consider_static_condensation = self.mesh_tying_settings["consider_static_condensation"].GetBool()

    def ExecuteInitialize(self):
        """ This method is executed at the begining to initialize the process

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # We call to the base process
        super().ExecuteInitialize()

        # Scale factor settings
        if self.mesh_tying_settings["scale_factor_parameters"]["manual_scale_factor"].GetBool():
            self.scale_factor = self.mesh_tying_settings["scale_factor_parameters"]["scale_factor"].GetDouble()
        else:
            scale_factor_var_parameters = KM.Parameters("""{
                "compute_penalty"      : false
            }""")
            scale_factor_var_parameters.AddValue("stiffness_factor", self.mesh_tying_settings["scale_factor_parameters"]["stiffness_factor"])
            self.scale_factor_process = CSMA.ALMVariablesCalculationProcess(self._get_process_model_part(), KM.NODAL_H, scale_factor_var_parameters)
            self.scale_factor_process.Execute()

        # If we consider static condensation
        if self.consider_static_condensation:
            computing_contact_model_part = self.main_model_part.GetSubModelPart("ComputingContact")
            self.assign_elements_conditions_process = CSMA.AssignParentElementConditionsProcess(computing_contact_model_part, self.main_model_part)
            self.assign_elements_conditions_process.ExecuteInitialize()

    def ExecuteBeforeSolutionLoop(self):
        """ This method is executed before starting the time loop

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        # We call to the base process
        super().ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):
        """ This method is executed in order to initialize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        # We call to the base process
        super().ExecuteInitializeSolutionStep()

        # If we consider static condensation
        if self.consider_static_condensation:
            self.assign_elements_conditions_process.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        """ This method is executed in order to finalize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        # We call to the base process
        super().ExecuteFinalizeSolutionStep()

    def ExecuteBeforeOutputStep(self):
        """ This method is executed right before the ouput process computation

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        # We call to the base process
        super().ExecuteBeforeOutputStep()

    def ExecuteAfterOutputStep(self):
        """ This method is executed right after the ouput process computation

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        # We call to the base process
        super().ExecuteAfterOutputStep()

    def ExecuteFinalize(self):
        """ This method is executed in order to finalize the current computation

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        # We call to the base process
        super().ExecuteFinalize()

    def _get_condition_name(self):
        """ This method returns the condition name

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # We define the condition name to be used
        return "MeshTyingMortar"

    def _initialize_search_conditions(self):
        """ This method initializes some conditions values

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # We call to the base process
        super()._initialize_search_conditions()

        # Setting tying variable
        for prop in self._get_process_model_part().GetProperties():
            prop[CSMA.TYING_VARIABLE] = self.variable_name

        # Initializing some values
        zero_vector = KM.Vector(3)
        zero_vector[0] = 0.0
        zero_vector[1] = 0.0
        zero_vector[2] = 0.0

        # Initilialize weighted variables and LM
        if self.type_variable == "Scalar":
            KM.VariableUtils().SetVariable(CSMA.WEIGHTED_SCALAR_RESIDUAL, 0.0, self._get_process_model_part().Nodes)
        else:
            KM.VariableUtils().SetVariable(CSMA.WEIGHTED_VECTOR_RESIDUAL, zero_vector, self._get_process_model_part().Nodes)

        # Setting the conditions
        KM.VariableUtils().SetNonHistoricalVariable(KM.NORMAL, zero_vector, self._get_process_model_part().Conditions)
