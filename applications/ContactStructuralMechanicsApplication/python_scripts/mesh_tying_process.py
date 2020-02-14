from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
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
            "help"                        : "This class is used in order to compute the a mortar mesh tying formulation. This class constructs the model parts containing the mesh tying conditions and initializes parameters and variables related with the mesh tying. The class creates search utilities to be used to create the tying pairs",
            "mesh_id"                     : 0,
            "model_part_name"             : "Structure",
            "mesh_tying_model_part"       : {"0":[],"1":[],"2":[],"3":[],"4":[],"5":[],"6":[],"7":[],"8":[],"9":[]},
            "assume_master_slave"         : {"0":[],"1":[],"2":[],"3":[],"4":[],"5":[],"6":[],"7":[],"8":[],"9":[]},
            "mesh_tying_property_ids"     : {"0": 0,"1": 0,"2": 0,"3": 0,"4": 0,"5": 0,"6": 0,"7": 0,"8": 0,"9": 0},
            "interval"                    : [0.0,"End"],
            "variable_name"               : "DISPLACEMENT",
            "zero_tolerance_factor"       : 1.0,
            "search_parameters" : {
                "type_search"                 : "in_radius_with_obb",
                "search_factor"               : 3.5,
                "active_check_factor"         : 0.01,
                "max_number_results"          : 1000,
                "bucket_size"                 : 4,
                "dynamic_search"              : false,
                "database_step_update"        : 999999999,
                "debug_mode"                  : false,
                "check_gap"                   : "check_mapping",
                "octree_search_parameters" : {
                    "bounding_box_factor"             : 0.1,
                    "debug_obb"                       : false,
                    "OBB_intersection_type"           : "SeparatingAxisTheorem",
                    "lower_bounding_box_coefficient"  : 0.0,
                    "higher_bounding_box_coefficient" : 1.0
                }
            },
            "integration_order"           : 2
        }
        """)

        # Overwrite the default settings with user-provided parameters
        self.mesh_tying_settings = settings
        self.mesh_tying_settings.ValidateAndAssignDefaults(default_parameters)

        # We transfer the parameters to the base class
        base_process_settings = KM.Parameters("""{}""")
        base_process_settings.AddValue("mesh_id", self.mesh_tying_settings["mesh_id"])
        base_process_settings.AddValue("model_part_name", self.mesh_tying_settings["model_part_name"])
        base_process_settings.AddValue("search_model_part", self.mesh_tying_settings["mesh_tying_model_part"])
        base_process_settings.AddValue("assume_master_slave", self.mesh_tying_settings["assume_master_slave"])
        base_process_settings.AddValue("search_property_ids", self.mesh_tying_settings["mesh_tying_property_ids"])
        base_process_settings.AddValue("interval", self.mesh_tying_settings["interval"])
        base_process_settings.AddValue("zero_tolerance_factor", self.mesh_tying_settings["zero_tolerance_factor"])
        base_process_settings.AddValue("integration_order", self.mesh_tying_settings["integration_order"])
        base_process_settings.AddValue("search_parameters", self.mesh_tying_settings["search_parameters"])

        # Construct the base process.
        super(MeshTyingProcess, self).__init__(Model, base_process_settings)

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

    def ExecuteInitialize(self):
        """ This method is executed at the begining to initialize the process

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # We call to the base process
        super(MeshTyingProcess, self).ExecuteInitialize()

    def ExecuteBeforeSolutionLoop(self):
        """ This method is executed before starting the time loop

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        # We call to the base process
        super(MeshTyingProcess, self).ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):
        """ This method is executed in order to initialize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        # We call to the base process
        super(MeshTyingProcess, self).ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        """ This method is executed in order to finalize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        # We call to the base process
        super(MeshTyingProcess, self).ExecuteFinalizeSolutionStep()

    def ExecuteBeforeOutputStep(self):
        """ This method is executed right before the ouput process computation

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        # We call to the base process
        super(MeshTyingProcess, self).ExecuteBeforeOutputStep()

    def ExecuteAfterOutputStep(self):
        """ This method is executed right after the ouput process computation

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        # We call to the base process
        super(MeshTyingProcess, self).ExecuteAfterOutputStep()

    def ExecuteFinalize(self):
        """ This method is executed in order to finalize the current computation

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        # We call to the base process
        super(MeshTyingProcess, self).ExecuteFinalize()

    def _get_condition_name(self):
        """ This method returns the condition name

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # We define the condition name to be used
        return "MeshTyingMortar"

    def _get_final_string(self, key = "0"):
        """ This method returns the final string of the condition name

        Keyword arguments:
        self -- It signifies an instance of a class.
        key -- The key to identify the current pair
        """
        # Determine the geometry of the element
        super(MeshTyingProcess, self)._get_final_string(key)
        number_nodes, number_nodes_master = self._compute_number_nodes_elements()
        geometry_element = self.__type_element(number_nodes)
        geometry_element_master = self.__type_element(number_nodes_master)
        # We compute the number of nodes of the conditions
        number_nodes, number_nodes_master = super(MeshTyingProcess, self)._compute_number_nodes()
        if number_nodes != number_nodes_master:
            return geometry_element + str(number_nodes_master) + "N" + geometry_element_master
        else:
            return geometry_element

    def _initialize_search_conditions(self):
        """ This method initializes some conditions values

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # We call to the base process
        super(MeshTyingProcess, self)._initialize_search_conditions()

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

    def _compute_number_nodes_elements(self):
        """ This method computes the  number of nodes per element

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        # We compute the number of nodes of the geometry
        if self.predefined_master_slave and self.dimension == 3:
            slave_defined = False
            master_defined = False
            for elem in self.main_model_part.Elements:
                nodes = elem.GetNodes()
                for node in nodes:
                    if node.Is(KM.SLAVE):
                        number_nodes = len(elem.GetNodes())
                        slave_defined = True
                    if node.Is(KM.MASTER):
                        number_nodes_master = len(elem.GetNodes())
                        master_defined = True
                if slave_defined and master_defined:
                    break
        else:
            number_nodes = len(self.main_model_part.Elements[1].GetNodes())
            number_nodes_master = number_nodes

        return number_nodes, number_nodes_master

    def __type_element(self, number_nodes):
        """ This method computes the type of element considered

        Keyword arguments:
        self -- It signifies an instance of a class.
        number_nodes -- The number of nodes of the element
        """

        if self.dimension == 2:
            if number_nodes == 3:
                return "Triangle"
            elif number_nodes == 4:
                return "Quadrilateral"
        else:
            if number_nodes == 4:
                return "Tetrahedron"
            elif number_nodes == 8:
                return "Hexahedron"
