from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics as KM

# Check that applications were imported in the main script
KM.CheckRegisteredApplications("StructuralMechanicsApplication")
KM.CheckRegisteredApplications("ContactStructuralMechanicsApplication")

import KratosMultiphysics.StructuralMechanicsApplication as SMA
import KratosMultiphysics.ContactStructuralMechanicsApplication as CSMA

def Factory(settings, Model):
    if(type(settings) != KM.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return MeshTyingProcess(Model, settings["Parameters"])


# All the processes python should be derived from "Process"


class MeshTyingProcess(KM.Process):
    """This class is used in order to compute the a mortar mesh tying formulation

    This class constructs the model parts containing the mesh tying conditions and
    initializes parameters and variables related with the mesh tying. The class creates
    search utilities to be used to create the tying pairs

    Only the member variables listed below should be accessed directly.

    Public member variables:
    model_part -- the model part used to construct the process.
    settings -- Kratos parameters containing solver settings.
    """

    def __init__(self, model_part, settings):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        model_part -- the model part used to construct the process.
        settings -- Kratos parameters containing solver settings.
        """

        # Settings string in json format
        default_parameters = KM.Parameters("""
        {
            "mesh_id"                     : 0,
            "model_part_name"             : "Structure",
            "computing_model_part_name"   : "computing_domain",
            "mesh_tying_model_part"       : "Tying_Part",
            "assume_master_slave"         : "",
            "type_variable"               : "Components",
            "geometry_element"            : "Quadrilateral",
            "search_parameters" : {
                "type_search"                 : "in_radius",
                "search_factor"               : 3.5,
                "active_check_factor"         : 0.01,
                "max_number_results"          : 1000,
                "bucket_size"                 : 4,
                "dynamic_search"              : false,
                "debug_mode"                  : false,
                "check_gap"                   : "check_mapping"
            },
            "integration_order"           : 2
        }
        """)

        # Overwrite the default settings with user-provided parameters
        self.settings = settings
        self.settings.ValidateAndAssignDefaults(default_parameters)

        self.main_model_part = model_part[self.settings["model_part_name"].GetString()]
        self.computing_model_part_name = self.settings["computing_model_part_name"].GetString()

        self.dimension = self.main_model_part.ProcessInfo[KM.DOMAIN_SIZE]

        self.mesh_tying_model_part = model_part[self.settings["mesh_tying_model_part"].GetString()]

        if (self.settings["assume_master_slave"].GetString() != ""):
            KM.VariableUtils().SetFlag(KM.SLAVE, False, self.mesh_tying_model_part.Nodes)
            KM.VariableUtils().SetFlag(KM.MASTER, True, self.mesh_tying_model_part.Nodes)
            model_part_slave = self.main_model_part.GetSubModelPart(self.settings["assume_master_slave"].GetString())
            KM.VariableUtils().SetFlag(KM.SLAVE, True, model_part_slave.Nodes)
            KM.VariableUtils().SetFlag(KM.MASTER, False, model_part_slave.Nodes)

        self.type_variable = self.settings["type_variable"].GetString()
        self.geometry_element = self.settings["geometry_element"].GetString()

    def ExecuteInitialize(self):
        """ This method is executed at the begining to initialize the process

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # The computing model part
        computing_model_part = self.main_model_part.GetSubModelPart(self.computing_model_part_name)

        # We compute NODAL_H that can be used in the search and some values computation
        self.find_nodal_h = KM.FindNodalHProcess(computing_model_part)
        self.find_nodal_h.Execute()

        # Appending the conditions created to the self.main_model_part
        if (computing_model_part.HasSubModelPart("Contact")):
            interface_model_part = computing_model_part.GetSubModelPart("Contact")
        else:
            interface_model_part = computing_model_part.CreateSubModelPart("Contact")

        # Setting the integration order and active check factor
        for prop in computing_model_part.GetProperties():
            prop[CSMA.INTEGRATION_ORDER_CONTACT] = self.settings["integration_order"].GetInt()

        self.main_model_part.ProcessInfo[CSMA.ACTIVE_CHECK_FACTOR] = self.settings["search_parameters"]["active_check_factor"].GetDouble()

        # We set the interface flag
        KM.VariableUtils().SetFlag(KM.INTERFACE, True, self.mesh_tying_model_part.Nodes)

        self.interface_preprocess = CSMA.InterfacePreprocessCondition(self.main_model_part)

        # It should create the conditions automatically
        interface_parameters = KM.Parameters("""{"simplify_geometry": false}""")
        if (self.dimension == 2):
            self.interface_preprocess.GenerateInterfacePart2D(self.mesh_tying_model_part, interface_parameters)
        else:
            self.interface_preprocess.GenerateInterfacePart3D(self.mesh_tying_model_part, interface_parameters)

        # When all conditions are simultaneously master and slave
        if (self.settings["assume_master_slave"].GetString() == ""):
            KM.VariableUtils().SetFlag(KM.SLAVE, True, self.mesh_tying_model_part.Conditions)

        # We copy the conditions to the ContactSubModelPart
        for cond in self.mesh_tying_model_part.Conditions:
            interface_model_part.AddCondition(cond)
        del(cond)
        for node in self.mesh_tying_model_part.Nodes:
            interface_model_part.AddNode(node, 0)
        del(node)

        # Creating the search
        condition_name = "MeshTyingMortar"
        search_parameters = KM.Parameters("""{"condition_name": "", "final_string": ""}""")
        search_parameters.AddValue("type_search",self.settings["search_parameters"]["type_search"])
        search_parameters.AddValue("check_gap",self.settings["search_parameters"]["check_gap"])
        search_parameters.AddValue("allocation_size",self.settings["search_parameters"]["max_number_results"])
        search_parameters.AddValue("bucket_size",self.settings["search_parameters"]["bucket_size"])
        search_parameters.AddValue("search_factor",self.settings["search_parameters"]["search_factor"])
        search_parameters["condition_name"].SetString(condition_name)
        search_parameters["final_string"].SetString(self.geometry_element + self.type_variable)

        # We compute the number of nodes of the geometry
        number_nodes = len(computing_model_part.Conditions[1].GetNodes())
        if (self.dimension == 2):
            self.contact_search = CSMA.TreeContactSearch2D2N(computing_model_part, search_parameters)
        else:
            if (number_nodes == 3):
                self.contact_search = CSMA.TreeContactSearch3D3N(computing_model_part, search_parameters)
            else:
                self.contact_search = CSMA.TreeContactSearch3D4N(computing_model_part, search_parameters)

        zero_vector = KM.Vector(3)
        zero_vector[0] = 0.0
        zero_vector[1] = 0.0
        zero_vector[2] = 0.0

        # Initilialize weighted variables and LM
        if (self.type_variable == "Scalar"):
            KM.VariableUtils().SetVariable(CSMA.WEIGHTED_SCALAR_RESIDUAL, 0.0, self.mesh_tying_model_part.Nodes)
        else:
            KM.VariableUtils().SetVariable(CSMA.WEIGHTED_VECTOR_RESIDUAL, zero_vector, self.mesh_tying_model_part.Nodes)

        # Setting the conditions
        KM.VariableUtils().SetNonHistoricalVariable(KM.NORMAL, zero_vector, self.mesh_tying_model_part.Conditions)

        self.contact_search.CreatePointListMortar()
        self.contact_search.InitializeMortarConditions()

    def ExecuteBeforeSolutionLoop(self):
        """ This method is executed before starting the time loop

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        self.contact_search.ClearMortarConditions()

        self.contact_search.UpdateMortarConditions()
        #self.contact_search.CheckMortarConditions()

        # We initialize the conditions (just in case, technically already done)
        for cond in self.mesh_tying_model_part.Conditions:
            cond.Initialize()
        del cond

    def ExecuteInitializeSolutionStep(self):
        """ This method is executed in order to initialize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        #self.contact_search.CheckMortarConditions()
        pass

    def ExecuteFinalizeSolutionStep(self):
        """ This method is executed in order to finalize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        pass

    def ExecuteBeforeOutputStep(self):
        """ This method is executed right before the ouput process computation

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        pass

    def ExecuteAfterOutputStep(self):
        """ This method is executed right after the ouput process computation

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        pass

    def ExecuteFinalize(self):
        """ This method is executed in order to finalize the current computation

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        pass

    def __get_enum_flag(self, param, label, dictionary):
        """ Parse enums settings using an auxiliary dictionary of acceptable values.

        Keyword arguments:
        self -- It signifies an instance of a class.
        param -- The label to add to the postprocess file
        label -- The label used to get the string
        dictionary -- The dictionary containing the list of possible candidates
        """

        keystring = param[label].GetString()
        try:
            value = dictionary[keystring]
        except KeyError:
            msg = "{0} Error: Unknown value \"{1}\" read for parameter \"{2}\"".format(self.__class__.__name__, value, label)
            raise Exception(msg)

        return value
