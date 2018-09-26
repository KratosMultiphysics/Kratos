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
    return SearchBaseProcess(Model, settings["Parameters"])

import sys

# All the processes python processes should be derived from "Process"

class SearchBaseProcess(KM.Process):
    """This class is a base class used to perform the search for contact and mesh tying

    This class constructs the model parts containing the conditions. The class creates
    search utilities to be used to create the pairs

    Only the member variables listed below should be accessed directly.

    Public member variables:
    Model -- the container of the different model parts.
    settings -- Kratos parameters containing solver settings.
    """

    def __init__(self, Model, settings):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        Model -- the container of the different model parts.
        settings -- Kratos parameters containing solver settings.
        """

        # NOTE: Due to recursive check "search_model_part" and "assume_master_slave" requires to pre-define configurations, if more that 10 pairs of contact are required, just add. I assume nobody needs that much
        # Settings string in json format
        default_parameters = KM.Parameters("""
        {
            "help"                        : "This class is a base class used to perform the search for contact and mesh tying. This class constructs the model parts containing the conditions. The class creates search utilities to be used to create the pairs",
            "mesh_id"                     : 0,
            "model_part_name"             : "Structure",
            "computing_model_part_name"   : "computing_domain",
            "search_model_part"           : {"0":[],"1":[],"2":[],"3":[],"4":[],"5":[],"6":[],"7":[],"8":[],"9":[]},
            "assume_master_slave"         : {"0":[],"1":[],"2":[],"3":[],"4":[],"5":[],"6":[],"7":[],"8":[],"9":[]},
            "search_property_ids"         : {"0": 0,"1": 0,"2": 0,"3": 0,"4": 0,"5": 0,"6": 0,"7": 0,"8": 0,"9": 0},
            "interval"                    : [0.0,"End"],
            "integration_order"           : 2,
            "search_parameters" : {
                "type_search"                 : "in_radius",
                "adapt_search"                : false,
                "search_factor"               : 3.5,
                "active_check_factor"         : 0.01,
                "max_number_results"          : 1000,
                "bucket_size"                 : 4,
                "dynamic_search"              : false,
                "database_step_update"        : 1,
                "debug_mode"                  : false,
                "check_gap"                   : "check_mapping"
            }
        }
        """)

        # Overwrite the default settings with user-provided parameters
        self.settings = settings
        self.settings.RecursivelyValidateAndAssignDefaults(default_parameters)

        # The main model part
        self.model = Model
        self.main_model_part = self.model[self.settings["model_part_name"].GetString()]
        # The computing model part
        computing_model_part_name = self.settings["computing_model_part_name"].GetString()
        self.computing_model_part = self.main_model_part.GetSubModelPart(computing_model_part_name)

        self.dimension = self.main_model_part.ProcessInfo[KM.DOMAIN_SIZE]

        self.database_step = 0

        # Detect "End" as a tag and replace it by a large number
        if(self.settings.Has("interval")):
            if(self.settings["interval"][1].IsString() ):
                if(self.settings["interval"][1].GetString() == "End"):
                    self.settings["interval"][1].SetDouble(1e30) # = default_settings["interval"][1]
                else:
                    raise Exception("the second value of interval can be \"End\" or a number, interval currently:"+settings["interval"].PrettyPrintJsonString())

        # Assign this here since it will change the "interval" prior to validation
        self.interval = KM.IntervalUtility(self.settings)

    def ExecuteInitialize(self):
        """ This method is executed at the begining to initialize the process

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # First we generate or identify the different model parts
        if (self.computing_model_part.HasSubModelPart("Contact")):
            self.preprocess = False
            # We get the submodelpart
            self.search_model_part = self.computing_model_part.GetSubModelPart("Contact")
        else:
            self.preprocess = True
            # We create the submodelpart
            self.search_model_part = self.computing_model_part.CreateSubModelPart("Contact")

        # In case of no "Contact" model part we create it
        if (self.preprocess is True):
            # In case no model part is assigned we detect the skin
            self.count_search_model_part = 0
            for key in self.settings["search_model_part"].keys():
                if (self.settings["search_model_part"][key].size() > 0):
                    self.count_search_model_part += 1
            if (self.count_search_model_part == 0):
                self.__detect_skin(self.main_model_part)
            else:
                for key in self.settings["search_model_part"].keys():
                    if (self.settings["search_model_part"][key].size() > 0):
                        self.__generate_search_model_part_from_input_list(self.settings["search_model_part"][key], key)

        # We compute NODAL_H that can be used in the search and some values computation
        self.find_nodal_h = KM.FindNodalHProcess(self.computing_model_part)
        self.find_nodal_h.Execute()

        ## We recompute the search factor and the check in function of the relative size of the mesh
        if (self.settings["search_parameters"]["adapt_search"].GetBool() is True):
            factor = CSMA.ContactUtilities.CalculateRelativeSizeMesh(self.computing_model_part)
            KM.Logger.PrintWarning("SEARCH ADAPT FACTOR: ", "{:.2e}".format(factor))
            search_factor = self.settings["search_parameters"]["search_factor"].GetDouble() * factor
            self.settings["search_parameters"]["search_factor"].SetDouble(search_factor)
            active_check_factor = self.settings["search_parameters"]["active_check_factor"].GetDouble() * factor
            self.settings["search_parameters"]["active_check_factor"].SetDouble(active_check_factor)

        # We call the process info
        self._initialize_process_info()

        #If the conditions doesn't exist we create them
        if (self.preprocess is False):
            master_slave_process = CSMA.MasterSlaveProcess(self.computing_model_part)
            master_slave_process.Execute()

        # Setting the integration order and active check factor
        for prop in self._get_process_model_part().GetProperties():
            prop[CSMA.INTEGRATION_ORDER_CONTACT] = self.settings["integration_order"].GetInt()

        # We initialize the contact values
        self._initialize_search_values()

        # We initialize some other process values
        self._initialize_problem_parameters()

        # Creating the search
        self.search_search = {}
        for key in self.settings["search_model_part"].keys():
            if (self.settings["search_model_part"][key].size() > 0):
                self._create_main_search(key)

        # We initialize conditions
        self._initialize_search_conditions()

        # We perform search initializations
        for key in self.settings["search_model_part"].keys():
            if (self.settings["search_model_part"][key].size() > 0):
                # We initialize the search utility
                self.search_search[key].CheckContactModelParts()
                self.search_search[key].CreatePointListMortar()
                self.search_search[key].InitializeMortarConditions()

    def ExecuteBeforeSolutionLoop(self):
        """ This method is executed before starting the time loop

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def ExecuteInitializeSolutionStep(self):
        """ This method is executed in order to initialize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        if(self._get_if_is_interval() is True):
            self.database_step += 1
            global_step = self.main_model_part.ProcessInfo[KM.STEP]
            database_step_update = self.settings["search_parameters"]["database_step_update"].GetInt()
            if (self.database_step >= database_step_update or global_step == 1):
                # We unset the flag MARKER (used in the nodes to not deactivate it)
                KM.VariableUtils().SetFlag(KM.MARKER, False, self._get_process_model_part().Nodes)
                # We solve one linear step with a linear strategy if needed
                for key in self.settings["search_model_part"].keys():
                    if (self.settings["search_model_part"][key].size() > 0):
                        # Clear current pairs
                        self.search_search[key].ClearMortarConditions()
                        # Update database
                        self.search_search[key].UpdateMortarConditions()
                        #self.search_search[key].CheckMortarConditions()

                # We unset the flag MARKER (used in the nodes to not deactivate it)
                KM.VariableUtils().SetFlag(KM.MARKER, False, self._get_process_model_part().Nodes)

                # Debug
                if (self.settings["search_parameters"]["debug_mode"].GetBool() is True):
                    self._debug_output(global_step, "", self._get_problem_name())
                    # We compute the total integrated area, for debugging
                    self.__get_integration_area()
        else:
            # We deactivate the conditions and nodes
            KM.VariableUtils().SetFlag(KM.ACTIVE, False, self._get_process_model_part().Nodes)
            KM.VariableUtils().SetFlag(KM.ACTIVE, False, self._get_process_model_part().Conditions)

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
        global_step = self.main_model_part.ProcessInfo[KM.STEP]
        if(self._get_if_is_interval() is True):
            modified = self.main_model_part.Is(KM.MODIFIED)
            database_step_update = self.settings["search_parameters"]["database_step_update"].GetInt()
            if (modified is False and (self.database_step >= database_step_update or global_step == 1)):
                for key in self.settings["search_model_part"].keys():
                    if (self.settings["search_model_part"][key].size() > 0):
                        self.search_search[key].ClearMortarConditions()
                self.database_step = 0

    def ExecuteFinalize(self):
        """ This method is executed in order to finalize the current computation

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def _get_condition_name(self):
        """ This method returns the condition name

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        return ""

    def _get_final_string(self, key = "0"):
        """ This method returns the final string of the condition name

        Keyword arguments:
        self -- It signifies an instance of a class.
        key -- The key to identify the current pair
        """
        self.__assume_master_slave(key)
        return ""

    def _get_problem_name(self):
        """ This method returns the problem name to be solved

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        return ""

    def _assign_master_flags(self, partial_model_part):
        """ This method initializes assigment of the master nodes and conditions

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        if (self.predefined_master_slave is True):
            KM.VariableUtils().SetFlag(KM.SLAVE, False, partial_model_part.Nodes)
            KM.VariableUtils().SetFlag(KM.MASTER, True, partial_model_part.Nodes)

            if (len(partial_model_part.Conditions) > 0):
                KM.VariableUtils().SetFlag(KM.SLAVE, False, partial_model_part.Conditions)
                KM.VariableUtils().SetFlag(KM.MASTER, True, partial_model_part.Conditions)

    def _assign_slave_flags(self, key = "0"):
        """ This method initializes assigment of the slave nodes and conditions

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        if (self.predefined_master_slave is True):
            if not self.settings["assume_master_slave"][key].IsArray():
                raise Exception("{0} Error: Model part list is unreadable".format(self.__class__.__name__))
            for i in range(0, self.settings["assume_master_slave"][key].size()):
                model_part_slave_name = self.settings["assume_master_slave"][key][i].GetString()
                model_part_slave = self.main_model_part.GetSubModelPart(model_part_slave_name)
                KM.VariableUtils().SetFlag(KM.SLAVE, True, model_part_slave.Nodes)
                KM.VariableUtils().SetFlag(KM.MASTER, False, model_part_slave.Nodes)

                if (len(model_part_slave.Conditions) > 0):
                    KM.VariableUtils().SetFlag(KM.SLAVE, True, model_part_slave.Conditions)
                    KM.VariableUtils().SetFlag(KM.MASTER, False, model_part_slave.Conditions)

    def _interface_preprocess(self, partial_model_part, key = "0"):
        """ This method creates the process used to compute the contact interface

        Keyword arguments:
        self -- It signifies an instance of a class.
        partial_model_part -- The partial model part that contains the structural problem to be solved
        """

        # We create the process for creating the interface
        self.interface_preprocess = CSMA.InterfacePreprocessCondition(self.computing_model_part)

        # It should create the conditions automatically
        interface_parameters = KM.Parameters("""{"simplify_geometry": false, "contact_property_id": 0}""")
        interface_parameters["contact_property_id"].SetInt(self.settings["search_property_ids"][key].GetInt())
        if (self.dimension == 2):
            self.interface_preprocess.GenerateInterfacePart2D(partial_model_part, interface_parameters)
        else:
            self.interface_preprocess.GenerateInterfacePart3D(partial_model_part, interface_parameters)

    def _initialize_process_info(self):
        """ This method initializes some values from the process info
        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # We call the process info
        process_info = self.main_model_part.ProcessInfo
        process_info[CSMA.ACTIVE_CHECK_FACTOR] = self.settings["search_parameters"]["active_check_factor"].GetDouble()

    def _initialize_search_values(self):
        """ This method initializes some values and variables used during search computations

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # We call the process info
        process_info = self.main_model_part.ProcessInfo

        # We set the distance threshold by default
        process_info[CSMA.DISTANCE_THRESHOLD] = 1.0e24

        # Setting the integration order and active check factor
        for prop in self._get_process_model_part().GetProperties():
            prop[CSMA.INTEGRATION_ORDER_CONTACT] = self.settings["integration_order"].GetInt()
            prop[CSMA.ACTIVE_CHECK_FACTOR] = self.settings["search_parameters"]["active_check_factor"].GetDouble()

    def _initialize_problem_parameters(self):
        """ This method initializes some values and variables used during problem computations

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def _initialize_search_conditions(self):
        """ This method initializes some conditions values

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def _create_main_search(self, key = "0"):
        """ This method creates the search process that will be use during contact search

        Keyword arguments:
        self -- It signifies an instance of a class.
        key -- The key to identify the current pair
        """

        search_parameters = KM.Parameters("""{"condition_name": "", "final_string": "", "predefined_master_slave" : true, "id_name" : ""}""")
        search_parameters.AddValue("type_search", self.settings["search_parameters"]["type_search"])
        search_parameters.AddValue("check_gap", self.settings["search_parameters"]["check_gap"])
        search_parameters.AddValue("allocation_size", self.settings["search_parameters"]["max_number_results"])
        search_parameters.AddValue("bucket_size", self.settings["search_parameters"]["bucket_size"])
        search_parameters.AddValue("search_factor", self.settings["search_parameters"]["search_factor"])
        search_parameters.AddValue("dynamic_search", self.settings["search_parameters"]["dynamic_search"])
        search_parameters["condition_name"].SetString(self._get_condition_name())
        search_parameters["final_string"].SetString(self._get_final_string())
        self.__assume_master_slave(key)
        search_parameters["predefined_master_slave"].SetBool(self.predefined_master_slave)
        search_parameters["id_name"].SetString(key)

        # We compute the number of nodes of the geometry
        number_nodes, number_nodes_master = self._compute_number_nodes()

        # We create the search process
        if (self.dimension == 2):
            self.search_search[key] = CSMA.TreeContactSearch2D2N(self.computing_model_part, search_parameters)
        else:
            if (number_nodes == 3):
                if (number_nodes_master == 3):
                    self.search_search[key] = CSMA.TreeContactSearch3D3N(self.computing_model_part, search_parameters)
                else:
                    self.search_search[key] = CSMA.TreeContactSearch3D3N4N(self.computing_model_part, search_parameters)

            elif (number_nodes == 4):
                if (number_nodes_master == 3):
                    self.search_search[key] = CSMA.TreeContactSearch3D4N3N(self.computing_model_part, search_parameters)
                else:
                    self.search_search[key] = CSMA.TreeContactSearch3D4N(self.computing_model_part, search_parameters)
            else:
                raise Exception("Geometries not compatible. Check all the geometries are linear")

    def _get_enum_flag(self, param, label, dictionary):
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

    def _debug_output(self, label, name, problem_name = ""):
        """ This method is used for debugging pourposes, it creates a postprocess file when called, in sucha  way that it can.

        Keyword arguments:
        self -- It signifies an instance of a class.
        label -- The label to add to the postprocess file
        name -- The name to append to the file
        """
        output_file = "POSTSEARCH"
        gid_mode = KM.GiDPostMode.GiD_PostBinary
        singlefile = KM.MultiFileFlag.SingleFile
        deformed_mesh_flag = KM.WriteDeformedMeshFlag.WriteUndeformed
        write_conditions = KM.WriteConditionsFlag.WriteConditions

        gid_io = KM.GidIO(output_file + name + "_STEP_" + str(label), gid_mode, singlefile, deformed_mesh_flag, write_conditions)

        gid_io.InitializeMesh(label)
        gid_io.WriteMesh(self.main_model_part.GetMesh())
        gid_io.FinalizeMesh()
        gid_io.InitializeResults(label, self.main_model_part.GetMesh())

        gid_io.WriteNodalFlags(KM.INTERFACE, "INTERFACE", self.main_model_part.Nodes, label)
        gid_io.WriteNodalFlags(KM.ACTIVE, "ACTIVE", self.main_model_part.Nodes, label)
        gid_io.WriteNodalFlags(KM.ISOLATED, "ISOLATED", self.main_model_part.Nodes, label)
        gid_io.WriteNodalFlags(KM.SLAVE, "SLAVE", self.main_model_part.Nodes, label)
        gid_io.WriteNodalResults(KM.NORMAL, self.main_model_part.Nodes, label, 0)
        gid_io.WriteNodalResultsNonHistorical(KM.NODAL_AREA, self.main_model_part.Nodes, label)
        gid_io.WriteNodalResults(KM.DISPLACEMENT, self.main_model_part.Nodes, label, 0)
        if (self.main_model_part.Nodes[1].SolutionStepsDataHas(KM.VELOCITY_X) is True):
            gid_io.WriteNodalResults(KM.VELOCITY, self.main_model_part.Nodes, label, 0)
            gid_io.WriteNodalResults(KM.ACCELERATION, self.main_model_part.Nodes, label, 0)
        gid_io.WriteNodalResultsNonHistorical(CSMA.NORMAL_GAP, self.main_model_part.Nodes, label)
        gid_io.WriteNodalResultsNonHistorical(CSMA.AUXILIAR_COORDINATES, self.main_model_part.Nodes, label)
        gid_io.WriteNodalResultsNonHistorical(CSMA.DELTA_COORDINATES, self.main_model_part.Nodes, label)

        if (problem_name == "Contact"):
            gid_io.WriteNodalResults(CSMA.WEIGHTED_GAP, self.main_model_part.Nodes, label, 0)
            gid_io.WriteNodalResultsNonHistorical(CSMA.AUGMENTED_NORMAL_CONTACT_PRESSURE, self.main_model_part.Nodes, label)
            if (self.is_frictional is True):
                gid_io.WriteNodalFlags(KM.SLIP, "SLIP", self.main_model_part.Nodes, label)
                gid_io.WriteNodalResultsNonHistorical(CSMA.AUGMENTED_TANGENT_CONTACT_PRESSURE, self.main_model_part.Nodes, label)
            if (self.main_model_part.Nodes[1].SolutionStepsDataHas(CSMA.LAGRANGE_MULTIPLIER_CONTACT_PRESSURE) is True):
                gid_io.WriteNodalResults(CSMA.LAGRANGE_MULTIPLIER_CONTACT_PRESSURE, self.main_model_part.Nodes, label, 0)
            else:
                gid_io.WriteNodalResults(KM.VECTOR_LAGRANGE_MULTIPLIER, self.main_model_part.Nodes, label, 0)

        gid_io.FinalizeResults()

        #raise NameError("DEBUG")

    def _get_process_model_part(self):
        """ Returns the search model part considered

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        return self.search_model_part

    def _get_if_is_interval(self):
        """ Returns if we are inside the time interval or not

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        current_time = self.main_model_part.ProcessInfo[KM.TIME]
        if(self.interval.IsInInterval(current_time)):
            return True
        else:
            return False

    def _compute_number_nodes(self):
        # We compute the number of nodes of the geometry
        if (self.predefined_master_slave is True and self.dimension == 3):
            slave_defined = False
            master_defined = False
            for cond in self.computing_model_part.Conditions:
                if (cond.Is(KM.SLAVE)):
                    number_nodes = len(cond.GetNodes())
                    slave_defined = True
                if (cond.Is(KM.MASTER)):
                    number_nodes_master = len(cond.GetNodes())
                    master_defined = True
                if (slave_defined and master_defined):
                    break
        else:
            number_nodes = len(self.computing_model_part.Conditions[1].GetNodes())
            number_nodes_master = number_nodes

        return number_nodes, number_nodes_master

    def __get_integration_area(self):
        """ Computes auxiliarly the current integration area

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        # We compute the total integrated area, for debugging
        total_area = 0.0
        if (self.dimension == 2):
            exact_integration = KM.ExactMortarIntegrationUtility2D2N(3)
        else:
            number_nodes, number_nodes_master = self._compute_number_nodes()
            if (number_nodes == 3):
                if (number_nodes_master == 3):
                    exact_integration = KM.ExactMortarIntegrationUtility3D3N(3)
                else:
                    exact_integration = KM.ExactMortarIntegrationUtility3D3N4N(3)
            else:
                if (number_nodes_master == 3):
                    exact_integration = KM.ExactMortarIntegrationUtility3D4N3N(3)
                else:
                    exact_integration = KM.ExactMortarIntegrationUtility3D4N(3)

        # We iterate over the conditions
        for cond in self._get_process_model_part().Conditions:
            if cond.Is(KM.SLAVE):
                area = exact_integration.TestGetExactAreaIntegration(self._get_process_model_part(), cond)
                total_area += area

        KM.Logger.PrintWarning("TOTAL INTEGRATED AREA: ", "{:.2e}".format(total_area))

        #exact_integration.TestGiDDebug(self._get_process_model_part())

    def __generate_search_model_part_from_input_list(self, param, key = "0"):
        """ Generates a contact model part from a list of model parts

        Keyword arguments:
        self -- It signifies an instance of a class.
        param -- The configuration parameters
        key -- The key to identify the current pair
        """

        sub_search_model_part_name = "ContactSub"+key
        sub_search_model_part = self._get_process_model_part().CreateSubModelPart(sub_search_model_part_name)

        # If we assume master/slave pairs or not
        self.__assume_master_slave(key)

        # At least verify that the input is a string
        if not param.IsArray():
            raise Exception("{0} Error: Model part list is unreadable".format(self.__class__.__name__))

        # We transfer the list of submodelparts to the contact model part
        for i in range(0, param.size()):
            model_part_name = param[i].GetString()
            partial_model_part = self.main_model_part.GetSubModelPart(model_part_name)

            # Assigning master and slave sides
            self._assign_master_flags(partial_model_part)
            self._assign_slave_flags(key)

            # We set the interface flag
            KM.VariableUtils().SetFlag(KM.INTERFACE, True, partial_model_part.Nodes)
            if (len(partial_model_part.Conditions) == 0):
                KM.Logger.PrintInfo("Contact Process", "Using nodes for interface. We recommend to use conditions instead")
            else:
                KM.VariableUtils().SetFlag(KM.INTERFACE, True, partial_model_part.Conditions)

        if (self.preprocess is True):
            id_prop = self.settings["search_property_ids"][key].GetInt()
            if (id_prop != 0):
                sub_search_model_part.SetProperties(self.main_model_part.GetProperties(id_prop))
            else:
                sub_search_model_part.SetProperties(self.main_model_part.GetProperties())

            # We transfer the list of submodelparts to the contact model part
            for i in range(0, param.size()):
                partial_model_part = self.main_model_part.GetSubModelPart(param[i].GetString())

                if (self.computing_model_part.Is(KM.MODIFIED)):
                    KM.VariableUtils().SetFlag(KM.TO_ERASE, True, partial_model_part.Conditions)
                    partial_model_part.RemoveConditions(KM.TO_ERASE)

                # We generate the conditions
                self._interface_preprocess(partial_model_part, key)

                # We copy the conditions to the contact model part
                transfer_process = KM.FastTransferBetweenModelPartsProcess(sub_search_model_part, partial_model_part, KM.FastTransferBetweenModelPartsProcess.EntityTransfered.NODESANDCONDITIONS)
                transfer_process.Execute()

    def __detect_skin(self, model_part, key = "0"):
        """ Generates a contact model part from the skin
        Keyword arguments:
        self -- It signifies an instance of a class
        key -- The key to identify the current pair.
        """
        detect_skin_parameters = KM.Parameters("""{"name_auxiliar_model_part": "Contact"}""")
        sub_search_model_part_name = "ContactSub"+key
        self._get_process_model_part().CreateSubModelPart(sub_search_model_part_name)
        detect_skin_parameters["name_auxiliar_model_part"].SetString(sub_search_model_part_name)
        if (self.dimension == 2):
            detect_skin = KM.SkinDetectionProcess2D(model_part, detect_skin_parameters)
        else:
            detect_skin = KM.SkinDetectionProcess3D(model_part, detect_skin_parameters)
        detect_skin.Execute()
        self.settings["search_model_part"][key].Append(sub_search_model_part_name)
        # Assigning master and slave sides
        self.__assume_master_slave(key)
        self._assign_master_flags(self._get_process_model_part())
        self._assign_slave_flags(key)

    def __assume_master_slave(self, key = "0"):
        """ Assigns as true or false if we assume master or slave
        self -- It signifies an instance of a class.
        key -- The key to identify the current pair
        """
        # When all conditions are simultaneously master and slave
        if (self.settings["assume_master_slave"][key].size() > 0):
            self.predefined_master_slave = True
        else:
            self.predefined_master_slave = False

