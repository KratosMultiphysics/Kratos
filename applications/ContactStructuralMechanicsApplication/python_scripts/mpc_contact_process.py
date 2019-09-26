from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics as KM

import KratosMultiphysics.StructuralMechanicsApplication as SMA
import KratosMultiphysics.ContactStructuralMechanicsApplication as CSMA

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return MPCContactProcess(Model, settings["Parameters"])

import sys

import KratosMultiphysics.ContactStructuralMechanicsApplication.search_base_process as search_base_process

class MPCContactProcess(search_base_process.SearchBaseProcess):
    """This class is used in order to compute the contact using a mortar MPC formulation

    This class constructs the model parts containing the contact conditions and
    initializes parameters and variables related with the contact. The class creates
    search utilities to be used to create the contact pairs

    Only the member variables listed below should be accessed directly.

    Public member variables:
    Model -- the container of the different model parts.
    settings -- Kratos parameters containing solver settings.
    """

    __normal_computation = {
        # JSON input
        "NO_DERIVATIVES_COMPUTATION": CSMA.NormalDerivativesComputation.NO_DERIVATIVES_COMPUTATION,
        "no_derivatives_computation": CSMA.NormalDerivativesComputation.NO_DERIVATIVES_COMPUTATION,
        "ELEMENTAL_DERIVATIVES":  CSMA.NormalDerivativesComputation.ELEMENTAL_DERIVATIVES,
        "elemental_derivatives":  CSMA.NormalDerivativesComputation.ELEMENTAL_DERIVATIVES,
        "NODAL_ELEMENTAL_DERIVATIVES": CSMA.NormalDerivativesComputation.NODAL_ELEMENTAL_DERIVATIVES,
        "nodal_elemental_derivatives": CSMA.NormalDerivativesComputation.NODAL_ELEMENTAL_DERIVATIVES,
        "NO_DERIVATIVES_COMPUTATION_WITH_NORMAL_UPDATE": CSMA.NormalDerivativesComputation.NO_DERIVATIVES_COMPUTATION_WITH_NORMAL_UPDATE,
        "no_derivatives_computation_with_normal_update": CSMA.NormalDerivativesComputation.NO_DERIVATIVES_COMPUTATION_WITH_NORMAL_UPDATE
        }

    def __init__(self, Model, settings):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        Model -- the container of the different model parts.
        settings -- Kratos parameters containing solver settings.
        """

        # NOTE: Due to recursive check "contact_model_part" and "assume_master_slave" requires to pre-define configurations, if more that 10 pairs of contact are required, just add. I assume nobody needs that much
        # Settings string in json format
        default_parameters = KM.Parameters("""
        {
            "help"                        : "This class is used in order to compute the contact using a mortar MPC formulation. This class constructs the model parts containing the contact conditions and initializes parameters and variables related with the contact. The class creates search utilities to be used to create the contact pairs",
            "mesh_id"                         : 0,
            "model_part_name"                 : "Structure",
            "computing_model_part_name"       : "computing_domain",
            "contact_model_part"              : {"0":[],"1":[],"2":[],"3":[],"4":[],"5":[],"6":[],"7":[],"8":[],"9":[]},
            "assume_master_slave"             : {"0":[],"1":[],"2":[],"3":[],"4":[],"5":[],"6":[],"7":[],"8":[],"9":[]},
            "contact_property_ids"            : {"0": 0,"1": 0,"2": 0,"3": 0,"4": 0,"5": 0,"6": 0,"7": 0,"8": 0,"9": 0},
            "contact_type"                    : "Frictionless",
            "not_normal_update_frictional"    : false,
            "interval"                        : [0.0,"End"],
            "normal_variation"                : "no_derivatives_computation",
            "frictional_law"                  : "Coulomb",
            "tangent_factor"                  : 1.0e-1,
            "zero_tolerance_factor"           : 1.0e2,
            "reaction_check_stiffness_factor" : 1.0e-10,
            "integration_order"               : 2,
            "clear_inactive_for_post"         : true,
            "update_condition_relation_step"  : false,
            "search_parameters" : {
                "type_search"                         : "in_radius_with_obb",
                "simple_search"                       : false,
                "adapt_search"                        : false,
                "search_factor"                       : 3.5,
                "active_check_factor"                 : 0.01,
                "max_number_results"                  : 1000,
                "bucket_size"                         : 4,
                "dynamic_search"                      : false,
                "static_check_movement"               : false,
                "database_step_update"                : 1,
                "normal_orientation_threshold"        : 1.0e-1,
                "consider_gap_threshold"              : false,
                "debug_mode"                          : false,
                "predict_correct_lagrange_multiplier" : false,
                "check_gap"                           : "check_mapping",
                "octree_search_parameters" : {
                    "bounding_box_factor"             : 0.1,
                    "debug_obb"                       : false,
                    "OBB_intersection_type"           : "SeparatingAxisTheorem",
                    "build_from_bounding_box"         : true,
                    "lower_bounding_box_coefficient"  : 0.0,
                    "higher_bounding_box_coefficient" : 1.0
                }
            }
        }
        """)

        # Overwrite the default settings with user-provided parameters
        self.contact_settings = settings
        self.contact_settings.RecursivelyValidateAndAssignDefaults(default_parameters)

        # We transfer the parameters to the base class
        base_process_settings = KM.Parameters("""{}""")
        base_process_settings.AddValue("mesh_id", self.contact_settings["mesh_id"])
        base_process_settings.AddValue("model_part_name", self.contact_settings["model_part_name"])
        base_process_settings.AddValue("computing_model_part_name", self.contact_settings["computing_model_part_name"])
        base_process_settings.AddValue("search_model_part", self.contact_settings["contact_model_part"])
        base_process_settings.AddValue("assume_master_slave", self.contact_settings["assume_master_slave"])
        base_process_settings.AddValue("search_property_ids", self.contact_settings["contact_property_ids"])
        base_process_settings.AddValue("interval", self.contact_settings["interval"])
        base_process_settings.AddValue("zero_tolerance_factor", self.contact_settings["zero_tolerance_factor"])
        base_process_settings.AddValue("integration_order", self.contact_settings["integration_order"])
        base_process_settings.AddValue("search_parameters", self.contact_settings["search_parameters"])

        # Construct the base process.
        super(MPCContactProcess, self).__init__(Model, base_process_settings)

        # Getting the normal variation flag
        self.normal_variation = super(MPCContactProcess, self)._get_enum_flag(self.contact_settings, "normal_variation", self.__normal_computation)

        # Name of the frictional law
        self.frictional_law = self.contact_settings["frictional_law"].GetString()

        # If we compute a frictional contact simulation
        contact_type = self.contact_settings["contact_type"].GetString()
        if "Frictional" in contact_type:
            not_normal_update_frictional = self.contact_settings["not_normal_update_frictional"].GetBool()
            if not not_normal_update_frictional and not "WithNormalUpdate" in contact_type:
                contact_type += "WithNormalUpdate"
            self.is_frictional = True
            if "PureSlip" in contact_type:
                self.pure_slip = True
            else:
                self.pure_slip = False
            self.mesh_tying = False
        elif "MeshTying" in contact_type:
            self.is_frictional = False
            self.pure_slip = False
            self.mesh_tying = True
        else:
            self.is_frictional = False
            self.pure_slip = False
            self.mesh_tying = False

        # In case we want a normal update
        if "WithNormalUpdate" in contact_type:
            if self.normal_variation == CSMA.NormalDerivativesComputation.NO_DERIVATIVES_COMPUTATION:
                self.normal_variation = CSMA.NormalDerivativesComputation.NO_DERIVATIVES_COMPUTATION_WITH_NORMAL_UPDATE

        # We consider frictional contact (We use the SLIP flag because was the easiest way)
        if self.is_frictional:
            self.computing_model_part.Set(KM.SLIP, True)
        else:
            self.computing_model_part.Set(KM.SLIP, False)

        # We consider mesh tying (We use the RIGID flag because was the easiest way)
        if self.mesh_tying:
            self.computing_model_part.Set(KM.RIGID, True)
        else:
            self.computing_model_part.Set(KM.RIGID, False)

    def ExecuteInitialize(self):
        """ This method is executed at the begining to initialize the process

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # We call to the base process
        super(MPCContactProcess, self).ExecuteInitialize()

    def ExecuteBeforeSolutionLoop(self):
        """ This method is executed before starting the time loop

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        # We call to the base process
        super(MPCContactProcess, self).ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):
        """ This method is executed in order to initialize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # We call to the base process
        super(MPCContactProcess, self).ExecuteInitializeSolutionStep()

        # We set the flag SLIP on frictional conditions
        if self.is_frictional:
            KM.VariableUtils().SetFlag(KM.SLIP, True, self.computing_model_part.Conditions)

    def ExecuteFinalizeSolutionStep(self):
        """ This method is executed in order to finalize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        # We call to the base process
        super(MPCContactProcess, self).ExecuteFinalizeSolutionStep()

        # Debug we compute if the total load corresponds with the total contact force and the reactions
        if self.settings["search_parameters"]["debug_mode"].GetBool():
            total_load = KM.Vector(3)
            total_load[0] = 0.0
            total_load[1] = 0.0
            total_load[2] = 0.0
            total_reaction = KM.Vector(3)
            total_reaction[0] = 0.0
            total_reaction[1] = 0.0
            total_reaction[2] = 0.0
            total_contact_force = 0

            # Computing total load applied (I will consider only surface loads for now)
            for cond in self.computing_model_part.Conditions:
                geom = cond.GetGeometry()
                if cond.Has(SMA.LINE_LOAD):
                    total_load += geom.Length() * cond.GetValue(SMA.LINE_LOAD)
                if cond.Has(SMA.SURFACE_LOAD):
                    total_load += geom.Area() * cond.GetValue(SMA.SURFACE_LOAD)

            for node in self.computing_model_part.Nodes:
                if node.Has(KM.NODAL_AREA) and node.Has(CSMA.AUGMENTED_NORMAL_CONTACT_PRESSURE):
                    total_contact_force += node.GetValue(KM.NODAL_AREA) * node.GetValue(CSMA.AUGMENTED_NORMAL_CONTACT_PRESSURE)

                total_reaction += node.GetSolutionStepValue(KM.REACTION)

            KM.Logger.PrintWarning("TOTAL LOAD: ", "X: {:.2e}".format(total_load[0]), "\t Y: {:.2e}".format(total_load[1]), "\tZ: {:.2e}".format(total_load[2]))
            KM.Logger.PrintWarning("TOTAL REACTION: ", "X: {:.2e}".format(total_reaction[0]), "\t Y: {:.2e}".format(total_reaction[1]), "\tZ: {:.2e}".format(total_reaction[2]))
            KM.Logger.PrintWarning("TOTAL CONTACT FORCE: ", "{:.2e}".format(total_contact_force))

    def ExecuteBeforeOutputStep(self):
        """ This method is executed right before the ouput process computation

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # We call to the base process
        super(MPCContactProcess, self).ExecuteBeforeOutputStep()

    def ExecuteAfterOutputStep(self):
        """ This method is executed right after the ouput process computation

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # We call to the base process
        super(MPCContactProcess, self).ExecuteAfterOutputStep()

    def ExecuteFinalizeSolutionStep(self):
        """ This method is executed in order to finalize the current time step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # We call to the base process
        super(MPCContactProcess, self).ExecuteFinalizeSolutionStep()

        # Blocking the conditions
        if not self.contact_settings["update_condition_relation_step"].GetBool():
            KM.VariableUtils().SetFlag(KM.BLOCKED, True, self.computing_model_part.Conditions)

    def ExecuteFinalize(self):
        """ This method is executed in order to finalize the current computation

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # We call to the base process
        super(MPCContactProcess, self).ExecuteFinalize()

    def _set_additional_parameters(self, param):
        """ This sets additional parameters for the search

        Keyword arguments:
        self -- It signifies an instance of a class.
        param -- The parameters where to set additional values
        """
        param.AddEmptyValue("pure_slip")
        param["pure_slip"].SetBool(self.pure_slip)

    def _get_condition_name(self):
        """ This method returns the condition name

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # We define the condition name to be used
        condition_name = "MPCMortarContact"

        return condition_name

    def _get_final_string(self, key = "0"):
        """ This method returns the final string of the condition name

        Keyword arguments:
        self -- It signifies an instance of a class.
        key -- The key to identify the current pair
        """
        # Determine the geometry of the element
        super(MPCContactProcess, self)._get_final_string(key)
        # We compute the number of nodes of the conditions
        number_nodes, number_nodes_master = super(MPCContactProcess, self)._compute_number_nodes()
        if number_nodes != number_nodes_master:
            return str(number_nodes_master) + "N"
        else:
            return ""

    def _get_problem_name(self):
        """ This method returns the problem name to be solved

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        return "Contact"

    def _initialize_process_info(self):
        """ This method initializes some values from the process info
        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # We call to the base process
        super(MPCContactProcess, self)._initialize_process_info()

        # We call the process info
        process_info = self.main_model_part.ProcessInfo
        # Initialize ACTIVE_SET_CONVERGED
        process_info[CSMA.ACTIVE_SET_CONVERGED] = True
        process_info[CSMA.REACTION_CHECK_STIFFNESS_FACTOR] = self.contact_settings["reaction_check_stiffness_factor"].GetDouble()

    def _initialize_search_values(self):
        """ This method initializes some values and variables used during contact computations

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # We call to the base process
        super(MPCContactProcess, self)._initialize_search_values()

        # We set the CONTACT flag
        self.computing_model_part.Set(KM.CONTACT, True)
        self._get_process_model_part().Set(KM.CONTACT, True)
        # We consider frictional contact (We use the SLIP flag because was the easiest way)
        if self.is_frictional:
            self.computing_model_part.Set(KM.SLIP, True)
            self._get_process_model_part().Set(KM.SLIP, True)
        else:
            self.computing_model_part.Set(KM.SLIP, False)
            self._get_process_model_part().Set(KM.SLIP, False)

        # We call the process info
        process_info = self.main_model_part.ProcessInfo

        # We set the value that scales in the tangent direction the penalty and scale parameter
        if self.is_frictional:
            process_info[KM.TANGENT_FACTOR] = self.contact_settings["tangent_factor"].GetDouble()

    def _initialize_problem_parameters(self):
        """ This method initializes the MPC parameters from the process info

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # We call to the base process
        super(MPCContactProcess, self)._initialize_problem_parameters()

        # We call the process info
        process_info = self.main_model_part.ProcessInfo

    def _initialize_search_conditions(self):
        """ This method initializes some conditions values

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # We call to the base process
        super(MPCContactProcess, self)._initialize_search_conditions()

        alm_init_var = CSMA.ALMFastInit(self._get_process_model_part())
        alm_init_var.Execute()

    def _create_main_search(self, key = "0"):
        """ This method creates the search process that will be use during contact search

        Keyword arguments:
        self -- It signifies an instance of a class.
        key -- The key to identify the current pair
        """
        # Create main parameters
        search_parameters = self._create_search_parameters(key)

        # We create the search process
        self.search_utility_list[key] = CSMA.MPCContactSearchProcess(self.computing_model_part, search_parameters)
