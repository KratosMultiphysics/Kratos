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
    return ALMContactProcess(Model, settings["Parameters"])

import sys

import search_base_process

class ALMContactProcess(search_base_process.SearchBaseProcess):
    """This class is used in order to compute the contact using a mortar ALM formulation

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
            "help"                        : "This class is used in order to compute the contact using a mortar ALM formulation. This class constructs the model parts containing the contact conditions and initializes parameters and variables related with the contact. The class creates search utilities to be used to create the contact pairs",
            "mesh_id"                     : 0,
            "model_part_name"             : "Structure",
            "computing_model_part_name"   : "computing_domain",
            "contact_model_part"          : {"0":[],"1":[],"2":[],"3":[],"4":[],"5":[],"6":[],"7":[],"8":[],"9":[]},
            "assume_master_slave"         : {"0":[],"1":[],"2":[],"3":[],"4":[],"5":[],"6":[],"7":[],"8":[],"9":[]},
            "contact_property_ids"        : {"0": 0,"1": 0,"2": 0,"3": 0,"4": 0,"5": 0,"6": 0,"7": 0,"8": 0,"9": 0},
            "contact_type"                : "Frictionless",
            "interval"                    : [0.0,"End"],
            "normal_variation"            : "no_derivatives_computation",
            "frictional_law"              : "Coulomb",
            "tangent_factor"              : 1.0e-1,
            "integration_order"           : 2,
            "clear_inactive_for_post"     : true,
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
            },
            "advance_ALM_parameters" : {
                "manual_ALM"                  : false,
                "stiffness_factor"            : 1.0,
                "penalty_scale_factor"        : 1.0,
                "use_scale_factor"            : true,
                "penalty"                     : 1.0e-12,
                "scale_factor"                : 1.0e0,
                "adapt_penalty"               : false,
                "max_gap_factor"              : 1.0e-3
            },
            "alternative_formulations" : {
                "axisymmetric"                : false
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
        base_process_settings.AddValue("integration_order", self.contact_settings["integration_order"])
        base_process_settings.AddValue("search_parameters", self.contact_settings["search_parameters"])

        # Construct the base process.
        super(ALMContactProcess, self).__init__(Model, base_process_settings)

        # A check necessary for axisymmetric cases (the domain can not be 3D)
        if (self.contact_settings["alternative_formulations"]["axisymmetric"].GetBool() is True) and (self.dimension == 3):
            raise NameError("3D and axisymmetric makes no sense")

        # Getting the normal variation flag
        self.normal_variation = super(ALMContactProcess, self)._get_enum_flag(self.contact_settings, "normal_variation", self.__normal_computation)

        # Name of the frictional law
        self.frictional_law = self.contact_settings["frictional_law"].GetString()

        # If we compute a frictional contact simulation
        if self.contact_settings["contact_type"].GetString() == "Frictional":
            self.is_frictional = True
            if self.normal_variation == CSMA.NormalDerivativesComputation.NO_DERIVATIVES_COMPUTATION:
                self.normal_variation = CSMA.NormalDerivativesComputation.NO_DERIVATIVES_COMPUTATION_WITH_NORMAL_UPDATE
        else:
            self.is_frictional = False

    def ExecuteInitialize(self):
        """ This method is executed at the begining to initialize the process

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # We call to the base process
        super(ALMContactProcess, self).ExecuteInitialize()

    def ExecuteBeforeSolutionLoop(self):
        """ This method is executed before starting the time loop

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        # We call to the base process
        super(ALMContactProcess, self).ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):
        """ This method is executed in order to initialize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # We call to the base process
        super(ALMContactProcess, self).ExecuteInitializeSolutionStep()

        # Before computing we reset the flags of slip
        if(self._get_if_is_interval() is True):
            if (self.is_frictional is True):
                KM.VariableUtils().SetFlag(KM.SLIP, False, self._get_process_model_part().Nodes)

    def ExecuteFinalizeSolutionStep(self):
        """ This method is executed in order to finalize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        # We call to the base process
        super(ALMContactProcess, self).ExecuteFinalizeSolutionStep()

    def ExecuteBeforeOutputStep(self):
        """ This method is executed right before the ouput process computation

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # We call to the base process
        super(ALMContactProcess, self).ExecuteBeforeOutputStep()

        if (self.contact_settings["clear_inactive_for_post"].GetBool()):
            zero_vector = KM.Array3()
            zero_vector[0] = 0.0
            zero_vector[1] = 0.0
            zero_vector[2] = 0.0
            KM.VariableUtils().SetNonHistoricalVariable(CSMA.AUGMENTED_NORMAL_CONTACT_PRESSURE, 0.0, self.computing_model_part.Nodes, KM.ACTIVE, False)
            KM.VariableUtils().SetNonHistoricalVariable(CSMA.AUGMENTED_TANGENT_CONTACT_PRESSURE, zero_vector, self.computing_model_part.Nodes, KM.ACTIVE, False)

    def ExecuteAfterOutputStep(self):
        """ This method is executed right after the ouput process computation

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # We call to the base process
        super(ALMContactProcess, self).ExecuteAfterOutputStep()

    def ExecuteFinalize(self):
        """ This method is executed in order to finalize the current computation

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # We call to the base process
        super(ALMContactProcess, self).ExecuteFinalize()

    def _get_condition_name(self):
        """ This method returns the condition name

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # We define the condition name to be used
        if self.contact_settings["contact_type"].GetString() == "Frictionless":
            if self.normal_variation == CSMA.NormalDerivativesComputation.NODAL_ELEMENTAL_DERIVATIVES:
                if self.contact_settings["alternative_formulations"]["axisymmetric"].GetBool() is True:
                    condition_name = "ALMNVFrictionlessAxisymMortarContact"
                else:
                    condition_name = "ALMNVFrictionlessMortarContact"
            else:
                if self.contact_settings["alternative_formulations"]["axisymmetric"].GetBool() is True:
                    condition_name = "ALMFrictionlessAxisymMortarContact"
                else:
                    condition_name = "ALMFrictionlessMortarContact"
        elif self.contact_settings["contact_type"].GetString() == "FrictionlessComponents":
            if self.normal_variation == CSMA.NormalDerivativesComputation.NODAL_ELEMENTAL_DERIVATIVES:
                condition_name = "ALMNVFrictionlessComponentsMortarContact"
            else:
                condition_name = "ALMFrictionlessComponentsMortarContact"
        elif self.is_frictional is True:
            if self.normal_variation == CSMA.NormalDerivativesComputation.NODAL_ELEMENTAL_DERIVATIVES:
                if self.contact_settings["alternative_formulations"]["axisymmetric"].GetBool() is True:
                    condition_name = "ALMNVFrictionalAxisymMortarContact"
                else:
                    condition_name = "ALMNVFrictionalMortarContact"
            else:
                if self.contact_settings["alternative_formulations"]["axisymmetric"].GetBool() is True:
                    condition_name = "ALMFrictionalAxisymMortarContact"
                else:
                    condition_name = "ALMFrictionalMortarContact"

        return condition_name

    def _get_final_string(self, key = "0"):
        """ This method returns the final string of the condition name

        Keyword arguments:
        self -- It signifies an instance of a class.
        key -- The key to identify the current pair
        """
        # Determine the geometry of the element
        super(ALMContactProcess, self)._get_final_string(key)
        # We compute the number of nodes of the conditions
        number_nodes, number_nodes_master = super(ALMContactProcess, self)._compute_number_nodes()
        if (number_nodes != number_nodes_master):
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
        super(ALMContactProcess, self)._initialize_process_info()

        # We call the process info
        process_info = self.main_model_part.ProcessInfo
        # We recompute the normal at each iteration (false by default)
        process_info[CSMA.CONSIDER_NORMAL_VARIATION] = self.normal_variation
        # Initialize ACTIVE_SET_CONVERGED
        process_info[CSMA.ACTIVE_SET_CONVERGED] = True
        # We set the max gap factor for the gap adaptation
        max_gap_factor = self.contact_settings["advance_ALM_parameters"]["max_gap_factor"].GetDouble()
        process_info[CSMA.ADAPT_PENALTY] = self.contact_settings["advance_ALM_parameters"]["adapt_penalty"].GetBool()
        process_info[CSMA.MAX_GAP_FACTOR] = max_gap_factor

    def _initialize_search_values(self):
        """ This method initializes some values and variables used during contact computations

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # We call to the base process
        super(ALMContactProcess, self)._initialize_search_values()

        # We set the CONTACT flag
        self.computing_model_part.Set(KM.CONTACT, True)
        self._get_process_model_part().Set(KM.CONTACT, True)
        # We consider frictional contact (We use the SLIP flag because was the easiest way)
        if self.is_frictional is True:
            self.computing_model_part.Set(KM.SLIP, True)
            self._get_process_model_part().Set(KM.SLIP, True)
        else:
            self.computing_model_part.Set(KM.SLIP, False)
            self._get_process_model_part().Set(KM.SLIP, False)

        # We call the process info
        process_info = self.main_model_part.ProcessInfo

        # We recompute the normal at each iteration (false by default)
        process_info[CSMA.CONSIDER_NORMAL_VARIATION] = self.normal_variation

        # We set the value that scales in the tangent direction the penalty and scale parameter
        if self.is_frictional is True:
            process_info[KM.TANGENT_FACTOR] = self.contact_settings["tangent_factor"].GetDouble()

    def _initialize_problem_parameters(self):
        """ This method initializes the ALM parameters from the process info

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # We call to the base process
        super(ALMContactProcess, self)._initialize_problem_parameters()

        # We call the process info
        process_info = self.main_model_part.ProcessInfo

        if (self.contact_settings["advance_ALM_parameters"]["manual_ALM"].GetBool() is False):
            # Computing the scale factors or the penalty parameters (StiffnessFactor * E_mean/h_mean)
            alm_var_parameters = KM.Parameters("""{}""")
            alm_var_parameters.AddValue("stiffness_factor", self.contact_settings["advance_ALM_parameters"]["stiffness_factor"])
            alm_var_parameters.AddValue("penalty_scale_factor", self.contact_settings["advance_ALM_parameters"]["penalty_scale_factor"])
            self.alm_var_process = CSMA.ALMVariablesCalculationProcess(self._get_process_model_part(), KM.NODAL_H, alm_var_parameters)
            self.alm_var_process.Execute()
            # We don't consider scale factor
            if (self.contact_settings["advance_ALM_parameters"]["use_scale_factor"].GetBool() is False):
                process_info[KM.SCALE_FACTOR] = 1.0
        else:
            # We set the values in the process info
            process_info[KM.INITIAL_PENALTY] = self.contact_settings["advance_ALM_parameters"]["penalty"].GetDouble()
            process_info[KM.SCALE_FACTOR] = self.contact_settings["advance_ALM_parameters"]["scale_factor"].GetDouble()

        # We set a minimum value
        if (process_info[KM.INITIAL_PENALTY] < sys.float_info.epsilon):
            process_info[KM.INITIAL_PENALTY] = 1.0e0
        if (process_info[KM.SCALE_FACTOR] < sys.float_info.epsilon):
            process_info[KM.SCALE_FACTOR] = 1.0e0

        # We print the parameters considered
        KM.Logger.PrintInfo("SCALE_FACTOR: ", "{:.2e}".format(process_info[KM.SCALE_FACTOR]))
        KM.Logger.PrintInfo("INITIAL_PENALTY: ", "{:.2e}".format(process_info[KM.INITIAL_PENALTY]))

    def _initialize_search_conditions(self):
        """ This method initializes some conditions values

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # We call to the base process
        super(ALMContactProcess, self)._initialize_search_conditions()

        alm_init_var = CSMA.ALMFastInit(self._get_process_model_part())
        alm_init_var.Execute()
