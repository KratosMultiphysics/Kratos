from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics as KM

import KratosMultiphysics.StructuralMechanicsApplication as SMA
import KratosMultiphysics.ContactStructuralMechanicsApplication as CSMA

def Factory(settings, Model):
    if(type(settings) != KM.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return PenaltyContactProcess(Model, settings["Parameters"])

import sys

import alm_contact_process

class PenaltyContactProcess(alm_contact_process.ALMContactProcess):
    """This class is used in order to compute the contact using a mortar penalty formulation

    This class constructs the model parts containing the contact conditions and
    initializes parameters and variables related with the contact. The class creates
    search utilities to be used to create the contact pairs

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
                "type_search"                         : "in_radius",
                "simple_search"                       : false,
                "adapt_search"                        : false,
                "search_factor"                       : 3.5,
                "active_check_factor"                 : 0.01,
                "max_number_results"                  : 1000,
                "bucket_size"                         : 4,
                "dynamic_search"                      : false,
                "static_check_movement"               : false,
                "database_step_update"                : 1,
                "consider_gap_threshold"              : false,
                "debug_mode"                          : false,
                "predict_correct_lagrange_multiplier" : false,
                "check_gap"                           : "check_mapping"
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
        super(PenaltyContactProcess, self).__init__(Model, base_process_settings)

        # A check necessary for axisymmetric cases (the domain can not be 3D)
        if (self.contact_settings["alternative_formulations"]["axisymmetric"].GetBool() is True) and (self.dimension == 3):
            raise NameError("3D and axisymmetric makes no sense")

        # Getting the normal variation flag
        self.normal_variation = super(PenaltyContactProcess, self)._get_enum_flag(self.contact_settings, "normal_variation", self.__normal_computation)

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
        super(PenaltyContactProcess, self).ExecuteInitialize()

    def ExecuteBeforeSolutionLoop(self):
        """ This method is executed before starting the time loop

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        # We call to the base process
        super(PenaltyContactProcess, self).ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):
        """ This method is executed in order to initialize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # We call to the base process
        super(PenaltyContactProcess, self).ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        """ This method is executed in order to finalize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        # We call to the base process
        super(PenaltyContactProcess, self).ExecuteFinalizeSolutionStep()

    def ExecuteBeforeOutputStep(self):
        """ This method is executed right before the ouput process computation

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # We call to the base process
        super(PenaltyContactProcess, self).ExecuteBeforeOutputStep()

    def ExecuteAfterOutputStep(self):
        """ This method is executed right after the ouput process computation

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # We call to the base process
        super(PenaltyContactProcess, self).ExecuteAfterOutputStep()

    def ExecuteFinalize(self):
        """ This method is executed in order to finalize the current computation

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # We call to the base process
        super(PenaltyContactProcess, self).ExecuteFinalize()

    def _get_condition_name(self):
        """ This method returns the condition name

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # We define the condition name to be used
        if self.contact_settings["contact_type"].GetString() == "Frictionless":
            if self.normal_variation == CSMA.NormalDerivativesComputation.NODAL_ELEMENTAL_DERIVATIVES:
                if self.contact_settings["alternative_formulations"]["axisymmetric"].GetBool() is True:
                    condition_name = "PenaltyNVFrictionlessAxisymMortarContact"
                else:
                    condition_name = "PenaltyNVFrictionlessMortarContact"
            else:
                if self.contact_settings["alternative_formulations"]["axisymmetric"].GetBool() is True:
                    condition_name = "PenaltyFrictionlessAxisymMortarContact"
                else:
                    condition_name = "PenaltyFrictionlessMortarContact"
        elif self.is_frictional is True:
            if self.normal_variation == CSMA.NormalDerivativesComputation.NODAL_ELEMENTAL_DERIVATIVES:
                if self.contact_settings["alternative_formulations"]["axisymmetric"].GetBool() is True:
                    condition_name = "PenaltyNVFrictionalAxisymMortarContact"
                else:
                    condition_name = "PenaltyNVFrictionalMortarContact"
            else:
                if self.contact_settings["alternative_formulations"]["axisymmetric"].GetBool() is True:
                    condition_name = "PenaltyFrictionalAxisymMortarContact"
                else:
                    condition_name = "PenaltyFrictionalMortarContact"

        return condition_name
