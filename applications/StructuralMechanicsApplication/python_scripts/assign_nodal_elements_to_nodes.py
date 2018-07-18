import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
from math import *

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignNodalElementsToNodes(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class AssignNodalElementsToNodes(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
        {
            "mesh_id"                        : 0,
            "main_model_part"                : "Structure",
            "sub_model_part_name"            : "",
            "rayleigh_damping"               : false,
            "assign_active_flag_node"        : true,
            "constitutive_law_name"          : "SpringConstitutiveLaw",
            "additional_dependence_variables": [],
            "interval"                       : [0.0, 1e30]
        }
        """
        )

        to_validate_parameters = KratosMultiphysics.Parameters("""{}""")
        if (settings.Has("mesh_id")):
            to_validate_parameters.AddValue("mesh_id", settings["mesh_id"])
        if (settings.Has("model_part_name")):
            to_validate_parameters.AddValue("model_part_name", settings["model_part_name"])
        if (settings.Has("sub_model_part_name")):
            to_validate_parameters.AddValue("sub_model_part_name", settings["sub_model_part_name"])
        if (settings.Has("rayleigh_damping")):
            to_validate_parameters.AddValue("rayleigh_damping", settings["rayleigh_damping"])
        if (settings.Has("assign_active_flag_node")):
            to_validate_parameters.AddValue("assign_active_flag_node", settings["assign_active_flag_node"])
        if (settings.Has("constitutive_law_name")):
            to_validate_parameters.AddValue("constitutive_law_name", settings["constitutive_law_name"])
        if (settings.Has("additional_dependence_variables")):
            to_validate_parameters.AddValue("additional_dependence_variables", settings["additional_dependence_variables"])
        if (settings.Has("interval")):
            to_validate_parameters.AddValue("interval", settings["interval"])

        # Overwrite the default settings with user-provided parameters
        to_validate_parameters.RecursivelyValidateAndAssignDefaults(default_settings)
        if (settings.Has("mesh_id")):
            settings.SetValue("mesh_id", to_validate_parameters["mesh_id"])
        else:
            settings.AddValue("mesh_id", to_validate_parameters["mesh_id"])
        if (settings.Has("model_part_name")):
            settings.SetValue("model_part_name", to_validate_parameters["model_part_name"])
        else:
            settings.AddValue("model_part_name", to_validate_parameters["model_part_name"])
        if (settings.Has("sub_model_part_name")):
            settings.SetValue("sub_model_part_name", to_validate_parameters["sub_model_part_name"])
        else:
            settings.AddValue("sub_model_part_name", to_validate_parameters["sub_model_part_name"])
        if (settings.Has("rayleigh_damping")):
            settings.SetValue("rayleigh_damping", to_validate_parameters["rayleigh_damping"])
        else:
            settings.AddValue("rayleigh_damping", to_validate_parameters["rayleigh_damping"])
        if (settings.Has("assign_active_flag_node")):
            settings.SetValue("assign_active_flag_node", to_validate_parameters["assign_active_flag_node"])
        else:
            settings.AddValue("assign_active_flag_node", to_validate_parameters["assign_active_flag_node"])
        if (settings.Has("constitutive_law_name")):
            settings.SetValue("constitutive_law_name", to_validate_parameters["constitutive_law_name"])
        else:
            settings.AddValue("constitutive_law_name", to_validate_parameters["constitutive_law_name"])
        if (settings.Has("additional_dependence_variables")):
            settings.SetValue("additional_dependence_variables", to_validate_parameters["additional_dependence_variables"])
        else:
            settings.AddValue("additional_dependence_variables", to_validate_parameters["additional_dependence_variables"])
        if (settings.Has("interval")):
            settings.SetValue("interval", to_validate_parameters["interval"])
        else:
            settings.AddValue("interval", to_validate_parameters["interval"])

        # List of auxiliar parameters to assign in case not defined
        auxiliar_parameters = KratosMultiphysics.Parameters("""
        {
            "nodal_mass"                     : null,
            "nodal_inertia"                  : [null, null, null],
            "nodal_stiffness"                : [null, null, null],
            "nodal_rotational_stiffness"     : [null, null, null],
            "nodal_damping_ratio"            : [null, null, null],
            "nodal_rotational_damping_ratio" : [null, null, null]
        }
        """
        )

        if not (settings.Has("nodal_mass")):
            settings.AddValue("nodal_mass", auxiliar_parameters["nodal_mass"])
        if not (settings.Has("nodal_inertia")):
            settings.AddValue("nodal_inertia", auxiliar_parameters["nodal_inertia"])
        if not (settings.Has("nodal_stiffness")):
            settings.AddValue("nodal_stiffness", auxiliar_parameters["nodal_stiffness"])
        if not (settings.Has("nodal_rotational_stiffness")):
            settings.AddValue("nodal_rotational_stiffness", auxiliar_parameters["nodal_rotational_stiffness"])
        if not (settings.Has("nodal_damping_ratio")):
            settings.AddValue("nodal_damping_ratio", auxiliar_parameters["nodal_damping_ratio"])
        if not (settings.Has("nodal_rotational_damping_ratio")):
            settings.AddValue("nodal_rotational_damping_ratio", auxiliar_parameters["nodal_rotational_damping_ratio"])

        # Finally we assign to the instance
        self.settings = settings

        # The main model part
        self.model = Model
        self.main_model_part = self.model[self.settings["model_part_name"].GetString()]

        # The creation of the process
        process_parameters = KratosMultiphysics.Parameters("""{}""")
        process_parameters.AddValue("main_model_part", self.settings["sub_model_part_name"])
        process_parameters.AddValue("rayleigh_damping", self.settings["rayleigh_damping"])
        process_parameters.AddValue("assign_active_flag_node", self.settings["assign_active_flag_node"])
        process_parameters.AddValue("nodal_mass", self.settings["nodal_mass"])
        process_parameters.AddValue("nodal_inertia", self.settings["nodal_inertia"])
        process_parameters.AddValue("nodal_stiffness", self.settings["nodal_stiffness"])
        process_parameters.AddValue("nodal_rotational_stiffness", self.settings["nodal_rotational_stiffness"])
        process_parameters.AddValue("nodal_damping_ratio", self.settings["nodal_damping_ratio"])
        process_parameters.AddValue("nodal_rotational_damping_ratio", self.settings["nodal_rotational_damping_ratio"])
        process_parameters.AddValue("additional_dependence_variables", self.settings["additional_dependence_variables"])
        process_parameters.AddValue("interval", self.settings["interval"])
        self.assign_nodal_elements_to_nodes = StructuralMechanicsApplication.AssignNodalElementsToNodes(self.main_model_part, process_parameters)

    def ExecuteInitialize(self):
        self.assign_nodal_elements_to_nodes.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):
        self.assign_nodal_elements_to_nodes.ExecuteInitializeSolutionStep()
