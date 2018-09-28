from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics as KM

# Check that applications were imported in the main script
KM.CheckRegisteredApplications("StructuralMechanicsApplication")

import KratosMultiphysics.StructuralMechanicsApplication as SMA

def Factory(settings, Model):
    if(type(settings) != KM.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return SPRISMProcess(Model, settings["Parameters"])

# All the processes python processes should be derived from "Process"

class SPRISMProcess(KM.Process):
    """This class is used in order to compute some pre and post process on the SPRISM solid shell elements

    Only the member variables listed below should be accessed directly.

    Public member variables:
    Model -- the container of the different model parts.
    settings -- Kratos parameters containing the settings.
    """

    def __init__(self, Model, settings):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        Model -- the container of the different model parts.
        settings -- Kratos parameters containing solver settings.
        """

        # Settings string in json format
        default_parameters = KM.Parameters("""
        {
            "help"                           :"This class is used in order to compute some pre and post process on the SPRISM solid shell elements",
            "mesh_id"                        : 0,
            "model_part_name"                : "Structure",
            "preprocess_shell_to_solidshell" : false,
            "parameters_shell_to_solidshell" : {
                "element_name"                         : "SolidShellElementSprism3D6N",
                "new_constitutive_law_name"            : "LinearElastic3DLaw",
                "number_of_layers"                     : 1,
                "export_to_mdpa"                       : false,
                "output_name"                          : "output",
                "computing_model_part_name"            : "computing_domain",
                "create_submodelparts_external_layers" : false,
                "append_submodelparts_external_layers" : false,
                "initialize_elements"                  : false
            }
        }
        """)

        # Overwrite the default settings with user-provided parameters
        self.settings = settings
        self.settings.RecursivelyValidateAndAssignDefaults(default_parameters)

        # We define the model parts
        self.solid_shell_model_part = Model[self.settings["model_part_name"].GetString()]
        self.main_model_part = self.solid_shell_model_part.GetRootModelPart()

        # We create the process to compute the neighbours (should be run each time we recompute connectivity)
        self.sprism_neighbour_search = SMA.PrismNeighboursProcess(self.solid_shell_model_part)

        # We create the process to compute the thickness of the solid shell (post-process info)
        self.thickness_compute_process = SMA.SolidShellThickComputeProcess(self.solid_shell_model_part)

    def ExecuteInitialize(self):
        """ This method is executed at the begining to initialize the process

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # We preprocess from triangle shells to SPRISM solid-shells
        if (self.settings["preprocess_shell_to_solidshell"].GetBool() is True):
            parameters_shell_to_solidshell = KM.Parameters("""{}""")
            parameters_shell_to_solidshell.AddValue("element_name", self.settings["parameters_shell_to_solidshell"]["element_name"])
            parameters_shell_to_solidshell.AddValue("new_constitutive_law_name", self.settings["parameters_shell_to_solidshell"]["new_constitutive_law_name"])
            parameters_shell_to_solidshell.AddValue("model_part_name", self.settings["model_part_name"])
            parameters_shell_to_solidshell.AddValue("number_of_layers", self.settings["parameters_shell_to_solidshell"]["number_of_layers"])
            parameters_shell_to_solidshell.AddValue("export_to_mdpa", self.settings["parameters_shell_to_solidshell"]["export_to_mdpa"])
            parameters_shell_to_solidshell.AddValue("output_name", self.settings["parameters_shell_to_solidshell"]["output_name"])
            parameters_shell_to_solidshell.AddValue("computing_model_part_name", self.settings["parameters_shell_to_solidshell"]["computing_model_part_name"])
            parameters_shell_to_solidshell.AddValue("create_submodelparts_external_layers", self.settings["parameters_shell_to_solidshell"]["create_submodelparts_external_layers"])
            parameters_shell_to_solidshell.AddValue("append_submodelparts_external_layers", self.settings["parameters_shell_to_solidshell"]["append_submodelparts_external_layers"])
            parameters_shell_to_solidshell.AddValue("initialize_elements", self.settings["parameters_shell_to_solidshell"]["initialize_elements"])

            preprocess_shell_to_solidshell = SMA.TriangleShellToSolidShellProcess(self.main_model_part, parameters_shell_to_solidshell)
            preprocess_shell_to_solidshell.Execute()

        # We compute the neighbours
        self.sprism_neighbour_search.Execute()

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
        # We compute the neighbours if we have remeshed the problem
        if (self.main_model_part.Is(KM.MODIFIED)):
            self.sprism_neighbour_search.Execute()

    def ExecuteFinalizeSolutionStep(self):
        """ This method is executed in order to finalize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        # We compute the thickness of the solid shell element
        self.thickness_compute_process.Execute()

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
