from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("StructuralMechanicsApplication")

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication


def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return FormfindingIO(Model,settings)

class FormfindingIO(KratosMultiphysics.Process):
    """
    This class is responsible for the input and output of prestress and modelpart data for formfinding/ membrane analysis.
    """
    def __init__(self, Model, settings):

        default_settings = KratosMultiphysics.Parameters(
            """
            {
                "help"              : "This class is responsible for the input and output of prestress and modelpart data for formfinding/ membrane analysis.",
                "model_part_name"   : "Structure",
                "print_mdpa"        : false,
                "print_prestress"   : false,
                "read_prestress"    : false

            }
            """
        );
        self.settings = settings["Parameters"]
        self.settings.ValidateAndAssignDefaults(default_settings)
        KratosMultiphysics.Process.__init__(self)
        model = Model[self.settings["model_part_name"].GetString()]
        self.print_mdpa = self.settings["print_mdpa"].GetBool()
        self.print_prestress = self.settings["print_prestress"].GetBool()
        self.read_prestress = self.settings["read_prestress"].GetBool()
        self.formfinding_io = StructuralMechanicsApplication.FormfindingIOUtility(model, settings)

    def ExecuteInitialize(self):
        if (self.read_prestress):
            self.formfinding_io.ReadPrestressData()
            KratosMultiphysics.Logger.PrintInfo("FormfindingIO", "Read prestress")

    def ExecuteBeforeSolutionLoop(self):
        pass

    def ExecuteInitializeSolutionStep(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        pass

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        pass


    def ExecuteFinalize(self):
        if (self.print_mdpa):
            self.formfinding_io.PrintModelPart()
            KratosMultiphysics.Logger.PrintInfo("FormfindingIO", "Print mdpa")
        if (self.print_prestress):
            self.formfinding_io.PrintPrestressData()
            KratosMultiphysics.Logger.PrintInfo("FormfindingIO", "Print prestress")
