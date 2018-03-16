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
    return PrintAfterFormfinding(Model,settings)

class PrintAfterFormfinding(KratosMultiphysics.Process):
    def __init__(self, Model, settings):

        default_settings = KratosMultiphysics.Parameters(
            """
            {
                "python_module"   : "print_after_formfinding_process",
                "kratos_module"   : "KratosMultiphysics.StructuralMechanicsApplication",
                "help"                  : "This process writes the mesh resulting from the formfinding in a .mdpa-file",
                "process_name"          : "PrintAfterFormfindingProcess",
                "Parameters"            : {
                    "model_part_name"   : "Structure"
                }
            }
            """
        );

        settings.ValidateAndAssignDefaults(default_settings)
        KratosMultiphysics.Process.__init__(self)
        model = Model[settings["Parameters"]["model_part_name"].GetString()]
        self.print_prestress = StructuralMechanicsApplication.FormfindingPrintUtility(model, settings)
                                                                              
    def ExecuteInitialize(self):
        pass
    
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
        self.print_prestress.PrintModelPart()
