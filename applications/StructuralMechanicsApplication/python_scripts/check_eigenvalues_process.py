# Importing the Kratos Library
import KratosMultiphysics
# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return CheckEigenvaluesProcess(Model, settings["Parameters"])

class CheckEigenvaluesProcess(KratosMultiphysics.Process, KratosUnittest.TestCase):
    def __init__(self, Model, settings):

        default_settings = KratosMultiphysics.Parameters(
            """
            {
                "help"            :"This process checks the solution obtained in a eigenvalue problem. Can be used to create tests",
                "model_part_name" : "Structure",
                "variable_name"   : "EIGENVALUE_VECTOR",
                "reference_values": "[1.,2.,3.]"
            }
            """
        );

        settings.ValidateAndAssignDefaults(default_settings)

        KratosMultiphysics.Process.__init__(self)
        self.model_part = Model[settings["model_part_name"].GetString()]
        self.variable = KratosMultiphysics.KratosGlobals.GetVariable(settings["variable_name"].GetString())
        self.reference_values = []
        reference_values = settings["reference_values"].GetString()
        for ev in reference_values.strip('[]').split(','):
            self.reference_values.append(float(ev))

    def ExecuteFinalizeSolutionStep(self):
        current_values = [ev for ev in self.model_part.ProcessInfo[self.variable]]
        for evs in zip(current_values,self.reference_values):
            self.assertAlmostEqual(evs[0],evs[1],5)
