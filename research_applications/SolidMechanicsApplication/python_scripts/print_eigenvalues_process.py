import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid

def Factory(settings, Model):
    if( not isinstance(settings,KratosMultiphysics.Parameters) ):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return PrintEigenvaluesProcess(Model, settings["Parameters"])

class PrintEigenvaluesProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings):

        default_settings = KratosMultiphysics.Parameters(
            """
            {
                "model_part_name" : "Structure",
                "variable_name"   : "EIGENVALUE_VECTOR"
            }
            """
        );

        settings.ValidateAndAssignDefaults(default_settings)

        KratosMultiphysics.Process.__init__(self)
        self.model_part = Model[settings["model_part_name"].GetString()]
        self.variable = getattr(KratosSolid, settings["variable_name"].GetString())

    def ExecuteFinalizeSolutionStep(self):

        eigen_values = [ev for ev in self.model_part.ProcessInfo[self.variable]]
        output_name = 'Eigenvalues.txt'
        with open(output_name, 'w') as output:
            output.write ("The Eigenvalues are:")
            output.write (str(eigen_values))
