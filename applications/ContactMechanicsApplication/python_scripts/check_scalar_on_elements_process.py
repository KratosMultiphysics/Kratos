import KratosMultiphysics
# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

def Factory(settings, Model):
    if( not isinstance(settings,KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return CheckScalarOnElementsProcess(Model, settings["Parameters"])


## All the processes python should be derived from "Process"
class CheckScalarOnElementsProcess(KratosMultiphysics.Process, KratosUnittest.TestCase):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
            {
                "mesh_id"         : 0,
                "model_part_name" : "please_specify_model_part_name",
                "variable_name"   : "SPECIFY_VARIABLE_NAME",
                "tolerance_rank"  : 3
            }
            """
            )


        settings.ValidateAndAssignDefaults(default_settings)

        self.model    = Model
        self.settings = settings
        self.variable = getattr(KratosMultiphysics, settings["variable_name"].GetString())
        self.tol = settings["tolerance_rank"].GetInt()

    def ExecuteInitialize(self):
        self.model_part = self.model[self.settings["model_part_name"].GetString()]
        self.mesh = self.model_part.GetMesh(self.settings["mesh_id"].GetInt())

    def ExecuteInitializeSolutionStep(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        FirstElement = self.mesh.Elements[1]

        SomeProcessInfo = KratosMultiphysics.ProcessInfo()
        VariableValue = FirstElement.CalculateOnIntegrationPoints( self.variable, SomeProcessInfo)
        VariableValue = VariableValue[0]


        for elem in self.mesh.Elements:
            ThisValue = elem.CalculateOnIntegrationPoints( self.variable, SomeProcessInfo)
            ThisValue = ThisValue[0]
            for i in range(0, VariableValue.Size1()):
                for j in range(0, VariableValue.Size2()):
                    self.assertAlmostEqual( VariableValue[i,j], ThisValue[i,j], self.tol)
