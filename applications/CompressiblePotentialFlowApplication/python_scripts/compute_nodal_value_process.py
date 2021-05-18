import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp

def Factory(settings, Model):
    if(not isinstance(settings, KratosMultiphysics.Parameters)):
        raise Exception(
            "expected input shall be a Parameters object, encapsulating a json string")
    return ComputeNodalValueProcess(Model, settings["Parameters"])

# all the processes python processes should be derived from "python_process"

class ComputeNodalValueProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings):
        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters(r'''{
            "model_part_name": "",
            "elemental_variables_list_to_project": []
        }''')

        settings.ValidateAndAssignDefaults(default_parameters)

        self.main_model_part = Model[settings["model_part_name"].GetString()]
        self.elemental_variables_list_to_project = settings["elemental_variables_list_to_project"].GetStringArray()

    def ExecuteFinalizeSolutionStep(self):
        if not self.elemental_variables_list_to_project:
            KratosMultiphysics.Logger.PrintWarning('ComputeNodalValueProcess', 'No elemental variables were specified to compute its nodal projection')
        else:
            CPFApp.ComputeNodalValueProcess(self.main_model_part, self.elemental_variables_list_to_project).Execute()


