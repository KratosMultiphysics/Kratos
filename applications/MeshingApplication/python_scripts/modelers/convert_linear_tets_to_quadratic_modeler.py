import KratosMultiphysics
from KratosMultiphysics import MeshingApplication as Meshing

class ConvertLinearTetsToQuadraticModeler(KratosMultiphysics.Modeler):
    def __init__(self, model, settings):
        super().__init__(model, settings)

        settings.ValidateAndAssignDefaults(self.__GetDefaultSettings())
        self.model = model
        self.settings=settings

    # Converts the input model part from linear tets to quadratic
    def SetupModelPart(self):
        super().SetupModelPart()

        if(KratosMultiphysics.IsDistributedRun()):
            raise Exception("ConvertLinearTetsToQuadraticModeler is not MPI-ready")

        model_part = self.model[self.settings["model_part_name"].GetString()]
        #quick check
        for elem in model_part.Elements:
            if elem.GetGeometry().PointsNumber()!=4:
                msg = "Only linear tets are accepted as input (4 nodes) \n"
                raise Exception(msg)

        #flag elems and conds
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(
            KratosMultiphysics.SPLIT_ELEMENT,
            True,
            model_part.Elements)
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(
            KratosMultiphysics.SPLIT_ELEMENT,
            True,
            model_part.Conditions)

        #calling utility to replace the tets:
        Meshing.LinearToQuadraticTetrahedraMeshConverter(model_part.GetRootModelPart()).LocalConvertLinearToQuadraticTetrahedraMesh(False,False)

        element_name = self.settings["element_name"].GetString()

        if element_name:
            replace_settings = KratosMultiphysics.Parameters("{}")
            replace_settings.AddString("element_name", element_name)
            replace_settings.AddString("condition_name", "")
            KratosMultiphysics.ReplaceElementsAndConditionsProcess(model_part,replace_settings).Execute()

    def __GetDefaultSettings(self):
        default_settings = KratosMultiphysics.Parameters('''{
            "model_part_name" : "main_model_part",
            "element_name" : ""
        }''')
        return default_settings

    def __str__(self):
        return "ConvertLinearTetsToQuadraticModeler"

def Factory(model, settings):
    return ConvertLinearTetsToQuadraticModeler(model, settings)