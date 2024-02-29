import KratosMultiphysics
from KratosMultiphysics import MeshingApplication as Meshing

class ConvertLinearTetrahedraToQuadraticModeler(KratosMultiphysics.Modeler):
    """
    This class is a modeler in KratosMultiphysics that converts linear tetrahedral (tet) elements to quadratic tet elements within a specified model part.

    Attributes:
        model (KratosMultiphysics.Model): The Kratos model containing the model part to be modified.
        settings (KratosMultiphysics.Parameters): Settings for the conversion process.

    Methods:
        __init__(model, settings): Constructor for the class.
        SetupModelPart(): Performs the conversion from linear to quadratic tet elements.
        __GetDefaultSettings(): Provides default settings for the modeler.
        __str__(): String representation of the class.
        Factory(model, settings): Factory method to create an instance of the class.
    """

    def __init__(self, model, settings):
        """
        Constructor for the ConvertLinearTetrahedraToQuadraticModeler class.

        Args:
            model (KratosMultiphysics.Model): The Kratos model containing the model part to be modified.
            settings (KratosMultiphysics.Parameters): Configuration settings for the modeler.

        The constructor initializes the modeler with the provided model and settings, and validates and assigns default settings.
        """
        super().__init__(model, settings)
        settings.ValidateAndAssignDefaults(self.__GetDefaultSettings())
        self.model = model
        self.settings = settings

    def SetupModelPart(self):
        """
        Converts the specified model part's linear tet elements to quadratic tets.

        Raises:
            Exception: If the operation is attempted in a distributed (MPI) run or if the input elements are not linear tets (4 nodes).

        This method checks the model part for compatibility, sets necessary variables, and calls the utility function for conversion. It also handles the replacement of elements and conditions based on provided settings.
        """
        super().SetupModelPart()

        if KratosMultiphysics.IsDistributedRun():
            raise Exception("ConvertLinearTetrahedraToQuadraticModeler is not MPI-ready")

        model_part = self.model[self.settings["model_part_name"].GetString()]
        if model_part != model_part.GetRootModelPart():
            raise Exception(f"A root ModelPart is required, got SubModelPart {model_part.FullName()}")

        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(
            KratosMultiphysics.SPLIT_ELEMENT,
            True,
            model_part.Elements)
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(
            KratosMultiphysics.SPLIT_ELEMENT,
            True,
            model_part.Conditions)

        Meshing.LinearToQuadraticTetrahedraMeshConverter(model_part).LocalConvertLinearToQuadraticTetrahedraMesh(False, False)

        element_name = self.settings["element_name"].GetString()
        condition_name = self.settings["condition_name"].GetString()
        if element_name or condition_name:
            replace_settings = KratosMultiphysics.Parameters("{}")
            replace_settings.AddString("element_name", element_name)
            replace_settings.AddString("condition_name", condition_name)
            KratosMultiphysics.ReplaceElementsAndConditionsProcess(model_part, replace_settings).Execute()

    def __GetDefaultSettings(self):
        """
        Provides the default settings for the modeler.

        Returns:
            KratosMultiphysics.Parameters: The default settings as a Parameters object.

        These settings include the name of the model part to be modified and the name of the element to use for replacement, if specified.
        """
        return KratosMultiphysics.Parameters('''{
            "model_part_name" : "main_model_part",
            "element_name" : "",
            "condition_name": ""
        }''')

    def __str__(self):
        """
        String representation of the ConvertLinearTetrahedraToQuadraticModeler class.

        Returns:
            str: A string describing the class.
        """
        return "ConvertLinearTetrahedraToQuadraticModeler"
