import KratosMultiphysics
from pathlib import Path

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception(
            "expected input shall be a Parameters object, encapsulating a json string"
        )
    return SetUpPreStressedOrientedCompositeMaterials(Model, settings["Parameters"])

class SetUpPreStressedOrientedCompositeMaterials(KratosMultiphysics.Process):

    """This process sets a proper orientation of the local axes of the elements intersected by line elements (steel tendons). Besides it also computes and sets a volumetric participation of steel within the concrete FE as well as an indicated pre-stressing strain. It also creates a submodelpart for each steel tendon intersected FE.

    Format of the file:
    # Intersection points (D is the diameter of the tendon and Ep the imposed pre-stressing strain):

    Begin Tendon - Hexahedra intersection: Tendon_Inf     D=1.0e-3    Ep=0.00001
	    4807	0.01000	     0.05000	   0.01500	    0.02000	     0.05000   0.01500
        ...
    End Tendon - Hexahedra intersection

    # FE inside each tendon:

    Begin Tendon - Hexahedra: Tendon_Inf
        4807    4812    4817    4822    4827    4832    4837    4842    4847    4852   ...
    End Tendon - Hexahedra


    Public member variables:
    Model -- the container of the different model parts.
    settings -- Kratos parameters containing solver settings.
    """

    def __init__(self, Model, settings):
        """The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        Model -- the container of the different model parts.
        settings -- Kratos parameters containing solver settings.
        """
        KratosMultiphysics.Process.__init__(self)

        # The value can be a double or a string (function)
        default_settings = KratosMultiphysics.Parameters(
            """
        {
            "help" : "This sets the initial conditions in terms of imposed pre-stressing strain, local axes and % participation of fiber",
            "model_part_name" : "please_specify_model_part_name",
            "intersection_file_name" : "please_include_directory_and_full_name_with_extension"
        }
        """
        )
        settings.ValidateAndAssignDefaults(default_settings)

        self.intersection_file_name = settings["intersection_file_name"].GetString()
        self.model_part = Model[settings["model_part_name"].GetString()]

    def ExecuteInitializeSolutionStep(self):
        """This method is executed in order to initialize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        KratosMultiphysics.Logger.PrintInfo("SetUpPreStressedOrientedCompositeMaterials ", "Reading intersections file " + self.intersection_file_name + "...")
        intersections_file = open(self.intersection_file_name, "r")












        intersections_file.close()