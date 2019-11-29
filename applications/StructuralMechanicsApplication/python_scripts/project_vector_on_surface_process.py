from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ProjectVectorOnSurfaceProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ProjectVectorOnSurfaceProcess(KratosMultiphysics.Process):

    def __init__(self, Model, settings ):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        Model -- the container of the different model parts.
        settings -- Kratos parameters containing solver settings.
        """
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""{
            "model_part_name"  : "Structure",
            "echo_level"       : 0,
            "projection_type"  : "planar",
            "global_direction" : [1,0,0],
            "variable_name"    : "PLEASE_SPECIFY",
            "method_specific_settings" : { }
        }""")

        StructuralMechanicsApplication.ProjectVectorOnSurfaceUtility.Execute(Model[settings["model_part_name"].GetString()],settings)
