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

        default_settings = KratosMultiphysics.Parameters("""
        {
            "main_model_part_name"        : "Structure",
            "echo_level"      : 1,
            "method_specific_settings" : {
                "projection_type"        : "planar",
                "global_fiber_direction" : [1,0,0],
                "local_variable_name"    : "LOCAL_MATERIAL_AXIS_1"
            }
        }
        """)

        # The main model part
        main_model_part = Model[settings["main_model_part_name"].GetString()]
        StructuralMechanicsApplication.ProjectVectorOnSurfaceUtility.Execute(main_model_part,settings)
